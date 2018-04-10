#!/usr/bin/env python

#===============================================================================
# Pileup 2 nucleotide rate
# written by : Ofer Isakov, Shomron Lab, Tel Aviv University
#
# takes in a pileup file (produced by samtools) and returns the nucleotide
# rate in each position, ignores start and end markers and counts reverse
# strand nucleotides as sense
#===============================================================================
VERSION = '0.51'
import optparse
import re
import sys
import time
import scipy
import numpy
from numpy import array,log
from scipy.stats import binom
from scipy.stats import mannwhitneyu
from scipy.stats import fisher_exact
from scipy.stats import poisson
from Queue import Queue
import threading

mean = lambda(x): sum(x) *1./len(x)

def getFileLen(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

## compute the log of the likelihood for a parameter f with data yis,pis
## works with arrays for optimal performance
def log_lik_f(f,yis,one_minus_yis,pis,one_minus_pis):
    f = max(0,min(f,1))
    q =  (1-f)*pis + one_minus_pis*f
    a = log(yis*q + one_minus_yis*(1-q))
    return -numpy.sum(a)

def log_lik_f_backup(f,yis,pis):
    f = max(0,min(f,1))
    a = log( (1-f)*pis + (1-pis)*f)
    b = log( (1-f)*(1-pis) + f*pis)
    return -sum(yis*a + (1-yis)*b)

## find the optimal f for this data using numerical optimization
def find_f(yis,pis):
    omy = 1-yis
    omp = 1-pis
    if len(yis)==0:
        return 0 
    if mean(yis)==0 or mean(yis)==1:
        return mean(yis)
    wrap = lambda f: log_lik_f(f,yis,omy,pis,omp)
    res = scipy.optimize.brute(wrap,[[0,1]])[0]
    return res

## find a confidence interval by finding the points [bottom,top] which
## are furthest from the MLE f, for which -2(log likelihood diff) is < 3.84
## this is an asymptotically correct CI.
def find_ci(f,yis,pis):
    omy = 1-yis
    omp = 1-pis
    step_size = max(1/10000,f/100)
    max_log_lik = log_lik_f(f,yis,omy,pis,omp)
    # find bottom -
        # our goal is to find the point for which -2*log lik diff is 3.84.
        # we do so by defining a function that returns the squared distance from
        # 3.84 and minimizing it.
    def wrap_bottom(b):
        b = max(0,min(b,f))
        b_log_lik = log_lik_f(b,yis,omy,pis,omp)
        return ((-2*(max_log_lik - b_log_lik)) - 3.84)**2

    bottom = scipy.optimize.brute(wrap_bottom,[[0,f]])[0]
    
    # find top
    def wrap_top(t):
        t = min(1,max(f,t))
        t_log_lik = log_lik_f(t,yis,omy,pis,omp)
        return ((-2*(max_log_lik - t_log_lik)) - 3.84)**2

    top = scipy.optimize.brute(wrap_top,[[f,1]])[0]

    return [bottom,top]      

class Feature:
    def __init__(self,bedLine):
        self.featureLine = bedLine.strip()
        self.featureArgs = bedLine.strip().split('\t')
        self.region = self.featureArgs[0]
        self.start = int(self.featureArgs[1])
        self.end = int(self.featureArgs[2])
        self.name = self.region
        self.cdsStartSites=None
        self.cdsEndSites=None
        if len(self.featureArgs)>3:
            name = self.featureArgs[3].strip()
            if name:
                self.name = self.featureArgs[3].strip()
        if len(self.featureArgs)>4:
            cdsStartSites = self.featureArgs[4].strip()
            if cdsStartSites:
                self.cdsStartSites = [int(startSite) for startSite in cdsStartSites.split(',')]
            cdsEndSites = self.featureArgs[5].strip()
            if cdsEndSites:
                self.cdsEndSites = [int(endSite) for endSite in cdsEndSites.split(',')]
            if cdsStartSites:
                for i,cdsStartSite in enumerate(self.cdsStartSites[1:]):
                    if cdsStartSite<self.cdsEndSites[i]:
                        raise Exception('Feature error: there are 2 start sites found before the same end site:\n%s\n'%bedLine)
                        
                    
    def getCDS(self,fullSeq):            
        cds = ''
        if self.cdsStartSites:
            for i,cdsstart in enumerate(self.cdsStartSites):
                cdsend = self.cdsEndSites[i]
                cds+=fullSeq[cdsstart-1:cdsend]
        else:
            cds+=fullSeq[self.start-1:self.end]
        return cds
    
    def isPosInFeature(self,region,pos):
        if region==self.region:
            if self.cdsStartSites:
                for i,cdsstart in enumerate(self.cdsStartSites):
                    cdsend = self.cdsEndSites[i]
                    if int(pos)>=cdsstart and int(pos)<=cdsend:
                        return True
                else:
                    return False
            if int(pos)>=self.start and int(pos)<=self.end:
                return True
        return False
    
def parseFeatureBed(bedFile):
    bedFile = open(bedFile,'r').xreadlines()
    features = []
    for line in bedFile:
        if line.strip():
            feature = Feature(line)
            features.append(feature)
    return features

def convert2phred(qualLine,type='sanger'):
    # get a string of ascii per-base qualities and return a list of the numerical
    # error probabilities 
    probs = []
    qualLine = re.sub('[\s\n]','',qualLine)
    for q in qualLine:
        asc2num = ord(q)
        if type=='illumina':
            phred = asc2num-64
        if type=='sanger':
            phred=asc2num-33
        
        prob = 10**(float(phred)/(-10))
        if prob>1:
            raise Exception('probability higher than 1 (%s=%s), probably %s is not the right quality type'%(q,prob,type))
        probs.append(prob)
    return probs


def testQualBias(baseRef,qualsDict,alleles):

    # calculate qual bias
    qualsBiasDict = dict((allele,[]) for allele in alleles if (allele!=baseRef and len(qualsDict[allele])!=0))
#    print 'DEBUG: REF: %s QUAL AVERAGE: %s'%(baseRef,','.join(['%s:%s:%s'%(allele,len(qualsDict[allele]),str(sum(qualsDict[allele])/len(qualsDict[allele]))) for allele in qualsDict.keys() if len(qualsDict[allele])!=0 ]))
    baseRefQualAvr = 0
    if len(qualsDict[baseRef])>0:
        baseRefQualAvr = sum(qualsDict[baseRef])/len(qualsDict[baseRef])
        
    for varAllele in qualsBiasDict.keys():
        if len(set(qualsDict[baseRef]+qualsDict[varAllele]))==1:
            # only one number in all of them, no point in test
            qualsBiasDict[varAllele]=1
        elif sum(qualsDict[varAllele])/len(qualsDict[varAllele])<baseRefQualAvr:
            # if the variant average error rate is lower than the reference, no point in test
            qualsBiasDict[varAllele]=1
        else:
            qualsBiasDict[varAllele] = mannwhitneyu(qualsDict[varAllele],qualsDict[baseRef])[1]
#    print 'DEBUG: BIAS pvalue : %s'%(','.join(['%s:%s'%(allele,qualsBiasDict[allele]) for allele in qualsBiasDict.keys()]))

    return qualsBiasDict

def testStrandBias(baseRef,strandDict,alleles):

    # calculate strand bias
    strandBiasDict = dict((allele,[]) for allele in alleles if (allele!=baseRef))
    for varAllele in strandBiasDict.keys():
        totalReads = sum(strandDict[varAllele].values())
        baseBinom = binom(totalReads,0.5)
        strandBiasDict[varAllele] = baseBinom.cdf(min(strandDict[varAllele].values()))
#     if 'T' in strandDict:
#         print 'DEBG:StrandBias : %s : %s\n'%(strandDict['T'],strandBiasDict)
#     print 'DEBUG: STRAND BIAS pvalue : %s\n%s'%(','.join(['%s:%s'%(allele,strandBiasDict[allele]) for allele in strandBiasDict.keys()]),strandDict)
    return strandBiasDict

def testStrandBiasFisher(baseRef,strandDict,alleles):

    # calculate strand bias
    strandBiasDict = dict((allele,[]) for allele in alleles if (allele!=baseRef))
    refSense = strandDict[baseRef]['sense']
    refAnti = strandDict[baseRef]['anti']
    refs = [strandDict[baseRef]['sense'],strandDict[baseRef]['anti']]
    for varAllele in strandBiasDict.keys():
        vars = [strandDict[varAllele]['sense'],strandDict[varAllele]['anti']]
        if refs==[0,0] or vars==[0,0] or (refs[0]==0 and vars[0]==0) or(refs[1]==0 and vars[1]==0):
            strandBiasDict[varAllele]=1
            continue
        fisherOdd,fisherPval = fisher_exact([refs,vars])
        strandBiasDict[varAllele]=fisherPval
#         if fisherPval<0.05:
#             print '%s : %s , %s: %s'%(baseRef,refs,varAllele,vars)
#             print 'DEBG:StrandBias : %s : odd:%s pval:%s\n'%(strandDict,fisherOdd,fisherPval)
#     print 'DEBUG: STRAND BIAS pvalue : %s\n%s'%(','.join(['%s:%s'%(allele,strandBiasDict[allele]) for allele in strandBiasDict.keys()]),strandDict)
    return strandBiasDict

def getAlleleStats(nucs,qualsDict):
    allQuals = []
    alleles = sorted(qualsDict.keys())
    for allele in alleles:
        allQuals+=qualsDict[allele]
    alleleStats = dict([(allele,{}) for allele in alleles])
    # the total number of reads in this position (n)
    totalBases = sum(nucs.values())
    
    # for each of the possible alleles
    for allele in alleles:
        numOallele= len(qualsDict[allele])
        # the original allele rate
        alleleRate = float(numOallele)/totalBases
        
        cis = []
        fis = []
        pis = []
        yis = []
        # go through the qualsDict dictionary and collect cis and fis
        for qallele in qualsDict.keys():
            yi=0
            if qallele==allele:
                yi=1
            qalleleProbs = qualsDict[qallele]
            for pi in qalleleProbs:
                pis.append(pi)
                yis.append(yi)
                if pi==1:
                    ci = 0
                else:
                    ci = (yi-pi)/(pi*(1-pi))**0.5
                cis.append(ci)

        # calculate Z
        z = sum(cis)/totalBases**0.5
        # calculate survival function from z (two-tailed pval)
        pVal = scipy.stats.norm.sf(z)
        
        yis = array(yis)
        pis = array(pis)
        ml_f = find_f(yis,pis)
        ci = find_ci(ml_f,yis,pis)
        wilks = -2*(log_lik_f_backup(ml_f,yis,pis) - log_lik_f_backup(0,yis,pis)) 
        pVal = 1 - scipy.stats.chi2.cdf(wilks, 1)
        
        alleleStats[allele]['pVal'] = pVal
        alleleStats[allele]['alleleRate']=alleleRate
        alleleStats[allele]['z']=z
        alleleStats[allele]['estimatedRate']=ml_f
        alleleStats[allele]['conf95']=ci[1]
        alleleStats[allele]['conf5']=ci[0]
    return alleleStats
            
def isBaseWithSNPs(nucs):
    numOalleles = len(nucs.keys())
    return (nucs.values().count(0)<(numOalleles-2) or (nucs.values().count(0)<(numOalleles-3) and nucs['N']!=0))     
def isBaseWithIndels(indels):
    return indels.keys()

def getAlleleRate(alleleCount,baseCoverage):
    return float(alleleCount)/baseCoverage

def getRMSD(nucs,baseRef,baseCov):
    numOalleles = len(nucs.keys())
    sumOdeviations = 0
    if baseCov==0:
        return 0
    for allele in nucs.keys():
        if allele==baseRef:
            expected=1
        else:
            expected=0
        alleleRate = getAlleleRate(nucs[allele], baseCov)
        alleleVarDeviation = (expected-alleleRate)**2
        sumOdeviations+=alleleVarDeviation
    return (sumOdeviations/numOalleles)**0.5

def getShannon(perAlleleStats,baseCov):  
    shannon = 0.0     
    if baseCov==0 or not perAlleleStats:
        # if there is no perAlleleStats it means that there was no analysis done because there are no variant alleles
        return shannon
    for allele in perAlleleStats.keys():
        alleleRate = perAlleleStats[allele]['estimatedRate']
        if alleleRate==0:
            continue
        shannon+=alleleRate*log(alleleRate)
    return -shannon

def getNucRate(baseRef,basePile,baseQual,baseCov,alleles,PVAL):
    nucs = dict((allele,0) for allele in alleles)
    indels = {}
    SNP=False
    INDEL =False
    qualsDict = dict((allele,[]) for allele in alleles)
    strandDict = dict((allele,{'sense':0,'anti':0}) for allele in alleles)
    # split the pileup line into single reads
    iterNuc = iter(basePile)
    
    readNum = 0
    # for read nucleotide in the pileup
    for nuc in iterNuc:    
        #if the read supports an indel 
        if re.match('[+-]', nuc):
           # readNum+=1
            # get the indel sign (+/-)
            indel = nuc
            # get the indel size
            sizeNuc = iterNuc.next()
            indelSize =''
            
            while re.match('\d',sizeNuc):
                indelSize+=sizeNuc
                sizeNuc = iterNuc.next()
            
            indel+=sizeNuc    
            # get the indel allele
            for i in range(int(indelSize)-1):
                indelNuc = iterNuc.next()
                indel+=indelNuc
            
            # add the indel to the indels dict
            try:
                indels[indel]+=1
            except KeyError:
                indels[indel]=1
            continue
        
        #if the same base as the reference
        elif nuc=='.' or nuc==',':
            base = baseRef.upper()
            if nuc=='.':
                strandDict[base]['sense']+=1
            else:
                strandDict[base]['anti']+=1
            qualsDict[base].append(baseQual[readNum])
            readNum+=1
        #if SNP
        elif nuc.upper() in nucs.keys():
            #test strand
            if nuc.upper()==nuc:
                strandDict[nuc.upper()]['sense']+=1
            else:
                strandDict[nuc.upper()]['anti']+=1
            base = nuc.upper()
            qualsDict[nuc.upper()].append(baseQual[readNum])
            readNum+=1
            
        #if ^ or $ (start or end of read)
        else:
            if nuc=='^':
                nuc = iterNuc.next()
            continue
        # add the read base to the dictionary
        nucs[base]+=1
    alleleStats=None
    strandBiasDict=None
    if re.search('[ACGT\*]',basePile,re.IGNORECASE):
        alleleStats = getAlleleStats(nucs, qualsDict)
        #qualsBiasDict = testQualBias(baseRef.upper(), qualsDict,alleles)
        
#         strandBiasDict = testStrandBias(baseRef.upper(), strandDict, alleles)
        strandBiasDict = testStrandBiasFisher(baseRef.upper(), strandDict, [sigAllele for sigAllele in alleleStats if alleleStats[sigAllele]['pVal']<PVAL])
            
    
    totalBases = sum(nucs.values())

    return nucs,indels,alleleStats,strandBiasDict

def collectBaseArgs(alleles,base,qualType,Q,PVAL):
    baseDict = {}
    baseArgs = base.split('\t')
    baseRegion = baseArgs[0]
    basePos = int(baseArgs[1])
    baseRef = baseArgs[2].upper()
    baseCov = int(baseArgs[3])
    basePile = baseArgs[4]
    baseQual = baseArgs[5]
    probs = convert2phred(baseQual, qualType)
    # get the nucleotide rates for the base
    baseDict['nucs'],baseDict['indels'],baseDict['alleleStats'],baseDict['strandBias'] = getNucRate(baseRef, basePile,probs,baseCov,alleles,PVAL)
    baseDict['baseRef']=baseRef
    baseDict['baseCov']=baseCov
    baseDict['nonVar']=False
    Q.put([baseRegion,basePos,baseDict])
    Q.task_done()
    return 

def collectBaseArgsNoThread(alleles,base,qualType,PVAL):
    baseDict = {}
    baseArgs = base.split('\t')
    baseRegion = baseArgs[0]
    basePos = int(baseArgs[1])
    baseRef = baseArgs[2].upper()
    baseCov = int(baseArgs[3])
    basePile = baseArgs[4]
    baseQual = baseArgs[5]
    probs = convert2phred(baseQual, qualType)
    # get the nucleotide rates for the base
    baseDict['nucs'],baseDict['indels'],baseDict['alleleStats'],baseDict['strandBias'] = getNucRate(baseRef, basePile,probs,baseCov,alleles,PVAL)
    baseDict['baseRef']=baseRef
    baseDict['baseCov']=baseCov
    baseDict['nonVar']=False
    return baseRegion,basePos,baseDict

def getNucsAndPvals(inputFile,alleles,qualType,totalBases,PVAL,THREAD):
    # take in a pileup file, go over it base by base
    inputPile = open(inputFile).xreadlines()
    nonVarPos = {}
    posDict = {}
    regions = []
    i=0
    if THREAD==1:
        print 'VIVAN: Number of threads selected = 1, will run without threads\n'
        totalStartTime = time.time()
        for base in inputPile:
            # count the number of bases
            i+=1
            baseStartTime = time.time()
            region,pos,baseDict = collectBaseArgsNoThread(alleles, base, qualType, PVAL)
            if region not in posDict:
                posDict[region]={}
                regions.append(region)
            posDict[region][pos]=baseDict
            print '\rfinished %s/%s (base: %s; total: %s)'%(i,totalBases,time.time()-baseStartTime,time.time()-totalStartTime),
            sys.stdout.flush()
    else:
        tCount = 0
        Q = Queue()
        totalThreadTime = 0
        for base in inputPile:
            # count the number of bases
            i+=1
            # collect the base arguments
            tCount+=1
            baseThreads=[]
            bT = threading.Thread(target=collectBaseArgs,args=[alleles,base,qualType,Q,PVAL])
            bT.daemon = True
            bT.start()
            time.sleep(0.001)   
            baseThreads.append(bT)
            if tCount==THREAD or i==totalBases:
                threadStart = time.time()
                for baseThread in baseThreads:
                    baseThread.join()
                for j in range(tCount):
                   # print j
                    region,pos,baseDict = Q.get()
                    #print baseDict
                  #  print 'collected base %s'%pos
                    if region not in posDict:
                        posDict[region]={}
                        regions.append(region)
                    posDict[region][pos]=baseDict
    #                print posDict[region].keys()
    
                            
                tCount=0
                baseThreads = []
                threadTime = time.time()-threadStart
                totalThreadTime +=threadTime
                print 'finished %s/%s (%.2f%%) thread time : %s ; total thread time : %s'%(i,totalBases,(float(i)*100/totalBases),threadTime,totalThreadTime)
            
    return posDict,regions
        
def getBenjaminiHochberg(pvals,PVAL=0.05):
    
    sortedPvals = sorted(pvals)
    for i,pval in enumerate(sortedPvals):
        if pval>((i+1)*PVAL)/len(pvals):
            if i==0:
                # this means that the lowest p-value is still not 
                # significant after correction. and should return p-value
                # where no value passes:
                return 0
            return sortedPvals[i-1]
        

def getAllPvals(posDict):
    pvals = []
    for region in posDict.keys():
        for pos in posDict[region].keys():
            if posDict[region][pos]['alleleStats']:
                pvals=pvals+[posDict[region][pos]['alleleStats'][nuc]['pVal'] for nuc in posDict[region][pos]['alleleStats'].keys() if nuc!=posDict[region][pos]['baseRef']]
    return pvals

def getAllStrandPvals(posDict):
    pvals = []
    for region in posDict.keys():
        for pos in posDict[region].keys():
            if posDict[region][pos]['strandBias']:
                pvals=pvals+[posDict[region][pos]['strandBias'][nuc] for nuc in posDict[region][pos]['strandBias'].keys() if nuc!=posDict[region][pos]['baseRef']]
    return pvals

def getPosFeatures(region,pos,features):
    posFeats = []
    for feature in features:
        if feature.isPosInFeature(region,pos):
            posFeats.append(feature.name)
    if len(posFeats)==0:
        posFeats=['.']
    return posFeats
    
def stats():
    print 'Pileup2NucleotideRate v%s'%VERSION
    numOPosWithSNP = 0
    numOPosWithIndel = 0
    numOPosWithRef = 0     
    posDict = {}
    startTime = time.time()
    #===========================================================================
    # Input
    #===========================================================================
    parser = optparse.OptionParser()
    parser.add_option('-i','--inputPile')
    parser.add_option('-q','--qualType',default='sanger')
    parser.add_option('-p','--pval',default=0.05)
    parser.add_option('-N','--disregardN',default=False,action='store_true')
    parser.add_option('-f','--featuresFile',default=None)
    parser.add_option('-t','--threads',default=3,help='the number of threads you wish to use')
    
    args = parser.parse_args(sys.argv)[0]
    
    qualType = args.qualType.strip()
    pval = float(args.pval)
    inputPrefix = re.match('[^\.]+', args.inputPile).group(0)
    outFile = open(inputPrefix+'_nucleotideRate.csv','w')
    DISREGARD_N = args.disregardN
    THREAD = int(args.threads)
    
    featuresFile = args.featuresFile
    if featuresFile:
        features = parseFeatureBed(featuresFile)
    
    totalNumObases = getFileLen(args.inputPile)
    
    if DISREGARD_N:
        alleles = ['A','T','C','G','*'] 
    else:
        alleles = ['A','T','C','G','N','*'] 
        
    print 'there are %s bases in %s'%(totalNumObases,args.inputPile)
    
    numOvarsPass = 0
    
    headers = ['#region',
               'position',
               'reference',
               'coverage',
               '%s'%'\t'.join(alleles),
               '%s'%'\t'.join(['raw_%s_rate'%allele for allele in alleles]),
               '%s'%'\t'.join(['%s_rate'%allele for allele in alleles]),
               'significantVars',
               'strandBiasVars',
               'RMSD',
               'ShannonEntropy',
               'PValues']
    
    if featuresFile:
        headers.append('features')

    outFile.write('%s\n'%'\t'.join(headers))
        
    posDict,regions = getNucsAndPvals(args.inputPile, alleles,qualType,totalNumObases,pval,THREAD) 
    allPvals = getAllPvals(posDict)
#     allStrandPvals = getAllStrandPvals(posDict)
    collected = 0
    for region in posDict.keys():
        collected +=len(posDict[region].keys())

    print 'collected %s positions'%collected
    bh_corrected_pval_threshold = getBenjaminiHochberg(allPvals,pval)
    print 'Corrected Pval : %s'%bh_corrected_pval_threshold
    if len(allPvals)>0:
        print 'Strand Pval : %s'%(pval/len(allPvals))
    for region in regions:
        
        for pos in sorted(posDict[region].keys()):
            # check if the position is a non-variant (has no alternate alleles)
            baseRef = posDict[region][pos]['baseRef']
            baseCov = posDict[region][pos]['baseCov']
            nucs = posDict[region][pos]['nucs']
            indels = posDict[region][pos]['indels']
            baseRegion = region
#            qualsBias = posDict[region][pos]['qualBias']
            strandBias = posDict[region][pos]['strandBias']
            perAlleleStats = posDict[region][pos]['alleleStats']
          #  baseRegion = posDict[region][pos]['baseRegion']
            
            if featuresFile:
                posFeatures = getPosFeatures(region,pos,features)
                posFeaturesString = ';'.join(posFeatures)
                
            # if there are reads supporting indels in this base, prepare the indel output in the form of:
            # indel <tab> number of supporting reads
            indelsString = '';
            if isBaseWithIndels(indels):
                numOPosWithIndel+=1
                indelItems =  ["%s\t%s"%(indel,numOsupportingReads) for  indel,numOsupportingReads in indels.iteritems()]
                indelsString = '\t'.join(indelItems)
            
            # if there are no reads supporting indels in this position, check if there was any SNP, if so, count
            # it as a location with a SNP
            elif isBaseWithSNPs(nucs):
                numOPosWithSNP+=1
            # if not, count it as same as reference position
            else:
                numOPosWithRef+=1
            # return the output in the nucleotide number of reads format
            
            counts = []
            rawRatios = []
            estRatios = []
            for allele in alleles:
                counts.append(str(nucs[allele]))
                if nucs[allele]==0:
                    rawRatios.append('0')
                    estRatios.append('0;(0,0)')
                elif not perAlleleStats:
                    rawRatios.append('1')
                    estRatios.append('1;(1,1)')
                else:
                    rawRatios.append(str(perAlleleStats[allele]['alleleRate']))
                    estRatios.append('%s;(%s,%s)'%(perAlleleStats[allele]['estimatedRate'],perAlleleStats[allele]['conf5'],perAlleleStats[allele]['conf95']))
            
            baseCounts = '\t'.join(counts)
            rawRatios = '\t'.join(rawRatios)
            estRatios = '\t'.join(estRatios)
            
            RMSD = getRMSD(nucs, baseRef, baseCov)
            shannonEntropy = getShannon(perAlleleStats, baseCov)
            
            significantVars = []
            strandBiasVars = []
            pvals = ['%s:1'%allele for allele in sorted(alleles) if allele!=baseRef]
            if perAlleleStats:
                pvals = ['%s:%s'%(nuc,perAlleleStats[nuc]['pVal']) for nuc in sorted(perAlleleStats.keys()) if nuc!=baseRef]
                significantVars = ['%s'%(nuc) for nuc in sorted(perAlleleStats.keys()) if (perAlleleStats[nuc]['pVal']<=bh_corrected_pval_threshold and nuc!=baseRef)]
                strandBiasVars = ['%s'%(nuc) for nuc in significantVars if strandBias[nuc]>0.05/len(allPvals)]
                
            if len(significantVars)>0:
                significantVars = ';'.join(significantVars)
            else:
                significantVars = '.'
            if len(strandBiasVars)>0:
                strandVars = ';'.join(strandBiasVars)
            else:
                strandVars = '.'
            if len(pvals)>0:
                pvals = ';'.join(pvals)
            else:
                pvals='.'
          
            outputLineArgs = [baseRegion,
                              pos,
                              baseRef,
                              baseCov,
                              baseCounts,
                              rawRatios,
                              estRatios,
                              significantVars,
                              strandVars,
                              RMSD,
                              shannonEntropy,
                              pvals] 
            
            if featuresFile:
                outputLineArgs.append(posFeaturesString)      
            outputLineArgs.append(indelsString)
            
            outFile.write('%s\n'%'\t'.join([str(arg) for arg in outputLineArgs]))
            
            # Sanity check to see that the number of nucleotides equals the overall coverage
            sumNucs = sum([nucs[allele] for allele in alleles])
            if int(baseCov)!=sumNucs:
                raise Exception("error in position %s, the number of bases counted (%s) is not equal to the base coverage (%s)"%(pos,sumNucs,baseCov))
               
    print "out of %s bases:\nthere were %s bases without variants, %s bases with SNPs and %s positions with Indels\n"%(totalNumObases,numOPosWithRef,numOPosWithSNP,numOPosWithIndel)
    print 'done.. (%s)\n'%(time.time()-startTime)    
if __name__=='__main__':
    stats()
        
