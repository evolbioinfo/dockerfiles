#!/usr/bin/env python

import optparse
import re
import sys
import time
import os

VERSION = 0.7
#===============================================================================
# Produce Variance Metrics 
#===============================================================================
print 'VIVAN: ProduceVarMetrics v%s'%VERSION

def getHeaderArgs(nucFile):
    header = open(nucFile,'r').readline()
    headerArgs = [arg.strip() for arg in header.split('\t')]
    return headerArgs

def parseAnnotationFile(annotationFileName):
    print 'Parsing %s for metrics production..'%annotationFileName
    annotationFile = open(annotationFileName,'r')
    annotationsDict = {}
    for line in annotationFile:
        lineArgs = line.strip().split('\t')
        if line.startswith('#'):
            headers = lineArgs
            continue
        region = lineArgs[0]
        if region not in annotationsDict:
            annotationsDict[region]={}
        feature = lineArgs[1]
        position = lineArgs[2]
        alternate = lineArgs[headers.index('alternate')]
        if feature not in annotationsDict[region]:
            annotationsDict[region][feature]={}
        if position not in annotationsDict[region][feature]:
            annotationsDict[region][feature][position]={}
        if alternate not in annotationsDict[region][feature][position]:
            annotationsDict[region][feature][position][alternate]={}
        for i,header in enumerate(headers):
            annotationsDict[region][feature][position][alternate][header]=lineArgs[i]
    print 'done..'
    return annotationsDict
    
def parseNucLine(posNucLine,headerArgs):    
    posDict = {} 
    posArgs = posNucLine.split('\t')
    for i,arg in enumerate(headerArgs):
        posDict[arg]=posArgs[i]
    
    return posDict

def getFiltersFromConfigFile(configFile):
    conf = open(configFile.strip(),'r')
    line = conf.readline()
    filters = {}
    while line:
        if line.startswith('#filters_start'):
            line = conf.readline()
            while line and not line.startswith('#filters_end'):
                filterName,filterVal = line.strip().split('\t')
                filters[filterName]=float(filterVal)
                line = conf.readline()
        else:
            line=conf.readline()
    return filters

def getTransitionRate(nucMatrix,metricsFile):
    
    nucs = ['A','C','T','G']
    finalTransitionRate=0
    sumOchanges = 0
    for nuc in nucs:
        if nuc=='A':
            transitionRate = nucMatrix['A']['G']
        elif nuc=='G':
            transitionRate = nucMatrix['G']['A']
        elif nuc=='T':
            transitionRate = nucMatrix['T']['C']
        elif nuc=='C':
            transitionRate = nucMatrix['C']['T']
        sumOfNucVars = sum([nucMatrix[nuc][varNuc] for varNuc in nucs if varNuc!=nuc]) 
        sumOchanges+=sumOfNucVars
        finalTransitionRate+=transitionRate
        if sumOfNucVars==0:
            transRate=0
            transVrate = 0
        else:
            transRate = transitionRate/sumOfNucVars
            transVrate = 1-transRate
        metricsFile.write( '%s transition rate : %.4f , transversion rate : %.4f\n'%(nuc,transRate,transVrate))
    if sumOchanges==0:
        finalTransitionRate=0
    else:  
        finalTransitionRate = finalTransitionRate/sumOchanges
    if finalTransitionRate==1:
        finalTransitonTransversionRate=0
    else:
        finalTransitonTransversionRate = finalTransitionRate/(1-finalTransitionRate)
    metricsFile.write( 'Transition rate : %s\nTransversion Rate : %s\nTransitions/Transversions : %s\n'%(finalTransitionRate,(1-finalTransitionRate),finalTransitonTransversionRate))

def getPassingAlleles(posDict,tests=[]):
#         tests = [re.match('(.+)Vars',header).group(1) for header in posDict.keys() if 'Vars' in header]
#         print tests
        sigAlleles = set(posDict['significantVars'].split(';'))
        for test in tests:
            sigAlleles = sigAlleles.intersection(set(posDict['%sVars'%test].split(';')))
        if sigAlleles!=set(['.']):
            return sigAlleles
        else:
            return []
    
def getMinRateAlleles(posDict,alleles,minRate):
    minRateAlleles=  []
    for allele in alleles:
        rate = float(posDict['%s_rate'%allele].split(';')[0])
        if rate>=minRate:
            minRateAlleles.append(allele)
    return minRateAlleles
    
def getMaxAllele(posDict,significantAlleles):
    maxAllele = ''
    maxAlleleRate = 0
    if len(significantAlleles)>0:
        refAndVars = significantAlleles+[posDict['reference']]
        for allele in refAndVars:
            alleleRate = float(posDict['%s_rate'%allele].split(';')[0])
            if alleleRate>maxAlleleRate:
                maxAlleleRate=alleleRate
                maxAllele = allele
    else:
        return posDict['reference']
    return maxAllele 

def getRefSeqs(referenceFile):
    refSeqs = {}
    refFile = open(referenceFile,'r')
    line = refFile.readline()
    header = None
    headers = []
    seq = ''
    while line:
        if re.match('[>@]', line):
            if seq!='':
                refSeqs[header]=list(seq)
                headers.append(header)
            seq=''
            header = re.sub('[>@]','',line).strip()
            line = refFile.readline()
            continue
        seq+=line.strip().upper()
        line = refFile.readline()
    if seq!='':
        refSeqs[header]=list(seq)
        headers.append(header)

    return refSeqs,headers

def split_len(seq, length):
    return [seq[i:i+length] for i in range(0, len(seq), length)]

def outputConsensusSeqs(consSeqs,prefix,headers):
    outFasta = open('%s.consensus.fa'%prefix,'w')
    seqChunkSize = 50
#     print consSeqs.keys()
    for header in headers:
        seq = ''.join(consSeqs[header])
        seqChunks = split_len(seq,seqChunkSize)
        headerName = '%s_%s'%(prefix,header)
        outFasta.write('>%s\n%s\n'%(headerName,'\n'.join(seqChunks)))

def getInputPrefix(inputFile):
    prefix = re.match('(.+)_nucleotideRate',inputFile).group(1)
    return prefix

def outputNucMatrix(nucs,nucMatrix,prefix):
    outMatrix = open('%s.matrix.csv'%prefix,'w')
    outMatrix.write('\t%s\n'%('\t'.join(nucs)))
    for nuc in nucs:
        rates = '\t'.join(['0' for varNuc in nucs])
        rateSum = sum(nucMatrix[nuc].values())
        if rateSum>0:
            outMatrix.write('%s\t%s\n'%(nuc,'\t'.join(str(nucMatrix[nuc][varNuc]/rateSum) for varNuc in nucs)))

def isSynonymous(region,feature,position,allele,annotationsDict):
    if region not in annotationsDict:
        return 'non-coding'
    if feature not in annotationsDict[region]:
        return 'non-coding'
    if allele=='*':
        return 'nonsynonymous'
    if position in annotationsDict[region][feature]:
        if allele not in annotationsDict[region][feature][position]:
            raise Exception('The given allele (%s) is not found in the tested position : %s'%(allele,annotationsDict[region][feature][position]))
        if 'synonymous' in annotationsDict[region][feature][position][allele]['AA_change']:
            return 'synonymous'
        else:
            return 'nonsynonymous'
    else:
        return 'non-coding'

def updateVariantAlleleRateDict(variantAlleleRate,posDict,significantAlleles,feature,annotationsDict=None): 
    region = posDict['#region']
    if feature not in variantAlleleRate[region]:
        variantAlleleRate[region][feature]={'totalMismatches':0,'synonymous':0,'nonsynonymous':0,'non-coding':0}
    for allele in significantAlleles:
        alleleRate = float(posDict['%s_rate'%allele].split(';')[0])
        variantAlleleRate[region][feature]['totalMismatches']+=alleleRate
        if annotationsDict:
            if feature=='all':
                continue
            syn = isSynonymous(region,feature, posDict['position'], allele, annotationsDict)
            variantAlleleRate[region]['all'][syn]+=alleleRate
            variantAlleleRate[region][feature][syn]+=alleleRate
        
    return variantAlleleRate

def updateMutationFrequencyDict(mutationFrequency,posDict,significantAlleles,feature,annotationsDict=None): 
    region = posDict['#region']
    if feature not in mutationFrequency[region]:
        mutationFrequency[region][feature]={'totalMismatches':0,'synonymous':0,'nonsynonymous':0,'non-coding':0}
    
    mutationFrequency[region][feature]['totalMismatches']+=len(significantAlleles) 
    if annotationsDict:
        if feature=='all':
            return mutationFrequency
        for allele in significantAlleles:
            syn = isSynonymous(region,feature, posDict['position'], allele, annotationsDict)
            mutationFrequency[region][feature][syn]+=1
            mutationFrequency[region]['all'][syn]+=1

    return mutationFrequency
 
    
def outputMutationFrequency(mutationFrequency,metricsFile,variantAlleleOrMutationFrequency,totalBases,perBases=10000): 
    
    metricsFile.write('%s\n'%variantAlleleOrMutationFrequency)
    
    perBasesText = 'per %s bases'%perBases if perBases>1 else 'per base'
    metricsFile.write('%s %s\n'%(variantAlleleOrMutationFrequency,perBasesText))
    
    headers = ['Region',
               'Feature',
               'covered bases',
               'synonymous',
               'non-synonymous',
               'non-coding',
               'total mismatch',
               'synonymous-freq',
               'nonsynonymous-freq',
               'non-coding-frequency',
               'total-freq']
    
    metricsFile.write('%s\n'%'\t'.join(headers))
    for region in mutationFrequency.keys():
        for feature in mutationFrequency[region].keys():
    
            types4freq = ['synonymous','nonsynonymous','non-coding','totalMismatches']
            freqs4metrics = []
            for type in types4freq:
                if totalBases[region][feature]>0:
                    freqs4metrics.append(str(float(mutationFrequency[region][feature][type])*perBases/totalBases[region][feature]))
                else:
                    freqs4metrics.append('Not-covered')
            outMetrics = [region,
                          feature,
                          str(totalBases[region][feature]),
                          '\t'.join([str(mutationFrequency[region][feature][type]) for type in types4freq]),
                          '\t'.join(freqs4metrics)]                        
           
            if feature=='all':
                allOut = outMetrics
                continue
            metricsFile.write('%s\n'%'\t'.join(outMetrics))
        metricsFile.write('%s\n'%'\t'.join(allOut))

def passFilters(posDict,filters,significantAlleles):
    if 'coverage' in filters:
        if int(posDict['coverage'])<filters['coverage']:
            return False
    if 'minRate' in filters:
        for allele in significantAlleles:
            rate = float(posDict['%s_rate'%allele].split(';')[0])
            if rate>=filters['minRate']:
                break
        else:
            return False
    return True

def passCoverageFilter(posDict,filters):
    if 'coverage' in filters:
        if int(posDict['coverage'])<filters['coverage']:
            return False
    return True

def parseIgnorePositions(ignoreArg):
    if not ignoreArg:
        return None
    ignorePos = {}
    ignoreStartEnd = ignoreArg.split(',')
    for regionStartEnd in ignoreStartEnd:
        region,startEnd = regionStartEnd.split(':')
        start,end = [int(pos) for pos in startEnd.split('-')]
        if region not in ignorePos:
            ignorePos[region]=[]
        ignorePos[region].append((start,end))
    return ignorePos
    
def isIgnorePos(posDict,ignorePos):
    # if no ignore positions are given
    if not ignorePos:
        return False    
    region = posDict['#region']
    if region in ignorePos:
        position = int(posDict['position'])
        for (start,end) in ignorePos[region]:
            if position>=start and position<=end:
                return True
    return False

def outputShannonEntropyMetrics(metricsFile,shannonEntropy):
    metricsFile.write('Region\tCoveredBases\tTotalEntropy\tAverageEntropy\n')
    for region in shannonEntropy.keys():
        for feature in shannonEntropy[region]:
            if feature=='all':
                continue
            featureSE = sum(shannonEntropy[region][feature])
            featureAverageSE = featureSE/len(shannonEntropy[region][feature])
            metricsFile.write('%s\t%s\t%s\t%s\n'%(feature,len(shannonEntropy[region][feature]),featureSE,featureAverageSE))
        regionSE = sum(shannonEntropy[region]['all'])
        regionAverageSE = regionSE/len(shannonEntropy[region]['all'])
        metricsFile.write('%s\t%s\t%s\t%s\n'%(region,len(shannonEntropy[region]['all']),regionSE,regionAverageSE))
    
    
def stats():
    parser = optparse.OptionParser()
    parser.add_option('-i','--inputNucleotideVariance')   
    parser.add_option('-r','--reference',default=None)
    parser.add_option('-f','--configFile')
    parser.add_option('-p','--ignorePositions',default=None,help='add comma delimited region:start-end positions in which metrics should not be calculated')
    parser.add_option('-S','--strandBias',default=False,action='store_true',help='add this if metrics should be calculated for quality variants passing strand bias test')
    args = parser.parse_args(sys.argv)[0]
    nucVarFile = open(args.inputNucleotideVariance,'r').xreadlines()
    headerArgs = getHeaderArgs(args.inputNucleotideVariance)
    
    STRAND_BIAS = args.strandBias
    tests = []
    if STRAND_BIAS:
        tests = ['strandBias']
    
    ignorePos = parseIgnorePositions(args.ignorePositions)
    
    annotationsDict = None
    annotationsFileName = re.sub('.csv','.annotation',args.inputNucleotideVariance)
    if os.path.exists(annotationsFileName):
        annotationsDict = parseAnnotationFile(annotationsFileName)
    
    referenceFile = args.reference    
    prefix = getInputPrefix(args.inputNucleotideVariance)
    
    metricsFile = open('%s.metrics.txt'%prefix,'w')
    filters = getFiltersFromConfigFile(args.configFile)
    
    MIN_RATE = 0
    if 'minRate' in filters:
        MIN_RATE = float(filters['minRate'])
    
    # get the sequences from the reference file
    consSeqs,headers = getRefSeqs(referenceFile)
    nucs = ['A','C','T','G']
    
    nucMatrix = dict([(nuc,dict([(nuc,0) for nuc in nucs])) for nuc in nucs])
    changedNucs = {}
    
    numOchangedNucs = 0
    ignoredPositions = 0
    totalBases = {}
    mutationFrequency={}
    variantAlleleRate={}
    shannonEntropy = {}
    for pos in nucVarFile:
        if pos.startswith('#'):
            continue
        posDict = parseNucLine(pos, headerArgs)
        
        if isIgnorePos(posDict, ignorePos):
            ignoredPositions+=1
            continue
        
        region = posDict['#region']
        
        if region not in totalBases:
            totalBases[region] = {'all':0}
        
        if region not in mutationFrequency:
            mutationFrequency[region] = {'all':{'totalCoverage':0,'totalMismatches':0,'synonymous':0,'nonsynonymous':0,'non-coding':0}}
            variantAlleleRate[region] = {'all':{'totalCoverage':0,'totalMismatches':0,'synonymous':0,'nonsynonymous':0,'non-coding':0}}
        
        position = int(posDict['position'])
        refBase = consSeqs[region][position-1]
        features = [posDict['#region']]
        if 'features' in headerArgs:
            features = posDict['features'].split(';')

        # validate that the base in the nucleotide rate file is the same as in
        # the reference file
        if posDict['reference']!=refBase:
            raise Exception('base mismatch in position %s consensus : %s, nucVar : %s'%(posDict['position'],refBase,posDict['reference']))
        
        if not passCoverageFilter(posDict, filters):
            continue
        
        # Shannon Entropy 
        if region not in shannonEntropy:
            shannonEntropy[region]={'all':[]}
        for feature in features:
            if feature not in shannonEntropy[region]:
                shannonEntropy[region][feature]=[]
            shannonEntropy[region][feature].append(float(posDict['ShannonEntropy']))
        shannonEntropy[region]['all'].append(float(posDict['ShannonEntropy']))
        
        # get the significantly changed alleles in the given positiongetSignificantAlleles
        significantAlleles = getPassingAlleles(posDict,tests)
        
        # get the significant alleles that pass the minimal rate filter
        finalAlleles = getMinRateAlleles(posDict, significantAlleles, MIN_RATE)
        
        for feature in features:
            if feature not in totalBases[region]:
                totalBases[region][feature]=0
            totalBases[region][feature]+=1
        totalBases[region]['all']+=1
        
        for feature in features:
            variantAlleleRate = updateVariantAlleleRateDict(variantAlleleRate, posDict, finalAlleles, feature,annotationsDict)
            mutationFrequency = updateMutationFrequencyDict(mutationFrequency, posDict, finalAlleles, feature,annotationsDict)
        
        variantAlleleRate = updateVariantAlleleRateDict(variantAlleleRate, posDict, finalAlleles, 'all',annotationsDict)
        mutationFrequency = updateMutationFrequencyDict(mutationFrequency, posDict, finalAlleles, 'all',annotationsDict)

        
        # if the significant allele is also the maximal, change the consensus
        maxAllele = getMaxAllele(posDict, finalAlleles) 
        if maxAllele!=refBase:
            consSeqs[posDict['#region']][int(posDict['position'])-1]=maxAllele
            if region not in changedNucs:
                changedNucs[region]=[]
            changedNucs[region].append('%s:%s>%s'%(posDict['position'],posDict['reference'],maxAllele))
            numOchangedNucs+=1
        for nuc in nucs:
            refAndVars = [posDict['reference']]+finalAlleles
            if nuc in refAndVars:
                nucRate = float(posDict['%s_rate'%nuc].split(';')[0])
                if posDict['reference']=='N':
                    continue
                nucMatrix[posDict['reference']][nuc]+=nucRate
            
    outputNucMatrix(nucs, nucMatrix,prefix)
    metricsFile.write('\nTransitions / Transversions\n===================\n')
    getTransitionRate(nucMatrix,metricsFile)
    
    metricsFile.write('\nVariation Metrics\n================\n')
    outputMutationFrequency(mutationFrequency, metricsFile,'Region Heterogeneity',totalBases)
    outputMutationFrequency(variantAlleleRate, metricsFile,'Region Variation Rate',totalBases)

    metricsFile.write('\nShannon Entropy\n================\n')
    outputShannonEntropyMetrics(metricsFile, shannonEntropy)

    outputConsensusSeqs(consSeqs, prefix,headers)
    
    changedNucsLine = ''
    for regionName in changedNucs.keys():
        changedNucsLine = '%s>%s\n%s\n'%(changedNucsLine,regionName,','.join(changedNucs[regionName]))
        

    metricsFile.write('\nConsensus sequences had %s base changes:\n%s'%(numOchangedNucs,changedNucsLine))
    print 'Ignored %s positions across the genome'%ignoredPositions
if __name__=='__main__':
    stats()
