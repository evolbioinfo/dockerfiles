#!/usr/bin/env python

#===============================================================================
# Parse Variant File
# takes in a variant file (from pileup2nucleotideRate) and a 
# bed file of the different features (CDS start and stop) and a possible FASTA file?
#===============================================================================
VALIDATE = True
VERSION = 0.64

print 'VIVAN: ParseVariantFile v%s'%VERSION

import optparse
import re
import sys
import time
from Bio.Data.CodonTable import TranslationError
from Bio import Seq, SeqIO
from Bio import pairwise2
import difflib

WARNINGS = []

def getFileLen(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

class BaseVar:
    def __init__(self,baseLine,headerArgs):
        self.baseArgs = baseLine.strip().split('\t')
        baseDict = {}
        for i,header in enumerate(headerArgs):
            baseDict[header]=self.baseArgs[i]
        self.baseDict = baseDict
        self.headerArgs = headerArgs
        self.ref = baseDict['reference'].upper()
        self.pos = baseDict['position']
        self.region = baseDict['#region']
        self.coverage = int(baseDict['coverage'])
        self.features = baseDict['features'].split(';')
        self.bases = {'A':0,'C':0,'T':0,'G':0}
        self.ratios = {'A':0,'C':0,'T':0,'G':0}
        
        for allele in self.bases.keys():
            self.bases[allele]=int(baseDict[allele])
        
        for base in self.bases.keys():
#            self.ratios[base]=float(self.bases[base])/self.coverage
            self.ratios[base]=float(baseDict['%s_rate'%base].split(';')[0]) 
               
    def getMajorAF(self):
        maf = 0
        for base in self.bases.keys():
            baseDepth = self.bases[base]
            baseFreq = float(baseDepth)/self.coverage
            if baseFreq>maf:
                maf=baseFreq
        return maf
    
    def getNonRefMajorAllele(self):
        maxAlleleCount = 0
        maxAllele = ''
        for allele in self.bases.keys():
            if allele==self.ref:
                continue
            else:
                if maxAlleleCount<self.bases[allele]:
                    maxAllele=allele
                    maxAlleleCount=self.bases[allele]
        return (maxAllele,maxAlleleCount)
    
    def getPassingAlleles(self,tests=[]):
#         tests = [re.match('(.+)Vars',header).group(1) for header in self.headerArgs if 'Vars' in header]
        sigAlleles = set(self.baseDict['significantVars'].split(';'))
        for test in tests:
            sigAlleles = sigAlleles.intersection(set(self.baseDict['%sVars'%test].split(';')))
        if sigAlleles!=set(['.']):
            return sigAlleles
        else:
            return []
        
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
        self.cds = None
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
            if VALIDATE and cdsStartSites:
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


def getHeaderArgs(header):
    headerArgs = [arg.strip() for arg in header.split('\t')]
    return headerArgs
  
def parseVariantFile(varFile,tests):
    print 'VIVAN: Parsing %s'%varFile
    varDict = {}
    fileLen =getFileLen(varFile)
    print 'VIVAN: There are %s lines in %s'%(fileLen,varFile)
    varFile = open(varFile,'r').xreadlines()
    regions = []
    bases = 0
    for base in varFile:
        if base.startswith('#'):
            headerArgs = getHeaderArgs(base)
            continue
        baseVar = BaseVar(base,headerArgs)
        bases+=1
        if baseVar.region not in varDict:
            varDict[baseVar.region]={}
            regions.append(baseVar.region)
        if divmod(bases,1000)[1]==0:
            print '\rVIVAN: Parsed %s/%s bases'%(bases,fileLen),
            sys.stdout.flush()
        if len(baseVar.getPassingAlleles(tests))==0:
            continue
        varDict[baseVar.region][baseVar.pos]=baseVar
       
    return varDict,regions
        
def parseFeatureBed(bedFile,regionSeqs):
    print 'VIVAN: Parsing %s'%bedFile
    bedFile = open(bedFile,'r').xreadlines()
    features = {}
    for line in bedFile:
        if line.strip():
            feature = Feature(line)
            regionSeq = regionSeqs[feature.region]
            feature.cds = feature.getCDS(regionSeq)
            try:
                translated = Seq.translate(feature.cds,cds=True)
            except TranslationError,e:
                translated = Seq.translate(feature.cds)
                WARNINGS.append('Translation error in feature : %s\nCDS : %s\nProtein : %s\n%s\n'%(feature.featureLine,feature.cds,translated,e))
            features[feature.name]=feature
    print 'VIVAN: Done, collected %s features'%(len(features))
    return features

def parseReference(referenceFileName):
    regionSeqs = {}
    ref = SeqIO.parse(open(referenceFileName), 'fasta')
    for record in ref:
        regionSeqs[record.id]=record.seq
    return regionSeqs

def getFeatureCDS(feature,referenceFileName):
    ref = SeqIO.parse(open(refernceFileName), 'fasta')
    for record in ref:
        if record.id==feature.region:
            cds = feature.getCDS(record.seq)
            break
    else:
        raise Exception('could not find a region %s in %s'%(feature.region,referenceFileName))
    return cds

def alnSequences(seq1,seq2):
    if len(seq1)!=len(seq2):
        raise Exception('Alignment exception : the two sequences have different lengths: \n%s\n%s\n'%(seq1,seq2))
    seq1 = list(seq1.upper())
    seq2 = list(seq2.upper())
    changes = []
    for i in xrange(len(seq1)):
        if seq1[i]!=seq2[i]:
            changes.append([i+1,seq1[i].upper(),seq2[i].upper()])
    return changes

def getAAnum(cdsNum):
    div,mod = divmod(cdsNum, 3)
    if mod!=0:
        div = div+1
    return div

def getFiltersFromConfigFile(configFile):
    conf = open(configFile,'r')
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

def isPassFilter(base,filters):
    for filter in filters.keys():
        PASS=True
        if filter=='coverage':
            PASS = base.coverage>=filters[filter]
        if not PASS:
            return False
    else:
        return True
        
def isSignificantAllele(allele,base,tests):
    return allele in base.getPassingAlleles(tests)

def translateSeq(cds):
    senseOrAnti = 'sense'
    finalCDS = cds
    try:
        translated = Seq.translate(cds,cds=True)
#         finalCDS = cds
    except TranslationError,e:
        try:
            reverseCDS = Seq.reverse_complement(cds)
            translated = Seq.translate(reverseCDS,cds=True)
            finalCDS = reverseCDS
            senseOrAnti = 'anti'
        except TranslationError,e:
            print 'Translation failed in %s'%cds
    return Seq.translate(finalCDS),senseOrAnti
    
def stats():

    startTime = time.time()
    #===========================================================================
    # Input
    #===========================================================================
    parser = optparse.OptionParser()
    parser.add_option('-f','--configFile')
    parser.add_option('-v','--varFile')
    parser.add_option('-b','--bedFile')
    parser.add_option('-r','--reference',default=None)
    parser.add_option('-c','--validateCDS',default=True,action='store_false')
    parser.add_option('-S','--strandBias',default=False,action='store_true')
    
    args = parser.parse_args(sys.argv)[0]
    varFileName = args.varFile
    bedFileName = args.bedFile
    referenceFileName = args.reference
    
    VALIDATE_CDS = args.validateCDS
    STRAND_BIAS = args.strandBias
    
    tests = []
    if STRAND_BIAS:
        tests = ['strandBias']
    
    filters = getFiltersFromConfigFile(args.configFile)
    # minimal rate below which variants are excluded. use this to reduce noise
    MIN_RATE = 0
    if 'minRate' in filters:
        MIN_RATE = float(filters['minRate'])
    
#     regionSeqs = getSequenceFromVarFile(varFileName)
    regionSeqs = parseReference(referenceFileName)
    bedFeatures = parseFeatureBed(bedFileName,regionSeqs)
    varDict,regions = parseVariantFile(varFileName,tests)
    
    varPrefix = re.match('[^.]+',varFileName).group(0)
    
    outFile = open('%s.annotation'%varPrefix,'w')
    warningsFile = open('%s.warnings'%varPrefix,'w')
    
    outFile.write('#region\tfeature\tposition\tcoverage\treference\talternate\talt_frequency\tCDS_position\tAA_position\tAA_change\n')
    for region in regions:
        print 'VIVAN: Analyzing bases in %s'%region
        def sortByPosition(position):
            return int(position)
        
        regionSeq = regionSeqs[region]
        numOpos=0
        positions = sorted(varDict[region].keys(),key=sortByPosition)
        for position in positions:
            numOpos+=1
            baseVar = varDict[region][position]
            if VALIDATE:
                if baseVar.ref!=regionSeq[int(position)-1].upper():
                    raise Exception('Position : %s, the reference allele in the variant file (%s) is not the same as the one in the reference sequence (%s).'%(position,baseVar.ref,regionSeq[int(pos)-1]))            
            def sortByRatio(base):
                return baseVar.ratios[base]
            if isPassFilter(baseVar, filters):
                sortedVars = sorted(baseVar.ratios.keys(),key=sortByRatio,reverse=True)
                
                for base in sortedVars: 
                    # check if base rate is above minRate threshold
                    if baseVar.ratios[base]<MIN_RATE:
                        continue
                    # check if allele variance found significant
                    if not isSignificantAllele(base, baseVar,tests):
                        continue
                    # check if the base is not the reference base
                    if base.upper()==baseVar.ref.upper():
                        continue
                    
                    if baseVar.features!=['.']:
                        for featureName in baseVar.features:
                            feature = bedFeatures[featureName]
                            cds = feature.cds
#                             translated = Seq.translate(cds)
                            translated,senseOrAnti = translateSeq(cds)
                            
                            modified_seq = list(regionSeq)
                            modified_seq[int(position)-1]=base
                            
                            modified_cds = feature.getCDS(''.join(modified_seq))
                            cds_len = len(modified_cds)
                            AA_change = '.\t.'
                            # check if the sequence change is found affects the coding
                            # sequence (for example if the change is found outside of the 
                            # spliced CDS)
                            if str(cds)!=str(modified_cds):
                                cds_changes = alnSequences(cds, modified_cds)
                                for cds_pos,ref,alt in cds_changes:
                                    if senseOrAnti=='anti':
                                        modified_translated = Seq.translate(Seq.reverse_complement(modified_cds))
                                    else:
                                        modified_translated = Seq.translate(modified_cds)
                                    if str(translated)!=str(modified_translated):
                                        trans_changes =  alnSequences(translated, modified_translated)
                                        for trans_pos,ref,alt in trans_changes:
                                            AA_change = '%s\t%s>%s'%(trans_pos,ref,alt)
                                    else:
                                        try:
                                            AA_pos = getAAnum(cds_pos)
                                            if AA_pos>(cds_len/3):
                                                warningsFile.write('Coding sequence error on %s : the number of nucleotides (%s) cant be divided by 3.\nthere is a qualifying mutation at cds position %s which cannot be translated\n'%(feature.name,cds_len,cds_pos))
                                                continue
                                            AA = modified_translated[AA_pos-1]
                                            AA_change = '%s\tsynonymous>%s'%(AA_pos,AA) 
                                        except IndexError:
                                            raise IndexError('Index error : %s\ncds pos : %s\ncds len : %s\nAA pos: %s\nAA length %s\n'%(feature.name,cds_pos,len(cds),AA_pos,len(modified_translated)))     
                            else:
                                cds_pos = '.'
                            outArgs = [baseVar.region,
                                       feature.name,
                                       baseVar.pos,
                                       baseVar.coverage,
                                       baseVar.ref,
                                       base,
                                       baseVar.ratios[base],
                                       cds_pos,
                                       AA_change]        
                            outFile.write('%s\n'%'\t'.join([str(arg) for arg in outArgs]))                     
                    else:
                        outArgs = [baseVar.region,
                                   'non-coding',
                                   baseVar.pos,
                                   baseVar.coverage,
                                   baseVar.ref,
                                   base,
                                   baseVar.ratios[base],
                                   'NA',
                                   'NA',
                                   'NA']
                        outFile.write('%s\n'%'\t'.join([str(arg) for arg in outArgs]))   
       
                     #       print outLine
            if divmod(numOpos,500)[1]==0:
                print '\rVIVAN: Finished %s/%s positions'%(numOpos,len(positions)),
                sys.stdout.flush()
    for warning in WARNINGS:
      #  print warning
        warningsFile.write('%s\n'%warning)
        
if __name__=='__main__':
    stats()
        
