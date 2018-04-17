#!/usr/bin/env python

import optparse
import re
import sys
import time
from os import chdir
import itertools
from os import popen
from os import path
from scipy import stats

VERSION = 0.61
#===============================================================================
# Compare nucleotide variance files
#===============================================================================
print 'VIVAN: CompareNucVarFiles v%s'%VERSION

def defineGroupDict(groupDict,group):
    if group not in groupDict:
        groupDict[group]={'samples':{},'reference':None,'features':None}
    return groupDict

def parseConfigFile(confFile,delim='\t'):
    confFile = open(confFile,'r')
    confLine = confFile.readline()
    groupDict = {}
    filters = {}
    while confLine:
        if confLine.startswith('#samples_start'):
            confLine = confFile.readline()
            while not(confLine.startswith('#samples_end')):
                sampleArgs = confLine.strip().split(delim)
                sample = sampleArgs[0]
                group = sampleArgs[1]
                fileNames = sampleArgs[2:]
                defineGroupDict(groupDict, group)
                groupDict[group]['samples'][sample]=fileNames
                confLine = confFile.readline()
        elif confLine.startswith('#referencesNfeatures_start'):
            confLine = confFile.readline()    
            while not(confLine.startswith('#referencesNfeatures_end')):    
                groupNreferenceNfeature = confLine.strip().split(delim)
                group = groupNreferenceNfeature[0]
                defineGroupDict(groupDict, group)
                groupDict[group]['reference']= groupNreferenceNfeature[1]
                # if feature file is also given
                if len(groupNreferenceNfeature)>2:
                    groupDict[group]['features']=groupNreferenceNfeature[2]
                confLine = confFile.readline()
        elif confLine.startswith('#filters_start'):
            confLine = confFile.readline()    
            while not(confLine.startswith('#filters_end')):    
                filter,threshold = confLine.strip().split(delim)
                filters[filter]=float(threshold)
                confLine = confFile.readline()
        else:
            confLine = confFile.readline()
            
    return groupDict,filters

def getFileNames(fileList):
    # take in the file names from the file list
    fileLines = open(fileList).xreadlines()
    ls = popen('ls').read()
    fileNames = []
    for file in fileLines:
        fileName,groupName = file.split('\t')  
        fileName = re.search('(%s[^\s]*nucleotideRate.csv)'%fileName.strip(),ls).group(1)
        fileNames.append((fileName,groupName.strip()))
    return fileNames

def getColumnNum(fileName,columnTitle):
    file = open(fileName,'r')
    firstLine = file.readline()
    for i,title in enumerate(firstLine.split('\t')):
        if columnTitle==title.strip():
            #print '%s : %s : %s'%(fileName,columnTitle,i)
            return i
    else:
        raise Exception('could not find %s in %s header :\n%s'%(columnTitle,fileName,firstLine))

def isRMSDvariant(RMSD,RMSDthreshold):
    return RMSD>=RMSDthreshold

def isSufficientCoverage(baseCov,baseCovThreshold):
    return baseCov>=baseCovThreshold

def isLowRefAllele(baseRefRate,baseRefThreshold):
    return baseRefRate<=baseRefThreshold

def isAnyVarAboveMinRate(varRates,minRate):
    minRate = float(minRate)
#    varsAboveMinRate = []
    for varRate in varRates:
        if varRate>=minRate:
            return True
    return False

def getSignificantAlleles(significance):
    sigAllels = significance.split(';')
    if sigAllels!=['.']:
        return sigAllels
    else:
        return []

def getPassingAlleles(baseDict,tests=[]):
#         tests = [re.match('(.+)Vars',header).group(1) for header in baseDict.keys() if 'Vars' in header]
        sigAlleles = set(baseDict['significantVars'].split(';'))
        for test in tests:
            sigAlleles = sigAlleles.intersection(set(baseDict['%sVars'%test].split(';')))
        if sigAlleles!=set(['.']):
            return sigAlleles
        else:
            return []
    

def getNumOFilesInGroup(fileNames):
    groupDict = {}
    for fileName,groupName in fileNames:
        if groupName not in groupDict:
            groupDict[groupName]=[fileName]
        else:
            groupDict[groupName].append(fileName)
    
    for groupName in groupDict.keys():
        print '%s has %s files'%(groupName,len(groupDict[groupName]))
    
    return groupDict

def getAlleleFreqFromBaseLine(baseDict,sortedAlleles,significantAlleles,AAchanges,region,MIN_RATE):

    if AAchanges:
        alleleFreqLine = {}
        for feature in AAchanges[region].keys():
#            print '%s AAChange : '%feature,AAchanges[region][feature]
            alleleFreqLine[feature]= ';'.join(['%s:%s:%.5f'%(allele,AAchanges[region][feature][allele],float(baseDict['%s_rate'%allele].split(';')[0])) for allele in sortedAlleles if (allele in significantAlleles and float(baseDict['%s_rate'%allele].split(';')[0])>=MIN_RATE)])
    else:
        alleleFreqLine=';'.join(['%s:%.5f'%(allele,float(baseDict['%s_rate'%allele].split(';')[0])) for allele in sortedAlleles if allele in significantAlleles])
    return alleleFreqLine

def compGroups(groupNames,groupDict,posDict, INCLUDE_FEATURES):
    print 'now testing : %s'%' vs '.join(groupNames)
    combName = 'vs'.join(groupNames)
    if len(groupNames)==len(groupDict.keys()):
        combName= 'allGroupsComparison'
    popen('mkdir %s'%combName)
    chdir('./%s'%combName)
    print 'curr dir : %s'%popen('pwd').read()
    #popen('cd %s'%combName)
    for groupName in groupNames:
        
        numOcommon = 0
        numOfilter = 0
        numOcommonFilter = 0
        
        commonPositions = open('%s.common.csv'%groupName,'w')
        filterPositions  = open('%s.unique.csv'%groupName,'w')
        commonFilterPositions = open('%s.common.unique.csv'%groupName,'w')
   
        outFiles = [commonPositions,filterPositions,commonFilterPositions]
        for outFile in outFiles:
            if INCLUDE_FEATURES:
                outFile.write('region\tfeature\tposition\treference\tcoverage\t%s\n'%('\t'.join(groupDict[groupName]['samples'].keys())))
            else:
                outFile.write('region\tposition\treference\tcoverage\t%s\n'%('\t'.join(groupDict[groupName]['samples'].keys())))    
   
        numOgroupFiles = len(groupDict[groupName]['samples'].keys())
        for region in posDict[groupName].keys():
            for position in sorted(posDict[groupName][region].keys()):
                
                numOpositionFiles = len([fileName for fileName in posDict[groupName][region][position]['files'].keys() if posDict[groupName][region][position]['files'][fileName]!='.'])
                # check if the position is common in all the files in the group
                common =  numOpositionFiles==numOgroupFiles
                
                # check if the position is not found in any of the other group files
                for otherGroup in groupNames:
                    if otherGroup!=groupName:
                        if region in posDict[otherGroup]:
                            if position in posDict[otherGroup][region]:
                                filter=False
                                break
                else:
                    filter=True
                
                # create the position line to write in the appropriate files
                positionLines = []
                if common or filter:
                    features = posDict[groupName][region][position]['features']
                    avrBaseCov = float(posDict[groupName][region][position]['baseCov'])/numOpositionFiles
                    if features and features!=['.']:
                        for feature in features:
                            positionLine = ['%s\t%s\t%s\t%s\t%.2f'%(region,feature,position,posDict[groupName][region][position]['reference'],avrBaseCov)]
                            for fileName in groupDict[groupName]['samples'].keys():
                                if feature in posDict[groupName][region][position]['files'][fileName]:
#                                    print feature
#                                    print posDict[groupName][region][position]['files'][fileName]
                                    positionLine.append(posDict[groupName][region][position]['files'][fileName][feature])
                                else:
                                    positionLine.append(posDict[groupName][region][position]['files'][fileName])
                            positionLine = '\t'.join(positionLine)
                            positionLines.append(positionLine)
                        
                    else:
                        positionLine = ['%s\t.\t%s\t%s\t%.2f'%(region,position,posDict[groupName][region][position]['reference'],avrBaseCov)]
                        for fileName in groupDict[groupName]['samples'].keys():
                            positionLine.append(posDict[groupName][region][position]['files'][fileName])
                        positionLine = '\t'.join(positionLine)
                        positionLines.append(positionLine)
#                print positionLines
                positionLines = '\n'.join(positionLines)
                
                if common:
                    commonPositions.write('%s\n'%positionLines)
                    numOcommon+=1
                if filter:
                    filterPositions.write('%s\n'%positionLines)
                    numOfilter+=1
                if common and filter:
                    commonFilterPositions.write('%s\n'%positionLines)
                    numOcommonFilter+=1
        
        for outFile in outFiles:
            outFile.close()
        print '%s had:\n%s common variant positions\n%s unique variant positions\n%s common unique variant positions'%(groupName,numOcommon,numOfilter,numOcommonFilter)
    #popen('cd ..')    
    chdir('..')     
    popen('mv %s ./comparisons'%combName)      

def parseAnnotationFile(annotationFileName):
    #feature    position    coverage    reference    alternate    alt_frequency    CDS_position    AA_position    AA_change
    annotationFile = open(annotationFileName,'r')
    featuresDict= {}
    for line in annotationFile:
        if line.startswith('#'):
            continue
        lineArgs = line.strip().split('\t')
        region = lineArgs[0]
        if region not in featuresDict:
            featuresDict[region]={}
        featureName = lineArgs[1]
        if featureName not in featuresDict[region]:
            featuresDict[region][featureName]={}
        alternate = lineArgs[5]
        AApos = lineArgs[8]
        AAchange = lineArgs[9]
        position = int(lineArgs[2])
        if position not in featuresDict[region][featureName]:
            featuresDict[region][featureName][position] = {}
            featuresDict[region][featureName][position]['*']=''
        featuresDict[region][featureName][position][alternate]='%s:%s'%(AApos,AAchange)
        
    return featuresDict

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
        
def getGroupConsenus(posDict,groupName,reference,refHeaders):
    print 'VIVAN: Producing %s consenus sequence'%groupName
    consChanges = []
    newCons = {}
    seqChunkSize = 50
    for region in posDict[groupName].keys():
        newRef = list(reference[region])
        for position in posDict[groupName][region].keys():
            #Validation
            refBase = newRef[int(position)-1]
            if posDict[groupName][region][position]['reference']!=refBase:
                raise Exception('VIVAN: compareNucVarFiles: The reference sequence did not have the same base as the one specified in the position dict (%s instead of %s)'%(posDict[groupName][region][position]['reference'],refBase))
            consensusVars = posDict[groupName][region][position]['consensus']
            if len(consensusVars)>0:
                newConsBase = list(consensusVars)[0]
                consChanges.append('%s:%s:%s>%s'%(region,position,posDict[groupName][region][position]['reference'],newConsBase))
                newRef[int(position)-1]=newConsBase
        newCons[region]=''.join(newRef)
    outFile = open('%s.consensus.fa'%(groupName),'w')
    for header in refHeaders:
        if header not in newCons:
            continue
        seq = newCons[header]
        seqChunks = split_len(seq,seqChunkSize)
        outFile.write('>%s\n%s\n'%(header,'\n'.join(seqChunks)))
    print 'VIVAN: Consensus Changes : %s'%','.join(consChanges)

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
        refSeqs[header]=seq
        headers.append(header)

    return refSeqs,headers

                
def compare():

    startTime = time.time()
    #===========================================================================
    # Input
    #===========================================================================
    parser = optparse.OptionParser()
    parser.add_option('-f','--confFile')
    parser.add_option('-N','--disregardN',default=False,action='store_true')
    parser.add_option('-v','--verbose',default=False,action='store_true')
    parser.add_option('-F','--includeFeatures',default=False,action='store_true')
    parser.add_option('-S','--strandBias',default=False,action='store_true',help='add this if only variants passing strand bias should be compared')
    
    args = parser.parse_args(sys.argv)[0]
    
    VERBOSE = args.verbose
    DISREGARD_N = args.disregardN
    confFile = args.confFile
    INCLUDE_FEATURES = args.includeFeatures
    STRAND_BIAS = args.strandBias
    
    tests = []
    if STRAND_BIAS:
        tests=['strandBias']
    
    groupDict,filters = parseConfigFile(confFile)
    
    posDict = {}
    
   # filters = args.filters.split(',')
    for filter in filters.keys():
        if filter=='coverage':
            print 'Coverage filter : only positions with coverage higher than %sx will be included'%filters[filter]
        if filter=='RMSD':
            print 'RMSD filter : only positions with RMSD higher than %s will be included'%filters[filter]
        if filter=='baseRef':
            print 'Reference allele rate filter : only positions where reference allele is lower than %s will be included'%filters[filter]
        if filter=='pval':
            print 'Significane filter : only positions where there was at least one significant (p<%s) allele variance'%filters[filter]
        if filter=='minRate':
            print 'Minimal Rate filter : only variants with a rate higher or equal to %s will be included'%filters[filter]
    for groupName in groupDict.keys():
        for sample in groupDict[groupName]['samples'].keys():
            fileName = '%s_nucleotideRate.csv'%sample
            annotationFile = '%s_nucleotideRate.annotation'%sample
            annotationDict = {}
            if path.exists(annotationFile):
                annotationDict = parseAnnotationFile(annotationFile)
            
            numOpositions = 0
            numOpassPos = 0
            numOfailPos = 0
            
            if groupName not in posDict:
                posDict[groupName]={}
            
            alleles = ['A','C','T','G','N','*']
            if DISREGARD_N:
                alleles = ['A','C','T','G','*']
            
            file = open(fileName,'r').xreadlines()
            for baseLine in file:
                baseDict = {}
                baseArgs = baseLine.strip().split('\t')
                if baseLine.startswith('#'):
                    headers = baseArgs
                    continue
                for i,header in enumerate(headers):
                    baseDict[header]=baseArgs[i]
                
                numOpositions+=1
                
                region = baseDict['#region']
                position = int(baseDict['position'])
                baseCov = int(baseDict['coverage'])
                baseRef = baseDict['reference']
                
                if INCLUDE_FEATURES:
                    features = baseDict['features'].strip().split(',')
                else:
                    features = None
                
                baseRefRate = float(baseDict['%s_rate'%baseRef].split(';')[0])
                
                RMSD = float(baseDict['RMSD'])
                PASS = False
                
                significantAlleles = getPassingAlleles(baseDict,tests)
                baseConsensus = [var for var in significantAlleles if float(baseDict['%s_rate'%var].split(';')[0])>0.5]
                
                varRates = [float(baseDict['%s_rate'%var].split(';')[0]) for var in significantAlleles]
                MIN_RATE = 0
                for filter2check in filters.keys():
                    if filter2check=='coverage':
                        PASS = isSufficientCoverage(baseCov, filters[filter2check])
                    if filter2check=='RMSD':
                        PASS = isRMSDvariant(RMSD,filters[filter2check])
                    if filter2check=='baseRef':
                        PASS = isLowRefAllele(baseRefRate,filters[filter2check])
                    if filter2check=='pval':
                        PASS = len(significantAlleles)>0
                    if filter2check=='minRate':
                        MIN_RATE = float(filters[filter2check])
                        PASS = isAnyVarAboveMinRate(varRates, filters[filter2check])    
                    
                    if not PASS:
                        numOfailPos+=1
                        break
                
                if PASS:
                    AAchanges = None
                    numOpassPos+=1
                    if features:
                        AAchanges={}
                        for feature in features:
                            if region in annotationDict: 
                                if feature in annotationDict[region]:
                                    if position in annotationDict[region][feature]:
                                        if region not in AAchanges:
                                            AAchanges[region]={}
                                        AAchanges[region][feature]=annotationDict[region][feature][position]
                                    else: 
                                        # the position passes the filters,
                                        # but is not found in the annotation
                                        # file. should only occur when
                                        # a significant deletion occurs
                                        if '*' in significantAlleles:
                                            AAchanges=None
                                        else:
                                            raise Exception('Comparison script : %s : could not find position %s in feature %s in the annotations file although the position passed the filters'%(sample,position,feature))

                    if region not in posDict[groupName]:
                        posDict[groupName][region]={}
                    
                    if position not in posDict[groupName][region]:
                        posDict[groupName][region][position]={}
                        posDict[groupName][region][position]['baseCov']=baseCov
                        posDict[groupName][region][position]['reference']=baseRef
                        posDict[groupName][region][position]['files']=dict([(groupSamp,'.') for groupSamp in groupDict[groupName]['samples'].keys()])
                        posDict[groupName][region][position]['features']=features
                        posDict[groupName][region][position]['consensus']=set(baseConsensus)
                        
                    posDict[groupName][region][position]['files'][sample]=getAlleleFreqFromBaseLine(baseDict,alleles,significantAlleles,AAchanges,region,MIN_RATE)
                    posDict[groupName][region][position]['baseCov']+=baseCov
                    posDict[groupName][region][position]['consensus']=posDict[groupName][region][position]['consensus'].intersection(set(baseConsensus))
                    
            print '%s had %s positions: %s passed %s filters and %s did not'%(fileName,numOpositions,numOpassPos,' and '.join(filters.keys()),numOfailPos)
    
    for groupName in groupDict.keys():
        reference,headers= getRefSeqs(groupDict[groupName]['reference'])
        getGroupConsenus(posDict, groupName, reference,headers)
    
    popen('mkdir comparisons')
    numOgroups = len(groupDict.keys())
    # this can result in thousands of comparisons, 
    # therefore im limiting it to 2 for now
    
    for i in range(2,3):
        groupCombs = itertools.combinations(groupDict.keys(),i)
        # compare all groups combinations
        for groupComb in groupCombs:
            compGroups(groupComb,groupDict,posDict, INCLUDE_FEATURES)
    if len(groupDict.keys())>2:    
        #compare all  
        compGroups(groupDict.keys(), groupDict, posDict, INCLUDE_FEATURES)
    
                
if __name__=='__main__':
    compare()
    
