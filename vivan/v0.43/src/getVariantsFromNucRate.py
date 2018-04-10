#!/usr/bin/env python

#===============================================================================
# Get variants from nucRate
# takes in a nucRate file and returns a simple list of all the variants detected
#===============================================================================

VERSION = 0.1

print 'VIVAN: getVariantsFromNucRate v%s'%VERSION

import optparse
import re
import sys
import time
from Bio.Data.CodonTable import TranslationError
from Bio import Seq, SeqIO
from Bio import pairwise2
import difflib


def parseNucRateFile(nucRateFileName):
    print 'VIVAN: parsing %s'%nucRateFileName
    nucRates = []
    with open(nucRateFileName,'r') as nucRateFile:
        for line in nucRateFile:
            posDict = {}
            lineArgs = line.split('\t')
            if line.startswith('#'):
                headers = lineArgs
                continue
            for i,header in enumerate(headers):
                posDict[header]=lineArgs[i]
            nucRates.append(posDict)
    return nucRates,headers

def getVarTests(headers):
    tests = []
    for header in headers:
        if 'Vars' in header:
            tests.append(re.match('(.+)Vars',header).group(1))
    return tests

def getPassingAlleles(posDict,tests):
    passingAlleles = set()
    for test in tests:
        if posDict['%sVars'%test]!='.':
            passingAlleles = passingAlleles.union(set(posDict['%sVars'%test].split(';')))
    return passingAlleles

def printVariants(nucRateFile,headers,nucRates,tests):
    prefix = re.match('(.+)_nucleotideRate.csv',nucRateFile).group(1)
    outFileName = '%s.variants'%prefix
    outFile = open(outFileName,'w')
    print 'VIVAN: Printing variants to %s'%outFileName
    headerArgs = ['#Region',
               'Position',
               'Allele',
               'AlleleRate',
               'CI',
               'FailedTests']
    if 'features' in headers:
        headerArgs.append('features')
    outFile.write('%s\n'%'\t'.join(headerArgs))
    for posDict in nucRates:
        passingAlleles = getPassingAlleles(posDict, tests)
        failedTests = [test for test in tests if posDict['%sVars'%test]=='.']
        if len(failedTests)==0:
            failedTests = ['.']
        for allele in passingAlleles:
            outArgs = [posDict['#region'],
                       posDict['position'],
                       allele,
                       '\t'.join(posDict['%s_rate'%allele].split(';')),
                       ','.join(failedTests)]
            if 'features' in headers:
                outArgs.append(posDict['features'])
            outFile.write('%s\n'%'\t'.join(outArgs))
            
def stats():

    startTime = time.time()
    #===========================================================================
    # Input
    #===========================================================================
    parser = optparse.OptionParser()
    parser.add_option('-i','--nucRateFile')
    
    args = parser.parse_args(sys.argv)[0]
    nucRates,headers = parseNucRateFile(args.nucRateFile)
    
    tests = getVarTests(headers)
    printVariants(args.nucRateFile, headers, nucRates, tests)
if __name__=='__main__':
    stats()
        
