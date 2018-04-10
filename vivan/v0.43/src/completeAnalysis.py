#!/usr/bin/env python
#===============================================================================
# ViVan complete analysis script.

#===============================================================================
# PLEASE UPDATE THE PATHES BELOW BEFORE RUNNING THIS SCRIPT
#===============================================================================

# Usage : python completeAnalysis.py [arguments]
# -f : the configuration file (see example and fill accordingly)
# -c : add this argument if adapter clipping and quality trimming should are required prior to alignment
# -a : add this argument if the input is sequence data and requires both alignment and pileup
# -p : Skip alignment and only perform pileup
# -P : add this if pileup 2 nucleotide rate has already been performed
# -N : add this argument if sequence reads with unknown (N) alleles should be discarded
# -V : add this if variant annotation has already been performed
# -M : Use this flag if your only interested in producing variant metrics for all the samples in the configuration file
# -S : add this if you want only variants that pass strand bias to be analyzed
# -I : internal do not use.
# -Q : internal do not use.
# -q : internal do not use.
# -l : internal do not use.
# -W : internal do not use.
#===============================================================================

# Dependancies:
# SAMtools is expected to be installed and in the path
# BWA (v0.7.8) is expected to be installed and found inside the ViVan scripts directory
# ea-utils (https://code.google.com/p/ea-utils/) folder must be installed inside the ViVan scripts directory
# Python v>2.6 is expected to be installed with these modules:
# Biopython, scipy and numpy

# The scripts directory, please change this to match the ViVan scripts directory
SCRIPTS_DIR = '/usr/local/ViVan/'

# The bwa path, should be inside the ViVan scripts directory
BWA_PATH = '/usr/local/ViVan/bwa/'

#===============================================================================
import optparse
import sys
import re
import time
import sys
import os
from os import system
from os import popen
import subprocess
import random
import datetime

VERSION = '0.4'
SUBVERSION= '3'


summaryReport = []

def getFileLen(fName):
    f = open(fName,'r')
    for i, l in enumerate(f):
        pass
    f.close()
    return i + 1

def verifyJobsComplete(jobNames):
    for jobScriptName in jobNames:
        errFile = open('%s.ERR'%jobScriptName,'r').read()
        isErr = re.search('Traceback.+',re.sub('\n','xxxx',errFile))
        if isErr:
            raise Exception('ViVan: Exception: %s:\n%s'%(jobScriptName,re.sub('xxxx','\n',isErr.group(0))))
    
def waitForJobs(jobNums):
    # wait for jobs, will only be used when QSUB flag is true
    SLEEP_TIME = 5
    timeWait = 0
    while jobNums!=[]:
        qstat = popen('qstat').read()
        for job in jobNums:
            if job not in qstat:
                jobNums.remove(job)
                print '\njob %s done'%job
                if jobNums==[]:
                    return
        time.sleep(SLEEP_TIME)
        timeWait+=SLEEP_TIME
        print ('\rwaited %s seconds'%timeWait),
        sys.stdout.flush()

def defineGroupDict(groupDict,group):
    # a function that defines a group dictionary
    if group not in groupDict:
        groupDict[group]={'samples':{},'reference':None,'features':None}
    return groupDict

def validateGroupDict(groupDict):
    # a validation that should be ran for each group after config file parsing
    for group in groupDict.keys():
        # validate that each group has at least one sample assigned to it
        if len(groupDict[group]['samples'])==0:
            raise Exception('Config File Exception : %s has no samples'%group)

def testQual(fileName):
    qualOut = popen('python %s/getQualFormat.py -f %s'%(SCRIPTS_DIR,fileName)).read()
    if 'FORMAT:sanger' in qualOut:
        return 'sanger'
    elif 'FORMAT:illumina' in qualOut:
        return 'illumina'
    else:
        raise Exception('could not determine the quality format in %s : %s'%(fileName,qualOut))

def fixQualFormats(groupDict,QSUB,ILLUMINA_FORMAT,VERBOSE,LIGHT,QUEUE):
    
    format2fix = 'illumina'
    outFormat = 'sanger'
    if ILLUMINA_FORMAT:
        format2fix = 'sanger'
        outFormat = 'illumina'
    
    # Align a given fastq file
    startTime = time.time()
    print 'Fixing File Quality Formats..'
    jobNums = []
    jobNames = []
    for group in groupDict.keys():
        for sample in groupDict[group]['samples'].keys():
            fileNames = groupDict[group]['samples'][sample]
            for fileName in fileNames:
                qual = testQual(fileName)
                if qual==format2fix:
                    groupDict[group]['samples'][sample]=[]
                    fixCommand = 'python %s/convertFileFormat.py -i %s -p %s -f fastq-%s -t fastq-%s'%(SCRIPTS_DIR,fileName,sample,format2fix,outFormat)
                    commmands = [fixCommand]
                    # if QSUB is true, submit a job with all the commands
                    if QSUB:
                        jobName = 'FixQual'
                        jobNum = qsubJob(sample, commmands,jobName,LIGHT,QUEUE)
                        jobNames.append('%s%s'%(sample,jobName))
                        jobNums.append(jobNum)
                    # if QSUB is false, run the commands one by one
                    else:
                        output = popen(fixCommand).read()
                        if VERBOSE:
                            print output
                    # fix the file name after conversion
                    fixFileName = '%s.fastq-%s'%(fileName,outFormat)
                    groupDict[group]['samples'][sample].append(fixFileName)
                
    if QSUB:
        waitForJobs(jobNums)
        verifyJobsComplete(jobNames)
    print 'Quality Fix done..(%s)\n'%(time.time()-startTime)
    

def parseConfigFile(confFile,delim='\t'):
    # a function for parsing the config file
    confFile = open(confFile,'r')
    confLine = confFile.readline()
    groupDict = {}
    filters = {}
    while confLine:
        # get the samples , groups and file names 
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
        # get the reference and feature files for each group
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
        # get the filters that should be used in the group comparison section
        elif confLine.startswith('#filters_start'):
            confLine = confFile.readline()    
            while not(confLine.startswith('#filters_end')):    
                filter,threshold = confLine.strip().split(delim)
                filters[filter]=float(threshold)
                confLine = confFile.readline()
        else:
            confLine = confFile.readline()
    # validate that the created group dict is ok
    validateGroupDict(groupDict)
    return groupDict,filters

def printGroups(groupDict):
    # print out the data parsed from the config file
    for group in sorted(groupDict.keys()):
        summaryReport.append('%s\n=================\n'%group)
        summaryReport.append( 'Reference\t:\t%s'%groupDict[group]['reference'])
        if groupDict[group]['features']:
            summaryReport.append( 'Features\t:\t%s'%groupDict[group]['features'])
        for sample in groupDict[group]['samples']:
            summaryReport.append( '%s\t:\t%s'%(sample,groupDict[group]['samples'][sample]))
        summaryReport.append('')
    
    print '\n'.join(summaryReport)
    
def qsubJob(sample,commands,jobSuffix=None,LIGHT=False,QUEUE=None):
    # a function that submits a job and returns the job number
    shscript="#!/bin/csh\n#$ -N %s%s\n#$ -S /bin/csh\n#$ -cwd\n#$ -o $JOB_NAME.OUT\n#$ -e $JOB_NAME.ERR\n#$ -l eran\n\nmodule load  gcc/gcc461\nmodule load  python/python-2.7.2\nmodule load  python/python-2.7.6\n"%(sample,jobSuffix)
    if QUEUE:
        shscript="#!/bin/csh\n#$ -N %s%s\n#$ -S /bin/csh\n#$ -cwd\n#$ -o $JOB_NAME.OUT\n#$ -e $JOB_NAME.ERR\n#$ -l %s\n\n\nmodule load  gcc/gcc461\nmodule load  python/python-2.7.2\nmodule load  python/python-2.7.6\n"%(sample,jobSuffix,QUEUE)
    if LIGHT:
        shscript="#!/bin/csh\n#$ -N %s%s\n#$ -S /bin/csh\n#$ -cwd\n#$ -o $JOB_NAME.OUT\n#$ -e $JOB_NAME.ERR\n#$ -l light\n#$ -ckpt BLCR\n\nmodule load blcr/blcr-0.8.2\nmodule load  gcc/gcc461\nmodule load  python/python-2.7.2\nmodule load  python/python-2.7.6\n"%(sample,jobSuffix)
        for i,command in enumerate(commands):
            commands[i]='cr_run %s'%command
    comScript = '%s%s\n'%(shscript,'\n'.join(commands))
    
    jobFileName = '%s%s.sh'%(sample,jobSuffix)
    
    jobFile = open(jobFileName,'w')
    jobFile.write(comScript)
    jobFile.close()
    out = popen('qsub %s'%(jobFileName)).read()
    jobNum = re.search('\d+', out).group(0)
    
    print '%s : submitted %s job %s'%(sample,jobSuffix,jobNum)
    return jobNum


def getAvrReadLength(fileName):
    avrReadLen = popen("head -100 %s"%fileName+"| awk '{if(NR%4==2) sum+=length($1)} END {print sum/25}'").read()
    return float(avrReadLen)

def Align(groupDict,NUM_O_GAPS,QSUB,ILLUMINA_FORMAT,VERBOSE,LIGHT,QUEUE):
    # Align a given fastq file
    startTime = time.time()
    print 'ViVAN: Alignment starting..'
    summaryReport.append('Alignment starting..')
    jobNums = []
    jobNames = []
    for group in groupDict.keys():
        reference = groupDict[group]['reference']
        for sample in groupDict[group]['samples'].keys():
            fileNames = groupDict[group]['samples'][sample]
            avrReadLen = getAvrReadLength(fileNames[0])
            if avrReadLen<70:
                # Single End
                if len(fileNames)==1:
                    alnCommand = '%s aln -o %s -t 10 %s %s | %s samse -r "@RG\tID:%s\tPL:ILLUMINA\tSM:viVan" %s - %s > %s.sam'%(BWA_PATH,NUM_O_GAPS,reference,fileNames[0],BWA_PATH,sample,reference,fileNames[0],sample)
                # Paired End
                else:
                    aln1Command = '%s aln -o %s -t 10 %s %s > %s_1.sai'%(BWA_PATH,NUM_O_GAPS,reference,fileNames[0],sample)
                    aln2Command = '%s aln -o %s -t 10 %s %s > %s_2.sai'%(BWA_PATH,NUM_O_GAPS,reference,fileNames[1],sample)
                    sampeCommand = '%s sampe -r "@RG\tID:%s\tPL:ILLUMINA\tSM:viVan" %s %s_1.sai %s_2.sai %s > %s.sam'%(BWA_PATH,sample,reference,sample,sample,' '.join(fileNames),sample)
                    alnCommand = '\n'.join([aln1Command,aln2Command,sampeCommand])
            else:
                print 'ViVAN: Read length longer than 70 (%s) will use bwa mem'%avrReadLen
                alnCommand = '%s mem -M -t 5 -R "@RG\tID:%s\tPL:ILLUMINA\tSM:viVan" %s %s > %s.sam'%(BWA_PATH,sample,reference,' '.join(fileNames),sample)

            bamNsortCommand = 'samtools view -bS %s.sam | samtools sort - %s_sorted'%(sample,sample)
            indexCommand = 'samtools index %s_sorted.bam'%sample
            
            commands = [alnCommand,bamNsortCommand,indexCommand]
            
            # if QSUB is true, submit a job with all the commands
            if QSUB:
                jobName = 'Aln'
                jobNum = qsubJob(sample, commands,jobName,LIGHT,QUEUE)
                jobNames.append('%s%s'%(sample,jobName))
                jobNums.append(jobNum)
            # if QSUB is false, run the commands one by one
            else:
                for command in commands:
                    output = popen(command).read()
                    if VERBOSE:
                        print output
    if QSUB:
        waitForJobs(jobNums)
        verifyJobsComplete(jobNames)
    print 'ViVAN: Alignment done..(%s)\n'%(time.time()-startTime)
    summaryReport.append('Alignment done..(%s)\n'%(time.time()-startTime))

def Pile(groupDict,QSUB,ILLUMINA_FORMAT,VERBOSE,LIGHT,QUEUE):
    # Perform pileup on the aligned file
    startTime = time.time()
    print 'ViVAN: Pileup starting..'
    summaryReport.append('Pileup starting..')
    jobNums = []
    jobNames = []
    for group in groupDict.keys():
        reference = groupDict[group]['reference']
        for sample in groupDict[group]['samples'].keys():

            pileCommand = 'samtools mpileup -B -Q 30 -f %s -d 10000000 %s_sorted.bam > %s.pileup' %(reference,sample,sample)
            
            commands = [pileCommand]
            
            # if QSUB is true, submit a job with all the commands
            if QSUB:
                jobName = 'Pile'
                jobNum = qsubJob(sample, commands,jobName,LIGHT,QUEUE)
                jobNames.append('%s%s'%(sample,jobName))
                jobNums.append(jobNum)
            # if QSUB is false, run the commands one by one
            else:
                for command in commands:
                    output = popen(command).read()
                    if VERBOSE:
                        print output
    if QSUB:
        waitForJobs(jobNums)
        verifyJobsComplete(jobNames)
    print 'ViVAN: Pileup done..(%s)\n'%(time.time()-startTime)
    summaryReport.append('Pileup done..(%s)\n'%(time.time()-startTime))


def removeNreads(groupDict,QSUB,VERBOSE,LIGHT,QUEUE):
    startTime = time.time()
    
    REMOVE_N_SCRIPT = '%s/removeNreads.py'%SCRIPTS_DIR
    
    print 'ViVAN: Removal of reads with N starting..'
    summaryReport.append('Clipping and Trimming starting..')
    jobNums = []
    jobNames = []
    for group in groupDict.keys():
        for sample in groupDict[group]['samples'].keys():
            fileNames = groupDict[group]['samples'][sample]
            if len(fileNames)!=1:
                print 'ViVAN: Reads containing unknown bases (N) will not be removed in sample %s since it is a paired end sample'%sample
                continue
            fileName = fileNames[0]
            removeNs = 'python %s -i %s -o %s'%(REMOVE_N_SCRIPT,fileName,sample)
            commands = [removeNs]
        
            if QSUB:
                jobName = 'removeN'
                jobNum = qsubJob(sample, commands,jobName,LIGHT,QUEUE)
                jobNames.append('%s%s'%(sample,jobName))
                jobNums.append(jobNum)
            else:
                for command in commands:
                    output = popen(command).read()
                    summaryReport.append(output)
                    if VERBOSE:
                        print output
    if QSUB:
        waitForJobs(jobNums)
        verifyJobsComplete(jobNames)
    print 'ViVAN: Removal of reads with N done..(%s)\n'%(time.time()-startTime)
    summaryReport.append('Removal of reads with N done..(%s)\n'%(time.time()-startTime))
    return groupDict

def clipAndTrim(groupDict,QSUB,VERBOSE,LIGHT,QUEUE):
    # Clip and Trim the given fastq files using a default adapters file 
    # while keeping only reads >15 nt and bases with quality>30
    startTime = time.time()
    print 'ViVAN: Clipping and Trimming starting..'
    summaryReport.append('Clipping and Trimming starting..')
    jobNums = []
    jobNames = []
    CLIP_AND_TRIM = True
    for group in groupDict.keys():
        reference = groupDict[group]['reference']
        for sample in groupDict[group]['samples'].keys():
            fileNames = groupDict[group]['samples'][sample]
            if len(fileNames)==1:
                outFiles = ['%s.trimmed.fq'%sample]
                clipNtrim = '%s/ea-utils/fastq-mcf -C 200000 -l 16 -q 30 -o %s %s/adapters.fa %s'%(SCRIPTS_DIR,' '.join(outFiles),SCRIPTS_DIR,fileNames[0])
            else:
                outFiles = ['%s.r1.trimmed.fq'%sample,'%s.r2.trimmed.fq'%sample]
                clipNtrim = '%s/ea-utils/fastq-mcf -C 200000 -l 16 -q 30 -o %s %s/adapters.fa %s'%(SCRIPTS_DIR,' -o '.join(outFiles),SCRIPTS_DIR,' '.join(fileNames))
            commands = [clipNtrim]
            
            # change the name of the file
            groupDict[group]['samples'][sample]=outFiles
            
            if QSUB:
                jobName = 'ClipNTrim'
                jobNum = qsubJob(sample, commands,jobName,LIGHT,QUEUE)
                jobNames.append('%s%s'%(sample,jobName))
                jobNums.append(jobNum)
            else:
                for command in commands:
                    output = popen(command).read()
                    print output
                    summaryReport.append(output)
                    if VERBOSE:
                        print output
    if QSUB:
        waitForJobs(jobNums)
        verifyJobsComplete(jobNames)
    
    print 'ViVAN: Clipping and Trimming done..(%s)\n'%(time.time()-startTime)
    summaryReport.append('Clipping and Trimming done..(%s)\n'%(time.time()-startTime))
    return groupDict

def pileup2nucleotideRate(groupDict,filters,QSUB,DISREGARD_N,VERBOSE,LIGHT,QUEUE,SIZE_LIMIT):
    # Pileup 2 Nucleotide rate - takes in a pileup file and creates a nucleotide
    # variance file with the allele frequency for each allele and the significance
    # of the change
    startTime = time.time()
    print 'ViVAN: Pileup 2 Nucleotide Rate starting..'
    summaryReport.append('Pileup 2 Nucleotide Rate starting..')
    
    
    PILE2NUC_SCRIPT = '%s/pileup2nucleotideRate.py'%SCRIPTS_DIR
    
    jobNums = []
    jobNames = []
    for group in groupDict.keys():
        for sample in groupDict[group]['samples'].keys():
            pileName = '%s.pileup'%sample
            nucRatePrefix = re.match('[^\.]+', sample).group(0)
            nucRateFileName = '%s_nucleotideRate.csv'%nucRatePrefix
            print 'checking %s'%nucRateFileName
            if os.path.exists('./%s'%nucRateFileName):
                print '%s: %s already exists, will not run pile2nucleotideRate analysis on it again'%(sample,nucRateFileName)
                continue
                
            pile2nucCommand = 'python %s -i %s'%(PILE2NUC_SCRIPT,pileName)
            
            if 'pval' in filters:
                pile2nucCommand = '%s -p %s'%(pile2nucCommand,filters['pval'])
            
            if DISREGARD_N:
                pile2nucCommand = '%s -N'%pile2nucCommand

            if groupDict[group]['features']:
                pile2nucCommand = '%s -f %s'%(pile2nucCommand,groupDict[group]['features'])
            
            commands = [pile2nucCommand]
    
            if QSUB:
                jobName = 'Pile2Nuc'
                jobNum = qsubJob(sample, commands,jobName,LIGHT,QUEUE)
                jobNames.append('%s%s'%(sample,jobName))
                jobNums.append(jobNum)
            else:
                for command in commands:
                    print 'running %s'%command
                    output = popen(command).read()
                    summaryReport.append(output)
                    if VERBOSE:
                        print output

    if QSUB:
        waitForJobs(jobNums)
        verifyJobsComplete(jobNames)

    print 'ViVAN: Pileup 2 Nucleotide Rate done..(%s)\n'%(time.time()-startTime)
    summaryReport.append('Pileup 2 Nucleotide Rate done..(%s)\n'%(time.time()-startTime))
    
def varAnnotation(groupDict,confFile,STRAND_BIAS,QSUB,VERBOSE,LIGHT,QUEUE):
    # Variation annotation - takes in a pile2nuc file and returns all the variant
    # positions and the predicted AA change according to the given reading
    # frame in the bed feature file 
    startTime = time.time()
    print 'Variant annotation starting..'
    summaryReport.append('Variant annotation starting..')
    
    VARANNO_SCRIPT = '%s/parseVariantFile.py'%SCRIPTS_DIR
    
    jobNums = []
    jobNames = []
    for group in groupDict.keys():
        reference = groupDict[group]['reference']
        featureFile = groupDict[group]['features']
        if not featureFile:
            print '%s has no features bed file, will not annotate it..'%group
            continue
        for sample in groupDict[group]['samples'].keys():
            fileName = groupDict[group]['samples'][sample]
            varAnnoCommand = 'python %s -v %s_nucleotideRate.csv -f %s -b %s -r %s -c'%(VARANNO_SCRIPT,sample,confFile,featureFile,reference)
            if STRAND_BIAS:
                varAnnoCommand+=' -S'
            commands = [varAnnoCommand]
    
            if QSUB:
                jobName = 'VarAnno'
                jobNum = qsubJob(sample, commands,jobName,LIGHT,QUEUE)
                jobNames.append('%s%s'%(sample,jobName))
                jobNums.append(jobNum)
            else:
                for command in commands:
                    print 'VIVAN:SHELLCOM: %s'%command
                    output = popen(command).read()
                    summaryReport.append(output)
                    if VERBOSE:
                        print output
    if QSUB and len(jobNums)>0:
        waitForJobs(jobNums)
        verifyJobsComplete(jobNames)
    print 'Variation annotation done..(%s)\n'%(time.time()-startTime)
    summaryReport.append('Variation annotation done..(%s)\n'%(time.time()-startTime))

def produceVarMetrics(groupDict,confFile, STRAND_BIAS,QSUB,VERBOSE,LIGHT,QUEUE):
    # Variance metrics production, takes in a nucleotide rate file and produces
    # 1 - a consensus sequence fasta
    # 2 - a nucleotide matrix file
    # 3 - some metrics on transition/transversion rates
    startTime = time.time()
    print 'ViVAN: Producing Variance metrics..'
    summaryReport.append('Producing Variance metrics..')
    
    
    VARMETRIC_SCRIPT = '%s/produceVarMetrics.py'%SCRIPTS_DIR
    
    jobNums = []
    jobNames = []
    for group in groupDict.keys():
        reference = groupDict[group]['reference']
        featureFile = groupDict[group]['features']
        if not featureFile:
            print '%s has no features bed file, will not annotate it..'%group
            continue
        for sample in groupDict[group]['samples'].keys():
            varMetCom = 'python %s -i %s_nucleotideRate.csv -r %s -f %s'%(VARMETRIC_SCRIPT,sample,reference,confFile)
            if STRAND_BIAS:
                varMetCom+=' -S'
            
            commands = [varMetCom]
    
            if QSUB:
                jobName = 'VarMetric-%s'%sample
                jobNum = qsubJob(sample, commands,jobName,LIGHT,QUEUE)
                jobNames.append('%s%s'%(sample,jobName))
                jobNums.append(jobNum)
                
            else:
                for command in commands:
                    output = popen(command).read()
                    summaryReport.append(output)
                    if VERBOSE:
                        print output
    if QSUB and len(jobNums)>0:
        waitForJobs(jobNums)
        verifyJobsComplete(jobNames)
    print 'ViVAN: Variance metrics production done..(%s)\n'%(time.time()-startTime)
    summaryReport.append('Variance metrics production done..(%s)\n'%(time.time()-startTime))

def compareGroups(confFile,STRAND_BIAS,QSUB,DISREGARD_N,INCLUDE_FEATURES,VERBOSE,LIGHT,QUEUE):
    # Compare groups - run a comparison between the different groups. 
    # first step - check for each sample in each group what are the positions 
    # that pass the given filters (coverage / pval / AF)
    # second step - for each group, output three files : 
    #    - the positions that pass the filter in all the samples in the group
    #    - the unique positions that pass the filter in at least one sample in the group and not 
    startTime = time.time()
    print 'ViVAN: Group comparison starting..'
    summaryReport.append('Group comparison starting..')
    jobNames = []
    COMP_GROUPS_SCRIPT = '%s/compareNucVarFiles.py'%SCRIPTS_DIR
    if DISREGARD_N:
        COMP_GROUPS_SCRIPT ='%s -N'%COMP_GROUPS_SCRIPT
    if INCLUDE_FEATURES:
        COMP_GROUPS_SCRIPT = '%s -F'%COMP_GROUPS_SCRIPT
    if STRAND_BIAS:
        COMP_GROUPS_SCRIPT = '%s -S'%COMP_GROUPS_SCRIPT
    jobNums = []
    
    commands = ['python %s -f %s'%(COMP_GROUPS_SCRIPT,confFile)]
    
    if QSUB:
        jobName = 'GroupCompare'
        jobNum = qsubJob('final', commands,jobName,LIGHT,QUEUE)
        jobNames.append('%s%s'%('final',jobName))
        jobNums.append(jobNum)
    else:
        for command in commands:
            output = popen(command).read()
            summaryReport.append(output)
            if VERBOSE:
                print output
    if QSUB and len(jobNums)>0:
        waitForJobs(jobNums)
        verifyJobsComplete(jobNames)
    summaryReport.append('Group Comparison done..(%s)\n'%(time.time()-startTime))
    print 'ViVAN: Group Comparison done..(%s)\n'%(time.time()-startTime)

def cleanup(groupDict,WEBSERVER):
    print 'VIVAN: Cleanup..'

    scriptsDir = './scripts'
    if not os.path.exists(scriptsDir):
        os.mkdir(scriptsDir)
    popen('mv *.ER* *.OU* *.sh %s 2>/dev/null'%scriptsDir)

    fastqDir = './fastq'
    alignmentDir = './alignAndPile'
    if not os.path.exists(alignmentDir):
        os.mkdir(alignmentDir)
    if not os.path.exists(fastqDir):
        os.mkdir(fastqDir)
    popen('mv *.fq *.filter.fq *.fastq %s 2>/dev/null'%(fastqDir))
    popen('mv *.sam *.bam* *.pileup %s 2>/dev/null'%(alignmentDir))
    
        
    for group in groupDict.keys():
        groupDir = './%s'%group
        if not os.path.exists(groupDir):
            os.mkdir(groupDir)
        popen('mv %s.consensus.fa %s'%(group,groupDir))
        for sample in groupDict[group]['samples'].keys():
            sampleDir = '%s/%s'%(groupDir,sample)
            if not os.path.exists(sampleDir):
                os.mkdir(sampleDir)       
            popen('mv %s_* %s.* %s 2>/dev/null'%(sample,sample,sampleDir))

    if WEBSERVER:
        # delete the alignment and fastq files
        popen('rm -r %s %s'%(fastqDir,alignmentDir))
        refDir= './REF'
        # move all reference files to the reference dir
        if not os.path.exists(refDir):
            os.mkdir(refDir)
        popen('mv *.fa* %s'%refDir)
        # create the results folder and copy results to it
        resFolder = 'ViVan-Results'
        if not os.path.exists(resFolder):
            os.mkdir(resFolder)
        for group in groupDict.keys():
            groupDir = './%s'%group
            popen('mv %s %s'%(groupDir,resFolder))
            popen('mv comparisons %s'%resFolder)
        
    print 'VIVAN: Done'
    
def finalize():
    print 'VIVAN: Finalizing.'
    finished = open('finished.txt','w')

    timestamp = datetime.datetime.now().strftime('ViVan-%y%m%d-%H%M%S')
    finished.write(timestamp)
    finished.close()
    
def verifyInput(groupDict):
    print 'ViVAN: Verifying input..'
    
    def getRefRegions(referenceName):
        regions = set()
        referenceFile = open(referenceName,'r')
        for line in referenceFile:
            if line.startswith('>'):
                regions.add(re.match('>(.+)',line.strip()).group(1))
        referenceFile.close()
        return regions
    
    def getFeaturesRegions(featuresName):
        fregions = set()
        featuresFile = open(featuresName,'r')
        for line in featuresFile:
            lineArgs = line.split('\t')
            fregions.add(lineArgs[0])
        return fregions
    
    fixedRefsAndFeatures = set()
    for group in groupDict:
        #=======================================================================
        # Verify sample and group names do not start with numbers
        #=======================================================================
        if re.match('\d',group):
            raise Exception('VIVAN: Input Exception: Config file : Group names cannot start with numebrs (BAD NAME: %s)'%group)
        for sample in groupDict[group]['samples']:
            if re.match('\d',sample): 
                raise Exception('VIVAN: Input Exception: Config file : Sample names cannot start with numebrs (BAD NAME: %s)'%sample)
        
        #=======================================================================
        # File existance validation
        #=======================================================================
        groupReference = groupDict[group]['reference'].strip()
        if not os.path.isfile(groupReference):
            raise Exception('Config file : the given reference file does not exist :\n%s'%groupReference)
        if groupDict[group]['features']:
            groupFeatures = groupDict[group]['features'].strip()
            if not os.path.isfile(groupFeatures):
                raise Exception('VIVAN: Input Exception: Config file : the given features file does not exist :\n%s'%groupFeatures)
        else:
            print 'WARNING : %s : no features file given..'%group
        
        #=======================================================================
        # Fix reference and features header names 
        #=======================================================================
        refName = groupDict[group]['reference'].strip()
        if refName not in fixedRefsAndFeatures:
            print 'VIVAN: Fixing %s headers'%refName
            tmpRefName = '%s.tmp'%refName
            with open(tmpRefName,'w') as tmpRef:
                with open(refName,'r') as ref:
                    for line in ref:
                        tmpRef.write('%s\n'%line.strip().replace(' ','_'))
            os.popen('mv -f %s %s'%(tmpRefName,refName))
            fixedRefsAndFeatures.add(refName)
        if groupDict[group]['features']:
            featName = groupDict[group]['features']
            if featName not in fixedRefsAndFeatures:
                print 'VIVAN: Fixing %s headers'%featName
                tmpFeatName = '%s.tmp'%featName
                with open(tmpFeatName,'w') as tmpFeat:
                    with open(featName,'r') as feat:
                        for line in feat:
                            lineArgs = line.split('\t')
                            lineArgs[0]=lineArgs[0].strip().replace(' ','_')
                            tmpFeat.write('\t'.join(lineArgs))
                os.popen('mv -f %s %s'%(tmpFeatName,featName))
                fixedRefsAndFeatures.add(featName)
        
        #=======================================================================
        # Validate reference and features regions match
        #=======================================================================
        if groupDict[group]['features']:
            print '%s : Validating the reference and features file regions match..'%group
            refRegions = getRefRegions(groupDict[group]['reference'])
            featRegions = getFeaturesRegions(groupDict[group]['features'])
            regionsUnion = refRegions.union(featRegions)
            if len(regionsUnion)>len(refRegions):
                extraRegions = regionsUnion.difference(refRegions)
                raise Exception('VIVAN: Input Exception: Reference (%s) and features (%s) mismatch :\nthere is a region (%s) in the features not found in the reference'%(groupDict[group]['reference'],groupDict[group]['features'],' ; '.join(extraRegions)))
        #=======================================================================
        # Reference indexing validation
        #=======================================================================
        referencePathArgs = groupReference.split('/')
        referencePath = '/'.join(referencePathArgs[:-1])
        referenceFile = referencePathArgs[-1]
        refereceDirFiles = popen('ls %s'%referencePath).read().strip().split('\n')
        suffixes = ['.amb','.ann','.bwt','.fai','.pac','.rbwt','.rpac','.rsa','.sa']
        for suffix in suffixes:
            if not os.path.exists('%s%s'%(groupDict[group]['reference'],suffix)):
                print 'The Database is not properly indexed (missing %s%s)..\nwill run db indexing..'%(groupDict[group]['reference'],suffix)
                popen('%s index %s'%(BWA_PATH,groupDict[group]['reference']))
                popen('samtools faidx %s'%groupDict[group]['reference'])
                break
    print 'done..'     
    
def analyze():
    parser = optparse.OptionParser(description='ViVan v%s-%s ; Run this file in order to perform the complete ViVan analysis.\nRequirements :\nBWA and SAMTools should both be installed in the PATH\n ; Python version >= 2.7'%(VERSION,SUBVERSION))
    parser.add_option('-f','--configFile',help='a configuration file containing all the analysis parameters')
    parser.add_option('-c','--clipNtrim',default=False,action='store_true',help='add this if adapter clipping and quality trimming should are required prior to alignment')
    parser.add_option('-a','--alignNpile',default=False,action='store_true',help='add this if the input is sequence data and requires both alignment and pileup')
    parser.add_option('-p','--onlyPileup',default=False,action='store_true',help='Skip alignment and only perform pileup')
    parser.add_option('-P','--skipPile2Nuc',default=False,action='store_true',help='add this if pileup 2 nucleotide rate has already been performed')
    parser.add_option('-V','--skipVarAnnotation',default=False,action='store_true',help='add this if variant annotation has already been performed')
    parser.add_option('-q','--useQsub',default=False,action='store_true',help='if a q manager is installed, add this argument to submit jobs through it')
    parser.add_option('-v','--verbose',default=False,action='store_true',help='print full output')
    parser.add_option('-I','--illuminaFormat',default=False,action='store_true',help='add this argument if the data is in illumina format')
    parser.add_option('-N','--removeNreads',default=False,action='store_true',help='add this argument if sequence reads with unknown (N) alleles should be discarded')
    parser.add_option('-l','--lightQue',default=False,action='store_true',help='run on light queue (private)')
    parser.add_option('-M','--onlyMetrics',default=False,action='store_true',help='Use this flag if your only interested in producing variant metrics for all the samples in the configuration file')
    parser.add_option('-Q','--queue',default=None,help='use this if you want to designate a specific queue (comp0.eran.q / comp1.eran.q). also accepts comma seperated queues')
    parser.add_option('-g','--numOgaps',default=1,help='max number of gaps when comparing to reference')
    parser.add_option('-C','--cleanupOnly',default=False,action='store_true',help='add this argument if only cleanup is necessary')
    parser.add_option('-S','--strandBias',default=False,action='store_true',help='add this if you want only variants that pass strand bias to be analyzed')
    parser.add_option('-L','--sizeLimit',default=10**7,help='change this value to limit the size of the input pileup, larger pileups will be split into a 100 jobs')
    parser.add_option('-W','--webserver',default=None,help='set this argument if the script is ran through the webserver')
    args = parser.parse_args(sys.argv)[0]
    
    startTime = time.time()
    
    CLIP_AND_TRIM = args.clipNtrim
    ALIGN_AND_PILE = args.alignNpile
    ONLY_PILEUP = args.onlyPileup
    ILLUMINA_FORMAT = args.illuminaFormat
    REMOVE_N_READS = args.removeNreads
    SKIP_PILE2NUC = args.skipPile2Nuc
    SKIP_VARANNO = args.skipVarAnnotation
    ONLY_METRICS = args.onlyMetrics
    QUEUE = args.queue
    NUM_O_GAPS = int(args.numOgaps)
    CLEANUP_ONLY = args.cleanupOnly
    STRAND_BIAS = args.strandBias
    SIZE_LIMIT = int(args.sizeLimit)
    WEBSERVER = args.webserver
    if WEBSERVER:
        if not os.path.exists(WEBSERVER):
            raise Exception('VIVAN:Exception : Could not find directory %s'%WEBSERVER)
        print 'VIVAN: Opening %s'%WEBSERVER
        os.chdir(WEBSERVER)
    LIGHT = args.lightQue
    QSUB = args.useQsub
    VERBOSE = args.verbose
    confFile = args.configFile
    groupDict,filters = parseConfigFile(confFile)
    
    # Verify Input
    verifyInput(groupDict)
    
    # only include features if all the groups have a feature file in them
    for group in groupDict.keys():
        if groupDict[group]['features']:
            continue
        else:
            INCLUDE_FEATURES=False
            break
    else:
        INCLUDE_FEATURES=True
    
    printGroups(groupDict)
    if CLEANUP_ONLY:
        cleanup(groupDict,WEBSERVER)
        return
    # if alignment is required, check the qualities of the input files
    if ALIGN_AND_PILE:
        fixQualFormats(groupDict,QSUB,ILLUMINA_FORMAT,VERBOSE,LIGHT,QUEUE)
    # only produce the analysis metrics (consensus, matrix, metrics)
    if ONLY_METRICS:
        produceVarMetrics(groupDict,confFile, STRAND_BIAS,QSUB, VERBOSE,LIGHT,QUEUE)
        return
    # clip and trim the sequence files
    if CLIP_AND_TRIM:
        groupDict = clipAndTrim(groupDict, QSUB, VERBOSE,LIGHT,QUEUE)
        #N reade will be removed only if clip and trim are also running
        if REMOVE_N_READS:
            removeNreads(groupDict, QSUB,VERBOSE,LIGHT,QUEUE)
            # now change all the names of the sample files

    if REMOVE_N_READS:    
        for group in groupDict.keys():
            for sample in groupDict[group]['samples'].keys():
                # only for single end
                if len(groupDict[group]['samples'][sample])==1:
                    groupDict[group]['samples'][sample]=['%s.filter.fq'%sample]
    
    # Either perform alignment and pileup or only pileup                
    if ALIGN_AND_PILE:
        Align(groupDict,NUM_O_GAPS, QSUB, ILLUMINA_FORMAT, VERBOSE, LIGHT, QUEUE)
        Pile(groupDict, QSUB, ILLUMINA_FORMAT, VERBOSE, LIGHT, QUEUE)
    elif ONLY_PILEUP:
        Pile(groupDict, QSUB, ILLUMINA_FORMAT, VERBOSE, LIGHT, QUEUE)
    
    # run the pileup 2 nucleotide rate script (pileup2nucleotideRate.py)
    if not SKIP_PILE2NUC:
        pileup2nucleotideRate(groupDict, filters,QSUB,REMOVE_N_READS,VERBOSE,LIGHT,QUEUE,SIZE_LIMIT)
    # run the variant annotation script (parseVariantFile.py)
    if not SKIP_VARANNO: 
        varAnnotation(groupDict, confFile,STRAND_BIAS,QSUB,VERBOSE,LIGHT,QUEUE)
    
    # produce run metrics
    produceVarMetrics(groupDict, confFile, STRAND_BIAS, QSUB, VERBOSE,LIGHT,QUEUE)
    
    # compare samples between and within groups
    compareGroups(confFile,STRAND_BIAS,QSUB,REMOVE_N_READS,INCLUDE_FEATURES,VERBOSE,LIGHT,QUEUE)
    print 'Analysis done..(%s)\n'%(time.time()-startTime)
    summaryReport.append('Analysis done..(%s)\n'%(time.time()-startTime))
    
    summaryFile = open('summary.txt','w')
    summaryFile.write('\n'.join(summaryReport))
    
    cleanup(groupDict,WEBSERVER)
    finalize()
    
if __name__=='__main__':
    try:
        analyze()
    except SystemExit:
        pass
    except Exception,e:
        import datetime
        import traceback
        errorFileName = 'viVanError.log'
#         errorFileName = datetime.datetime.now().strftime('viVanErrorLog-%y%m%d-%H%M%S.txt')
        errorFile = open(errorFileName,'w')
        errorFile.write(traceback.format_exc())
        print 'EXCEPTION: Analysis Failed. Exception printed into %s.\nTraceback:\n'%errorFileName
        print traceback.format_exc()
         
        
        
