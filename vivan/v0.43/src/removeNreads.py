#!/usr/bin/env python

import optparse
import re
import sys
import time

def checkFormat(fasFileName):
    firstLine = open(fasFileName,'r').readline()
    if firstLine.startswith('@'):
        return 'fastq'
    elif firstLine.startswith('>'):
        return 'fasta'
    else:
        raise "could not detect the file format"
            
def stats():
    parser = optparse.OptionParser()
    parser.add_option('-i','--inputFast') 
    parser.add_option('-o','--outputPrefix')
    args = parser.parse_args(sys.argv)[0]
    inputFast = args.inputFast

    outFast = open('%s.filter.fq'%args.outputPrefix,'w')
    FORMAT = checkFormat(inputFast)
    filteredFast = open('%s.Nreads.fq'%args.outputPrefix,'w')
    if FORMAT=='fastq':
        skipLines=3
    else:
        skipLines=1
    
    inFile = open(inputFast,'r')
    inLine = inFile.readline()
    filtered = 0
    total = 0
    while inLine:
        total+=1
        lines = []
        FILTER = False
        if re.match('[>@]',inLine):
            lines.append(inLine)
            inLine = inFile.readline()
            if 'N' in inLine:
                filtered+=1
                FILTER=True
            for i in range(skipLines):
                lines.append(inLine)
                inLine = inFile.readline()
            
        if FILTER:
            if divmod(filtered,100)[1]==0:
                print '\rFiltered %s reads out of %s'%(filtered,total),
                sys.stdout.flush()
            
            filteredFast.write(''.join(lines))
        else:
            outFast.write(''.join(lines))
    print 'filtered %s reads out of %s'%(filtered,total)
if __name__=='__main__':
    stats()
