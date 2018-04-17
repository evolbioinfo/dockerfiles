#!/usr/bin/env python

from Bio import SeqIO
import optparse
import sys
import re

#===============================================================================
# Convert file format v0.1
# using biopython, converts the file -i that has the format -f into format -t
#===============================================================================

def convert():
    parser = optparse.OptionParser(description='get a file and its format and convert it to the given format')
    parser.add_option('-i','--inputFile')
    parser.add_option('-p','--prefix',default=None)
    parser.add_option('-f','--inputFormat')
    parser.add_option('-t','--outputFormat')
    
    args = parser.parse_args(sys.argv)[0]
    
    if args.prefix:
        inputPrefix = args.prefix
    else:
        inputPrefix = re.search('[^\.]+',args.inputFile).group(0)
    
    outputFile = open(inputPrefix+'.%s'%args.outputFormat,'w')
    
    SeqIO.convert(args.inputFile, args.inputFormat, outputFile, args.outputFormat)
    
    

if __name__=='__main__':
    convert()
