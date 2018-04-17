#!/usr/bin/env python

import optparse
import re
import sys
import time
from os import chdir
import os
import itertools
from os import popen

VERSION = 0.1

def getQFormat():
    
    print 'MAIN:NAME: Get Quality Format v%s'%VERSION
    parser = optparse.OptionParser()
    
    parser.add_option('-f','--fileName')
    parser.add_option('-n','--numOQuals2check',default=10000)
    args = parser.parse_args(sys.argv)[0]
    
    inFile = open(args.fileName,'r')
    
    numOquals2check = int(args.numOQuals2check)
    numOqualsChecked = 0
    
#    my $sanger_regexp = qr/[!"#$%&'()*+,-.\/0123456789:]/;
#    my $solexa_regexp = qr/[\;<=>\?]/;
#    my $solill_regexp = qr/[JKLMNOPQRSTUVWXYZ\[\]\^\_\`abcdefgh]/;
#    my $all_regexp = qr/[\@ABCDEFGHI]/;

    
    format = ""
    # set regular expressions
    sanger_regexp = re.escape('''!"#$%&'()*+,-./0123456789:''')
    solexa_regexp = re.escape(''';<=>?''')
    solill_regexp = re.escape('''JKLMNOPQRSTUVWXYZ[]^_`abcdefgh''')
    all_regexp = re.escape('''@ABCDEFGHI''')
    
    # set counters
    sanger_counter = 0
    solexa_counter = 0
    solill_counter = 0
    
    i =0
    print 'MAIN:PROC: Collecting Quals'
    for line in inFile:
        i+=1
        # retrieve qualities
        if divmod(i,4)[1]!=0:
            continue
        
        numOqualsChecked+=1
        
        # check qualities
        if( re.search('[%s]'%sanger_regexp,line) ):
            sanger_counter = 1
            break
        if( re.search('[%s]'%solexa_regexp,line )):
            solexa_counter = 1
        if( re.search('[%s]'%solill_regexp,line) ):
            solill_counter = 1;
        
        if numOqualsChecked==numOquals2check:
            break
        
    # determine format
    if( sanger_counter ):
        format = "sanger"

    elif(  not (sanger_counter) and solexa_counter ):
        format = "solexa"

    elif( not (sanger_counter) and  not (solexa_counter) and solill_counter ):
        format = "illumina"

    
    print "MAIN:OUT: FORMAT:%s"%format
if __name__=='__main__':
    getQFormat()
