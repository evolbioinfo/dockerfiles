#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

cranrepo=args[1]
rlib=args[2]
packages=args[3]

packageslist = read.table(packages,header=F)

for (l in packageslist) {

    install.packages(l, dependencies=TRUE,lib=rlib,repo=cranrepo);
    #if ( ! library(l, character.only=TRUE, logical.return=TRUE) ) {
    #    quit(status=1, save='no')
    #}
}
