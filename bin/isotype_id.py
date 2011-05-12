#! /usr/bin/env python

import sys
import warnings
import optparse

import seqtools
import vdj
import vdj.pipeline

parser = optparse.OptionParser()
parser.add_option('-i','--IGHC',dest='ighc_fasta')
(options, args) = parser.parse_args()

if len(args) == 2:
    inhandle = open(args[0],'r')
    outhandle = open(args[1],'w')
elif len(args) == 1:
    inhandle = open(args[0],'r')
    outhandle = sys.stdout
elif len(args) == 0:
    inhandle = sys.stdin
    outhandle = sys.stdout

# load isotypes
isotypes = vdj.pipeline.load_isotypes(options.ighc_fasta)

for chain in vdj.parse_imgt(inhandle):
    vdj.pipeline.id_isotype(chain,isotypes)
    print >>outhandle, chain
