#! /usr/bin/env python

import sys
import optparse

import seqtools

import vdj

parser = optparse.OptionParser()
parser.add_option('-b','--basename')
(options, args) = parser.parse_args()

if len(args) == 1:
    inhandle = open(args[0],'r')
elif len(args) == 0:
    inhandle = sys.stdin

# NOTE: this script ignores the allele numbers

partition_files = vdj.pipeline.partition_VJ(inhandle,options.basename)
for filename in partition_files:
    print filename
