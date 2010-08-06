#! /usr/bin/env python

import sys
import optparse

import vdj
import vdj.alignment

parser = optparse.OptionParser()
parser.add_option('-L','--locus',action='append',dest='loci')
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

aligner = vdj.alignment.vdj_aligner_combined(loci=options.loci)
for chain in vdj.parse_VDJXML(inhandle):
    aligner.coding_chain(chain)
    print >>outhandle, chain
