#! /usr/bin/env python

import sys
import optparse

import vdj
import vdj.alignment

parser = optparse.OptionParser()
parser.add_option('-L','--locus',action='append',dest='loci')
parser.add_option('-D','--debug',action='store_true',default=False)
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
else:
    raise Exception, "Wrong number of arguments."

if options.debug:
    import pdb
    pdb.set_trace()

aligner = vdj.alignment.vdj_aligner_combined(loci=options.loci)
for chain in vdj.parse_imgt(inhandle):
    aligner.align_chain(chain)
    print >>outhandle, chain
