#! /usr/bin/env python
"""sample_vdjxml_with_replacement.py

Subsample an input vdjxml stream. Sample with replacement.
"""

import sys
import random
import optparse
import subprocess

import statstools
import vdj

option_parser = optparse.OptionParser()
option_parser.add_option('-n','--num',type='int')
(options,args) = option_parser.parse_args()

if len(args) == 2:
    inhandle = open(args[0],'r')
    outhandle = open(args[1],'w')
elif len(args) == 1:
    inhandle = open(args[0],'r')
    outhandle = sys.stdout
elif len(args) == 0:
    raise ValueError, "must provide at least an input file"

# determine the total number of chains in the file (using unix grep and wc)
p = subprocess.Popen('cat %s | grep "<ImmuneChain>" | wc -l' % args[0],shell=True,stdout=subprocess.PIPE)
total_chains = int(p.stdout.read().strip())

print >>outhandle, "<root>"

random.seed()
idxs = sorted(statstools.sample_with_replacement(xrange(total_chains),options.num))
for (i,chain) in enumerate(vdj.parse_VDJXML(inhandle)):
    if len(idxs) == 0:
        break
    if i == idxs[0]:
        while i == idxs[0]:
            print >>outhandle, chain
            idxs.pop(0)

print >>outhandle, "</root>"
