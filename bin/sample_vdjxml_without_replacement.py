#! /usr/bin/env python
"""sample_vdjxml_without_replacement.py

Subsample an input vdjxml stream. Sample without replacement. If the sampling
level is larger than the number of chains in the file, it simply copies the
file and issues a warning. 
"""

import sys
import random
import optparse
import subprocess
import warnings

import vdj

warnings.simplefilter('always')

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

# check if subsampling should include the entire file
if options.num >= total_chains:
    warnings.warn("Subsampling level is greater than or equal to number of chains in file: printing whole file to output.")
    for chain in vdj.parse_VDJXML(inhandle):
        print >>outhandle, chain
else:
    # choose a random set of indices to select for the output file
    random.seed()
    idxs = sorted(random.sample(xrange(total_chains),options.num))
    for (i,chain) in enumerate(vdj.parse_VDJXML(inhandle)):
        if len(idxs) == 0:
            break
        if i == idxs[0]:
            print >>outhandle, chain
            idxs.pop(0)

print >>outhandle, "</root>"
