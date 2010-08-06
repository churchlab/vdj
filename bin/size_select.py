#! /usr/bin/env python

import sys
import optparse

import vdj

parser = optparse.OptionParser()
parser.add_option('-m','--min',type='int',default=0)
parser.add_option('-M','--max',type='int',default=float('inf'))
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

for chain in vdj.parse_VDJXML(inhandle):
    if len(chain) >= options.min and len(chain) <= options.max:
        print >>outhandle, chain
