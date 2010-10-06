#! /usr/bin/env python

import sys
import optparse

import vdj

parser = optparse.OptionParser()
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

print >>outhandle, "<root>"
for chain in vdj.parse_VDJXML(inhandle):
    if hasattr(chain,'v') and hasattr(chain,'j'):
        print >>outhandle, chain
print >>outhandle, "</root>"
