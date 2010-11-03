#! /usr/bin/env python

import sys

import vdj
import vdj.legacy

if len(sys.argv) == 3:
    inhandle = open(sys.argv[1],'r')
    outhandle = open(sys.argv[2],'w')
elif len(sys.argv) == 2:
    inhandle = open(sys.argv[1],'r')
    outhandle = sys.stdout
elif len(sys.argv) == 1:
    inhandle = sys.stdin
    outhandle = sys.stdout

for chain in vdj.legacy.parse_VDJXML_old(inhandle):
    print >>outhandle, chain
