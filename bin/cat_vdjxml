#! /usr/bin/env python

import sys
import optparse
import glob

import vdj

parser = optparse.OptionParser()
(options, args) = parser.parse_args()

files = []
for arg in args:
    files.extend(glob.glob(arg))

print >>sys.stdout, "<root>"
for f in files:
    inhandle = open(f,'r')
    for chain in vdj.parse_VDJXML(inhandle):
        print >>sys.stdout, chain
print >>sys.stdout, "</root>"
