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

counts_dict = vdj.get_clone_counts(inhandle)
counts = count_dict_clone_counts(counts_dict)

for count in counts:
    print >>outhandle, count
