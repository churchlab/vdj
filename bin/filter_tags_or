#! /usr/bin/env python

import sys
import optparse

import vdj

parser = optparse.OptionParser()
parser.add_option('-t','--tag',action='append',dest='tags')
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

empty_set = set()
query_tags = set(options.tags)
for chain in vdj.parse_VDJXML(inhandle):
    if query_tags & chain.all_tags != empty_set:
        print >>outhandle, chain
