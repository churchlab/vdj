#! /usr/bin/env python

import sys
import optparse
import os

import vdj

parser = optparse.OptionParser()
parser.add_option('-t','--tag',action='append',dest='tags')
(options, args) = parser.parse_args()

if len(args) == 1:
    inname = args[0]
    inhandle = open(args[0],'r')
else:
    raise Exception, "Need a single input file."

(basename,ext) = os.path.splitext(inname)
basename = os.path.basename(basename)

outhandles = {}
for tag in options.tags:
    outname = basename+'.'+vdj.sequtils.cleanup_id(tag)+ext
    outhandles[tag] = open(outname,'w')

query_tags = set(options.tags)
for chain in vdj.parse_VDJXML(inhandle):
    try:
        tag = (query_tags & chain.all_tags).pop()
        print >>outhandles[tag], chain
    except KeyError:
        continue

for handle in outhandles.itervalues():
    handle.close()
