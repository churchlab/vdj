#! /usr/bin/env python

import sys
import optparse

import vdj
import seqtools

option_parser = optparse.OptionParser()
option_parser.add_option('-f','--feature')
option_parser.add_option('-o','--outputbasename')
(options,args) = option_parser.parse_args()

if len(args) == 1:
    inhandle = open(args[0],'r')
elif len(args) == 0:
    inhandle = sys.stdin
else:
    raise ValueError, "must give a single input file as an argument or stdin"

feature = options.feature
cleanup_id = seqtools.cleanup_id

outhandles = {}
for chain in vdj.parse_imgt(inhandle):
    try: curr_feature = chain.__getattribute__(feature)
    except AttributeError: continue
    
    try: print >>outhandles[curr_feature], chain
    except KeyError:
        outname = '.'.join([options.outputbasename,cleanup_id(curr_feature),'imgt'])
        outhandles[curr_feature] = open(outname,'w')
        print >>outhandles[curr_feature], chain

for outhandle in outhandles.itervalues():
    outhandle.close()
