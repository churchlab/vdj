#! /usr/bin/env python

import sys
import optparse

import vdj
import vdj.pipeline

parser = optparse.OptionParser()
parser.add_option('-b','--basename')
parser.add_option('-p','--packetsize',type='int')
(options, args) = parser.parse_args()

if len(args) == 1:
    inhandle = open(args[0],'r')
elif len(args) == 0:
    inhandle = sys.stdin
else:
    raise Exception, "Too many arguments."

parts = vdj.pipeline.iterator2parts( vdj.parse_VDJXML(inhandle),
                                     options.basename,
                                     options.packetsize,
                                     prefix='<root>',
                                     suffix='</root>')

for part in parts:
    print part
