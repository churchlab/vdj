#! /usr/bin/env python

import sys
import optparse

import vdj

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

parts = vdj.split_vdjxml_into_parts(options.packetsize,inhandle,options.basename)

for part in parts:
    print part
