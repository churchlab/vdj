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

parts = []
chains_processed = 0
file_num = 0

curr_outname = options.basename+'.'+str(file_num)
for chain in vdj.parse_VDJXML(inhandle):
    if chains_processed == 0:
        op = open(curr_outname,'w')
        print >>op, "<root>"
        parts.append(curr_outname)
    
    print >>op, chain
    chains_processed += 1
    
    if chains_processed == options.packetsize:
        print >>op, "</root>"
        op.close()
        chains_processed = 0
        file_num += 1
        curr_outname = options.basename+'.'+str(file_num)
if not op.closed:
    print >>op, "</root>"
    op.close()

for part in parts:
    print part
