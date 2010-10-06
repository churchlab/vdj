#! /usr/bin/env python

import sys
import optparse

import seqtools

import vdj

parser = optparse.OptionParser()
parser.add_option('-b','--basename')
(options, args) = parser.parse_args()

if len(args) == 1:
    inhandle = open(args[0],'r')
elif len(args) == 0:
    inhandle = sys.stdin

# NOTE: this script ignores the allele numbers

def vj_id_no_allele(chain):
    return seqtools.cleanup_id(chain.v.split('*')[0]) + '_' + seqtools.cleanup_id(chain.j.split('*')[0])

outhandles = {}
for chain in vdj.parse_VDJXML(inhandle):
    curr_vj_id = vj_id_no_allele(chain)
    try:
        print >>outhandles[curr_vj_id], chain
    except KeyError:
        outhandles[curr_vj_id] = open("%s.%s.vdjxml" % (options.basename,curr_vj_id),'w')
        print >>outhandles[curr_vj_id], "<root>"
        print >>outhandles[curr_vj_id], chain

for outhandle in outhandles.itervalues():
    print >>outhandle, "</root>"
