#! /usr/bin/env python

import sys

import vdj

infile = sys.argv[1]
outfile = sys.argv[2]

with open(outfile,'w') as op:
    for chain in vdj.parse_imgt(infile):
        try:
            op.write( ">%s\n%s\n" % (chain.id,chain.junction_nt) )
        except AttributeError:
            pass

