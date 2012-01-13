#! /usr/bin/env python

import sys
import argparse

import vdj

argparser = argparse.ArgumentParser(description=None)
argparser.add_argument('positional',nargs='+')
args = argparser.parse_args()

if len(args.positional) == 2:
    inhandle = open(args.positional[0],'r')
    outhandle = open(args.positional[1],'w')
elif len(args.positional) == 1:
    inhandle = open(args.positional[0],'r')
    outhandle = sys.stdout
elif len(args.positional) == 0:
    inhandle = sys.stdin
    outhandle = sys.stdout

for chain in vdj.parse_imgt(inhandle):
    try:
        outhandle.write( ">%s\n%s\n" % (chain.id,chain.junction_nt) )
    except AttributeError:
        pass

