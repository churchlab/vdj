#! /usr/bin/env python

import sys
import argparse
from collections import defaultdict

import vdj

argparser = argparse.ArgumentParser(description=None)
argparser.add_argument('input',nargs='?',type=argparse.FileType('r'),default=sys.stdin)
argparser.add_argument('output',nargs='?',type=argparse.FileType('w'),default=sys.stdout)
args = argparser.parse_args()

# read in all the junctions
junctions = defaultdict(list)
for chain in vdj.parse_imgt(args.input):
    try:
        junctions[chain.junction_nt].append(chain.id)
    except AttributeError:
        pass

for junction in sorted(junctions.iterkeys(), key=lambda k: len(junctions[k]), reverse=True):
    for id_ in junctions[junction]:
        args.output.write(">%s\n%s\n" % (id_, junction))
