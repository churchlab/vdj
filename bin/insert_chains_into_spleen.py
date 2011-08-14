#! /usr/bin/env python

import sys
import argparse

import pymongo

import vdj
import vdj.mongo

argparser = argparse.ArgumentParser(description=None)
argparser.add_argument('-d','--db',required=True)
argparser.add_argument('-c','--collection',default='chains')
argparser.add_argument('-i','--input')
# argparser.add_argument('--option',dest='xxx',action='store_const',default=5)
args = argparser.parse_args()

inputfile = args.input
db = vdj.mongo.connect_to_spleen(connect_to=args.db)
chains = db[args.collection]
for (i,chain) in enumerate(vdj.parse_imgt(inputfile)):
    if i%1000 == 0:
        sys.stdout.write("%i "%i)
        sys.stdout.flush()
    chains.insert(vdj.mongo.encode_chain(chain))
