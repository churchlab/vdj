#! /usr/bin/env python

import sys
import argparse

import pymongo

import vdj
import vdj.mongo

argparser = argparse.ArgumentParser(description=None)
argparser.add_argument('-d','--db',required=True)
argparser.add_argument('-c','--collection',required=True)
argparser.add_argument('-i','--input',nargs='?',type=argparse.FileType('r'),default=sys.stdin)
argparser.add_argument('-p','--padding',type=int,default=0)
# argparser.add_argument('--option',dest='xxx',action='store_const',default=5)
args = argparser.parse_args()

inputfile = args.input
db = vdj.mongo.connect_to_lymph(connect_to=args.db)
chains = db[args.collection]
for (i,chain) in enumerate(vdj.parse_imgt(inputfile)):
    if i%1000 == 0:
        sys.stdout.write("%i "%i)
        sys.stdout.flush()
    doc = vdj.mongo.encode_chain(chain)
    if args.padding > 0:
        doc['__padding__'] = '0' * args.padding
    chains.insert(doc, safe=True)

# delete padding from elements
if args.padding > 0:
    chains.update({}, {"$unset" : {"__padding__" : 1}}, upsert=False, safe=True, multi=True)