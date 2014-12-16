#! /usr/bin/env python
# Copyright 2014 Uri Laserson
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

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