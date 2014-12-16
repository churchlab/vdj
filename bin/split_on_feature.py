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
import optparse

import vdj
import seqtools

option_parser = optparse.OptionParser()
option_parser.add_option('-f','--feature')
option_parser.add_option('-o','--outputbasename')
(options,args) = option_parser.parse_args()

if len(args) == 1:
    inhandle = open(args[0],'r')
elif len(args) == 0:
    inhandle = sys.stdin
else:
    raise ValueError, "must give a single input file as an argument or stdin"

feature = options.feature
cleanup_id = seqtools.cleanup_id

outhandles = {}
for chain in vdj.parse_imgt(inhandle):
    try: curr_feature = chain.__getattribute__(feature)
    except AttributeError: continue
    
    try: print >>outhandles[curr_feature], chain
    except KeyError:
        outname = '.'.join([options.outputbasename,cleanup_id(curr_feature),'imgt'])
        outhandles[curr_feature] = open(outname,'w')
        print >>outhandles[curr_feature], chain

for outhandle in outhandles.itervalues():
    outhandle.close()
