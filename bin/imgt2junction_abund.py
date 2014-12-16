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
