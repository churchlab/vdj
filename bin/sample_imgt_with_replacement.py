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

"""sample_vdjxml_with_replacement.py

Subsample an input vdjxml stream. Sample with replacement.
"""

import sys
import random
import optparse
import subprocess

import statstools
import vdj

option_parser = optparse.OptionParser()
option_parser.add_option('-n','--num',type='int')
(options,args) = option_parser.parse_args()

if len(args) == 2:
    inhandle = open(args[0],'r')
    outhandle = open(args[1],'w')
elif len(args) == 1:
    inhandle = open(args[0],'r')
    outhandle = sys.stdout
elif len(args) == 0:
    raise ValueError, "must provide at least an input file"

# determine the total number of chains in the file (using unix grep and wc)
p = subprocess.Popen('cat %s | grep ^ID | wc -l' % args[0],shell=True,stdout=subprocess.PIPE)
total_chains = int(p.stdout.read().strip())

random.seed()
idxs = sorted(statstools.sample_with_replacement(xrange(total_chains),options.num))
for (i,chain) in enumerate(vdj.parse_imgt(inhandle)):
    if len(idxs) == 0:
        break
    if i == idxs[0]:
        while len(idxs) > 0 and i == idxs[0]:
            print >>outhandle, chain
            idxs.pop(0)

