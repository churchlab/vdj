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

import os
import argparse

import vdj
from pyutils import cleanup_id

argparser = argparse.ArgumentParser(description=None)
argparser.add_argument('input_file')
argparser.add_argument('output_dir',default=os.getcwd())
args = argparser.parse_args()

for chain in vdj.parse_imgt(args.input_file):
    output_file = os.path.join(args.output_dir,'%s.imgt' % cleanup_id(chain.id))
    with open(output_file,'w') as op:
        print >>op, chain
