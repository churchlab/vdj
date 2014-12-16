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
import urllib2
from StringIO import StringIO

from Bio import SeqIO
from BeautifulSoup import BeautifulSoup

url = sys.argv[1]
page = urllib2.urlopen( url )
soup = BeautifulSoup( page )
pre = soup('pre')[-1]
data = pre.contents[-1].strip()+'\n'
record = SeqIO.parse(StringIO(data),'fasta').next()
group = record.description.split('|')[1][:4]
with open(group+'.fasta','w') as op:
    op.write(data)
