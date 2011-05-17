#! /usr/bin/env python

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
