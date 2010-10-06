#! /usr/bin/env python

import sys
import optparse
import glob

import vdj.pipeline

parser = optparse.OptionParser()
(options, args) = parser.parse_args()

files = []
for arg in args:
    files.extend(glob.glob(arg))

vdj.pipeline.cat_vdjxml(files,sys.stdout)
