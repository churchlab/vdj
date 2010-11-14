#! /usr/bin/env python

import optparse
import sys

import numpy as np

import vdj
import vdj.analysis
import timeseries

option_parser = optparse.OptionParser()
option_parser.add_option('-t','--times')
option_parser.add_option('-q','--quantify',choices=['clone','junction','v','j','vj','vdj'])
(options,args) = option_parser.parse_args()

if len(args) == 2:
    inhandle = open(args[0],'r')
    outhandle = open(args[1],'w')
elif len(args) == 1:
    inhandle = open(args[0],'r')
    outhandle = sys.stdout
else:
    raise ValueError, "Must give a single argument to vdjxml file"

# determine mapping between barcodes and times
timedict = {}
ip = open(options.times,'r')
for line in ip:
    timedict[line.split()[0]] = float(line.split()[1])
ip.close()

# compute counts
features = [options.quantify,'barcode']
(uniq_feature_values,countdict) = vdj.analysis.vdjxml2countdict(inhandle,features)
times = np.array([timedict[bc] for bc in uniq_feature_values['barcode']])
timesort = np.argsort(times)
times = times[timesort]
countmatrix = vdj.analysis.countdict2matrix(features,uniq_feature_values,countdict)
countmatrix = countmatrix[:,timesort]  # sort in increasing times

# write time series data to file
timeseries.write_timeseries(outhandle,labels=uniq_feature_values[options.quantify],times=times,matrix=countmatrix)
