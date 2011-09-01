#! /usr/bin/env python

import sys
import optparse

# HACK
# importing vdj.clustering will import scipy.cluster, which will import pylab,
# if it's available.  By default, pylab will try to open some kind of
# windowing device.  If I dispatch jobs to the batch system, this causes an
# error.  But if I preempt this import by calling matplotlib and connecting it
# to the non-windowing Agg backend, it should avoid the problem.
import matplotlib as mpl
mpl.use('agg')

import vdj
import vdj.clustering

parser = optparse.OptionParser()
parser.add_option('-c','--cutoff',type='float',default=4.5)
parser.add_option('-t','--tag',default='')
parser.add_option('-l','--linkage',type='choice',choices=['single','complete'],default='single')
(options, args) = parser.parse_args()

if len(args) == 2:
    inhandle = open(args[0],'r')
    outhandle = open(args[1],'w')
elif len(args) == 1:
    inhandle = open(args[0],'r')
    outhandle = sys.stdout
elif len(args) == 0:
    inhandle = sys.stdin
    outhandle = sys.stdout

# NOTE: this script requires there to be a well-defined junction
#       sequence.  It raises an exception if not.  Therefore, seqs
#       must be pre-filtered for having legit junctions
# NOTE: this script must hold all chains in memory in order to 
#       perform the clustering and then assign cluster names

# load data
chains = []
junctions = []
for chain in vdj.parse_imgt(inhandle):
    # check for presence of V, J, and non-trivial junction
    if not hasattr(chain,'v') or not hasattr(chain,'j') or not hasattr(chain,'junction'):
        raise ValueError, "Chain %s has no junction of V-J aln." % chain.descr
    chains.append(chain)
    junctions.append(chain.junction)

# perform the sequence clustering
(T,seq_idxs) = vdj.clustering.cluster_seqs(junctions,options.cutoff,options.linkage)

# tag chains with unique cluster IDs
if options.tag == '':
    tag = ''
else:
    tag = options.tag+'|'

for (i,chain) in enumerate(chains):
    cloneID = '%s%s' % (tag,T[seq_idxs[chain.junction]])
    chain.annotations['clone'] = cloneID
    print >>outhandle, chain
