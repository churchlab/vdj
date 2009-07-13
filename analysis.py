import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import vdj

def scatter_repertoires(reps,info='VJCDR3',gooddata=False,measurement='proportions'):
	"""Create a grid of scatter plots showing corelations between all pairs of repertoires.
	
	reps -- list of Repertoire objects
	
	"""
	numreps = len(reps)
	
	datalist = []
	for rep in reps:
		datalist.append( vdj.counts_ontology_1D(rep,info,gooddata) )
	
	print "finished computing counts"
	
	if measurement == 'proportions':
		for data in datalist:
			data = np.float_(data) / np.sum(data)
	
	print "finished normalizing data"
	
	min_nonzero = np.min([np.min(data[data>0]) for data in datalist])
	max_nonzero = np.max([np.max(data[data>0]) for data in datalist])
	axislo = 10**np.floor( np.frexp(min_nonzero)[1] * np.log10(2) )
	axishi = 10**np.ceil(  np.frexp(max_nonzero)[1] * np.log10(2) )
	
	fig = plt.figure()
	
	for row in xrange(1,numreps):
		for col in xrange(row,numreps):
			plotnum = (numreps-1)*(row-1) + col
			ax = fig.add_subplot(numreps-1,numreps-1,plotnum)
			ax.scatter(datalist[row-1],datalist[col-1],c='k',marker='o',s=2,edgecolors=None)
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.axis([axislo,axishi,axislo,axishi])
	
	return fig