import numpy as np
#import matplotlib
import matplotlib.pyplot as plt

import vdj

def scatter_repertoires_ontology(reps,info='VJCDR3',gooddata=False,measurement='proportions'):
	"""Create a grid of scatter plots showing corelations between all pairs of repertoires.
	
	reps -- list of Repertoire objects
	
	"""
	numreps = len(reps)
	
	datalist = []
	for rep in reps:
		datalist.append( vdj.counts_ontology_1D(rep,info,gooddata) )
	
	if measurement == 'proportions':
		for i in xrange(len(datalist)):
			datalist[i] = np.float_(datalist[i]) / np.sum(datalist[i])
	
	min_nonzero = np.min([np.min(data[data>0]) for data in datalist])
	max_nonzero = np.max([np.max(data[data>0]) for data in datalist])
	axislo = 10**np.floor( np.frexp(min_nonzero)[1] * np.log10(2) )
	axishi = 10**np.ceil(  np.frexp(max_nonzero)[1] * np.log10(2) )
	
	fig = plt.figure()
	
	hist_axs = []
	for row in xrange(numreps):
		col = row
		plotnum = numreps*row + col + 1
		ax = fig.add_subplot(numreps,numreps,plotnum)
		ax.hist(datalist[row],bins=100,log=True,facecolor='k')
		hist_axs.append(ax)
	
	scatter_axs = []
	for row in xrange(numreps-1):
		for col in xrange(row+1,numreps):
			plotnum = numreps*row + col + 1
			ax = fig.add_subplot(numreps,numreps,plotnum)
			ax.scatter(datalist[row],datalist[col],c='k',marker='o',s=2,edgecolors=None)
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.axis([axislo,axishi,axislo,axishi])
			scatter_axs.append(ax)
	
	return fig

def scatter_repertoires_clusters(reps,refclusters,measurement='proportions'):
	"""Create a grid of scatter plots showing corelations between all pairs of repertoires.
	
	reps -- list of Repertoire objects
	
	"""
	numreps = len(reps)
	
	datalist = []
	for rep in reps:
		clusters = vdj.getClusters(rep)
		datalist.append( vdj.countsClusters(clusters,refclusters) )
	
	if measurement == 'proportions':
		for i in xrange(len(datalist)):
			datalist[i] = np.float_(datalist[i]) / np.sum(datalist[i])
	
	min_nonzero = np.min([np.min(data[data>0]) for data in datalist])
	max_nonzero = np.max([np.max(data[data>0]) for data in datalist])
	axislo = 10**np.floor( np.frexp(min_nonzero)[1] * np.log10(2) )
	axishi = 10**np.ceil(  np.frexp(max_nonzero)[1] * np.log10(2) )
	
	fig = plt.figure()
	
	hist_axs = []
	for row in xrange(numreps):
		col = row
		plotnum = numreps*row + col + 1
		ax = fig.add_subplot(numreps,numreps,plotnum)
		ax.hist(datalist[row],bins=100,log=True,facecolor='k')
		hist_axs.append(ax)
	
	scatter_axs = []
	for row in xrange(numreps-1):
		for col in xrange(row+1,numreps):
			plotnum = numreps*row + col + 1
			ax = fig.add_subplot(numreps,numreps,plotnum)
			ax.scatter(datalist[row],datalist[col],c='k',marker='o',s=0.5,edgecolors=None)
			ax.set_xscale('log')
			ax.set_yscale('log')
			ax.axis([axislo,axishi,axislo,axishi])
			scatter_axs.append(ax)
	
	return fig

def timeseries_repertoires(times,reps,refclusters,allpositive=False):
	"""Create a time-series of the different clusters in refclusters.
	
	If allpositive is True, then it will limit itself to drawing timeseries
	only for those clusters that are non-zero at all timepoints.
	
	"""
	ax = plt.gca()
	
	numreps = len(reps)
	numclusters = len(refclusters)
	
	countdata = np.zeros((numclusters,numreps))
	for (i,rep) in enumerate(reps):
		clusters = vdj.getClusters(rep)
		countdata[:,i] = vdj.countsClusters(clusters,refclusters)
	
	sums = countdata.sum(0)
	proportions = np.float_(countdata) / sums
	
	ax.plot(times,countdata.transpose(),'k-',linewidth=0.2)
	
	plt.draw_if_interactive()
	return ax