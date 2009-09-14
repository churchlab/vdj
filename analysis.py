import numpy as np
import numpy.ma as ma
import scipy as sp
import scipy.stats
import scipy.special
import matplotlib as mpl
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

def reps2timeseries(reps,refclusters):
	"""Return time series matrix from list or Repertoire objs in chron order.
	
	reps is list of Repertoire objects
	refclusters is the master list of reference clusters
	"""
	numreps = len(reps)
	numclusters = len(refclusters)
	
	countdata = np.zeros((numclusters,numreps))
	for (i,rep) in enumerate(reps):
		clusters = vdj.getClusters(rep)
		countdata[:,i] = vdj.countsClusters(clusters,refclusters)
	
	return countdata

def timeseries_repertoires(times,reps,refclusters,idxsbool=None,allpositive=False):
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
	
	if idxsbool == None:
		if allpositive == True:
			idxsbool = np.sum(proportions,axis=1) > 0
		else:
			idxsbool = np.array([True]*proportions.shape[0])
	
	ax.plot(times,countdata[idxsbool,:].transpose(),'k-',linewidth=0.2)
	
	plt.draw_if_interactive()
	return ax

def rep2spectratype(rep):
	"""Compute spectratype curves from Repertoire object."""
	
	cdr3s = np.array([c.cdr3 for c in rep if c.junction != ''])
	min_raw_cdr3 = np.min(cdr3s)
	max_raw_cdr3 = np.max(cdr3s)
	min_cdr3 = np.int(np.ceil( min_raw_cdr3 / 3.) * 3)	# will be a nonzero mult of 3
	max_cdr3 = np.int(np.floor(max_raw_cdr3 / 3.) * 3)	# will be a mult of 3
	
	# bin the CDR3s lengths.  The first elt is rep zero len (and should be zero)
	# and the last bin always represents one greater than the biggest mult of 3
	binnedcdr3s = np.histogram(cdr3s,bins=np.arange(0,max_cdr3+2))[0]	# the +2 is due to the pecul. of np.hist.
	
	gaussians = []
	for cdr3len in np.arange(min_cdr3,max_raw_cdr3,3):
		totalcdr3s = np.sum(binnedcdr3s[cdr3len-1:cdr3len+2])
		goodcdr3s  = binnedcdr3s[cdr3len]
		if totalcdr3s == 0:
			continue
		mu = cdr3len
		x = cdr3len-0.5
		tail = (1 - (np.float(goodcdr3s)/totalcdr3s)) / 2.
		sigma = (x-mu) / (np.sqrt(2.)*sp.special.erfinv(2*tail-1))
		rv = sp.stats.norm(loc=mu,scale=sigma)
		gaussians.append( (totalcdr3s,rv) )
	
	t = np.linspace(0,max_cdr3+1,1000)
	y = np.zeros(len(t))
	for (s,rv) in gaussians:
		y += s*rv.pdf(t)
	return (t,y)

def circlemapVJ(ax,counts,rowlabels=None,collabels=None,log=False):
    numV = counts.shape[0]
    numJ = counts.shape[1]
    X,Y = np.meshgrid(range(numJ),range(numV))
    
    # mask zero positions
    X,Y = ma.array(X), ma.array(Y)
    C = ma.array(counts)
    
    zeromask = (counts == 0)
    X.mask = zeromask
    Y.mask = zeromask
    C.mask = zeromask
    
    # ravel nonzero elts (deletes zero-positions)
    x = ma.compressed(X)
    y = ma.compressed(Y)
    c = ma.compressed(C)
    
    # log normalize counts if requested
    if log == True:
        c = ma.log10(c)
    
    # normalize counts to desired size-range
    max_counts = ma.max(c)
    min_counts = ma.min(c)
    counts_range = max_counts - min_counts
    
    max_size = 100
    min_size = 5
    size_range = max_size - min_size
    
    sizes = (np.float(size_range) / counts_range) * (c - min_counts) + min_size
    
    collection = mpl.collections.CircleCollection(
                                        sizes,
                                        offsets = zip(x,y),
                                        transOffset = ax.transData) # i may need to explicitly set the xlim and ylim info for this to work correctly
    
    ax.add_collection(collection)
    
    return

# define colormap for -1 to 1 (green-black-red) like gene expression
redgreencdict = {'red':	[(0.0,   0.0,   0.0),
						 (0.5,   0.0,   0.0),
						 (1.0,   1.0,   0.0)],
						
				'green':[(0.0,   0.0,   1.0),
						 (0.5,   0.0,   0.0),
						 (1.0,   0.0,   0.0)],
						
				'blue': [(0.0,   0.0,   0.0),
						 (0.5,   0.0,   0.0),
						 (1.0,   0.0,   0.0)]}

redgreen = mpl.colors.LinearSegmentedColormap('redgreen',redgreencdict,256)
redgreen.set_bad(color='w')