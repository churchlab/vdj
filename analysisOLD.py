# analysis.py

# VDJ package for analysis and visualization functions

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

#####################################################################
#####################################################################
##
##
##		CLEAN UP TO HERE
##
##
#####################################################################
#####################################################################





# define hot colormap for counts.  it can be log-normalized
hotcounts = matplotlib.colors.LinearSegmentedColormap('hotcounts',matplotlib.cm.datad['hot'],256)
#hotsafe.set_bad(color='#3E623E')
hotcounts.set_bad(color='#000060')
hotcounts.set_over(color='w')
hotcounts.set_under(color='k')
#hotsafe.set_under(color='#3333FF')

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

redgreen = matplotlib.colors.LinearSegmentedColormap('redgreen',redgreencdict,256)
redgreen.set_bad(color='w')

def heatmap(repertoire,info='VJCDR3',scale='log',measurement='proportions'):
	if info == 'VJ':
		data = ma.array(np.float_(repertoire.countsVJ()))
	elif info == 'VDJ':
		data = ma.array(np.float_(repertoire.countsVDJ_2D()))
	elif info == 'VJCDR3':
		data = ma.array(np.float_(repertoire.countsVJCDR3_2D()))
	
	if measurement == 'proportions':
		data = data / np.sum(data)
	
	if scale == 'log':
		data = np.log10(data)
	
	data.mask = np.logical_not(np.isfinite(data))
	
	return pylab.imshow(data,interpolation='nearest',cmap=hotcounts)

# def heatmapVDJ(repertoire,scale='log',measurement='counts'):
# 	data = ma.array(np.float_(repertoire.countsVDJ_2D()))
# 	if measurement == 'proportions':
# 		data = data / np.sum(data)
# 	if scale == 'log':
# 		data = np.log10(data)
# 	data.mask = np.logical_not(np.isfinite(data))
# 	return pylab.imshow(data,interpolation='nearest',cmap=hotcounts)
# 
# def heatmapVJ(repertoire,scale='log'):
# 	if scale == 'linear':
# 		scaling = None
# 	elif scale == 'log':
# 		scaling = matplotlib.colors.LogNorm(vmin=1)
# 	data = repertoire.countsVJ()
# 	return pylab.imshow(data,interpolation='nearest',cmap=hotcounts,norm=scaling)

def scatter(rep1,rep2,p=0,info='VJCDR3',gooddata=True,measurement='proportions'):
	if info == 'VJ':
		if gooddata == True:
			counts1 = rep1.get_rep_fullVJ().countsVJ_1D()
			counts2 = rep2.get_rep_fullVJ().countsVJ_1D()
		else:
			counts1 = rep1.countsVJ_1D()
			counts2 = rep2.countsVJ_1D()
	elif info == 'VDJ':
		if gooddata == True:
			counts1 = rep1.get_rep_fullVDJ().countsVDJ_1D()
			counts2 = rep2.get_rep_fullVDJ().countsVDJ_1D()
		else:
			counts1 = rep1.countsVDJ_1D()
			counts2 = rep2.countsVDJ_1D()
	elif info == 'VJCDR3':
		if gooddata == True:
			counts1 = rep1.get_rep_fullVJCDR3().countsVJCDR3_1D()
			counts2 = rep2.get_rep_fullVJCDR3().countsVJCDR3_1D()
		else:
			counts1 = rep1.countsVJCDR3_1D()
			counts2 = rep2.countsVJCDR3_1D()
	
	data1 = ma.array(np.float_(counts1))
	data2 = ma.array(np.float_(counts2))
	
	if measurement == 'proportions':
		data1 = data1 / np.sum(data1)
		data2 = data2 / np.sum(data2)
	
	if p==0:	# no signif calc
		return pylab.scatter(data1,data2)
	
	pvals = countdata.pvalsCounts(counts1,counts2)
	signif = (pvals <= p)
	notsignif = np.logical_not( signif )

	ax = pylab.scatter(data1[notsignif],data2[notsignif],c='b')
	pylab.scatter(data1[signif],data2[signif],c='r')
	
	return ax

def scatter_repertoires(reps,info='VJCDR3',gooddata=False,measurement='proportions'):
	"""Create grid of scatter plots showing correlations between all pairs of repertoires.
	
	reps -- list of Repertoire objects
	
	"""
	numreps = len(reps)
	
	datalist = []
	for rep in reps:
		datalist.append( vdj.counts_ontology_1D(rep,info,gooddata) )
	
	if measurement == 'proportions':
		for data in datalist:
			data = np.float_(data) / np.sum(data)

# def scatterVJCDR3(rep1,rep2,p=0,properCDR3=True):
# 	# see if I should only get chains with good CDR3s
# 	if properCDR3 == True:
# 		counts1 = rep1.get_rep_fullVJCDR3().countsVJCDR3_1D()
# 		counts2 = rep2.get_rep_fullVJCDR3().countsVJCDR3_1D()
# 	else:
# 		counts1 = rep1.countsVJCDR3_1D()
# 		counts2 = rep2.countsVJCDR3_1D()
# 	
# 	if p == 0:	# no significance calc
# 		return pylab.scatter(counts1,counts2)
# 	
# 	# compute p-values
# 	pvals = countdata.pvalsCounts(counts1,counts2)
# 	signif = (pvals <= p)
# 	notsignif = np.logical_not( signif )
# 	
# 	fig = pylab.scatter(counts1[notsignif],counts2[notsignif],c='b')
# 	pylab.scatter(counts1[signif],counts2[signif],c='r')
# 	
# 	return fig
# 
# def scatterVDJ(rep1,rep2,p=0,fullVDJ=True):
# 	# see if I should only use chains that have full V, D, and J alignment
# 	if fullVDJ == True:
# 		counts1 = rep1.get_rep_fullVDJ().countsVDJ_1D()
# 		counts2 = rep2.get_rep_fullVDJ().countsVDJ_1D()
# 	else:
# 		counts1 = rep1.countsVDJ_1D()
# 		counts2 = rep2.countsVDJ_1D()
# 	
# 	if p == 0:	# no significance calc
# 		return pylab.scatter(counts1,counts2)
# 	
# 	# compute p-values
# 	pvals = countdata.pvalsCounts(counts1,counts2)
# 	signif = (pvals <= p)
# 	notsignif = np.logical_not( signif )
# 	
# 	fig = pylab.scatter(counts1[notsignif],counts2[notsignif],c='b')
# 	pylab.scatter(counts1[signif],counts2[signif],c='r')
# 	
# 	return fig
# 
# def scatterVJ(rep1,rep2,p=0,fullVJ=True):
# 	# see if I should only use chains that have full V and J alignment
# 	if fullVJ == True:
# 		counts1 = rep1.get_rep_fullVJ().countsVJ_1D()
# 		counts2 = rep2.get_rep_fullVJ().countsVJ_1D()
# 	else:
# 		counts1 = rep1.countsVJ_1D()
# 		counts2 = rep2.countsVJ_1D()
# 	
# 	if p == 0:	# no significance calc
# 		return pylab.scatter(counts1,counts2)
# 	
# 	# compute p-values
# 	pvals = countdata.pvalsCounts(counts1,counts2)
# 	signif = (pvals <= p)
# 	notsignif = np.logical_not( signif )
# 	
# 	fig = pylab.scatter(counts1[notsignif],counts2[notsignif],c='b')
# 	pylab.scatter(counts1[signif],counts2[signif],c='r')
# 	
# 	return fig

def sortedpvals(rep1,rep2,info='VDJ'):
	'''
	takes two repertoire objects, and computes the pvals of the
	components, when counts are measured according to
	info:
		VJ
		VDJ
		VJCDR3
	it returns a sorted list (lowest pval first) of 2-tuples: (component,pval)
	where component is a list of either:
		[V,J]
		[V,D,J]
		[V,J,CDR3 len]
	'''
	
	if info == 'VJ':
		counts1 = rep1.countsVJ_1D()
		counts2 = rep2.countsVJ_1D()
		componentref = refseq.VJref
	elif info == 'VDJ':
		counts1 = rep1.countsVDJ_1D()
		counts2 = rep2.countsVDJ_1D()
		componentref = refseq.VDJref
	elif info == 'VJCDR3':
		counts1 = rep1.countsVJCDR3_1D()
		counts2 = rep2.countsVJCDR3_1D()
		componentref = refseq.VJCDR3ref(counts1)
	
	pvals = countdata.pvalsCounts(counts1,counts2)
	
	if len(pvals) != len(componentref):
		raise ValueError, 'number of components in pval is different from reference raveled array'
	
	sortedpvals = zip(componentref,pvals)
	sortedpvals.sort(key=lambda x:x[1])
	
	return sortedpvals

def heatmappvals(rep1,rep2,info='VJCDR3'):
	
	if info == 'VJ':
		counts1 = rep1.countsVJ()
		counts2 = rep2.countsVJ()
	elif info == 'VDJ':
		counts1 = rep1.countsVDJ_2D()
		counts2 = rep2.countsVDJ_2D()
	elif info == 'VJCDR3':
		counts1 = rep1.countsVJCDR3_2D()
		counts2 = rep2.countsVJCDR3_2D()
	
	pvals = countdata.pvalsCounts(counts1.ravel(),counts2.ravel())
	pvals.resize(counts1.shape)
	
	return pylab.imshow(1-pvals,interpolation='nearest',cmap=pylab.hot())

def heatmaplogratio(rep1,rep2,info='VJCDR3'):
	if info == 'VJ':
		counts1 = np.float_(rep1.countsVJ())
		counts2 = np.float_(rep2.countsVJ())
	elif info == 'VDJ':
		counts1 = np.float_(rep1.countsVDJ_2D())
		counts2 = np.float_(rep2.countsVDJ_2D())
	elif info == 'VJCDR3':
		counts1 = np.float_(rep1.countsVJCDR3_2D())
		counts2 = np.float_(rep2.countsVJCDR3_2D())
	
	proportions1 = counts1 / np.sum(counts1)
	proportions2 = counts2 / np.sum(counts2)	
	
	logratios = ma.array( np.log10(proportions2/proportions1) )
	logratios.mask = np.logical_not(np.isfinite(logratios))
	
	
	return pylab.imshow(logratios,interpolation='nearest',cmap=redgreen)

def circlemapVJ(rep):
	counts = vdj.countsVJ(rep)
	logcounts = np.log10(counts)
	logcounts[logcounts==-np.Inf] = 0
	norm_counts = counts / np.float_(np.max(counts))
	norm_logcounts = logcounts / np.float_(np.max(logcounts))
	
	num_rows = counts.shape[0]
	num_cols = counts.shape[1]
	fig_unit_frac = 1./(num_rows+2)
	axes_unit_frac = 1./num_rows
	axes_width  = num_cols*fig_unit_frac
	axes_height = num_rows*fig_unit_frac
	max_circle_radius = 0.98*axes_unit_frac / 2
	max_circle_area = np.pi * (max_circle_radius**2)
	
	circle_centers = np.empty(counts.shape,dtype=np.object_)
	for row in range(num_rows):
		for col in range(num_cols):
			x = (col*axes_unit_frac) + axes_unit_frac/2
			y = ((num_rows-row-1)*axes_unit_frac) + axes_unit_frac/2
			circle_centers[row,col] = (x,y)
	
	fig = plt.figure()
	ax = fig.add_axes([ (1-axes_width)/2, fig_unit_frac, axes_width, axes_height ])
	
	for row in range(num_rows):
		for col in range(num_cols):
			center = circle_centers[row,col]
			radius = norm_logcounts[row,col] * max_circle_radius
			circ = matplotlib.patches.Circle( center, radius, transform = ax.transAxes, linewidth=0, facecolor='b' )
			ax.add_patch(circ)
	
	#ax.draw()
	return ax