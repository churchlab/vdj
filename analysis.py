import numpy as np
import scipy as sp
import scipy.stats
import scipy.special
import matplotlib as mpl
import matplotlib.pyplot as plt

import vdj

# import numpy.ma as ma


def clone_timeseries(inhandle,time_tags,reference_clones=None):
    # get count data
    time_tags_set = set(time_tags)
    clone_counts = {}
    for tag in time_tags:
        clone_counts[tag]={}
    for chain in vdj.parse_VDJXML(inhandle):
        try:
            curr_time_tag = (chain.tags&time_tags_set).pop()
        except KeyError:
            continue    
        
        try:
            clone_counts[curr_time_tag][vdj.get_clone(chain)] += 1
        except KeyError:
            clone_counts[curr_time_tag][vdj.get_clone(chain)] = 1
    
    # set up reference clones
    if reference_clones == None:
        reference_clones = set()
        for counts in clone_counts.itervalues():
            reference_clones.update(counts.keys())
        reference_clones = list(reference_clones)
    
    # build timeseries matrix
    num_clones = len(reference_clones)
    num_times = len(time_tags)
    countdata = np.zeros((num_clones,num_times))
    for (i,tag) in enumerate(time_tags):
        countdata[:,i] = vdj.count_dict_clone_counts(clone_counts[tag],reference_clones)
    
    return countdata


def cdr3s2spectratype(cdr3s):
    min_raw_cdr3 = np.min(cdr3s)
    max_raw_cdr3 = np.max(cdr3s)
    min_cdr3 = np.int(np.ceil( min_raw_cdr3 / 3.) * 3)  # will be a nonzero mult of 3
    max_cdr3 = np.int(np.floor(max_raw_cdr3 / 3.) * 3)  # will be a mult of 3
    
    # bin the CDR3s lengths.  The first elt is rep zero len (and should be zero)
    # and the last bin always represents one greater than the biggest mult of 3
    binnedcdr3s = np.histogram(cdr3s,bins=np.arange(0,max_cdr3+2))[0]   # the +2 is due to the pecul. of np.hist.
    
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


def spectratype_curves(inhandle):
    # NOTE: requires chains with V and J alns
    
    # init data structure
    cdr3s = {}
    for v_seg in vdj.refseq.IGHV_seqs.keys():
        for j_seg in vdj.refseq.IGHJ_seqs.keys():
            cdr3s[vdj.vj_id(v_seg,j_seg)] = []
    
    # load data
    for chain in vdj.parse_VDJXML(inhandle):
        if chain.v == '' or chain.j == '' or chain.junction == '':
            continue
        cdr3s[vdj.vj_id(chain.v,chain.j)].append(chain.cdr3)
    
    spectras = {}
    for v_seg in vdj.refseq.IGHV_seqs.keys():
        for j_seg in vdj.refseq.IGHJ_seqs.keys():
            if len(cdr3s[vdj.vj_id(v_seg,j_seg)]) == 0:
                # empty VJ combo:
                spectras[vdj.vj_id(v_seg,j_seg)] = (np.array([0,150]),np.array([0,0]))
            else:
                spectras[vdj.vj_id(v_seg,j_seg)] = cdr3s2spectratype(cdr3s[vdj.vj_id(v_seg,j_seg)])
    
    return spectras



# ========================
# = Diversity estimation =
# ========================


def estimator_chao1(counts):
    """Bias corrected.  See EstimateS doc (Colwell)"""
    Sobs = len(counts)
    F1 = np.float_(np.sum(np.int_(counts)==1))
    F2 = np.float_(np.sum(np.int_(counts)==2))
    chao1 = Sobs + F1*(F1-1)/(2*(F2+1))
    return chao1


def estimator_chao1_variance(counts):
    F1 = np.float_(np.sum(np.int_(counts)==1))
    F2 = np.float_(np.sum(np.int_(counts)==2))
    if F1 > 0 and F2 > 0:
        chao1_var = (F1*(F1-1)/(2*(F2+1))) + (F1*(2*F1-1)*(2*F1-1)/(4*(F2+1)*(F2+1))) + (F1*F1*F2*(F1-1)*(F1-1)/(4*(F2+1)*(F2+1)*(F2+1)*(F2+1)))
    elif F1 > 0 and F2 == 0:
        Schao1 = estimator_chao1(counts)
        chao1_var = (F1*(F1-1)/2) + (F1*(2*F1-1)*(2*F1-1)/4) - (F1*F1*F1*F1/(4*Schao1))
    elif F1 == 0:
        N = np.float_(np.sum(counts))
        Sobs = np.float_(len(counts))
        chao1_var = Sobs*np.exp(-1*N*Sobs) * (1-np.exp(-1*N*Sobs))
    return chao1_var


def estimator_ace(counts,rare_cutoff=10):
    Sobs = np.float_(len(counts))
    Srare = np.float_(np.sum(np.int_(counts)<=rare_cutoff))
    Sabund = Sobs - Srare
    F1 = np.float_(np.sum(np.int_(counts)==1))
    F = lambda i: np.sum(np.int_(counts)==i)
    Nrare = np.float_(np.sum([i*F(i) for i in range(1,rare_cutoff+1)]))
    #Nrare = np.float_(np.sum(counts[np.int_(counts)<=rare_cutoff]))
    Cace = 1 - (F1/Nrare)
    gamma_squared = Srare*np.sum([i*(i-1)*F(i) for i in range(1,rare_cutoff+1)])/(Cace*Nrare*(Nrare-1))
    if gamma_squared < 0:
        gamma_squared = 0
    Sace = Sabund + (Srare/Cace) + (F1/Cace)*gamma_squared
    return Sace



# =========================
# = Statistical utilities =
# =========================

def counts2sample(counts):
    """Computes a consistent sample from a vector of counts.
    
    Takes a vector of counts and returns a vector of indices x
       such that len(x) = sum(c) and each elt of x is the index of
       a corresponding elt in c
    """
    x = np.ones(np.sum(counts),dtype=np.int_)
    
    start_idx = 0
    end_idx = 0
    for i in xrange(len(counts)):
        start_idx = end_idx
        end_idx = end_idx + counts[i]
        x[start_idx:end_idx] = x[start_idx:end_idx] * i 
    return x


def sample2counts(sample, categories=0):
    """Return count vector from list of samples.
    
    Take vector of samples and return a vector of counts.  The elts
       refer to indices in something that would ultimately map to the
       originating category (like from a multinomial).  Therefore, if there
       are, say, 8 categories, then valid values in sample should be 0-7.
       If categories is not given, then i compute it from the highest value
       present in sample (+1).
    """
    counts = np.bincount(sample)
    if (categories > 0) and (categories > len(counts)):
        counts = np.append( counts, np.zeros(categories-len(counts)) )
    return counts


def scoreatpercentile(values,rank):
    return sp.stats.scoreatpercentile(values,rank)


def percentileofscore(values,score):
    values.sort()
    return values.searchsorted(score) / np.float_(len(values))


def bootstrap(x, nboot, theta, *args):
    '''return n bootstrap replications of theta from x'''
    N = len(x)
    th_star = np.zeros(nboot)
    
    for i in xrange(nboot):
        th_star[i] = theta( x[ np.random.randint(0,N,N) ], *args )    # bootstrap repl from x
    
    return th_star


def subsample(x, num_samples, sample_size, theta, *args):
    """return num_samples evaluations of the statistic theta
    on subsamples of size sample_size"""
    N = len(x)
    th_star = np.zeros(num_samples)
    
    for i in xrange(num_samples):
        th_star[i] = theta( x[ np.random.randint(0,N,sample_size) ], *args )    # subsample from from x
    
    return th_star


def randint_without_replacement(low,high=None,size=None):
    if high == None:
        high = low
        low = 0
    if size == None:
        size = 1
    urn = range(low,high)
    N = len(urn)
    flip = False
    if size > N/2:
        flip = True
        size = N - size
    sample = []
    for i in xrange(size):
        draw = np.random.randint(0,N-i)
        sample.append(urn.pop(draw))
    if not flip:
        return np.asarray(sample)
    else:
        return np.asarray(urn)


def subsample_without_replacement(x, num_samples, sample_size, theta, *args):
    """return num_samples evaluations of the statistic theta
    on subsamples of size sample_size"""
    N = len(x)
    th_star = np.zeros(num_samples)
    
    for i in xrange(num_samples):
        th_star[i] = theta( x[ randint_without_replacement(0,N,sample_size) ], *args )    # subsample from from x
    
    return th_star










#==============================================================================
#==============================================================================
#==============================================================================

def rep2spectratype(rep):
    """Compute spectratype curves from Repertoire object."""
    
    cdr3s = np.array([c.cdr3 for c in rep if c.junction != ''])
    min_raw_cdr3 = np.min(cdr3s)
    max_raw_cdr3 = np.max(cdr3s)
    min_cdr3 = np.int(np.ceil( min_raw_cdr3 / 3.) * 3)  # will be a nonzero mult of 3
    max_cdr3 = np.int(np.floor(max_raw_cdr3 / 3.) * 3)  # will be a mult of 3
    
    # bin the CDR3s lengths.  The first elt is rep zero len (and should be zero)
    # and the last bin always represents one greater than the biggest mult of 3
    binnedcdr3s = np.histogram(cdr3s,bins=np.arange(0,max_cdr3+2))[0]   # the +2 is due to the pecul. of np.hist.
    
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
    min_cdr3 = np.int(np.ceil( min_raw_cdr3 / 3.) * 3)  # will be a nonzero mult of 3
    max_cdr3 = np.int(np.floor(max_raw_cdr3 / 3.) * 3)  # will be a mult of 3
    
    # bin the CDR3s lengths.  The first elt is rep zero len (and should be zero)
    # and the last bin always represents one greater than the biggest mult of 3
    binnedcdr3s = np.histogram(cdr3s,bins=np.arange(0,max_cdr3+2))[0]   # the +2 is due to the pecul. of np.hist.
    
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
    
    ax.set_aspect('equal')
    ax.autoscale_view()
    
    ax.xaxis.set_major_locator(mpl.ticker.FixedLocator(range(counts.shape[1])))
    ax.yaxis.set_major_locator(mpl.ticker.FixedLocator(range(counts.shape[0])))
    
    if rowlabels != None:
        ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter(collabels))
    if collabels != None:
        ax.yaxis.set_major_formatter(mpl.ticker.FixedFormatter(rowlabels))
    
    for ticklabel in ax.xaxis.get_ticklabels():
        ticklabel.set_horizontalalignment('left')
        ticklabel.set_rotation(-45)
        ticklabel.set_size(8)
    
    for ticklabel in ax.yaxis.get_ticklabels():
        ticklabel.set_size(8)
    
    return

# define colormap for -1 to 1 (green-black-red) like gene expression
redgreencdict = {'red': [(0.0,   0.0,   0.0),
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