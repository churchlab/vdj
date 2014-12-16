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

import itertools

import numpy as np
# from numpy import ma
# import scipy as sp
# import scipy.stats
# import scipy.special
# import matplotlib as mpl
# import matplotlib.pyplot as plt

import pyutils
import vdj

def iterator2countdict(iterable,features,count='read'):
    counts = pyutils.nesteddict()
    uniq_feature_values = dict([(f,set()) for f in features])
    for chain in iterable:
        try:    # get the feature tuple
            feature_list = [chain.__getattribute__(f) for f in features]
            for (feature,value) in zip(features,feature_list): uniq_feature_values[feature].add(value)
        except AttributeError:  # chain is missing a feature; abandon it
            continue
        
        # update the dictionary
        if count == 'read':
            counts.nested_increment(feature_list)
        elif count in ['junction','clone']:
            counts.nested_add(feature_list,chain.__getattribute__(count))
        else:
            raise ValueError, "'count' must be 'read', 'junction', or 'clone'"
        
    # if counting clones/junctions, convert sets to numbers
    if count in ['junction','clone']:
        for tup in counts.walk():
            (keylist,val) = (tup[:-1],tup[-1])
            counts.nested_assign(keylist,len(val))
    
    counts.lock()
    for feature in features: uniq_feature_values[feature] = list(uniq_feature_values[feature])
    
    return (uniq_feature_values,counts)

def imgt2countdict(inhandle,features,count='read'):
    return iterator2countdict(vdj.parse_imgt(inhandle),features,count)

def countdict2matrix(features,feature_values,countdict):
    # feature_values is a dict where keys are the features and the values are
    # the list of specific values I should process for that feature.
    
    dim = tuple([len(feature_values[f]) for f in features])
    matrix = np.zeros(dim,dtype=np.int)
    
    for posvals in itertools.product( *[list(enumerate(feature_values[f])) for f in features] ):
        (pos,vals) = zip(*posvals)
        count = countdict
        for val in vals:
            try:
                count = count[val]
            except KeyError:
                count = 0
                break
        matrix[pos] = count
    
    return matrix

def barcode_clone_counts(inhandle):
    """Return count dict from vdjxml file with counts[barcode][clone]"""
    counts = dict()
    for chain in vdj.parse_VDJXML(inhandle):
        try:    # chain may not have barcode
            counts_barcode = counts.setdefault(chain.barcode,dict())
        except AttributeError:
            continue
        counts_barcode[chain.clone] = counts_barcode.get(chain.clone,0) + 1
    return counts


def barcode_junction_counts(inhandle):
    """Return count dict from vdjxml file with counts[barcode][junction]"""
    counts = dict()
    for chain in vdj.parse_VDJXML(inhandle):
        try:    # chain may not have barcode
            counts_barcode = counts.setdefault(chain.barcode,dict())
        except AttributeError:
            continue
        counts_barcode[chain.junction] = counts_barcode.get(chain.junction,0) + 1
    return counts


def barcode_clone_counts2matrix(counts,barcodes=None,clones=None):
    """Generates matrix from count dict"""
    if barcodes == None:
        barcodes = counts.keys()
    if clones == None:
        clones = list( reduce( lambda x,y: x|y, [set(c.keys()) for c in counts.itervalues()] ) )
    matrix = np.zeros((len(clones),len(barcodes)))
    for (col,barcode) in enumerate(barcodes):
        for (row,clone) in enumerate(clones):
            matrix[row,col] = counts.get(barcode,dict()).get(clone,0)
    return (clones,barcodes,matrix)

barcode_junction_counts2matrix = barcode_clone_counts2matrix


# ====================================================================
# = OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD OLD  =
# ====================================================================


# # ===============
# # = Time series =
# # ===============
# 
# def clone_timeseries(inhandle, barcodes, reference_clones=None):
#     # generate count data
#     for chain in vdj.parse_VDJXML(inhandle):
#         counts
# 
# def clone_timeseries(inhandle,time_tags,reference_clones=None):
#     # get count data
#     time_tags_set = set(time_tags)
#     clone_counts = {}
#     for tag in time_tags:
#         clone_counts[tag]={}
#     for chain in vdj.parse_VDJXML(inhandle):
#         try:
#             curr_time_tag = (chain.tags&time_tags_set).pop()
#         except KeyError:
#             continue    
#         
#         try:
#             clone_counts[curr_time_tag][vdj.get_clone(chain)] += 1
#         except KeyError:
#             clone_counts[curr_time_tag][vdj.get_clone(chain)] = 1
#     
#     # set up reference clones
#     if reference_clones == None:
#         reference_clones = set()
#         for counts in clone_counts.itervalues():
#             reference_clones.update(counts.keys())
#         reference_clones = list(reference_clones)
#     
#     # build timeseries matrix
#     num_clones = len(reference_clones)
#     num_times = len(time_tags)
#     countdata = np.zeros((num_clones,num_times))
#     for (i,tag) in enumerate(time_tags):
#         countdata[:,i] = vdj.count_dict_clone_counts(clone_counts[tag],reference_clones)
#     
#     return countdata,reference_clones
# 
# 
# def timeseries2proportions(countdata,freq=True,log=True,pseudocount=1e-1):
#     num_time_series, num_times = countdata.shape
#     num_transitions = num_times - 1
#     if pseudocount != 0:
#         proportions = np.zeros((num_time_series,num_transitions))
#         countdata_pseudo = countdata + np.float_(pseudocount)
#         if freq == True:
#             countdata_pseudo = countdata_pseudo / countdata_pseudo.sum(axis=0)
#         for i in range(num_transitions):
#             proportions[:,i] = countdata_pseudo[:,i+1] / countdata_pseudo[:,i]
#     else:   # only look at time series that are non-zero the whole way through
#         idxs = np.sum(countdata>0,axis=1)==num_times
#         proportions = np.zeros((np.sum(idxs),num_transitions))
#         if freq == True:
#             countdata_modified = np.float_(countdata) / countdata.sum(axis=0)
#         else:
#             countdata_modified = np.float_(countdata)
#         for i in range(num_transitions):
#             proportions[:,i] = countdata_modified[idxs,i+1] / countdata_modified[idxs,i]
#     if log==True:
#         return np.log10(proportions)
#     else:
#         return proportions
# 
# 
# def timeseries2autocorrelation(timeseries):
#     ac = [1.] + [sp.stats.pearsonr(timeseries[:-i],timeseries[i:])[0] for i in range(1,len(timeseries)-2)]
#     return ac
# 
# # ================
# # = Spectratypes =
# # ================
# 
# def cdr3s2spectratype(cdr3s):
#     min_raw_cdr3 = np.min(cdr3s)
#     max_raw_cdr3 = np.max(cdr3s)
#     min_cdr3 = np.int(np.ceil( min_raw_cdr3 / 3.) * 3)  # will be a nonzero mult of 3
#     max_cdr3 = np.int(np.floor(max_raw_cdr3 / 3.) * 3)  # will be a mult of 3
#     
#     # bin the CDR3s lengths.  The first elt is rep zero len (and should be zero)
#     # and the last bin always represents one greater than the biggest mult of 3
#     binnedcdr3s = np.histogram(cdr3s,bins=np.arange(0,max_cdr3+2))[0]   # the +2 is due to the pecul. of np.hist.
#     
#     gaussians = []
#     for cdr3len in np.arange(min_cdr3,max_raw_cdr3,3):
#         totalcdr3s = np.sum(binnedcdr3s[cdr3len-1:cdr3len+2])
#         goodcdr3s  = binnedcdr3s[cdr3len]
#         if totalcdr3s == 0:
#             continue
#         mu = cdr3len
#         x = cdr3len-0.5
#         tail = (1 - (np.float(goodcdr3s)/totalcdr3s)) / 2.
#         sigma = (x-mu) / (np.sqrt(2.)*sp.special.erfinv(2*tail-1))
#         rv = sp.stats.norm(loc=mu,scale=sigma)
#         gaussians.append( (totalcdr3s,rv) )
#     
#     t = np.linspace(0,max_cdr3+1,1000)
#     y = np.zeros(len(t))
#     for (s,rv) in gaussians:
#         y += s*rv.pdf(t)
#     return (t,y)
# 
# 
# def spectratype_curves(inhandle):
#     # NOTE: requires chains with V and J alns
#     
#     # init data structure
#     cdr3s = {}
#     for v_seg in vdj.refseq.IGHV_seqs.keys():
#         for j_seg in vdj.refseq.IGHJ_seqs.keys():
#             cdr3s[vdj.vj_id(v_seg,j_seg)] = []
#     
#     # load data
#     for chain in vdj.parse_VDJXML(inhandle):
#         if chain.v == '' or chain.j == '' or chain.junction == '':
#             continue
#         cdr3s[vdj.vj_id(chain.v,chain.j)].append(chain.cdr3)
#     
#     spectras = {}
#     for v_seg in vdj.refseq.IGHV_seqs.keys():
#         for j_seg in vdj.refseq.IGHJ_seqs.keys():
#             if len(cdr3s[vdj.vj_id(v_seg,j_seg)]) == 0:
#                 # empty VJ combo:
#                 spectras[vdj.vj_id(v_seg,j_seg)] = (np.array([0,150]),np.array([0,0]))
#             else:
#                 spectras[vdj.vj_id(v_seg,j_seg)] = cdr3s2spectratype(cdr3s[vdj.vj_id(v_seg,j_seg)])
#     
#     return spectras
# 
# 
# 
# # ========================
# # = Diversity estimation =
# # ========================
# 
# 
# def estimator_chao1(counts):
#     """Bias corrected.  See EstimateS doc (Colwell)"""
#     Sobs = len(counts)
#     F1 = np.float_(np.sum(np.int_(counts)==1))
#     F2 = np.float_(np.sum(np.int_(counts)==2))
#     chao1 = Sobs + F1*(F1-1)/(2*(F2+1))
#     return chao1
# 
# 
# def estimator_chao1_variance(counts):
#     F1 = np.float_(np.sum(np.int_(counts)==1))
#     F2 = np.float_(np.sum(np.int_(counts)==2))
#     if F1 > 0 and F2 > 0:
#         chao1_var = (F1*(F1-1)/(2*(F2+1))) + (F1*(2*F1-1)*(2*F1-1)/(4*(F2+1)*(F2+1))) + (F1*F1*F2*(F1-1)*(F1-1)/(4*(F2+1)*(F2+1)*(F2+1)*(F2+1)))
#     elif F1 > 0 and F2 == 0:
#         Schao1 = estimator_chao1(counts)
#         chao1_var = (F1*(F1-1)/2) + (F1*(2*F1-1)*(2*F1-1)/4) - (F1*F1*F1*F1/(4*Schao1))
#     elif F1 == 0:
#         N = np.float_(np.sum(counts))
#         Sobs = np.float_(len(counts))
#         chao1_var = Sobs*np.exp(-1*N*Sobs) * (1-np.exp(-1*N*Sobs))
#     return chao1_var
# 
# 
# def estimator_ace(counts,rare_cutoff=10):
#     Sobs = np.float_(len(counts))
#     Srare = np.float_(np.sum(np.int_(counts)<=rare_cutoff))
#     Sabund = Sobs - Srare
#     F1 = np.float_(np.sum(np.int_(counts)==1))
#     F = lambda i: np.float_(np.sum(np.int_(counts)==i))
#     Nrare = np.float_(np.sum([i*F(i) for i in range(1,rare_cutoff+1)]))
#     #Nrare = np.float_(np.sum(counts[np.int_(counts)<=rare_cutoff]))
#     if Nrare == F1: # in accordance with EstimateS
#         return estimator_chao1(counts)
#     Cace = 1 - (F1/Nrare)
#     gamma_squared = Srare*np.sum([i*(i-1)*F(i) for i in range(1,rare_cutoff+1)])/(Cace*Nrare*(Nrare-1))
#     if gamma_squared < 0:
#         gamma_squared = 0
#     Sace = Sabund + (Srare/Cace) + (F1/Cace)*gamma_squared
#     return Sace
# 
# 
# def accumulation_curve(sample,sampling_levels):
#     pass
# 
# 
# # =========================
# # = Statistical utilities =
# # =========================
# 
# def counts2sample(counts):
#     """Computes a consistent sample from a vector of counts.
#     
#     Takes a vector of counts and returns a vector of indices x
#        such that len(x) = sum(c) and each elt of x is the index of
#        a corresponding elt in c
#     """
#     x = np.ones(np.sum(counts),dtype=np.int_)
#     
#     start_idx = 0
#     end_idx = 0
#     for i in xrange(len(counts)):
#         start_idx = end_idx
#         end_idx = end_idx + counts[i]
#         x[start_idx:end_idx] = x[start_idx:end_idx] * i 
#     return x
# 
# 
# def sample2counts(sample):
#     """Return count vector from list of samples.
#     
#     The ordering etc is ignored; only the uniqueness
#     of the objects is considered.
#     """
#     num_categories = len(set(sample))
#     count_dict = {}
#     for elt in sample:
#         try: count_dict[elt] += 1
#         except KeyError: count_dict[elt] = 1
#     return count_dict.values()
# 
# 
# # def sample2counts(sample, categories=0):
# #     """Return count vector from list of samples.
# #     
# #     Take vector of samples and return a vector of counts.  The elts
# #        refer to indices in something that would ultimately map to the
# #        originating category (like from a multinomial).  Therefore, if there
# #        are, say, 8 categories, then valid values in sample should be 0-7.
# #        If categories is not given, then i compute it from the highest value
# #        present in sample (+1).
# #     """
# #     counts = np.bincount(sample)
# #     if (categories > 0) and (categories > len(counts)):
# #         counts = np.append( counts, np.zeros(categories-len(counts)) )
# #     return counts
# 
# 
# def scoreatpercentile(values,rank):
#     return sp.stats.scoreatpercentile(values,rank)
# 
# 
# def percentileofscore(values,score):
#     values.sort()
#     return values.searchsorted(score) / np.float_(len(values))
# 
# 
# def bootstrap(x, nboot, theta, *args):
#     '''return n bootstrap replications of theta from x'''
#     N = len(x)
#     th_star = np.zeros(nboot)
#     
#     for i in xrange(nboot):
#         th_star[i] = theta( x[ np.random.randint(0,N,N) ], *args )    # bootstrap repl from x
#     
#     return th_star
# 
# 
# def subsample(x, num_samples, sample_size, theta, *args):
#     """return num_samples evaluations of the statistic theta
#     on subsamples of size sample_size"""
#     N = len(x)
#     th_star = np.zeros(num_samples)
#     
#     for i in xrange(num_samples):
#         th_star[i] = theta( x[ np.random.randint(0,N,sample_size) ], *args )    # subsample from from x
#     
#     return th_star
# 
# 
# def randint_without_replacement(low,high=None,size=None):
#     if high == None:
#         high = low
#         low = 0
#     if size == None:
#         size = 1
#     urn = range(low,high)
#     N = len(urn)
#     flip = False
#     if size > N/2:
#         flip = True
#         size = N - size
#     sample = []
#     for i in xrange(size):
#         draw = np.random.randint(0,N-i)
#         sample.append(urn.pop(draw))
#     if not flip:
#         return np.asarray(sample)
#     else:
#         return np.asarray(urn)
# 
# 
# def subsample_without_replacement(x, num_samples, sample_size, theta, *args):
#     """return num_samples evaluations of the statistic theta
#     on subsamples of size sample_size"""
#     N = len(x)
#     th_star = np.zeros(num_samples)
#     
#     for i in xrange(num_samples):
#         th_star[i] = theta( x[ randint_without_replacement(0,N,sample_size) ], *args )    # subsample from from x
#     
#     return th_star
# 
# 
# 
# # =================
# # = Visualization =
# # =================
# 
# class ConstWidthRectangle(mpl.patches.Patch):
#     
#     def __init__(self, x, y1, y2, w, **kwargs):
#         self.x  = x
#         self.y1 = y1
#         self.y2 = y2
#         self.w  = w
#         mpl.patches.Patch.__init__(self,**kwargs)
#     
#     def get_path(self):
#         return mpl.path.Path.unit_rectangle()
#     
#     def get_transform(self):
#         box = np.array([[self.x,self.y1],
#                         [self.x,self.y2]])
#         box = self.axes.transData.transform(box)
#         
#         w = self.w * self.axes.bbox.width / 2.0
#         
#         box[0,0] -= w
#         box[1,0] += w
#         
#         return mpl.transforms.BboxTransformTo(mpl.transforms.Bbox(box))
# 
# class ConstWidthLine(mpl.lines.Line2D):
#     
#     def __init__(self,x,y,w,**kwargs):
#         self.x = x
#         self.y = y
#         self.w = w
#         mpl.lines.Line2D.__init__(self,[0,1],[0,0],**kwargs) # init to unit line
#     
#     def get_transform(self):
#         # define transform that takes unit horiz line seg
#         # and places it in correct position using display
#         # coords
#         
#         box = np.array([[self.x,self.y],
#                         [self.x,self.y+1]])
#         box = self.axes.transData.transform(box)
#         
#         w = self.w * self.axes.bbox.width / 2.0
#         
#         box[0,0] -= w
#         box[1,0] += w
#         
#         #xdisp,ydisp = self.axes.transData.transform_point([self.x,self.y])
#         #xdisp -= w
#         #xleft  = xdisp - w
#         #xright = xdisp + w
#         
#         return mpl.transforms.BboxTransformTo(mpl.transforms.Bbox(box))
#         #return mpl.transforms.Affine2D().scale(w,1).translate(xdisp,ydisp)
#     
#     def draw(self,renderer):
#         # the ONLY purpose of redefining this function is to force the Line2D
#         # object to execute recache().  Otherwise, certain changes in the scale
#         # do not invalidate the Line2D object, and the transform will not be
#         # recomputed (and so the Axes coords computed earlier will be obsolete)
#         self.recache()
#         return mpl.lines.Line2D.draw(self,renderer)
# 
# 
# class ConstHeightRectangle(mpl.patches.Patch):
#     
#     def __init__(self, x1, x2, y, h, **kwargs):
#         self.x1 = x1
#         self.x2 = x2
#         self.y  = y
#         self.h  = h
#         mpl.patches.Patch.__init__(self,**kwargs)
#     
#     def get_path(self):
#         return mpl.path.Path.unit_rectangle()
#     
#     def get_transform(self):
#         box = np.array([[self.x1,self.y],
#                         [self.x2,self.y]])
#         box = self.axes.transData.transform(box)
#         
#         h = self.h * self.axes.bbox.height / 2.0
#         
#         box[0,1] -= h
#         box[1,1] += h
#         
#         return mpl.transforms.BboxTransformTo(mpl.transforms.Bbox(box))
# 
# class ConstHeightLine(mpl.lines.Line2D):
#     
#     def __init__(self,x,y,h,**kwargs):
#         self.x = x
#         self.y = y
#         self.h = h
#         mpl.lines.Line2D.__init__(self,[0,0],[0,1],**kwargs) # init to unit line
#         
#         # self.x = x
#         # self.y = y
#         # self.w = w
#         # mpl.lines.Line2D.__init__(self,[0,1],[0,0],**kwargs) # init to unit line
#     
#     def get_transform(self):
#         # define transform that takes unit horiz line seg
#         # and places it in correct position using display
#         # coords
#         
#         box = np.array([[self.x,self.y],
#                         [self.x+1,self.y]])
#         box = self.axes.transData.transform(box)
#         
#         h = self.h * self.axes.bbox.height / 2.0
#         
#         box[0,1] -= h
#         box[1,1] += h
#         
#         #xdisp,ydisp = self.axes.transData.transform_point([self.x,self.y])
#         #xdisp -= w
#         #xleft  = xdisp - w
#         #xright = xdisp + w
#         
#         return mpl.transforms.BboxTransformTo(mpl.transforms.Bbox(box))
#         #return mpl.transforms.Affine2D().scale(w,1).translate(xdisp,ydisp)
#     
#     def draw(self,renderer):
#         # the ONLY purpose of redefining this function is to force the Line2D
#         # object to execute recache().  Otherwise, certain changes in the scale
#         # do not invalidate the Line2D object, and the transform will not be
#         # recomputed (and so the Axes coords computed earlier will be obsolete)
#         self.recache()
#         return mpl.lines.Line2D.draw(self,renderer)
# 
# 
# def boxplot(ax, x, positions=None, widths=None, vert=1):
#     # adapted from matplotlib
#     
#     # convert x to a list of vectors
#     if hasattr(x, 'shape'):
#         if len(x.shape) == 1:
#             if hasattr(x[0], 'shape'):
#                 x = list(x)
#             else:
#                 x = [x,]
#         elif len(x.shape) == 2:
#             nr, nc = x.shape
#             if nr == 1:
#                 x = [x]
#             elif nc == 1:
#                 x = [x.ravel()]
#             else:
#                 x = [x[:,i] for i in xrange(nc)]
#         else:
#             raise ValueError, "input x can have no more than 2 dimensions"
#     if not hasattr(x[0], '__len__'):
#         x = [x]
#     col = len(x)
#     
#     # get some plot info
#     if positions is None:
#         positions = range(1, col + 1)
#     if widths is None:
#         widths = min(0.3/len(positions),0.05)
#     if isinstance(widths, float) or isinstance(widths, int):
#         widths = np.ones((col,), float) * widths
#     
#     # loop through columns, adding each to plot
#     for i,pos in enumerate(positions):
#         d = np.ravel(x[i])
#         row = len(d)
#         if row==0:
#             # no data, skip this position
#             continue
#         # get distrib info
#         q1, med, q3 = mpl.mlab.prctile(d,[25,50,75])
#         dmax = np.max(d)
#         dmin = np.min(d)
#         
#         line_color = '#074687'
#         face_color = '#96B7EC'
#         if vert == 1:
#             medline = ConstWidthLine(pos,med,widths[i],color=line_color,zorder=3)
#             box = ConstWidthRectangle(pos,q1,q3,widths[i],facecolor=face_color,edgecolor=line_color,zorder=2)
#             vertline = mpl.lines.Line2D([pos,pos],[dmin,dmax],color=line_color,zorder=1)
#         else:
#             medline = ConstHeightLine(med,pos,widths[i],color=line_color,zorder=3)
#             box = ConstHeightRectangle(q1,q3,pos,widths[i],facecolor=face_color,edgecolor=line_color,zorder=2)
#             vertline = mpl.lines.Line2D([dmin,dmax],[pos,pos],color=line_color,zorder=1)
#         
#         ax.add_line(vertline)
#         ax.add_patch(box)
#         ax.add_line(medline)
# 
# 
# 
# 
# 
# 
# #==============================================================================
# #==============================================================================
# #==============================================================================
# 
# def rep2spectratype(rep):
#     """Compute spectratype curves from Repertoire object."""
#     
#     cdr3s = np.array([c.cdr3 for c in rep if c.junction != ''])
#     min_raw_cdr3 = np.min(cdr3s)
#     max_raw_cdr3 = np.max(cdr3s)
#     min_cdr3 = np.int(np.ceil( min_raw_cdr3 / 3.) * 3)  # will be a nonzero mult of 3
#     max_cdr3 = np.int(np.floor(max_raw_cdr3 / 3.) * 3)  # will be a mult of 3
#     
#     # bin the CDR3s lengths.  The first elt is rep zero len (and should be zero)
#     # and the last bin always represents one greater than the biggest mult of 3
#     binnedcdr3s = np.histogram(cdr3s,bins=np.arange(0,max_cdr3+2))[0]   # the +2 is due to the pecul. of np.hist.
#     
#     gaussians = []
#     for cdr3len in np.arange(min_cdr3,max_raw_cdr3,3):
#         totalcdr3s = np.sum(binnedcdr3s[cdr3len-1:cdr3len+2])
#         goodcdr3s  = binnedcdr3s[cdr3len]
#         if totalcdr3s == 0:
#             continue
#         mu = cdr3len
#         x = cdr3len-0.5
#         tail = (1 - (np.float(goodcdr3s)/totalcdr3s)) / 2.
#         sigma = (x-mu) / (np.sqrt(2.)*sp.special.erfinv(2*tail-1))
#         rv = sp.stats.norm(loc=mu,scale=sigma)
#         gaussians.append( (totalcdr3s,rv) )
#     
#     t = np.linspace(0,max_cdr3+1,1000)
#     y = np.zeros(len(t))
#     for (s,rv) in gaussians:
#         y += s*rv.pdf(t)
#     return (t,y)
# 
# 
# 
# 
# def scatter_repertoires_ontology(reps,info='VJCDR3',gooddata=False,measurement='proportions'):
#     """Create a grid of scatter plots showing corelations between all pairs of repertoires.
#     
#     reps -- list of Repertoire objects
#     
#     """
#     numreps = len(reps)
#     
#     datalist = []
#     for rep in reps:
#         datalist.append( vdj.counts_ontology_1D(rep,info,gooddata) )
#     
#     if measurement == 'proportions':
#         for i in xrange(len(datalist)):
#             datalist[i] = np.float_(datalist[i]) / np.sum(datalist[i])
#     
#     min_nonzero = np.min([np.min(data[data>0]) for data in datalist])
#     max_nonzero = np.max([np.max(data[data>0]) for data in datalist])
#     axislo = 10**np.floor( np.frexp(min_nonzero)[1] * np.log10(2) )
#     axishi = 10**np.ceil(  np.frexp(max_nonzero)[1] * np.log10(2) )
#     
#     fig = plt.figure()
#     
#     hist_axs = []
#     for row in xrange(numreps):
#         col = row
#         plotnum = numreps*row + col + 1
#         ax = fig.add_subplot(numreps,numreps,plotnum)
#         ax.hist(datalist[row],bins=100,log=True,facecolor='k')
#         hist_axs.append(ax)
#     
#     scatter_axs = []
#     for row in xrange(numreps-1):
#         for col in xrange(row+1,numreps):
#             plotnum = numreps*row + col + 1
#             ax = fig.add_subplot(numreps,numreps,plotnum)
#             ax.scatter(datalist[row],datalist[col],c='k',marker='o',s=2,edgecolors=None)
#             ax.set_xscale('log')
#             ax.set_yscale('log')
#             ax.axis([axislo,axishi,axislo,axishi])
#             scatter_axs.append(ax)
#     
#     return fig
# 
# def scatter_repertoires_clusters(reps,refclusters,measurement='proportions'):
#     """Create a grid of scatter plots showing corelations between all pairs of repertoires.
#     
#     reps -- list of Repertoire objects
#     
#     """
#     numreps = len(reps)
#     
#     datalist = []
#     for rep in reps:
#         clusters = vdj.getClusters(rep)
#         datalist.append( vdj.countsClusters(clusters,refclusters) )
#     
#     if measurement == 'proportions':
#         for i in xrange(len(datalist)):
#             datalist[i] = np.float_(datalist[i]) / np.sum(datalist[i])
#     
#     min_nonzero = np.min([np.min(data[data>0]) for data in datalist])
#     max_nonzero = np.max([np.max(data[data>0]) for data in datalist])
#     axislo = 10**np.floor( np.frexp(min_nonzero)[1] * np.log10(2) )
#     axishi = 10**np.ceil(  np.frexp(max_nonzero)[1] * np.log10(2) )
#     
#     fig = plt.figure()
#     
#     hist_axs = []
#     for row in xrange(numreps):
#         col = row
#         plotnum = numreps*row + col + 1
#         ax = fig.add_subplot(numreps,numreps,plotnum)
#         ax.hist(datalist[row],bins=100,log=True,facecolor='k')
#         hist_axs.append(ax)
#     
#     scatter_axs = []
#     for row in xrange(numreps-1):
#         for col in xrange(row+1,numreps):
#             plotnum = numreps*row + col + 1
#             ax = fig.add_subplot(numreps,numreps,plotnum)
#             ax.scatter(datalist[row],datalist[col],c='k',marker='o',s=0.5,edgecolors=None)
#             ax.set_xscale('log')
#             ax.set_yscale('log')
#             ax.axis([axislo,axishi,axislo,axishi])
#             scatter_axs.append(ax)
#     
#     return fig
# 
# def reps2timeseries(reps,refclusters):
#     """Return time series matrix from list or Repertoire objs in chron order.
#     
#     reps is list of Repertoire objects
#     refclusters is the master list of reference clusters
#     """
#     numreps = len(reps)
#     numclusters = len(refclusters)
#     
#     countdata = np.zeros((numclusters,numreps))
#     for (i,rep) in enumerate(reps):
#         clusters = vdj.getClusters(rep)
#         countdata[:,i] = vdj.countsClusters(clusters,refclusters)
#     
#     return countdata
# 
# def timeseries_repertoires(times,reps,refclusters,idxsbool=None,allpositive=False):
#     """Create a time-series of the different clusters in refclusters.
#     
#     If allpositive is True, then it will limit itself to drawing timeseries
#     only for those clusters that are non-zero at all timepoints.
#     
#     """
#     ax = plt.gca()
#     
#     numreps = len(reps)
#     numclusters = len(refclusters)
#     
#     countdata = np.zeros((numclusters,numreps))
#     for (i,rep) in enumerate(reps):
#         clusters = vdj.getClusters(rep)
#         countdata[:,i] = vdj.countsClusters(clusters,refclusters)
#     
#     sums = countdata.sum(0)
#     proportions = np.float_(countdata) / sums
#     
#     if idxsbool == None:
#         if allpositive == True:
#             idxsbool = np.sum(proportions,axis=1) > 0
#         else:
#             idxsbool = np.array([True]*proportions.shape[0])
#     
#     ax.plot(times,countdata[idxsbool,:].transpose(),'k-',linewidth=0.2)
#     
#     plt.draw_if_interactive()
#     return ax
# 
# 
# def rep2spectratype(rep):
#     """Compute spectratype curves from Repertoire object."""
#     
#     cdr3s = np.array([c.cdr3 for c in rep if c.junction != ''])
#     min_raw_cdr3 = np.min(cdr3s)
#     max_raw_cdr3 = np.max(cdr3s)
#     min_cdr3 = np.int(np.ceil( min_raw_cdr3 / 3.) * 3)  # will be a nonzero mult of 3
#     max_cdr3 = np.int(np.floor(max_raw_cdr3 / 3.) * 3)  # will be a mult of 3
#     
#     # bin the CDR3s lengths.  The first elt is rep zero len (and should be zero)
#     # and the last bin always represents one greater than the biggest mult of 3
#     binnedcdr3s = np.histogram(cdr3s,bins=np.arange(0,max_cdr3+2))[0]   # the +2 is due to the pecul. of np.hist.
#     
#     gaussians = []
#     for cdr3len in np.arange(min_cdr3,max_raw_cdr3,3):
#         totalcdr3s = np.sum(binnedcdr3s[cdr3len-1:cdr3len+2])
#         goodcdr3s  = binnedcdr3s[cdr3len]
#         if totalcdr3s == 0:
#             continue
#         mu = cdr3len
#         x = cdr3len-0.5
#         tail = (1 - (np.float(goodcdr3s)/totalcdr3s)) / 2.
#         sigma = (x-mu) / (np.sqrt(2.)*sp.special.erfinv(2*tail-1))
#         rv = sp.stats.norm(loc=mu,scale=sigma)
#         gaussians.append( (totalcdr3s,rv) )
#     
#     t = np.linspace(0,max_cdr3+1,1000)
#     y = np.zeros(len(t))
#     for (s,rv) in gaussians:
#         y += s*rv.pdf(t)
#     return (t,y)
# 
# def circlemapVJ(ax,counts,rowlabels=None,collabels=None,scale='linear'):
#     numV = counts.shape[0]
#     numJ = counts.shape[1]
#     X,Y = np.meshgrid(range(numJ),range(numV))
#     
#     # mask zero positions
#     X,Y = ma.array(X), ma.array(Y)
#     C = ma.array(counts)
#     
#     zeromask = (counts == 0)
#     X.mask = zeromask
#     Y.mask = zeromask
#     C.mask = zeromask
#     
#     # ravel nonzero elts (deletes zero-positions)
#     x = ma.compressed(X)
#     y = ma.compressed(Y)
#     c = ma.compressed(C)
#     
#     # log normalize counts if requested
#     if scale == 'log':
#         c = ma.log10(c)
#     
#     if scale == 'linear' or scale == 'log':
#         # normalize counts to desired size-range
#         max_counts = ma.max(c)
#         min_counts = ma.min(c)
#         counts_range = max_counts - min_counts
#     
#         max_size = 100
#         min_size = 5
#         size_range = max_size - min_size
#     
#         sizes = (np.float(size_range) / counts_range) * (c - min_counts) + min_size
#     
#     if scale == 'custom':
#         trans_counts = 1000
#         linear_positions = c >= trans_counts
#         log_positions = c < trans_counts
#         
#         min_size = 3
#         trans_size = 40 # 30
#         max_size = 200 # 150
#         log_size_range = trans_size - min_size
#         linear_size_range = max_size - trans_size
#         
#         linear_max_counts = ma.max(c[linear_positions])
#         linear_min_counts = ma.min(c[linear_positions])
#         linear_counts_range = linear_max_counts - linear_min_counts
#         log_max_counts = ma.max(c[log_positions])
#         log_min_counts = ma.min(c[log_positions])
#         log_counts_range = np.log10(log_max_counts) - np.log10(log_min_counts)
#         
#         sizes = np.zeros(len(c))
#         sizes[linear_positions] = (np.float(linear_size_range) / linear_counts_range) * (c[linear_positions] - linear_min_counts) + trans_size
#         sizes[log_positions] = (np.float(log_size_range) / log_counts_range) * (ma.log10(c[log_positions]) - ma.log10(log_min_counts)) + min_size
#     
#     collection = mpl.collections.CircleCollection(
#                                         sizes,
#                                         offsets = zip(x,y),
#                                         transOffset = ax.transData, # i may need to explicitly set the xlim and ylim info for this to work correctly
#                                         facecolors = '#1873C1',
#                                         linewidths = 0.25,
#                                         clip_on = False)
#     
#     ax.add_collection(collection)
#     
#     ax.set_aspect('equal')
#     ax.autoscale_view()
#     
#     ax.xaxis.set_major_locator(mpl.ticker.FixedLocator(range(counts.shape[1])))
#     ax.yaxis.set_major_locator(mpl.ticker.FixedLocator(range(counts.shape[0])))
#     
#     if rowlabels != None:
#         ax.xaxis.set_major_formatter(mpl.ticker.FixedFormatter(collabels))
#     if collabels != None:
#         ax.yaxis.set_major_formatter(mpl.ticker.FixedFormatter(rowlabels))
#     
#     for ticklabel in ax.xaxis.get_ticklabels():
#         ticklabel.set_horizontalalignment('left')
#         ticklabel.set_rotation(-45)
#         ticklabel.set_size(8)
#     
#     for ticklabel in ax.yaxis.get_ticklabels():
#         ticklabel.set_size(8)
#     
#     if scale == 'linear' or scale == 'log':
#         return (min_counts,max_counts),(min_size,max_size)
#     else:
#         return (linear_min_counts,trans_counts,log_max_counts),(min_size,trans_size,max_size)
# 
# # define colormap for -1 to 1 (green-black-red) like gene expression
# redgreencdict = {'red': [(0.0,   0.0,   0.0),
#                          (0.5,   0.0,   0.0),
#                          (1.0,   1.0,   0.0)],
#                         
#                 'green':[(0.0,   0.0,   1.0),
#                          (0.5,   0.0,   0.0),
#                          (1.0,   0.0,   0.0)],
#                         
#                 'blue': [(0.0,   0.0,   0.0),
#                          (0.5,   0.0,   0.0),
#                          (1.0,   0.0,   0.0)]}
# 
# redgreen = mpl.colors.LinearSegmentedColormap('redgreen',redgreencdict,256)
# redgreen.set_bad(color='w')
