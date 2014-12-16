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

import numpy as np
import scipy as sp
import scipy.cluster

import vdj
import clusteringcore


def pdist(X,metric):
    m = len(X)
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)
    k = 0
    for i in xrange(0, m - 1):
        for j in xrange(i+1, m):
            dm[k] = metric(X[i], X[j])
            k += 1
    return dm


def cluster_seqs(seqs,cutoff=4.5,linkage='single'):
    # check trivial cases
    if len(seqs) == 0:
        return (np.array([]),{})
        # raise Exception, "chains has nothing it"
    
    # collapse identical seqs into each other
    unique_seqs = list(set(seqs))
    seq_idxs = dict( [(j,i) for (i,j) in enumerate(unique_seqs)] )
    
    # check trivial case
    if len(unique_seqs) == 1:
        T = np.array([1]*len(seqs))
        return (T,seq_idxs)
    
    # compute the distance matrix
    Y = pdist( unique_seqs, clusteringcore.levenshtein )
    
    # compute the linkage
    Z = sp.cluster.hierarchy.linkage(Y,method=linkage)
    
    # determine the clusters at level cutoff
    T = sp.cluster.hierarchy.fcluster(Z,cutoff,criterion='distance')
    
    return (T,seq_idxs)
