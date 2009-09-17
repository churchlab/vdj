import numpy as np

def pdist(X,metric):
    m = len(X)
    dm = np.zeros((m * (m - 1) / 2,), dtype=np.double)
    k = 0
    for i in xrange(0, m - 1):
        for j in xrange(i+1, m):
            dm[k] = metric(X[i], X[j])
            k += 1
    return dm

def clusterChains(chains,cutoff=4.5,tag_chains=False,tag=''):
    """Cluster list of ImmuneChains using Levenshtein edit distance between junctions.
    
    Returns a vector T where T[i] is the cluster number for chains[i]
    
    Note: to improve the performance, this method first collapses all identical chains
          and records the weight.  It then computes the initial distance matrix taking
          the weights into account.  The rest of the UPGMA method will propagate this
          information.
    
    """
    # check trivial cases
    if len(chains) == 0:
        raise Exception, "chains has nothing it"
    
    # collapse identical junctions into each other
    unique_junctions = list( set( [c.junction for c in chains] ) )
    junction_idx = dict( [(j,i) for i,j in enumerate(unique_junctions)] )
    
    if len(unique_junctions) == 1:
        T = np.array([1]*len(chains))
        if tag_chains == True:
            for chain in chains:
                chain.add_tags('cluster|'+tag+'|'+str(1))
        return T
    
    # compute the distance matrix
    Y = pdist( unique_junctions, clusteringcore.levenshtein )
    
    # compute the linkage
    Z = sp.cluster.hierarchy.linkage(Y,method='single')
    
    # determine the clusters at level cutoff
    T = sp.cluster.hierarchy.fcluster(Z,cutoff,criterion='distance')
    
    # perform chain tagging
    if tag_chains == True:
        for (i,chain) in enumerate(chains):
            chain.add_tags('cluster|'+tag+'|'+str(T[junction_idx[chain.junction]]))
    
    return T

def clusterRepertoire(rep,cutoff=4.5,tag_rep=False,tag_chains=False,tag=''):
    """Cluster the chains in a Repertoire object.
    
    First the algorithm partitions the chains according to V-J combo.
    Then clusters based on junction sequences for each partition
    separately.
    
    Only clusters the set of chains that has both a V and J and CDR3.
    
    If tag_chains is True, it will add a cluster tag to each chain.
    tag will be incorporated as well.
    
    If tag_chains is False, the fn will return a list of lists,
    each one of which represents a cluster, and which is composed
    of the descriptions of the corresponding chains in the cluster.
    
    """
    if tag == '':
        reptag = ''
    else:
        reptag = tag+'|'
    
    repgood = rep.get_chains_fullVJCDR3()
    clusters = []
    for vseg in refseq.IGHV_seqs.keys():
        for jseg in refseq.IGHJ_seqs.keys():
            currtag = reptag+vseg+'|'+jseg
            currchains = repgood.get_chains_AND([vseg,jseg]).chains
            if len(currchains) == 0:
                continue
            T = clusterChains(repgood.get_chains_AND([vseg,jseg]).chains,cutoff,tag_chains,currtag)
            numclusters = len(set(T))
            currclusters = [ [] for i in np.arange(numclusters)]
            for (i,clust) in enumerate(T):
                currclusters[clust-1].append(currchains[i].descr)
            clusters.extend(currclusters)
    if tag_rep == True:
        rep.add_metatags("Clustering|" + tag + "levenshtein|single_linkage|cutoff="+str(cutoff)+"|"+timestamp())
    return clusters

def getClusters(rep):
    clusters = {}
    for (tag,idxs) in rep.tags.iteritems():
        if tag.startswith('cluster'):
            # error checking: make sure every new cluster is unique
            if clusters.has_key(tag):
                raise Exception, "repertoire object's tags has multiple copies of the same tag"
            clusters[tag]=idxs
    return clusters