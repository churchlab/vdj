import types
import string

import numpy as np
import Bio.SeqIO

import refseq
import sequtils
import alignment

# import xml.sax
# import xml.sax.handler
# import time
# import datetime
# import operator
# import os
# 
# import numpy as np
# import scipy as sp
# import scipy.cluster
# 
# import seqtools
# import refseq
# import alignmentcore
# import clusteringcore



# ===================
# = DATA STRUCTURES =
# ===================

class ImmuneChain(object):
    """Data structure to represent an immune chain."""
    
    def __init__(self, seq='', func='', v='', d='', j='', ighc='', junction='', descr='', tags=set([])):
        """Initialize ImmuneChain
        
        seq is 5'->3'
        """
        self.seq = seq.upper()
        self.descr = descr
        if isinstance(tags,types.StringTypes): tags = [tags]
        self.tags = set(tags)   # tag for sample number/experiment etc
        self.v = v
        self.d = d
        self.j = j
        self.ighc = ighc
        self.junction = junction
        self.func = func
    
    @property
    def cdr3(self):
        return len(self.junction)
    
    @property
    def all_tags(self):
        """Return set object containing all non-empty tags and identifiers.
        
        This includes all tags, v, d, j, ighc, and descr
        """
        return (self.tags | set([self.v,self.d,self.j,self.ighc,self.descr])).discard('')
    
    def add_tags(self,tagset):
        if isinstance(tagset,types.StringTypes): tagset = [tagset]
        self.tags.update(tagset)
    
    def remove_tags(self,tagset):
        if isinstance(tagset,types.StringTypes): tagset = [tagset]
        for tag in tagset: self.tags.remove(tag)
    
    def __len__(self):
        return len(self.seq)
    
    def __str__(self):
        return self.__repr__()
    
    def __repr__(self):
        return self.get_XML()
    
    def get_XML(self):
        xmlstring = ''
        xmlstring += '<ImmuneChain>\n'
        xmlstring += '\t<descr>'    + self.descr +      '</descr>\n'
        xmlstring += '\t<seq>'      + self.seq +        '</seq>\n'
        xmlstring += '\t<v>'        + self.v +          '</v>\n' 
        xmlstring += '\t<d>'        + self.d +          '</d>\n'
        xmlstring += '\t<j>'        + self.j +          '</j>\n'
        xmlstring += '\t<ighc>'     + self.ighc +       '</ighc>\n'
        xmlstring += '\t<cdr3>'     + str(self.cdr3) +  '</cdr3>\n' # measured in nt
        xmlstring += '\t<junction>' + self.junction +   '</junction>\n'
        xmlstring += '\t<func>'     + self.func +       '</func>\n'
        for tag in self.tags:
            xmlstring += '\t<tag>' + tag + '</tag>\n'
        xmlstring += '</ImmuneChain>\n'
        return xmlstring



# ================
# = INPUT/OUTPUT =
# ================

def parse_VDJXML(inputfile):
    """Load a data from a VDJXML file as a Repertoire or list of ImmuneChains
    
    NOTE: this fn does NOT utilize the XML libraries; it implements a manual parser
    that takes input line by line.
    
    THIS ASSUMES THAT EVERY XML ELEMENT TAKES ONE AND ONLY ONE LINE
    
    """
    
    if isinstance(inputfile,types.StringTypes):
        ip = open(inputfile,'r')
    elif isinstance(inputfile,file):
        ip = inputfile
    
    numChains = 0
    
    possible_elements = [
                'descr',
                'seq',
                'v',
                'd',
                'j',
                'ighc',
                'cdr3',
                'junction',
                'func',
                'tag'
                ]
    
    for line in ip:
        line = line.strip()
        endelementpos = line.find('>') + 1
        xmlelement = line[0:endelementpos]
        element = xmlelement[1:-1]
        
        if xmlelement == '<ImmuneChain>':
            chain = ImmuneChain()
        elif xmlelement == '</ImmuneChain>':
            numChains += 1
            yield chain
        elif element in possible_elements:
            if element == 'cdr3':
                chain.cdr3 = eval(line[endelementpos:-1*(endelementpos+1)])
            elif element == 'tag':
                chain.add_tags(line[endelementpos:-1*(endelementpos+1)])
            else:
                chain.__setattr__(element,line[endelementpos:-1*(endelementpos+1)])
    
    if isinstance(inputfile,types.StringTypes):
        ip.close()

def write_VDJXML(data, outputfile):
    """Write list of ImmuneChains to a file in VDJXML format"""
    
    if isinstance(outputfile,types.StringTypes):
        op = open(outputfile,'w')
    elif isinstance(outputfile,file):
        op = outputfile
    
    if isinstance(data,list) and isinstance(data[0],ImmuneChain):
        for (i,chain) in enumerate(data):
            print >>op, chain
    else:
        raise Exception, "Must supply a list of ImmuneChain objects."
    
    if isinstance(outputfile,types.StringTypes):
        op.close()



# ============
# = Counting =
# ============

def counts_VJ(inputfile):
    if isinstance(inputfile,types.StringTypes):
        ip = open(inputfile,'r')
    elif isinstance(inputfile,file):
        ip = inputfile
    
    counts = np.zeros( (len(refseq.IGHV_list),len(refseq.IGHJ_list)) )
    for chain in parseVDJXML(ip):
        counts[refseq.IGHV_idx[chain.v],refseq.IGHJ_idx[chain.j]] += 1
    
    if isinstance(inputfile,types.StringTypes):
        ip.close()
    
    return cn


def counts_VDJ(rep):
    cn = np.zeros( (len(refseq.IGHV_list),len(refseq.IGHD_list),len(refseq.IGHJ_list)) )
    for chain in rep.chains:
        cn[refseq.IGHV_idx[chain.v],refseq.IGHD_idx[chain.d],refseq.IGHJ_idx[chain.j]] += 1
    return cn


def reshape_counts_VDJ_2D(counts):
    return counts.reshape(len(refseq.IGHV_list),len(refseq.IGHD_list)*len(refseq.IGHJ_list))


def counts_clones(clone_idxs,reference_clones=None):
    """Takes a dictionary of cluster names mapped to a sequence of indices into an ImmuneChain list.
    
    Returns an np array of the same length as reference_clones with the counts of each
    cluster in reference_clones.
    
    The need for reference_clones is due to the fact that splitting a given repertoire
    may result in some parts not observing any of a given clone, so there needs to be a common way
    to compare two clone sets.
    
    If reference_clones is left out, then the set of clones present in clone_idxs is used.
    
    """
    if reference_clones == None:
        reference_clones = clone_idxs.keys()
    counts = np.zeros(len(reference_clones))
    for (i,name) in enumerate(reference_clones):
        counts[i] = len(clone_idxs.get(name,[]))
    return counts



# =================================
# = Retrieving tags and filtering =
# =================================

def get_tag_with_prefix(chain,prefix):
    for tag in chain.tags:
        if tag.startswith(prefix):
            return tag
    raise ValueError, "Tag that starts with " + prefix + " not found."


def get_clone(chain):
    return get_tag_with_prefix(chain,'clone')


def get_barcode(chain):
    try:
        return get_tag_with_prefix(chain,'barcode')
    except ValueError:
        return ''


def filter_tags_and(tags,inhandle,outhandle):
    if isinstance(tags,types.StringTypes): tags = [tags]
    tags = set(tags)
    for chain in parse_VDJXML(inhandle):
        if tags <= chain.all_tags:    # test that everything in tags is in all_tags
            print >>outhandle, chain


def filter_tags_or(tags,inhandle,outhandle):
    if isinstance(tags,types.StringTypes): tags = [tags]
    tags = set(tags)
    empty_set = set()
    for chain in parse_VDJXML(inhandle):
        if tags & chain.all_tags != empty_set:    # test that tags and all_tags share something
            print >>outhandle, chain


def is_full_VJ(chain):
    if (chain.v in refseq.IGHV_seqs.keys()) and (chain.j in refseq.IGHJ_seqs.keys()):
        return True
    else:
        return False



# ======================
# = Pipeline functions =
# ======================

def fasta2vdjxml(inhandle,outhandle):
    multiple_fields = False
    
    for record in Bio.SeqIO.parse(inhandle,'fasta'):
        description = record.description.split()
        sequence = record.seq.tostring()   # SeqRecord object
        if not multiple_fields and len(description) > 1:
            multiple_fields = True
        chain = ImmuneChain(descr=description[0],seq=sequence)
        print >>outhandle, chain
    
    if multiple_fields == True:
        print "WARNING: input fasta file has descriptions with multiple fields"


def size_select(min_=None,max_=None,inhandle,outhandle):
    if min_ == None:
        min_ = 0
    if max_ == None:
        max_ = float('inf')
    for chain in parse_VDJXML(inhandle):
        if len(chain) >= min_ and len(chain) <= max_:
            print >>outhandle, chain


def barcode_id(barcode_fasta,inhandle,outhandle):
    # NOTE: all barcodes must be the same length
    # NOTE: all barcode names must start with 'barcode'
    if isinstance(barcode_fasta,types.StringTypes):
        bcip = open(barcode_fasta,'r')
    elif isinstance(barcode_fasta,file):
        bcip = barcode_fasta
    
    barcodes = {}
    for record in Bio.SeqIO.parse(bcip,'fasta'):
        barcodes[record.seq.tostring()] = record.id
    
    barcode_len = len(barcodes.keys()[0])
    for bc in barcodes.keys():
        if len(bc) != barcode_len:
            raise Exception, "ERROR: All barcode lengths must be equal."
        if not barcodes[bc].startswith('barcode'):
            raise Exception, "ERROR: All barcode names must start with 'barcode'"
    
    for chain in parse_VDJXML(inhandle):
        curr_barcode = barcodes.get(chain.seq[:barcode_len],'')
        if curr_barcode != '':
            chain.seq = chain.seq[barcode_len:] # remove barcode from seq
            chain.add_tags(curr_barcode)
            print >>op, chain
        else:   # no barcode found; print chain unchanged
            print >>op, chain
    
    if isinstance(barcode_fasta,types.StringTypes):
        bcip.close()


def isotype_id(ighc_fasta,inhandle,outhandle):
    if isinstance(ighc_fasta,types.StringTypes):
        ighcip = open(ighc_fasta,'r')
    elif isinstance(ighc_fasta,file):
        ighcip = ighc_fasta
    
    isotypes = {}
    for record in Bio.SeqIO.parse(ighcip.'fasta'):
        isotypes[record.seq.reverse_complement().tostring()] = record.id
    
    for chain in parse_VDJXML(inhandle):
        get_tag_with_prefix(chain,'positive')   # will throw ValueError if finds non-positive chain
        for iso in isotypes.iteritems():
            if iso[0] in chain.seq[-50:]:   # arbitrary cutoff from 3' end
                chain.ighc = iso[1]
        print >>outhandle, chain
    
    if isinstance(ighc_fasta,types.StringTypes):
        ighcip.close()


def positive_strand(inhandle,outhandle):
    aligner = alignment.vdj_aligner()
    for chain in parse_VDJXML(inhandle):
        strand = aligner.seq2posstrand(chain.seq)
        chain.add_tags('positive')
        if strand == -1:
            chain.add_tags('revcomp')
            chain.seq = sequtils.reverse_complement(chain.seq)
        print >>outhandle, chain


def align_vdj(inhandle,outhandle):
    aligner = alignment.vdj_aligner()
    for chain in parse_VDJXML(inhandle):
        aln = aligner.align_chain(chain)
        print >>outhandle, chain





















# =============================
# = Clustering/CDR3 functions =
# =============================

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








#===============================================================================

# =============
# = Utilities =
# =============

def split_rep(rep,IDs):
    """Split rep into multiple Repertoire objects based on identifiers in IDs."""
    reps = np.empty(len(IDs),dtype=np.object)
    for (i,ID) in enumerate(IDs):
        reps[i] = rep.get_chains_AND(ID)
    return reps

def reps2timeseries(reps,refclones):
    """Generate time series of clones from list of Repertoire objects, using refclones as reference."""
    numClones = len(refclones)
    numRepertoires = len(reps)
    countdata = np.zeros((numClones,numRepertoires))
    for (i,rep) in enumerate(reps):
        clones = getClusters(rep)
        countdata[:,i] = countsClusters(clones,refclones)
    return countdata




#===============================================================================

# ===========================================
# = Generating specificities reference data =
# ===========================================

if not os.path.exists(os.path.join(refseq.refdatadir,refseq.imgtspecfasta)) or not os.path.exists(os.path.join(refseq.refdatadir,refseq.imgtspecvdjxml)):
    refseq.get_LIGM_with_specificities(refseq.refdatadir,refseq.imgtdat,refseq.imgtfasta,refseq.imgtspecfasta,refseq.imgtspecvdjxml)

