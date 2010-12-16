import types
import xml.etree.cElementTree as ElementTree

import numpy as np
import Bio.SeqIO

import seqtools

# import refseq
# import sequtils
# import alignment
# import clustering
# import LSF
# import params
# import analysis



# ===================
# = DATA STRUCTURES =
# ===================

class ImmuneChain(object):
    """Data structure to represent an immune chain."""
    
    def __init__(self,**kw):
        """Initialize ImmuneChain
        
        seq is 5'->3'
        """
        def kw_init(attrib):
            if kw.has_key(attrib):
                self.__setattr__(attrib,kw[attrib])
        
        kw_init('seq')
        kw_init('descr')
        kw_init('locus')
        kw_init('v')
        kw_init('d')
        kw_init('j')
        kw_init('c')
        kw_init('junction')
        
        if kw.has_key('tags'):
            tags = kw['tags']
            if isinstance(tags,types.StringTypes): tags = [tags]
            self.tags = set(tags)
        else:
            self.tags = set([])
    
    def get_cdr3(self):
        return len(self.junction)
    def set_cdr3(self,value):
        pass
    cdr3 = property(fget=get_cdr3,fset=set_cdr3)
    
    def get_vj(self):
        return '|'.join([self.v,self.j])
    def set_vj(self):
        pass
    vj = property(fget=get_vj,fset=set_vj)
    
    def get_vdj(self):
        return '|'.join([self.v,self.d,self.j])
    def set_vdj(self):
        pass
    vdj = property(fget=get_vdj,fset=set_vdj)
    
    def add_tags(self,tagset):
        if isinstance(tagset,types.StringTypes): tagset = [tagset]
        self.tags.update(tagset)
    
    def add_tag(self,tag):
        self.add_tags(tag)
    
    def remove_tags(self,tagset):
        if isinstance(tagset,types.StringTypes): tagset = [tagset]
        for tag in tagset: self.tags.remove(tag)
    
    def remove_tag(self,tag):
        self.remove_tags(tag)
    
    def has_tag(self,tag):
        if tag in self.tags:
            return True
        else:
            return False
    
    def __len__(self):
        return len(self.seq)
    
    def __str__(self):
        return self.__repr__()
    
    def __repr__(self):
        return self.get_XML()
    
    def get_XML(self):
        format_xml = lambda attrib,value: "\t<%(attrib)s>%(value)s</%(attrib)s>\n" % {'attrib':attrib,'value':value}
        xmlstring = '<ImmuneChain>\n'
        for (attrib,value) in self.__dict__.iteritems():
            if attrib == 'tags':
                for tag in self.tags:
                    xmlstring += format_xml('tag',tag)
            else:
                xmlstring += format_xml(attrib,value)
        xmlstring += '</ImmuneChain>\n'
        return xmlstring



# ================
# = INPUT/OUTPUT =
# ================

class ParserVDJXML(object):
    """Parser for VDJXML"""
    def __init__(self):
        self.chain = None
    
    def start_handler(self,elem):
        if elem.tag == 'ImmuneChain':
            self.chain = ImmuneChain()
    
    def end_handler(self,elem):
        if elem.tag == 'tag':
            self.chain.add_tags(elem.text)
        elif elem.tag == 'v_end_idx' or elem.tag == 'j_start_idx':
            self.chain.__setattr__(elem.tag,int(elem.text))
        else:
            self.chain.__setattr__(elem.tag,elem.text)
    
    def parse(self,inputfile):
        for event, elem in ElementTree.iterparse(inputfile,events=('start','end')):
            if event == 'start':
                if elem.tag == 'root':  # to ensure non-incorp of <root> obj in chain
                    pass
                else:
                    self.start_handler(elem)
            elif event == 'end':
                if elem.tag == 'ImmuneChain':
                    yield self.chain
                elif elem.tag == 'root':    # to ensure clean exit at end of file
                    pass
                else:
                    self.end_handler(elem)


class PredicateParserVDJXML(ParserVDJXML):
    """VDJXML Parser that takes a predicate function"""
    def __init__(self,predicate):
        ParserVDJXML.__init__(self)
        self.predicate = predicate
    
    def parse(self,inputfile):
        for event, elem in ElementTree.iterparse(inputfile,events=('start','end')):
            if event == 'start':
                self.start_handler(elem)
            elif event == 'end':
                if elem.tag == 'ImmuneChain':
                    if self.predicate(self.chain) == True:
                        yield self.chain
                else:
                    self.end_handler(elem)


def parse_VDJXML(inputfile):
    vdjxmlparser = ParserVDJXML()
    return vdjxmlparser.parse(inputfile)


def filter_parse_VDJXML(inputfile,predicate):
    vdjxmlparser = PredicateParserVDJXML(predicate)
    return vdjxmlparser.parse(inputfile)



# ==========================================================
# ==========================================================
# = DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED =
# ==========================================================
# = DEPRECATED DEPRECATED DEPRECATED DEPRECATED DEPRECATED =
# ==========================================================
# ==========================================================

# # ============
# # = Counting =
# # ============
# 
# def counts_VJ(inputfile):
#     if isinstance(inputfile,types.StringTypes):
#         ip = open(inputfile,'r')
#     elif isinstance(inputfile,file):
#         ip = inputfile
#     
#     counts = np.zeros( (len(refseq.IGHV_list),len(refseq.IGHJ_list)) )
#     for chain in parseVDJXML(ip):
#         counts[refseq.IGHV_idx[chain.v],refseq.IGHJ_idx[chain.j]] += 1
#     
#     if isinstance(inputfile,types.StringTypes):
#         ip.close()
#     
#     return counts
# 
# 
# def counts_VDJ(rep):
#     cn = np.zeros( (len(refseq.IGHV_list),len(refseq.IGHD_list),len(refseq.IGHJ_list)) )
#     for chain in rep.chains:
#         cn[refseq.IGHV_idx[chain.v],refseq.IGHD_idx[chain.d],refseq.IGHJ_idx[chain.j]] += 1
#     return cn
# 
# 
# def reshape_counts_VDJ_2D(counts):
#     return counts.reshape(len(refseq.IGHV_list),len(refseq.IGHD_list)*len(refseq.IGHJ_list))
# 
# 
# def count_dict_clone_idxs(clone_idxs,reference_clones=None):
#     """Takes a dictionary of cluster names mapped to a sequence of indices into an ImmuneChain list.
#     
#     Returns an np array of the same length as reference_clones with the counts of each
#     cluster in reference_clones.
#     
#     The need for reference_clones is due to the fact that splitting a given repertoire
#     may result in some parts not observing any of a given clone, so there needs to be a common way
#     to compare two clone sets.
#     
#     If reference_clones is left out, then the set of clones present in clone_idxs is used.
#     
#     """
#     if reference_clones == None:
#         reference_clones = clone_idxs.keys()
#     counts = np.zeros(len(reference_clones))
#     for (i,name) in enumerate(reference_clones):
#         counts[i] = len(clone_idxs.get(name,[]))
#     return counts
# 
# 
# def count_dict_clone_counts(clone_counts,reference_clones=None):
#     if reference_clones == None:
#         reference_clones = clone_counts.keys()
#     counts = np.zeros(len(reference_clones))
#     for (i,name) in enumerate(reference_clones):
#         counts[i] = clone_counts.get(name,0)
#     return counts
# 
# 
# # =================================
# # = Retrieving tags and filtering =
# # =================================
# 
# def get_tag_with_prefix(chain,prefix):
#     for tag in chain.tags:
#         if tag.startswith(prefix):
#             return tag
#     raise ValueError, "Tag that starts with " + prefix + " not found."
# 
# 
# def get_clone(chain):
#     return get_tag_with_prefix(chain,'clone')
# 
# 
# def get_barcode(chain):
#     try:
#         return get_tag_with_prefix(chain,'barcode')
#     except ValueError:
#         return ''
# 
# 
# def filter_tags_and(tags,inhandle,outhandle):
#     if isinstance(tags,types.StringTypes): tags = [tags]
#     tags = set(tags)
#     for chain in parse_VDJXML(inhandle):
#         if tags <= chain.all_tags:    # test that everything in tags is in all_tags
#             print >>outhandle, chain
# 
# 
# def filter_tags_or(tags,inhandle,outhandle):
#     if isinstance(tags,types.StringTypes): tags = [tags]
#     tags = set(tags)
#     empty_set = set()
#     for chain in parse_VDJXML(inhandle):
#         if tags & chain.all_tags != empty_set:    # test that tags and all_tags share something
#             print >>outhandle, chain
# 
# 
# def is_full_VJ(chain):
#     if (chain.v in refseq.IGHV_seqs.keys()) and (chain.j in refseq.IGHJ_seqs.keys()):
#         return True
#     else:
#         return False
# 
# 
# def get_clone_idxs(inhandle):
#     clusters = {}
#     i = 0
#     for chain in parse_VDJXML(inhandle):
#         try: clusters[get_clone(chain)] += [i]
#         except KeyError: clusters[get_clone(chain)] = [i]
#         i += 1
#     return clusters
# 
# 
# def get_clone_counts(inhandle):
#     clusters = {}
#     for chain in parse_VDJXML(inhandle):
#         try: clusters[get_clone(chain)] += 1
#         except KeyError: clusters[get_clone(chain)] = 1
#     return clusters
# 
# 
# # ======================
# # = Pipeline functions =
# # ======================
# 
# def vdjxml2fasta(inhandle,outhandle):
#     for chain in parse_VDJXML(inhandle):
#         print >>outhandle, '>'+chain.descr
#         print >>outhandle, chain.seq
# 
# 
# # for generating identifiers from VJ combos
# def vj_id(v_seg,j_seg):
#     return seqtools.cleanup_id(v_seg)+'_'+seqtools.cleanup_id(j_seg)
# 
# 
# def split_vdjxml_into_VJ_parts(inhandle,outname):
#     parts = []
#     vj_ids = []
#     outhandles = {}
#     
#     # open output files for all VJ combos
#     i = 0
#     for v_seg in refseq.IGHV_seqs.keys():
#         for j_seg in refseq.IGHJ_seqs.keys():
#             curr_outname = outname + '.' + str(i)
#             curr_vj_id = vj_id(v_seg,j_seg)
#             parts.append(curr_outname)
#             vj_ids.append(curr_vj_id)
#             outhandles[curr_vj_id] = open(curr_outname,'w')
#             i += 1
#     
#     for chain in parse_VDJXML(inhandle):
#         curr_vj_id = vj_id(chain.v,chain.j)
#         print >>outhandles[curr_vj_id], chain
#     
#     for handle in outhandles.itervalues():
#         handle.close()
#     
#     return (parts,vj_ids)
# 
# 
# def parse_VDJXML_parts(parts):
#     for part in parts:
#         for chain in parse_VDJXML(part):
#             yield chain
# 
# 
# def wait_for_subprocesses(process_list,interval=30):
#     finished = False
#     while not finished:
#         finished = True
#         time.sleep(interval)
#         for p in process_list:
#             if p.poll() == None:
#                 finished = False
#                 break
