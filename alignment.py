import time
import operator

import numpy as np

import subprocess

import vdj
import refseq
import sequtils
import maligner
import hasher

class vdj_aligner(object):
    
    def __init__(self,hashed=True,verbose=False):
        
        self.hashed = hashed
        
        # instantiate multialigner objects
        self.Vma = maligner.MAligner()
        self.Dma = maligner.MAligner()
        self.Jma = maligner.MAligner()
        
        # load reference sequences
        vdj_aligner.load_references(self.Vma,refseq.IGHV_seqs)
        vdj_aligner.load_references(self.Dma,refseq.IGHD_seqs)
        vdj_aligner.load_references(self.Jma,refseq.IGHJ_seqs)
        
        if hashed:
            # instantiate hasher objects
            self.Vhasher = vdj.hasher.SequenceHasher()
            self.Dhasher = vdj.hasher.SequenceHasher()
            self.Jhasher = vdj.hasher.SequenceHasher()
            
            # load reference sequences
            vdj_aligner.load_references(self.Vhasher,refseq.IGHV_seqs)
            vdj_aligner.load_references(self.Dhasher,refseq.IGHD_seqs)
            vdj_aligner.load_references(self.Jhasher,refseq.IGHJ_seqs)
            
            self.num_V_candidates = 5
            self.num_D_candidates = 10
            self.num_J_candidates = 2
    
    @staticmethod
    def load_references(objinst,seqdict):
        for (name,seq) in seqdict.iteritems():
            objinst.addReference(name,seq)
        return
    
    def hash_V(self,seq,verbose=False):
        return self.Vhasher.hash(seq)
    
    def hash_D(self,seq,verbose=False):
        return self.Dhasher.hash(seq)
    
    def hash_J(self,seq,verbose=False):
        return self.Jhasher.hash(seq)
    
    def align_V(self,seq,reference_segments=None,verbose=False):
        return self.Vma.align(seq,reference_segments)
    
    def align_D(self,seq,reference_segments=None,verbose=False):
        return self.Dma.align(seq,reference_segments)
    
    def align_J(self,seq,reference_segments=None,verbose=False):
        return self.Jma.align(seq,reference_segments)
    
    def align_seq(self,seq,verbose=False):
        chain = vdj.ImmuneChain(descr='sequence',seq=seq)
        self.align_chain(chain,verbose)
        return chain
    
    def align_chain(self,chain,verbose=False):
        
        query = chain.seq
        
        # align V
        if self.hashed:
            candidateVsegs = self.hash_V(query,verbose)[:self.num_V_candidates]
        else:
            candidateVsegs = refseq.IGHV_list
        (bestV,ref_gapped,query_gapped) = self.align_V(query,candidateVsegs,verbose)
        chain.v = bestV
        
        # prune V
        v_end_idx = self.compute_v_end_idx(bestV,ref_gapped,query_gapped)
        chain.add_tags('v_end_idx|%d'%v_end_idx)
        query = query[v_end_idx:]
        
        # align J
        if self.hashed:
            candidateJsegs = self.hash_J(query,verbose)[:self.num_J_candidates]
        else:
            candidateJsegs = refseq.IGHJ_list
        (bestJ,ref_gapped,query_gapped) = self.align_J(query,candidateJsegs,verbose)
        chain.j = bestJ
        
        # prune J
        j_start_idx = self.compute_j_start_idx(bestJ,ref_gapped,query_gapped)
        chain.add_tags('j_start_idx|%d'%(j_start_idx))
        query = query[:j_start_idx]
        
        if bestV != '' and bestJ != '':
            chain.junction = query
            
            # align D, only if both V and J were successfully pruned
            # TODO: NOT IMPLEMENTED
        
        return
    
    @staticmethod
    def compute_v_end_idx(ref_id,ref_gapped,query_gapped):
        FR3_end = refseq.IGHV_offset[ref_id]
        ref_gaps = ref_gapped[:FR3_end].count('-') # count gaps up to putative CYS pos
        seen_gaps = 0
        while ref_gaps > 0:
            seen_gaps += ref_gaps
            FR3_end += ref_gaps
            ref_gaps = ref_gapped[:FR3_end].count('-') - seen_gaps
        v_end_idx = FR3_end - query_gapped[:FR3_end].count('-')
        return v_end_idx
    
    @staticmethod
    def compute_j_start_idx(ref_id,ref_gapped,query_gapped):
        FR4_start = refseq.IGHJ_offset[ref_id] + 3  # extra 3 for J-TRP
        ref_gaps = ref_gapped[:FR4_start].count('-')
        seen_gaps = 0
        while ref_gaps > 0:
            seen_gaps += ref_gaps
            FR4_start += ref_gaps
            ref_gaps = ref_gapped[:FR4_start].count('-') - seen_gaps
        j_start_idx = -1 * (len(ref_gapped) - FR4_start - query_gapped[FR4_start:].count('-'))
        return j_start_idx
