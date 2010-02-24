import time
import operator

import numpy as np

import subprocess

import vdj
import refseq
import sequtils
import maligner

class vdj_full_aligner(object):
    
    def __init__(self,verbose=False):
        
        # instantiate multialigner objects
        self.Vma = maligner.MAligner()
        self.Dma = maligner.MAligner()
        self.Jma = maligner.MAligner()
        
        # load reference sequences
        vdj_full_aligner.load_references(self.Vma,refseq.IGHV_seqs)
        vdj_full_aligner.load_references(self.Dma,refseq.IGHD_seqs)
        vdj_full_aligner.load_references(self.Jma,refseq.IGHJ_seqs)
    
    @staticmethod
    def load_references(aligner,seqdict):
        for (name,seq) in seqdict.iteritems():
            aligner.addEntry(name,seq)
        return
    
    def align_V(self,seq,verbose=False):
        return self.Vma.align(seq)
    
    def align_D(self,seq,verbose=False):
        return self.Dma.align(seq)
    
    def align_J(self,seq,verbose=False):
        return self.Jma.align(seq)
    
    def align_seq(self,seq,verbose=False):
        chain = vdj.ImmuneChain(descr='sequence',seq=seq)
        self.align_chain(chain,verbose)
        return chain
    
    def align_chain(self,chain,verbose=False):
        
        query = chain.seq
        
        # align V
        # TODO: this should also return the full alignment for pruning afterwards
        bestV = self.align_V(query,verbose)
        chain.v = bestV
        
        
        # prune V
        # modify query variable
        # chain.add_tags('v_end_idx|%d'%v_end_idx)
        
        
        # align J
        bestJ = self.align_J(query,verbose)
        chain.j = bestJ
        
        # prune J
        # modify query variable
        # chain.add_tags('j_start_idx|%d'%(v_end_idx+j_start_idx))
        
        
        # align D, only if both V and J were successfully pruned
        if chain.v != '' and chain.j != '':
            bestD = self.align_D(query,verbose)
            chain.d = bestD
    
