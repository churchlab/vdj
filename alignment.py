import warnings
import copy

import numpy as np

import seqtools

import vdj
import refseq
import sequtils
import alignmentcore

warnings.simplefilter('always')

class vdj_aligner(object):
    
    def __init__(self,**kw):
        
        self.numCrudeVCandidates = 5
        self.numCrudeDCandidates = 10
        self.numCrudeJCandidates = 2
        
        # Define seed patterns
        patternA='111011001011010111'
        patternB='1111000100010011010111'
        patternC='111111111111'
        patternD='110100001100010101111'
        patternE='1110111010001111'
        self.seedpatterns = [patternA,patternB,patternC,patternD,patternE]
        self.miniseedpatterns = ['111011','110111']
        self.patternPos = '111111111111'
        
        self.minVscore = 100    # derived from calibration data 20090710
        self.minDscore = 4
        self.minJscore = 13
        
        # set reference sequences (locus) and generate hashes from ref data
        self.locus = kw['locus']
        
        self.refV_seqs = refseq.__getattribute__(self.locus+'V_seqs')
        self.Vseqlistkeys = vdj_aligner.seqdict2kmers( self.refV_seqs, self.seedpatterns )
        
        self.refJ_seqs = refseq.__getattribute__(self.locus+'J_seqs')
        self.Jseqlistkeys = vdj_aligner.seqdict2kmers( self.refJ_seqs, self.seedpatterns )
        
        try:    # this locus may not have D segments
            self.refD_seqs = refseq.__getattribute__(self.locus+'D_seqs')
            self.Dseqlistkeysmini = vdj_aligner.seqdict2kmers( self.refD_seqs, self.miniseedpatterns )
        except AttributeError:
            pass
        
        self.refV_offset = refseq.__getattribute__(self.locus+'V_offset')
        self.refJ_offset = refseq.__getattribute__(self.locus+'J_offset')
        
        # Generate reference data for positive sequence ID
        posVseqlistkeys = vdj_aligner.seqdict2kmers( self.refV_seqs, [self.patternPos] )
        posJseqlistkeys = vdj_aligner.seqdict2kmers( self.refJ_seqs, [self.patternPos] )
        negVseqlistkeys = vdj_aligner.seqdict2kmers( vdj_aligner.seqdict2revcompseqdict(self.refV_seqs), [self.patternPos] )
        negJseqlistkeys = vdj_aligner.seqdict2kmers( vdj_aligner.seqdict2revcompseqdict(self.refJ_seqs), [self.patternPos] )
        
        # collect possible keys
        posset = set([])
        for key in posVseqlistkeys.keys():
            posset.update(posVseqlistkeys[key][self.patternPos])
        for key in posJseqlistkeys.keys():
            posset.update(posJseqlistkeys[key][self.patternPos])
        
        negset = set([])
        for key in negVseqlistkeys.keys():
            negset.update(negVseqlistkeys[key][self.patternPos])
        for key in negJseqlistkeys.keys():
            negset.update(negJseqlistkeys[key][self.patternPos])
        
        # get keys unique to positive or negative versions of reference set
        possetnew = posset - negset
        negsetnew = negset - posset
        
        self.posset = possetnew
        self.negset = negsetnew
    
    
    def Valign_chain(self,chain,verbose=False):
        query = chain.seq
        
        # compute hashes from query seq
        querykeys = vdj_aligner.seq2kmers(query,self.seedpatterns)
        
        # for each reference V segment and each pattern, how many shared k-mers are there?
        Vscores_hash = vdj_aligner.hashscore(self.Vseqlistkeys,querykeys)
        
        # get numCrudeVCandidates highest scores in Vscores and store their names in descending order
        goodVseglist = sorted(self.refV_seqs.keys(),key=lambda k: Vscores_hash[k],reverse=True)[0:self.numCrudeVCandidates]
        goodVsegdict = dict([(seg,self.refV_seqs[seg]) for seg in goodVseglist])
        
        # Needleman-Wunsch of V segment
        (bestVseg,bestVscore,bestVscoremat,bestVtracemat) = vdj_aligner.bestalignNW(goodVsegdict,query,self.minVscore)
        
        # if successful alignment
        if bestVseg is not None:
            chain.v = bestVseg
            
            # find CDR3 boundary
            # construct alignment
            Valnref,Valnrefcoords,Valnquery,Valnquerycoords = vdj_aligner.construct_alignment( self.refV_seqs[bestVseg], query, bestVscoremat, bestVtracemat )
            
            # find CDR3 boundary
            chain.v_end_idx = vdj_aligner.pruneVregion( Valnref, Valnrefcoords, Valnquery, Valnquerycoords, self.refV_offset[bestVseg] )
        
        return bestVscore
    
    
    def Jalign_chain(self,chain,verbose=False):
        # try pruning off V region for J alignment
        try:
            query = chain.seq[chain.v_end_idx:]
        except AttributeError:
            query = chain.seq
        
        # compute hashes from query seq
        querykeys = vdj_aligner.seq2kmers(query,self.seedpatterns)
        
        # for each reference J segment and each pattern, how many shared k-mers are there?
        Jscores_hash = vdj_aligner.hashscore(self.Jseqlistkeys,querykeys)
        
        # get numCrudeJCandidates highest scores in Jscores and store their names in descending order
        goodJseglist = sorted(self.refJ_seqs.keys(),key=lambda k: Jscores_hash[k],reverse=True)[0:self.numCrudeJCandidates]
        goodJsegdict = dict([(seg,self.refJ_seqs[seg]) for seg in goodJseglist])
        
        # Needleman-Wunsch of J segment
        (bestJseg,bestJscore,bestJscoremat,bestJtracemat) = vdj_aligner.bestalignNW(goodJsegdict,query,self.minJscore)
        
        # if successful alignment
        if bestJseg is not None:
            chain.j = bestJseg
            
            # find CDR3 boundary
            # construct alignment
            Jalnref,Jalnrefcoords,Jalnquery,Jalnquerycoords = vdj_aligner.construct_alignment( self.refJ_seqs[bestJseg], query, bestJscoremat, bestJtracemat )
            
            # find CDR3 boundary
            j_start_offset = vdj_aligner.pruneJregion( Jalnref, Jalnrefcoords, Jalnquery, Jalnquerycoords, self.refJ_offset[bestJseg] )
            try:
                chain.j_start_idx = chain.v_end_idx + j_start_offset
            except AttributeError:
                chain.j_start_idx = j_start_offset
        
        return bestJscore
    
    
    def Dalign_chain(self,chain,verbose=False):
        # prune off V and J regions for D alignment
        # we should not be attempting D alignment unless we have
        # a well-defined CDR3
        query = chain.junction
        
        # compute hashes from query seq
        querykeys = vdj_aligner.seq2kmers(query,self.miniseedpatterns)
        
        # for each reference D segment and each pattern, how many shared k-mers are there?
        Dscores_hash = vdj_aligner.hashscore(self.Dseqlistkeysmini,querykeys)
        
        # get numCrudeJCandidates highest scores in Jscores and store their names in descending order
        goodDseglist = sorted(self.refD_seqs.keys(),key=lambda k: Dscores_hash[k],reverse=True)[0:self.numCrudeDCandidates]
        goodDsegdict = dict([(seg,self.refD_seqs[seg]) for seg in goodDseglist])
        
        # Needleman-Wunsch of J segment
        (bestDseg,bestDscore,bestDscoremat,bestDtracemat) = vdj_aligner.bestalignSW(goodDsegdict,query,self.minDscore)
        
        # if successful alignment
        if bestDseg is not None:
            chain.d = bestDseg
        
        return bestDscore
    
    
    def align_chain(self,chain,verbose=False):
        
        if not chain.has_tag('positive') and not chain.has_tag('coding'):
            warnings.warn('chain %s may not be the correct strand' % chain.descr)
        
        scores = {}
        
        scores['v'] = self.Valign_chain(chain,verbose)
        
        scores['j'] = self.Jalign_chain(chain,verbose)
        
        # only process junction if V and J successful
        if hasattr(chain,'v') and hasattr(chain,'j'):
            chain.junction = chain.seq[chain.v_end_idx:chain.j_start_idx]
            
            # only align D if I am in a locus that has D chains
            if self.locus in ['IGH','TRB','TRD']:
                scores['d'] = self.Dalign_chain(chain,verbose)
        
        return scores
    
    
    def align_seq(self,seq,verbose=False):
        chain = vdj.ImmuneChain(descr='sequence',seq=seq)
        self.align_chain(chain,verbose)
        return chain
    
    
    def coding_chain(self,chain,verbose=False):
        strand = self.seq2coding(chain.seq)
        if strand == -1:
            chain.seq = seqtools.reverse_complement(chain.seq)
            chain.add_tag('revcomp')
        chain.add_tag('coding')
    
    
    def seq2coding(self,seq):        
        seqkeys = vdj_aligner.seq2kmers(seq,[self.patternPos])
        seqwords = seqkeys[self.patternPos]
        strandid = 1
        if len(self.negset & seqwords) > len(self.posset & seqwords):
            strandid = -1
        return strandid
    
    
    @staticmethod
    def seq2kmers(seq,patterns):
        """Given sequence and patterns, for each pattern, compute all corresponding k-mers from sequence.
        
        The result is seqannot[pattern][key]=[pos1,pos2,...,posN] in seq
                      seqkeys[pattern] = set([kmers])
        
        """
        seqkeys = {}
        patlens = []
        for pattern in patterns:
            patlens.append(len(pattern))
            seqkeys[pattern] = set()
        
        maxpatlen = max(patlens)
        
        for i in xrange(len(seq)):
            word = seq[i:i+maxpatlen]
            for pattern in patterns:
                patlen = len(pattern)
                if len(word) >= patlen:
                    key = ''
                    for j in xrange(patlen):
                        if pattern[j] == '1':
                            key += word[j]
                    seqkeys[pattern].add(key)
        
        return seqkeys
    
    
    @staticmethod
    def seqdict2kmers(seqdict,patterns):
        seqlistkeys  = {}
        for seq in seqdict.iteritems():
            seqlistkeys[seq[0]] = vdj_aligner.seq2kmers(seq[1],patterns)
        return seqlistkeys
    
    
    @staticmethod
    def hashscore(refkeys,querykeys):
        """Compute number of common keys for each reference sequence.
    
        querykeys is dict of sets, where dict keys are patterns
        reference keys is dict of ref seqs, where each elt is a
        dict of patterns with sets as values.  the patterns must be
        the same
        """
        scores = {}
        for seg in refkeys.iterkeys():
            score = 0
            for pattern in querykeys.iterkeys():
                score += len( refkeys[seg][pattern] & querykeys[pattern] )
            scores[seg] = score
        return scores
    
    
    @staticmethod
    def bestalignNW(candidatedict,query,minscore):
        bestseg = None
        bestscore = minscore
        bestscoremat = None
        besttracemat = None
        
        seq2 = query
        for (seg,seq1) in candidatedict.iteritems():
            # C implementation:
            # carve out memory
            # note that we are using zero initial conditions, so matrices are initialized too
            # notation is like Durbin p.29
            scores  = np.zeros( [len(seq1)+1, len(seq2)+1] )
            Ix = np.zeros( [len(seq1)+1, len(seq2)+1] )
            Iy = np.zeros( [len(seq1)+1, len(seq2)+1] )
            trace = np.zeros( [len(seq1)+1, len(seq2)+1], dtype=np.int)
            alignmentcore.alignNW( scores, Ix, Iy, trace, seq1, seq2 )

            currscore = vdj_aligner.scoreVJalign(scores)
            if currscore > bestscore:
                bestscore = currscore
                bestseg = seg
                bestscoremat = scores
                besttracemat = trace
        
        return (bestseg,bestscore,bestscoremat,besttracemat)
    
    
    @staticmethod
    def bestalignSW(candidatedict,query,minscore):
        bestseg = None
        bestscore = minscore
        bestscoremat = None
        besttracemat = None
        
        seq2 = query
        for (seg,seq1) in candidatedict.iteritems():
            # C implementation:
            # carve out memory
            # note that we are using zero initial conditions, so matrices are initialized too
            # notation is like Durbin p.29
            scores  = np.zeros( [len(seq1)+1, len(seq2)+1] )
            trace = np.zeros( [len(seq1)+1, len(seq2)+1], dtype=np.int)
            alignmentcore.alignSW( scores, trace, seq1, seq2 )
            
            currscore = vdj_aligner.scoreDalign(scores)
            if currscore > bestscore:
                bestscore = currscore
                bestseg = seg
                bestscoremat = scores
                besttracemat = trace
        
        return (bestseg,bestscore,bestscoremat,besttracemat)
    
    
    @staticmethod
    def pruneVregion( alnref, alnrefcoords, alnquery, alnquerycoords, offset ):
        """Prune V region out of query sequence based on alignment.
        
        Given ref and query alignments of V region, refID, and the original
        query sequence, return a sequence with the V region cut out, leaving
        the 2nd-CYS.  Also needs query alignment coords.
        
        """
        # FR3end = self.refV_offset[refID] - alnrefcoords[0]        # first candidate position  
        FR3end = offset - alnrefcoords[0]        # first candidate position  
        refgaps = alnref[:FR3end].count('-')    # count gaps up to putative CYS pos
        seengaps = 0
        while refgaps > 0:     # iteratively find all gaps up to the CYS
            seengaps += refgaps
            FR3end   += refgaps     # adjust if for gaps in ref alignment
            refgaps   = alnref[:FR3end].count('-') - seengaps   # any add'l gaps?
        
        querygaps = alnquery[:FR3end].count('-')
        
        # v_end_idx = idx of start of aln of query + distance into aln - # of gaps
        v_end_idx = alnquerycoords[0] + FR3end - querygaps
        
        return v_end_idx
    
    
    @staticmethod
    def pruneJregion( alnref, alnrefcoords, alnquery, alnquerycoords, offset ):
        """Prune J region out of query sequence based on alignment.
        
        Given ref and query alignments of J region, refID, and the original
        query sequence, return a sequence with the J region cut out, leaving
        the J-TRP.  Also needs query alignment coords.
        
        """
        # FR4start = self.refJ_offset[refID] - alnrefcoords[0]  # first candidate position of J-TRP start
        FR4start = offset - alnrefcoords[0]  # first candidate position of J-TRP start
        refgaps = alnref[:FR4start].count('-')  # count gaps up to putative TRP pos
        seengaps = 0
        while refgaps > 0:     # iteratively find all gaps up to the TRP
            seengaps += refgaps
            FR4start += refgaps     # adjust for gaps in ref alignment
            refgaps   = alnref[:FR4start].count('-') - seengaps # any add'l gaps?
        
        querygaps = alnquery[:FR4start].count('-')
        
        # j_start_offset = idx of start of aln of query + distance into aln - # of gaps
        # note: j_start_offset is from the pruned query seq
        j_start_offset = alnquerycoords[0] + FR4start - querygaps
        
        return j_start_offset
    
    
    @staticmethod
    def construct_alignment(seq1,seq2,scoremat,tracemat):
        """Construct alignment of ref segment to query from score and trace matrices."""
        nrows,ncols = scoremat.shape
        
        # do some error checking
        if len(seq1)+1 != nrows or len(seq2)+1 != ncols:
            raise Exception, "nrows and ncols must be equal to len(seq1)+1 and len(seq2)+1"
        
        # translate integer traces to coords
        deltas = {
            0 : (1,1),
            1 : (1,0),
            2 : (0,1),
            3 : (0,0)
        }
        
        # compute col where alignment should start
        if nrows <= ncols:
            col = np.argmax( scoremat[nrows-1,:] )
            row = nrows-1
        else:
            col = ncols-1
            row = np.argmax( scoremat[:,ncols-1] )
        
        # if row is coord in matrix, row-1 is coord in seq (b/c of init conditions)
        
        aln1 = seq1[row-1]
        aln2 = seq2[col-1]
        
        aln1end = row
        aln2end = col
        
        while (row-1 > 0) and (col-1 > 0):
            # compute direction of moves
            rowchange,colchange = deltas[ tracemat[row,col] ]
            
            # emit appropriate symbols
            if rowchange == 1:
                row -= 1
                aln1 = seq1[row-1] + aln1
            elif rowchange == 0:
                aln1 = '-' + aln1
            else:
                raise Exception, "Trace matrix contained jump of greater than one row/col."
            
            if colchange == 1:
                col -= 1
                aln2 = seq2[col-1] + aln2
            elif colchange == 0:
                aln2 = '-' + aln2
            else:
                raise Exception, "Trace matrix contained jump of greater than one row/col."
        
        aln1start = row-1
        aln2start = col-1
        
        # the coords refer to coords in the sequence (pythonic)
        return aln1, (aln1start,aln1end), aln2, (aln2start,aln2end)
    
    
    @staticmethod
    def scoreVJalign(scorematrix):
        """Computes score of V alignment given Needleman-Wunsch score matrix
        
        ASSUMES num rows < num cols, i.e., refseq V seg is on vertical axis
        
        """
        nrows,ncols = scorematrix.shape
        
        if nrows <= ncols:
            return np.max( scorematrix[nrows-1,:] )
        else:
            return np.max( scorematrix[:,ncols-1] )
    
    
    @staticmethod
    def scoreDalign(scorematrix):
        """Computes score of D alignment given Smith-Waterman score matrix
        
        """
        return np.max( scorematrix )
    
    
    @staticmethod
    def seqdict2revcompseqdict(seqdict):
        revcompdict = {}
        for item in seqdict.iteritems():
            revcompdict[item[0]] = sequtils.reverse_complement(item[1])
        return revcompdict


class vdj_aligner_combined(object):
    """vdj aligner for 'light' chains
    
    this class will perform alignment for both loci, e.g., IGK and IGL
    and pick the one with the better V score
    """
    def __init__(self,**kw):
        self.loci = kw['loci']
        self.aligners = [vdj_aligner(locus=locus) for locus in self.loci]
        
        self.patternPos = '111111111111'
        self.posset = set()
        self.negset = set()
        for aligner in self.aligners:
            self.posset.update(aligner.posset)
            self.negset.update(aligner.negset)
    
    def align_chain(self,chain,verbose=False):
        alignments = []
        for aligner in self.aligners:
            curr_chain = copy.deepcopy(chain)
            curr_score = aligner.align_chain(curr_chain)
            alignments.append((curr_chain,curr_score))
        alignments = sorted(filter(lambda a: hasattr(a[0],'v'),alignments),key=lambda a:a[1]['v'],reverse=True)
        if len(alignments) > 0:
            bestchain = alignments[0][0]
            if hasattr(bestchain,'v'):
                chain.v = bestchain.v
                chain.v_end_idx = bestchain.v_end_idx
            if hasattr(bestchain,'j'):
                chain.j = bestchain.j
                chain.j_start_idx = bestchain.j_start_idx
            if hasattr(bestchain,'junction'):
                chain.junction = bestchain.junction
            if hasattr(bestchain,'d'):
                chain.d = bestchain.d
            return alignments[0][1]     # NOTE: I only return the scores upon successful aln
    
    def align_seq(self,seq,verbose=False):
        chain = vdj.ImmuneChain(descr='sequence',seq=seq)
        self.align_chain(chain,verbose)
        return chain
    
    def coding_chain(self,chain):
        strand = self.seq2coding(chain.seq)
        if strand == -1:
            chain.seq = seqtools.reverse_complement(chain.seq)
            chain.add_tag('revcomp')
        chain.add_tag('coding')
    
    def seq2coding(self,seq):
        seqkeys = vdj_aligner.seq2kmers(seq,[self.patternPos])
        seqwords = seqkeys[self.patternPos]
        strandid = 1
        if len(self.negset & seqwords) > len(self.posset & seqwords):
            strandid = -1
        return strandid

def igh_aligner():
    return vdj_aligner(locus='IGH')

def igk_aligner():
    return vdj_aligner(locus='IGK')

def igl_aligner():
    return vdj_aligner(locus='IGL')

def igkl_aligner():
    return vdj_aligner_combined(loci=['IGK','IGL'])

def trb_aligner():
    return vdj_aligner(locus='TRB')

def tra_aligner():
    return vdj_aligner(locus='TRA')

def trd_aligner():
    return vdj_aligner(locus='TRD')

def trg_aligner():
    return vdj_aligner(locus='TRG')
