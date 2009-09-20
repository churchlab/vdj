import time
import operator

import numpy as np

import vdj
import refseq
import sequtils
import alignmentcore

class vdj_aligner(object):
    
    def __init__(self, verbose=False):
        
        t0 = time.time()
        
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
        
        # Generate hashes from reference data for sequence alignment
        self.Vseqlistannot,self.Vseqlistkeys = vdj_aligner.seqdict2kmerannot( refseq.IGHV_seqs, self.seedpatterns )
        self.Dseqlistannotmini,self.Dseqlistkeysmini = vdj_aligner.seqdict2kmerannot( refseq.IGHD_seqs, self.miniseedpatterns )
        self.Jseqlistannot,self.Jseqlistkeys = vdj_aligner.seqdict2kmerannot( refseq.IGHJ_seqs, self.seedpatterns )
        
        # Generate reference data for positive sequence ID
        posVseqlistannot,posVseqlistkeys = vdj_aligner.seqdict2kmerannot( refseq.IGHV_seqs, [self.patternPos] )
        posJseqlistannot,posJseqlistkeys = vdj_aligner.seqdict2kmerannot( refseq.IGHJ_seqs, [self.patternPos] )
        negVseqlistannot,negVseqlistkeys = vdj_aligner.seqdict2kmerannot( vdj_aligner.seqdict2revcompseqdict(refseq.IGHV_seqs), [self.patternPos] )
        negJseqlistannot,negJseqlistkeys = vdj_aligner.seqdict2kmerannot( vdj_aligner.seqdict2revcompseqdict(refseq.IGHJ_seqs), [self.patternPos] )
        
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
        
        t1 = time.time()
        
        if verbose: print "Database init:", t1-t0
    
    def align_seq(self,seq,verbose=False):
        chain = vdj.ImmuneChain(descr='sequence',seq=seq)
        self.align_chain(chain,verbose)
        return chain
    
    def align_chain(self,chain,verbose=False):
        
        query = chain.seq
        
        t0 = time.time()
        
        # compute hashes from query seq
        queryannot,querykeys = vdj_aligner.seq2kmerannot(query,self.seedpatterns)
        
        t1 = time.time()
        
        Vscores_hash = {}
        
        # for each reference V segment and each pattern, how many shared k-mers are there?
        for Vseg in refseq.IGHV_seqs.keys():
            score = 0
            for pattern in self.seedpatterns:
                score += len( self.Vseqlistkeys[Vseg][pattern] & querykeys[pattern] )
            Vscores_hash[Vseg] = score
        
        # get numCrudeVCandidates highest scores in Vscores and store their names in descending order
        goodVseglist = [ seg[0] for seg in vdj_aligner.dict2sorteddecreasingitemlist(Vscores_hash,'value')[0:self.numCrudeVCandidates] ]
        
        if verbose:
            print goodVseglist
        
        t2 = time.time()
        
        bestVseg = ''
        bestVscore = 100    # derived from calibration data 20090710
        bestVscoremat = []
        bestVtracemat = []
        
        # perform Needleman-Wunsch on top V seg candidates and remember which had the highest score
        for goodVseg in goodVseglist:
            # C implementation:
            # carve out memory
            # note that we are using zero initial conditions, so matrices are initialized too
            # notation is like Durbin p.29
            seq1 = refseq.IGHV_seqs[goodVseg]
            seq2 = query
            scores  = np.zeros( [len(seq1)+1, len(seq2)+1] )
            Ix = np.zeros( [len(seq1)+1, len(seq2)+1] )
            Iy = np.zeros( [len(seq1)+1, len(seq2)+1] )
            trace = np.zeros( [len(seq1)+1, len(seq2)+1], dtype=np.int)
            alignmentcore.alignNW( scores, Ix, Iy, trace, seq1, seq2 )
            
            currscore = vdj_aligner.scoreVJalign(scores)
            if currscore > bestVscore:
                bestVscore = currscore
                bestVseg = goodVseg
                bestVscoremat = scores
                bestVtracemat = trace
        
        chain.v = bestVseg
        
        t3 = time.time()
        
        # reconstruct the alignment and chop off V region through beginning of CDR3 (IMGT)
        v_end_idx = 0   # to ensure it gets defined (for processing j_start_idx below)
        if bestVseg != '':
            Valnref,Valnrefcoords,Valnquery,Valnquerycoords = vdj_aligner.construct_alignment( refseq.IGHV_seqs[bestVseg], query, bestVscoremat, bestVtracemat )
            query,v_end_idx = vdj_aligner.pruneVregion( Valnref, Valnrefcoords, Valnquery, Valnquerycoords, bestVseg, query )
            chain.add_tags('v_end_idx|%d'%v_end_idx)
        
        t4 = time.time()
        
        # compute hashes from (pruned) query seq (junction + J)
        queryannot,querykeys = vdj_aligner.seq2kmerannot(query,self.seedpatterns)
        
        t5 = time.time()
        
        Jscores_hash = {}
        
        # for each reference J segment and each pattern, how many shared k-mers are there?
        for Jseg in refseq.IGHJ_seqs.keys():
            score = 0
            for pattern in self.seedpatterns:
                score += len( self.Jseqlistkeys[Jseg][pattern] & querykeys[pattern] )
            Jscores_hash[Jseg] = score
        
        # get numCrudeJCandidates highest scores in Jscores and store their names in descending order
        goodJseglist = [ seg[0] for seg in vdj_aligner.dict2sorteddecreasingitemlist(Jscores_hash,'value')[0:self.numCrudeJCandidates] ]
        
        if verbose:
            print goodJseglist
        
        t6 = time.time()
        
        bestJseg = ''
        bestJscore = 13     # derived from calibration data 20090710
        bestJscoremat = []
        bestJtracemat = []
        
        # perform Needleman-Wunsch on top J seg candidates and remember which had the highest score
        for goodJseg in goodJseglist:
            # C implementation:
            # carve out memory
            # note that we are using zero initial conditions, so matrices are initialize too
            # notation is like Durbin p.29
            seq1 = refseq.IGHJ_seqs[goodJseg]
            seq2 = query
            scores  = np.zeros( [len(seq1)+1, len(seq2)+1] )
            Ix = np.zeros( [len(seq1)+1, len(seq2)+1] )
            Iy = np.zeros( [len(seq1)+1, len(seq2)+1] )
            trace = np.zeros( [len(seq1)+1, len(seq2)+1], dtype=np.int)
            alignmentcore.alignNW( scores, Ix, Iy, trace, seq1, seq2 )
            
            # pure python:
            #scores,trace = vdj_aligner.alignNW( refseq.IGHJ_seqs[goodJseg], query )
            
            currscore = vdj_aligner.scoreVJalign(scores)
            if currscore > bestJscore:
                bestJscore = currscore
                bestJseg = goodJseg
                bestJscoremat = scores
                bestJtracemat = trace
        
        chain.j = bestJseg
        
        t7 = time.time()
        
        # reconstruct the alignment and chop off J region after the TRP (IMGT)
        if bestJseg != '':
            Jalnref,Jalnrefcoords,Jalnquery,Jalnquerycoords = vdj_aligner.construct_alignment( refseq.IGHJ_seqs[bestJseg], query, bestJscoremat, bestJtracemat )
            (query,j_start_idx) = vdj_aligner.pruneJregion( Jalnref, Jalnrefcoords, Jalnquery, Jalnquerycoords, bestJseg, query )
            chain.add_tags('j_start_idx|%d'%(v_end_idx+j_start_idx))
        
        t8 = time.time()
        
        # only attempt D alignment if both V and J were successful
        if bestVseg != '' and bestJseg != '':
            chain.junction = query
            chain.cdr3 = len(query)
            
            # compute hashes from junction sequence using mini seed patterns
            queryannot,querykeys = vdj_aligner.seq2kmerannot(query,self.miniseedpatterns)
            
            t9 = time.time()
            
            Dscores_hash = {}
            
            # for each reference D segment and each pattern, how many shared k-mers are there?
            for Dseg in refseq.IGHD_seqs.keys():
                score = 0
                for pattern in self.miniseedpatterns:
                    score += len( self.Dseqlistkeysmini[Dseg][pattern] & querykeys[pattern] )
                Dscores_hash[Dseg] = score
            
            # get numCrudeDCandidates highest scores in Dscores and store their names in descending order
            goodDseglist = [ seg[0] for seg in vdj_aligner.dict2sorteddecreasingitemlist(Dscores_hash,'value')[0:self.numCrudeDCandidates] ]
            
            if verbose:
                print goodDseglist
            
            t10 = time.time()
            
            bestDseg = ''
            bestDscore = 4      # derived from calibration data 20090710
            bestDscoremat = []
            bestDtracemat = []
            
            # perform Smith-Waterman on top D seg candidates and remember which had the highest score
            for goodDseg in goodDseglist:
                # C implementation:
                # carve out memory
                # note that we are using zero initial conditions, so matrices are initialize too
                # notation is like Durbin p.29
                seq1 = refseq.IGHD_seqs[goodDseg]
                seq2 = query
                scores  = np.zeros( [len(seq1)+1, len(seq2)+1] )
                trace = np.zeros( [len(seq1)+1, len(seq2)+1], dtype=np.int)
                alignmentcore.alignSW( scores, trace, seq1, seq2 )
                
                # pure python:
                #scores,trace = vdj_aligner.alignSW( refseq.IGHD_seqs[goodDseg], query )
                
                currscore = vdj_aligner.scoreDalign(scores)
                if currscore > bestDscore:
                    bestDscore = currscore
                    bestDseg = goodDseg
                    bestDscoremat = scores
                    bestDtracemat = trace
            
            t11 = time.time()
            
            chain.d = bestDseg
            
        else:
            bestDseg = ''
            t9 = t8
            t10 = t8
            t11 = t8
        
        if verbose:
            print t1-t0, "Full query hashes"
            print t2-t1, "Comparing hashes to V database"
            print t3-t2, "V seg NW alignment"
            print t4-t3, "Construct alignment and prune V region off"
            print t5-t4, "Pruned query hashes"
            print t6-t5, "Comparing hashes to J database"
            print t7-t6, "J seg NW alignment"
            print t8-t7, "Construct alignment and prune J region off"
            print t9-t8, "Pruned query hashes (junction only)"
            print t10-t9, "Comparing hashes to D database"
            print t11-t10, "D seg SW alignment"
            print t11-t0, "Total time"
        
        return chain.v, chain.d, chain.j, chain.junction
    
    def seq2posstrand(self,seq):        
        seqannot,seqkeys = vdj_aligner.seq2kmerannot(seq,[self.patternPos])
        seqwords = seqkeys[self.patternPos]
        strandid = 1
        if len(self.negset & seqwords) > len(self.posset & seqwords):
            strandid = -1
        return strandid
    
    @staticmethod
    def seq2kmerannot(seq,patterns):
        """Given sequence and patterns, for each pattern, compute all corresponding k-mers from sequence.
        
        The result is seqannot[pattern][key]=[pos1,pos2,...,posN] in seq
                      seqkeys[pattern] = set([kmers])
        
        """
        seqannot = {}
        patlens = []
        for pattern in patterns:
            patlens.append(len(pattern))
            seqannot[pattern] = {}
        
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
                    prevkmers = seqannot[pattern].get(key,[])
                    seqannot[pattern][key] = prevkmers + [i]
        
        seqkeys = {}
        for pattern in patterns:
            seqkeys[pattern] = set( seqannot[pattern].keys() )
        
        return seqannot,seqkeys
    
    @staticmethod
    def seqdict2kmerannot(seqdict,patterns):
        seqlistannot = {}
        seqlistkeys  = {}
        for seq in seqdict.iteritems():
            seqannot,seqkeys = vdj_aligner.seq2kmerannot(seq[1],patterns)
            seqlistannot[seq[0]] = {}
            seqlistkeys[seq[0]]  = {}
            for pattern in patterns:
                seqlistannot[seq[0]][pattern] = seqannot[pattern]
                seqlistkeys[seq[0]][pattern]  = seqkeys[pattern]
        return seqlistannot,seqlistkeys
    
    @staticmethod
    def pruneVregion( alnref, alnrefcoords, alnquery, alnquerycoords, refID, queryseq ):
        """Prune V region out of query sequence based on alignment.
        
        Given ref and query alignments of V region, refID, and the original
        query sequence, return a sequence with the V region cut out, leaving
        the 2nd-CYS.  Also needs query alignment coords.
        
        """
        
        #DEBUG
        # # check that alnref actually has the whole reference segment
        #       # otherwise, I would need to pass something like alnrefcoords
        #       if alnref.replace('-','') != refseq.IGHV_seqs[refID]:
        #           raise Exception, "Aligned reference segment is not equal to vdj.refseq reference segment."
        
        FR3end = refseq.IGHV_offset[refID] - alnrefcoords[0]        # first candidate position  
        #FR3end = refseq.IGHV_offset[refID]     # first candidate position  
        refgaps = alnref[:FR3end].count('-')    # count gaps up to putative CYS pos
        seengaps = 0
        while refgaps != 0:     # iteratively find all gaps up to the CYS
            seengaps += refgaps
            FR3end   += refgaps     # adjust if for gaps in ref alignment
            refgaps   = alnref[:FR3end].count('-') - seengaps   # any add'l gaps?
        
        querygaps = alnquery[:FR3end].count('-')
        
        # v_end_idx = idx of start of aln of query + distance into aln - # of gaps
        v_end_idx = alnquerycoords[0] + FR3end - querygaps
        
        return (queryseq[v_end_idx:], v_end_idx)
    
    @staticmethod
    def pruneJregion( alnref, alnrefcoords, alnquery, alnquerycoords, refID, queryseq ):
        """Prune J region out of query sequence based on alignment.
        
        Given ref and query alignments of J region, refID, and the original
        query sequence, return a sequence with the J region cut out, leaving
        the J-TRP.  Also needs query alignment coords.
        
        """
        #DEBUG
        # # check that alnref actually has the whole reference segment
        #       # otherwise, I would need to pass something like alnrefcoords
        #       if alnref.replace('-','') != refseq.IGHJ_seqs[refID]:
        #           raise Exception, "Aligned reference segment is not equal to vdj.refseq reference segment."
        
        FR4start = refseq.IGHJ_offset[refID] - alnrefcoords[0]  # first candidate position of J-TRP start   
        refgaps = alnref[:FR4start].count('-')  # count gaps up to putative TRP pos
        seengaps = 0
        while refgaps != 0:     # iteratively find all gaps up to the TRP
            seengaps += refgaps
            FR4start += refgaps     # adjust for gaps in ref alignment
            refgaps   = alnref[:FR4start].count('-') - seengaps # any add'l gaps?
        
        querygaps = alnquery[:FR4start].count('-')
        
        # v_end_idx = idx of start of aln of query + distance into aln - # of gaps + 3 nt for J-TRP
        j_start_idx = alnquerycoords[0] + FR4start - querygaps + 3
        
        return (queryseq[:j_start_idx],j_start_idx)
    
    @staticmethod
    def construct_alignment(seq1,seq2,scoremat,tracemat):
        """Construct alignment of ref segment to query from score and trace matrices."""
        nrows,ncols = scoremat.shape
        
        # do some error checking
        if len(seq1)+1 != nrows or len(seq2)+1 != ncols:
            raise Exception, "nrows and ncols must be equal to len(seq1)+1 and len(seq2)+1"
        
        #DEBUG
        # if not nrows <= ncols:
        #           raise Exception, "score matrix must have nrows < ncols"
        #       if not len(seq1) <= len(seq2):
        #           raise Exception, "len of seq1 must be smaller than seq2"
        
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
        #DEBUG
        # col = np.argmax( scoremat[nrows-1,:] )
        # row = nrows-1
        
        # if row is coord in matrix, row-1 is coord in seq
        
        aln1 = seq1[row-1]
        aln2 = seq2[col-1]
        
        aln1end = row
        aln2end = col
        
        #DEBUG
        #while row-1 > 0:
        while (row-1 > 0) and (col-1 > 0):
            # compute direction of moves
            rowchange,colchange = deltas[ tracemat[row,col] ]
            
            # WORKS WITH PURE PYTHON alignNW trace return
            #rowchange = row-tracemat[row,col][0]
            #colchange = col-tracemat[row,col][1]
            
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
        
        #DEBUG
        #OLD WAY
        # nrows,ncols = scorematrix.shape
        #       
        #       if not nrows < ncols:
        #           raise Exception, "score matrix must have nrows < ncols"
        #       
        #       return np.max( scorematrix[nrows-1,:] )
        
    @staticmethod
    def scoreDalign(scorematrix):
        """Computes score of D alignment given Smith-Waterman score matrix
        
        """
        return np.max( scorematrix )
    
    @staticmethod
    def dict2sorteddecreasingitemlist(dictionary,keyorvalue='value'):
        pos = {'key':0, 'value':1}
        di = dictionary.items()
        di.sort(key=operator.itemgetter(pos[keyorvalue]))
        di.reverse()
        return di
    
    @staticmethod
    def seqdict2revcompseqdict(seqdict):
        revcompdict = {}
        for item in seqdict.iteritems():
            revcompdict[item[0]] = sequtils.reverse_complement(item[1])
        return revcompdict