import warnings
import copy

import numpy as np
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import seqtools

import vdj
import refseq
import seqtools
import alignmentcore

warnings.simplefilter('always')

class vdj_aligner(object):
    
    def __init__(self,**kw):
        
        self.numCrudeVCandidates = 5
        self.numCrudeDCandidates = 10
        self.numCrudeJCandidates = 2
        
        self.minVscore = 100    # derived from calibration data 20090710
        self.minDscore = 4
        self.minJscore = 13
        
        if kw.has_key('rigorous') and kw['rigorous'] == True:
            self.numCrudeVCandidates = 10000
            self.numCrudeDCandidates = 10000
            self.numCrudeJCandidates = 10000

            self.minVscore = 20
            self.minDscore = 1
            self.minJscore = 5
        
        # Define seed patterns
        patternA='111011001011010111'
        patternB='1111000100010011010111'
        patternC='111111111111'
        patternD='110100001100010101111'
        patternE='1110111010001111'
        self.seedpatterns = [patternA,patternB,patternC,patternD,patternE]
        self.miniseedpatterns = ['111011','110111']
        self.patternPos = '111111111111'
        
        # set reference sequences (locus) and generate hashes from ref data
        self.locus = kw['locus']
        
        self.refV = refseq.__getattribute__(self.locus+'V')
        refV_seqs = dict([(allele,record.seq.tostring()) for (allele,record) in self.refV.iteritems()])
        self.Vseqlistkeys = vdj_aligner.seqdict2kmers( refV_seqs, self.seedpatterns )
        
        self.refJ = refseq.__getattribute__(self.locus+'J')
        refJ_seqs = dict([(allele,record.seq.tostring()) for (allele,record) in self.refJ.iteritems()])
        self.Jseqlistkeys = vdj_aligner.seqdict2kmers( refJ_seqs, self.seedpatterns )
        
        try:    # this locus may not have D segments
            self.refD = refseq.__getattribute__(self.locus+'D')
            refD_seqs = dict([(allele,record.seq.tostring()) for (allele,record) in self.refD.iteritems()])
            self.Dseqlistkeysmini = vdj_aligner.seqdict2kmers( refD_seqs, self.miniseedpatterns )
        except AttributeError:
            pass
        
        # Generate reference data for positive sequence ID
        posVseqlistkeys = vdj_aligner.seqdict2kmers( refV_seqs, [self.patternPos] )
        posJseqlistkeys = vdj_aligner.seqdict2kmers( refJ_seqs, [self.patternPos] )
        negVseqlistkeys = vdj_aligner.seqdict2kmers( vdj_aligner.seqdict2revcompseqdict(refV_seqs), [self.patternPos] )
        negJseqlistkeys = vdj_aligner.seqdict2kmers( vdj_aligner.seqdict2revcompseqdict(refJ_seqs), [self.patternPos] )
        
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
        # compute hashes from query seq
        querykeys = vdj_aligner.seq2kmers(chain.seq.tostring(),self.seedpatterns)
        
        # for each reference V segment and each pattern, how many shared k-mers are there?
        Vscores_hash = vdj_aligner.hashscore(self.Vseqlistkeys,querykeys)
        
        # get numCrudeVCandidates highest scores in Vscores and store their names in descending order
        goodVseglist = sorted(self.refV.keys(),key=lambda k: Vscores_hash[k],reverse=True)[0:self.numCrudeVCandidates]
        goodVsegdict = dict([(seg,self.refV[seg].seq.tostring()) for seg in goodVseglist])
        
        # Needleman-Wunsch of V segment
        (bestVseg,bestVscore,bestVscoremat,bestVtracemat) = vdj_aligner.bestalignNW(goodVsegdict,chain.seq.tostring(),self.minVscore)
        
        # if successful alignment
        if bestVseg is not None:
            # copy features from ref to query
            Vrefaln,Vqueryaln = vdj_aligner.construct_alignment( self.refV[bestVseg].seq.tostring(), chain.seq.tostring(), bestVscoremat, bestVtracemat )
            coord_mapping = vdj_aligner.ungapped_coord_mapping(Vrefaln, Vqueryaln)
            seqtools.copy_features(self.refV[bestVseg], chain, coord_mapping, erase=['translation'], replace=False)
            
            # store gapped aln
            chain.annotations['gapped_query'] = Vqueryaln
            chain.annotations['gapped_reference'] = Vrefaln
            
            # annotate mutations
            curr_annot = chain.letter_annotations['alignment']
            aln_annot = vdj_aligner.alignment_annotation(Vrefaln,Vqueryaln)
            aln_annot = aln_annot.translate(None,'D')
            lNER = len(aln_annot) - len(aln_annot.lstrip('I'))
            rNER = len(aln_annot.rstrip('I'))
            chain.letter_annotations['alignment'] = curr_annot[:lNER] + aln_annot[lNER:rNER] + curr_annot[rNER:]
            
            # perform some curating; esp, CDR3-IMGT is annotated in V
            # references, though it's not complete. I will recreate that
            # annotation manually.
            chain._update_feature_dict()
            try:    # some reference entries do not have CDR3 annotations
                chain.features.pop(chain._features['CDR3-IMGT'][0])
                chain._features.pop('CDR3-IMGT')
                chain._update_feature_dict()
            except KeyError:
                pass
            
            # update codon_start of V-REGION anchored to the CDR3 2nd-CYS
            cys = chain.features[ chain._features['2nd-CYS'][0] ]
            v_reg = chain.features[ chain._features['V-REGION'][0] ]
            v_reg.qualifiers['codon_start'] = [cys.location.start.position % 3 + 1]
        
        return bestVscore
    
    
    def Jalign_chain(self,chain,verbose=False):
        # try pruning off V region for J alignment
        try:
            second_cys = chain.__getattribute__('2nd-CYS')
            second_cys_offset = second_cys.location.end.position
            query = chain.seq.tostring()[second_cys_offset:]
        except AttributeError:
            query = chain.seq.tostring()
            second_cys_offset = 0
        
        # compute hashes from query seq
        querykeys = vdj_aligner.seq2kmers(query,self.seedpatterns)
        
        # for each reference J segment and each pattern, how many shared k-mers are there?
        Jscores_hash = vdj_aligner.hashscore(self.Jseqlistkeys,querykeys)
        
        # get numCrudeJCandidates highest scores in Jscores and store their names in descending order
        goodJseglist = sorted(self.refJ.keys(),key=lambda k: Jscores_hash[k],reverse=True)[0:self.numCrudeJCandidates]
        goodJsegdict = dict([(seg,self.refJ[seg].seq.tostring()) for seg in goodJseglist])
        
        # Needleman-Wunsch of J segment
        (bestJseg,bestJscore,bestJscoremat,bestJtracemat) = vdj_aligner.bestalignNW(goodJsegdict,query,self.minJscore)
        
        # if successful alignment
        if bestJseg is not None:
            # copy features from ref to query
            Jrefaln,Jqueryaln = vdj_aligner.construct_alignment( self.refJ[bestJseg].seq.tostring(), query, bestJscoremat, bestJtracemat )
            coord_mapping = vdj_aligner.ungapped_coord_mapping(Jrefaln, Jqueryaln)
            seqtools.copy_features(self.refJ[bestJseg], chain, coord_mapping, offset=second_cys_offset, erase=['translation'], replace=False)
            chain._update_feature_dict()
            
            # update gapped aln
            gapped_query = chain.annotations.get('gapped_query','')
            gapped_reference = chain.annotations.get('gapped_reference','')
            gapped_CDR3_offset = vdj_aligner.ungapped2gapped_coord(chain.seq.tostring(),gapped_query,second_cys_offset)
            gapped_Vref_aln_end = len(gapped_reference.rstrip('-'))
            chain.annotations['gapped_query'] = gapped_query[:gapped_Vref_aln_end] + Jqueryaln[gapped_Vref_aln_end-gapped_CDR3_offset:]
            chain.annotations['gapped_reference'] = gapped_reference[:gapped_Vref_aln_end] + Jrefaln[gapped_Vref_aln_end-gapped_CDR3_offset:]
            
            # annotate mutations
            curr_annot = chain.letter_annotations['alignment']
            aln_annot = vdj_aligner.alignment_annotation(Jrefaln,Jqueryaln)
            aln_annot = aln_annot.translate(None,'D')
            lNER = len(aln_annot) - len(aln_annot.lstrip('I'))
            rNER = len(aln_annot.rstrip('I'))
            chain.letter_annotations['alignment'] = curr_annot[:second_cys_offset+lNER] + aln_annot[lNER:rNER] + curr_annot[second_cys_offset+rNER:]
        
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
        goodDseglist = sorted(self.refD.keys(),key=lambda k: Dscores_hash[k],reverse=True)[0:self.numCrudeDCandidates]
        goodDsegdict = dict([(seg,self.refD[seg].seq.tostring()) for seg in goodDseglist])
        
        # Needleman-Wunsch of J segment
        (bestDseg,bestDscore,bestDscoremat,bestDtracemat) = vdj_aligner.bestalignSW(goodDsegdict,query,self.minDscore)
        
        # if successful alignment
        if bestDseg is not None:
            # TEMPORARY SOLUTION
            chain.annotations['D-REGION'] = bestDseg
        
        return bestDscore
    
    
    def align_chain(self,chain,verbose=False,debug=False):
        
        # DEBUG
        # import vdj
        # import vdj.alignment
        # from Bio import SeqIO
        # from Bio.Alphabet import generic_dna
        # iter = SeqIO.parse('smallset.fasta','fasta',generic_dna)
        # iter = SeqIO.parse('donor12_cd8_memory_raw_reads.fasta','fasta',generic_dna)
        # aligner = vdj.alignment.igh_aligner()
        # aligner = vdj.alignment.trb_aligner()
        # a = iter.next()
        # a = vdj.ImmuneChain(a)
        # aligner.coding_chain(a)
        # aligner.align_chain(a)
        # print a
        # 
        if debug:
            import pdb
            pdb.set_trace()
        
        if chain.seq.tostring() != chain.seq.tostring().upper():
            raise ValueError, "aligner requires all uppercase alphabet."
        
        if not chain.has_tag('positive') and not chain.has_tag('coding'):
            warnings.warn('chain %s may not be the correct strand' % chain.id)
        
        # insert letter annotations for alignment annotation
        chain.letter_annotations["alignment"] = '_' * len(chain)
        
        scores = {}
        
        scores['v'] = self.Valign_chain(chain,verbose)
        
        scores['j'] = self.Jalign_chain(chain,verbose)
        
        # manually annotate CD3-IMGT, only if V and J alns are successful
        try:
            if chain.v and chain.j:
                cdr3_start = chain.__getattribute__('2nd-CYS').location.end.position
                try:
                    cdr3_end = chain.__getattribute__('J-PHE').location.start.position
                except AttributeError:
                    cdr3_end = chain.__getattribute__('J-TRP').location.start.position
                cdr3_feature = SeqFeature(location=FeatureLocation(cdr3_start,cdr3_end),type='CDR3-IMGT',strand=1)
                chain.features.append(cdr3_feature)
                chain._update_feature_dict()
                
                # erase alignment annotations in CDR3.  can't tell SHM from TdT at this point
                curr_annot = chain.letter_annotations['alignment']
                chain.letter_annotations['alignment'] = curr_annot[:cdr3_start] + '3' * (cdr3_end-cdr3_start) + curr_annot[cdr3_end:]
                
                # if I am in a locus with D segments, try aligning that as well
                if self.locus in ['IGH','TRB','TRD']:
                    scores['d'] = self.Dalign_chain(chain,verbose)
        except AttributeError:    # chain.v or chain.j raised an error
            pass
        
        return scores
    
    
    def coding_chain(self,chain,verbose=False):
        strand = self.seq2coding(chain.seq.tostring())
        if strand == -1:
            chain.seq = chain.seq.reverse_complement()
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
    def alignment_annotation(aln_ref,aln_query):
        # should be given equivalenced region
        assert len(aln_query) == len(aln_ref)
        annot = ''
        for (ref_letter,query_letter) in zip(aln_ref,aln_query):
            if query_letter == '-':
                annot += 'D'
            elif ref_letter == '-':
                annot += 'I'
            elif query_letter == ref_letter:
                annot += '.'
            else:
                annot += 'S'
        
        return annot
    
    
    @staticmethod
    def ungapped_coord_mapping(aln_from, aln_to):
        if len(aln_from) != len(aln_to):
            raise ValueError, "from and to strings must be same length"
        
        coord_from = 0
        coord_to = 0
        mapping = {}
        for coord_gapped in range(len(aln_from)):
            
            if aln_from[coord_gapped-1:coord_gapped+1] == '--':
                coord_to += 1
                continue
            
            mapping.setdefault(coord_from,[]).append(coord_to)
            
            if aln_from[coord_gapped] != '-':
                coord_from += 1
            
            if aln_to[coord_gapped] != '-':
                coord_to += 1
        
        mapping.setdefault(coord_from,[]).append(coord_to)
        
        return mapping
    
    
    @staticmethod
    def ungapped2gapped_coord(ungapped,gapped,ungapped_coord):
        left_gaps = len(gapped) - len(gapped.lstrip('-'))
        gapped_coord = ungapped_coord + left_gaps
        gaps = gapped.count('-',0,gapped_coord)
        while gapped_coord - gaps < ungapped_coord:
            gapped_coord += gaps
            gaps = gapped.count('-',0,gapped_coord)
        return gapped_coord
    
    
    @staticmethod
    def construct_alignment(seq1,seq2,scoremat,tracemat):
        """Construct alignment of ref segment to query from score and trace
        matrices.
        """
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
        # if nrows <= ncols:
        #     col = np.argmax( scoremat[nrows-1,:] )
        #     row = nrows-1
        # else:
        #     col = ncols-1
        #     row = np.argmax( scoremat[:,ncols-1] )
        col = np.argmax( scoremat[nrows-1,:] )
        row = nrows-1
        
        # if row is coord in matrix, row-1 is coord in seq (b/c of init conditions)
        aln1 = seq1[row-1:] + '-'*(ncols-col-1)
        aln2 = seq2[col-1:] + '-'*(nrows-row-1)
        
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
        
        aln1 = seq1[:row-1]+ '-'*(col-1) + aln1
        aln2 = seq2[:col-1]+ '-'*(row-1) + aln2
        
        return aln1, aln2
    
    
    @staticmethod
    def scoreVJalign(scorematrix):
        """Computes score of V alignment given Needleman-Wunsch score matrix
        
        ASSUMES num rows < num cols, i.e., refseq V seg is on vertical axis
        
        """
        nrows,ncols = scorematrix.shape
        
        # if nrows <= ncols:
        #     return np.max( scorematrix[nrows-1,:] )
        # else:
        #     return np.max( scorematrix[:,ncols-1] )
        return np.max( scorematrix[nrows-1,:] )
    
    
    @staticmethod
    def scoreDalign(scorematrix):
        """Computes score of D alignment given Smith-Waterman score matrix
        
        """
        return np.max( scorematrix )
    
    
    @staticmethod
    def seqdict2revcompseqdict(seqdict):
        revcompdict = {}
        for item in seqdict.iteritems():
            revcompdict[item[0]] = seqtools.reverse_complement(item[1])
        return revcompdict


class vdj_aligner_combined(object):
    """vdj aligner for 'light' chains
    
    this class will perform alignment for both loci, e.g., IGK and IGL
    and pick the one with the better V score
    """
    def __init__(self,**kw):
        self.loci = kw['loci']
        self.aligners = [vdj_aligner(locus=locus,**kw) for locus in self.loci]
        
        self.patternPos = '111111111111'
        self.posset = set()
        self.negset = set()
        for aligner in self.aligners:
            self.posset.update(aligner.posset)
            self.negset.update(aligner.negset)
    
    def align_chain(self,chain,verbose=False,debug=False):
        alignments = []
        for aligner in self.aligners:
            curr_chain = copy.deepcopy(chain)
            curr_score = aligner.align_chain(curr_chain,debug=debug)
            alignments.append((curr_chain,curr_score))
        alignments = sorted(filter(lambda a: hasattr(a[0],'v'),alignments),key=lambda a:a[1]['v'],reverse=True)
        if len(alignments) > 0:
            bestchain = alignments[0][0]
            chain.__init__(bestchain)
            return alignments[0][1]     # NOTE: I only return the scores upon successful aln
    
    def coding_chain(self,chain,verbose=False):
        strand = self.seq2coding(chain.seq.tostring())
        if strand == -1:
            chain.seq = chain.seq.reverse_complement()
            chain.add_tag('revcomp')
        chain.add_tag('coding')
    
    def seq2coding(self,seq):        
        seqkeys = vdj_aligner.seq2kmers(seq,[self.patternPos])
        seqwords = seqkeys[self.patternPos]
        strandid = 1
        if len(self.negset & seqwords) > len(self.posset & seqwords):
            strandid = -1
        return strandid

def igh_aligner(**kw):
    return vdj_aligner(locus='IGH',**kw)

def igk_aligner(**kw):
    return vdj_aligner(locus='IGK',**kw)

def igl_aligner(**kw):
    return vdj_aligner(locus='IGL',**kw)

def igkl_aligner(**kw):
    return vdj_aligner_combined(loci=['IGK','IGL'],**kw)

def trb_aligner(**kw):
    return vdj_aligner(locus='TRB',**kw)

def tra_aligner(**kw):
    return vdj_aligner(locus='TRA',**kw)

def trd_aligner(**kw):
    return vdj_aligner(locus='TRD',**kw)

def trg_aligner(**kw):
    return vdj_aligner(locus='TRG',**kw)
