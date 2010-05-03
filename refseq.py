






IGHV = []
IGHD = []
IGHJ = []

IGKV = []
IGKJ = []

IGLV = []
IGLJ = []

TRBV = []
TRBD = []
TRBJ = []

TRAV = []
TRAJ = []

TRDV = []
TRDD = []
TRDJ = []

TRGV = []
TRGJ = []





# refseq.py

"""Contains data and functions for dealing with reference IMGT data.

refdatadir must be set to the directory where the IMGT flat file release
is present.  This must include imgt.dat, imgt.fasta, and imgtrefseq.fasta.

IGHn_list -- list of the refseq identifiers (e.g., IGHJ1, IGHV6-1, etc.)
IGHn_acc -- dict where keys are refseq IDs and values are the IMGT accession numbers
IGHn_coords -- dict where keys are refseq IDs and values are pairs of coords that are
                extracted from the refseq database.  These numbers are unmodified seq
                coords of the reference element within the LIGM entry
IGHn_seqs -- dict where the keys are refseq IDs and the values are the actual sequences
                from refseq.
IGHn_idx -- dict where the keys are refseq IDs and the values are the indices into the
                IGHn_list list.
IGHV_FR3_IMGT_end_coord -- dict where the keys are refseq IDs and the values are the coord of
                            the end of FR3-IMGT (incl the 2nd-CYS) unmodified
IGHJ_J_TRP_start_coord -- dict where the keys are refseq IDs and the values are are the
                            coord of the start of the J-TRP site unmodified

"""

import os
import cPickle

from Bio import SeqIO

import vdj
import params
import sequtils

# ===============================================
# = UTILITY FNs for parsing reference databases =
# ===============================================

def get_refseq_elements(locus,alleles,func,species,getcoords,refdatadir,imgtrefseqfasta):
    """Extract all identifiers that meet certain criteria from IMGT/GENE-DB (refseq).
    
    locus -- identifier such as 'IGH', 'TRBV', etc
    alleles -- can be either '01' or 'all'
    func -- list of func IDs, e.g., ['F','ORF','P'].  If empty, includes everything
    species -- species of animal (e.g., 'Homo+sapiens')
    getcoords -- determines whether to process and return coords (boolean)
    
    refdatadir -- directory where IMGT resides
    imgtrefseqfasta -- filename of imgtrefseq fasta file
    
    """
    
    refIDs = []
    refacc = {}
    refseqs = {}
    if getcoords:
        refcoords = {}
    allrefseqs = sequtils.load_fasta( os.path.join(refdatadir,imgtrefseqfasta) )
    
    if func == []:
        func = ['F','ORF','P']
    
    for seq in allrefseqs:
        # get info on curr seq
        currdescr = seq.description
        currlocus = currdescr.split(',')[0].split('*')[0]
        currallele = currdescr.split(',')[0].split('*')[1]
        currspecies = currdescr.split(',')[1]
        currfunc = currdescr.split(',')[2].split()[-1].strip(' ()')
        currID = currdescr.split(',')[0].lstrip('>')
        curracc = currdescr.split(',')[3]
        if getcoords:
            try:
                currcoords = ( eval(currdescr.split(',')[4].split('.')[0]), eval(currdescr.split(',')[4].split('.')[2]) )
            except:
                #print "Coordinates for", currlocus, curracc, "are not interpretable"
                pass
        
        # perform tests to exclude current seq
        if locus not in currlocus:
            continue
        
        if alleles != 'all':
            if currallele != alleles:
                continue
        
        if currspecies != species:
            continue
                
        if currfunc not in func:
            continue
        
        # check for "partial"
        if len(currdescr.split(',')) > 5 and 'partial' in ' '.join(currdescr.split(',')[5:]):
            continue
        
        refIDs.append( currID )
        refacc[currID] = curracc
        refseqs[currID] = seq.seq.tostring().upper()
        if getcoords:
            refcoords[currID] = currcoords
        
    refIDs.append( '' )
    refIDs.sort()
        
    if not getcoords:
        return refIDs,refacc,refseqs
    return refIDs,refacc,refseqs,refcoords

def get_FR3_IMGT_end(IGHV_acc, refdatadir, imgtdat, verbose=False):
    """Get coord of end of FR3 from IMGT LIGM database.
    
    Requires list of IGHV elts and corresponding accession numbers.
    Returns the coord as-is in the IMGT file.  Must be compared with
    unmodified coordinate in IGHV_coords.
    
    NOTE: FR3-IMGT includes the 2nd-CYS
    
    IGHV_acc -- dict where keys are refseq IDs and values are accessions
                where the refseq resides
                
    refdatadir -- directory where IMGT resides
    imgtdat -- filename of IMGT LIGM flat file
    
    """
    
    LIGMflat = open( os.path.join(refdatadir,imgtdat), 'r' )
    
    IGHV_junc_start = {}
    
    numRecordsSeen = 0
    numAccFound = 0
    numCoordsFound = 0
    inTargetRecord = False
    
    for line in LIGMflat:
        splitline = line.split()
        
        if splitline[0] == 'ID':
            AC = splitline[1]
            if AC in IGHV_acc.values():
                inTargetRecord = True
                matchedID = False
                numAccFound += 1
            numRecordsSeen += 1
        elif splitline[0] == 'FT':
            if inTargetRecord:
                if splitline[1].startswith('/allele='):
                    if matchedID == True:
                        if verbose: print "Found", ID, "in", AC, "correctly but failed to find the coords."
                        matchedID = False
                    ID = splitline[1].split('/allele=')[1].strip('"')
                    if ID in IGHV_acc.keys() and IGHV_acc[ID] == AC:
                        matchedID = True
                    else:
                        if verbose:
                            if ID not in IGHV_acc.keys():
                                print "Found", ID, "in", AC, "but it's not a target."
                            elif IGHV_acc[ID] != AC:
                                print "Found", ID, "in", AC, "but I'm supposed to find it elsewhere."
                elif splitline[1] == 'FR3-IMGT':
                    if matchedID == True:
                        if '<' not in splitline[2] and '>' not in splitline[2]:
                            IGHV_junc_start[ID] = eval(splitline[2].split('.')[2])
                            matchedID = False
                            numCoordsFound += 1
                            if verbose: print "Successfully extracted coords for", ID, "in", AC
        elif splitline[0] == 'SQ':
            if inTargetRecord == True:  # if we missed the annotation
                inTargetRecord = False
                if verbose: print "Reached the end of target record", AC
    
    if verbose:
        print "Number of LIGM records read:", numRecordsSeen
        print "Number of target accessions searched:", len(set(IGHV_acc.values()))
        print "Number of target accessions found:", numAccFound
        print "Number of coords searched:", len(IGHV_acc.keys())
        print "Number of coords successfully found:", numCoordsFound
    
    LIGMflat.close()
    
    return IGHV_junc_start

def get_J_TRP_start(IGHJ_acc, refdatadir, imgtdat, verbose=False):
    """Get coord of start of FR4 from IMGT LIGM database.
    
    Requires list of IGHJ elts and corresponding accession numbers.
    Returns the coord as-is in the IMGT file.  Must be compared with
    unmodified coordinate in IGHJ_coords.
    
    NOTE: returns coord of start of conserved J-TRP
    
    IGHJ_acc -- dict where keys are refseq IDs and values are accessions
                where the refseq resides
                
    refdatadir -- directory where IMGT resides
    imgtdat -- filename of IMGT LIGM flat file
    
    """
    
    LIGMflat = open( os.path.join(refdatadir,imgtdat), 'r' )
    
    IGHJ_junc_start = {}
    
    numRecordsSeen = 0
    numAccFound = 0
    numCoordsFound = 0
    inTargetRecord = False
    
    for line in LIGMflat:
        splitline = line.split()
        
        if splitline[0] == 'ID':
            AC = splitline[1]
            if AC in IGHJ_acc.values():
                inTargetRecord = True
                matchedID = False
                numAccFound += 1
            numRecordsSeen += 1
        elif splitline[0] == 'FT':
            if inTargetRecord:
                if splitline[1].startswith('/allele='):
                    if matchedID == True:
                        if verbose: print "Found", ID, "in", AC, "correctly but failed to find the coords."
                        matchedID = False
                    ID = splitline[1].split('/allele=')[1].strip('"')
                    if ID in IGHJ_acc.keys() and IGHJ_acc[ID] == AC:
                        matchedID = True
                    else:
                        if verbose:
                            if ID not in IGHJ_acc.keys():
                                print "Found", ID, "in", AC, "but it's not a target."
                            elif IGHJ_acc[ID] != AC:
                                print "Found", ID, "in", AC, "but I'm supposed to find it elsewhere."
                elif splitline[1] == 'J-TRP':
                    if matchedID == True:
                        if '<' not in splitline[2] and '>' not in splitline[2]:
                            IGHJ_junc_start[ID] = eval(splitline[2].split('.')[0])
                            matchedID = False
                            numCoordsFound += 1
                            if verbose: print "Successfully extracted coords for", ID, "in", AC
        elif splitline[0] == 'SQ':
            if inTargetRecord == True:  # if we missed the annotation
                inTargetRecord = False
                if verbose: print "Reached the end of target record", AC
    
    if verbose:
        print "Number of LIGM records read:", numRecordsSeen
        print "Number of target accessions searched:", len(set(IGHJ_acc.values()))
        print "Number of target accessions found:", numAccFound
        print "Number of coords searched:", len(IGHJ_acc.keys())
        print "Number of coords successfully found:", numCoordsFound
        
    LIGMflat.close()
    
    return IGHJ_junc_start

# ================================
# = Definition of reference data =
# ================================

LOCI = ['IGH','IGK','IGL','TRA','TRB','TRD','TRG']

IGHV_list,IGHV_acc,IGHV_seqs,IGHV_coords = get_refseq_elements(locus='IGHV',alleles='01',func=['F','ORF'],species='Homo+sapiens',getcoords=True,refdatadir=params.refdatadir,imgtrefseqfasta=params.imgtrefseqfasta)
IGHD_list,IGHD_acc,IGHD_seqs = get_refseq_elements(locus='IGHD',alleles='01',func=['F','ORF'],species='Homo+sapiens',getcoords=False,refdatadir=params.refdatadir,imgtrefseqfasta=params.imgtrefseqfasta)
IGHJ_list,IGHJ_acc,IGHJ_seqs,IGHJ_coords = get_refseq_elements(locus='IGHJ',alleles='01',func=['F','ORF'],species='Homo+sapiens',getcoords=True,refdatadir=params.refdatadir,imgtrefseqfasta=params.imgtrefseqfasta)

IGHV_idx = dict([(g,i) for i,g in enumerate(IGHV_list)])
IGHD_idx = dict([(g,i) for i,g in enumerate(IGHD_list)])
IGHJ_idx = dict([(g,i) for i,g in enumerate(IGHJ_list)])

IGHC_list = [
            '',
            'IGHA1',
            'IGHA2',
            'IGHD',
            'IGHE',
            'IGHG1',
            'IGHG2',
            'IGHG3',
            'IGHG4',
            'IGHM'
            ]
IGHC_idx = dict([(g,i) for i,g in enumerate(IGHC_list)])

ALL_IDs = list(set(IGHV_list + IGHD_list + IGHJ_list + IGHC_list))
ALL_IDs.sort()

imgtrefseqFR3endcoords = 'imgtrefseqFR3endcoords.dat'
imgtrefseqJTRPstartcoords = 'imgtrefseqJTRPstartcoords.dat'

if os.path.exists(os.path.join(params.refdatadir,params.imgtrefseqFR3endcoords)):
    pickle_in = open(os.path.join(params.refdatadir,params.imgtrefseqFR3endcoords),'r')
    IGHV_FR3_IMGT_end_coord = cPickle.load(pickle_in)
    pickle_in.close()
else:
    IGHV_FR3_IMGT_end_coord = get_FR3_IMGT_end(IGHV_acc=IGHV_acc, refdatadir=params.refdatadir, imgtdat=params.imgtdat, verbose=False)
    pickle_out = open(os.path.join(params.refdatadir,params.imgtrefseqFR3endcoords),'w')
    cPickle.dump(IGHV_FR3_IMGT_end_coord,pickle_out)
    pickle_out.close()

if os.path.exists(os.path.join(params.refdatadir,params.imgtrefseqJTRPstartcoords)):
    pickle_in = open(os.path.join(params.refdatadir,params.imgtrefseqJTRPstartcoords),'r')
    IGHJ_J_TRP_start_coord = cPickle.load(pickle_in)
    pickle_in.close()
else:
    IGHJ_J_TRP_start_coord  = get_J_TRP_start(IGHJ_acc=IGHJ_acc, refdatadir=params.refdatadir, imgtdat=params.imgtdat, verbose=False)
    pickle_out = open(os.path.join(params.refdatadir,params.imgtrefseqJTRPstartcoords),'w')
    cPickle.dump(IGHJ_J_TRP_start_coord,pickle_out)
    pickle_out.close()

IGHV_offset = {}
IGHJ_offset = {}
for seg in IGHV_list[1:]:
    # offset to end before the 2nd-CYS
    IGHV_offset[seg] = IGHV_FR3_IMGT_end_coord[seg] - IGHV_coords[seg][0] - 2
for seg in IGHJ_list[1:]:
    # offset to start before the J-TRP
    IGHJ_offset[seg] = IGHJ_J_TRP_start_coord[seg]  - IGHJ_coords[seg][0]

# ====================
# = Specificity data =
# ====================

def get_LIGM_with_specificities(refdatadir,imgtdat,imgtfasta,outputfasta,outputvdjxml):
    '''
    Take the IMGT LIGM flat and fasta files, and return a fasta with only those seqs
    that have a specificity associated with them in a fasta file.  The header contains
    the accession and the specificity.
    
    vdj.refseq.get_LIGM_with_specificities("imgt.dat","imgt.fasta","imgt.specificities.fasta","imgt.specificities.vdjxml")
    '''
    specificities = {}      # list of 2-tuples: (ACCESION,specificity)
    LIGMflat = open(os.path.join(refdatadir,imgtdat),'r')
    LIGMfasta = open(os.path.join(refdatadir,imgtfasta),'r')
    opSpecificityfasta = open(os.path.join(refdatadir,outputfasta),'w')
    
    numRecords = 0
    numRecordswithSpec = 0
    
    ID = ''
    DE = ''
    for line in LIGMflat:
        splitline = line.split()
        
        if splitline[0] == 'ID':
            inRecord = True
            ID = splitline[1]
            numRecords += 1
        elif splitline[0] == 'DE':
            if inRecord:
                DE += ' '.join(splitline[1:]) + ' '
        elif splitline[0] == 'XX':
            if DE == '':    # if i haven't stored the description yet
                continue
            else:   # finished record
                if 'specificity' in DE:
                    specidx = DE.rfind('specificity')
                    spec = DE[specidx+len('specificity'):].strip().rstrip('.')
                    if spec.startswith('anti'):
                        specificities[ID] = spec
                        numRecordswithSpec += 1
                ID = ''
                DE = ''
                inRecord = False
    
    print "Number of LIGM records read: " + str(numRecords)
    print "Number of LIGM records that have specificities: " + str(numRecordswithSpec)
    
    numFasta = 0
    
    for seq in SeqIO.parse(LIGMfasta,'fasta'):
        spec = specificities.get(seq.id,'')
        if spec == '':
            continue
        else:
            print >>opSpecificityfasta, ">" + seq.id + " | " + spec
            print >>opSpecificityfasta, seq.seq.tostring().upper()
            numFasta += 1
    
    print "Number of Fasta records with specificities found and printed: " + str(numFasta)
    
    LIGMflat.close()
    LIGMfasta.close()
    opSpecificityfasta.close()
    
    # VDJXML and alignment
    rep = vdj.initial_import([os.path.join(refdatadir,outputfasta)],os.path.join(refdatadir,outputvdjxml),metatags=['specificity_reference : '+vdj.timestamp()])
    rep = vdj.positive_strand(rep)
    rep = vdj.align_rep(rep)
    
    #repfiltered = rep.get_chains_fullVJ()
    repfiltered = rep
    for chain in repfiltered:
        spec = specificities.get(chain.descr,'')
        if spec == '':
            print "Reference chain has empty specificity: " + chain.descr
            continue
        else:
            chain.add_tags( 'specificity|'+spec )
    
    vdj.writeVDJ(repfiltered,os.path.join(refdatadir,outputvdjxml)) 
    
    return


if not os.path.exists(os.path.join(params.refdatadir,params.imgtspecfasta)) or not os.path.exists(os.path.join(params.refdatadir,params.imgtspecvdjxml)):
    get_LIGM_with_specificities(params.refdatadir,params.imgtdat,params.imgtfasta,params.imgtspecfasta,params.imgtspecvdjxml)


if os.path.exists(os.path.join(params.refdatadir,params.imgtspecfasta)) and os.path.exists(os.path.join(params.refdatadir,params.imgtspecvdjxml)):
    ipspecfasta = open(os.path.join(params.refdatadir,params.imgtspecfasta),'r')
    
    SPEC_list = set()
    for line in ipspecfasta:
        if line[0] == '>':
            currspec = '|'.join(line.split('|')[1:]).strip()
            SPEC_list.add(currspec)
    ipspecfasta.close()
    SPEC_list.add('')
    SPEC_list = list(SPEC_list)
    SPEC_list.sort()
