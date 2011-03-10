import os
import sys
import types
import cPickle as pickle

import params
# import refsequtils

import vdj

# There are two ways to initialize the reference database:
# 
#     1.  Reference IMGT set (used for IMGT/V-QUEST)
#         
#         The set should be a FASTA file from the V-QUEST set, or formatted
#         similarly.
#     
#     2.  From your own set.  The data should be formatted as IMGT-flavored
#         INSDC (i.e., Genbank/EMBL).
#         
#         The set will trivially just transfer any annotation in the
#         reference sequences.  This should, minimally, include the IMGT
#         annotations for the CDR3 (2nd-CYS for V region and J-TRP or
#         J-PHE for J region).
# 
# The package will load the IMGT set by default.  The first time it loads the
# IMGT set it will save an IMGT-INSDC-type file for each locus.
# 
# If you want to load your own set, you can than do so afterwards and it will
# override the IMGT set.  These are premade.

# ==================
# = Load IMGT data =
# ==================

# check if we already processed the IMGT reference data before

ref_dir_full_path = os.path.join(params.vdj_dir,params.data_dir)
ligm_index = None   # LIGM not loaded by default

# figure out which files to process

processed_dir_full_path = os.path.join(params.vdj_dir,params.processed_dir,vdj.organism)
if not os.path.exists(processed_dir_full_path):
    os.mkdir(processed_dir_full_path)

reference_files = glob.glob( os.path.join(ref_dir_full_path,'*.fasta') )
processed_files = glob.glob( os.path.join(processed_dir_full_path,'*.imgt') )

get_group = lambda f: os.path.splitext(os.path.basename(f))[0]
loci_to_process = set(map(get_group,reference_files)) - set(map(get_group,processed_files))
files_to_process = filter(lambda f: get_group(f) in loci_to_process,reference_files)

# process those files

for reference_fasta in files_to_process:
    if ligm_index == None: ligm_index = SeqIO.index( os.path.join(params.imgt_dir,ligm_filename), 'imgt')
    group = get_group(reference_file)
    reference_imgt = os.path.join(processed_dir_full_path,group+'.imgt')
    
    reference_records = []
    for record in SeqIO.parse(reference_fasta,'fasta'):
        header_data = parse_imgt_fasta_header(record)
        imgt_record = ligm_index[header_data['accession']]
        
        # prune down to the V-REGION annotated in V-QUEST fasta file
        reference_record = imgt_record[header_data['coords'][0]:header_data['coords'][1]]
        reference_records.append(reference_record)

def parse_imgt_fasta_header(record):
    raw_data = record.description.strip().split('|')
    data = {}
    data['accession'] = raw_data[0]
    data['allele'] = raw_data[1]
    data['group'] = data['allele'][0:4]
    data['gene'] = data['allele'].split('*')[0]
    data['species'] = raw_data[2]
    data['functionality'] = raw_data[3]
    data['imgt_label'] = raw_data[4]
    data['frame'] = int(raw_data[7]) - 1   # note change to python numbering (0-based)
    data['partial'] = raw_data[13]
    raw_coords_start = int(raw_data[5].split('.')[0])
    raw_coords_end   = int(raw_data[5].split('.')[-1])
    coords = (raw_coords_start - 1,raw_coords_end)    # note the conversion to python coord system
    data['coords'] = coords
    return data


group = os.splitext(os.path.basename(reference_file))[0]



from Bio import SeqIO

            
        
        

def process_imgt_reference_fasta( ref_fasta ):
    pass

def process_fasta_reference_dir( ref_dir ):
    pass



# ==============================
# = First-time initializations =
# ==============================

# does the pickle directory exist?
if not os.path.exists( os.path.join(params.vdj_dir,params.pickle_dir) ):
    os.mkdir( os.path.join(params.vdj_dir,params.pickle_dir) )

# test for each gene type pickle file
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.IGHV_pickle)): refsequtils.process_IMGT_references(refsequtils.VReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.IGHV_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.IGHV_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.IGHD_pickle)): refsequtils.process_IMGT_references( refsequtils.ReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.IGHD_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.IGHD_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.IGHJ_pickle)): refsequtils.process_IMGT_references(refsequtils.JReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.IGHJ_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.IGHJ_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.IGKV_pickle)): refsequtils.process_IMGT_references(refsequtils.VReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.IGKV_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.IGKV_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.IGKJ_pickle)): refsequtils.process_IMGT_references(refsequtils.JReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.IGKJ_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.IGKJ_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.IGLV_pickle)): refsequtils.process_IMGT_references(refsequtils.VReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.IGLV_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.IGLV_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.IGLJ_pickle)): refsequtils.process_IMGT_references(refsequtils.JReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.IGLJ_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.IGLJ_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.TRBV_pickle)): refsequtils.process_IMGT_references(refsequtils.VReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.TRBV_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.TRBV_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.TRBD_pickle)): refsequtils.process_IMGT_references( refsequtils.ReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.TRBD_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.TRBD_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.TRBJ_pickle)): refsequtils.process_IMGT_references(refsequtils.JReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.TRBJ_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.TRBJ_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.TRAV_pickle)): refsequtils.process_IMGT_references(refsequtils.VReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.TRAV_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.TRAV_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.TRAJ_pickle)): refsequtils.process_IMGT_references(refsequtils.JReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.TRAJ_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.TRAJ_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.TRDV_pickle)): refsequtils.process_IMGT_references(refsequtils.VReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.TRDV_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.TRDV_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.TRDD_pickle)): refsequtils.process_IMGT_references( refsequtils.ReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.TRDD_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.TRDD_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.TRDJ_pickle)): refsequtils.process_IMGT_references(refsequtils.JReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.TRDJ_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.TRDJ_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.TRGV_pickle)): refsequtils.process_IMGT_references(refsequtils.VReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.TRGV_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.TRGV_pickle),verbose=True)
if not os.path.exists(os.path.join(params.vdj_dir,params.pickle_dir,params.TRGJ_pickle)): refsequtils.process_IMGT_references(refsequtils.JReferenceEntry,os.path.join(params.vdj_dir,params.data_dir,params.TRGJ_fasta),os.path.join(params.vdj_dir,params.pickle_dir,params.TRGJ_pickle),verbose=True)

# at this point, there should be pickle files with fully processed reference genes
# (including the LIGM-dependent parts)


# =======================
# = Load reference data =
# =======================

IGHV = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.IGHV_pickle)))
IGHD = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.IGHD_pickle)))
IGHJ = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.IGHJ_pickle)))
IGKV = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.IGKV_pickle)))
IGKJ = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.IGKJ_pickle)))
IGLV = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.IGLV_pickle)))
IGLJ = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.IGLJ_pickle)))
# TRBV = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.TRBV_pickle)))
# TRBD = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.TRBD_pickle)))
# TRBJ = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.TRBJ_pickle)))
# TRAV = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.TRAV_pickle)))
# TRAJ = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.TRAJ_pickle)))
# TRDV = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.TRDV_pickle)))
# TRDD = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.TRDD_pickle)))
# TRDJ = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.TRDJ_pickle)))
# TRGV = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.TRGV_pickle)))
# TRGJ = pickle.load(open(os.path.join(params.vdj_dir,params.pickle_dir,params.TRGJ_pickle)))


# =========================================================
# = Define data structures for compatibility with aligner =
# =========================================================

# IGHn_list -- list of the refseq identifiers (e.g., IGHJ1, IGHV6-1, etc.)
# IGHn_seqs -- dict where the keys are refseq IDs and the values are the actual sequences
#                 from refseq.
# IGHn_idx -- dict where the keys are refseq IDs and the values are the indices into the
#                 IGHn_list list.
# IGHV_offset

def legacy_data(ref_data):
    """Generate data structure used by aln tools etc."""
    locus_list = [elt.allele for elt in ref_data]
    locus_seqs = dict([(elt.allele,elt.seq) for elt in ref_data])
    locus_idx  = dict([(g,i) for i,g in enumerate(locus_list)])
    try:
        locus_offset = dict([(elt.allele,elt.cdr3_boundary-elt.accession_coords[0]) for elt in ref_data])
        return (locus_list,locus_seqs,locus_idx,locus_offset)
    except AttributeError:
        return (locus_list,locus_seqs,locus_idx)

(IGHV_list,IGHV_seqs,IGHV_idx,IGHV_offset) = legacy_data(IGHV)
(IGHD_list,IGHD_seqs,IGHD_idx)             = legacy_data(IGHD)
(IGHJ_list,IGHJ_seqs,IGHJ_idx,IGHJ_offset) = legacy_data(IGHJ)
(IGKV_list,IGKV_seqs,IGKV_idx,IGKV_offset) = legacy_data(IGKV)
(IGKJ_list,IGKJ_seqs,IGKJ_idx,IGKJ_offset) = legacy_data(IGKJ)
(IGLV_list,IGLV_seqs,IGLV_idx,IGLV_offset) = legacy_data(IGLV)
(IGLJ_list,IGLJ_seqs,IGLJ_idx,IGLJ_offset) = legacy_data(IGLJ)
# (TRBV_list,TRBV_seqs,TRBV_idx,TRBV_offset) = legacy_data(TRBV)
# (TRBD_list,TRBD_seqs,TRBD_idx)             = legacy_data(TRBD)
# (TRBJ_list,TRBJ_seqs,TRBJ_idx,TRBJ_offset) = legacy_data(TRBJ)
# (TRAV_list,TRAV_seqs,TRAV_idx,TRAV_offset) = legacy_data(TRAV)
# (TRAJ_list,TRAJ_seqs,TRAJ_idx,TRAJ_offset) = legacy_data(TRAJ)
# (TRDV_list,TRDV_seqs,TRDV_idx,TRDV_offset) = legacy_data(TRDV)
# (TRDD_list,TRDD_seqs,TRDD_idx)             = legacy_data(TRDD)
# (TRDJ_list,TRDJ_seqs,TRDJ_idx,TRDJ_offset) = legacy_data(TRDJ)
# (TRGV_list,TRGV_seqs,TRGV_idx,TRGV_offset) = legacy_data(TRGV)
# (TRGJ_list,TRGJ_seqs,TRGJ_idx,TRGJ_offset) = legacy_data(TRGJ)





# # refseq.py
# 
# """Contains data and functions for dealing with reference IMGT data.
# 
# refdatadir must be set to the directory where the IMGT flat file release
# is present.  This must include imgt.dat, imgt.fasta, and imgtrefseq.fasta.
# 
# IGHn_list -- list of the refseq identifiers (e.g., IGHJ1, IGHV6-1, etc.)
# IGHn_acc -- dict where keys are refseq IDs and values are the IMGT accession numbers
# IGHn_coords -- dict where keys are refseq IDs and values are pairs of coords that are
#                 extracted from the refseq database.  These numbers are unmodified seq
#                 coords of the reference element within the LIGM entry
# IGHn_seqs -- dict where the keys are refseq IDs and the values are the actual sequences
#                 from refseq.
# IGHn_idx -- dict where the keys are refseq IDs and the values are the indices into the
#                 IGHn_list list.
# IGHV_FR3_IMGT_end_coord -- dict where the keys are refseq IDs and the values are the coord of
#                             the end of FR3-IMGT (incl the 2nd-CYS) unmodified
# IGHJ_J_TRP_start_coord -- dict where the keys are refseq IDs and the values are are the
#                             coord of the start of the J-TRP site unmodified
# 
# """
# 
# import os
# import cPickle
# 
# from Bio import SeqIO
# 
# import vdj
# import params
# import sequtils
# 
# # ===============================================
# # = UTILITY FNs for parsing reference databases =
# # ===============================================
# 
# def get_refseq_elements(locus,alleles,func,species,getcoords,refdatadir,imgtrefseqfasta):
#     """Extract all identifiers that meet certain criteria from IMGT/GENE-DB (refseq).
#     
#     locus -- identifier such as 'IGH', 'TRBV', etc
#     alleles -- can be either '01' or 'all'
#     func -- list of func IDs, e.g., ['F','ORF','P'].  If empty, includes everything
#     species -- species of animal (e.g., 'Homo+sapiens')
#     getcoords -- determines whether to process and return coords (boolean)
#     
#     refdatadir -- directory where IMGT resides
#     imgtrefseqfasta -- filename of imgtrefseq fasta file
#     
#     """
#     
#     refIDs = []
#     refacc = {}
#     refseqs = {}
#     if getcoords:
#         refcoords = {}
#     allrefseqs = sequtils.load_fasta( os.path.join(refdatadir,imgtrefseqfasta) )
#     
#     if func == []:
#         func = ['F','ORF','P']
#     
#     for seq in allrefseqs:
#         # get info on curr seq
#         currdescr = seq.description
#         currlocus = currdescr.split(',')[0].split('*')[0]
#         currallele = currdescr.split(',')[0].split('*')[1]
#         currspecies = currdescr.split(',')[1]
#         currfunc = currdescr.split(',')[2].split()[-1].strip(' ()')
#         currID = currdescr.split(',')[0].lstrip('>')
#         curracc = currdescr.split(',')[3]
#         if getcoords:
#             try:
#                 currcoords = ( eval(currdescr.split(',')[4].split('.')[0]), eval(currdescr.split(',')[4].split('.')[2]) )
#             except:
#                 #print "Coordinates for", currlocus, curracc, "are not interpretable"
#                 pass
#         
#         # perform tests to exclude current seq
#         if locus not in currlocus:
#             continue
#         
#         if alleles != 'all':
#             if currallele != alleles:
#                 continue
#         
#         if currspecies != species:
#             continue
#                 
#         if currfunc not in func:
#             continue
#         
#         # check for "partial"
#         if len(currdescr.split(',')) > 5 and 'partial' in ' '.join(currdescr.split(',')[5:]):
#             continue
#         
#         refIDs.append( currID )
#         refacc[currID] = curracc
#         refseqs[currID] = seq.seq.tostring().upper()
#         if getcoords:
#             refcoords[currID] = currcoords
#         
#     refIDs.append( '' )
#     refIDs.sort()
#         
#     if not getcoords:
#         return refIDs,refacc,refseqs
#     return refIDs,refacc,refseqs,refcoords
# 
# 
# 
# # ================================
# # = Definition of reference data =
# # ================================
# 
# LOCI = ['IGH','IGK','IGL','TRA','TRB','TRD','TRG']
# 
# IGHV_list,IGHV_acc,IGHV_seqs,IGHV_coords = get_refseq_elements(locus='IGHV',alleles='01',func=['F','ORF'],species='Homo+sapiens',getcoords=True,refdatadir=params.refdatadir,imgtrefseqfasta=params.imgtrefseqfasta)
# IGHD_list,IGHD_acc,IGHD_seqs = get_refseq_elements(locus='IGHD',alleles='01',func=['F','ORF'],species='Homo+sapiens',getcoords=False,refdatadir=params.refdatadir,imgtrefseqfasta=params.imgtrefseqfasta)
# IGHJ_list,IGHJ_acc,IGHJ_seqs,IGHJ_coords = get_refseq_elements(locus='IGHJ',alleles='01',func=['F','ORF'],species='Homo+sapiens',getcoords=True,refdatadir=params.refdatadir,imgtrefseqfasta=params.imgtrefseqfasta)
# 
# IGHV_idx = dict([(g,i) for i,g in enumerate(IGHV_list)])
# IGHD_idx = dict([(g,i) for i,g in enumerate(IGHD_list)])
# IGHJ_idx = dict([(g,i) for i,g in enumerate(IGHJ_list)])
# 
# IGHC_list = [
#             '',
#             'IGHA1',
#             'IGHA2',
#             'IGHD',
#             'IGHE',
#             'IGHG1',
#             'IGHG2',
#             'IGHG3',
#             'IGHG4',
#             'IGHM'
#             ]
# IGHC_idx = dict([(g,i) for i,g in enumerate(IGHC_list)])
# 
# ALL_IDs = list(set(IGHV_list + IGHD_list + IGHJ_list + IGHC_list))
# ALL_IDs.sort()
# 
# 
# 
# # ====================
# # = Specificity data =
# # ====================
# 
# def get_LIGM_with_specificities(refdatadir,imgtdat,imgtfasta,outputfasta,outputvdjxml):
#     '''
#     Take the IMGT LIGM flat and fasta files, and return a fasta with only those seqs
#     that have a specificity associated with them in a fasta file.  The header contains
#     the accession and the specificity.
#     
#     vdj.refseq.get_LIGM_with_specificities("imgt.dat","imgt.fasta","imgt.specificities.fasta","imgt.specificities.vdjxml")
#     '''
#     specificities = {}      # list of 2-tuples: (ACCESION,specificity)
#     LIGMflat = open(os.path.join(refdatadir,imgtdat),'r')
#     LIGMfasta = open(os.path.join(refdatadir,imgtfasta),'r')
#     opSpecificityfasta = open(os.path.join(refdatadir,outputfasta),'w')
#     
#     numRecords = 0
#     numRecordswithSpec = 0
#     
#     ID = ''
#     DE = ''
#     for line in LIGMflat:
#         splitline = line.split()
#         
#         if splitline[0] == 'ID':
#             inRecord = True
#             ID = splitline[1]
#             numRecords += 1
#         elif splitline[0] == 'DE':
#             if inRecord:
#                 DE += ' '.join(splitline[1:]) + ' '
#         elif splitline[0] == 'XX':
#             if DE == '':    # if i haven't stored the description yet
#                 continue
#             else:   # finished record
#                 if 'specificity' in DE:
#                     specidx = DE.rfind('specificity')
#                     spec = DE[specidx+len('specificity'):].strip().rstrip('.')
#                     if spec.startswith('anti'):
#                         specificities[ID] = spec
#                         numRecordswithSpec += 1
#                 ID = ''
#                 DE = ''
#                 inRecord = False
#     
#     print "Number of LIGM records read: " + str(numRecords)
#     print "Number of LIGM records that have specificities: " + str(numRecordswithSpec)
#     
#     numFasta = 0
#     
#     for seq in SeqIO.parse(LIGMfasta,'fasta'):
#         spec = specificities.get(seq.id,'')
#         if spec == '':
#             continue
#         else:
#             print >>opSpecificityfasta, ">" + seq.id + " | " + spec
#             print >>opSpecificityfasta, seq.seq.tostring().upper()
#             numFasta += 1
#     
#     print "Number of Fasta records with specificities found and printed: " + str(numFasta)
#     
#     LIGMflat.close()
#     LIGMfasta.close()
#     opSpecificityfasta.close()
#     
#     # VDJXML and alignment
#     rep = vdj.initial_import([os.path.join(refdatadir,outputfasta)],os.path.join(refdatadir,outputvdjxml),metatags=['specificity_reference : '+vdj.timestamp()])
#     rep = vdj.positive_strand(rep)
#     rep = vdj.align_rep(rep)
#     
#     #repfiltered = rep.get_chains_fullVJ()
#     repfiltered = rep
#     for chain in repfiltered:
#         spec = specificities.get(chain.descr,'')
#         if spec == '':
#             print "Reference chain has empty specificity: " + chain.descr
#             continue
#         else:
#             chain.add_tags( 'specificity|'+spec )
#     
#     vdj.writeVDJ(repfiltered,os.path.join(refdatadir,outputvdjxml)) 
#     
#     return
# 
# 
# if not os.path.exists(os.path.join(params.refdatadir,params.imgtspecfasta)) or not os.path.exists(os.path.join(params.refdatadir,params.imgtspecvdjxml)):
#     get_LIGM_with_specificities(params.refdatadir,params.imgtdat,params.imgtfasta,params.imgtspecfasta,params.imgtspecvdjxml)
# 
# 
# if os.path.exists(os.path.join(params.refdatadir,params.imgtspecfasta)) and os.path.exists(os.path.join(params.refdatadir,params.imgtspecvdjxml)):
#     ipspecfasta = open(os.path.join(params.refdatadir,params.imgtspecfasta),'r')
#     
#     SPEC_list = set()
#     for line in ipspecfasta:
#         if line[0] == '>':
#             currspec = '|'.join(line.split('|')[1:]).strip()
#             SPEC_list.add(currspec)
#     ipspecfasta.close()
#     SPEC_list.add('')
#     SPEC_list = list(SPEC_list)
#     SPEC_list.sort()
