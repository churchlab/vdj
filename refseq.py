# Copyright 2014 Uri Laserson
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
import sys
import types
import glob
import warnings

from Bio import SeqIO

from seqtools import FastaIterator

import params

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
# 
# Note that you should ensure params.organism is properly set before calling
# refseq.

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

def get_X_REGION_feature(record,X):
    for feature in record.features:
        if feature.type == ('%s-REGION'%X):
            return feature

def get_any_X_REGION_feature(record):
    for feature in record.features:
        if feature.type in ['V-REGION','D-REGION','J-REGION']:
            return feature

def has_CYS(record):
    for feature in record.features:
        if feature.type == '2nd-CYS':
            return True
    return False

def has_TRP_PHE(record):
    for feature in record.features:
        if feature.type in ['J-TRP','J-PHE']:
            return True
    return False


ref_dir = os.path.join(params.vdj_dir,params.data_dir)
ligm_index = None   # LIGM not loaded by default

# figure out which files to process
processed_dir = os.path.join(params.vdj_dir,params.processed_dir,params.organism)
if not os.path.exists(processed_dir):
    os.makedirs(processed_dir,mode=0755)

reference_files = glob.glob( os.path.join(ref_dir,'*.fasta') )
processed_files = glob.glob( os.path.join(processed_dir,'*.imgt') )

get_group = lambda f: os.path.splitext(os.path.basename(f))[0]
groups_to_process = set(map(get_group,reference_files)) - set(map(get_group,processed_files))
files_to_process = filter(lambda f: get_group(f) in groups_to_process,reference_files)

# process those files
for reference_fasta in files_to_process:
    if ligm_index == None: ligm_index = SeqIO.index( os.path.join(params.imgt_dir,'imgt.dat'), 'imgt')
    group = get_group(reference_fasta)
    reference_imgt = os.path.join(processed_dir,group+'.imgt')
    
    reference_records = []
    for record in SeqIO.parse(reference_fasta,'fasta'):
        #DEBUG
        # print record.id; sys.stdout.flush()
        # if record.id == 'M62379|TRBV2*02|Homo':
        #     import pdb
        #     pdb.set_trace()
        
        try:
            header_data = parse_imgt_fasta_header(record)
        except ValueError, e:
            print >>sys.stderr, "Problem with parsing fasta header\n%s\n%s" % (record.description,e)
            continue
        
        imgt_record = ligm_index[header_data['accession']]
        
        # prune down to the V-REGION annotated in V-QUEST fasta file
        reference_record = imgt_record[header_data['coords'][0]:header_data['coords'][1]]
        reference_feature = get_X_REGION_feature(reference_record,group[-1])
        if reference_feature == None:
            print >>sys.stderr, "Couldn't find a %s-REGION feature in %s when clipped to fasta coords.  Skipping." % (group[-1],reference_record.id)
            continue
        if not reference_feature.qualifiers.has_key('allele'):
            warnings.warn("Reference record %s's %s is missing an `allele` qualifier; it (%s) was added manually." % (reference_record.id,reference_feature.type,header_data['allele']))
            reference_feature.qualifiers['allele'] = [header_data['allele']]
        elif reference_feature.qualifiers['allele'][0] != header_data['allele']:
            warnings.warn("Reference record %s's %s is annotated %s while the corres. fasta header annotates it as %s.  I picked the LIGM version." % (reference_record.id,reference_feature.type,reference_feature.qualifiers['allele'][0],header_data['allele']))
        
        # ensure that reference sequence has either 2nd-CYS, J-TRP, or J_PHE
        if group[-1] == 'V' and not has_CYS(reference_record):
            continue
        elif group[-1] == 'J' and not has_TRP_PHE(reference_record):
            continue
        
        reference_records.append(reference_record)
    
    SeqIO.write(reference_records,reference_imgt,'imgt')

get_allele = lambda record: get_any_X_REGION_feature(record).qualifiers['allele'][0]

# load the reference data (as dict: key=allele, value=SeqRecord of reference)
IGHV = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'IGHV.imgt'),'imgt')])
IGHD = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'IGHD.imgt'),'imgt')])
IGHJ = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'IGHJ.imgt'),'imgt')])
IGKV = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'IGKV.imgt'),'imgt')])
IGKJ = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'IGKJ.imgt'),'imgt')])
IGLV = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'IGLV.imgt'),'imgt')])
IGLJ = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'IGLJ.imgt'),'imgt')])
TRBV = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'TRBV.imgt'),'imgt')])
TRBD = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'TRBD.imgt'),'imgt')])
TRBJ = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'TRBJ.imgt'),'imgt')])
TRAV = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'TRAV.imgt'),'imgt')])
TRAJ = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'TRAJ.imgt'),'imgt')])
TRDV = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'TRDV.imgt'),'imgt')])
TRDD = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'TRDD.imgt'),'imgt')])
TRDJ = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'TRDJ.imgt'),'imgt')])
TRGV = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'TRGV.imgt'),'imgt')])
TRGJ = dict([(get_allele(r),r) for r in SeqIO.parse(os.path.join(processed_dir,'TRGJ.imgt'),'imgt')])
ALL = reduce(lambda a,b: dict(a.items()+b.items()),[IGHV,IGHD,IGHJ,IGKV,IGKJ,IGLV,IGLJ,TRBV,TRBD,TRBJ,TRAV,TRAJ,TRDV,TRDD,TRDJ,TRGV,TRGJ])

# load gapped sequences
data_dir = os.path.join(params.vdj_dir,params.data_dir)
IGHV_gapped = dict(FastaIterator(os.path.join(data_dir,'IGHV.fasta'),lambda s: s.split('|')[1]))
IGHD_gapped = dict(FastaIterator(os.path.join(data_dir,'IGHD.fasta'),lambda s: s.split('|')[1]))
IGHJ_gapped = dict(FastaIterator(os.path.join(data_dir,'IGHJ.fasta'),lambda s: s.split('|')[1]))
IGKV_gapped = dict(FastaIterator(os.path.join(data_dir,'IGKV.fasta'),lambda s: s.split('|')[1]))
IGKJ_gapped = dict(FastaIterator(os.path.join(data_dir,'IGKJ.fasta'),lambda s: s.split('|')[1]))
IGLV_gapped = dict(FastaIterator(os.path.join(data_dir,'IGLV.fasta'),lambda s: s.split('|')[1]))
IGLJ_gapped = dict(FastaIterator(os.path.join(data_dir,'IGLJ.fasta'),lambda s: s.split('|')[1]))
TRBV_gapped = dict(FastaIterator(os.path.join(data_dir,'TRBV.fasta'),lambda s: s.split('|')[1]))
TRBD_gapped = dict(FastaIterator(os.path.join(data_dir,'TRBD.fasta'),lambda s: s.split('|')[1]))
TRBJ_gapped = dict(FastaIterator(os.path.join(data_dir,'TRBJ.fasta'),lambda s: s.split('|')[1]))
TRAV_gapped = dict(FastaIterator(os.path.join(data_dir,'TRAV.fasta'),lambda s: s.split('|')[1]))
TRAJ_gapped = dict(FastaIterator(os.path.join(data_dir,'TRAJ.fasta'),lambda s: s.split('|')[1]))
TRDV_gapped = dict(FastaIterator(os.path.join(data_dir,'TRDV.fasta'),lambda s: s.split('|')[1]))
TRDD_gapped = dict(FastaIterator(os.path.join(data_dir,'TRDD.fasta'),lambda s: s.split('|')[1]))
TRDJ_gapped = dict(FastaIterator(os.path.join(data_dir,'TRDJ.fasta'),lambda s: s.split('|')[1]))
TRGV_gapped = dict(FastaIterator(os.path.join(data_dir,'TRGV.fasta'),lambda s: s.split('|')[1]))
TRGJ_gapped = dict(FastaIterator(os.path.join(data_dir,'TRGJ.fasta'),lambda s: s.split('|')[1]))
ALL_gapped = reduce(lambda a,b: dict(a.items()+b.items()),[IGHV_gapped,IGHD_gapped,IGHJ_gapped,IGKV_gapped,IGKJ_gapped,IGLV_gapped,IGLJ_gapped,TRBV_gapped,TRBD_gapped,TRBJ_gapped,TRAV_gapped,TRAJ_gapped,TRDV_gapped,TRDD_gapped,TRDJ_gapped,TRGV_gapped,TRGJ_gapped])




# ****************************************************************************


# DEPRECATED FUNCTIONS
# 
# The functions below are deprecated, but are saved, since they implemented
# some functionality that was complicated, and that I may reuse in the future.


# Pull LIGM records from the IMGT web server.

# def pull_LIGM_record(self):
#     """Get SeqRecord object for LIGM record from IMGT server"""
#     
#     # NOTE: this can potentially be significantly simplified by accessing the URL
#     # interface to LIGM, through:
#     # http://imgt.cines.fr/cgi-bin/IMGTlect.jv?query=5+numacc
#     # where numacc is the accession number
#     
#     request = urllib2.Request('http://imgt.cines.fr/cgi-bin/IMGTlect.jv?livret=0')
#     # LIGM page
#     response = urllib2.urlopen(request)
#     forms = ClientForm.ParseResponse(response,
#                                      form_parser_class=ClientForm.XHTMLCompatibleFormParser,
#                                      backwards_compat=False)
#     form = forms[1]
#     form['l01p01c02'] = self.accession
#     request2 = form.click()
#     # data format page
#     response2 = urllib2.urlopen(request2)
#     forms2 = ClientForm.ParseResponse(response2,
#                                      form_parser_class=ClientForm.XHTMLCompatibleFormParser,
#                                      backwards_compat=False)
#     form2 = forms2[0]
#     assert( form2.controls[8].attrs['value'] == '2 IMGT flat-file' )
#     form2.controls[8].id = 'flatfile'
#     request3 = form2.click(id='flatfile')
#     # LIGM record results
#     response3 = urllib2.urlopen(request3)
# 
#     # ghetto parse of the results.  the text of the LIGM record is in <pre>...</pre> tags
#     rawdata1 = response3.read()
#     rawdata2 = rawdata1.split('<pre>')[1].split('</pre>')[0].lstrip()
#     rawdata3 = StringIO(rawdata2)
#     self.record = SeqIO.read(rawdata3,'imgt')


# Get specificity information from LIGM.  This was written before I
# implemented the biopython IMGT parser, so it's extra-ugly.  It is probably
# much easier to implement now.

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
