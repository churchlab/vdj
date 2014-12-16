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

# functions for pipeline operations

import warnings

from Bio.Seq import Seq
from Bio.Alphabet import DNAAlphabet

import vdj
import alignment
import seqtools

def parse_jobfile(filename):
    parameters = {}
    ip = open(filename,'r')
    for line in ip:
        data = line.split('#')[0].strip()
        if data == '': continue
        data = data.split('\t')
        name = data[0].strip()
        value = data[1].strip()
        if value == '': continue
        
        if name == 'locus':
            parameters.setdefault(name,[]).append(value)
            continue
        if name == 'rigorous':
            parameters[name] = bool(value)
            continue
        parameters[name] = value
    ip.close()
    return parameters

def iterator2parts(iterator,basename,packetsize,prefix='',suffix=''):
    """Split data from iterator into multiple files"""
    parts = []
    num_processed = 0
    file_num = 1

    curr_outname = basename+'.'+str(file_num)
    
    for obj in iterator:
        if num_processed == 0:
            op = open(curr_outname,'w')
            print >>op, prefix
            parts.append(curr_outname)

        print >>op, obj
        num_processed += 1

        if num_processed == packetsize:
            print >>op, suffix
            op.close()
            num_processed = 0
            file_num += 1
            curr_outname = basename+'.'+str(file_num)
    if not op.closed:
        print >>op, suffix
        op.close()
    
    return parts

def load_barcodes(barcode_file):
    bcip = open(barcode_file,'r')
    barcodes = {}
    for (descr,seq) in seqtools.FastaIterator(bcip):
        barcodes[seq.upper()] = descr
    bcip.close()
    
    # check that barcodes meet necessary criteria
    barcode_len = len(barcodes.keys()[0])
    for bc in barcodes.keys():
        if len(bc) != barcode_len:
            raise Exception, "ERROR: All barcode lengths must be equal."
    
    return barcodes

def id_barcode(chain,barcodes):
    # barcodes assumed to be single length
    barcode_len = len(barcodes.keys()[0])
    try:
        curr_barcode = barcodes[chain.seq.tostring()[:barcode_len].upper()]
    except KeyError:    # barcode not found; chain unchanged
        return    # chain remains unchanged
    chain.__init__(chain[barcode_len:])
    chain.annotations['barcode'] = curr_barcode

def load_isotypes(isotype_file):
    ighcip = open(isotype_file,'r')
    isotypes = {}
    for (descr,seq) in seqtools.FastaIterator(ighcip):
        isotypes[seq.upper()] = descr
    ighcip.close()
    
    return isotypes

def id_isotype(chain,isotypes):
    if not chain.has_tag('positive') and not chain.has_tag('coding'):
        warnings.warn('chain %s may not be the correct strand' % chain.descr)
    
    for iso in isotypes.iteritems():
        if iso[0] in chain.seq[-50:]:   # arbitrary cutoff from 3' end
            chain.annotations['c'] = iso[1]

def fasta2imgt(inhandle,outhandle):
    for seq in SeqIO.parse(inhandle,'fasta'):
        chain = vdj.ImmuneChain(seq)
        print >>outhandle, chain

def imgt2fasta(inhandle,outhandle):
    for chain in vdj.parse_imgt(inhandle):
        outhandle.write( chain.format('fasta') )

def partition_VJ(inhandle,basename):
    # ignores allele numbers
    def vj_id_no_allele(chain):
        return seqtools.cleanup_id(chain.v.split('*')[0]) + '_' + seqtools.cleanup_id(chain.j.split('*')[0])
    
    def outname(basename,vj_id):
        return "%s.%s.imgt" % (basename,vj_id)
    
    outhandles = {}
    for chain in vdj.parse_imgt(inhandle):
        curr_vj_id = vj_id_no_allele(chain)
        try:
            print >>outhandles[curr_vj_id], chain
        except KeyError:
            outhandles[curr_vj_id] = open( outname(basename,curr_vj_id), 'w' )
            print >>outhandles[curr_vj_id], chain
    
    for outhandle in outhandles.itervalues():
        outhandle.close()
    
    return [outname(basename,vj_id) for vj_id in outhandles.iterkeys()]

def translate_chain( chain ):
    chain.annotations['translation'] = chain.seq.translate().tostring()
    for feature in chain.features:
        offset = int(feature.qualifiers.get('codon_start',[1])[0]) - 1
        feature.qualifiers['translation'] = feature.extract(chain.seq)[offset:].translate().tostring()

def sequence_force_in_frame(chain,replace=False):
    nt = ''
    cdr3_start = alignment.vdj_aligner.ungapped2gapped_coord(chain.seq.tostring(),chain.annotations['gapped_query'],chain.__getattribute__('CDR3-IMGT').location.nofuzzy_start)
    cdr3_len = len(chain.junction_nt)
    extra_junction = cdr3_len % 3
    for (i,(r,q)) in enumerate(zip(chain.annotations['gapped_reference'],chain.annotations['gapped_query'])):
        # am I in the CDR3?
        if i >= cdr3_start and i < cdr3_start + cdr3_len:
            in_cdr3 = True
        elif i == cdr3_start + cdr3_len:
            in_cdr3 = False
            nt += '-' * ((3 - (cdr3_len % 3)) % 3)
        else:
            in_cdr3 = False
        
        if r == '-':    # insertion in the query
            nt += '' if not in_cdr3 else q.upper()
        elif q == '-':  # deletion in query
            nt += r.lower() if replace else q
        else:
            nt += q.upper()
    return nt

def translate_chain_force_in_frame(chain):
    nt = sequence_force_in_frame(chain,replace=False)
    return Seq(nt.replace('-','N'),DNAAlphabet()).translate().tostring()
