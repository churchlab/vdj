# functions for pipeline operations

import vdj
import vdj.clustering
import seqtools

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
        curr_barcode = barcodes[chain.seq[:barcode_len].upper()]
    except KeyError:    # barcode not found; chain unchanged
        return    # chain remains unchanged
    chain.seq = chain.seq[barcode_len:] # prune off barcode from seq
    chain.barcode = curr_barcode

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
            chain.c = iso[1]

def fasta2vdjxml(inhandle,outhandle):
    print >>outhandle, "<root>"
    for (descr,seq) in seqtools.FastaIterator(inhandle,lambda d: d.split()[0]):
        chain = vdj.ImmuneChain(descr=descr,seq=seq)
        print >>outhandle, chain
    print >>outhandle, "</root>"

def size_select(inhandle,outhandle,min_size,max_size):
    print >>outhandle, "<root>"
    for chain in vdj.parse_VDJXML(inhandle):
        if len(chain) >= min_size and len(chain) <= max_size:
            print >>outhandle, chain
    print >>outhandle, "</root>"

def filter_VJ(inhandle,outhandle):
    print >>outhandle, "<root>"
    for chain in vdj.parse_VDJXML(inhandle):
        if hasattr(chain,'v') and hasattr(chain,'j'):
            print >>outhandle, chain
    print >>outhandle, "</root>"

def cat_vdjxml(files,outhandle):
    print >>outhandle, "<root>"
    for f in files:
        inhandle = open(f,'r')
        for chain in vdj.parse_VDJXML(inhandle):
            print >>outhandle, chain
    print >>outhandle, "</root>"

def partition_VJ(inhandle,basename):
    # ignores allele numbers
    def vj_id_no_allele(chain):
        return seqtools.cleanup_id(chain.v.split('*')[0]) + '_' + seqtools.cleanup_id(chain.j.split('*')[0])
    
    def outname(basename,vj_id):
        return "%s.%s.vdjxml" % (basename,vj_id)
    
    outhandles = {}
    for chain in vdj.parse_VDJXML(inhandle):
        curr_vj_id = vj_id_no_allele(chain)
        try:
            print >>outhandles[curr_vj_id], chain
        except KeyError:
            outhandles[curr_vj_id] = open( outname(basename,curr_vj_id), 'w' )
            print >>outhandles[curr_vj_id], "<root>"
            print >>outhandles[curr_vj_id], chain
    
    for outhandle in outhandles.itervalues():
        print >>outhandle, "</root>"
    
    return [outname(basename,vj_id) for vj_id in outhandles.iterkeys()]
