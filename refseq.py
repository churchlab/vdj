# refseq.py

# contain definitions from ref IMGT database
# accessed around 1 Oct 2008

from Bio import SeqIO
import numpy as np

LOCI = ['IGH','IGK','IGL','TRA','TRB','TRD','TRG']

IGHV = [
        'IGHV1-18',
        'IGHV1-2',
        'IGHV1-24',
        'IGHV1-3',
        'IGHV1-45',
        'IGHV1-46',
        'IGHV1-58',
        'IGHV1-69',
        'IGHV1-8',
        'IGHV1-f',
        'IGHV1/OR15-1',
        'IGHV1/OR15-3',
        'IGHV1/OR15-5',
        'IGHV1/OR15-9',
        'IGHV2-26',
        'IGHV2-5',
        'IGHV2-70',
        'IGHV2/OR16-5',
        'IGHV3-11',
        'IGHV3-13',
        'IGHV3-15',
        'IGHV3-20',
        'IGHV3-21',
        'IGHV3-23',
        'IGHV3-30',
        'IGHV3-30-3',
        'IGHV3-30-5',
        'IGHV3-33',
        'IGHV3-43',
        'IGHV3-48',
        'IGHV3-49',
        'IGHV3-53',
        'IGHV3-64',
        'IGHV3-66',
        'IGHV3-7',
        'IGHV3-72',
        'IGHV3-73',
        'IGHV3-74',
        'IGHV3-9',
        'IGHV3-d',
        'IGHV3-16',
        'IGHV3-35',
        'IGHV3-38',
        'IGHV3/OR15-7',
        'IGHV3/OR16-10',
        'IGHV3/OR16-12',
        'IGHV3/OR16-13',
        'IGHV3/OR16-8',
        'IGHV3/OR16-9',
        'IGHV4-28',
        'IGHV4-30-1',
        'IGHV4-30-2',
        'IGHV4-30-4',
        'IGHV4-31',
        'IGHV4-34',
        'IGHV4-39',
        'IGHV4-4',
        'IGHV4-59',
        'IGHV4-61',
        'IGHV4-b',
        'IGHV5-51',
        'IGHV5-a',
        'IGHV6-1',
        'IGHV7-4-1',
        'IGHV4-61',
        'IGHV4/OR15-8',
        'IGHV7-81'
       ]
IGHV.sort()
IGHV = [''] + IGHV  # for unclassifiable

# change to dictionary comprehension in python 3
IGHVdict = dict([(g,i) for i,g in enumerate(IGHV)])

IGHD = [
        'IGHD1-1',
        'IGHD1-20',
        'IGHD1-26',
        'IGHD1-7',
        'IGHD2-15',
        'IGHD2-2',
        'IGHD2-21',
        'IGHD2-8',
        'IGHD3-10',
        'IGHD3-16',
        'IGHD3-22',
        'IGHD3-3',
        'IGHD3-9',
        'IGHD4-17',
        'IGHD4-4',
        'IGHD5-12',
        'IGHD5-18',
        'IGHD5-5',
        'IGHD6-13',
        'IGHD6-19',
        'IGHD6-25',
        'IGHD6-6',
        'IGHD7-27',
        'IGHD1-14',
        'IGHD1/OR15-1a',
        'IGHD1/OR15-1b',
        'IGHD2/OR15-2a',
        'IGHD2/OR15-2b',
        'IGHD3/OR15-3a',
        'IGHD3/OR15-3b',
        'IGHD4-11',
        'IGHD4-23',
        'IGHD4/OR15-4a',
        'IGHD4/OR15-4b',
        'IGHD5-24',
        'IGHD5/OR15-5a',
        'IGHD5/OR15-5b'
       ]
IGHD.sort()
IGHD = [''] + IGHD  # for unclassifiable

IGHDdict = dict([(g,i) for i,g in enumerate(IGHD)])

IGHJ = [
        'IGHJ1',
        'IGHJ2',
        'IGHJ3',
        'IGHJ4',
        'IGHJ5',
        'IGHJ6'
       ]
IGHJ.sort()
IGHJ = [''] + IGHJ  # for unclassifiable

IGHJdict = dict([(g,i) for i,g in enumerate(IGHJ)])

IGHC = [
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
IGHC.sort()
IGHC = [''] + IGHC  # for unclassifiable

IGHCdict = dict([(g,i) for i,g in enumerate(IGHC)])


# structures to find VDJ etc combos from raveled numpy arrays
VJref = np.empty( shape=(len(IGHV),len(IGHJ)), dtype=np.object )
for (i,Vgene) in enumerate(IGHV):
	for (j,Jgene) in enumerate(IGHJ):
		VJref[i,j] = [Vgene,Jgene]
VJref = VJref.ravel()


VDJref = np.empty( shape=(len(IGHV),len(IGHD),len(IGHJ)), dtype=np.object )
for (i,Vgene) in enumerate(IGHV):
	for (j,Dgene) in enumerate(IGHD):
		for (k,Jgene) in enumerate(IGHJ):
			VDJref[i,j,k] = [Vgene,Dgene,Jgene]
VDJref = VDJref.ravel()

def VJCDR3ref(data):
	# data is array of counts
	num_elts = 1
	for d in data.shape:
		num_elts *= d
	
	if num_elts % len(IGHV) != 0 or (num_elts/len(IGHV)) % len(IGHJ) != 0:
		raise ValueError, 'dimensions of data array do not correspond with numbers of V, D, and J elts'
	
	# NOTE: first CDR3 len has to be 1aa
	num_CDR3_vals = num_elts / (len(IGHV)*len(IGHJ))
	
	VJCDR3refdata = np.empty( shape=(len(IGHV),len(IGHJ),num_CDR3_vals), dtype=np.object )
	for (i,Vgene) in enumerate(IGHV):
		for (j,Jgene) in enumerate(IGHJ):
			for cdrlen in xrange(num_CDR3_vals):
				VJCDR3refdata[i,j,cdrlen] = [Vgene,Jgene,3*(cdrlen+1)]
	VJCDR3refdata = VJCDR3refdata.ravel()
	
	return VJCDR3refdata



# after downloading fresh imgtrefseq.fasta, process it into
# files for use by ABACUS
# ftp://ftp.cines.fr/IMGT/imgtrefseq.fasta

def splitrefseqs(directory,infile):
    ip  = open(directory+infile,'r')
    opV = open(directory+'V3.txt','w')
    opD = open(directory+'D3.txt','w')
    opJ = open(directory+'J3.txt','w')

    for seq in SeqIO.parse(ip,'fasta'):
        ident = seq.id.split(',')[0]
        ident2 = ident.split('*')
        species = seq.id.split(',')[1]
        seq.id = '|'+ident
        seq.description = ''
        if (ident[0:3] == 'IGH') and (species=='Homo+sapiens'):
            if ident[3] == 'V':
                if (ident2[0] in refseq.IGHV) and (ident2[1] == '01'):
                    SeqIO.write([seq],opV,'fasta')
            elif ident[3] == 'D':
                if (ident2[0] in refseq.IGHD) and (ident2[1] == '01'):
                    SeqIO.write([seq],opD,'fasta')
            elif ident[3] == 'J':
                if (ident2[0] in refseq.IGHJ) and (ident2[1] == '01'):
                    SeqIO.write([seq],opJ,'fasta')

    ip.close()
    opV.close()
    opD.close()
    opJ.close()


# All IGHV FR3 seqs
# http://imgt.cines.fr/IMGT_GENE-DB/GENElect?query=8.1+IGHV&species=Homo+sapiens&IMGTlabel=FR3-IMGT
#
# All IGHJ coding regions
# http://imgt.cines.fr/IMGT_GENE-DB/GENElect?query=8.1+IGHJ&species=Homo+sapiens&IMGTlabel=FR4-IMGT
#
# All IGHJ conserved TRP
# http://imgt.cines.fr/IMGT_GENE-DB/GENElect?query=8.1+IGHJ&species=Homo+sapiens&IMGTlabel=J-TRP
#
#
# All IGHV full GAPPED seqs:
# http://imgt.cines.fr/IMGT_GENE-DB/GENElect?query=7.1+IGHV&species=Homo+sapiens

# Load dictionaries (referenced by gene name)
# the IGHV dictionary contains the FR3 list, so that the next codon should be CYS
# the IGHJ dictionary returns a pair (pos of J-TRP,sequence of J-REGION)

import os
if os.path.exists('/Users/laserson/'):
	align_ref_dir = '/Users/laserson/research/church/vdj-ome/ref-data/align/'
elif os.path.exists('/home/ul2/'):
	align_ref_dir = '/home/ul2/vdj-ome/ref-data/align/'

ref_FR3_file      = 'IGHV.FR3-IMGT.20090205.fasta'
ref_J_REGION_file = 'IGHJ.J-REGION.20090213.fasta'
ref_J_TRP_file    = 'IGHJ.J-TRP.20090213.fasta'

# IGHV_FR3
ip = open(align_ref_dir+ref_FR3_file,'r')
IGHV_FR3 = {}
for seq in SeqIO.parse(ip,'fasta'):
    data = seq.description.split('|')
    gene,allele = data[1].split('*')[0],data[1].split('*')[1]
    func = data[3]
    if gene in IGHV and allele == '01':
        IGHV_FR3[gene] = seq.seq.data.upper()

ip.close()

# IGHJ_FR4
ipREGION = open(align_ref_dir+ref_J_REGION_file,'r')
ipTRP = open(align_ref_dir+ref_J_TRP_file,'r')

j_region = list(SeqIO.parse(ipREGION,'fasta'))
j_trp = list(SeqIO.parse(ipTRP,'fasta'))

IGHJ_FR4 = {}
for trp in j_trp:
    data1 = trp.description.split('|')
    gene1,allele1 = data1[1].split('*')[0],data1[1].split('*')[1]
    if gene1 in IGHJ and allele1 == '01':
        abs_trp_pos = eval(data1[5].split('.')[0])
        for region in j_region:
            data2 = region.description.split('|')
            gene2,allele2 = data2[1].split('*')[0],data2[1].split('*')[1]
            if gene2 == gene1 and allele2 == '01':
                abs_reg_start = eval(data2[5].split('.')[0])
                rel_pos = abs_trp_pos - abs_reg_start
                IGHJ_FR4[gene2] = (rel_pos,region.seq.data.upper())

ipREGION.close()
ipTRP.close()
