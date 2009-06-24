# refseq.py

# contain definitions from ref IMGT database
# accessed around 1 Oct 2008

from Bio import SeqIO
import numpy as np
import seqtools
import vdj

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

# ====================================
# = Methods to manage reference sets =
# ====================================

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


def get_LIGM_with_specificities(inputfileLIGMflat,inputfileLIGMfasta,outputfasta,outputvdjxml):
	'''
	Take the IMGT LIGM flat and fasta files, and return a fasta with only those seqs
	that have a specificity associated with them in a fasta file.  The header contains
	the accession and the specificity.
	
	vdj.refseq.get_LIGM_with_specificities("imgt.dat","imgt.fasta","imgt.specificities.fasta","imgt.specificities.vdjxml")
	'''
	specificities = {}		# list of 2-tuples: (ACCESION,specificity)
	LIGMflat = open(inputfileLIGMflat,'r')
	LIGMfasta = open(inputfileLIGMfasta,'r')
	opSpecificityfasta = open(outputfasta,'w')
	
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
			if DE == '':	# if i haven't stored the description yet
				continue
			else:	# finished record
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
			print >>opSpecificityfasta, seqtools.seqString(seq)
			numFasta += 1
	
	print "Number of Fasta records with specificities found and printed: " + str(numFasta)
	
	LIGMflat.close()
	LIGMfasta.close()
	opSpecificityfasta.close()
	
	# VDJXML and alignment
	rep = vdj.initial_import([outputfasta],outputvdjxml,metatags=['specificity_reference : '+vdj.timestamp()])
	rep = vdj.positive_strand(rep)
	rep = vdj.align_rep(rep)
	
	repfiltered = rep.get_chains_fullVJ()
	for chain in repfiltered:
		spec = specificities.get(chain.descr,'')
		if spec == '':
			print "Reference chain has empty specificity: " + chain.descr
			continue
		else:
			chain.add_tags( 'specificity|'+spec )
	
	vdj.writeVDJ(repfiltered,outputvdjxml)	
	
	return


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
else:
	align_ref_dir = './'

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

SPEC = \
	['anti-(6-4) photoproduct',
	 'anti-(T,G)-A--L (poly-(Tyr,Glu)-Ala-Lys)',
	 'anti-11-dehydro-thromboxane B2 (11-dehydro TXB2)',
	 'anti-17 beta-estradiol',
	 'anti-2-phenyl oxazolone',
	 'anti-2-phenyl-5-oxazolone (phOx)',
	 'anti-4-1BB (CD137) (TNFR superfamily member, T cell activation molecule) [human]',
	 'anti-A (blood group A)',
	 'anti-A33 (transmembrane glycoprotein) [human]',
	 'anti-B (blood group B)',
	 'anti-B cell [human fetal spleen]',
	 'anti-B cell [lymphoma]',
	 'anti-B cell [lymphoma] [human]',
	 'anti-BSA',
	 'anti-Bluetongue virus VP7',
	 'anti-BrMRBC (bromelain-treated mouse red blood cells)',
	 'anti-Brucella',
	 'anti-Bt-LC-Z2Ph (Diphenyl N-[6-(biotinamido) hexanoyl] amino (4- amidinophenyl) methanephosphonate)',
	 'anti-Burkholderia pseudomallei exotoxin',
	 'anti-Burkholderia pseudomallei protease',
	 'anti-C3 [guinea pig]',
	 'anti-C3 nephritic factor',
	 'anti-C3a receptor (C3aR) [human]',
	 'anti-C5',
	 'anti-C5 [guinea pig]',
	 'anti-C5a',
	 'anti-CD10 (cALL)',
	 'anti-CD16 [human]',
	 'anti-CD18',
	 'anti-CD19',
	 'anti-CD19, anti-CD16',
	 'anti-CD1a',
	 'anti-CD20 [human]',
	 'anti-CD22',
	 'anti-CD22 [human]',
	 'anti-CD23',
	 'anti-CD25',
	 'anti-CD28',
	 'anti-CD28 [human]',
	 'anti-CD28/alpha-1 antitrypsine',
	 'anti-CD29',
	 'anti-CD3',
	 'anti-CD3 [human]',
	 'anti-CD3 [monkey]',
	 'anti-CD30',
	 'anti-CD34',
	 'anti-CD34 [human]',
	 'anti-CD37 [human]',
	 'anti-CD4',
	 'anti-CD4 (Leu3a) [murine]',
	 'anti-CD4 [human]',
	 'anti-CD40',
	 'anti-CD40 [human]',
	 'anti-CD46',
	 'anti-CD5',
	 'anti-CD5 [human]',
	 'anti-CD7 [human]',
	 'anti-CD72',
	 'anti-CD79B (Ig-beta/B29) [human]',
	 'anti-CD8',
	 'anti-CD8 [bovine]',
	 'anti-CD95 [human]',
	 'anti-CD98',
	 'anti-CHL1 (neural adhesion molecule) [mouse]',
	 'anti-Candida albicans surface antigens',
	 'anti-CasBrE (murine neurotropic ecotropic retrovirus) envelope glycoprotein VRA motif',
	 'anti-Citrus tristeza virus (CTV) coat protein',
	 'anti-DNA',
	 'anti-DNA double-stranded (ds)',
	 'anti-DNA double-stranded (ds), anti-HP8 protein',
	 'anti-DNA double-stranded (ds), anti-Sm (ribonucleoprotein RNP)',
	 'anti-DNA double-stranded (ds), anti-Sm (ribonucleoprotein RNP), anti-cardiolipin',
	 'anti-DNA double-stranded (ds), anti-Sm (ribonucleoprotein RNP), anti-cardiolipin, anti-fibronectin',
	 'anti-DNA double-stranded (ds), anti-cardiolipin',
	 'anti-DNA double-stranded (ds), anti-fibronectin',
	 'anti-DNA double-stranded (ds), anti-glomerule (glomerular-binding, glomerulotropic)',
	 'anti-DNA quadruplex',
	 'anti-DNA single-stranded (ss)',
	 'anti-DNA single-stranded (ss) [bacterial]',
	 'anti-DNA single-stranded (ss) [human]',
	 'anti-DNA single-stranded (ss), anti-DNA double-stranded (ds)',
	 'anti-DNA triplex',
	 'anti-DNA, anti-phosphocholine',
	 'anti-DNP (2,4-dinitrophenyl)',
	 'anti-Dns (anti-5-dimethylaminophtalene-1-sulfonyl)',
	 'anti-E-selectin',
	 'anti-Ebola virus',
	 'anti-Echinococcus granulosus antigen B',
	 'anti-Eimeria acervulina apical complex protein',
	 'anti-Eimeria acervulina surface antigens',
	 'anti-Entamoeba histolytica',
	 'anti-EpCAM (epithelial cell adhesion molecule)',
	 'anti-EpCAM (epithelial cell adhesion molecule) [human], anti-CD28',
	 'anti-EpCAM (epithelial cell adhesion molecule), anti-CD46',
	 'anti-Epstein-Barr virus (EBV)',
	 'anti-Epstein-Barr virus (EBV) nuclear antigen 3, epitope FLRGRAYGL, HLA-B8-restricted',
	 "anti-F(ab')2",
	 'anti-F4 alpha (patent)',
	 'anti-FB5',
	 'anti-Fas',
	 'anti-Fas [human]',
	 'anti-Fc epsilon RI (HUGO: FCERI) (receptor for Fc fragment of IgE)',
	 'anti-Fc epsilon RI (HUGO: FCERI) (receptor for Fc fragment of IgE, high affinity I)',
	 'anti-Fc gamma RI',
	 'anti-Fc gamma RIII [human]',
	 'anti-Fusarium',
	 'anti-GAT (poly(Glu60 Ala30 Tyr10))',
	 'anti-GMP-140 [human]',
	 'anti-GRP (glycin-rich cell wall protein) [cereal]',
	 'anti-Gal (Gal alpha1-3gal beta1-4glcNAc-R)',
	 'anti-Gal (beta-(1,6)-galactan)',
	 'anti-Gal (terminal galactose alpha1,3-galactose carbohydrate structures)',
	 'anti-HIV-1 (human immunodeficiency virus type 1)',
	 'anti-HIV-1 Rev protein',
	 'anti-HIV-1 Tat protein',
	 'anti-HIV-1 Vif protein',
	 'anti-HIV-1 gp120 (envelope protein)',
	 'anti-HIV-1 gp120 (envelope protein) (CD4 binding site)',
	 'anti-HIV-1 gp41 (envelope protein)',
	 'anti-HIV-1 p17 (matrix protein)',
	 'anti-HIV-1 p24 (capsid protein)',
	 'anti-HIV-1 p25',
	 'anti-HIV-1 p31 (integrase)',
	 'anti-HIV-1 p34',
	 'anti-HIV-1 reverse transcriptase (RT)',
	 'anti-HLA',
	 'anti-HLA class II',
	 'anti-HLA-A2',
	 'anti-HLA-A2, anti-HLA-A28',
	 'anti-HLA-B7',
	 'anti-HLA-DQ3',
	 'anti-HLA-DR',
	 'anti-HM1.24 (patent)',
	 'anti-HMWG (high-molecular-weight glycoprotein) [human], anti-CD28',
	 'anti-HMWG (high-molecular-weight glycoprotein) [human], anti-CD95 [human]',
	 'anti-HMWG (high-molecular-weight glycoprotein), anti-CD28',
	 'anti-HTLV-1 (human T leukaemia virus type 1)',
	 'anti-Haemophilus influenzae type a (Hia) [human]',
	 'anti-Haemophilus influenzae type b (Hib) capsular polysaccharide',
	 'anti-I (carbohydrate antigens on red blood cells)',
	 'anti-ICOS (inducible T-cell co-stimulator, alias AILIM)',
	 'anti-IgD [mouse] allotype a',
	 'anti-IgE',
	 'anti-IgE [human]',
	 'anti-IgE receptor',
	 'anti-IgG [human]',
	 'anti-IgG rheumatoid factor',
	 'anti-InlB (internalins B) [Listeria monocytogenes]',
	 'anti-L-selectin',
	 'anti-L1CAM (neural cell adhesion molecule L1) [mouse]',
	 'anti-LDL receptor',
	 'anti-La (RNP, ribonucleoprotein)',
	 'anti-Lewis b [human]',
	 'anti-MN (tumor-associated-antigen)',
	 "anti-N(alpha)-(5'-phosphopyridoxyl)-L-lysine",
	 'anti-N-glycolyl GM3 (ganglioside GM3)',
	 'anti-NP',
	 'anti-NP (4-hydroxy-3-nitrophenyl (phenolate + phenolic form))acetyl',
	 'anti-NP (4-hydroxy-3-nitrophenyl (phenolate form)) acetyl',
	 'anti-NP (4-hydroxy-3-nitrophenylacetic acid)',
	 'anti-NP (4-hydroxy-3-nitrophenylacetic acid), anti-NIP (4-hydroxy-3-iodo-5-nitrophenylacetic acid)',
	 'anti-NP (4-hydroxy-3-nitrophenylacetic acid)-ficoll',
	 'anti-Neisseria meningitidis PorA protein P1.16 epitope',
	 'anti-Neisseria meningitidis PorA protein P1.7 epitope',
	 'anti-Neisseria meningitidis outer membrane protein P1.16 epitope',
	 'anti-P-selectin',
	 'anti-P7 Paragloboside',
	 'anti-PAF receptor',
	 'anti-PAH (polycyclic aromatic hydrocarbons)',
	 'anti-Pasteurella pneumotropica',
	 'anti-Plasmodium berghei Pbs21 antigen',
	 'anti-Plasmodium berghei ookinete protein Pbs21',
	 'anti-Plasmodium circumsporozoite surface protein',
	 'anti-Plasmodium falciparum apical membrane antigen 1 (AMA1)',
	 'anti-Plasmodium falciparum merozoite surface protein-1 (MSP1) Block 2 region',
	 'anti-Plasmodium vivax',
	 'anti-Pseudomonas aeruginosa (B type) lipopolysaccharide (LPS)',
	 'anti-Pseudomonas aeruginosa (E type) lipopolysaccharide(LPS)',
	 'anti-Pseudomonas aeruginosa (F4 type)',
	 'anti-Pseudomonas aeruginosa (H1 type)',
	 'anti-Pseudomonas aeruginosa (I type) lipopolysaccharide(LPS) [human]',
	 'anti-Pseudomonas aeruginosa exotoxin',
	 'anti-Pseudomonas aeruginosa inner core oligosaccharide',
	 'anti-Pseudomonas aeruginosa lipid-A',
	 'anti-Pseudomonas aeruginosa lipopolysaccharide (LPS) A-band',
	 'anti-Pseudomonas aeruginosa lipopolysaccharide (LPS) O-antigen',
	 'anti-Pseudomonas aeruginosa lipopolysaccharide (LPS) O-antigen, serogroup 06ad',
	 'anti-Pseudomonas aeruginosa lipopolysaccharide (LPS) serotype IATS 05 O-antigen',
	 'anti-Pseudomonas aeruginosa lipopolysaccharide (LPS) serotype IATS 06 O-antigen',
	 'anti-Pseudomonas aeruginosa lipoprotein I (outer membrane protein I)',
	 'anti-Pseudomonas aeruginosa outer core oligosaccharide',
	 'anti-Pseudomonas lipopolysaccharide (LPS)',
	 'anti-RBC (red blood cells)',
	 'anti-RBC (red blood cells) [human]',
	 'anti-RBC (red blood cells) [mouse]',
	 'anti-RBC (red blood cells) [sheep]',
	 'anti-RNA',
	 'anti-RNA double-strand (ds) [rice dwarf virus (RDV)]',
	 'anti-RNA-bacteriophage fr',
	 'anti-Ralstonia solanacearum race 3 lipopolysaccharide (LPS)',
	 'anti-Rh(D) (blood group)',
	 'anti-SSEA-1',
	 'anti-SSEA-1 fucosyl',
	 'anti-Schistosoma japonicum membrane protein',
	 'anti-Shigella dysenteriae type 1 polysaccharide',
	 'anti-Sia-1b (sialylated carbohydrates) (formerly anti-Gd) cold agglutinin',
	 'anti-Siglec-9 (dendritic cell (DC) membrane molecule)',
	 'anti-Sm',
	 'anti-Sm (ribonucleoprotein RNP)',
	 'anti-Sm (ribonucleoprotein RNP), anti-DNA double-stranded (ds)',
	 'anti-Sm (ribonucleoprotein RNP), anti-DNA double-stranded (ds), anti-glomerule (glomerular-binding, glomerulotropic)',
	 'anti-Sm (ribonucleoprotein RNP), anti-DNA single-stranded (ss)',
	 'anti-Streptococcus pneumonia (serotype 23F) (capsular polysaccharide)',
	 'anti-Streptococcus pneumonia (type 6B) (capsular polysaccharide)',
	 'anti-Syk (tyrosine kinase)',
	 'anti-T cell receptor beta chain',
	 'anti-T3 (triiodo-L-thyronine) (thyroid hormone)',
	 'anti-TNF alpha (tumor necrosis factor)',
	 'anti-TNF alpha (tumor necrosis factor) [human]',
	 'anti-TNP (2,4,6-trinitrophenyl)',
	 'anti-TNP (2,4,6-trinitrophenyl) H-2b class I-restricted',
	 'anti-TRAIL-R (TNF-related apoptosis-inducing ligand receptors )',
	 'anti-TRAIL-R (TNF-related apoptosis-inducing ligand receptors)',
	 'anti-Tn (tumor-associated-antigen)',
	 'anti-Toxoplasma gondii',
	 'anti-Toxoplasma gondii P30 (major surface antigen)',
	 'anti-Toxoplasma gondii P66 antigen',
	 'anti-VSG (variant surface glycoprotein)',
	 'anti-Venezuelan equine encephalitis virus (VEE)',
	 'anti-Western equine encephalitis virus (WEE)',
	 'anti-abscisic acid',
	 'anti-acetylaminofluorene',
	 'anti-acetylcholine receptor',
	 'anti-acetyllysine',
	 'anti-acid phosphatase [Japanese radish]',
	 'anti-acidic isoferritin',
	 'anti-acidic isoferritin [Escherichia coli]',
	 'anti-activating transcription factor 1',
	 'anti-adenocarcinoma',
	 'anti-adenocarcinoma [gastric human] (GA733-2)',
	 'anti-adseverin [human]',
	 'anti-aflatoxin',
	 'anti-aldolase',
	 'anti-aldolase [human]',
	 'anti-alkaline phosphatase [human placental]',
	 'anti-alkylphenol',
	 'anti-alkylphenol ethoxylate',
	 'anti-allergen Bet v1 (birch pollen allergen)',
	 'anti-allergen Der p1 (house dust mite antigen)',
	 'anti-allergen Poa p9 (Poa pratensis) peptide 26',
	 'anti-alpha-2-antiplasmin',
	 'anti-alpha-fetoprotein',
	 'anti-alpha-fetoprotein [human]',
	 'anti-amiloride-binding protein',
	 'anti-angiotensin II',
	 'anti-anthrax toxin',
	 'anti-antigen BA46 [breast epithelial]',
	 'anti-apolipoprotein A [human] kringle domain',
	 'anti-apolipoprotein A [human] protease domain',
	 'anti-apolipoprotein A-I [human plasma]',
	 'anti-apolipoprotein B-100 [human plasma]',
	 'anti-apolipoprotein B-100 [human]',
	 'anti-apolipoprotein E [human]',
	 'anti-arsonate',
	 'anti-atrazine (pesticide)',
	 'anti-bacterial levan',
	 'anti-band 3 protein [human erythrocytes]',
	 'anti-bax inhibitor protein',
	 'anti-beet necrotic yellow vein virus (BNYVV)',
	 'anti-beta-2-glycoprotein I (beta2GPI)',
	 'anti-beta-amyloid 40',
	 'anti-beta-amyloid [human peripheral blood]',
	 'anti-beta-amyloid peptide',
	 'anti-beta-galactosidase',
	 'anti-beta-lactamase',
	 'anti-biotin',
	 'anti-bisphenol A',
	 'anti-blood coagulation factor VIII (antihemophilic factor)',
	 'anti-blood coagulation factor VIII A3-C1 domains',
	 'anti-blood coagulation factor VIII C2 domain',
	 'anti-blood group antigen [human]',
	 'anti-botulinum neurotoxin type A (BoNT/A)',
	 'anti-c-erbB-2 (HUGO: ERBB2)',
	 'anti-c-myc (HUGO: MYC) peptide',
	 'anti-cancer',
	 'anti-cancer [colorectal]',
	 'anti-cancer [gastric human]',
	 'anti-cancer-associated antigen CA125',
	 'anti-canine parvovirus capsid',
	 'anti-capsular glucuronoxylomannan (GXM) [Cryptococcus neoformans]',
	 'anti-carbohydrate antigen (Lewis Y)',
	 'anti-carbohydrate antigen (Lewis Y) (embryonic)',
	 'anti-carbohydrate antigen (Lewis Y) (tumor-associated)',
	 'anti-carboxypeptidase [bovine]',
	 'anti-carcinoembryonic antigen (CEA) (HUGO:CEACAM5)',
	 'anti-carcinoembryonic antigen (CEA71)',
	 'anti-carcinoma',
	 'anti-carcinoma [colorectal]',
	 'anti-carcinoma [human breast]',
	 'anti-carcinoma [human ovarian]',
	 'anti-carcinoma [human pulmonary]',
	 'anti-carcinoma [human]',
	 'anti-carcinoma [pulmonary]',
	 'anti-carcinoma surface antigen',
	 'anti-carcinoma surface antigen (GA733-2)',
	 'anti-cardiolipin',
	 'anti-cardiolipin, anti-beta-2-glycoprotein I (beta2GPI)',
	 'anti-cardiolipin/beta 2 glycoprotein I complex',
	 'anti-cathepsin L [human]',
	 'anti-chemokine receptor (CCR2)',
	 'anti-chloramphenicol (CAM) monoester (hydrolytic)',
	 'anti-chloramphenicol derivative',
	 'anti-chlorpyrifos-ethyl',
	 'anti-chondroitin sulphate proteoglycan (NG2)',
	 'anti-chorionic gonadotropin [human] (hcG)',
	 'anti-chorionic somatotropin [human] (HCS)',
	 'anti-chrysanthemate (insecticide)',
	 'anti-cocaine',
	 'anti-collagenase IV',
	 'anti-connective tissue growth factor',
	 'anti-core-streptavidin',
	 'anti-coronavirus',
	 'anti-cortisol',
	 'anti-crotoxin',
	 'anti-cruzipain',
	 'anti-cucumber mosaic virus (CMV)',
	 'anti-cucumber mosaic virus (CMV) (coat protein)',
	 'anti-cyclobutane type thymine dimer',
	 'anti-cyclosporine (Cs)',
	 'anti-cytidine, anti-guanosine',
	 'anti-cytochrome c',
	 'anti-cytochrome c [horse]',
	 'anti-cytochrome c [mammalian]',
	 'anti-cytochrome c [mouse]',
	 'anti-cytochrome c [pigeon]',
	 'anti-cytochrome c oxidase [mammalian] subunits VIac (nuclear encoded)',
	 'anti-cytokines [human]',
	 'anti-cytotoxic T lymphocyte-associated antigen 4 (CTLA-4)',
	 'anti-dansyl',
	 'anti-dendritic cells',
	 'anti-dengue virus',
	 'anti-deoxynivalenol (mycotoxin)',
	 'anti-desipramine (tricylic antidepressant)',
	 'anti-desmoglein 3 [human]',
	 'anti-desmosome',
	 'anti-dextran (alpha (1->6))',
	 'anti-dextran B1355',
	 'anti-dextran B512',
	 'anti-digoxin',
	 'anti-dihydrolipoamide acetyltransferase (E2 subunit of pyruvate dehydrogenase complex)',
	 'anti-dihydrotestosterone',
	 'anti-disialoganglioside 2',
	 'anti-diuron (a phenylurea herbicide)',
	 'anti-dyskerin',
	 'anti-eEF-1 alpha (translation Elongation Factor 1 alpha)',
	 'anti-eEF1A-1 (eukaryotic Elongation Factor 1A-1)',
	 'anti-epidermal growth factor receptor',
	 'anti-epithelial transmembrane glycoprotein-2 (EGP-2) [human]',
	 'anti-erbB-2',
	 'anti-estradiol',
	 'anti-estrogen receptor [calf uterus]',
	 'anti-feline calici virus',
	 'anti-feline herpes virus-1 (FHV-1)',
	 'anti-ferritin [human]',
	 'anti-fibrin [human]',
	 'anti-fibroblast activation protein, alpha (FAP)',
	 'anti-fibroblast growth factor 8 (androgen-induced) (fgf-8)',
	 'anti-fibronectin ED-B domain (a marker of angiogenesis)',
	 'anti-fibronectin [human]',
	 'anti-flavocytochrome b2',
	 "anti-flavocytochrome b2 [baker's yeast L-lactate dehydrogenase]",
	 'anti-fluorescein',
	 'anti-fluorescein (FITC)',
	 'anti-folate binding protein',
	 'anti-foot-and-mouth disease virus (FMDV)',
	 'anti-ganglioside GA1 (asialo GM1)',
	 'anti-ganglioside GD2',
	 'anti-ganglioside GD3',
	 'anti-ganglioside GM2',
	 'anti-ganglioside GM2 [human type]',
	 'anti-gibberellin A24',
	 'anti-gibberellin A4 (GA4)',
	 'anti-gingipain K protease (KGP) [Porphyromonas gingivalis]',
	 'anti-gliadin',
	 'anti-glioma [mouse]',
	 'anti-glioma surface antigen (cell cycle independent)',
	 'anti-glomerular basement membrane (GBM)',
	 'anti-glutathione',
	 'anti-glutathione S-transferase (GST)',
	 'anti-glycophorin A',
	 'anti-glycophorin A type M',
	 'anti-glycophorin A type M and type N [human]',
	 'anti-glycophorin A type N',
	 'anti-glycoprotein VI (GPVI)',
	 'anti-glycyrrhetic acid',
	 'anti-gp39 [human]',
	 'anti-granzyme B (GrB)',
	 'anti-growth hormone [human] (HGH)',
	 'anti-heat shock protein 70 (HSP70) [recombinant] [human]',
	 'anti-heat shock protein gp96 (grp94)',
	 'anti-hemagglutinin H3 [influenza virus]',
	 'anti-hemagglutinin [Porphyromonas gingivalis]',
	 'anti-hemagglutinin [influenza virus A/PR/8/34]',
	 'anti-hemagglutinin [influenza virus A]',
	 'anti-hemagglutinin [influenza virus]',
	 'anti-hepatitis A virus (HAV)capsid protein',
	 'anti-hepatitis B virus (HBV) X protein (HBx)',
	 'anti-hepatitis B virus (HBV) core protein',
	 'anti-hepatitis B virus (HBV) pre-S1 surface antigen',
	 'anti-hepatitis B virus (HBV) pre-S2 surface antigen',
	 'anti-hepatitis B virus (HBV) surface antigen (HBsAg)',
	 'anti-hepatitis B virus (HBV) surface antigen (HBsAg) (human)',
	 'anti-hepatitis C virus (HCV)',
	 'anti-hepatitis C virus (HCV) RNA-dependent RNA polymerase (RdRp)',
	 'anti-hepatitis C virus (HCV) core protein',
	 'anti-hepatitis C virus (HCV) envelope glycoprotein (E2)',
	 'anti-hepatitis C virus (HCV) non-structural 4A protein',
	 'anti-hepatitis C virus (HCV) non-structural 5A protein',
	 'anti-hepatitis C virus (HCV) non-structural protein 3 (NS3)',
	 'anti-hepatitis C virus (HCV) serine protease',
	 'anti-hepatitis E virus (HEV) ORF2 protein',
	 'anti-herpes simplex virus (HSV)',
	 'anti-herpes simplex virus (HSV) type 1',
	 'anti-herpes simplex virus (HSV) type 1 and type 2',
	 'anti-herpes simplex virus (HSV) type 1 glycoprotein B peptide (SSIEFARL)/H2-Kb class-I restricted',
	 'anti-herpes virus-1 (bovine)',
	 'anti-histamine',
	 'anti-histone',
	 'anti-histone H1',
	 'anti-histone H2A',
	 'anti-histone H2A, H4',
	 'anti-histone H2A/H2B',
	 'anti-histone H2B',
	 'anti-histone H3',
	 'anti-histone H3A, H4',
	 'anti-histone H4',
	 'anti-histone, anti-nucleosome',
	 'anti-human cerebromicrovascular endothelial cells (HCEC)',
	 'anti-human cytomegalovirus (CMV)',
	 'anti-human cytomegalovirus (CMV) 65kD antigen',
	 'anti-human cytomegalovirus (CMV) glycoprotein B',
	 'anti-human cytomegalovirus (CMV) glycoprotein B antigenic domain 2 (AD-2)',
	 'anti-human cytomegalovirus (CMV) glycoprotein H',
	 'anti-human hepatitis A virus',
	 'anti-human milk fat globule',
	 'anti-i (carbohydrate antigens on red blood cells)',
	 'anti-idiotype',
	 'anti-idiotype (A48: levan-specific BALB/c myeloma protein ABPC48)',
	 'anti-idiotype (Ab3) (anti-Neisseria meningitidis polysaccharide C)',
	 'anti-idiotype (anti-sigma receptor)',
	 'anti-idiotype (anti-tumor effect)',
	 'anti-idiotype > cocaine',
	 'anti-idiotype > cyclosporine (Cs)',
	 'anti-idiotype HLA class II',
	 'anti-idiotype [human]',
	 'anti-idiotype, anti-cancer [human]',
	 'anti-idiotype, anti-epidermal growth factor receptor',
	 'anti-imidazolinone (herbicide)',
	 'anti-insulin',
	 'anti-insulin [beef]',
	 'anti-insulin [human]',
	 'anti-integrin IIIa (GPIIIa, glycoprotein IIIa) *L33 (P1A1, HPA-1a) [platelet]',
	 'anti-integrin IIb (GPIIb, glycoprotein IIb, integrin alphaIIb chain) [human platelet]',
	 'anti-integrin IIb/IIIa (GPIIb/IIIa, glycoprotein IIb/IIIa, integrin alphaIIb-beta3, fibrinogen receptor)',
	 'anti-integrin IIb/IIIa (GPIIb/IIIa, glycoprotein IIb/IIIa, integrin alphaIIb-beta3, fibrinogen receptor) [human platelet]',
	 'anti-integrin IIb/IIIa (GPIIb/IIIa, glycoprotein IIb/IIIa, integrin alphaIIb-beta3, fibrinogen receptor) [platelet]',
	 'anti-integrin Ia/IIa (GPIa/IIa, glycoprotein Ia/IIa, integrin alpha2-beta1, collagen receptor, VLA-2) [human platelet]',
	 'anti-integrin Ib (GPIb, glycoprotein Ib) [human platelet]',
	 'anti-integrin alpha4 chain',
	 'anti-integrin alpha4-beta7',
	 'anti-intercellular adhesion molecule (ICAM)',
	 'anti-intercellular adhesion molecule-1 (ICAM-1)',
	 'anti-intercellular adhesion molecule-4 (ICAM-4)',
	 'anti-interferon gamma (IFNG)',
	 'anti-interferon gamma (IFNG) [human]',
	 'anti-interferon gamma receptor [human] (HUGO: IFNGR1)',
	 'anti-interleukin-1 (IL-1)',
	 'anti-interleukin-1 receptor (IL-1R) [human]',
	 'anti-interleukin-11 (IL-11)',
	 'anti-interleukin-18 precursor (IL-18)',
	 'anti-interleukin-18 receptor (IL-18R) [human]',
	 'anti-interleukin-2 receptor (IL-12R) [human]',
	 'anti-interleukin-2 receptor (IL-2R)',
	 'anti-interleukin-2 receptor (IL-2R) [human]',
	 'anti-interleukin-2 receptor (IL-2R) beta1 chain',
	 'anti-interleukin-2 receptor (IL-2R) p55 [human]',
	 'anti-interleukin-4 (IL-4) [human]',
	 'anti-interleukin-5 (IL-5)',
	 'anti-interleukin-5 receptor (IL-5R) [human]',
	 'anti-interleukin-6 (IL-6) [human]',
	 'anti-interleukin-6 (IL-6) [human] helix D region',
	 'anti-interleukin-6 receptor (IL-6R) [human]',
	 'anti-interleukin-8 (IL-8)',
	 'anti-interleukin-8 (IL-8) [human]',
	 'anti-interleukin-8 receptor (IL-8R) [human]',
	 'anti-keratin',
	 'anti-lactococcal bacteriophage p2 baseplate protein [Lactococcus lactis]',
	 'anti-laminin',
	 'anti-lipid A [bacterial]',
	 'anti-lipoprotein [human]',
	 'anti-liver-specific membrane lipoprotein [human]',
	 'anti-louping ill virus',
	 'anti-lysozyme',
	 'anti-lysozyme [hen egg white]',
	 'anti-lysozyme [hen egg white] [chicken]',
	 'anti-lysozyme [turkey egg white]',
	 'anti-macrophage migration inhibitory factor (MIF)',
	 'anti-malathion (pesticide)',
	 'anti-malonyl-CoA decarboxylase (MCD) [recombinant] [rat]',
	 'anti-medulloblastoma [human]',
	 'anti-medulloblastoma cell [human]',
	 'anti-melanocyte differentiation antigen TRP-2 [mouse]',
	 'anti-melanoma',
	 'anti-melanoma [human',
	 'anti-melanoma [human]',
	 'anti-melanoma [human] (anti-MAGE1 peptide/HLA-A1)',
	 'anti-melanoma [human] (anti-Pmel17/gp100 (HUGO: SILV) peptide/HLA-A2)',
	 'anti-mesothelin',
	 'anti-metal-EDTA (chelate) complex',
	 'anti-methamphetamine',
	 'anti-microcystin-LR [Microcystis aeruginosa] (cyanobacterial hepatotoxin)',
	 'anti-mitochondrial translocase (TOM70) [human]',
	 'anti-morphine',
	 'anti-mucin (MUC1)',
	 'anti-mucin (MUC1) [breast cancer associated]',
	 'anti-mucin (MUC1) core protein',
	 'anti-mucin [human breast]',
	 'anti-mucin type synthetic glycolipid',
	 'anti-multidrug resistance protein (MRP1)',
	 'anti-myelin',
	 'anti-myelin basic protein (MBP) [human]',
	 'anti-myelin basic protein (MBP) cross-reactive idiotope (CRI)',
	 'anti-myelin oligodendrocyte glycoprotein (MOG)',
	 'anti-myelin-associated glycoprotein (MAG)',
	 'anti-myelin-associated neurite growth inhibitors',
	 'anti-myeloperoxidase (MPO)',
	 'anti-myosin',
	 'anti-myosin, anti-N-acetyl-glucosamine (streptococcal polysaccharide)',
	 'anti-nerve growth factor (NGF) 30 (NGF30)',
	 'anti-neuraminidase',
	 'anti-neuroblastoma [human]',
	 'anti-neuropeptide substance P',
	 'anti-neurotensin G protein-coupled receptor NTS1 [human]',
	 'anti-neurotoxin AahI [Androctonus australis hector scorpion] (venom)',
	 'anti-nicotinic acetylcholine receptor (AChR) [human]',
	 'anti-nuclear',
	 'anti-nucleosome',
	 'anti-nucleosome, anti-DNA single-stranded (ss)',
	 'anti-nucleosome, anti-glomerule (glomerular-binding, glomerulotropic)',
	 'anti-nucleotide',
	 'anti-ovalbumin peptide (SIINFEKL)/H2-Kb class-I restricted',
	 'anti-ovomucoid',
	 'anti-oxidized LDL',
	 'anti-p-azobenzenearsonate',
	 'anti-p-azophenylarsonate',
	 'anti-p-azophenylarsonate-KLH (keyhole limpet hemocyanin)',
	 'anti-p-nitrophenyl phosphonate, esterolytic',
	 'anti-p21 RAS',
	 'anti-p53 (tumor-suppressor protein)',
	 'anti-paraquat (herbicide)',
	 'anti-parathion',
	 'anti-parathormone [human] related peptide',
	 'anti-parathyrin (PTHrP, parathyroid hormone-related protein) [human]',
	 'anti-parathyroid hormone',
	 'anti-peripherin',
	 'anti-phosphatidylcholine',
	 'anti-phosphatidylinositol-3,4,5-triphosphate',
	 'anti-phosphatidylserine',
	 'anti-phosphocholine',
	 'anti-phosphocholine [Proteus morganii]',
	 'anti-phosphocholine-KLH (keyhole limpet hemocyanin)',
	 'anti-phospholipid',
	 'anti-phospholipid [human platelet]',
	 'anti-phosphonate esters, catalytic, proteolytic',
	 'anti-phosphorylcholine',
	 'anti-phosphorylcholine (low-affinity)',
	 'anti-phthalate',
	 'anti-picloram',
	 'anti-plasminogen activator inhibitor type-1 (PAI-1)',
	 'anti-plum pox virus (PPV) RNA replicase NIb',
	 'anti-pneumococcal',
	 'anti-pneumococcal [human]',
	 'anti-pneumococcal capsular polysaccharide',
	 'anti-pneumococcal capsular polysaccharide serotype 6B',
	 'anti-pneumococcal, anti-DNA double-stranded (ds)',
	 'anti-pneumolysin',
	 'anti-poliovirus type 1',
	 'anti-poly(A)/poly(U)',
	 'anti-poly(A)/poly(dT)',
	 'anti-poly(dC)',
	 'anti-polycation',
	 'anti-porphyrin',
	 'anti-potato cyst nematode',
	 'anti-potato virus V coat protein',
	 'anti-potato virus Y (PVY) nuclear inclusion (NIa)',
	 'anti-procollagenase [human]',
	 'anti-prostate',
	 'anti-proteolipid protein (PLP)',
	 'anti-pseudouridine',
	 'anti-rabies virus',
	 'anti-respiratory syncytial virus (RSV)',
	 'anti-rhinovirus HRV2 [human]',
	 'anti-ribosomal P (phospho-) protein (systemic lupus erythematosus) (SLE)',
	 'anti-ribosomal P (phospho-) protein [eukaryote]',
	 'anti-rice stripe virus P20 protein',
	 'anti-rotavirus VP7',
	 'anti-seminoprotein [human]',
	 'anti-serum albumin [bovine]',
	 'anti-serum albumin [human]',
	 'anti-severe acute respiratory syndrome coronavirus (SARS)',
	 'anti-shiga toxin 1 (Stx1)',
	 'anti-sialoprotein protein Fv (pFv) 175 kDa [human] gut-associated',
	 'anti-signaling lymphocytic activation molecule (SLAM)',
	 'anti-simetryn (herbicide)',
	 'anti-sm',
	 'anti-solamargine',
	 'anti-sperm [human] (immobilising)',
	 'anti-sperm agglutination antigen-1 (SAGA-1) (sperm glycoform of CD25) [human] carbohydrate epitope',
	 'anti-staphylococcal nuclease',
	 'anti-steroids',
	 'anti-streptococcal polyglycerol phosphate',
	 'anti-streptococcal, anti-myosin',
	 'anti-tenascin-C (TN-C) [human]',
	 'anti-testosterone',
	 'anti-testosterone [bovine]',
	 'anti-tetanus',
	 'anti-tetanus toxoid (TeTox)',
	 'anti-thrombospondin',
	 'anti-thymocyte [mouse]',
	 'anti-thyroglobulin',
	 'anti-thyroglobulin [bovine]',
	 'anti-thyroglobulin [human]',
	 'anti-thyroid [human] (stimulating autoantibody)',
	 'anti-thyroid hormone',
	 'anti-thyroid peroxidase (TPO)',
	 'anti-thyroid peroxidase (TPO) (degraded)',
	 'anti-thyroid stimulating hormone (TSH) [human]',
	 'anti-thyrotropin receptor [human]',
	 'anti-tissue factor [human]',
	 'anti-tobacco mosaic virus (TMV)',
	 'anti-transferrin [human]',
	 'anti-transferrin receptor [human]',
	 'anti-transferrin receptor [rat]',
	 'anti-transforming growth factor beta II receptor (TGF-beta-II receptor)',
	 'anti-transition state analogue',
	 'anti-transition state analogue (nitrophenyl phosphonate)',
	 'anti-transmembrane antigen [breast human]',
	 'anti-transmissible gastroenteritis coronavirus (TGEV)',
	 'anti-triazine (herbicide)',
	 'anti-trypsin',
	 'anti-tryptophan synthase beta-subunit [Escherichia coli]',
	 'anti-tumor',
	 'anti-tumor associated antigen',
	 'anti-tumor surface antigen',
	 'anti-tumor surface antigen p185',
	 'anti-tumor-associated glycoprotein-72 (TAG-72)',
	 'anti-tumor-associated glycoprotein-72 (TAG-72) [human]',
	 'anti-type II collagen (CII)',
	 'anti-urokinase plasminogen activator receptor (uPAR) [human]',
	 'anti-urokinase receptor',
	 'anti-varicella-zoster virus (VZV)',
	 'anti-vascular cell adhesion molecule (VCAM) [porcine]',
	 'anti-vascular endothelial growth factor (VEGF)',
	 'anti-vascular endothelial growth factor (VEGF) [human]',
	 'anti-verotoxin II',
	 'anti-vesicular stomatitis virus (VSV)',
	 'anti-viral haemorrhagic septicaemia virus (vhsv)',
	 'anti-von Willebrand Factor (vWF)',
	 'anti-white pine blister rust',
	 'anti-yellow grouper nervous necrosis virus (YGNNV) coat protein',
	 'anti-zearalenone (mycotoxin)']
SPEC.sort()
SPEC = [''] + SPEC