# clones.py
# Methods to deal with sequencing clones (incl clustering)

import vdj
import cluster
import os
from seqtools import editdist
import seqtools


# ====================
# = Distance metrics =
# ====================

def chain_Levenshtein(x,y): return editdist(x.seq,y.seq)

# DEPRECATE
def chain_NGLD(x,y):
	'''
	Normalized Generalized Levenshtein Distance
	(generalization is trivial case here; alpha=1)
	based on IEEE Trans Pattern Analys Mach Intel 29(6):1091
	'''
	GLD = editdist(x.junction,y.junction)
	Nx	= float(len(x.junction))
	Ny	= float(len(y.junction))
	return (2.*GLD)/(Nx+Ny+GLD)

def NGLD(x,y):
	'''
	Normalized Generalized Levenshtein Distance
	(generalization is trivial case here; alpha=1)
	based on IEEE Trans Pattern Analys Mach Intel 29(6):1091
	'''
	try:
		GLD = editdist(x,y)
		Nx	= float(len(x))
		Ny	= float(len(y))
		return (2.*GLD)/(Nx+Ny+GLD)
	except ZeroDivisionError, e:
		return 0.

# ==============
# = Clustering =
# ==============

def repertoireMapper(repertoire,dist_fn):
	return lambda n,m: dist_fn(repertoire[int(n)].junction,repertoire[int(m)].junction)

def clusterRepertoire(repertoire, linkage='single'):
	clusters = []
	trees = {}
	for Vgene in vdj.refseq.IGHV:
		trees[Vgene] = {}
		for Jgene in vdj.refseq.IGHJ:
			trees[Vgene][Jgene] = cluster.HierarchicalClustering( repertoire.get_idxs_AND(Vgene,Jgene), repertoireMapper(repertoire,NGLD), linkage=linkage )
			clusters += trees[Vgene][Jgene].getlevel(0.05)
		
	return (clusters,trees)

# USE OF hcluster:
#Xrep14 = rep[vdj.refseq.IGHV[1],vdj.refseq.IGHJ[4]]
#Yrep14 = hcluster.pdist(np.arange(len(Xrep14)).reshape((len(Xrep14),1)), vdj.clones.repertoireMapper(Xrep14,vdj.clones.NGLD) )
#Zrep14 = hcluster.linkage(Yrep14,'single')
#Trep14 = hcluster.fcluster(Zrep14,0.8)


# ==============
# = CDR3 manip =
# ==============

def extractCDR3(chains):
	'''
	takes list of ImmuneChain objects.	They must have a sequence, a V, and J region alignment
	'''
	
	if isinstance(chains,vdj.Repertoire): chain_list = chains.chains
	if isinstance(chains,list): chain_list = chains
	if isinstance(chains,vdj.ImmuneChain): chain_list = [chains]

	curdir = os.getcwd()
	if curdir[-1] != '/': curdir += '/'	   
	
	query_file_name = 'query'+seqtools.randalphanum(8)+'.fasta'
	while os.path.exists(query_file_name):
		query_file_name = 'query'+seqtools.randalphanum(8)+'.fasta'
	
	ref_FR3_file_name = 'ref_FR3_'+seqtools.randalphanum(8)+'.fasta'
	while os.path.exists(ref_FR3_file_name):
		ref_FR3_file_name = 'ref_FR3_'+seqtools.randalphanum(8)+'.fasta'
	
	ref_J_file_name	  = 'ref_J_'+seqtools.randalphanum(8)+'.fasta'
	while os.path.exists(ref_J_file_name):
		ref_J_file_name	  = 'ref_J_'+seqtools.randalphanum(8)+'.fasta'

	cdr3s = []
	for chain in chain_list:
		try:
			# init indexes
			v_end_idx = len(chain.seq)
			j_start_idx = 0
		
			# can't get the CDR3
			if chain.v == '' or chain.j == '':
				cdr3s.append('')
				continue
	
			# write sequence and references to fasta file
			query_file	 = open(curdir+query_file_name, 'w')
			ref_FR3_file = open(curdir+ref_FR3_file_name,'w')
			ref_J_file	 = open(curdir+ref_J_file_name, 'w')
	
			print >>query_file, '>query_V'
			print >>query_file, chain.seq
	
			print >>ref_FR3_file, '>ref_V_FR3'
			print >>ref_FR3_file, vdj.refseq.IGHV_FR3[chain.v]
	
			print >>ref_J_file, '>ref_J_FR4'
			print >>ref_J_file, vdj.refseq.IGHJ_FR4[chain.j][1]
	
			query_file.close()
			ref_FR3_file.close()
			ref_J_file.close()
	
			v_align = os.popen('bl2seq -i ' + query_file_name + ' -j ' + ref_FR3_file_name + ' -p blastn -g T -F F -e 1.0 -D 1 -r 3 -G 3 -E 6 -W 9')
	
			for line in v_align:
				if line[0]=='#':
					continue
				linedata = line.split()
			
				refend   = eval(linedata[-3])
				reflen   = len(vdj.refseq.IGHV_FR3[chain.v])
				ext_deficit = reflen - refend
				
				# TODO: sometimes bl2seq will not extend to the end of FR3
				# to resolve this, I need to do the same thing I do with the
				# J region, and move the 2nd-CYS to an internal site.  the prob I
				# forsee here is that there is not much add'l conserved seq
				# after the 2nd-CYS
				if ext_deficit > 3:
					raise Exception, 'V segment FR3 alignment for sequence ' + chain.descr + ' was not extended through the end'
				queryend = eval(linedata[-5])
				v_end_idx = eval(line.split()[-5]) + ext_deficit - 3		# FR3 includes the CYS
				break
	
			query_file	 = open(curdir+query_file_name, 'w')
	
			print >>query_file, '>query_J'
			print >>query_file, chain.seq[v_end_idx:]
	
			query_file.close()
	
			j_align = os.popen('bl2seq -i ' + query_file_name + ' -j ' + ref_J_file_name + ' -p blastn -g T -F F -e 1.0 -D 0 -r 3 -G 3 -E 6 -W 9')
	
			# as the conserved TRP is internal to the ref sequence, we have to deal with possible gapped alignments
			for line in j_align:
				linedata = line.split()
				if linedata == [] or linedata[0] != 'Query:':
					continue
		
				querystart = eval(linedata[1])
				querymatch = linedata[2]
		
				j_align.next()	# burn a line
				line = j_align.next()
				linedata = line.split()
		
				refstart   = eval(linedata[1])
				refmatch   = linedata[2]
		
				# count the number of gaps before the position of the J-TRP
				TRPstart = vdj.refseq.IGHJ_FR4[chain.j][0] - (refstart - 1)	# first candidate pos
				if TRPstart < 0 or TRPstart > len(refmatch):
					raise Exception, 'J segment FR4 alignment for sequence ' + chain.descr + ' did not include the J-TRP site'
				refgaps  = refmatch[:TRPstart].count('-')	# count gaps up to putative TRP pos
				seengaps = 0
				while refgaps != 0:	# iteratively find all gaps up to TRP
					seengaps += refgaps
					TRPstart += refgaps	# adjust it for any gaps in ref seq
					refgaps   = refmatch[:TRPstart].count('-') - seengaps	# any add'l gaps?
		
				querygaps = querymatch[:TRPstart].count('-')	# num gaps in corres query segment
		
				# j_start_idx = pos i left off for V + distance into what's left + distance to TRP - adj for gaps + 3 for TRP
				j_start_idx = v_end_idx + (querystart - 1) + TRPstart - querygaps + 3
			
				break	# in case there are less high-scoring matches
			
			cdr3s.append(chain.seq[v_end_idx:j_start_idx])
		except Exception, e:
			print e
			cdr3s.append('')
			
	
	if os.path.exists(curdir+query_file_name): os.remove(curdir+query_file_name)
	if os.path.exists(curdir+ref_FR3_file_name): os.remove(curdir+ref_FR3_file_name)
	if os.path.exists(curdir+ref_J_file_name): os.remove(curdir+ref_J_file_name)

	if isinstance(chains,vdj.Repertoire):
		if len(chains.chains) != len(cdr3s):
			raise Exception, 'wrong number of CDR3s extracted'
		for (chain,cdr3) in zip(chains.chains,cdr3s):
			chain.junction = cdr3
			chain.cdr3 = len(cdr3)
	if isinstance(chains,list):
		if len(chains) != len(cdr3s):
			raise Exception, 'wrong number of CDR3s extracted'
		for (chain,cdr3) in zip(chains,cdr3s):
			chain.junction = cdr3
			chain.cdr3 = len(cdr3)
	if isinstance(chains,vdj.ImmuneChain):
		if len(cdr3s) != 1:
			raise Exception, 'wrong number of CDR3s extracted'
		chains.junction = cdr3s[0]
		chains.cdr3 = len(cdr3s[0])
	
	return chains



# def new_dict_partition():
# 	'''
# 	returns a dictionary that indexes by V region, each elt of which
# 	is a dictionary which indexes by J, each elt of which is an
# 	empty list
# 	'''
# 	partitioned = {}
# 	for Vgene in vdj.refseq.IGHV+['']:	# [''] is for no alignment
# 		partitioned[Vgene] = {}
# 		for Jgene in vdj.refseq.IGHJ+['']:
# 			partitioned[Vgene][Jgene]=[]
# 	return partitioned
# 
# 
# 
# 
# def partition_chains( chains ):
# 	'''
# 	take a list of ImmuneChain objects, and return a dict
# 	that first indexes by V region and then by J region,
# 	that contains all the assoc chains.	 There is also an
# 	unclassified bin
# 	'''
# 	# define data structure
# 	partitioned = new_dict_partition()
# 	
# 	# bin the objects
# 	for chain in chains:
# 		v = chain.v
# 		
# 		partitioned[chain.v][chain.j].append(chain)
# 	
# 	return partitioned
# 
# 
# 
# def cluster_partitioned_chains( partitioned_chains, linkage='single' ):
# 	'''
# 	return a dict partition obj (see above), where each
# 	elt is a cluster.HierarchicalClustering object for the
# 	data that needs to be clustered in that partition
# 	'''
# 	cluster_objs = new_dict_partition()
# 	for Vgene in vdj.refseq.IGHV+['']:
# 		for Jgene in vdj.refseq.IGHJ+['']:
# 			cl = cluster.HierarchicalClustering( partitioned_chains[Vgene][Jgene], chain_NGLD, linkage=linkage )
# 			cluster_objs[Vgene][Jgene] = cl
# 			
# 	return cluster_objs
