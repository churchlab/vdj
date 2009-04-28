# alignment

import os
import urllib2
import csv
import ClientForm
import vdj
from vdj import refseq
from vdj import vdjexcept
import abacuscore
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet  import IUPAC
import numpy as np
import seqtools

################################
##
## Interface to ABACUS aligner
##
################################

# DEFAULTS.	 CAN BE SET IN abacuscore.set_params()
if os.path.exists('/Users/laserson/'):
	indirectory = '/Users/laserson/research/church/vdj-ome/ref-data/align/'
elif os.path.exists('/home/ul2/'):
	indirectory = '/home/ul2/vdj-ome/ref-data/align/'
outdirectory1=os.getcwd()
vfile='V3.txt'
dfile='D3.txt'
jfile='J3.txt'

class ABACUSer(object):
		
		def __init__(self, refdirectory=indirectory,outputdir=outdirectory1,vref=vfile,dref=dfile,jref=jfile):
				if refdirectory[-1] != '/': refdirectory += '/'
				if outputdir[-1] != '/': outputdir += '/'
				abacuscore.set_params(refdirectory,outputdir,vref,dref,jref)
		
		def alignseq(self,seq):
				return self.alignseqlist([seq])[0]
		
		def alignchainlist(self,chainlist):
			seqlist = [SeqRecord(Seq(chain.seq,IUPAC.unambiguous_dna),description=chain.descr) for chain in chainlist]
			alignedchains = self.alignseqlist(seqlist)
			if len(chainlist) != len(alignedchains):
				raise Exception, "Alignment error."
			for (aligned,chain) in zip(alignedchains,chainlist):
				chain.v = aligned.v
				chain.d = aligned.d
				chain.j = aligned.j
			return
		
		def alignseqlist(self,seqlist):
				# diagnostic
				print "number seqs to align:",len(seqlist)
				
				# put sequences in a temporary file
				tempdir = os.getcwd()
				if tempdir[-1] != '/': tempdir += '/'
				tempfasta = 'tempseqlist'+seqtools.randalphanum(8)+'.fasta'
				while os.path.exists(tempfasta):
					tempfasta = 'tempseqlist'+seqtools.randalphanum(8)+'.fasta'
				op = open(tempdir+tempfasta,'w')
				#seqrecs = False
				if isinstance(seqlist[0],Bio.SeqRecord.SeqRecord):
						#seqrecs = True
						SeqIO.write(seqlist, op, 'fasta')
				else:
						places = int(np.ceil(np.log10(len(seqlist))))
						for i,seq in enumerate(seqlist):
								formatstring = ">%0"+str(places)+"i"
								print >>op, formatstring % i
								print seqtools.seqString(seq)
				op.close()

				# perform the alignment
				temptag = seqtools.randalphanum(8)
				outputfile = 'trialclass'+temptag+'.class'
				while os.path.exists(outputfile):
					temptag = seqtools.randalphanum(8)
					outputfile = 'trialclass'+temptag+'.class'
				abacuscore.fasta2classT2(tempdir,tempfasta,temptag)

				# remove temp input file
				if os.path.exists(tempfasta): os.remove(tempfasta)
				
				# open and parse output
				ip = open(outputfile,'r')
				chains = []
				for k,line in enumerate(ip):
						name = line.split('|')[0].strip()
						rawvdj = line.split('|')[1:]
						v = ''
						d = ''
						j = ''
						for segment in rawvdj:
								if segment.split('*')[0][3] == 'V':
										v = segment.split('*')[0]
								elif segment.split('*')[0][3] == 'D':
										d = segment.split('*')[0]
								elif segment.split('*')[0][3] == 'J':
										j = segment.split('*')[0]
						parsedvdj = ( v, d, j )
						curseq = seqlist[k]
						if curseq.description.strip() != name.strip():
								print curseq.description, 'is different from', name
								raise Exception
						#if seqrecs:
						#		 curseq = seqlist[eval(name.lstrip('0'))]
						#else:
						#		 curseq = seqlist[eval(name.lstrip('0'))]
						chain = vdj.ImmuneChain(seq=seqtools.seqString(curseq), descr = name, vdj = parsedvdj)
						chains.append(chain)
				print "num alignments obtained:",len(chains)
				ip.close()

				# cleanup
				if os.path.exists(outputfile): os.remove(outputfile)
				if os.path.exists('newclass'+temptag+'.class'): os.remove('newclass'+temptag+'.class')

				return chains

		def alignfasta(self,fasta):
				seqlist = seqtools.getFasta(fasta)
				return self.alignseqlist(seqlist)
				

# define aligner for use
abacus = ABACUSer()


def seqlist2positiveseqlist(seqlist):
		'''
		Take a list of SeqRecords and determine whether
		VDJ is on the plus or minus strand.	 If minus, then
		take revcomp of strand and write it to ImmuneChain
		obj, and add a '_revcomp' to the descr (which should
		be just one 'word' with no whitespace)
		'''
		strands = abacuscore.seqlist2positiveseqids(abacuscore.uriseqlist2erezseqlist(seqlist))
		neg = 0
		for (strand,seq) in zip(strands,seqlist):
				if strand == -1:
						seq.description = seq.description.strip() + '_revcomp'
						seq.seq = seq.seq.reverse_complement()
						neg += 1
		print 'Total:', str(len(seqlist)), 'Revcomped:', str(neg)								 
		return seqlist

def chainlist2strandlist(chainlist):
	'''
	Takes a list of ImmuneChains or a Repertoire object, and returns
	a list of +1 or -1 for positive or negative strand at the same
	position
	'''
	seqlist = [SeqRecord(Seq(chain.seq,IUPAC.unambiguous_dna),description=chain.descr) for chain in chainlist]
	strands = abacuscore.seqlist2positiveseqids(abacuscore.uriseqlist2erezseqlist(seqlist))
	return strands


################################
##
## Interface to other aligners
##
################################


class VDJaligner(object):
	'''Object that instantiates a VDJ aligner.	It will interface with
	   multiple alignment methods, chosen at instantiation.
	'''
	def __init__(self, algo='V-QUEST', **kwargs):
		
		self.algo = algo
		self.setAlgo(algo)
	
	def setAlgo(self,algo):
		if algo == 'V-QUEST':
			self.align = self.__VQUESTer
			self.__initVQUESTform()
		elif algo == 'iHMMuneAlign':
			self.align = self.__iHMMuneAligner
		elif algo == 'IgBLAST':
			self.align = self.__IgBLASTer
			self.__initIgBLASTform()

	def __VQUESTer(self,seq):
		'''Aligns sequence.	 If no data, raises exception.	Otherwise, interp
		   is left to the calling function
		'''
		# form already exists in object, with most parameters filled
		# here, we fill out the seq, send it, and parse the results
		
		# submit form
		inputfasta = '>INPUTSEQ\n' + seq
		self.form['l01p01c10'] = inputfasta
		
		request = self.form.click()
		response = urllib2.urlopen( request )
		
		# get data
		for line in response:
			if 'Sequence ID' in line:
				header = line
				rawdata = response.next()
				break
		else:
			raise vdjexcept.NoData("Alignment failed",seq)
		
		header = header.split(';')
		rawdata = rawdata.split(';')
		
		# parse data (no allele)
		data = {}
		data['func'] = rawdata[2]
		data['v'] = rawdata[1].split('*')[0]
		data['d'] = rawdata[6].split('*')[0]
		data['j'] = rawdata[5].split('*')[0]
		#NOTE: cdr3 and junction here refer to the standard IMGT annotations (of amino acids),
		#	   including the C and W nucleotides around the CDR3
		try:
			data['cdr3'] = eval( rawdata[8].split('.')[2][:-1] )
		except NameError, e:
			# no CDR3; probably 'X' in data
			data['cdr3'] = None
		data['junction'] = rawdata[9]
		
		return data
				
	
	def __initVQUESTform(self):
		# get form
		request = urllib2.Request('http://imgt.cines.fr/IMGT_vquest/vquest?livret=0&Option=humanIg')
		response = urllib2.urlopen(request)
		forms = ClientForm.ParseResponse(response,\
										 form_parser_class=ClientForm.XHTMLCompatibleFormParser,\
										 backwards_compat=False)
		response.close()
		form = forms[0]
		
		# fill out base part of form - Synthesis view with no extra options - TEXT
		form['l01p01c03'] = ['inline']
		form['l01p01c07'] = ['2. Synthesis']
		form['l01p01c05'] = ['TEXT'] # may need to be 'TEXT'
		form['l01p01c09'] = ['60']
		form['l01p01c35'] = ['F+ORF+ in-frame P']
		form['l01p01c36'] = ['0']
		form['l01p01c40'] = ['1']	# ['1'] for searching with indels
		form['l01p01c25'] = ['default']
		form['l01p01c37'] = ['default']
		form['l01p01c38'] = ['default']
		form['l01p01c39'] = ['default']
		form['l01p01c27'] = 0
		form['l01p01c28'] = 0
		form['l01p01c29'] = 0
		form['l01p01c30'] = 0
		form['l01p01c31'] = 0
		form['l01p01c32'] = 0
		form['l01p01c33'] = 0
		form['l01p01c34'] = 0
		
		self.form = form
	
	
	def __iHMMuneAligner(self,seq):
		
		# place sequence in file for processing
		prevdir = os.getcwd()
		os.chdir('/Users/laserson/research/church/vdj-ome/align/iHMMuneAlign-full')
		seqf = open('VDJalignSEQtemp.fasta','w')
		print >>seqf, '>INPUTSEQ'
		print >>seqf, seq
		seqf.close()
		
		# call iHMMuneAlign on file
		#os.chdir('/Users/laserson/research/church/vdj-ome/align/iHMMuneAlign-full')
		os.system('perl iHMMuneAlignBATCH_linux_UL20081102.pl VDJalignSEQtemp.fasta')
		
		# read and parse csv format of raw output
		outf = open('VDJalignSEQtemp.fasta.out 1-1.txt','r')
		outfcsv = csv.reader(outf, delimiter=';')
		for rec in outfcsv:
			if rec != [] and rec[0] == 'INPUTSEQ':
				rawdata = rec
		outf.close()
		
		# check for a successful alignment
		if 'NA' in (rawdata[0] + rawdata[1] + rawdata[2]):
			# if the V, D, or J segment has a 'NA' in it
			# perform cleanup and raise the exception
			if os.path.exists('VDJalignSEQtemp.fasta.out 1-1.txt'): os.remove('VDJalignSEQtemp.fasta.out 1-1.txt')
			if os.path.exists('VDJalignSEQtemp.fasta.out 1-1.txt.fostream'): os.remove('VDJalignSEQtemp.fasta.out 1-1.txt.fostream')
			if os.path.exists('gene_nucleotide_mutation_probabilities.txt'): os.remove('gene_nucleotide_mutation_probabilities.txt')
			if os.path.exists('mutations in alignments.txt'): os.remove('mutations in alignments.txt')
			if os.path.exists('VDJalignSEQtemp.fasta'): os.remove('VDJalignSEQtemp.fasta')
			os.chdir(prevdir)
			raise vdjexcept.AlignmentError(seq,"Alignment failed")
		
		# parse output for VDJ info and length of CDR3 (no allele)
		data = {}
		data['id'] = rawdata[0]
		if (rawdata[13] == 'true') and (eval(rawdata[15]) == 0):	# if J in frame and 0 stop codons
			data['func'] = 'Productive'
		else:
			data['func'] = 'Unproductive'
		data['v'] = rawdata[1].split('*')[0]
		data['d'] = rawdata[2].split('*')[0]
		data['j'] = rawdata[3].split('*')[0].split('_')[1]			# iHMMuneAlign includes accession #
		# NOTE: cdr3 and junction refer to different things here:  is it the N1+D+N2 nucl.
		#		and the total length of this junction.	It is not standardized as in IMGT
		data['junction'] = rawdata[5].strip('.') + rawdata[6].strip('.') + rawdata[7].strip('.')
		data['cdr3'] = len(data['junction'])
		
		# delete all intermediate files
		if os.path.exists('VDJalignSEQtemp.fasta.out 1-1.txt'): os.remove('VDJalignSEQtemp.fasta.out 1-1.txt')
		if os.path.exists('VDJalignSEQtemp.fasta.out 1-1.txt.fostream'): os.remove('VDJalignSEQtemp.fasta.out 1-1.txt.fostream')
		if os.path.exists('gene_nucleotide_mutation_probabilities.txt'): os.remove('gene_nucleotide_mutation_probabilities.txt')
		if os.path.exists('mutations in alignments.txt'): os.remove('mutations in alignments.txt')
		if os.path.exists('VDJalignSEQtemp.fasta'): os.remove('VDJalignSEQtemp.fasta')
		os.chdir(prevdir)
		
		return data
	
	
	# NOT FINISHED
	def __IgBLASTer(self,seq):
		# form already exists in object, with most parameters filled
		# here, we fill out the seq, send it, and parse the results
		
		# submit form
		inputfasta = '>INPUTSEQ\n' + seq
		self.form['queryseq'] = inputfasta
		
		request = self.form.click()
		response = urllib2.urlopen( request )
	
		
	def __initIgBLASTform(self):
		request = urllib2.Request('http://www.ncbi.nlm.nih.gov/igblast/')
		response = urllib2.urlopen(request)
		forms = ClientForm.ParseResponse(response,\
										 form_parser_class=ClientForm.XHTMLCompatibleFormParser,\
										 backwards_compat=False)
		response.close()
		form = forms[0]
		
		form['list']   = ['Human']
		form['domain'] = ['2']
		
		self.form = form
	
