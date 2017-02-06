# get pyVCF

from abc import ABCMeta, abstractmethod
import argparse
import re
import logging
import sys
numbers = re.compile('[0-9\.]+')
alt_map = {'X':'X'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

##TODO: clean up start/stop codon scanning
##TODO: implement filters in polymorphism check
##TODO: implement logging
def openFRlist(f):
	#put file into memory, don't read SSD repeatedly
	file = open(f,'rU')
	out = list(file)
	file.close()
	return out

def validGene(seq):
	stopCodons = ['TAG','TAA','TGA']
	startCodons = ['ATG'] # add alternative starts here
	return {'start':seq[0:3].upper() in startCodons, 'stop':seq[-3:].upper() in stopCodons}

def gnCheck(gn,seq):
	if gn.strand == '+':
		seq = ''.join(seq)
		tempseq = seq
		firstcoord = gn.abscoords[0]
		firstpos = firstcoord[0]
		startseq = tempseq
		bpprepend = 0
		while not validGene(startseq)['start']:
		# fix the start codon
			bpprepend += 1
			pos = firstpos-bpprepend
			startseq = gn.contigseq[pos] + startseq
			if bpprepend > 20:
				break	
		if validGene(startseq)['start']:
			tempseq = startseq
			stopseq = tempseq
			lastcoord = gn.abscoords[len(gn.abscoords)-1]
			lastpos = lastcoord[1]
			bpappend = 0
			while not validGene(stopseq)['stop']:
			# fix the stop codon
				bpappend += 1
				pos = lastpos + bpappend
				stopseq += gn.contigseq[pos]
				if bpappend > 20:
					break
			if validGene(stopseq)['stop']:
				tempseq = stopseq
			return tempseq
		else:
			return False
	elif gn.strand == '-':
		seq = "".join(complement[base.upper()] for base in reversed(seq))
		tempseq = seq
		firstcoord = gn.abscoords[len(gn.abscoords)-1]
		firstpos = firstcoord[1]
		startseq = tempseq
		bpprepend = 0
		while not validGene(startseq)['start']:
			bpprepend += 1
			pos = firstpos + bpprepend
			startseq = complement[gn.contigseq[pos]] + startseq
			if bpprepend > 10:
				break
		if validGene(startseq)['start']:
			tempseq = startseq
			stopseq = tempseq
			lastcoord = gn.abscoords[0]
			lastpos = firstcoord[0]
			bpappend = 0
			while not validGene(stopseq)['stop']:
				bpappend += 1
				pos = lastpos - bpappend
				stopseq += complement[gn.contigseq[pos]]
				if bpappend > 10:
					break
			if validGene(stopseq)['stop']:
				tempseq = stopseq
			return tempseq
		else:
			return False

class Polymorphism (object):
	'''Parent class for polymorphism data

	attributes:
		pmFile: list of file contents
		coordinates: list of coordinates of SNPs
		strains: list of strains
		alleles: coordinates:strains:allele
	'''
	__metaclass__ = ABCMeta
	def __init__ (self, file):
		with open(file) as f:
			self.pmFile  = list(f)

	@abstractmethod
	def file_type(self):
		pass

class VCF (Polymorphism):
	'''Variant calling format
	inherits class Polymorphism

	attributes:
		pmFile: list of file contents
	'''
	def __init__(self,file):
		Polymorphism.__init__(self,file)

	def file_type(self):
		return 'VCF'

class SNP (Polymorphism)
	pass
class CDS (object):
	def __init__(self, name, contig, strand):
		self.name = name
		self.contig = contig
		self.abscoords = []
		self.genecoords = []
		self.genepolymorphs = {}
		self.abspolymorphs = {}
		self.strand = strand
		self.refseq = ''
		self.len = 0
		self.contigseq = ''
	def addCoords(self,abscoords,relcoords):
		#coords should be a tuple of tuples.  coord[][0] = position relative to the contig; coord[][1] = position relative to the transcript
		self.abscoords.append(abscoords)
		self.genecoords.append(relcoords)

	def addPolymorphicSite(self,site,seqs): #a dictionary - polymorphic[site][strain]=SNP
		self.genepolymorphs[site[0]] = seqs
		self.abspolymorphs[site[1]] = seqs

	def buildrefSeq(self, contigseq):
		tempseq = ''
		self.contigseq = contigseq
		for locs in self.abscoords:
			for i in range(locs[0],locs[1]):
				if i in self.abspolymorphs.keys():
					tempseq += 'X'
				else:
					tempseq += contigseq[i].lower()
		self.refseq = tempseq

def parseSourceGFF3(gene):
	gene = openFRlist(gene)

	#extract geneName
	firstLine = gene[0]
	geneName = ''.join(firstLine).strip('\n')

	#extract only coding sequence lines
	cdsLst = list(filter(lambda s: 'CDS' in s, gene))

	#get the first coding sequence line for contig,strand
	firstcds = cdsLst[0].split('\t')
	contig = firstcds[0]
	strand = firstcds[6]

	#instantiate CDS 
	gene = CDS(geneName, contig, strand)
	# this could probably be a filter too
	for entry in cdsLst:
		entries = entry.split('\t')
		#split each entry into its components, parse
		coordOne = int(entries[3])
		coordTwo = int(entries[4])
		cdslen = coordTwo-coordOne
		gene.addCoords((coordOne-1,coordTwo),(gene.len,cdslen))
		gene.len += cdslen
	return gene

def polymorphismCheck(pms,cds):
	## put a filter in here
	pms = openFRlist(pms)
	polymorphicSites = {}
	for line in pms:
		if line.startswith('#position'): #find the strains
			strains = line.split()
			strains = strains[1:]

		if numbers.findall(cds.contig)[0] in line: #identify if any SNPs exist in our contig
			snps = line.split() #split the line into <site>, strain 1, strain 2, strain 3...
			site = snps[0].split('_')[2]
			snps = snps[1:] #grab the SNPs

			for i in range(len(cds.abscoords)): 
				abscoords = cds.abscoords[i]
				gncoords = cds.genecoords[i]
				if int(site) in range(abscoords[0],abscoords[1]): #check to see if its in our coordinate range relative to the contig
					relpossite = int(site)-abscoords[0] #relative position within the coordinate range
					genepossite = relpossite + gncoords[0]+i # position within the transcribed gene
					cds.addPolymorphicSite((genepossite, int(site)), {strains[i]:snps[i].upper() for i in range(len(strains))}) #a dictionary - polymorphic[site][strain]=SNP; site is shifted by coords[0]
	return strains

def extractGenomeContig(genome, cds):
	genomeRecords = openFRlist(genome)
	append = False
	contig = ''
	# altered to not use biopython based filtering
	for line in genomeRecords:
		if line.startswith('>'+str(cds.contig)):
			append = True
		elif not line.startswith('>'+ str(cds.contig)) and line.startswith('>supercont'):
			append = False
		elif append and not line.startswith('>supercont'):
			contig += line.strip()
	
	return contig

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Output list of polymorphic FASTA sequences from GFF3 format of one coding sequence, a file of polymorphisms, and the source genome')
	parser.add_argument('-gn',type=str, nargs=1, help = 'source gene', required=True)
	parser.add_argument('-pm',type=str, nargs=1, help = 'polymorphism data', required=True)
	parser.add_argument('-gm',type=str, nargs=1, help = 'genome', required=True)
	parser.add_argument('-o', type=str, nargs=1, help = 'output folder', required=True)
	args = parser.parse_args()
	ofolder = args.o[0]

	if ofolder[-1] != '/': #check for valid folder path
		ofolder += '/'

	gn = parseSourceGFF3(args.gn[0])
	strains = polymorphismCheck(args.pm[0],gn)
	gm = extractGenomeContig(args.gm[0],gn)
	buf = ''

	built = gn.buildrefSeq(gm)
	# log = open(ofolder+'log.txt','a')
	# if built:
	# 	log.write(gn.name+' built gene successfully \n')
	# else:
	# 	log.write(gn.name + ' FAILED, start codon: ' + built[0] + ' - stop codon: ' + built[1])
	# log.close()

	for strain in strains:
	 	#write the name of each strain
		buf+=('>'+strain+'_'+gn.name + '\n')
		gene = list(gn.refseq)
		for i in range(len(gene)):
				if gene[i] == 'X':
					gene[i]=gn.genepolymorphs[i][strain]

		gene = gnCheck(gn,gene)
		if gene:
			buf += gene +'\n'
		else:
			buf += "gene failed"


	buf = buf[:-1]

	output = open(ofolder+gn.name+gn.strand+'_polymorphs.txt','w')
	output.write(buf)
	output.close()

