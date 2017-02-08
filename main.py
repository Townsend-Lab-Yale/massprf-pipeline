import pprint
import gffutils
import vcf
from Bio import SeqIO, Sequence
from pathlib import Path
import itertools
import csv
'''

top priority: unify indexing systems to linear sequences
1) implement full gff3 initialization of Genome class
2) implement FASTA coordinate/chromosome referencing to create CDS
5) implement genome wide composite class of genes
7) think about class hierarchy and inheritance
8) figure out how to make the adaptors & builders singletons
9) figure out an optimal time to export files so variants do not have to be reloaded continuously
10) make variants class iterable
11) add scaling algorithm
12) add subprocess spawning for MUSCLE, massprf
13) implement reverse complementation
14) implement codon scanning
'''

complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def allele(gt):
	'''OR gate heterozygous/homozygous dominant'''
	return gt in ['0/1','1/1']
class TwoWayDict(dict):
	'''TwoWayDict borrowed from http://stackoverflow.com/questions/1456373/two-way-reverse-map
	This mapping allows bidirectional dictionary type. used by homology map'''
	def __setitem__(self, key, value):
		if key in self:
			del self[key]
		if value in self:
			del self[value]
		dict.__setitem__(self, key, value)
		dict.__setitem__(self, value, key)

	def __delitem__(self, key):
		dict.__delitem__(self, self[key])
		dict.__delitem__(self, key)

	def __len__(self):
		return dict.__len__(self) // 2

class CDSAdaptorBuilder(object):
	'''CDSAdaptorBuilderis a factory for coding sequence adaptors classes
		To extend functionality to new file types, write a new class of Adaptor that creates objects of class Gene and accepts **kwargs
		return the new class using an elif statement in __new__ of CDSAdaptorFactory'''
	def __new__(self,cls, **kwargs):
		if cls == 'gff3':
			return GffUtilAdaptor(**kwargs)

class GffUtilAdaptor(object):
	'''to do:

	1) fix **kwargs checking to account for all scenarios - consider using paper logic diagram first
		- document logic as right now it is a bit obtuse
	2) fix createGene to initalize Gene classes
	3) integrate with Genome?
	'''
	def __init__(self, **kwargs):
		if not kwargs:
			raise IOError("No gffutils database or gff3 supplied")
		self.dbname = ''
		self.gff3 = ''
		if 'gff3' in kwargs:
			print(kwargs['gff3'])
			self.gff3 = Path(kwargs['gff3'])
			if not 'db' in kwargs:
				self.dbname = Path(str(kwargs['gff3']) + 'db')
		if 'db' in kwargs:
			print(kwargs['db'])
			if Path(kwargs['db']).is_file():
				self.dbname = Path(kwargs['db'])

		if not self.dbname.is_file():
			self.db = gffutils.create_db(str(self.gff3), str(self.dbname))
		self.db = gffutils.FeatureDB(str(self.dbname))

	def getGenes(self):
		for feature in self.db.all_features():
			if feature.featuretype=='gene':
				yield feature

	def createGene(self, feature, *kwargs):
		if not isinstance(feature, gffutils.feature.Feature):
			raise AttributeError("passed feature is not of type gffutils.feature.Feature")

class VariantsAdaptorBuilder(object):
	'''Interpret passed variant format file and return properly interfaced variant file'''
	def __new__(self,frmt, variants):
		if frmt == 'vcf':
			return PyVCFAdaptor(variants)

class PyVCFAdaptor(object):
	'''adapts VCF files into Variants class of Variant leafs
	skips all non-snps for now.  Fix!'''
	def __new__(self, variants):
		if not Path(variants).is_file():
			raise IOError("Invalid filename supplied for VCF")
		else:
			self.filename = variants
			self.file = open(self.filename)
			self.vcf = vcf.Reader(self.file)
			
			strains = self.vcf.samples
			variants = Variants(strains)
			for record in self.vcf:
				calls = []
				if record.is_snp:
					coordinate = Coordinate(record.CHROM, record.POS)
					for strain in variants.strains:
						calls.append(allele(record.genotype(strain)['GT']))

					ref = record.REF
					alt = record.ALT[0] #See above - SNP reading breaks >1 alternate allele

					alleles = [alt if gt else ref for gt in calls]

					for a, s in zip(alleles, variants.strains):
						variants.addVariant(Variant(coordinate, s, a))
			self.file.close()
		return variants

class GenomeAdaptorBuilder(object):
	'''similar to other AdaptorBuilders, this one builds genome adaptors based on passed formatting
	add support **kwargs to add variants'''
	def __new__(self, reference, species, frmt = "FASTA", variants = None, **kwargs):
		if not Path(reference).is_file():
			raise AttributeError("no valid reference genome supplied")
		if not species:
			raise AttributeError("No species supplied")
		if not frmt:
			raise IOError("No parse format supplied")
		if not strain:
			raise AttributeError("No strain supplied")
		if frmt == 'placeholder':
			#add new formats for parsing here
			pass
		elif frmt == "FASTA":
			#default to FASTA format for genome
			return FastaGenomeAdaptor(reference, species, variants)
		else:
			pass
class HomologAdaptorBuilder(object);
	'''Adapts input homolog formats to HomologAdaptor'''
	def __new__(self, homologyfile, frmt = "CSVa"):
		if not Path(homologyfile).is_file():
			raise AttributeError("No valid homology file specified")
		if frmt == "CSVa":
			#this file type is useful when you have two distinct divergent species with a CSV file mapping of the homologs
			return CSVaHomologs(homologyfile)
		elif frmt == "CSVb":
			return CSVbHomologs(homologyfile)
class CSVHomologs(object):
	def __init__(self, homologyfile):
		self.csvfilename = homologyfile
		self.csvfilehandle = open(self.csvfilename)
		self.reader = csv.reader(csvfilehandle)
	def close(self):
		self.csvfilehandle.close()
class CSVaHomologs(CSVHomologs):
	'''accepts CSV homology files where row = homolog, column = species'''
	def __init__(self):
		pass
class CSVbHomologs(CSVHomologs):
	'''this file format has two groupings of species/strains, and the species are using the same reference genome'''
	def __init__(self):
		pass

class HomologyMap(object):
	def __init__(self, species, homologs):
		pass

class FastaGenomeAdaptor(object):
	'''FastaGenomeAdaptor adapts biopython SeqIO genome into the genome class
	"species" here is a generic term, eg rice, banana, etc, that is useful for bluntly differentiating objects of type genome
	'''
	def __new__(self, reference, species, variants = None):
		#GenomeAdaptorBuilder already did the type checking of reference to insure it is a filetype
		self.reference_filename = reference
		self.reference_file = open(reference)

		#if variants are passed, we are going to build a bunch of new genome objects and return the genomes as a list
		#the genomes should build & insert their variants themselves. 

class Genome(object):
	''' in truth this should be a polymorphic class that can either be the fasta sequence or can be gff based.  the final output should combine both, and contain a list of Gene objects'''
	def __init__(self, species, strain, reference):

		self.reference = reference
		self.species = species
		self.strain = strain

	def __repr__(self):
		return repr(species, strain)
class Features(Genome);
	''' subclass of Genome
		features within a genome
	'''
	pass
class GenomeSequence(Genome):
	'''subclass of Genome
		genome w/ reference sequence'''
	pass

class Chromosome(Bio.Seq.Seq):
	pass

class CodingSequence(Bio.Seq.Seq):
	pass

class Variants(object):
	'''Variants, a composite class of Variant leafs
		attributes:
			strains: the strains that have been called
			locations: the Coordinate locations of each variant
			_snps_by_loc: private dictionary of dict[Coordinate] : [Variant,Variant,Variant]
			_snps_by_strain: private dictionary of dict[strain] : dict[chromosome]: [Variant, Variant, Variant]
		methods:
			initialization: provide strains or raise exception;
				initialize all attributes to empty
			addVariant: Given a variant of type Variant, 
				1) add its coordinate to self.locations
				2) add its coordinate to the keys of self._snps_by_loc
				3) append variants to that list
				4) append variants to _snps_by_strain[variant.strain]
			of_strain: 
				arguments: string(strain), string(chromosome)
				returns list of Variants of strain
			on_chromosome:
				returns hits, dict[hits] of variants on a chromosome
			in_range:
				returns variants on a chromosome within a range

	'''
	def __init__(self, strains):
		if not strains:
			raise AttributeError("No strains provided")
		self.strains = strains
		self.locations = []
		self._snps_by_loc = {}
		self._snps_by_strain = {strain:[] for strain in self.strains}

	def addVariant(self, variant):
		if variant.coordinate not in self.locations:
			self.locations.append(variant.coordinate)
		if variant.coordinate not in self._snps_by_loc.keys():
			self._snps_by_loc[variant.coordinate]=[]
		self._snps_by_loc[variant.coordinate].append(variant)
		self._snps_by_strain[variant.strain].append(variant)

	def of_strain(self, strain, chromosome = False):
		if strain not in self.strains:
			raise ValueError("strain not in Variants")
		else:
			if chromosome:
				return filter(lambda l: l.chromosome == str(chromosome), self._snps_by_strain[strain])
			else:
				return self._snps_by_strain[strain]

	def on_chromosome(self, chromosome):
		hits = filter(lambda l: l.chromosome == str(chromosome), self.locations)
		variants = {}
		for n in hits:
			variants[n] = self._snps_by_loc[n]
		return (hits, variants)

	def in_range(self, chromosome, x1,x2):
		if x2 <= x1:
			x2, x1 = x1, x2
		hits, variants = self.on_chromosome(chromosome)
		in_range = filter(lambda x: x1 <= x.start <= x2, hits)

		return (list(in_range), {x:variants[x] for x in in_range})
		
class Variant(object):
	def __init__(self, coordinate, strain, gt):
		if not strain:
			raise AttributeError("no strains supplied to variants")
		if not isinstance(coordinate, Coordinate):
			raise AttributeError("invalid coordinate supplied")

		self.coordinate = coordinate
		self.strain = strain
		self.gt = gt

	def __repr__(self):
		return repr((self.strain, self.coordinate, self.gt))

class Coordinate(object):
	'''Coordinate leaf
	chromosome is stored as a string for now to facilitate cases where the chromosome is oddly named, such as "supercont12.1"
	'''
	def __init__(self, chromosome, start, stop = False):
		if not chromosome:
			return AttributeError("No chromosome supplied for coordinate")
		if not start:
			return AttributeError("No start position supplied for coordinate")
		self.chromosome = chromosome
		self.start = int(start)
		if not stop:
			self.stop = start
		else:
			self.stop = stop
	def __repr__(self):
		return repr((self.chromosome, self.start))

class Gene(object):
	def __init__(self, coordinates, strand):
		pass
		'''else:
			self.geneid = self.gff3feature.id
			self.features = self.db.children(gff3feature)
			self.CDS = self.db.children(self.gff3feature,featuretype=='CDS')
#implement CDS coords

			self._strand = self.gff3feature.strand

	def __len__(self):
		pass 
	@property
	def strand(self):
		return self.strand
	@strand.setter
	def strand(self, watsoncrick):
		if watsoncrick in ['-','negative','minus','antisense','-1',-1,'0',0]: #default to positive strand
			self._strand = '-'
		elif watsoncrick is '.':
			self._strand = 'unknown'
		else:
			self._strand = '+'

	#implement CDS coordinate generator
'''

#filename = '../ricemassprf/Oryza_sativa.IRGSP-1.0.24.gtf.mRNA.clean.gff3'
#featuredb = '../ricemassprf/ricefeatureDB'

#kwargs = {"gff3":filename, "db": featuredb}

#x = CDSAdaptorBuilder("gff3",**kwargs)
#print(type(x))

