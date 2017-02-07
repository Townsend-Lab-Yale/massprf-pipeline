import pprint
import gffutils
import vcf
from Bio import SeqIO
from pathlib import Path

'''1) implement full gff3 initialization of Gene class
2) implement FASTA coordinate/chromosome referencing to create CDS
5) implement genome wide composite class of genes
7) think about class hierarchy and inheritance
8) figure out how to make the adaptors & builders singletons
9) figure out an optimal time to export files so variants do not have to be reloaded continuously
10) make variants class iterable
'''
def allele(gt):
	'''OR gate heterozygous/homozygous dominant'''
	return gt in ['0/1','1/1']

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
			self.vcf = vcf.Reader(open(variants))

		'''initialize the first row of the variants array - figure out the best way to do this initialization w/in one loop through entire file'''
		calls = []		
		firstrecord = next(self.vcf)
		strains = [s.sample for s in firstrecord.samples]
		variants = Variants(strains)
		if firstrecord.is_snp:
			coordinate = Coordinate(firstrecord.CHROM, firstrecord.POS)
			for strain in variants.strains:
				calls.append(allele(firstrecord.genotype(strain)['GT']))
		
			ref = firstrecord.REF
			alt = firstrecord.ALT[0] #What happens if there is more than one ALT? SNP reading breaks

			calls = [alt if gt else ref for gt in calls]

			for a, s in zip(calls, variants.strains):
				variants.addVariant(Variant(coordinate,s,a))
		
		# repeat the above for the rest of the variants
		for record in self.vcf:
			if record.is_snp:
				calls = []
				coordinate = Coordinate(record.CHROM, record.POS)
				for strain in variants.strains:
					calls.append(allele(record.genotype(strain)['GT']))

				ref = record.REF
				alt = record.ALT[0] #See above - SNP reading breaks >1 alternate allele

				calls = [alt if gt else ref for gt in calls]

				for a, s in zip(calls, variants.strains):
					variants.addVariant(Variant(coordinate, s, a))

		return variants

class GenomeAdaptorBuilder(object):
	'''similar to other AdaptorBuilders, this one builds genome adaptors based on passed formatting
	add support **kwargs to add variants'''
	def __new__(self, frmt, species, strain, reference, variants = None, **kwargs):
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
		else:
			#default to FASTA format for genome
			return Genome(species, strain, SeqIO.to_dict(SeqIO.parse(reference)), variants)

class Coordinate(object):
		def __init__(self, chromosome, start, end = False):
			if not chromosome:
				return AttributeError("No chromosome supplied for coordinate")
			if not start:
				return AttributeError("No start position supplied for coordinate")
			self.chromosome = chromosome
			self.start = start
			if not end:
				self.end = start
			else:
				self.end = end

class Genome(object):
	''' in truth this should be a polymorphic class that can either be the fasta sequence or can be gff based.  the final output should combine both, and contain a list of Gene objects'''
	def __init__(self, species, strain, reference, vcf = None):
		if not bool(reference):
			raise AttributeError("Error creating genome; valid file but empty dictionary")

		self.reference = reference
		self.species = species
		self.strain = strain

		self.chromosomes = len(reference)
class Variants(object):
	'''figure out how to cause this class to redirect instantiation'''
	def __init__(self, strains):
		self.strains = strains
		self.locations = []
		self.snps = {}

	def addVariant(self, variant):
		if variant.coordinate not in self.locations:
			self.locations.append(variant.coordinate)
		if not self.snps[variant.coordinate]:
			self.snps[variant.coordinate]=[]
		self.snps[variant.coordinate].append(variant.call)

		
class Variant(object):
	def __init__(self, coordinate, strain, call):
		if not strains:
			raise AttributeError("no strains supplied to variants")
		
		self.coordinate = coordinate
		self.strains = strains
		self.gt = call
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

filename = '../ricemassprf/Oryza_sativa.IRGSP-1.0.24.gtf.mRNA.clean.gff3'
featuredb = '../ricemassprf/ricefeatureDB'

kwargs = {"gff3":filename, "db": featuredb}

x = CDSAdaptorBuilder("gff3",**kwargs)
print(type(x))

