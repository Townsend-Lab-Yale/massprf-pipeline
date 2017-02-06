import pprint
import gffutils
from pathlib import Path
'''1) implement full gff3 initialization of Gene class
2) implement FASTA coordinate/chromosome referencing to create CDS
3) implement variant calling factory
4) implement variant calling format inserter
5) implement genome wide composite class of genes
6) implement coordinate class
7) think about class hierarchy and inheritance
'''
class GeneCoordinates(object):
		def __init__(self, start, end):
			pass

class CDSAdaptorFactory(object):
	'''CDSFactory is a factory for coding sequence adaptors classes
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

x = CDSAdaptorFactory("gff3",**kwargs)
print(type(x))

