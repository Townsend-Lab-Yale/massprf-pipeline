import pprint
import gffutils
import vcf
import subprocess
from functools import reduce
from Bio import SeqIO, Seq, SeqRecord
from Bio.Alphabet import generic_dna
from pathlib import Path
import itertools
import csv


'''

top priority: unify indexing systems to linear sequences
1) implement full gff3 initialization of Genome/Exon classes (UGH) (below 2,5 are same tasks)
    2) implement FASTA coordinate/chromosome referencing to create CDS
    5) implement genome wide composite class of genes (02/10/17 half way done, fasta reference to chromosome is done, add method to insert features)
7) think about class hierarchy and inheritance (see bio.seq.seq inheritance, dict inheritance,csvhomologs inheritance, genome inheritance)
8) figure out how to make the adaptors & builders singletons (02/10/17 this is probably not super necessary but would be nice; however appears to require making a new super class for all builders)
9) figure out an optimal time to export files so variants do not have to be reloaded continuously
10) make variants class iterable (02/10/17 is this still necessary? - .of_strain() method seems to handle all that is demanded
11) add scaling algorithm (02/10/17 scaler.py)
12) add subprocess spawning for MUSCLE, massprf (02/10/17 see aligntrim.py)
13) implement reverse complementation (02/10/17 by inheriting Bio.Seq.Seq, Chromosomes and CodingSequences should be capable of reverse_complement)
14) implement codon scanning (02/10/17 this can get messy really quickly - reconsider implementation; 
        if desire to include, see extractPoly.py)
'''

'''constants definitions'''

'''tool functions'''

def allele(gt):
    '''OR gate heterozygous/homozygous dominant'''
    return gt in ['0/1','1/1']


''' begin fundamental data structure definitions'''

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

class HomologyMap(object):
    def __init__(self, species, homologs):
        pass

class Genome(object):
    '''Parent class for ReferenceGenome and VariantGenome

    may change representation to abstract class, not sure how this would work
    '''
    def __init__(self, species, numchromosomes):
        self.species = species
        self.numchromosomes = numchromosomes
        self.CDS = {}
        self.chromosomes = {n:'' for n in range(1, numchromosomes+1)}

    def __repr__(self):
        return repr(self.species + "_" + self.strain)

    def __str__(self):
        return str(self.species + "_" + self.strain)

    def getChromosomes(self):
        for n in range(1, self.numchromosomes+1):
            yield(self.chromosomes[n])

    def addCDS(self, gene):
        if not isinstance(gene, CodingSequence):
            raise AttributeError("passed gene is not of type CodingSequence")
        self.CDS[gene.name] = gene

    def write(self, directory, filename = None):
        if not filename:
            filename = str(self)
        if str(directory):
            SeqIO.write([SeqRecord.SeqRecord(chromo, id=chromo.id, description=chromo.description) for chromo in self.getChromosomes()], str(directory) + filename, "fasta")
        
class ReferenceGenome(Genome):
    '''ReferenceGenome:
        extends Genome

        attributes:
            species: a blunt grouping name useful
            numchromosomes: the total number of chromosomes in the genome
            chromosomes: a dictionary mapping each chromosome number to a Chromosome object that holds the sequence
            substrains: dictionary mapping of VariantGenome's

        methods:
            __init__:
                initialize with no substrains, set strain = "reference"
                initialize numchromosomes
                initialize Chromosomes dictionary
            addStrain(strain):
                type(strain) = VariantGenome
                add a variant strain to the substrains dictionary
    '''
    def __init__(self, species, referenceChromosomes):
        #referenceChromosomes is a list of chromosomes w/ sequences
        Genome.__init__(self, species, max(referenceChromosomes))
        self.substrains = {}
        self.strain = "reference"
        for n in range(1, self.numchromosomes+1):
            self.chromosomes[n] = Chromosome(self, n, referenceChromosomes[n])

    def addStrain(self, strain):
        if not isinstance(strain, VariantGenome):
            raise AttributeError("error adding variantstrain, improper variantGenome supplied")
        self.substrains[strain.strain] = strain

class VariantGenome(Genome):
    '''VariantGenome
    extends: Genome
    
    attributes:
        species: a blunt grouping name useful
        numchromosomes: the total number of chromosomes in the genome
        chromosomes: a dictionary mapping each chromosome number to a Chromosome object that holds the sequence
        strain: substrains of species grouping
        variants: a list of Variant s, each with a chromosome number and index adjusted position.

    methods:
        __init__: initialized with a ReferenceGenome object, a strain name, and a list of Variant s
            1) initializes chromosomes to ReferenceGenome Chromosomes
            2) iterates over snps (Variants) and inserts the SNPs at corresponding positions in the strings held by
                 the chromosome dictionary
    '''
    def __init__(self, reference, strain, variants):
        Genome.__init__(self, reference.species, reference.numchromosomes)
        self.reference = reference
        self.strain = strain
        self.variants = variants

        self.chromosomes = reference.chromosomes

        for snp in self.variants:
            self.chromosomes[snp.chromosome].sequence[snp.pos] = snp.gt

class Chromosome(Seq.Seq):
    '''Chromosome
    extends Bio.Seq.Seq
    See biopython documentation for full documentation
    (new) attributes:
        reference: parent genome, can be a ReferenceGenome or VariantGenome
        number: chromosome number, cast as int
        sequence: str sequence of chromosome extracted 
    '''
    def __init__(self, reference, number, sequence):
        Seq.Seq.__init__(self, sequence, generic_dna)
        self.genome = reference
        self.number = int(number)
        self.id = str(self.number)
        self._description = str(self.genome) + " chromosome " + str(self.id)

    @property
    def description(self):
        return self._description

class CodingSequence(object):

    '''PLEASE make sure to adopt 0 indexing of coordinates from gff3'''
    def __init__(self, reference, name, coordinates, strand, homolog = None, gene_id = None):
        self.reference = reference
        self.species = reference.species
        self.genename = name
        self.coordinates = coordinates
        self.chromosome = coordinates[0].chromosome
        self.strand = strand
        self.homolog = homolog
        if gene_id:
            self.gene_id = gene_id
        else:
            self.gene_id = name
        self.length = len(self)

    def __repr__(self):
        return repr(str(self))

    def __str__(self):
        return str(self.species) + ' ' + str(self.genename)

    def __len__(self):
        return reduce(lambda x, y: x+y, map(lambda x: len(x), self.coordinates))

    def getSequence(self, strain):
        curstrain = self.reference.substrains[strain]
        chromosome = curstrain.chromosomes[self.chromosome]
        sequence = Seq(''.join([chromosome[coordinate.pos[0]:coordinate.pos[1]] for coordinate in self.coordinates]),generic_dna)
        if self.strand is '-':
            return sequence.reverse_complement()
        else:
            return str(sequence)



class Variants(object):s
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
                return self._snps_by_strain[strain].sort(key=lambda x: x.chromosome).sort(key=lambda x: x.pos)

    def on_chromosome(self, chromosome):
        hits = list(filter(lambda l: l.chromosome == str(chromosome), self.locations))
        variants = {}
        for n in hits:
            variants[n] = self._snps_by_loc[n]
        return (hits, variants)

    def in_range(self, chromosome, x1,x2):
        if x2 <= x1:
            x2, x1 = x1, x2
        hits, variants = self.on_chromosome(chromosome)
        in_range = filter(lambda x: x1 <= x.pos <= x2, hits)

        return (list(in_range), {x:variants[x] for x in in_range})
        
class Variant(object):
    '''Variant

    SNPs
    attributes:
        self.coordinate: type(coordinate) = Coordinate; 0 indexed position on 1-indexed chromosome
        self.strain: string representation of strain
        self.gt: string representation of allele

    methods:
        __init__: 
            type checking of strain, coordinate
            init coordinate, strain, gt
        __repr__:
            internal representation of Variant
            returns repr((self.strain, self.coordinate, self.gt))
    '''
    def __init__(self, coordinate, strain, gt, ref):
        if not strain:
            raise AttributeError("no strains supplied to variants")
        if not isinstance(coordinate, Coordinate):
            raise AttributeError("invalid coordinate supplied")
        self.coordinate = coordinate
        self.strain = strain
        self.gt = gt
        self.ref = ref

    def __repr__(self):
        return repr((self.strain, self.coordinate, self.gt))

    def ispolymorphic(self):
        return self.gt != self.ref

class Coordinate(object):
    '''
    Coordinate

    attributes:
        self.chromosome: integer representation of chromosome. 1 - indexed
        self.start: 0-indexed single identifier position - exchangeable for absolute position in the instance of no range
        self.stop: used if the coordinate exists over a range, as in exons/features.  otherwise, points to self.start
        self._pos: returned by @property self.pos, contains either just self.start or (self.start, self.stop)

    methods:
        __init__(chromosome, start, stop = False:
            initialize with chromosome position - casts chromosome to integer representation if not already
            initialize start/stop positions - cast to int
            default is single position
        __repr__:
            internal representation of Coordinate
            return repr((self.chromosome, self.start))
        @property
        pos:
            return self._pos, content of self._pos dependent on initialization (see above)

    '''
    def __init__(self, chromosome, start, stop = False):
        if not chromosome:
            raise AttributeError("No chromosome supplied for coordinate")
        if not start:
            raise AttributeError("No start position supplied for coordinate")
        self.chromosome = int(chromosome)

        self.start = int(start)
        if not stop:
            self.stop = start
            self._pos = self.start
        else:
            self.stop = int(stop)
            self._pos = (self.start, self.stop)

    def __len__(self):
        if self.stop != self.start:
            return self.start-self.stop
        else:
            return 1

    def __repr__(self):
        return repr((self.chromosome, self.pos))

    @property
    def pos(self):
        return self._pos

class Alignment(object):
    pass

class Trimmed(object):
    '''may be unnecessary'''
    pass

class PolyDiv(object):
    '''accepts massprf_preprocess output, accepts scaling, returns a reinstantiation descendant that cannot scale'''
    pass

'''adaptor builders'''

class CDSAdaptorBuilder(object):
    '''CDSAdaptorBuilderis a factory for coding sequence adaptors classes
        To extend functionality to new file types, write a new class of Adaptor that creates objects of class Gene and accepts **kwargs
        return the new class using an elif statement in __new__ of CDSAdaptorFactory'''
    def __new__(self,cls, **kwargs):
        if cls == 'gff3':
            return GffUtilAdaptor(**kwargs)

class VariantsAdaptorBuilder(object):
    '''Interpret passed variant format file and return properly interfaced variant file'''
    def __new__(self,frmt, variants):
        if frmt == 'vcf':
            return PyVCFAdaptor(variants)

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
        if frmt == 'placeholder':
            #add new formats for parsing here
            pass
        elif frmt == "FASTA":
            #default to FASTA format for genome
            return FastaGenomeAdaptor(reference, species, variants)
        else:
            pass

class HomologAdaptorBuilder(object):
    '''Adapts input homolog formats to HomologAdaptor'''
    def __new__(self, homologyfile, frmt = "CSVa"):
        if not Path(homologyfile).is_file():
            raise AttributeError("No valid homology file specified")
        if frmt == "CSVa":
            #this file type is useful when you have two distinct divergent species with a CSV file mapping of the homologs
            return CSVaHomologs(homologyfile)
        elif frmt == "CSVb":
            return CSVbHomologs(homologyfile)

'''adaptors'''

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

    def createCodingSequence(self, feature, reference, **kwargs):
        if not isinstance(feature, gffutils.feature.Feature):
            raise AttributeError("passed feature is not of type gffutils.feature.Feature")
        coordinates = list(map(lambda x: Coordinate(x, x['chrom'], x['start']-1, x['stop']), self.db.children(feature, featuretype = 'CDS', order_by = 'start')))
        name = feature['name']
        gene_id = feature['ID']
        strand = feature['strand']
        cds = CodingSequence(reference, name, coordinates, strand, gene_id = gene_id)
        reference.addCDS(cds)
        return cds


## PLEASE NOTE****** GFF3 is 1-indexed, inclusive on both sides.  Therefore, to translate to python string slicing, the lower coordinate must be -1.  
## if the range on the GFF is 1-10, then the python string indices will be 0-9.  Therefore, str[0:10].  Check this to make sure it is right.
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
                    coordinate = Coordinate(record.CHROM, record.POS-1)
                    for strain in variants.strains:
                        calls.append(allele(record.genotype(strain)['GT']))

                    ref = record.REF
                    alt = record.ALT[0] #See above - SNP reading breaks >1 alternate allele

                    alleles = [alt if gt else ref for gt in calls]

                    for a, s in zip(alleles, variants.strains):
                        variants.addVariant(Variant(coordinate, s, a, ref))
            self.file.close()
        return variants

class FastaGenomeAdaptor(object):
    '''FastaGenomeAdaptor adapts biopython SeqIO genome into the genome class
    "species" here is a generic term, eg rice, banana, etc, that is useful for bluntly differentiating objects of type genome

    Consider reworking the external representation of reference genome and its variants
    '''
    def __new__(self, reference, species, variants = None):
        #GenomeAdaptorBuilder already did the type checking of reference to insure it is a filetype
        reference_filename = reference
        reference_file = open(reference)

        parsedchromosomes = {int(n.id):str(n.seq) for n in SeqIO.parse(reference_file, "fasta")}
        referenceGenome = ReferenceGenome(species, parsedchromosomes)
        reference_file.close()

        if variants:
            for strain in variants.strains:
                referenceGenome.addStrain(VariantGenome(referenceGenome, strain, variants.of_strain(strain)))

        return referenceGenome
        #if variants are passed, we are going to build a bunch of new genome objects and return the genomes as a list
        #the genomes should build & insert their variants themselves. 

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

''' begin algorithm objects '''
'''these should all inherit the ability to spawn subprocesses and maybe do file writing'''

class SubProcessSpawner(object):
    '''parent of Aligner, MASSPRF_Pre'''
    def __init__(self, cores):
        self.numcores = cores-1 # cores - 1 because the allocating process consumes 1 core

class Aligner(SubProcessSpawner):
    '''Aligner: Requires MUSCLE, Biopython

    given a collection of coding sequences, pass them to MUSCLE and align them.  return them in a form analogous to their native structures
    '''
    pass

class Trimmer(object):
    '''Trimmer:
    Given an alignment, trim all '-' and return sequences in their native structure
    '''
    pass

class MASSPRF_Pre(SubProcessSpawner):
    '''MASSPRF_Pre:
    Given aligned & trimmed sequences, spawn a MASSPRF subprocess to preprocess them w/ annotation of polymorphism sites
    '''
    pass

class Scaler(object):
    '''Scaler:
    Given a MASSPRF_preprocess output file, scale the preprocessed genes to a computationally-friendly length
    '''
    pass

class MASSPRF_Queuer(object):
    '''MASSPRF_Queuer:
    Given a list of preprocessed & scaled MASSPRF output files, export a text file with corresponding MASSPRF commands for queueing system
    '''
    pass
    
class DirectoryTree(object):
    """docstring for DirectoryTree"""
    def __init__(self, rootdir):
        self._rootdir = Path(rootdir)
        if not self.rootdir.is_dir():
            raise IOError("Invalid root directory supplied")
        self.outdir = self.mksubdir(self.rootdir,"out")
        self.genomedir = self.mksubdir(self.outdir,"genomes")
        self.alignments = self.mksubdir(self.outdir,"alignments")
        self.prf_pre = self.mksubdir(self.outdir,"prf_pre")
        self.scaled = self.mksubdir(self.outdir,"scaled")
        self.csv = self.mksubdir(self.outdir, "csv")

    def mksubdir(self, root, directory):
        directory = root.joinpath(str(directory))
        if not directory.is_dir():
            directory.mkdir()
        return directory

class Program(object):
    """docstring for program"""
    def __init__(self):
        pass
        
