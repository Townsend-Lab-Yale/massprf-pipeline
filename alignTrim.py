import os
import itertools
import subprocess
import sys
import signal
import psutil
from argparse import ArgumentParser
from Bio.Align.Applications import MuscleCommandline
import numpy as np
from Bio import SeqIO,AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def is_valid_file(parser, arg):
	# used to check in argparse for valid file
    if not os.path.exists(arg):
    	if parser:
        	parser.error("The file/folder %s does not exist!" % arg)
    	else: 
    		return False
    else:
    	return os.path.abspath(arg)

def kill(proc_pid):
	#added as a workaround to kill MUSCLE process
    process = psutil.Process(proc_pid)
    for proc in process.children(recursive=True):
        proc.kill()
    process.kill()

class MUSCLE (object):
	def __init__(self):
		#opens a stdin/stdout pipe to MUSCLE
		#could be a bottleneck
		muscle_cline = MuscleCommandline(clwstrict=True)
		self.child = subprocess.Popen(str(muscle_cline), stdin = subprocess.PIPE, stdout = subprocess.PIPE, stderr = subprocess.PIPE, universal_newlines = True, shell=True, preexec_fn=os.setsid)
	
	def align(self, seqs):
		# uses MUSCLE pipe to send  sequences for alignment
		#could be a bottleneck
		SeqIO.write(seqs, self.child.stdin, "fasta")
		self.child.stdin.close()
		alignment = AlignIO.read(self.child.stdout,"clustal")

		#added as a work around to see if killing the muscle subprocess would result in functional parallel computing through slurm
		kill(self.child.pid)
		
		return alignment
	
	def trim(self, alignment):
		inserts = []
		# this could probably be converted to np.arrays and use logical indexing to delete gaps for (?) more efficiency
		for i in range(len(alignment)): # go through each sequence in the alignment
			inserts.append([c for c, x in enumerate(list(alignment[i].seq)) if x =='-']) # check the sequence for gaps ('-'), append a list of the indices of gaps for each alignment
		inserts=set([item for sublist in inserts for item in sublist]) # iterate through the nested list and condense into a set (removes duplicates)
		inserts = list(inserts) #flip it back to a list for sorting and indexing
		for i in range(len(alignment)): # go back through each alignment
			record = list(alignment[i].seq) # set to current alignment
			for z in sorted(inserts,reverse=True): # index of gaps, in reverse so that the indices do not change
				del record[z]

			alignment[i].seq = ''.join(record)
		return alignment

def script(args):
	#spawn alignment/trim class
	muscle = MUSCLE()

	#default output folder
	if not args.outfolder:
		outfolder = './ATout/'
		if not is_valid_file(None, outfolder):
			os.makedirs(outfolder)
		outfolder = os.path.abspath(outfolder)
	else:
		outfolder = args.outfolder[0]
	outfolder += '/'

	#progress file
	prog = open('progress.txt','a')
	prog.write(args.inputseqs[0] + '\n')
	prog.close()

	#seq i/o into dictionary, generator
	seqParser =  {os.path.basename(f)[:-4]:list(SeqIO.parse(f,"fasta")) for f in args.inputseqs}
	seqs = (r for r in itertools.chain.from_iterable(seqParser.values()))
	#align them, trim them, put into dictionary that facilitates output
	aligned = muscle.align(seqs)
	slicedaligned = muscle.trim(aligned)
	slicedaligned = {s.description:s.seq for s in slicedaligned}

	for f in seqParser.keys(): 
		#match aligned/trimmed sequences back to the original file they came from, write to a new file of similar name
		tempseqs = []
		for record in seqParser[f]:
			tempseqs.append(SeqRecord(Seq(slicedaligned[record.description]), id=record.description, description = ''))
		seqlen = len(tempseqs[0])
		SeqIO.write(tempseqs,outfolder+f+'_'+str(seqlen)+'.txt',"fasta")

if __name__== "__main__":
	parser = ArgumentParser(description='align two sequence sets, output trimmed alignment')
	parser.add_argument('-s', dest = "inputseqs", type=lambda x: is_valid_file(parser,x), nargs='+', help = 'sources', required=True)
	parser.add_argument('-o', dest="outfolder", type=lambda x: is_valid_file(parser,x), nargs = 1, help = 'output folder')
	args = parser.parse_args()
	script(args)