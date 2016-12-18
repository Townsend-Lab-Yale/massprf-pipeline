#!/Users/jaystanley/anaconda/bin/python
from Bio import SeqIO
import csv
import time
import argparse

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='from a csv list of orthologs -csv, map first argument (-gf) to a gff3, map second argument (-fas) to CDS in fasta')
	parser.add_argument('-csv',type=str, nargs=1, help = 'source csv', required=True)
	parser.add_argument('-gf',type=str, nargs=1, help = 'polymorphism data', required=True)
	parser.add_argument('-fas',type=str, nargs=1, help = 'genome', required=True)
	args = parser.parse_args()

	csvpath = args.csv[0]
	refgffpath = args.gf[0]
	orthofaspath = args.fas[0]
	
	#open the reference gff map of refgff
	f = open(refgffpath,'r')
	gff = list(f)
	f.close()
	# get cds of orthologs
	orthocds = list(SeqIO.parse(orthofaspath,"fasta"))

	t0 = time.time()
	with open(csvpath,'r') as csvfile:
		orthologs = csv.reader(csvfile)
		#ignore all w/o orthologs in both species
		orthologs = list(filter(lambda r: r[0] and r[7], orthologs))

		for n in range(len(orthologs)): # using n here in order to facilitate progress tracking
			percent = n/len(orthologs) * 100
			if percent % 1 == 0:
				time_elapsed = time.time()-t0
				print(str(time_elapsed) +'\n'+ str( n/len(orthologs) * 100 ) + '%')
			row = orthologs[n]
			ncras = row[0]
			ntet = row[7]

			fnc = open('/Users/jaystanley/Documents/School/Grad/townsend/ntetncras/output/'+ncras+'_polymap.txt','w')
			buf = ''
			for line in gff:
				if line.split()[9].replace('"','').replace(';','').startswith(ncras):
					buf += line
			if buf != '':
				buf = ncras + '\n' + buf
			else: 
				buf = ncras + '\nFAILED: no gene found in gff'
			fnc.write(buf)
			fnc.close()

			buf = ''
			fnt = open('/Users/jaystanley/Documents/School/Grad/townsend/ntetncras/output/'+ntet+'_'+ncras+'.txt','w')
			gn = list(filter(lambda s: s.name.startswith('jgi|Neute_matA2|'+ntet[5:]), orthocds))
			
			if not gn:
				buf += 'FAILED: no CDS found in fasta'
				buf += '>' + ntet + '\n'
			elif len(gn)>1:
				buf += 'WARNING: multiple CDS identified \n'
				for x in range(len(gn)):
					buf += '>' + ntet + '_'+ncras+'_'+str(x)+'\n'
					buf += str(gn[x].seq) + '\n'
			elif len(gn)==1:
				buf += '>' + ntet + '\n'
				buf += str(gn[0].seq)
			else:
				buf += 'failed, unidentified case'

			fnt.write(buf)
			fnt.close()
	print(time.time()-t0)