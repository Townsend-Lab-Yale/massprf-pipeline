import numpy as np
import os
from argparse import ArgumentParser
#build a text file of polymorph:divergence pairs to align and trim in alignTrim.py
def is_valid_file(parser, arg):
    if not os.path.exists(arg):
    	if parser:
        	parser.error("The file/folder %s does not exist!" % arg)
    	else: 
    		return False
    else:
    	return os.path.abspath(arg)

if __name__== "__main__":
	parser = ArgumentParser(description='build a text file of pairs for alignTrim')
	parser.add_argument('-s', dest = "directory", type=lambda x: is_valid_file(parser,x), nargs='+', help = 'sourcefolder', required=True)
	args = parser.parse_args()
	
	pmfiles = os.listdir(args.directory[0])
	divfiles = os.listdir(args.directory[1])
	polymorphs = list(filter(lambda x: x.endswith('polymorphs.txt'), pmfiles))
	orthologs = list(filter(lambda x: x.startswith('NtLA_'), divfiles))
	output = []
	for x in polymorphs:
		output.append((x, list(filter(lambda c: x[:-15] in c, orthologs))))
	ofile = open("pairs.txt", "w")
	error = open("aligntrimpreerror.txt",'w')
	for x in output:
		if not ''.join(x[1]) == '':
			ofile.write("source ~/.bashrc; cd /ysm-gpfs/home/jss245/scratch/20161216aligntrim; python alignTrim.py -s " + args.directory[0]+'/' +x[0] + ' '+ args.directory[1]+'/' + ''.join(x[1]) +';'+'\n')
		else:
			error.write(x[0])
	error.close()
	ofile.close()

