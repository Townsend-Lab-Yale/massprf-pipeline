from argparse import ArgumentParser
import csv
import glob
import os
def is_valid_file(parser, arg):
	# used to check in argparse for valid file
    if not os.path.exists(arg):
    	if parser:
        	parser.error("The file/folder %s does not exist!" % arg)
    	else: 
    		return False
    else:
    	return os.path.abspath(arg)

if __name__=="__main__":
	parser = ArgumentParser()
	parser = ArgumentParser(description='crawl folder, make CSV of consensus files w/ gene lengths and massprf command')
	parser.add_argument('-i', dest = "infolder", type=lambda x: is_valid_file(parser,x), nargs='+', help = 'sources', required=True)
	parser.add_argument('-o', dest="ofolder", type=lambda x: is_valid_file(parser,x), nargs = 1, help = 'ofolder')
	args = parser.parse_args()

	files = glob.glob(args.infolder[0]+'/*_prep.txt')
	fwlen = []

	outcsv = args.ofolder[0]+'/jobs.csv'
	jobsfile = args.ofolder[0] + '/jobs.list'
	with open(outcsv,'a') as outputfile:
		for f in files:
			file = open(f, 'r')
			lines = [line for line in file.read().splitlines() if line]
			file.close()

			genename = file.name.split('/')[-1].split('_')[0]
			try:
				poly = next(filter(lambda s: s.startswith('Polymorphism:'),lines)).split()[1]
				div = next(filter(lambda s: s.startswith('Divergence:'),lines)).split()[1]
				clength = len(poly)
				polyfile = args.ofolder[0] + '/' + genename + '_poly_consensus.txt'
				divfile = args.ofolder[0] + '/' + genename + '_div_consensus.txt'
				fwlen.append((clength, genename, polyfile, divfile))
				with open(polyfile,'w') as file:
					file.write('>Polymorphism \n' + poly)
				with open(divfile,'w') as file:
					file.write('>Divergence \n' + div)
			except:
				print(genename)
		fwlen.sort()
		string = '';
		for t in fwlen:
			string = t[1] +', ' + str(t[0]) + ', ' +  t[2] + ', ' + t[3] + ", source ~/.bashrc; source ~/.bash_profile; cd ../massprfout; massprf -p  " + t[2] + ' -d ' + t[3] + ' -ic 1 -sn 49 -ci_m 1 -r 1 -ci_r 1 -exact 0 -mn 10000 -t 1 -v 0 >> ' + t[1] + '.prfout \n'
			outputfile.write(string)
	with open(jobsfile, 'a') as jobsfile:
		parsefile = open(outcsv, 'r')
		print(parsefile.name)
		lines = [line for line in parsefile.read().splitlines() if line]
		for line in lines:
			splitline = line.split(',')
			if int(splitline[1]) <= 600:
				jobsfile.write(splitline[4]+'\n')
			
