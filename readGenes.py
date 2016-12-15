import argparse
if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Separate gene info from txt file list of genes')
	parser.add_argument('-s',type=str, nargs=1, help = 'source to separate', required=True)
	args = parser.parse_args()
	
	srcFile = open(args.s[0], 'rU')
	srcFile = list(srcFile)
	numLines = len(srcFile)

	i = 0
	for n in range(numLines):
		if srcFile[n] == '\n':
			geneName = srcFile[i].strip('\n')
			geneLst = srcFile[i:n]
			geneStr = ''.join(geneLst)

			output = open(geneName + '.txt','w')
			output.write(geneStr)
			output.close()
			
			i = n+1
