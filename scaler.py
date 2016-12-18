from argparse import ArgumentParser
import os

def is_valid_file(parser, arg):
	#http://stackoverflow.com/questions/11540854/file-as-command-line-argument-for-argparse-error-message-if-argument-is-not-va
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return open(arg, 'r')  # return an open file handle

class CONSENSUS(object):
	def __init__(self, f):
		lines = [line for line in f.read().splitlines() if line]
		if "Can't" in lines[0]:
			raise ValueError("file error due to silent clustering")
		if "Error" in lines[0]:
			raise ValueError("error in MASSPRF_preprocess input parameters")
		if "MAC-PRF" not in lines[0]:
			raise ValueError("Not a massprf file")
		if "Mission accomplished" not in lines[-1]:
			raise ValueError("file did not pass preprocessing")

		self.genename = f.name.split('/')[-1].split('_')[0]
		self.ogpolymorphism = next(filter(lambda s: s.startswith('Polymorphism:'),lines)).split()[1]
		self.ogdivergence = next(filter(lambda s: s.startswith('Divergence:'),lines)).split()[1]
		print(self.ogdivergence)
		f.close()
		self.ogplen = len(self.ogpolymorphism)
		self.ogdlen = len(self.ogdivergence)

		if self.ogplen != self.ogdlen:
			raise ValueError("Somehow, divergence and polymorphism sequence lengths are different")


if __name__ == "__main__":
	parser = ArgumentParser(description='extract and scale consensus sequences from MASSPRF_preprocess')
	parser.add_argument('-i', dest="filename", required=True, help="input consensus file from MASSPRF preprocess", metavar="FILE", type=lambda x: is_valid_file(parser,x))
	args = parser.parse_args()

	consensus = CONSENSUS(args.filename)