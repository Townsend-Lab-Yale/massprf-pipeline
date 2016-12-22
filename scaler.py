from argparse import ArgumentParser
import os

def is_valid_file(parser, arg, permission):
	#http://stackoverflow.com/questions/11540854/file-as-command-line-argument-for-argparse-error-message-if-argument-is-not-va
    if not os.path.exists(arg) and parser:
        parser.error("The file %s does not exist!" % arg)
    elif not os.path.exists(arg) and not parser:
    	raise IOError("The file %s does not exist!" % arg)
    else:    
		return open(arg, permission)  # return an open file handle
def write_error(file, string):
	errstring = ''
	try:
		f = is_valid_file(False, file, 'a')
		f.write(string)
	except:
		errstring = "Could not open", file , '\n', 'Attempted to write ', string
		f = open("logerror.txt",'a')
		f.write(errstring)
		f.close()

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

		self._genename = ''
		self._ogpol = ''
		self._ogdiv = ''
		self._scalex = 1

		self.genename = f.name.split('/')[-1].split('_')[0]	
		self.ogpol = next(filter(lambda s: s.startswith('Polymorphism:'),lines)).split()[1]
		self.ogdiv = next(filter(lambda s: s.startswith('Divergence:'),lines)).split()[1]
		f.close()
		
		self.scaledpol = self.ogpol
		self.scaleddiv = self.ogdiv

	def __len__(self):
		return len(self.scaledpol)

	@property
	def genename(self):
		return str(self._genename)
	@genename.setter
	def genename(self, name):
		self._genename = name

	@property
	def ogpol(self):
		return self._ogpol
	@ogpol.setter
	def ogpol(self, pol):
		self._ogpol = pol
	
	@property
	def ogdiv(self):
		return self._ogdiv
	@ogdiv.setter
	def ogdiv(self, div):
		self._ogdiv = div
		
	@property
	def scalex(self):
		return self._scalex
	@scalex.setter
	def scalex(self, scale):
		self._scalex = scale
		
	@property
	def scaledpol(self):
		return self._scaledpol
	@scaledpol.setter
	def scaledpol(self,sequence):
		self._scaledpol = self._scale(len(self), sequence)

	@property
	def scaleddiv(self):
		return self._scaleddiv
	@scaleddiv.setter
	def scaleddiv(self, sequence):
		self._scaleddiv = self._scale(len(self),sequence)

	def _scale(self,l,sequence):
		#logic for scaling
		if l <= 600:
			self.scalex = 1
		elif l > 600 and l <=1800:
			self.scalex = 3
		elif l > 1800 and l <= 3600:
			self.scalex = 6
		elif l > 3600 and l <= 5400:
			self.scalex = 9
		elif l > 5400 and l <= 7200:
			self.scalex = 12
		elif l > 7200 and l <= 9000:
			self.scalex = 15
		elif l > 9000 and l <= 10800:
			self.scalex = 18
		elif l > 10800 and l <= 12600:
			self.scalex = 21
		elif l > 12600 and l <= 14400:
			self.scalex = 24
		elif l > 14400 and l <= 16200:
			self.scalex = 27
		elif l > 16200 and l <= 18000:
			self.scalex = 30
		elif l > 18000 and l <= 30000:
			self.scalex = 50
		elif l > 30000 and l <= 70000:
			self.scalex = 117
		else:
			self.scalex = False
			write_error("output/scalingerrors.txt", self.genename + ' is of length ' + l + '\n  Did not scale' + self.genename)

		if self.scalex:
			#main body of scaling
if __name__ == "__main__":
	parser = ArgumentParser(description='extract and scale consensus sequences from MASSPRF_preprocess')
	parser.add_argument('-i', dest="filename", required=True, help="input consensus file from MASSPRF preprocess", metavar="FILE", type=lambda x: is_valid_file(parser,x ,'r'))
	args = parser.parse_args()

	consensus = CONSENSUS(args.filename)