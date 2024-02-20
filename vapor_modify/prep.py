from __future__ import print_function
vapor_version='vapor.py V1.0 12-21-2015'
def print_read_me():
	print (vapor_version)
	print ('Contact: Xuefang Zhao (xuefzhao@umich.edu)')
	print ('')
	print ('Usage: vapor.py [Options] [Parameters]')
	print ('Options: ')
	print ('	svelter')
	print ('	vcf')
	print ('	bed')
	print ('Parameters:')
	print ('	--sv-input:		input file in bed or vcf format')
	print ('	--output-path:		folder where the recurrence plot will be kept')
	print ('	--reference:		reference genome that pacbio files are aligned against')
	print ('	--pacbio-input:		absolute path of input pacbio file')

def readme_melt():
	print (vapor_version)
	print ('Contact: Xuefang Zhao (xuefzhao@umich.edu)')
	print ('')
	print ('Usage: vapor.py ins [Parameters]')
	print ('Parameters:')
	print (' --sv-input-prefix:	prefix of input file in vcf and fa format	')
	print ('	--output-path:		folder where the recurrence plot will be kept')
	print ('	--reference:		reference genome that pacbio files are aligned against')
	print ('	--pacbio-input:		absolute path of input pacbio file')

def readme_bed():
	print (vapor_version)
	print ('Contact: Xuefang Zhao (xuefzhao@umich.edu)')
	print ('')
	print ('Usage: vapor.py bed [Parameters]')
	print ('Parameters:')
	print ('	--sv-input:		input file in bed format with SV type labeled in the last column')
	print ('	--output-path:		folder where the recurrence plot will be kept')
	print ('	--output-file:		name of output file including vapor.py scores')
	print ('	--reference:		reference genome that pacbio files are aligned against')
	print ('	--pacbio-input:		absolute path of input pacbio file')


def readme_vcf():
	print (vapor_version)
	print ('Contact: Xuefang Zhao (xuefzhao@umich.edu)')
	print ('')
	print ('Usage: vapor.py vcf [Parameters]')
	print ('Parameters:')
	print ('	--sv-input:		input file in vcf format')
	print ('	--output-path:		folder where the recurrence plot will be kept')
	print ('	--reference:		reference genome that pacbio files are aligned against'		)
	print ('	--pacbio-input:		absolute path of input pacbio file')

