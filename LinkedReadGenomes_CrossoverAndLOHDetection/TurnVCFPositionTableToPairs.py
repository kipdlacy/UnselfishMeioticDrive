#! /usr/bin/env python

### 20220506
### This script is designed to read in a positions-only table output from a vcf using gatk variantstotable
### and output all pairs of SNPs.

### Import Modules
import sys
import numpy

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Passing Arguments
InFile=sys.argv[1]
List=read_csv(InFile)
if len(sys.argv)>2:
	print('List[0]:')
	print(List[0])
	print('List[1]:')
	print(List[1])
	print('len(List)')
	print(len(List))

### Opening the Output File
OutFile=open(InFile+'_Pairs','w')
OutFile.write('Chrom\tFirstSNP\tSecondSNP')

### Looping through the input file and writing pairs to the output file
### Skipping the header (index=0) and starting with the second SNP (index=2) since we are grabbing pairs
for i in range(2,len(List)):
	### Checking to make sure the line isn't an empty line
	if len(List[i])==len(List[1]):
		### If the current variant and the previous one are on the same chromosome
		if List[i][0]==List[i-1][0]:
			### Write out a line that contains the Chromosome, the previous variant, and the current variant
			OutFile.write('\n'+List[i][0]+'\t'+List[i-1][1]+'\t'+List[i][1])
		### if not
		else:
			### Check to see if the user has asked for the script to be run in a verbose manner for debugging
			### this is achieved by adding extra arguments 
			if len(sys.argv)>2:
				print('Starting '+List[i][0])

### Closing the output file and ending the script
OutFile.close()
if len(sys.argv)>2:
	print('Done :)')
