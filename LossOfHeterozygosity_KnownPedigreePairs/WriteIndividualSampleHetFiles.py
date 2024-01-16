#! /usr/bin/env python

# 20210713
# Script to write one outfile per sample in a VCF-style table and only sites where that sample is heterozygous

# Usage
# python ./WriteIndividualSampleHetFiles.py FileName.table

# Import Modules
import sys

# Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

# Read the header and body of the table from the input file
Head=read_csv(sys.argv[1])[0]
print('Head:')
print(Head)
Body=read_csv(sys.argv[1])[1:]
print('Body[0]:')
print(Body[0])

# Process each sample in the table
for i in range(5,len(Head)):
	# Create and open an output file for each sample with a suffix derived from the sample name
	HetsFile=open(sys.argv[1]+'_'+Head[i].split('.GT')[0]+'Het','w')
	# Write the header to the output file
	HetsFile.write(Head[0])
	for k in range(1,5):
		HetsFile.write('\t'+Head[k])
	# Write the sample name to the output file header
	HetsFile.write('\t'+Head[i])
	# Then loop through each variant in the table
	for Line in Body:
		# Check to make sure the line isn't empty (like there sometimes are at the end of files)
		if len(Line)==len(Body[0]):
			# Check whether the individual is heterozygous at that locus
			if Line[i][0]!=Line[i][2]:
				# If so, write the Variant to the output file
				# first the positional information and metadata
				HetsFile.write('\n'+Line[0])
				for k in range(1,5):
					HetsFile.write('\t'+Line[k])
				# then the genotype for the relevant sample
				HetsFile.write('\t'+Line[i])

