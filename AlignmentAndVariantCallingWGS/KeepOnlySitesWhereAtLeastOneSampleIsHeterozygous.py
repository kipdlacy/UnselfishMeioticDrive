#! /usr/bin/env python

# This script is designed to remove any sites at which no samples are heterozygous

# Usage
# python ./src/KeepSitesWhereAtLeastOneSampleIsHeterozygous.py VCFtableName

# Import Modules
import sys

# Defining function to read tabular data
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

# Passing the vcf table name as an arguments from the command line
SNPfile=sys.argv[1]

# Reading the header and the data from the VCF file
SNPhead=read_csv(SNPfile)[0]
SNPs=read_csv(SNPfile)[1:]

# Creating and writing the header of the output file
OutFile=open(SNPfile+"_AtLeastOneSampleHet",'w')
OutFile.write(SNPhead[0])
for i in range(1,len(SNPhead)):
	OutFile.write('\t'+SNPhead[i])

# Processing each line in the VCF file
for Locus in SNPs:
	Homs=0 # Count the number of samples that are homozygous
	SampleCounter=0
	# Loop across samples in the file
	for i in range(5,len(Locus)):
		# Determine whether the sample is homozygous.
		# Genotypes are recorded like "A/T" in these files, and we are dealing with them as strings
		SampleCounter+=1
		if Locus[i][0]==Locus[i][2]:
			# If the sample is homozygous add one to the homozygous sample counter
			Homs+=1
	# Lines for which not all samples are homozygous get written to the output file
	if Homs < SampleCounter:
		# Write the line to the output file
		OutFile.write('\n'+Locus[0])
		for i in range(1,len(Locus)):
			OutFile.write('\t'+Locus[i])

# Close the file
OutFile.close()
print('Done')
