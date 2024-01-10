#! /usr/bin/env python

# This script is designed to remove any sites at which not all samples are heterozygous
# and also to remove sites where any sample contains a non-nucleotide (or Null) genotype call

# Usage
# python ./src/CleanData_KeepHetSitesRemoveNullCalls.py VCFtableName

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
OutFile=open(SNPfile+"_HetsOnlyNoNull",'w')
OutFile.write(SNPhead[0])
for i in range(1,len(SNPhead)):
	OutFile.write('\t'+SNPhead[i])

# Processing each line in the VCF file
for Locus in SNPs:
	Homs=0 # Count the number of samples that are homozygous
	Nucleotide='Yes' # Flag for Null Genotype Calls.
	# Loop across samples in the file
	for i in range(5,len(Locus)):
		# Determine whether the sample is homozygous.
		# Genotypes are recorded like "A/T" in these files, and we are dealing with them as strings
		if Locus[i][0]==Locus[i][2]:
			# If the sample is homozygous add one to the homozygous sample counter
			Homs+=1
		# Check for non-nucleotide bases, and change the flag to "No" if we find any
		if Locus[i][0] not in ['A', 'T', 'G', 'C'] or Locus[i][2] not in ['A', 'T', 'G', 'C']:
				Nucleotide='No'
	# Lines for which no samples are homozygous and no samples contain null (or non nucleotide) genotype calls get written to the output file
	# We also remove mitochrondrial sites since we are interested in the nuclear genome.
	if Homs == 0 and Nucleotide == 'Yes' and "itochondrion" not in Locus[0]:
		# Write the line to the output file
		OutFile.write('\n'+Locus[0])
		for i in range(1,len(Locus)):
			OutFile.write('\t'+Locus[i])

# Close the file
OutFile.close()
print('Done')
