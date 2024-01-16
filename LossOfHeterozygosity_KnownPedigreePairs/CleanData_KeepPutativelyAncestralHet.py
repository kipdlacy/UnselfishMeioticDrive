#! /usr/bin/env python

#### Note that this must be run after all sites that are heterozygous in haploids are removed from the dataset.

# 20210227
# Writing script to read in VCF-style table and write only sites that were putatively ancestrally heterozygous in a line to an output file.
# Putative ancestral heterozygosity at any single locus is defined as follows:
# 	Heterozygous in at least one diploid
# 	All heterozygous diploid have the same genotype
# 	Any homozygous diploids bear one of the alles present in the heterozygous genotype
# 	Haploids should bear one of the alleles present in the heterozygous genotype

# Usage
# python ./src/CleanData_KeepPutativelyAncestralHet.py VCFtableName

# Import Modules
import sys

# Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

# Retrieve the VCF-style table from command line argument
SNPfile=sys.argv[1]

# Read the header and the data separately
SNPhead=read_csv(SNPfile)[0]
print('SNPhead')
print(SNPhead)
SNPs=read_csv(SNPfile)[1:]
print('SNPs[0]')
print(SNPs[0])

# Create output files for storing results

# File for Variants that meet the criteria of being putatively ancestrally heterozygous
OutFile=open(SNPfile+"_PutAncHet",'w')
OutFile.write(SNPhead[0])
for i in range(1,len(SNPhead)):
	OutFile.write('\t'+SNPhead[i])

# File for variants that do NOT meet those criteria
NotFile=open(SNPfile+"_NOTPutAncHet",'w')
NotFile.write(SNPhead[0])
for i in range(1,len(SNPhead)):
	NotFile.write('\t'+SNPhead[i])

# Assess each variant "Locus" in the dataset
for Locus in SNPs:
	# Checking to make sure that the line isn't an empty line at the end of the file
	if len(Locus)==len(SNPs[0]):
		Hets=[] # List for heterozygous genotypes
		Homs=[] # List for homozygous genotypes
		AllHetsSameGenotype='False' # Flag for whether it is true that for that locus all heterozygous genotypes are identical
		WriteOut='No' # Flag denoting to which file the Variant should be written
		# For each sample
		for k in range(5,len(Locus)):
			# Check whether the genotype is heterozygous (this works because genotypes are written as 'A/T' in these files)
			if Locus[k][0] != Locus[k][2]:
				# If heterozygous append to the hets list
				Hets.append(Locus[k])
			else:
				# Otherwise append to the homs list
				Homs.append(Locus[k])
		# If there aren't zero heterozygous sites
		if len(Hets) != 0:
			# if all of the heterozygous genotypes are the same
			if len(set(Hets))==1:
				AllHetsSameGenotype='True'
		# Assuming all heterozygous genotypes are the same, it is possible that the variant is putatively ancestrally heterozygous
		# but that depends on whether there are any homozygous genotypes, and if so, what those homozygous genotypes are.
		if AllHetsSameGenotype=='True':
			# If there are no homozygous samples then the variant automatically passes
			if len(Homs) == 0:
				WriteOut='Yes'
			# If there is at least one homozygous sample then we need to inspect the genotype
			elif len(Homs) >= 1:
				HomPassCount=0 # <- instantiating a counter for homozygous genotype inspection
				for Geno in Homs:
					# if the allele found in the homozygous genotype is the same as one of the alleles in the heterozygous genotype,
					# then this homozygous genotype is compatible with the variant being putatively ancestrally heterozygous
					if Geno[0]==Hets[0][0] or Geno[0]==Hets[0][2]:
						HomPassCount+=1
				# If all homozygous genotypes are compatible with the variant being putatively ancestrally heterozygous...
				if len(Homs)==HomPassCount:
					WriteOut='Yes'
		# Write the variant to the appropriate output file
		if WriteOut=='Yes':
			OutFile.write('\n'+Locus[0])
			for i in range(1,len(Locus)):
				OutFile.write('\t'+Locus[i])
		if WriteOut!='Yes':
			NotFile.write('\n'+Locus[0])
			for i in range(1,len(Locus)):
				NotFile.write('\t'+Locus[i])

# Close the output files
NotFile.close()
OutFile.close()
print('Done')
