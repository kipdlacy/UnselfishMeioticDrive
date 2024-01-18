#! /usr/bin/env python

# 20210227
# Writing script to read in VCF-style table and write only sites that were putatively ancestrally heterozygous in a line to an output file.
# Putative ancestral heterozygosity at any single locus is defined as follows:
# 	Heterozygous in at least one diploid
# 	All heterozygous diploid have the same genotype
# 	Any homozygous diploids bear one of the alles present in the heterozygous genotype
# 	All haploids must be hemizygous for either of the alleles present in the heterozygote(s)

# Usage
# python ./src/CleanData_KeepPutativelyAncestralHet.py VCFtableName

# Import Modules
import sys

# Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

# Passing Arguments
SNPfile=sys.argv[1]
SNPhead=read_csv(SNPfile)[0]
print('SNPhead')
print(SNPhead)
SNPs=read_csv(SNPfile)[1:]
print('SNPs[0]')
print(SNPs[0])

# OutFiles
OutFile=open(SNPfile+"_PutAncHet",'w')
OutFile.write(SNPhead[0])
for i in range(1,len(SNPhead)):
	OutFile.write('\t'+SNPhead[i])
NotFile=open(SNPfile+"_NOTPutAncHet",'w')
NotFile.write(SNPhead[0])
for i in range(1,len(SNPhead)):
	NotFile.write('\t'+SNPhead[i])

# Instantiate lists of diploids and haploids
DipInd=[]
HapInd=[]
# Distinguish diploids from haploids based on internal codes
for i in range(5,len(SNPhead)):
	if 'SM' in SNPhead[i] or 'KL' in SNPhead[i]:
		HapInd.append(i)
		print('Hap: '+SNPhead[i])
	else:
		DipInd.append(i)
		print('Dip: '+SNPhead[i])

# Loop through each SNP in the dataset
for Locus in SNPs:
	# make sure the line isn't empty
	if len(Locus)==len(SNPs[0]):
		Hets=[] # List for heterozygous genotypes
		Homs=[] # List for homozygous genotypes
		AllHetsSameGenotype='False' # Flag for whether it is true that for that locus all heterozygous genotypes are identical
		WriteOut='No' # Flag denoting to which file the Variant should be written
		# For each diploid sample
		for k in DipInd:
			# check if heterozygous and if so add to hets list
			if Locus[k][0] != Locus[k][2]:
				Hets.append(Locus[k])
			else:
				# otherwise add to homozygotes list
				Homs.append(Locus[k])
		# If not all samples are homozygous
		if len(Hets) != 0:
			# if there is only one heterozygous sample then all samples of course have the same genotype
			if len(Hets) == 1:
				AllHetsSameGenotype="True"
			elif len(Hets) > 1:
				# But if there's more you need to check if there is only one
				if len(set(Hets))==1:
					AllHetsSameGenotype="True"
		# If all hets have the same genotype then it's plausible that the locus could be ancestrally heterozygous
		if AllHetsSameGenotype=="True":
			# But we need to check the genotypes of the haploids
			HapPassCount=0
			for k in HapInd:
				# don't allow any heterozygous sites, only homozygous ones
				if Locus[k][0] == Locus[k][2]:
					# then check if the allele possessed by the haploid sample is found in the genotype of the heterozygous sample(s)
					if Locus[k][0] == Hets[0][0] or Locus[k][0] == Hets[0][2]:
						HapPassCount+=1
			# If this is true for all haploids
			if HapPassCount==len(HapInd):
				# Inspect the homozygotes. If there are none you're good to go.
				if len(Homs) == 0:
					WriteOut='Yes'
				# but if there are some we need to look further
				elif len(Homs) >= 1:
					HomPassCount=0
					# As with the haploids, check to see if the allele possessed by the homozygotes is one of the alleles in the heterozygous sample(s)'s genotypes
					for Geno in Homs:
						if Geno[0]==Hets[0][0] or Geno[0]==Hets[0][2]:
							HomPassCount+=1
					if len(Homs)==HomPassCount:
						WriteOut='Yes'
		# Write to the appropriate output file
		if WriteOut=='Yes':
			OutFile.write('\n'+Locus[0])
			for i in range(1,len(Locus)):
				OutFile.write('\t'+Locus[i])
		if WriteOut!='Yes':
			NotFile.write('\n'+Locus[0])
			for i in range(1,len(Locus)):
				NotFile.write('\t'+Locus[i])

NotFile.close()
OutFile.close()
print('Done')
