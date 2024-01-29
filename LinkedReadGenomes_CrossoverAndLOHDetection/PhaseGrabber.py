#! /usr/bin/env python

#### This file is designed to read in an input file of consecutive pairs of Pairs (i.e., at which crossovers have been inferred)
#### and grab the genotype/phase of single individuals from input VCF-style tab-delimited files generated with GATK VariantsToTable.
#### Need to also include as an argument the name of the sample.

### Import Modules
import sys

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Pass the crossover file
### This file is a bed style file that contains each crossover
### And over time, the haplotypic phase of each sample
### Pair stands not for pairs of samples, but for pairs of sites between which a crossover is inferred to have occurred.
PairFile=sys.argv[1]
PairHead=read_csv(PairFile)[0]
print('PairHead:')
print(PairHead)
Pairs=read_csv(PairFile)[1:]
print('Pairs[0]:')
print(Pairs[0])

### Pass the sample genotype file
GTFile=sys.argv[2]
GTHead=read_csv(GTFile)[0]
print('GTHead:')
print(GTHead)
GTs=read_csv(GTFile)[1:]
print('GTs[0]:')
print(GTs[0])

### Pass the name of the sample
Sample=sys.argv[3]

### Creating Output file
OutFile=open(PairFile+"_"+Sample,'w')
OutFile.write(str(PairHead[0]))
for i in range(1,len(PairHead)):
	OutFile.write('\t'+PairHead[i])
OutFile.write('\t'+Sample)

#### Looping across pairs of Sites between which crossovers occurred
for Pair in Pairs:
	### Checking to make sure the line isn't empty
	if len(Pair)==len(Pairs[0]):
		### Initialize a list
		Store=[]
		### For each SNP in the pair (1, and 2 index the positions of SNPs)
		for i in [1,2]:
			### Looping through all genotypes in the VCF-style table
			for GT in GTs:
				### Assuming the line isn't empty, the chromosome scaffold is the same, and the SNP position matches
				if len(GT)==len(GTs[0]) and Pair[0]==GT[0] and Pair[i]==GT[1]:
					### Store the last entry in the line in the genotype file, which contains the phased genotype
					Store.append(GT[-1])
		### Now write to the output file
		OutFile.write('\n'+Pair[0])
		for k in range(1,len(Pair)):
			OutFile.write('\t'+Pair[k])
		if len(Store)<2:
			OutFile.write('\tn/a')
		else:
			OutFile.write('\t'+str(Store[0])+','+str(Store[1]))
