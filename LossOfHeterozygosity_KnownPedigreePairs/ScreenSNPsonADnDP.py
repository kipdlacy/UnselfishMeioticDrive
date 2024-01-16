#! /usr/bin/env python

# This script reads in a genotype table and an allelic depth (AD) table.
# It discards variants with heterozygous genotypes where the minor allelic depth is less than 25% of the total read depth at that variant for that sample
# and also discards variants with homozygous genotypes where the minor allelic depth is greater than 0.

# Usage
# python ./MDpair_GetDiffs.py MDpair_1_SNPs_filtered.table

# Import Modules
import sys

# Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

Head=read_csv(sys.argv[1])[0]
print('Head:')
print(Head)
Body=read_csv(sys.argv[1])[1:]
print('Body[0]:')
print(Body[0])

ADtable=read_csv(sys.argv[2])[1:]
print('ADtable[0]:')
print(ADtable[0])

OutFile=open(sys.argv[1]+'_ADDPscreenNucs','w')
OutFile.write(Head[0])
for i in range(1,len(Head)):
	OutFile.write('\t'+Head[i])

# Create a dictionary for quick lookup in ADtable
# The key is a tuple (Site[0], Site[1]), which refers to the chromosome number and position
# and the value is the entire "Site" or line for each variant in the dataset
ADtableDictionary = {(Site[0], Site[1]): Site for Site in ADtable if len(Site) == len(ADtable[0])}

# Iterate through each Variant "Line" in in the dataset "Body"
for Line in Body:
	# Screening out empty lines at the end of the file
	if len(Line) != len(Body[0]):
		continue
	# Checking to see if every sample in the table has genotype calls that are nucleotides (rather than null calls)
	# This line uses all() with a generator expression for efficiency
	if not all(Line[k][0] in 'ATGC' and Line[k][2] in 'ATGC' for k in range(5, len(Line))):
		continue # Skip lines that don't meet the character criteria
	# Look up the corresponding Site in ADtableDictionary
	Site = ADtableDictionary.get((Line[0], Line[1]))
	if not Site:
		continue # Skip if there is no matching Site
	# Now that you have both the genotype line and the AD line for this variant, we will check to see if the genotype calls are valid
	# for each sample 
	for k in range(5, len(Line)):
		# AD values are the number of reads that match each allele, and are separated by a comma like "42,30"
		# So by splitting the values by the comma we can access the read depths for each allele for each sample
		SplitADs = [int(i) for i in Site[k].split(',')]
		SumADs = sum(SplitADs)
		if SumADs < 15 or SumADs > 80:
			# Skip further processing if the sum is zero to avoid division by zero
			break
		# Sorting the individual AD values so that you know which one is the one with lower read depth (the "minor allele")
		SortedADs = sorted(SplitADs)
		# For heterozygous genotypes
		if Line[k][0] != Line[k][2]:
			# Divide the "minor allele" (i.e., the one with lower allelic depth) by the total number of reads to get the "Minor Allelic Depth" proportion
			HetMAD = float(SortedADs[-2]) / float(SumADs)
			# If this proportion is less than 0.25, then this is likely a poor quality genotype call.
			# In this case, it is likely that this site is not truly heterozygous, but that
			# some technical artifact lead to this genotype call.
			# thus, we discard this variant if that is the case
			if HetMAD < 0.25:
				break  # Stop checking further if condition is not met
		# For homozygous genotypes
		elif Line[k][0] == Line[k][2]:
			# We discard homozygous genotype calls if one or more read actually possesses the minor allele
			# that is, if the minor allele read depth is greater than zero.
			# Such sites are likely falsely called homozygous, and may have technical artifacts biasing genotype calls
			if sum(SortedADs[:-1]) > 0:
				break
	else: # Note: The 'else' block of the for loop is executed only when the loop is not terminated by a 'break'
		# Write the Line to OutFile if all conditions are met
		OutFile.write('\n'+Line[0]+'\t'+'\t'.join(Line[1:]))

OutFile.close()
print('DOne')