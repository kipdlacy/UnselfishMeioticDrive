#! /usr/bin/env python

### 20220409
### Script designed to read in a list of pairs of SNPs in a bed file format from multiple chromosomes and spit out all SNPs

### Import the 'sys' module to access command-line arguments
import sys

### Define a function 'read_csv' to read a CSV file and return its contents as a list of lists.
### Each inner list represents a line in the file, split by tabs.
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Retrieve the input file name from command-line arguments
InFile=sys.argv[1]

### Read the input file using the read_csv function
List=read_csv(InFile)

### Open a new file for output, naming it by modifying the input file's name
OutFile=open(InFile.split('.bed')[0]+'_UniqueSNPs_gRanges.txt','w')

### Write the header line to the output file
OutFile.write('Chrom\tposition')

### Loop through each chromosome
for CHROM in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14']:
	print('Starting '+CHROM)
	### Initialize a set to store unique SNP positions for the current chromosome
	ChromSet=[]
	### Iterate over each line in the input file
	for Line in List:
		### Check if the line corresponds to the current chromosome
		if Line[0]==CHROM:
			### Add the start and end positions of the SNP pair to the set
			ChromSet.append(int(Line[1]))
			ChromSet.append(int(Line[2]))
	# Iterate over each unique SNP position in sorted order
	for SNP in sorted(set(ChromSet)):
		### Write the chromosome and SNP position to the output file
		OutFile.write('\n'+CHROM+'\t'+str(SNP))

### Close the output file
OutFile.close()
print('Done =~]')