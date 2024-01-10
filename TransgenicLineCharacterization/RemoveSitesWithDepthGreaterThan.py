#! /usr/bin/env python

### 20210907
### This is a script designed to crawl through a file with depth values
### and outputs a file with all lines except those that exceed a value
### input by the user as an argument from the command line

### Usage
### python ./src/RemoveSitesWithDepthGreaterThan.py InFile Number

### Import Modules
import sys

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Passing Arguments
SNPfile=sys.argv[1] # Input file name
SNPs=read_csv(SNPfile) # Read the SNP file
Cutoff=float(sys.argv[2]) # Cutoff depth value

### OutFile
# Note: The output file is named after the original input file with "_NonGreaterThan[Cutoff]" appended.
OutFile=open(SNPfile+"_NoneGreaterThan"+str(Cutoff),'w')

### Looping through each locus in the SNP data
for Locus in SNPs:
	# Check if depth value is less than or equal to the cutoff
	if float(Locus[2])<=Cutoff:
		# Writing the line to the output file
		OutFile.write('\n'+Locus[0]+'\t'+Locus[1]+'\t'+str(float(Locus[2])))

### Closing the output file
OutFile.close()
print('Done')

