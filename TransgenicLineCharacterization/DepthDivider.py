#! /usr/bin/env python

### 20210907
### This is a script designed to crawl through the output from samtools depth 
### and output a file with all depth values divided by a user input value 
### (typically the median or mean read depth for the entire genome)

### Usage
### python ./src/DepthDivider.py InFile Number

### Import Modules
import sys

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

# Passing command-line arguments
SNPfile=sys.argv[1] # Input file name
SNPs=read_csv(SNPfile) # Read the SNP file
Divisor=float(sys.argv[2]) # The divisor value (median or mean depth)

### OutFile
### Note: The output file is named after the original input file with "_DivBy[Divisor]" appended.
OutFile=open(SNPfile+"_DivBy"+str(Divisor),'w')

SignificantDigits=3 # Set the number of significant digits you want the normalized read depths to have.
FirstLine='True' # Just a flag to make sure there isn't an empty line at the beginning of the file
### Loop through each locus in the SNP data
for Locus in SNPs:
	# Perform the division and round the result to the desired number of significant digits
	NormalizedDepth=round(float(Locus[2])/Divisor,SignificantDigits)
	# Construct the line to be written, including dividing the depth value by the divisor 
	Line = Locus[0]+'\t'+Locus[1]+'\t'+str(NormalizedDepth)
	if FirstLine=='True':
		### Divide the depth value by the divisor and write the result to the output file
		OutFile.write(Line)
		FirstLine='False'
	elif FirstLine=='False':
		### Divide the depth value by the divisor and write the result to the output file
		### Add a new line character to the beginning
		OutFile.write('\n'+Line)

### Closing the output file
OutFile.close()
print('Done')

