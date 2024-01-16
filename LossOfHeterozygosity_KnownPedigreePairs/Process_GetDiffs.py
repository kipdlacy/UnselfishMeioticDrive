#! /usr/bin/env python

# Script designed to input a VCF-style
# table of genotypes with two samples,
# and output an analogous table containing all variants
# at which these two samples have different genotypes

### Usage
### python ./src/Process_GetDiffs.py VCFtable

### Import Modules
import sys

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Passing Arguments
SNPfile=sys.argv[1]
SNPhead=read_csv(SNPfile)[0]
print('SNPhead')
print(SNPhead)
SNPs=read_csv(SNPfile)[1:]
print('SNPs[0]')
print(SNPs[0])

### OutFile
OutFile=open(SNPfile+"_Diffs",'w')

### Loop
for Line in SNPs:
	# This style of table contains the genotype of the first sample in the column indexed by 5
	# and the genotype of the second sample in the column indexed by 6
	if Line[5]!=Line[6]:
		OutFile.write('\n'+Line[0]+'\t'+Line[1])
		for i in range(2,len(Line)):
			OutFile.write('\t'+Line[i])

### End
OutFile.close()
print('Done!')
