#! /usr/bin/env python

# 20210628
# Script designed to input a VCF-style
# table and output a bed file
# with one column for chrom
# and one for position

### Usage
### python ./src/Process_VCFtoGRangesInput.py VCFtable

### Import Modules
import sys

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Passing Arguments
SNPfile=sys.argv[1]

### Reading the SNP file header
SNPhead=read_csv(SNPfile)[0]
print('SNPhead')
print(SNPhead)
SNPs=read_csv(SNPfile)[1:]
print('SNPs[0]')
print(SNPs[0])

### Name the outfile and write the header
OutFile=open(SNPfile+".bed",'w')
OutFile.write(SNPs[0][0]+'\t'+str(int(SNPs[0][1])-1)+'\t'+SNPs[0][1])

### Loop through the file and write each line to the outfile
for Line in SNPs[1:]:
	OutFile.write('\n'+Line[0]+'\t'+str(int(Line[1])-1)+'\t'+Line[1])

### End
OutFile.close()
print('Done!')
