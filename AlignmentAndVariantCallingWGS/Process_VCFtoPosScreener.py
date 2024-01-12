#! /usr/bin/env python

# 20210628
# Modified for a bed output

# 20210628
# Script designed to input a VCF-style
# table and output a GRanges-compatible
# input table with one column for chrom
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
SNPhead=read_csv(SNPfile)[0]
print('SNPhead')
print(SNPhead)
SNPs=read_csv(SNPfile)[1:]
print('SNPs[0]')
print(SNPs[0])

### OutFile
OutFile=open(SNPfile+".PosScreener",'w')
OutFile.write('CHROM\tPOS')

### Loop
for Line in SNPs:
	OutFile.write('\n'+Line[0]+'\t'+Line[1])

### End
OutFile.close()
print('Done!')
