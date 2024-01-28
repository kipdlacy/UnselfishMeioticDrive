#! /usr/bin/env python

### Read in a file output from samtools depth, and write an output file that contains only sites at which the read depth is equal to zero.

### Import sys
import sys

OutFile=open(sys.argv[1]+'_ZeroesOnly.txt','w')
OutFile.write('Chrom\tPosition\tDepth')

with open(sys.argv[1]) as File:
	for Line in File:
		### If the line indicates that this site has a read depth of zero
		if Line.split('\t')[2]=='0\n':
			### Write the line to an output file
			OutFile.write('\n'+Line.split('\t')[0]+'\t'+Line.split('\t')[1]+'\t'+Line.split('\t')[2].split('\n')[0])

OutFile.close()
print('Done :~]')