#! /usr/bin/env python

### This script takes as input a file output from samtools depth, and writes out a new file including only sites where the read depth is equal to or greater than a user-input value

### import packages
import sys

### Pass the 2x genome-wide mean read depth variable in put as an argument from the command line
Num=float(sys.argv[2])

OutFile=open(sys.argv[1]+'_DepthAbove2xMean.txt','w')
OutFile.write('Chrom\tPosition\tDepth')

with open(sys.argv[1]) as File:
	for Line in File:
		### If the read depth is equal to or greater than the user-input argument
		if float(Line.split('\t')[2].split('\n')[0])>=Num:
			### Write to the output file
			OutFile.write('\n'+Line.split('\t')[0]+'\t'+Line.split('\t')[1]+'\t'+Line.split('\t')[2].split('\n')[0])

OutFile.close()
print('Done :~]')