#! /usr/bin/env python

### 20220407
### This script is designed to read in a list of pairs of all SNPs from a VCF with annotations
### and output a bed file with features corresponding to the boundaries of contiguous runs of "clean" pairs of SNPs

### Import Modules
import sys
import numpy

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Passing Arguments
InFile=sys.argv[1]
List=read_csv(InFile)
print('List[0]:')
print(List[0])
print('List[1]:')
print(List[1])

### Opening Output File
OutFile=open(InFile+'_ContiguousRuns.bed','w')
OutFile.write('Chrom\tStart\tEnd')
OutFile.write('\n'+List[1][0]+'\t'+List[1][1])

### Looping through input file
### Initialize a variable 'Run' to keep track of whether a contiguous run is ongoing
Run='On'
### Loop through each line in the input file, starting from the second line
for i in range(2,len(List)):
	### Check if the current line has the same number of columns as the second line
	### this is to make sure we aren't operating on an empty row, which sometimes occurs at the end of files
	if len(List[i])==len(List[1]):
		### If run is off
		if Run=='Off':
			### if the position of this site is one plus the position of the previous site (i.e., they are adjacent)
			if int(List[i][1]) == int(List[i-1][1])+1:
				### Start a new run
				Run='On'
				### and write the start position of the new contiguous block to the output file
				OutFile.write('\n'+List[i-1][0]+'\t'+List[i-1][1])
				### If this is the last line in the file,
				if i+1==len(List):
					### write the end position of the block
					OutFile.write('\t'+List[i][1])
		### If run is on
		elif Run=='On':
			### If this site is on a different chromosome than the previous, or the site isn't one plus the previous site
			if List[i][0] != List[i-1][0] or int(List[i][1]) != int(List[i-1][1])+1:
				### End the current run
				Run='Off'
				### and write the end position of the contiguous block
				OutFile.write('\t'+List[i-1][1])

OutFile.close()
print('Done =~]')