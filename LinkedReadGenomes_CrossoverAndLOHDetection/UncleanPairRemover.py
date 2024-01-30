#! /usr/bin/env python

### 20220407
### This script is designed to read in a list of pairs of all SNPs from a VCF with annotations
### and remove all but "clean" pairs of SNPs

### Import Modules
import sys

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
OutFile=open(InFile.split('.bed')[0]+'_CleanOnly.bed','w')

### Looping through input file
for i in range(1,len(List)):
	### Checking to be sure the line isn't empty
	if len(List[i])==len(List[1]):
		### If the pair has not been associated with any "bad" features
		if List[i][3]=='None' and List[i][4]=='None':
			### Write the pair to the output file
			OutFile.write('\n'+List[i][0]+'\t'+List[i][1]+'\t'+List[i][2])

### Finishing Up
OutFile.close()
print('Done =~]')