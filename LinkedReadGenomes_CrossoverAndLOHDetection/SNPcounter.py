#! /usr/bin/env python

### 20220409
### This script is designed to read in a bedfile style table of contiguous runs with the number of SNPs in the fourth column
### and sum all the fourth columns in the file
### and print out the sum.

### Import Modules
import sys

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Passing Arguments
InFile=sys.argv[1]
List=read_csv(InFile)
# print('List[0]:')
# print(List[0])
# print('List[1]:')
# print(List[1])

### Looping through input file
### Initiate the counter
Count=0
### For each line in the file
for i in range(1,len(List)):
	### Be sure that the line's not empty
	if len(List[i])==len(List[1]):
		### Add to the counter
		Count+=int(List[i][3])

### Finishing Up
print(InFile.split('/')[4]+' '+str(Count))
# print('Done =~]')