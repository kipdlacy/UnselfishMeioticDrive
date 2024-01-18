#! /usr/bin/env python

### 2023FEB10
### This script is designed to read in a list of inferred phase switches, and remove any phase blocks that are the first on their contig
### This is useful in cases where you don't want to count the beginning of the contig as a crossover

### Usage
# python ./src/PhaseCleaner_MergeAdjacents.py PStable

### Import Modules
import sys

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Passing Arguments
PSfile=sys.argv[1]

### Reading in FIles
PShead=read_csv(PSfile)[0]
PSs=read_csv(PSfile)[1:]
print('PShead:')
print(PShead)
print('PSs[0]:')
print(PSs[0])

### OutFile
OutFile=open(PSfile.split('.txt')[0]+"_FirstBlocksOnContigRemoved.txt",'w')
OutFile.write(PShead[0])
for i in range(1,len(PShead)):
	OutFile.write('\t'+PShead[i])

### Looping, identifying mergeable adjacents, writing out either unmerged or merged
for i in range(0,len(PSs)):
	# If the SNP at the start of the phase block is not the first on the contig
	if PSs[i][6] == 'No':
		# write it to the output file
		OutFile.write('\n'+PSs[i][0])
		for k in range(1,len(PSs[i])):
			OutFile.write('\t'+PSs[i][k])

### Finish
OutFile.close()
print('Done')