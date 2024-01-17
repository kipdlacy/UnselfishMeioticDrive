#! /usr/bin/env python

### 20210803
### This script is designed to read in a list of inferred phase switches, and merge any 
### adjacent PBs within a contig that do not differ in phase
### such might have previously been interrupted by zero length PBs.

### Usage
# python ./src/PhaseCleaner_MergeAdjacents.py PStable

### Import Modules
import sys

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Passing Arguments
PSfile=sys.argv[1]
NCOlen=int(sys.argv[2])

### Reading in FIles
PShead=read_csv(PSfile)[0]
PSs=read_csv(PSfile)[1:]
print('PShead:')
print(PShead)
print('PSs[0]:')
print(PSs[0])

### OutFile
OutFile=open(PSfile.split('.txt')[0]+"_NoPBsLessThan"+str(NCOlen)+"BP.txt",'w')
OutFile.write(PShead[0])
for i in range(1,len(PShead)):
	OutFile.write('\t'+PShead[i])

### Looping, identifying mergeable adjacents, writing out either unmerged or merged
for i in range(0,len(PSs)):
	if int(PSs[i][5]) >= NCOlen:
		OutFile.write('\n'+PSs[i][0])
		for k in range(1,len(PSs[i])):
			OutFile.write('\t'+PSs[i][k])

### Finish
OutFile.close()
print('Done')