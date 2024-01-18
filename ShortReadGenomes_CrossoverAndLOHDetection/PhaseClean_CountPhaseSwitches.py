#! /usr/bin/env python

### 20210803
### This script is designed to read in a list of inferred phase switches,
### And write to a separate file in GRanges input format only phase switches
### that is the beginnings of phase blocks that are not the first on the contig

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
OutFile=open(PSfile.split('.txt')[0]+"_PSforGRanges.txt",'w')
OutFile.write('chrom\tposition')

### Looping, identifying mergeable adjacents, writing out either unmerged or merged
for i in range(0,len(PSs)):
	# if not the first on the contig
	# we count the beginning of each phase block as a recombination event
	# and so the first SNP on a contig can't be one
	if PSs[i][6] == 'No':
		OutFile.write('\n'+PSs[i][0]+'\t'+PSs[i][2])

### Finish
OutFile.write('\n')
OutFile.close()
print('Done')