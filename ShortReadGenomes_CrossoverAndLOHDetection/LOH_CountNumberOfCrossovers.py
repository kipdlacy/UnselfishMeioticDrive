#! /usr/bin/env python

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
OutFile=open(PSfile.split('.txt')[0]+"_COsFromLOH_NonEndOfContigOnly_FullLines.txt",'w')
OutFile.write(PShead[0])
for k in range(1,len(PShead)):
	OutFile.write('\t'+PShead[k])

## Looping, identifying mergeable adjacents, writing out either unmerged or merged
for i in range(0,len(PSs)):
	# checking whether the breakpoint of the LOH is the first SNP on the contig or the last on the contig
	# if so, then we can't call it as a crossover, 
	# but if not (i.e., if PSs[i][6] or [7] == 'No') then it likely is
	if PSs[i][6] == 'No' or PSs[i][7] == 'No':
		OutFile.write('\n'+PSs[i][0])
		for k in range(1,len(PSs[i])):
			OutFile.write('\t'+PSs[i][k])

### Finish
OutFile.close()
print('Done')