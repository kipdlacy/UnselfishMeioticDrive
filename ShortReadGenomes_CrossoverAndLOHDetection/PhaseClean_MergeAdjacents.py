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

### Reading in FIles
PShead=read_csv(PSfile)[0]
PSs=read_csv(PSfile)[1:]
print('PShead:')
print(PShead)
print('PSs[0]:')
print(PSs[0])

### OutFile
OutFile=open(PSfile.split('.txt')[0]+"_MergeAdjacents.txt",'w')
OutFile.write(PShead[0])
for i in range(1,len(PShead)):
	OutFile.write('\t'+PShead[i])

### Looping, identifying mergeable adjacents, writing out either unmerged or merged
i=0
while i < len(PSs):
	if i+1 < len(PSs):
		if PSs[i][0]==PSs[i+1][0] and PSs[i][1]==PSs[i+1][1] and PSs[i][8]==PSs[i+1][8] and PSs[i][9]==PSs[i+1][9] and PSs[i][10]==PSs[i+1][10] and PSs[i][11]==PSs[i+1][11]:
			Run=[]
			Run.append(PSs[i])
			Run.append(PSs[i+1])
			p=i+2
			RunGood='True'
			while RunGood=='True':
				if p >= len(PSs):
					RunGood='False'
				elif PSs[i][0]==PSs[p][0] and PSs[i][1]==PSs[p][1] and PSs[i][8]==PSs[p][8] and PSs[i][9]==PSs[p][9] and PSs[i][10]==PSs[p][10] and PSs[i][11]==PSs[p][11]:
					Run.append(PSs[p])
					p+=1
				else:
					RunGood='False'
			OutFile.write('\n'+PSs[i][0]+'\t'+PSs[i][1]+'\t'+PSs[i][2]+'\t'+Run[-1][3])
			NumSNPs=0
			for y in Run:
				NumSNPs+=int(y[4])
			OutFile.write('\t'+str(NumSNPs)+'\t'+str(int(Run[-1][3])-int(PSs[i][2]))+'\t'+PSs[i][6]+'\t'+Run[-1][7])
			for k in range(8,len(PSs[i])):
				OutFile.write('\t'+PSs[i][k])
			i=p
		else:
			OutFile.write('\n'+PSs[i][0])
			for k in range(1,len(PSs[i])):
				OutFile.write('\t'+PSs[i][k])
			i+=1
	else:
		OutFile.write('\n'+PSs[i][0])
		for k in range(1,len(PSs[i])):
			OutFile.write('\t'+PSs[i][k])
		i+=1

### Finish
OutFile.close()
print('Done')