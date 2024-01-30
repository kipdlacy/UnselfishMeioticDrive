#! /usr/bin/env python

### 20211208
### Python script to count the number of sites at putative phase switches for which a "BAD SITE RANGE" of some kind overlaps with the two focal SNPs, or between either focal SNP and the
### SNP two upstream or two downstream.

### Importing sys to read command-line input arguments
import sys

### Defining a function to read files
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

#####Check if two ranges (start1, end1) and (start2, end2) overlap.
def overlap(start1, end1, start2, end2):
	### Returns: bool: True if the ranges overlap, False otherwise.
	return not (end1 < start2 or end2 < start1)

### Passing arguments
### A file of pairs of SNPs
PS=read_csv(sys.argv[1])
### A file of known sites indicative of poor genome assembly
Dirty1=read_csv(sys.argv[2])[1:]

### Deciding what to name the output file and then naming it
if '.bed' in sys.argv[1]:
	OutFile=open(sys.argv[1].split('.bed')[0]+'_'+sys.argv[3]+'InterveneOrFlank.bed','w')
else:
	OutFile=open(sys.argv[1]+'_'+sys.argv[3]+'InterveneOrFlank.bed','w')

### Writing the header of the output file
OutFile.write(PS[0][0])
for i in range(1,len(PS[0])):
	if PS[0][i]!='INDV':
		OutFile.write('\t'+PS[0][i])
OutFile.write('\t'+sys.argv[3]+'InterveneOrFlank')

### Initiating count variable
Count=0
### Loop through the table of pairs of SNPs
for p in range(2,len(PS)-1):
	### Checking if it's an empty line
	if len(PS[p])==len(PS[0]):
		### Initiating count variable
		Count=0
		### If it's the first iteration
		if p-1==0:
			### Loop through "Dirty" ranges
			for k in Dirty1:
				### If the line's not empty and the dirty range and the current SNP pair are on the same chromosome
				if len(k)==len(Dirty1[0]) and k[0]==PS[p][0]:
					### Ask if the current "Dirty" range overlaps with the current SNP pair
					##### WARNING, this logic is equivalent to just asking if the bad SNPs intervene the two focal SNPs
					##### so we are not asking if it flanks. If we truly wanted to ask if it flanked or intervened, we would want to do "if k[0]==PS[p][0] and int(k[1]) > int(PS[p-1][1]) and int(k[1]) < int(PS[p+1][2]):""
					if overlap(int(k[1]),int(k[2]),int(PS[p][1]),int(PS[p+1][1]))==True:
						### Add to the count
						Count+=1
		### If it's the last iteration
		elif p+1==len(PS):
			### Loop through "Dirty" ranges
			for k in Dirty1:
				### If the line's not empty and the dirty range and the current SNP pair are on the same chromosome
				if len(k)==len(Dirty1[0]) and k[0]==PS[p][0]:
					### Ask if the current "Dirty" range overlaps with the current SNP pair
					##### WARNING, this logic is equivalent to just asking if the bad SNPs intervene the two focal SNPs
					##### so we are not asking if it flanks. If we truly wanted to ask if it flanked or intervened, we would want to do "if k[0]==PS[p][0] and int(k[1]) > int(PS[p-1][1]) and int(k[1]) < int(PS[p+1][2]):""
					if overlap(int(k[1]),int(k[2]),int(PS[p-1][2]),int(PS[p][2]))==True:
						### Add to the count
						Count+=1
		### If it's not the first or the last iteration
		else:
			### Loop through "Dirty" ranges
			for k in Dirty1:
				### If the line's not empty and the dirty range and the current SNP pair are on the same chromosome
				if len(k)==len(Dirty1[0]) and k[0]==PS[p][0]:
					### Ask if the current "Dirty" range overlaps with the current SNP pair
					##### WARNING, this logic is equivalent to just asking if the bad SNPs intervene the two focal SNPs
					##### so we are not asking if it flanks. If we truly wanted to ask if it flanked or intervened, we would want to do "if k[0]==PS[p][0] and int(k[1]) > int(PS[p-1][1]) and int(k[1]) < int(PS[p+1][2]):""
					if overlap(int(k[1]),int(k[2]),int(PS[p-1][2]),int(PS[p+1][1]))==True:
						### Add to the count
						Count+=1
		### Write to the output file
		OutFile.write('\n'+PS[p][0])
		for j in range(1, len(PS[p])):
			if PS[0][j]!='INDV':
				OutFile.write('\t'+PS[p][j])
		if Count>0:
			OutFile.write('\t'+str(Count))
		else:
			OutFile.write('\tNone')
OutFile.close()
print('Done :~]')