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

def compare_elements(list1, list2, indices):
	"""
	Compares elements of two lists at specified indices.
	This function is designed to compare elements of list1 and list2 at the indices provided in the 'indices' list.
	It will only compare up to the length of the shorter list to avoid IndexError in case one list is shorter than the other.
	indices (list of int): The indices at which to compare elements in the lists.
	Returns:
	bool: True if all compared elements are equal, False otherwise.
	"""
	for index in indices:
		# Check if the current index is within the bounds of both lists
		if index < len(list1) and index < len(list2):
			# Compare elements at the current index; if any pair doesn't match, return False
			if list1[index] != list2[index]:
				return False
		else:
			# If the index is out of bounds for either list, stop comparing further
			break
	return True

### Looping, identifying mergeable adjacents, writing out either unmerged or merged
i = 0
while i < len(PSs):
	if i + 1 < len(PSs):
		# Use the compare_elements function to compare PSs[i] and PSs[i+1] at the specified indices
		if compare_elements(PSs[i], PSs[i+1], [0, 1, 8, 9, 10, 11, 12, 13]):
			Run = []
			Run.append(PSs[i])
			Run.append(PSs[i+1])
			p = i + 2
			RunGood = 'True'
			while RunGood == 'True':
				if p >= len(PSs):
					RunGood = 'False'
				else:
					# Again, use the compare_elements function for comparison with PSs[p]
					if compare_elements(PSs[i], PSs[p], [0, 1, 8, 9, 10, 11, 12, 13]):
						Run.append(PSs[p])
						p += 1
					else:
						RunGood = 'False'
			# Writing merged data to OutFile
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
