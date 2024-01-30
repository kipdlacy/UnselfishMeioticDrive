#! /usr/bin/env python

### 20220407
### This script is designed to read in a list of pairs of all SNPs from a VCF with annotations
### and output a bed file with features corresponding to the boundaries of contiguous runs of "clean" pairs of SNPs

### Import Modules
import sys

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Passing Arguments
InFile=sys.argv[1]
List=read_csv(InFile)


### Opening Output File
OutFile=open(InFile.split('.bed')[0]+'_ContiguousCleanPairs.bed','w')

### Looping through input file

### Instantiating "Run" variable which indicates whether a "Run" of clean contiguous pairs is currently ongoing
Run='Off'
### Instantiating a counter for the number of clean pairs in that run of pairs
Count=0
### Instantiating a counter for the total number of clean pairs
TotalCount=0

### Loop through all pairs in the file
for i in range(1,len(List)):
	### Checking if the list is empty
	if len(List[i])==len(List[1]):
		### If there is no ongoing run
		if Run=='Off':
			### If no bad features have been associated with this pair of SNPs
			if List[i][3]=='None' and List[i][4]=='None':
				### Turn the run 'ON'
				Run='On'
				### Update the count and total count
				Count+=1
				TotalCount+=1
				### Initiate the outfile entry for the run
				OutFile.write('\n'+List[i][0]+'\t'+List[i][1])
				### If this is the end of the file
				if i+1==len(List):
					### Finish writing the line
					OutFile.write('\t'+List[i][2]+'\t'+str(Count+1))
		### If there is an ongoing run
		elif Run=='On':
			### If bad features have been associated with this pair of SNPs
			if List[i][3]!='None' or List[i][4]!='None':
				### Turn off the run
				Run='Off'
				### Finish the outfile entry for the just-ended run
				OutFile.write('\t'+List[i-1][2]+'\t'+str(Count+1))
				### Reset the counter
				Count=0
			### If no bad features have been associated with this pair of SNPs
			elif List[i][3]=='None' and List[i][4]=='None':
				### Add to the counter and the total count
				Count+=1
				TotalCount+=1
				### If this is the end of the file
				if i+1==len(List):
					### Finish the outfile entry for the just-ended run
					OutFile.write('\t'+List[i][2]+'\t'+str(Count+1))


### Finishing Up
OutFile.close()
# print(str(TotalCount))
# print('Done =~]')