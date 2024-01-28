#! /usr/bin/env python

### 20220419
### This script is designed to read in a Table with AD, and output a file that contains instead the MAF, by calculating it

### Importing modules
import sys

### Defining functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Read in variables
Table=read_csv(sys.argv[1])

### Get index for AD
for i in range(2,len(Table[0])):
	if Table[0][i].split('.')[1]=='AD':
		Index=i

### Open OutFile and write header
OutFile=open(sys.argv[1]+'_MAF','w')
OutFile.write(Table[0][0]+'\t'+Table[0][1]+'\t'+'AD'+'\t'+'MAF')

### Loop through table and write to outfile if condition is satisfied
for i in range(1,len(Table)):
	if len(Table[i])==len(Table[0]):
		List=[]
		List.append(int(Table[i][Index].split(',')[0]))
		List.append(int(Table[i][Index].split(',')[1]))
		OutFile.write('\n'+Table[i][0]+'\t'+Table[i][1]+'\t'+Table[i][Index]+'\t'+str(round(float(min(List))/sum(List), 3)))

### Close the OutFile
OutFile.close
print('Done! :~]')