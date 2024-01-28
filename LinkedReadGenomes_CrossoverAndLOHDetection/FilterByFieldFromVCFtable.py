#! /usr/bin/env python

### 20220418
### This program is designed to for loop through a table output by GATK variants to table from a vcf file, and output only variants that pass a user-input criterion

### Import Packages
import sys

### Define functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Read in variables
Table=read_csv(sys.argv[1])
Field=sys.argv[2]
Operand=sys.argv[3]
Number=int(sys.argv[4])

### Get index for the Field of interest
for i in range(2,len(Table[0])):
	if Table[0][i].split('.')[1]==Field:
		Index=i

### Open OutFile and write header
OutFile=open(sys.argv[1]+'_'+Field+Operand+str(Number),'w')
OutFile.write(Table[0][0]+'\t'+Table[0][1])

### Loop through table and write to outfile if condition is satisfied
for i in range(1,len(Table)):
	if len(Table[i])==len(Table[0]):
		WriteOut='No'
		if Table[i][Index]!='NA':
			if Operand=="LessThan":
				if int(Table[i][Index]) < Number:
					WriteOut='Yes'
			elif Operand=="GreaterThan":
				if int(Table[i][Index]) > Number:
					WriteOut='Yes'
			elif Operand=="EqualTo":
				if int(Table[i][Index]) == Number:
					WriteOut='Yes'
		if WriteOut=='Yes':
			OutFile.write('\n'+Table[i][0]+'\t'+Table[i][1])

### Close the OutFile
OutFile.close
print('Done! :~]')