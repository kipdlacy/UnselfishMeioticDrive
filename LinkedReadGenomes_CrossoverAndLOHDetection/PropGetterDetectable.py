#! /usr/bin/env python

### 20220508
### Writing this script to read in a bed file, sum the number of bp included, and divide by the length of the genome composed of chromosomes.

### Import Modules
import sys

### Define a function 'read_csv' to read a CSV file and return its contents as a list of lists.
### Each inner list represents a line in the file, split by tabs.
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Retrieve the input file name from command-line arguments
InFile=sys.argv[1]
### Read the input file using the read_csv function
List=read_csv(InFile)

### Initialize a variable to count the number of base pairs
BpNum=0

### Iterate over each line in the input file
for Line in List:
	### Check if the line has the same number of elements as the header (to avoid processing malformed lines)
	if len(Line)==len(List[0]):
		### Calculate the length of the region described in this line and add it to the total count
		### The region length is calculated as (end position - start position + 1)
		BpNum+=((int(Line[2])-int(Line[1]))+1)

### Define the total size of the reference genome (in base pairs)
GenomeSize=221797443

### Open an output file to write the results
OutFile=open(InFile+"_PropGenomeDetectable.txt",'w')
OutFile.write('Number of Bp Detectable: '+str(BpNum))
OutFile.write('Number of Bp in Reference Genome (Chroms Only): '+str(GenomeSize))
OutFile.write('Proportion of Genome CO Detectable: '+str(float(BpNum)/float(GenomeSize)))

### Print the name of the input file and the proportion of the genome that is detectable
print(InFile.split('/')[4]+'     '+str(float(BpNum)/float(GenomeSize)))