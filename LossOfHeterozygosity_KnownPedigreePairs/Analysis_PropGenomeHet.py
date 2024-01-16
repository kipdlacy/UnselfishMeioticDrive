#! /usr/bin/env python

# 20210629
# Script to count the number and proportion of sliding windows
# for which one or more genomic features (in this case SNPs)
# are present

### Import Modules
import sys

### Defining Functions
def read_csv(filename, l_char='\n', spl_char='\t'):
	return [line.split(spl_char) for line in open(filename).read().split(l_char)]

### Passing Arguments
SNPfile=sys.argv[1]
SNPhead=read_csv(SNPfile)[0]
SNPs=read_csv(SNPfile)[1:]
print('SNPs[0]')
print(SNPs[0])
WindHead=read_csv(sys.argv[2])[0]
Windows=read_csv(sys.argv[2])[1:]
print('Windows[0]')
print(Windows[0])

### Creating Output file and writing header
OutFile=open(sys.argv[2]+'_PropHet_'+SNPhead[-1].split('.GT')[0],'w')
OutFile.write(WindHead[0])
for i in range(1,len(WindHead)):
	OutFile.write('\t'+WindHead[i])
OutFile.write('\tNumHetSNP')

###DoubleForLoop
Total=0 # instantiate counter of total number of windows
Het=0 # instantiate counter for number of windows that contain one or more heterozygous SNP
# loop through the windows
for Wind in Windows:
	NumHet=0 # instantiate counter for number of SNPs in window
	# loop through the VCF-style table
	for SNP in SNPs:
		# Check if the SNP falls within the window
		if Wind[1]==SNP[0] and SNP[1]>=Wind[3] and SNP[1]<=Wind[4]:
			# if yes, add to the tally of SNPs in this window
			# note -- since all SNPs within the window are counted, you have to input a file that already contains
			# SNPs that are heterozygous
			NumHet+=1
	# Once you finish looping through the SNPs, write the info to the output file
	OutFile.write('\n'+Wind[0])
	for i in range(1,len(Wind)):
		OutFile.write('\t'+Wind[i])
	OutFile.write('\t'+str(NumHet))
	# Add to the tally of number of windows that contain at least one heterozygous SNP, if appropriate
	if NumHet>0:
		Het+=1
	# Always add to the tally of total number of windows
	Total+=1

OutFile.close()

# Create a summary output file
PropOut=open(SNPfile+'_Prop_'+sys.argv[2].split('bp')[0].split('/')[-1]+'bpWindowsHet','w')
PropOut.write('SampleName\tWindowSize(bp)\tWindowsHet\tTotalWindows\tProportionWindowsHet')
PropOut.write('\n'+SNPfile.split('_')[-1].split('Het')[0]+'\t'+sys.argv[2].split('bp')[0].split('/')[-1]+'\t'+str(Het)+'\t'+str(Total)+'\t'+str(float(Het)/float(Total)))
print('Windows Het: '+str(Het))
print('Total Windows: '+str(Total))
print('Proportion Windows Het: '+str(float(Het)/float(Total)))
print('done')