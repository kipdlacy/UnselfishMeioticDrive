#! /usr/bin/bash

#run_tellread.sh \
#	-i <path/to/raw/data> \
#	-o <path/to/output> \
#	-f <path/to/reference> \
#	-s <comma separated sample list> \
#	-g <comma separated genome list>
	

# Had to create the genome directory according to the TellRead manual. 
# Had to port the raw data directory (entire output from Illumina run) over to the server
# Had to manually change a file titled runParameters.xml to RunParameters.xml

/home/kip/tools/tellread-release/run_tellread.sh \
	-i /path/to/raw/data/directory/220214_A00815_0328_BHVKVLDRXY/ \
	-o /path/to/output/directory/tellreadOUT/ \
	-f /path/to/GenomeReferenceDirectory/ \
	# I named by sample list by the TELLseq indexes that I used
	-s T505,T506,T507,T508,T509,T510,T511,T512,T513,T514,T515,T516 \
	# Specifying the genome used
	-g Obir,Obir,Obir,Obir,Obir,Obir,Obir,Obir,Obir,Obir,Obir,Obir


### For loop to run tell sort

# run_tellsort.sh \
# -r1 <R1 read>.fastq.gz \
# -r2 <R2 read>.fastq.gz \
# -i1 <I1 read>.fastq.gz \
# -r <genome reference file in fasta> \
# -o <path/to/output> \
# -p <prefix name> \
# -b <genome_variants>.bed \
# -v <genome_variants>.vcf.gz \
# -t number of threads

# Looping over samples
for CODE in T505 T506 T507 T508 T509 T510 T511 T512 T513 T514 T515 T516
do
echo 'starting '$CODE
/home/kip/tools/tellsort-release/run_tellsort.sh \
	-r1 /path/to/output/directory/tellreadOUT/Full/tellreadOUT_R1_$CODE.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
	-r2 /path/to/output/directory/tellreadOUT/Full/tellreadOUT_R2_$CODE.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
	-i1 /path/to/output/directory/tellreadOUT/Full/tellreadOUT_I1_$CODE.fastq.gz.corrected.fastq.err_barcode_removed.fastq.gz \
	-r /path/to/GenomeReferenceDirectory/Obir/GenomeAssemblyFilename.fasta \
	-o /path/to/output/directory/tellreadOUT/tellsortOUT/$CODE \
	-p $CODE \
	# A trusted SNP file is used for some preliminary QC. But SNPs were called separately in tellread.
	-v /path/to/KnownVariantFileForPrelimQC/TrustedSNPs.vcf \
	-t 38
done
