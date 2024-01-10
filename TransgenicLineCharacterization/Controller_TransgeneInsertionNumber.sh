#! /usr/bin/bash

# The goal of this script is to help determine the number of times a transgene inserted into the genome
# This is accomplished by sequencing the genome of an individual with a transgene using short read whole genome DNA sequencing
# And separately aligning that data to the species' reference genome, and the sequence of the transgene insert as a "reference genome"
# These alignments are done as shown in the "AlignmentAndVariantCallingWGS" directory using the "Commands_AlignAndProcess.sh" script.
# For the data aligned to the species' reference genome, it should be included in a variant calling run
# As described in the "AlignmentAndVariantCallingWGS" directory using the "Commands_VariantCallAndPrelimFilter.sh" script.

# Here, we illustrate processing downstream of those steps, which include filtering variant calls to identify trusted diploid heterozygous sites
# And then finding the read depth at those sites, and normalizing by the genome wide median read depth, which is found in the "AlignmentAndVariantCallingWGS" directory using "ObtainReadDepthDetails.sh".
# After also finding the read depth at all sites along the transgene insert, we produce a file that contains the read depth at all sites in the transgene insert,
# and the same number of randomly selected trusted diploid sites. This is then used downstream by an R script to plot the data and identify the number of insertions.

# Defining paths to programs and files at the beginning
GATK=/Path/To/gatk-4.2.0.0/gatk
REF=/Path/To/Reference/Genome/Obir.assembly.v5.4.fasta
AllSampleVCF=/Path/To/Variant/Call/File/For/AllSamples/AllIndividualAntsToDate.vcf

# First we need to align the raw DNA sequencing data (as shown in the "AlignmentAndVariantCallingWGS" directory using the "Commands_AlignAndProcess.sh" script)
# from a transgenic ant to both the reference genome for the species, and a "reference genome" of the transgene insert sequence

# The individual genome from the transgenic marker line that I used to determine the number of insertions (based on relative read depth)
# was "ShortReadPedigree1_Mother", which had an internal sample code of "MD-B-1-M". To avoid the confusion of these long sample names
# I record the process here in a for loop so that we can simply use the variable "SAMPLE".
for SAMPLE in MD-B-1-M; do
	### First keep only the relevant sample
	$GATK SelectVariants \
		-R $REF \
		-V $AllSampleVCF \
		--sample-name ${SAMPLE} \
		-O ./data/AllIndividualAntsToDate_${SAMPLE}.vcf

	### Exclude non-variant sites and keep only SNPs
	$GATK SelectVariants \
		-R $REF \
		-V ./data/AllIndividualAntsToDate_${SAMPLE}.vcf \
		--exclude-non-variants \
		--select-type-to-include SNP \
		-O ./data/AllIndividualAntsToDate_${SAMPLE}_VariantSNPs.vcf

	### Then write variants to table
	$GATK VariantsToTable \
		-V ./data/AllIndividualAntsToDate_${SAMPLE}_VariantSNPs.vcf \
		-F CHROM -F POS -F REF -F ALT -F FILTER \
		-GF GT \
		-O ./data/AllIndividualAntsToDate_${SAMPLE}_VariantSNPs.table

	### Remove all sites that are not heterozygous, and any sites at which there is a null genotype call
	python ./CleanData_KeepHetSitesRemoveNullCalls.py ./data/AllIndividualAntsToDate_${SAMPLE}_VariantSNPs.table

	### Convert VCF-style table to bedfile for samtools depth input
	python ./Process_VCFtoBed.py ./data/AllIndividualAntsToDate_${SAMPLE}_VariantSNPs.table_HetsOnlyNoNull

	### Extract read depth at het-non-zero sites using samtools depth
	samtools depth -b ./data/AllIndividualAntsToDate_${SAMPLE}_VariantSNPs.table_HetsOnlyNoNull.bed /Path/To/Aligned/Sample/File/${SAMPLE}_Trimmed_aligned_sorted_merged_dedup.bam >> ./data/AllIndividualAntsToDate_${SAMPLE}_VariantSNPs.table_HetsOnlyNoNull.bed_Depths

	### Then normalize that read depth by dividing by the genome wide MEDIAN read depth.
	### This value is obtained in the "AlignmentAndVariantCallingWGS" directory using ObtainReadDepthDetails.sh, and is entered here manually from the command line as the second argument
	python ./DepthDivider.py ./data/AllIndividualAntsToDate_${SAMPLE}_VariantSNPs.table_HetsOnlyNoNull.bed_Depths 41

	### Remove sites with normalized read depth greater than 2 to avoid erroneously assembled genomic regions
	python ./RemoveSitesWithDepthGreaterThan.py ./data/AllIndividualAntsToDate_${SAMPLE}_VariantSNPs.table_HetsOnlyNoNull.bed_Depths_DivBy41.0 2

	### NEXT for the data aligned to the transgene insert sequence, we will get the read depth at all sites using "samtools depth -aa" 
	samtools depth -aa /Path/To/Reads/Aligned/To/Transgene/Insert/Alignment/File/${SAMPLE}_Trimmed_aligned_sorted_merged_dedup.bam >> ./data/${SAMPLE}_Insert_RawDepth

	### and then normalize by dividing the raw depth by the genome-wide median read depth (as obtained in the "AlignmentAndVariantCallingWGS" directory using ObtainReadDepthDetails.sh)
	python ./DepthDivider.py ./data/${SAMPLE}_Insert_RawDepth.txt 41
	python ./SelectRandomSitesAndPrepPlottingFile.py ./data/${SAMPLE}_Insert_RawDepth_DivBy41.0 ./data/AllIndividualAntsToDate_MD-B-1-M_VariantSNPs.table_HetsOnlyNoNull.bed_Depths_DivBy41.0_NoneGreaterThan2.0

done

### After this, the plots can be produced using "TransgeneInsertionNumberViaDepthPlot.R"
