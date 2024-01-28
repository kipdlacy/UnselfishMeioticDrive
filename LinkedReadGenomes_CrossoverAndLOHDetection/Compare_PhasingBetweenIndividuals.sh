#! /usr/bin/bash

### This file is designed to compare phasing among different individuals of interest after filtering

#####################################################
###                                               ###
###   Known-pedigree pairs and negative control   ###
###                                               ###
#####################################################

### Define an array of sample pairs for comparison. Each pair is separated by a space.
### Sample pairs are formatted as 'Sample1_V_Sample2'.
declare -a SamplePairs=("NOV_T501_V_NOV_T503" "NOV_T504_V_NOV_T505" "NOV_T505_V_NOV_T506" "FEB_T506_V_FEB_T505" "FEB_T506_V_FEB_T507" "FEB_T508_V_FEB_T509" "FEB_T508_V_FEB_T510" "FEB_T511_V_FEB_T512" "FEB_T513_V_FEB_T514" "FEB_T515_V_FEB_T516")
### Iterate over each pair of samples in the array using a loop.
for Pair in "${SamplePairs[@]}"; do
	### The 'read' command splits each pair into separate variables based on the underscore delimiter.
	### 'SAMPLE1' and 'SAMPLE2' will store the names of the two samples to be compared.
	### Extract sample1 and sample2 from the pair string.
	SAMPLE1="${Pair%%_V_*}"  # Remove everything after "_V_" to get the first sample name
	SAMPLE2="${Pair##*_V_}"  # Remove everything before and including "_V_" to get the second sample name
	### Echo statement for monitoring progress.
	echo "Starting to compare $SAMPLE1 to $SAMPLE2"
	### A second loop iterates over chromosome numbers from 1 to 14.
	for CHROM in {1..14}; do
		### Construct file paths for the VCF files of each sample for the current chromosome.
		### The paths are based on a consistent directory structure and file naming convention.
		File1="./data/Samples/${SAMPLE1}/hapblock_chr${CHROM}.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf"
		File2="./data/Samples/${SAMPLE2}/hapblock_chr${CHROM}.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf"
		### Execute the vcftools command to compare the VCF files of the two samples.
		### The --diff-switch-error option looks for differences in phasing between the two files
		### each such difference is a putative crossover event
		mkdir -p "./results/Comparisons/${SAMPLE1}_v_${SAMPLE2}/"
		vcftools --vcf "$File1" --diff "$File2" --diff-switch-error --out "./results/Comparisons/${SAMPLE1}_v_${SAMPLE2}/chr${CHROM}_SwitchErrors"
	done
done

################
##
##
#### Performing all pairwise comparisons between representative samples from each known pedigree
##
##
################

### This is an array of all the sample names.
Samples=("FEB_T506" "FEB_T508" "NOV_T504" "FEB_T511" "FEB_T513" "NOV_T501" "FEB_T515")
### The outer loop iterates over all the sample names.
### 'i' is the index of the current sample in the samples array.
for i in "${!Samples[@]}"; do
	### The inner loop starts from the next sample (i+1) and goes till the last sample.
	### This way, we ensure we only do one-way comparisons (A vs B but not B vs A).
	for j in $(seq $((i+1)) $((${#Samples[@]}-1))); do
		### Fetch the sample names based on the current indices.
		SAMPLE1=${Samples[$i]}
		SAMPLE2=${Samples[$j]}
		### Print the sample names to be compared.
		echo "Comparing: $SAMPLE1 vs $SAMPLE2"
		### Uncomment the code below once you confirm the pairs are correct.
		mkdir -p "./results/Comparisons/${SAMPLE1}_v_${SAMPLE2}"
		### Loop through the chromosomes (1 to 14) for each pairwise comparison.
		for CHROM in {1..14}; do
		 vcftools \
		 --vcf ./data/Samples/$SAMPLE1/hapblock_chr$CHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode'_SNPsOnly.vcf'_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		 --diff ./data/Samples/$SAMPLE2/hapblock_chr$CHROM.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode'_SNPsOnly.vcf'_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		 --diff-switch-error \
		 --out ./results/Comparisons/${SAMPLE1}_v_${SAMPLE2}/chr$CHROM'_SwitchErrors'
		done
	done
done
