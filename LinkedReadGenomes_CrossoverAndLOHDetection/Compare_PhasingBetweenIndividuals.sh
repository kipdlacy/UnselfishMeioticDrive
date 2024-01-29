#! /usr/bin/bash

### This file is designed to compare phasing among different individuals of interest after filtering
### And then filter those "phase switches", or putative crossovers, to provide a list of not obviously false positive crossovers for manual screening

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
		#######################
		######### Now #########
		###### Filtering ######
		###### Putative #######
		###### Crossovers #####
		#######################
		### Define file paths
		SwitchErrorsFile="./results/Comparisons/${SAMPLE1}_v_${SAMPLE2}/chr$CHROM'_SwitchErrors.diff.switch'"
		ScreenedSitesSAMPLE1="./data/Samples/$SAMPLE1/chr$CHROM'_AllScreenedSitesUnordered_NoIndels.txt'"
		ScreenedSitesSAMPLE2="./data/Samples/$SAMPLE2/chr$CHROM'_AllScreenedSitesUnordered_NoIndels.txt'"
		VCFtableSAMPLE1="./data/Samples/$SAMPLE1/hapblock_chr$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table'"
		VCFtableSAMPLE2="./data/Samples/$SAMPLE2/hapblock_chr$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table'"
		### first remove all crossovers in which the focal SNPs are intervened by or flanked by (within two flanking SNPs on either side) previously screened SNPs
		### such sites are likely genome assembly errors, and so the phasing cannot be trusted
		python ./src/PSscreen_BadSNPsIntervenePS.py $SwitchErrorsFile $ScreenedSitesSAMPLE1 $ScreenedSitesSAMPLE2 BadSNP Remove
		python ./src/PSscreen_BadSNPsFlankPS.py ${SwitchErrorsFile}_NoBadSNPInter.bed $ScreenedSitesSAMPLE1 $ScreenedSitesSAMPLE2 $VCFtableSAMPLE1 $VCFtableSAMPLE2 BadSNP 1 Remove
		### Remove crossovers that occur at sites that are intervened by, or are flanked by (within two flanking SNPs on either side) sites with zero read depth
		### such sites are likely genome assembly errors, and so the phasing cannot be trusted
		python ./src/PSscreen_BadSNPsIntervenePS.py ${SwitchErrorsFile}_NoBadSNPInter_NoBadSNPWN_1.bed ${ScreenedSitesSAMPLE1}_ZeroesOnly.txt ${ScreenedSitesSAMPLE2}_ZeroesOnly.txt 0cov Remove
		python ./src/PSscreen_BadSNPsFlankPS.py ${SwitchErrorsFile}_NoBadSNPInter_NoBadSNPWN_1_No0covInter.bed ${ScreenedSitesSAMPLE1}_ZeroesOnly.txt ${ScreenedSitesSAMPLE2}_ZeroesOnly.txt $VCFtableSAMPLE1 $VCFtableSAMPLE2 0cov 1 Remove
		### Remove phase switches flanked or intervened by sites with 0 coverage if read depth counts only count reads with MAPQ>=1
		### This excludes sites with MAPQ<1, for which the phase cannot be trusted, because reads may have aligned to multiple sites in the genome
		python ./src/PSscreen_BadSNPsIntervenePS.py ${SwitchErrorsFile}_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1.bed ${ScreenedSitesSAMPLE1}_MAPQ0_ZeroesOnly.txt ${ScreenedSitesSAMPLE2}_MAPQ0_ZeroesOnly.txt MAPQ0 Remove
		python ./src/PSscreen_BadSNPsFlankPS.py ${SwitchErrorsFile}_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1_NoMAPQ0Inter.bed ${ScreenedSitesSAMPLE1}_MAPQ0_ZeroesOnly.txt ${ScreenedSitesSAMPLE2}_MAPQ0_ZeroesOnly.txt $VCFtableSAMPLE1 $VCFtableSAMPLE2 MAPQ0 1 Remove
		### Then count the number of intervening and flanking sites with coverage greater than 2x the genome-wide mean in either focal sample.
		python ./src/PSscreen_BadSNPsIntervenePS.py ${SwitchErrorsFile}_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1_NoMAPQ0Inter_NoMAPQ0WN_1.bed ${ScreenedSitesSAMPLE1}_DepthAbove2xMean.txt ${ScreenedSitesSAMPLE2}_DepthAbove2xMean.txt 2xMeanCov Count
		python ./src/PSscreen_BadSNPsFlankPS.py ${SwitchErrorsFile}_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1_NoMAPQ0Inter_NoMAPQ0WN_1_2xMeanCovInterCount.bed ${ScreenedSitesSAMPLE1}_DepthAbove2xMean.txt ${ScreenedSitesSAMPLE2}_DepthAbove2xMean.txt $VCFtableSAMPLE1 $VCFtableSAMPLE2 2xMeanCov 2 Count
		### Append minor allele frequencies to each putative crossover
		python ./src/PSscreen_FetchMAF.py ${SwitchErrorsFile}_final.bed $SAMPLE1 $VCFtableSAMPLE1 $SAMPLE2 $VCFtableSAMPLE2
	done
	### For each comparison, compile all of the screened (those that weren't screened out) phase switches from each chromosome
	head -n 1 ./results/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/chr1_SwitchErrors.diff.switch_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1_NoMAPQ0Inter_NoMAPQ0WN_1_2xMeanCovInterCount_2xMeanCovWN_2_Count_MAF.bed' > ./results/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/AllChromsAllSteps.bed'
	for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14
	do
		tail -n +2 ./results/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/chr'$i'_SwitchErrors.diff.switch_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1_NoMAPQ0Inter_NoMAPQ0WN_1_2xMeanCovInterCount_2xMeanCovWN_2_Count_MAF.bed' >> ./results/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/AllChromsAllSteps.bed'
		echo '' >> ./results/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/AllChromsAllSteps.bed'
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
			#######################
			######### Now #########
			###### Filtering ######
			###### Putative #######
			###### Crossovers #####
			#######################
			### first remove all phase switches in which the focal SNPs are intervened by previously screened SNPs.
			echo "remove all phase switches in which the focal SNPs are intervened by previously screened SNPs. ./src/PSscreen_BadSNPsIntervenePS.py"
			python ./src/PSscreen_BadSNPsIntervenePS.py ./results/Comparisons/$COMP/chr$CHROM'_SwitchErrors.diff.switch' ./data/Samples/$SAMPLE1/chr$CHROM'_AllScreenedSitesUnordered_NoIndels.txt' ./data/Samples/$SAMPLE2/chr$CHROM'_AllScreenedSitesUnordered_NoIndels.txt' BadSNP Remove
			### then remove all phase switches that have a previously screened SNP within two flanking SNPs on either side
			echo "remove all phase switches that have a previously screened SNP within two flanking SNPs on either side ./src/PSscreen_BadSNPsFlankPS.py"
			python ./src/PSscreen_BadSNPsFlankPS.py ./results/Comparisons/$COMP/chr$CHROM'_SwitchErrors.diff.switch_NoBadSNPInter.bed' ./data/Samples/$SAMPLE1/chr$CHROM'_AllScreenedSitesUnordered_NoIndels.txt' ./data/Samples/$SAMPLE2/chr$CHROM'_AllScreenedSitesUnordered_NoIndels.txt' ./data/Samples/$SAMPLE1/hapblock_chr$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table' ./data/Samples/$SAMPLE2/hapblock_chr$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table' BadSNP 1 Remove
			### Next remove phase switches that are flanked by or intervened by sites with zero coverage in either of the focal samples. Need to first identify such sites using samtools depth on the focal bam files.
			echo "remove all phase switches that have a previously screened SNP within two flanking SNPs on either side ./src/PSscreen_BadSNPsIntervenePS.py ./src/PSscreen_BadSNPsFlankPS.py"
			python ./src/PSscreen_BadSNPsIntervenePS.py ./results/Comparisons/$COMP/chr$CHROM'_SwitchErrors.diff.switch_NoBadSNPInter_NoBadSNPWN_1.bed' ./data/Samples/$SAMPLE1/AllDepth.txt_ZeroesOnly.txt ./data/Samples/$SAMPLE1/AllDepth.txt_ZeroesOnly.txt 0cov Remove
			python ./src/PSscreen_BadSNPsFlankPS.py ./results/Comparisons/$COMP/chr$CHROM'_SwitchErrors.diff.switch_NoBadSNPInter_NoBadSNPWN_1_No0covInter.bed' ./data/Samples/$SAMPLE1/AllDepth.txt_ZeroesOnly.txt ./data/Samples/$SAMPLE1/AllDepth.txt_ZeroesOnly.txt ./data/Samples/$SAMPLE1/hapblock_chr$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table' ./data/Samples/$SAMPLE2/hapblock_chr$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table' 0cov 1 Remove
			### Remove phase switches flanked or intervened by sites with 0 coverage if read depth counts only count reads with MAPQ>=1
			echo "Remove phase switches flanked or intervened by sites with 0 coverage if read depth counts only count reads with MAPQ>=1 ./src/PSscreen_BadSNPsIntervenePS.py ./src/PSscreen_BadSNPsFlankPS.py"
			python ./src/PSscreen_BadSNPsIntervenePS.py ./results/Comparisons/$COMP/chr$CHROM'_SwitchErrors.diff.switch_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1.bed' ./data/Samples/$SAMPLE1/AllQualDepth.txt_ZeroesOnly.txt ./data/Samples/$SAMPLE1/AllQualDepth.txt_ZeroesOnly.txt MAPQ0 Remove
			python ./src/PSscreen_BadSNPsFlankPS.py ./results/Comparisons/$COMP/chr$CHROM'_SwitchErrors.diff.switch_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1_NoMAPQ0Inter.bed' ./data/Samples/$SAMPLE1/AllQualDepth.txt_ZeroesOnly.txt ./data/Samples/$SAMPLE1/AllQualDepth.txt_ZeroesOnly.txt ./data/Samples/$SAMPLE1/hapblock_chr$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table' ./data/Samples/$SAMPLE2/hapblock_chr$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table' MAPQ0 1 Remove
			### Then count the number of intervening and flanking sites with coverage greater than 2x the genome-wide mean in either focal sample.
			echo "count the number of intervening and flanking sites with coverage greater than 2x the genome-wide mean in either focal sample. ./src/PSscreen_BadSNPsIntervenePS.py ./src/PSscreen_BadSNPsFlankPS.py"
			python ./src/PSscreen_BadSNPsIntervenePS.py ./results/Comparisons/$COMP/chr$CHROM'_SwitchErrors.diff.switch_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1_NoMAPQ0Inter_NoMAPQ0WN_1.bed' ./data/Samples/$SAMPLE1/AllDepth.txt_DepthAbove2xMean.txt ./data/Samples/$SAMPLE1/AllDepth.txt_DepthAbove2xMean.txt 2xMeanCov Count
			python ./src/PSscreen_BadSNPsFlankPS.py ./results/Comparisons/$COMP/chr$CHROM'_SwitchErrors.diff.switch_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1_NoMAPQ0Inter_NoMAPQ0WN_1_2xMeanCovInterCount.bed' ./data/Samples/$SAMPLE1/AllDepth.txt_DepthAbove2xMean.txt ./data/Samples/$SAMPLE1/AllDepth.txt_DepthAbove2xMean.txt ./data/Samples/$SAMPLE1/hapblock_chr$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table' ./data/Samples/$SAMPLE2/hapblock_chr$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table' 2xMeanCov 2 Count
			### Append minor allele frequencies to each putative phase switch
			echo "Append minor allele frequencies to each putative phase switch ./src/PSscreen_BadSNPsIntervenePS.py ./src/PSscreen_BadSNPsFlankPS.py"
			python ./src/PSscreen_FetchMAF.py ./results/Comparisons/$COMP/chr$CHROM'_SwitchErrors.diff.switch_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1_NoMAPQ0Inter_NoMAPQ0WN_1_2xMeanCovInterCount_2xMeanCovWN_2_Count.bed' ${SAMPLE1#*_} ./data/Samples/FEB_T508/hapblock_chr$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table_AD_MAF' ${SAMPLE2#*_} ./data/Samples/FEB_T513/hapblock_chr$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table_AD_MAF'
		done
		### For each comparison, compile all of the screened (those that weren't screened out) phase switches from each chromosome
		head -n 1 ./results/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/chr1_SwitchErrors.diff.switch_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1_NoMAPQ0Inter_NoMAPQ0WN_1_2xMeanCovInterCount_2xMeanCovWN_2_Count_MAF.bed' > ./results/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/AllChromsAllSteps.bed'
		for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14; do
			tail -n +2 ./results/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/chr'$i'_SwitchErrors.diff.switch_NoBadSNPInter_NoBadSNPWN_1_No0covInter_No0covWN_1_NoMAPQ0Inter_NoMAPQ0WN_1_2xMeanCovInterCount_2xMeanCovWN_2_Count_MAF.bed' >> ./results/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/AllChromsAllSteps.bed'
			echo '' >> ./results/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/AllChromsAllSteps.bed'
		done
	done
done
