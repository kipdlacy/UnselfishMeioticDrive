#! /usr/bin/bash

### 20220506
### This script is designed to call analyses used to obtain the proportion of the genome for each sample for which crossovers are detectable using the TELL-Seq method

### Defining an array of sample names (internal codes used)
declare -a Samples=("FEB_T505" "FEB_T506" "FEB_T507" "FEB_T508" "FEB_T509" "FEB_T510" "FEB_T511" "FEB_T512" "FEB_T513" "FEB_T514" "FEB_T515" "FEB_T516" "NOV_T501" "NOV_T503" "NOV_T504" "NOV_T505" "NOV_T506")

### For each sample
for SAMPLE in "${Samples[@]}"; do
	### Print that you're starting to process that sample
	echo 'Starting '$SAMPLE
	### Identifying runs of zero coverage sites (if you only consider reads with MAPQ>=1)
	python ./src/ContiguousZeroSites.py ./data/Samples/$SAMPLE'/AllQualDepth.txt_ZeroesOnly.txt'
	### Make a directory to hold the data for each sample
	mkdir -p ./results/Detectabilities/Samples/$SAMPLE
	### Write the header of a new file that will contain all pairs of SNPs at which we could possibly detect a crossover
	### without the result being obfuscated by genome assembly errors and errors with alignment
	# # rm ./results/Detectabilities/Samples/$SAMPLE/AllChroms_CleanContiguousPairs.bed # Optional line to remove such files from previous versions
	echo 'Chrom\tStart\tEnd' > ./results/Detectabilities/Samples/$SAMPLE/AllChroms_CleanContiguousPairs.bed
	### Also create one to hold all pairs of SNPs
	# # rm ./results/Detectabilities/Samples/$SAMPLE/AllChroms_AllPairs.bed # Optional line to remove such files from previous versions
	echo 'Chrom\tStart\tEnd' > ./results/Detectabilities/Samples/$SAMPLE/AllChroms_AllPairs.bed
	for CHROM in {1..14}; do
		echo 'chr'$CHROM
		### Turn all possible pairs of adjacent SNPs into start and end positions in bedfiles. Each pair of adjacent SNPs is a possible place at which a crossover could be detected using our method
		python ./src/TurnVCFPositionTableToPairs.py ./data/Samples/$SAMPLE'/hapblock_chr'$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table'
		### Copy that file to our new directory
		cp ./data/Samples/$SAMPLE'/hapblock_chr'$CHROM'.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf.table_Pairs' ./results/Detectabilities/Samples/$SAMPLE'/chr'$CHROM'.pairs'
		### Write all such pairs to an output file
		tail -n +2 ./results/Detectabilities/Samples/$SAMPLE'/chr'$CHROM'.pairs_BadSNPsInterveneOrFlank_0covMAPQ0InterveneOrFlank.bed' >> ./results/Detectabilities/Samples/$SAMPLE/AllChroms_AllPairs.bed
		echo '' >> ./results/Detectabilities/Samples/$SAMPLE/AllChroms_AllPairs.bed
		### Screen out pairs that are flanked by or intervened by SNPs previously identified as low quality ones that are likely associated with genomic regions for which there is no possibility of confidently identifying a crossover
		python ./src/ScreenPairs_FlankOrIntervene.py ./results/Detectabilities/Samples/$SAMPLE'/chr'$CHROM'.pairs' ./data/Samples/$SAMPLE'/chr'$CHROM'_AllScreenedSitesUnordered_NoIndels.txt' BadSNPs
		### Screen out pairs flanked by or intervened by sites with zero coverage of reads with MAPQ>=1
		python ./src/ScreenPairs_FlankOrInterveneRanges.py ./results/Detectabilities/Samples/$SAMPLE'/chr'$CHROM'.pairs_BadSNPsInterveneOrFlank.bed' ./data/Samples/$SAMPLE'/AllQualDepth.txt_ZeroesOnly.txt_ContiguousRuns.bed' 0covMAPQ0
		### Identify runs of contiguous pairs that have passed screening
		python ./src/ContiguousCleanPairGetter.py ./results/Detectabilities/Samples/$SAMPLE'/chr'$CHROM'.pairs_BadSNPsInterveneOrFlank_0covMAPQ0InterveneOrFlank.bed'
		### Write these contiguous clean pairs to an output file
		cat ./results/Detectabilities/Samples/$SAMPLE'/chr'$CHROM'.pairs_BadSNPsInterveneOrFlank_0covMAPQ0InterveneOrFlank_ContiguousCleanPairs.bed' >> ./results/Detectabilities/Samples/$SAMPLE/AllChroms_CleanContiguousPairs.bed
		echo '' >> ./results/Detectabilities/Samples/$SAMPLE/AllChroms_CleanContiguousPairs.bed
	done
	### Count the number of "trustworthy" SNPs
	python ./src/SNPcounter.py ./results/Detectabilities/Samples/$SAMPLE'/AllChroms_CleanContiguousPairs.bed'
	### Remove all unclean pairs from the "All Pairs" file
	python ./src/UncleanPairRemover.py ./results/Detectabilities/Samples/$SAMPLE'/AllChroms_AllPairs.bed'
	### Get the proportion of bp in the reference genome that, based on this analysis, have the possibility that we would be able to detect a crossover that occurred there
	python ./src/PropGetterPhased.py ./results/Detectabilities/Samples/$SAMPLE'/AllChroms_CleanContiguousPairs.bed'
done

#####################################################
###   Known-pedigree pairs and negative control   ###
#####################################################
declare - Comparisons=('FEB_T506_v_FEB_T505" "FEB_T506_v_FEB_T507" "FEB_T506_v_FEB_T511" "FEB_T508_v_FEB_T509" "FEB_T508_v_FEB_T510" "FEB_T508_v_FEB_T513" "FEB_T511_v_FEB_T512" "FEB_T513_v_FEB_T514" "FEB_T515_v_FEB_T516" "NOV_T501_v_NOV_T503" "NOV_T504_v_NOV_T501" "NOV_T504_v_NOV_T505" "NOV_T505_v_NOV_T506')
for COMP in "${Comparisons[@]}"; do
	mkdir -p ./results/Detectabilities/Comparisons/$COMP
	### Extract SAMPLE1 and SAMPLE2 from the COMP string.
	SAMPLE1="${COMP%%_V_*}"  # Remove everything after "_V_" to get the first sample name
	SAMPLE2="${COMP##*_V_}"  # Remove everything before and including "_V_" to get the second sample name
	### Intersect the "Detectable sites" for the two samples
	bedtools intersect \
		-a ./results/Detectabilities/Samples/$SAMPLE1'/AllChroms_CleanContiguousPairs.bed' \
		-b ./results/Detectabilities/Samples/$SAMPLE2'/AllChroms_CleanContiguousPairs.bed' \
		> ./results/Detectabilities/Comparisons/$COMP'/Intersection.bed'
	bedtools intersect \
		-a ./results/Detectabilities/Samples/$SAMPLE1'/AllChroms_AllPairs_CleanOnly.bed' \
		-b ./results/Detectabilities/Samples/$SAMPLE2'/AllChroms_AllPairs_CleanOnly.bed' \
		> ./results/Detectabilities/Comparisons/$COMP'/AllPairsIntersection.bed'
	### Count the proportion of the genome for which both samples could hypothetically be reliably phased
	python ./src/PropGetterDetectable.py ./results/Detectabilities/Comparisons/$COMP/Intersection.bed
	### Get SNPs for plotting histograms of coverage if wanted
	python ./src/GetUniqueSNPsFromPairsList.py ./results/Detectabilities/Comparisons/$COMP/AllPairsIntersection.bed
	wc ./results/Detectabilities/Comparisons/$COMP/AllPairsIntersection_UniqueSNPs_gRanges.txt
done

####### Now doing this for all pairwise comparisons
### This is an array of all the sample names.
samples=("FEB_T506" "FEB_T508" "NOV_T504" "FEB_T511" "FEB_T513" "NOV_T501" "FEB_T515")
### These are the combinations that are already done and should be skipped.
### The outer loop iterates over all the sample names.
### 'i' is the index of the current sample in the samples array.
for i in "${!samples[@]}"; do
	### The inner loop starts from the next sample (i+1) and goes till the last sample.
	### This way, we ensure we only do one-way comparisons (A vs B but not B vs A).
	for j in $(seq $((i+1)) $((${#samples[@]}-1))); do
		### Fetch the sample names based on the current indices.
		SAMPLE1=${samples[$i]}
		SAMPLE2=${samples[$j]}
		mkdir -p "./results/Detectabilities/Comparisons/${SAMPLE1}_v_${SAMPLE2}"
		bedtools intersect \
			-a ./results/Detectabilities/Samples/$SAMPLE1'/AllChroms_CleanContiguousPairs.bed' \
			-b ./results/Detectabilities/Samples/$SAMPLE2'/AllChroms_CleanContiguousPairs.bed' \
			> ./results/Detectabilities/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/Intersection.bed'
		bedtools intersect \
			-a ./results/Detectabilities/Samples/$SAMPLE1'/AllChroms_AllPairs_CleanOnly.bed' \
			-b ./results/Detectabilities/Samples/$SAMPLE2'/AllChroms_AllPairs_CleanOnly.bed' \
			> ./results/Detectabilities/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/AllPairsIntersection.bed'
		python ./src/PropGetterDetectable.py ./results/Detectabilities/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/Intersection.bed'
		python ./src/GetUniqueSNPsFromPairsList.py ./results/Detectabilities/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}'/AllPairsIntersection.bed'
		wc ./results/Detectabilities/Comparisons/${SAMPLE1}'_v_'${SAMPLE2}/AllPairsIntersection_UniqueSNPs_gRanges.txt
	done
done
