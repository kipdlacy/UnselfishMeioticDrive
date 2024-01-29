#! /usr/bin/bash

#### This script is designed to extract the haplotypic phase of each sample at all detected crossovers.
#### A file in Bed-like-format that contains all detected crossovers is read in as the "CrossoverFile"

### Pass path to crossover file
CrossoverFile="./data/AllCOsDetected.txt"
# New file path with modified name
FileToWrite="${CrossoverFile}_PhasedGTs.txt"
# Copy the file
cp "$CrossoverFile" "$FileToWrite"

GATK=/home/kip/tools/gatk-4.2.0.0/gatk

Samples=("FEB_T505" "FEB_T506" "FEB_T507" "FEB_T508" "FEB_T509" "FEB_T510" "FEB_T511" "FEB_T512" "FEB_T513" "FEB_T514" "FEB_T515" "FEB_T516" "NOV_T501" "NOV_T503" "NOV_T504" "NOV_T505" "NOV_T506")

for Sample in "${Samples[@]}"; do
	java -jar /home/sean/tools/picard/picard.jar GatherVcfs \
		INPUT=./data/Samples/$Sample/hapblock_chr1.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr2.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr3.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr4.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr5.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr6.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr7.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr8.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr9.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr10.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr11.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr12.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr13.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		INPUT=./data/Samples/$Sample/hapblock_chr14.phased_TEfiltered.recode_NoHapHet.recode_NonAbove2xMeanDepth.recode_SNPsOnly.vcf_NoGQLessThan99.recode_NoPQLessThan100.recode.vcf \
		OUTPUT=./data/Samples/$Sample/AllChroms_CleanPhasedGenotypes.vcf
	#### Then write to table
	$GATK VariantsToTable \
		-V ./data/Samples/$Sample/AllChroms_CleanPhasedGenotypes.vcf \
		-F CHROM -F POS -F REF -F ALT -F FILTER \
		-GF GT \
		-O ./data/Samples/$Sample/AllChroms_CleanPhasedGenotypes.vcf'.table'
	#### Then grab the phase at all crossovers
	python ./src/PhaseGrabber.py $FileToWrite ./data/Samples/$Sample/AllChroms_CleanPhasedGenotypes.vcf.table $Sample
	FileToWriteAndSample="${FileToWrite}_${Sample}"
	rm $FileToWrite
	mv $FileToWriteAndSample $FileToWrite
done
