#! /usr/bin/bash

### This script is designed to clean and filter VCF files for analyses of haplotypic recombination among individuals of unknown pedigree.
### For haploid males, haplotypic recombination can be directly inferred, while for diploid females, haplotypic recombination can be inferred from LOH

GATK=/home/kip/tools/gatk-4.2.0.0/gatk

mkdir ./data/

### First we will perform within-colony analyses of recombination

### We start by obtaining working files that contain just the samples we need for each analysis.
### C16
# $GATK SelectVariants \
# 	-V ../DataFiles/data/ShortReadGenomes.vcf \
# 	--sample-name C16B-3-1 \
# 	--sample-name C16B-4-2 \
# 	--sample-name SM65 \
# 	--sample-name SM76 \
# 	--sample-name SM83 \
# 	--sample-name SM87 \
# 	--sample-name W06 \
# 	--sample-name W07 \
# 	--select-type-to-include SNP \
# 	-O ./data/C16_SNPs.vcf
# 
for COL in C16
do
# $GATK VariantFiltration \
# 	-V ./data/$COL'_SNPs.vcf' \
# 	--filter-name "QD" \
# 	--filter-expression "QD < 2.0" \
# 	--filter-name "FS" \
# 	--filter-expression "FS > 55.0" \
# 	--filter-name "MQ" \
# 	--filter-expression "MQ < 55.0" \
# 	--filter-name "MQRankSum" \
# 	--filter-expression "MQRankSum < -8.5" \
# 	--filter-name "ReadPosRankSum" \
# 	--filter-expression "ReadPosRankSum < -6.0" \
# 	-O ./data/$COL'_SNPs_GATKfilters.vcf'
# $GATK SelectVariants \
# 	-V ./data/$COL'_SNPs_GATKfilters.vcf' \
# 	--exclude-filtered \
# 	-O ./data/$COL'_SNPs_GATKfilt_RmFilt.vcf'
# vcftools \
# 	--vcf ./data/$COL'_SNPs_GATKfilt_RmFilt.vcf' \
# 	--exclude-positions ../DataFiles/data/HaploidMaleHet.PosScreener  \
# 	--recode \
# 	--recode-INFO-all \
# 	--out ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet'
# $GATK VariantsToTable \
# 	-V ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table'
# $GATK VariantsToTable \
# 	-V ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF AD \
# 	-O ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.AD.table'
# python ../LossOfHeterozygosity_KnownPedigreePairs/CleanData_KeepPutativelyAncestralHet.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table'
# python ./CleanData_KeepPutativelyAncestralHet_ModForHaploids.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table'
# python ../LossOfHeterozygosity_KnownPedigreePairs/ScreenSNPsonADnDP.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet' ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.AD.table'
# python ../AlignmentAndVariantCallingWGS/Process_VCFtoPosScreener.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs'
# vcftools \
# 	--vcf ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf' \
# 	--positions ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.PosScreener' \
# 	--recode \
# 	--recode-INFO-all \
# 	--out ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs'
# done
# 
# 
# $GATK SelectVariants \
# 	-V ./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf \
# 	--sample-name C16B-3-1 \
# 	--sample-name C16B-4-2 \
# 	--sample-name W06 \
# 	--sample-name W07 \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf
# 
# $GATK SelectVariants \
# 	-V ./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf \
# 	--sample-name SM65 \
# 	--sample-name SM76 \
# 	--sample-name SM83 \
# 	--sample-name SM87 \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf
# 
# 
# for COL in C16
# do
# for TYPE in DipFem HapMale
# do
# $GATK VariantsToTable \
# 	-V ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_'$TYPE'.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_'$TYPE'.vcf.table'
# python ./CleanData_AddContigNames.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_'$TYPE'.vcf.table' ../UsefulFiles/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt
# done
# done

###################
#
#   Pick up here
#
###################

# for COL in C16
# do
# python ./src/LOH_GroupOfDiploidsHomRunFinder.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs' ./data/
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt' 1500
# python ./src/PhaseClean_RemoveCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt' 1500
# python ./src/PhaseClean_RemoveCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt' 5000
# ## create a file for HaplotypeBlock presentation
# python ./src/PhaseClean_RemoveFirstBlockForEachContig.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt'
# done
# 
# 
# for COL in  C16 STC6_3ind
# do
# python ./src/Recomb_HapPhaseSwitchScanner.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs'
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches.txt' 1500
# python ./src/PhaseClean_RemoveCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches.txt' 1500
# python ./src/PhaseClean_MergeAdjacents.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP.txt'
# python ./src/PhaseClean_CountPhaseSwitches.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents.txt'
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches.txt' 5000
# python ./src/PhaseClean_RemoveCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches.txt' 5000
# python ./src/PhaseClean_MergeAdjacents.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan5000BP.txt'
# python ./src/PhaseClean_CountPhaseSwitches.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan5000BP_MergeAdjacents.txt'
# ## create a file for HaplotypeBlock presentation
# python ./src/PhaseClean_RemoveFirstBlockForEachContig.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsGr8rThan1500BP.txt'
# python ./src/PhaseClean_RemoveFirstBlockForEachContig.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents.txt'
# cat ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved.txt' > ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt'
# echo '' >> ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt'
# tail -n +2 ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsGr8rThan1500BP_FirstBlocksOnContigRemoved.txt' >> ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt'
# done


###################################
#
#
#
#	Now multi-colony
#
#
#
###################################


# $GATK SelectVariants \
# 	-V /home/antqueen/genomics/experiments/analyses/KDL20211215_VarCall4AllIndieGenomesToDateGATK4.2/results/AllIndidualsAsOf_20211215_AllChroms_SNPs.vcf \
# 	--sample-name C16B-3-1 \
# 	--sample-name C16B-4-2 \
# 	--sample-name SM65 \
# 	--sample-name SM76 \
# 	--sample-name SM83 \
# 	--sample-name SM87 \
# 	--sample-name W06 \
# 	--sample-name W07 \
# 	--sample-name KL46 \
# 	--sample-name KL66 \
# 	--sample-name C18 \
# 	--sample-name C1A-1 \
# 	-O ./data/LineA_MultiCol_SNPs.vcf
# 
# $GATK SelectVariants \
# 	-V /home/antqueen/genomics/experiments/analyses/KDL20211215_VarCall4AllIndieGenomesToDateGATK4.2/results/AllIndidualsAsOf_20211215_AllChroms_SNPs.vcf \
# 	--sample-name SM03 \
# 	--sample-name SM74 \
# 	--sample-name SM75 \
# 	--sample-name STC6L \
# 	--sample-name W04 \
# 	--sample-name W05 \
# 	--sample-name KL70 \
# 	--sample-name SM56 \
# 	--sample-name STC1_W \
# 	--sample-name STC12_W \
# 	-O ./data/LineB_MultiCol_SNPs.vcf

# for COL in LineA_MultiCol LineB_MultiCol
# do
# $GATK VariantFiltration \
# 	-V ./data/$COL'_SNPs.vcf' \
# 	--filter-name "QD" \
# 	--filter-expression "QD < 2.0" \
# 	--filter-name "FS" \
# 	--filter-expression "FS > 55.0" \
# 	--filter-name "MQ" \
# 	--filter-expression "MQ < 55.0" \
# 	--filter-name "MQRankSum" \
# 	--filter-expression "MQRankSum < -8.5" \
# 	--filter-name "ReadPosRankSum" \
# 	--filter-expression "ReadPosRankSum < -6.0" \
# 	-O ./data/$COL'_SNPs_GATKfilters.vcf'
# $GATK SelectVariants \
# 	-V ./data/$COL'_SNPs_GATKfilters.vcf' \
# 	--exclude-filtered \
# 	-O ./data/$COL'_SNPs_GATKfilt_RmFilt.vcf'
# vcftools \
# 	--vcf ./data/$COL'_SNPs_GATKfilt_RmFilt.vcf' \
# 	--exclude-positions ./data/AllIndidualsAsOf_20211215_AllChroms_HapMalesOnly_HetOnly.table \
# 	--recode \
# 	--recode-INFO-all \
# 	--out ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet'
# $GATK VariantsToTable \
# 	-V ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table'
# python ./src/CleanData_KeepPutativelyAncestralHet.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table'
# python ./src/Process_VCFtoPosScreener.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet'
# vcftools \
# 	--vcf ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf' \
# 	--positions ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet.PosScreener' \
# 	--recode \
# 	--recode-INFO-all \
# 	--out ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet'
# $GATK VariantsToTable \
# 	-V ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf.table'
# $GATK VariantsToTable \
# 	-V ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF AD \
# 	-O ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf_AD.table'
# python ./src/ScreenSNPsonADnDP.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf.table' ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf_AD.table'
# python ./src/Process_VCFtoPosScreener.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf.table_ADDPscreenNucs'
# vcftools \
# 	--vcf ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf' \
# 	--positions ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode.vcf.table_ADDPscreenNucs.PosScreener' \
# 	--recode \
# 	--recode-INFO-all \
# 	--out ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs'
# done


# $GATK SelectVariants \
# 	-V ./data/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode.vcf \
# 	--sample-name C16B-3-1 \
# 	--sample-name C16B-4-2 \
# 	--sample-name W06 \
# 	--sample-name W07 \
# 	--sample-name C18 \
# 	--sample-name C1A-1 \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_DipFem.vcf
# 
# $GATK SelectVariants \
# 	-V ./data/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode.vcf \
# 	--sample-name STC6L \
# 	--sample-name W04 \
# 	--sample-name W05 \
# 	--sample-name STC1_W \
# 	--sample-name STC12_W \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_DipFem.vcf
# 
# $GATK SelectVariants \
# 	-V ./data/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode.vcf \
# 	--sample-name SM65 \
# 	--sample-name SM76 \
# 	--sample-name SM83 \
# 	--sample-name SM87 \
# 	--sample-name KL46 \
# 	--sample-name KL66 \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_HapMale.vcf
# 
# $GATK SelectVariants \
# 	-V ./data/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode.vcf \
# 	--sample-name SM03 \
# 	--sample-name SM74 \
# 	--sample-name SM75 \
# 	--sample-name KL70 \
# 	--sample-name SM56 \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_HapMale.vcf


# for COL in LineA_MultiCol LineB_MultiCol
# do
# for TYPE in DipFem HapMale
# do
# $GATK VariantsToTable \
# 	-V ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_'$TYPE'.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_'$TYPE'.vcf.table'
# python ./src/CleanData_AddContigNames.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_'$TYPE'.vcf.table' ./data/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt
# done
# done

# for COL in LineA_MultiCol LineB_MultiCol
# do
# python ./src/LOH_GroupOfDiploidsHomRunFinder.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs' ./data/
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs.txt' 1500
# done

# for COL in LineA_MultiCol LineB_MultiCol
# do
# python ./src/Recomb_HapPhaseSwitchScanner.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs'
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches.txt' 1500
# python ./src/PhaseClean_RemoveCOLengthPBs.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches.txt' 1500
# python ./src/PhaseClean_MergeAdjacents.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP.txt'
# python ./src/PhaseClean_RemoveFirstBlockForEachContig.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents.txt'
# python ./src/PhaseClean_RemoveFirstBlockForEachContig.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsGr8rThan1500BP.txt'
# python ./src/PhaseClean_CountPhaseSwitches.py ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved.txt'
# cat ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved.txt' > ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt'
# echo '' >> ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt'
# tail -n +2 ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsGr8rThan1500BP_FirstBlocksOnContigRemoved.txt' >> ./data/$COL'_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt'
# done
# 
# ####### Pairwise analyses
# # mkdir -p ./data/Pairwise/LineA/DipFemLOH/
# # mkdir -p ./data/Pairwise/LineA/HapMaleRec/
# # mkdir -p ./data/Pairwise/LineB/DipFemLOH/
# # mkdir -p ./data/Pairwise/LineB/HapMaleRec/
# # 
# for SAMPLE in SM76 SM83 SM87 KL46 KL66
# do
# GATK=/home/kip/tools/gatk-4.2.0.0/gatk
# $GATK SelectVariants \
# 	-V ./data/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode.vcf \
# 	--sample-name SM65 \
# 	--sample-name $SAMPLE \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/Pairwise/LineA/HapMaleRec/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE'.vcf'
# $GATK VariantsToTable \
# 	-V ./data/Pairwise/LineA/HapMaleRec/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE'.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/Pairwise/LineA/HapMaleRec/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE'.vcf.table'
# python ./src/CleanData_AddContigNames.py ./data/Pairwise/LineA/HapMaleRec/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE.'vcf.table' ./data/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt
# python ./src/Recomb_HapPhaseSwitchScanner.py ./data/Pairwise/LineA/HapMaleRec/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE.'vcf.table_Contigs'
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/Pairwise/LineA/HapMaleRec/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches.txt' 1500
# python ./src/PhaseClean_MergeAdjacents_ForPairwise.py ./data/Pairwise/LineA/HapMaleRec/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP.txt'
# python ./src/PhaseClean_RemoveFirstBlockForEachContig.py ./data/Pairwise/LineA/HapMaleRec/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents.txt'
# python ./src/PhaseClean_CountPhaseSwitches.py ./data/Pairwise/LineA/HapMaleRec/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM65_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved.txt'
# done
# 
# 
# for SAMPLE in SM74 SM75 SM56 KL70
# do
# GATK=/home/kip/tools/gatk-4.2.0.0/gatk
# $GATK SelectVariants \
# 	-V ./data/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode.vcf \
# 	--sample-name SM03 \
# 	--sample-name $SAMPLE \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/Pairwise/LineB/HapMaleRec/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE'.vcf'
# $GATK VariantsToTable \
# 	-V ./data/Pairwise/LineB/HapMaleRec/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE'.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/Pairwise/LineB/HapMaleRec/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE'.vcf.table'
# python ./src/CleanData_AddContigNames.py ./data/Pairwise/LineB/HapMaleRec/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE.'vcf.table' ./data/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt
# python ./src/Recomb_HapPhaseSwitchScanner.py ./data/Pairwise/LineB/HapMaleRec/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE.'vcf.table_Contigs'
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/Pairwise/LineB/HapMaleRec/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches.txt' 1500
# python ./src/PhaseClean_MergeAdjacents_ForPairwise.py ./data/Pairwise/LineB/HapMaleRec/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP.txt'
# python ./src/PhaseClean_RemoveFirstBlockForEachContig.py ./data/Pairwise/LineB/HapMaleRec/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents.txt'
# python ./src/PhaseClean_CountPhaseSwitches.py ./data/Pairwise/LineB/HapMaleRec/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_SM03_v_$SAMPLE.'vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved.txt'
# done
# 
# 
# for SAMPLE in W07 C16B-3-1 C16B-4-2 C18 C1A-1
# do
# GATK=/home/kip/tools/gatk-4.2.0.0/gatk
# $GATK SelectVariants \
# 	-V ./data/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode.vcf \
# 	--sample-name W06 \
# 	--sample-name $SAMPLE \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/Pairwise/LineA/DipFemLOH/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W06_v_$SAMPLE'.vcf'
# $GATK VariantsToTable \
# 	-V ./data/Pairwise/LineA/DipFemLOH/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W06_v_$SAMPLE'.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/Pairwise/LineA/DipFemLOH/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W06_v_$SAMPLE'.vcf.table'
# python ./src/CleanData_AddContigNames.py ./data/Pairwise/LineA/DipFemLOH/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W06_v_$SAMPLE.'vcf.table' ./data/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt
# python ./src/LOH_GroupOfDiploidsHomRunFinder.py ./data/Pairwise/LineA/DipFemLOH/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W06_v_$SAMPLE.'vcf.table_Contigs' ./data/Pairwise/LineA/DipFemLOH/
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/Pairwise/LineA/DipFemLOH/LineA_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W06_v_$SAMPLE.'vcf.table_Contigs_UniqLOHPresAbs.txt' 1500
# done
# 
# for SAMPLE in W05 STC6L STC1_W STC12_W
# do
# GATK=/home/kip/tools/gatk-4.2.0.0/gatk
# $GATK SelectVariants \
# 	-V ./data/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode.vcf \
# 	--sample-name W04 \
# 	--sample-name $SAMPLE \
# 	--exclude-non-variants \
# 	--remove-unused-alternates \
# 	-O ./data/Pairwise/LineB/DipFemLOH/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W04_v_$SAMPLE'.vcf'
# $GATK VariantsToTable \
# 	-V ./data/Pairwise/LineB/DipFemLOH/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W04_v_$SAMPLE'.vcf' \
# 	-F CHROM -F POS -F REF -F ALT -F FILTER \
# 	-GF GT \
# 	-O ./data/Pairwise/LineB/DipFemLOH/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W04_v_$SAMPLE'.vcf.table'
# python ./src/CleanData_AddContigNames.py ./data/Pairwise/LineB/DipFemLOH/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W04_v_$SAMPLE.'vcf.table' ./data/ContigPos_for_Obir.assembly.v5.4_ChromsOnly.txt
# python ./src/LOH_GroupOfDiploidsHomRunFinder.py ./data/Pairwise/LineB/DipFemLOH/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W04_v_$SAMPLE.'vcf.table_Contigs' ./data/Pairwise/LineB/DipFemLOH/
# python ./src/PhaseClean_RemoveNCOLengthPBs.py ./data/Pairwise/LineB/DipFemLOH/LineB_MultiCol_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet.recode_ADDPscreenNucs.recode_W04_v_$SAMPLE.'vcf.table_Contigs_UniqLOHPresAbs.txt' 1500
# done
# 
