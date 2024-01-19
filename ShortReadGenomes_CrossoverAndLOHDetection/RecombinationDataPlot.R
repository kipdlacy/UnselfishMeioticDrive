setwd("/Users/kiplacy/../../../../Volumes/antqueen/booster/KDL20240109_UploadCodeParthPaper/UnselfishMeioticDrive/ShortReadGenomes_CrossoverAndLOHDetection/")
getwd()

### This script is designed to plot recombination events observed among haploid male genomes, 
### and inferred among diploid female genomes using segmental losses of heterozygosity as a proxy.

###########################################
#                                         #
#                                         #
#     Plotting crossovers genome-wide     #
#                                         #
#                                         #
###########################################

### Load plotting packages
library(ggplot2)
library(karyoploteR)

### Import genome-specific files for KaryoploteR
Obir5.4withChr<-toGRanges("../UsefulFiles/GenomeForKaryoplotR_Obir5_4_withChr.txt")
Obir5.4withChr
ObirTigswithChr<-toGRanges("../UsefulFiles/ContigsForKaryoplotR_Obir5_4_withChr.txt")
ObirTigswithChr

### Import data files
# Crossovers inferred among C16 diploids using LOH as a proxy
C16_LOH_CO<-toGRanges('./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs_NoPBsLessThan1500BP_Curated_COonly.txt')
# Crossovers observed among C16 haploids
C16_HapRec<-toGRanges('./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_PSforGRanges_ManuallyScreened.txt')
# Putatively ancestrally heterozygous sites in colony C16
C16_PutAncHet<-toGRanges('./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.PosScreener')
# Crossovers inferred among STC6 diploids using LOH as a proxy
STC6_LOH_CO<-toGRanges('./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs_NoPBsLessThan1500BP_Curated_COonly.txt')
# Crossovers observed among STC6 haploids
STC6_HapRec<-toGRanges('./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_PSforGRanges_ManuallyScreened.txt')
# Putatively ancestrally heterozygous sites in colony STC6
STC6_PutAncHet<-toGRanges('./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.PosScreener')

### Plot this data genome-wide using Karyoploter
pdf("./results/Karyoplot_BothLines_UnknownPedigree_Border_LOHCOs.pdf",width = 15,height = 5)
kp<-plotKaryotype(genome = Obir5.4withChr, cytobands = ObirTigswithChr,plot.type = 4)
at <- autotrack(current.track = 6,total.tracks = 6)
kpPlotRegions(kp, r0=at$r0, r1=at$r1, data=C16_LOH_CO,col='black',border='black')
kpAddLabels(kp, labels="LOH among Diploid Females", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 5,total.tracks = 6)
kpPlotRegions(kp, r0=at$r0, r1=at$r1, data=C16_HapRec,col='black',border='black')
kpAddLabels(kp, labels="CO among Haploid Males", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 4,total.tracks = 6)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=C16_PutAncHet,window.size = 200000,col='grey',border=NA)
kpAxis(kp, ymax=65, r0=at$r0, r1=at$r1, cex=0.8)
at <- autotrack(current.track = 3,total.tracks = 6)
kpPlotRegions(kp, r0=at$r0, r1=at$r1, data=STC6_LOH_CO,col='black',border='black')
kpAddLabels(kp, labels="LOH among Diploid Females", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 2,total.tracks = 6)
kpPlotRegions(kp, r0=at$r0, r1=at$r1, data=STC6_HapRec,col='black',border='black')
kpAddLabels(kp, labels="CO among Haploid Males", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 1,total.tracks = 6)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=STC6_PutAncHet,window.size = 200000,col='grey',border=NA)
kpAxis(kp, ymax=151, r0=at$r0, r1=at$r1, cex=0.8)
dev.off()

### Get axes for density plots
kp<-kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=C16_PutAncHet,window.size = 200000,col='grey',border=NA)
kp$latest.plot$computed.values$max.density ## result: 65
kp<-kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=STC6_PutAncHet,window.size = 200000,col='grey',border=NA)
kp$latest.plot$computed.values$max.density ## result: 151

################################################
#                                              #
#                                              #
#     Haplotype block length distributions     #
#                                              #
#                                              #
################################################

library(ggplot2)
library(scales)
options(scipen=0)

### Load the haplotype blocks for recombination events observed among C16 haploids
C16HapRecHaplotypes<-read.table("./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions.txt",header = T)
C16blockLengths<-as.data.frame(C16HapRecHaplotypes$InfBlockLength)
C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`<-as.numeric(C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`)
### Load the LOH blocks for recombination events inferred via LOH among C16 diploids
C16DipFemHaplotypes<-read.table("./data/C16_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs_ManuallyScreened.txt",header = T)
C16DipFemblockLengths<-as.data.frame(C16DipFemHaplotypes$InfLOHLength)
C16DipFemblockLengths$`C16DipFemHaplotypes$InfLOHLength`<-as.numeric(C16DipFemblockLengths$`C16DipFemHaplotypes$InfLOHLength`)

### Plot haplotype block length histogram
pdf("./results/BlockLengths_C16haprec_NoLogX.pdf",width = 5,height = 4)
PICTURE <-ggplot(C16blockLengths, aes(x=C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 50) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",limits = c(-100000,4000000)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()
### C16 recombination among haploids--Plot haplotype block length vs SNP number scatterplot to show that long haplotype blocks are not supported by just a few SNPs
pdf("./results/NumSNPsVBlockLengths_C16haprec_NoLogX_inclGeneConv.pdf",width = 5,height = 4)
SCAT <-ggplot(C16HapRecHaplotypes, aes(x=C16HapRecHaplotypes$InfBlockLength,y=C16HapRecHaplotypes$NumSNPsInBlock)) +
  geom_point(alpha=0.5) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(0,1000000,2000000,3000000,4000000)) +
  scale_y_continuous(name="Number of SNPs in Haplotype Block") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()
### Plot LOH block length histogram
pdf("./results/BlockLengths_C16dipfem_NoLogX.pdf",width = 5,height = 4)
PICTURE <-ggplot(C16DipFemblockLengths, aes(x=C16DipFemblockLengths$`C16DipFemHaplotypes$InfLOHLength`)) + 
  geom_histogram(bins = 30) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)", breaks = c(0,500000,1000000,1500000,2000000,2500000)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()
### C16 LOH among diploids -- scatterplot of block length vs SNP number to show that long LOH tracts are supported by a corresponding number of SNPs
pdf("./results/NumSNPsVBlockLengths_C16dipFem_NoLogX_inclGeneConv.pdf",width = 5,height = 4)
SCAT <-ggplot(C16DipFemHaplotypes, aes(x=C16DipFemHaplotypes$InfLOHLength,y=C16DipFemHaplotypes$NumSNPsInLOH)) +
  geom_point(alpha=0.5, size=3) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(0,500000,1000000,1500000,2000000,2500000)) +
  scale_y_continuous(name="Number of SNPs in Haplotype Block") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()

### Load the haplotype blocks for recombination events observed among STC6 haploids
STC6HapRecHaplotypes<-read.table("./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents_FirstBlocksOnContigRemoved_MergedWithGeneConversions_ManuallyScreened.txt",header = T)
STC6blockLengths<-as.data.frame(STC6HapRecHaplotypes$InfBlockLength)
STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`<-as.numeric(STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`)
### Load the LOH blocks for recombination events inferred via LOH among STC6 diploids
STC6DipFemHaplotypes<-read.table("./data/STC6_SNPs_GATKfilt_RmFilt_NoHapHet.recode_PutAncHet_ADDPscreenNucs.recode_DipFem.vcf.table_Contigs_UniqLOHPresAbs_ManuallyScreened.txt",header = T)
STC6DipFemblockLengths<-as.data.frame(STC6DipFemHaplotypes$InfLOHLength)
STC6DipFemblockLengths$`STC6DipFemHaplotypes$InfLOHLength`<-as.numeric(STC6DipFemblockLengths$`STC6DipFemHaplotypes$InfLOHLength`)

### STC6 haploids -- Plot haplotype block length histogram
pdf("./results/BlockLengths_STC6haprec_NoLogX.pdf",width = 5,height = 4)
PICTURE <-ggplot(STC6blockLengths, aes(x=STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 50) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",limits = c(-100000,6000000),breaks = c(0,1000000,2000000,3000000,4000000,5000000,6000000)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()
### STC6 recombination among haploids--Plot haplotype block length vs SNP number scatterplot to show that long haplotype blocks are not supported by just a few SNPs
pdf("./results/NumSNPsVBlockLengths_STC6haprec_NoLogX_inclGeneConv.pdf",width = 5,height = 4)
SCAT <-ggplot(STC6HapRecHaplotypes, aes(x=STC6HapRecHaplotypes$InfBlockLength,y=STC6HapRecHaplotypes$NumSNPsInBlock)) +
  geom_point(alpha=0.5) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(0,1000000,2000000,3000000,4000000,5000000,6000000)) +
  scale_y_continuous(name="Number of SNPs in Haplotype Block") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()
### STC6 diploids -- Plot LOH block length histogram
pdf("./results/BlockLengths_STC6dipfem_NoLogX.pdf",width = 5,height = 4)
PICTURE <-ggplot(STC6DipFemblockLengths, aes(x=STC6DipFemblockLengths$`STC6DipFemHaplotypes$InfLOHLength`)) + 
  geom_histogram(bins = 30) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(0,500000,1000000,1500000), limits = c(-50000,1500000)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()
### STC6 LOH among diploids -- scatterplot of block length vs SNP number to show that long LOH tracts are supported by a corresponding number of SNPs
pdf("./results/NumSNPsVBlockLengths_STC6dipFem_NoLogX_inclGeneConv.pdf",width = 5,height = 4)
SCAT <-ggplot(STC6DipFemHaplotypes, aes(x=STC6DipFemHaplotypes$InfLOHLength,y=STC6DipFemHaplotypes$NumSNPsInLOH)) +
  geom_point(alpha=0.5, size=3) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(0,500000,1000000,1500000,2000000), limits=c(-50000,1500000)) +
  scale_y_continuous(name="Number of SNPs in Haplotype Block") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()




#####
#####
#  End Official, Prelim Below
#####
#####


### first for crossovers only
C16HapRecHaplotypes<-read.table("./data/C16_SNPs_GATKfilters_RmFilt_NoHapHet.recode_PutAncHet.recode_HapMale.vcf.table_Nucs_ADDPscreen_Intersected.recode.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1500BP_MergeAdjacents.txt",header = T)
head(C16HapRecHaplotypes)
C16blockLengths<-as.data.frame(C16HapRecHaplotypes$InfBlockLength)
head(C16blockLengths)
C16blockLengths(0,1)
str(C16blockLengths)
C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`<-as.numeric(C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`)
str(C16blockLengths)
head(C16blockLengths)


pdf("./results/UnknownPedigree/BlockLengths_C16haprec_LogX.pdf",width = 6,height = 5)
PICTURE <-ggplot(C16blockLengths, aes(x=C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 20) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",trans = 'log10',breaks = c(1,10,100,1000,10000,100000,1000000),limits = c(999,4000000)) +
  #scale_y_continuous(name="Counts",trans='log10',breaks = c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,150),limits = c(1,150)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()

pdf("./results/UnknownPedigree/BlockLengths_C16haprec_NoLogX.pdf",width = 6,height = 5)
PICTURE <-ggplot(C16blockLengths, aes(x=C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 50) +
  #geom_histogram(binwidth = 1000) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",limits = c(1,4000000)) +
  #scale_y_continuous(name="Counts",trans='log10',breaks = c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,150),limits = c(1,150)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()

### THEN INCLUDING GENE CONVERSIONS
C16HapRecHaplotypes<-read.table("./data/C16_SNPs_GATKfilters_RmFilt_NoHapHet.recode_PutAncHet.recode.vcf_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches.txt",header = T)
head(C16HapRecHaplotypes)
C16blockLengths<-as.data.frame(C16HapRecHaplotypes$InfBlockLength)
head(C16blockLengths)
C16blockLengths(0,1)
str(C16blockLengths)
C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`<-as.numeric(C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`)
str(C16blockLengths)
head(C16blockLengths)

pdf("./results/UnknownPedigree/BlockLengths_C16haprec_LogX_inclGeneConv.pdf",width = 6,height = 5)
PICTURE <-ggplot(C16blockLengths, aes(x=C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 50) +
  #geom_histogram(binwidth = 1000) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)", trans='log10', breaks = c(1,10,100,1000,10000,100000,1000000,10000000)) +
  #scale_y_continuous(name="Counts",trans='log10',breaks = c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,150),limits = c(1,150)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()

pdf("./results/UnknownPedigree/BlockLengths_C16haprec_NoLogX_inclGeneConv.pdf",width = 6,height = 5)
PICTURE <-ggplot(C16blockLengths, aes(x=C16blockLengths$`C16HapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 50) +
  #geom_histogram(binwidth = 1000) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(1000,1000000,2000000,3000000,4000000)) +
  #scale_y_continuous(name="Counts",trans='log10',breaks = c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,150),limits = c(1,150)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()

pdf("./results/UnknownPedigree/NumSNPsVBlockLengths_C16haprec_NoLogX_inclGeneConv.pdf",width = 6,height = 5)
SCAT <-ggplot(C16HapRecHaplotypes, aes(x=C16HapRecHaplotypes$InfBlockLength,y=C16HapRecHaplotypes$NumSNPsInBlock)) +
  geom_point(alpha=0.5) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(1000,1000000,2000000,3000000,4000000)) +
  #scale_x_continuous(name="Length of Haplotype Block (Bp)") +
  scale_y_continuous(name="Number of SNPs in Haplotype Block") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()


#### first for crossovers only 
STC6HapRecHaplotypes<-read.table("./data/STC6_3ind_SNPs_GATKfilters_RmFilt_NoHapHet.recode_PutAncHet.recode.vcf_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches_NoPBsLessThan1000BP_MergeAdjacents.txt",header = T)
head(STC6HapRecHaplotypes)
STC6blockLengths<-as.data.frame(STC6HapRecHaplotypes$InfBlockLength)
head(STC6blockLengths)
str(STC6blockLengths)
STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`<-as.numeric(STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`)
str(STC6blockLengths)
head(STC6blockLengths)


pdf("./results/UnknownPedigree/BlockLengths_STC6haprec_LogX.pdf",width = 6,height = 5)
PICTURE <-ggplot(STC6blockLengths, aes(x=STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 20) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",trans = 'log10',breaks = c(1,10,100,1000,10000,100000,1000000),limits = c(999,4000000)) +
  #scale_y_continuous(name="Counts",trans='log10',breaks = c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,150),limits = c(1,150)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()

pdf("./results/UnknownPedigree/BlockLengths_STC6haprec_NoLogX.pdf",width = 6,height = 5)
PICTURE <-ggplot(STC6blockLengths, aes(x=STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 20) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)") +
  #scale_y_continuous(name="Counts",trans='log10',breaks = c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,150),limits = c(1,150)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()


### THEN INCLUDING GENE CONVERSIONS
STC6HapRecHaplotypes<-read.table("./data/STC6_3ind_SNPs_GATKfilters_RmFilt_NoHapHet.recode_PutAncHet.recode.vcf_ADDPscreenNucs.recode_HapMale.vcf.table_Contigs_PhaseSwitches.txt",header = T)
head(STC6HapRecHaplotypes)
STC6blockLengths<-as.data.frame(STC6HapRecHaplotypes$InfBlockLength)
head(STC6blockLengths)
str(STC6blockLengths)
STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`<-as.numeric(STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`)
str(STC6blockLengths)
head(STC6blockLengths)

pdf("./results/UnknownPedigree/BlockLengths_STC6haprec_NoLogX_inclGeneConv.pdf",width = 6,height = 5)
PICTURE <-ggplot(STC6blockLengths, aes(x=STC6blockLengths$`STC6HapRecHaplotypes$InfBlockLength`)) + 
  geom_histogram(bins = 50) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(1000,1000000,2000000,3000000,4000000,5000000,6000000)) +
  #scale_y_continuous(name="Counts",trans='log10',breaks = c(1,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,150),limits = c(1,150)) +
  scale_y_continuous(name="Number of Blocks") +
  theme_bw(base_size = 15)
print(PICTURE)
dev.off()

pdf("./results/UnknownPedigree/NumSNPsVBlockLengths_STC6haprec_NoLogX_inclGeneConv.pdf",width = 6,height = 5)
SCAT <-ggplot(STC6HapRecHaplotypes, aes(x=STC6HapRecHaplotypes$InfBlockLength,y=STC6HapRecHaplotypes$NumSNPsInBlock)) +
  geom_point(alpha=0.5) +
  scale_x_continuous(name="Length of Haplotype Block (Bp)",breaks = c(1000,1000000,2000000,3000000,4000000,5000000,6000000)) +
  #scale_x_continuous(name="Length of Haplotype Block (Bp)") +
  scale_y_continuous(name="Number of SNPs in Haplotype Block") +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()
