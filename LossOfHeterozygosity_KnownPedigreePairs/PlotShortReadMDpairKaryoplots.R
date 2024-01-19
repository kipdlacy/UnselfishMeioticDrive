setwd("/Users/kiplacy/../../../../Volumes/antqueen/booster/KDL20240109_UploadCodeParthPaper/UnselfishMeioticDrive/LossOfHeterozygosity_KnownPedigreePairs/")

### Load plotting package
library(karyoploteR)

### Load O. biroi genome parameters for KaryoploteR
Obir5.4withChr<-toGRanges("../UsefulFiles/GenomeForKaryoplotR_Obir5_4_withChr.txt")

### Load data
# Heterozygous sites in the mother of pair 1 --- note that these are NOT manually curated
Pair1Mhet<-toGRanges('./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode.vcf.table_MD-B-1-MHet.PosScreener')
# Heterozygous sites in the daughter of pair 1 --- note that these are NOT manually curated
Pair1Dhet<-toGRanges('./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode.vcf.table_MD-B-1-DHet.PosScreener')
# Manually-curated sites that differ in heterozygosity between mother and daughter of pair 1 -- there was only a single SNP that was heterozygous in the mother but not in the daughter
Pair1LOH<-toGRanges('./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode_MD-B-1-M_v_MD-B-1-D.vcf.table_Diffs_Curated.PosScreener')
# Heterozygous sites in the mother of pair 2 --- note that these are NOT manually curated
Pair2Mhet<-toGRanges('./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode.vcf.table_MD-B-7-MHet.PosScreener')
# Heterozygous sites in the daughter of pair 2 --- note that these are NOT manually curated
Pair2Dhet<-toGRanges('./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode.vcf.table_MD-B-7-D1Het.PosScreener')
# Manually-curated sites that differ in heterozygosity between mother and daughter of pair 2 -- there was only a single SNP that was homozygous in the mother but heterozygous in the daughter
Pair2GOH<-toGRanges('./data/ShortRead_MDpairs_SNPs_GATKfilt_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode_MD-B-7-M_v_MD-B-7-D1.vcf.table_Diffs_Curated.PosScreener')

### Plot the data along the reference genome
pdf("./results/Karyoplot_GOHLOH_Borders_Dhet2.pdf",width = 15,height = 7) ## <- extra height to accomodate
kp<-plotKaryotype(genome = Obir5.4withChr, plot.type = 4)
at <- autotrack(current.track = 6,total.tracks = 6)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1LOH,window.size = 50000,col='magenta',border='magenta')
kpAddLabels(kp, labels="1", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 5,total.tracks = 6)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1Mhet,window.size = 250000,col='grey',border=NA)
kpAxis(kp, ymax=332, r0=at$r0, r1=at$r1, cex=0.8)
kpAddLabels(kp, labels="1Mhet", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 4,total.tracks = 6)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1Dhet,window.size = 250000,col='grey',border=NA)
kpAxis(kp, ymax=332, r0=at$r0, r1=at$r1, cex=0.8)
kpAddLabels(kp, labels="1Dhet", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 3,total.tracks = 6)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2GOH,window.size = 50000,col='darkgreen',border='darkgreen')
kpAddLabels(kp, labels="2", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 2,total.tracks = 6)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2Mhet,window.size = 250000,col='grey',border=NA)
kpAxis(kp, ymax=332, r0=at$r0, r1=at$r1, cex=0.8)
kpAddLabels(kp, labels="2Mhet", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 1,total.tracks = 6)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2Dhet,window.size = 250000,col='grey',border=NA)
kpAxis(kp, ymax=332, r0=at$r0, r1=at$r1, cex=0.8)
kpAddLabels(kp, labels="2Dhet", r0=at$r0, r1=at$r1, data.panel = 1)
dev.off()

### Get axes for density plots
# For pair 1 mother
kp<-kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1Mhet,window.size = 250000,col='grey',border=NA)
kp$latest.plot$computed.values$max.density # result: 332
# For pair 1 daughter
kp<-kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1Dhet,window.size = 250000,col='grey',border=NA)
kp$latest.plot$computed.values$max.density # result: 332
# For pair 2 mother
kp<-kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2Mhet,window.size = 250000,col='grey',border=NA)
kp$latest.plot$computed.values$max.density # result: 332
# For pair 2 daughter
kp<-kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2Dhet,window.size = 250000,col='grey',border=NA)
kp$latest.plot$computed.values$max.density # result: 332


