setwd("/Users/kiplacy/../../../../Volumes/antqueen/genomics/experiments/analyses/KDL20220509_FormalLOHnGOHAnalyses/")
getwd()

library(karyoploteR)
Obir5.4withChr<-toGRanges("../KDL20220407_TELLSeq_PhaseSwitchAutomatedDetectNScreen/data/KaryoplotR/TestObir5_4_withChr.txt")
Obir5.4withChr
ObirTigswithChr<-toGRanges("../KDL20220407_TELLSeq_PhaseSwitchAutomatedDetectNScreen/data/KaryoplotR/ContigFile_withChr.txt")
ObirTigswithChr

Pair1Mhet<-toGRanges('./results/NonTELL_MDpairs_SNPs_GATKfilters_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode_Pair1.vcf.table_MD-B-1-MHet.PosScreener')
Pair1Mhet
Pair1Dhet<-toGRanges('./results/NonTELL_MDpairs_SNPs_GATKfilters_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode_Pair1.vcf.table_MD-B-1-DHet.PosScreener')
Pair1Dhet
Pair1LOH<-toGRanges('./results/NonTELL_MDpairs_SNPs_GATKfilters_RmFilt_NoHapHet.recode_Pair1.vcf.table_Diffs_AllNucs_HetMADatLeast0.25_HomMADequalTo0_NoCovLessThn15_Curated.PosScreener')
Pair1LOH
Pair2Mhet<-toGRanges('./results/NonTELL_MDpairs_SNPs_GATKfilters_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode_Pair7.vcf.table_MD-B-7-MHet.PosScreener')
Pair2Mhet
Pair2Dhet<-toGRanges('./results/NonTELL_MDpairs_SNPs_GATKfilters_RmFilt_NoHapHet.recode.vcf.table_PutAncHet_ADDPscreenNucs.recode_Pair7.vcf.table_MD-B-7-D1Het.PosScreener')
Pair2Dhet
Pair2GOH<-toGRanges('./results/NonTELL_MDpairs_SNPs_GATKfilters_RmFilt_NoHapHet.recode_Pair7.vcf.table_Diffs_AllNucs_HetMADatLeast0.25_HomMADequalTo0_NoCovLessThn15_Curated.PosScreener')
Pair2GOH
MumCOH<-toGRanges('./results/NonTELL_MDpairs_SNPs_GATKfilters_RmFilt_NoHapHet.recode_MumsOnlyUnknownPedigree.vcf.table_Diffs_AllNucs_HetMADatLeast0.25_HomMADequalTo0_NoCovLessThn15_Screened.PosScreener_Merged')
MumCOH
MumLOH<-toGRanges('./results/NonTELL_MDpairs_SNPs_GATKfilters_RmFilt_NoHapHet.recode_MumsOnlyUnknownPedigree.vcf.table_Diffs_AllNucs_HetMADatLeast0.25_HomMADequalTo0_NoCovLessThn15_Screened.PosScreener_Merged_LOH')
MumLOH
MumGOH<-toGRanges('./results/NonTELL_MDpairs_SNPs_GATKfilters_RmFilt_NoHapHet.recode_MumsOnlyUnknownPedigree.vcf.table_Diffs_AllNucs_HetMADatLeast0.25_HomMADequalTo0_NoCovLessThn15_Screened.PosScreener_Merged_GOH')
MumGOH

?pdf()

#########New plot version with all daughter het in density plots too
pdf("./results/Plots/Karyoplot_GOHLOH_Borders_Dhet2.pdf",width = 15,height = 7) ## <- extra height to accomodate
pdf("/Users/kiplacy/Desktop/KaryoplotNonTELL_GOHLOH_Borders_Dhet2.pdf",width = 18,height = 7)
pdf("/Users/kiplacy/Desktop/KaryoplotJustIdeogram.pdf",width = 15,height = 7)
#pdf("/Users/kiplacy/Desktop/KaryoplotNonTELL_GOHLOH_Borders_Dhet2_Smaller.pdf",width = 7.08,height = 2.75)
#kp<-plotKaryotype(genome = Obir5.4withChr, cytobands = ObirTigswithChr,plot.type = 4)
kp<-plotKaryotype(genome = Obir5.4withChr, plot.type = 4)
#kpAddMainTitle(kp, main = 'Crossovers')
at <- autotrack(current.track = 7,total.tracks = 7)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1LOH,window.size = 50000,col='magenta',border='magenta')
kpAddLabels(kp, labels="1", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 6,total.tracks = 7)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1Mhet,window.size = 250000,col='grey',border=NA)
kpAxis(kp, ymax=332, r0=at$r0, r1=at$r1, cex=0.8)
kpAddLabels(kp, labels="1Mhet", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 5,total.tracks = 7)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1Dhet,window.size = 250000,col='grey',border=NA)
kpAxis(kp, ymax=332, r0=at$r0, r1=at$r1, cex=0.8)
kpAddLabels(kp, labels="1Dhet", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 4,total.tracks = 7)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2GOH,window.size = 50000,col='darkgreen',border='darkgreen')
kpAddLabels(kp, labels="2", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 3,total.tracks = 7)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2Mhet,window.size = 250000,col='grey',border=NA)
kpAxis(kp, ymax=332, r0=at$r0, r1=at$r1, cex=0.8)
kpAddLabels(kp, labels="2Mhet", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 2,total.tracks = 7)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2Dhet,window.size = 250000,col='grey',border=NA)
kpAxis(kp, ymax=332, r0=at$r0, r1=at$r1, cex=0.8)
kpAddLabels(kp, labels="2Dhet", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 1,total.tracks = 7)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=MumLOH,window.size = 50000,col='magenta',border='magenta')
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=MumGOH,window.size = 50000,col='darkgreen',border='darkgreen')
kpAddLabels(kp, labels="Mums", r0=at$r0, r1=at$r1, data.panel = 1)
dev.off()

print('Stop Kip You Went Too Far!')
print('Stop Kip You Went Too Far!')
print('Stop Kip You Went Too Far!')



pdf("./results/Plots/Karyoplot_GOHLOH_Borders.pdf",width = 15,height = 6)
pdf("/Users/kiplacy/Desktop/KaryoplotNonTELL_GOHLOH_Borders.pdf",width = 18,height = 6)
kp<-plotKaryotype(genome = Obir5.4withChr, cytobands = ObirTigswithChr,plot.type = 4)
#kpAddMainTitle(kp, main = 'Crossovers')
at <- autotrack(current.track = 5,total.tracks = 5)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1LOH,window.size = 80000,col='magenta',border='magenta')
kpAddLabels(kp, labels="1", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 4,total.tracks = 5)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1Mhet,window.size = 250000,col='grey',border=NA)
kpAxis(kp, ymax=332, r0=at$r0, r1=at$r1, cex=0.8)
kpAddLabels(kp, labels="1Mhet", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 3,total.tracks = 5)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2GOH,window.size = 80000,col='darkgreen',border='darkgreen')
kpAddLabels(kp, labels="2", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 2,total.tracks = 5)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2Mhet,window.size = 250000,col='grey',border=NA)
kpAxis(kp, ymax=332, r0=at$r0, r1=at$r1, cex=0.8)
kpAddLabels(kp, labels="2Mhet", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 1,total.tracks = 5)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=MumCOH,window.size = 2000,col='black',border='black')
kpAddLabels(kp, labels="Mums", r0=at$r0, r1=at$r1, data.panel = 1)
dev.off()

### Get axes for density plots
kp<-kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1Mhet,window.size = 250000,col='grey',border=NA)
kp$latest.plot$computed.values$max.density
## 332
kp<-kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2Mhet,window.size = 250000,col='grey',border=NA)
kp$latest.plot$computed.values$max.density
## 332

pdf("./results/Plots/Karyoplot_GOHLOH.pdf",width = 15,height = 6)
kp<-plotKaryotype(genome = Obir5.4withChr, cytobands = ObirTigswithChr,plot.type = 4)
#kpAddMainTitle(kp, main = 'Crossovers')
at <- autotrack(current.track = 5,total.tracks = 5)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1LOH,window.size = 80000,col='red',border=NA)
kpAddLabels(kp, labels="1", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 4,total.tracks = 5)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1Mhet,window.size = 120000,col='grey',border=NA)
kpAddLabels(kp, labels="1Mhet", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 3,total.tracks = 5)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2GOH,window.size = 80000,col='blue',border=NA)
kpAddLabels(kp, labels="2", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 2,total.tracks = 5)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2Mhet,window.size = 120000,col='grey',border=NA)
kpAddLabels(kp, labels="2Mhet", r0=at$r0, r1=at$r1, data.panel = 1)
at <- autotrack(current.track = 1,total.tracks = 5)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=MumCOH,window.size = 80000,col='black',border=NA)
kpAddLabels(kp, labels="Mums", r0=at$r0, r1=at$r1, data.panel = 1)
dev.off()


pdf("./results/Plots/Karyoplot_GOHLOH_ChromStack.pdf",width = 6,height = 10)
kp<-plotKaryotype(genome = Obir5.4withChr, cytobands = ObirTigswithChr,plot.type = 1)
at <- autotrack(current.track = 2,total.tracks = 2)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1LOH,window.size = 80000,col='red',border=NA)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair2GOH,window.size = 80000,col='blue',border=NA)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=MumCOH,window.size = 80000,col='black',border=NA)
at <- autotrack(current.track = 1,total.tracks = 2)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=Pair1Mhet,window.size = 200000,col='grey',border=NA)
dev.off()

