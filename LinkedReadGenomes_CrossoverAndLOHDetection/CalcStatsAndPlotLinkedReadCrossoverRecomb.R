setwd('/Volumes/antqueen/genomics/experiments/analyses/KDL20220407_TELLSeq_PhaseSwitchAutomatedDetectNScreen/results/Plots/')
getwd()

### Load packages
library(ggplot2)
library(ggbeeswarm)
library(dplyr)
library(tidyr)

### Load Data
dat<-read.csv('./NumCONumMeiosisPropGenomeCoPhased.txt', header = T, sep = '\t')
dat
### Creating a dataframe with additional columns
dat <- dat %>%
  mutate(
    COperMeiosis = NumCO / NumMeiosis,
    AdjustedNumCO = NumCO / PropGenomesCoPhased,
    AdjustedCOperMeiosis = AdjustedNumCO / NumMeiosis,
    AdjustedCOperMeiosisAdjusted4Uninherited = (4 * AdjustedCOperMeiosis) / 3
  )
### Print the updated dataframe
print(dat)
### Calculating the mean for each column
mean_values <- dat %>%
  summarise(
    mean_PropGenomesCoPhased = mean(PropGenomesCoPhased, na.rm = TRUE),
    mean_NumMeiosis = mean(NumMeiosis, na.rm = TRUE),
    mean_NumCO = mean(NumCO, na.rm = TRUE),
    mean_COperMeiosis = mean(COperMeiosis, na.rm = TRUE),
    sd_COperMeiosis = sd(COperMeiosis, na.rm = TRUE),
    n = n(),
    sem_mean_COperMeiosis = sd_COperMeiosis / sqrt(n), # standard error
    lower95ci_COperMeiosis = mean_COperMeiosis - qt(0.975, n-1) * sem_mean_COperMeiosis,
    upper95ci_COperMeiosis = mean_COperMeiosis + qt(0.975, n-1) * sem_mean_COperMeiosis,
    mean_AdjustedNumCO = mean(AdjustedNumCO, na.rm = TRUE),
    mean_AdjustedCOperMeiosis = mean(AdjustedCOperMeiosis, na.rm = TRUE),
    mean_AdjustedCOperMeiosisAdjusted4Uninherited = mean(AdjustedCOperMeiosisAdjusted4Uninherited, na.rm = TRUE)
  )
### Print the summary table
print(mean_values)

### Calculating the mean and 95% confidence interval for AdjustedCOperMeiosisAdjusted4Uninherited
summary_stats <- dat %>%
  summarise(
    mean_AdjustedCOperMeiosisAdjusted4Uninherited = mean(AdjustedCOperMeiosisAdjusted4Uninherited, na.rm = TRUE),
    sd_AdjustedCOperMeiosisAdjusted4Uninherited = sd(AdjustedCOperMeiosisAdjusted4Uninherited, na.rm = TRUE),
    n = n(),
    sem = sd_AdjustedCOperMeiosisAdjusted4Uninherited / sqrt(n), # standard error
    ci_lower = mean_AdjustedCOperMeiosisAdjusted4Uninherited - qt(0.975, n-1) * sem,
    ci_lowerPerChrom = ci_lower / 14,
    ci_upper = mean_AdjustedCOperMeiosisAdjusted4Uninherited + qt(0.975, n-1) * sem,
    ci_upperPerChrom = ci_upper / 14
  )
### Print the summary stats, including the 95% confidence interval
print(summary_stats)

### Plotting the phase block lengths
BlockLengths<-read.csv('../Detectabilities/SampleScatterPlot_FEBT505_AllChroms_CleanContiguousPairs.txt', header = T, sep = '\t')
head(BlockLengths)
pdf("./SCAT_numSNPsvsBlockLength_FEB_T505.pdf",width = 5,height = 5)
SCAT <-ggplot(BlockLengths, aes(x=NumBp,y=NumSNPs)) +
  geom_point(alpha=0.5) +
  scale_x_continuous(name="Length of Phased Block (Bp)") +
  scale_y_continuous(name="Number of SNPs in Phased Block") +
  geom_vline(xintercept=27332, linetype="dashed", color = "red",linewidth=1) +
  theme_bw(base_size = 15)
print(SCAT)
dev.off()

### Plotting proportion of genome phased
dat=read.table('../Detectabilities/PropGenomePhasedOrDetectableSummary.txt', header=T)
dat
summary_stats <- dat %>%
  filter(Type == "Pair") %>%
  summarise(
    mean = mean(PropGenomePhasedOrDetectable),
    median = median(PropGenomePhasedOrDetectable),
    min = min(PropGenomePhasedOrDetectable),
    max = max(PropGenomePhasedOrDetectable),
    IQR = IQR(PropGenomePhasedOrDetectable)
  )
print(summary_stats)
pdf("PropGenomePhasedOrDetectable.pdf",width = 4,height = 5)
p<-ggplot(data=dat, aes(x=Type, y=PropGenomePhasedOrDetectable)) +
  geom_beeswarm(shape=21,fill="black",color="black",alpha=1.0,size=4,cex = 4) +
  scale_y_continuous(name="Proportion of Genome",breaks = c(0.5,0.55,0.6,0.65,0.7), limits = c(0.5,0.7)) +
  theme_bw(base_size = 20)
p
dev.off()


### Plotting all data
dat=read.table('./AllCrossoverData.txt',header = T)
dat
pdf("AllCOs.pdf",width = 4,height = 5)
p<-ggplot(data=dat, aes(x=NumMeiosis, y=CO)) +
  geom_beeswarm(shape=21,fill="black",color="black",alpha=1.0,size=4,cex = 4) +
  scale_x_discrete(name="Number of Meioses",labels=c("D_3"=">=3","C_2"="2","B_1"="1","A_0"="0")) +
  scale_y_continuous(name="Number of Crossovers",breaks = c(0,5,10,15,20,25,30,35,40)) +
  theme_bw(base_size = 20)
p
dev.off()

### Calculating the correlation between the number of crossovers detected and the number of meioses
datNoPos=read.table('./AllCrossoverData_NoPosCTRL.txt',header = T)
datNoPos
datNoPos$NumMeiosis
datNoPosNum=read.table('./AllCrossoverData_NoPosCTRL_Numeric.txt',header = T)
cor(datNoPosNum$NumMeiosis,datNoPosNum$CO)
linearModel<-lm(CO ~ NumMeiosis, data=datNoPosNum)
linearModel
summary(linearModel)


##################################
##################################
##
##
##
##        Karyoplots
##
##
##
##################################
##################################

library(karyoploteR)

### Defining genome
Obir5.4withChr<-toGRanges("../../data/KaryoplotR/TestObir5_4_withChr.txt")
ObirTigswithChr<-toGRanges("../../data/KaryoplotR/ContigFile_withChr.txt")

### Loading in data
# first the 1 meiosis comparisons (Mother-Daughter Pairs)
NOVT504vNOVT505<-toGRanges('../Comparisons/NOV_T504_v_NOV_T505/AllChromsAllSteps_ManuallyScreened_gRangesInput.txt')
NOVT504vNOVT505
FEBT506vFEBT505<-toGRanges('../Comparisons/FEB_T506_v_FEB_T505/AllChromsAllSteps_ManuallyScreened_gRangesInput.txt')
FEBT506vFEBT505
FEBT506vFEBT507<-toGRanges('../Comparisons/FEB_T506_v_FEB_T507/AllChromsAllSteps_ManuallyScreened_gRangesInput.txt')
FEBT506vFEBT507
FEBT508vFEBT509<-toGRanges('../Comparisons/FEB_T508_v_FEB_T509/AllChromsAllSteps_ManuallyScreened_gRangesInput.txt')
FEBT508vFEBT509
FEBT508vFEBT510<-toGRanges('../Comparisons/FEB_T508_v_FEB_T510/AllChromsAllSteps_ManuallyScreened_gRangesInput.txt')
FEBT508vFEBT510
### then the 2 meiosis comparisons (sister-sister pairs)
NOVT501vNOVT503<-toGRanges('../Comparisons/NOV_T501_v_NOV_T503/AllChromsAllSteps_ManuallyScreened_gRangesInput.txt')
NOVT501vNOVT503
FEBT511vFEBT512<-toGRanges('../Comparisons/FEB_T511_v_FEB_T512/AllChromsAllSteps_ManuallyScreened_gRangesInput.txt')
FEBT511vFEBT512
FEBT513vFEBT514<-toGRanges('../Comparisons/FEB_T513_v_FEB_T514/AllChromsAllSteps_ManuallyScreened_gRangesInput.txt')
FEBT513vFEBT514
FEBT515vFEBT516<-toGRanges('../Comparisons/FEB_T515_v_FEB_T516/AllChromsAllSteps_ManuallyScreened_gRangesInput.txt')
FEBT515vFEBT516
### then the >=3 meiosis comparisons (All possible comparisons between pedigrees)
AllPairwise<-toGRanges('../AllCOsAllPairwiseUnknownPedigreePairs_ForGRanges.txt')
AllPairwise
### Then sites at which crossovers should be detectable
DetectableHet<-toGRanges('../Detectabilities/NOVT504vNOVT505_AllPairsIntersection_UniqueSNPs_gRanges.txt')
DetectableHet

### Plotting Whole Genome
pdf("/Users/kiplacy/Desktop/TELLkaryoplot_noLabels_AllPairwise.pdf",width=10,height=4.5)
kp<-plotKaryotype(genome = Obir5.4withChr, cytobands = ObirTigswithChr,plot.type = 4)
at <- autotrack(current.track = 11,total.tracks = 11)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=NOVT504vNOVT505,window.size = 40000,col='black',border='black')
at <- autotrack(current.track = 10,total.tracks = 11)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=FEBT506vFEBT505,window.size = 40000,col='black',border='black')
at <- autotrack(current.track = 9,total.tracks = 11)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=FEBT506vFEBT507,window.size = 40000,col='black',border='black')
at <- autotrack(current.track = 8,total.tracks = 11)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=FEBT508vFEBT509,window.size = 40000,col='black',border='black')
at <- autotrack(current.track = 7,total.tracks = 11)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=FEBT508vFEBT510,window.size = 40000,col='black',border='black')
at <- autotrack(current.track = 6,total.tracks = 11)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=NOVT501vNOVT503,window.size = 40000,col='black',border='black')
at <- autotrack(current.track = 5,total.tracks = 11)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=FEBT511vFEBT512,window.size = 40000,col='black',border='black')
at <- autotrack(current.track = 4,total.tracks = 11)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=FEBT513vFEBT514,window.size = 40000,col='black',border='black')
at <- autotrack(current.track = 3,total.tracks = 11)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=FEBT515vFEBT516,window.size = 40000,col='black',border='black')
at <- autotrack(current.track = 2,total.tracks = 11)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=AllPairwise,window.size = 4000,col='black',border='black')
at <- autotrack(current.track = 1,total.tracks = 11)
kpPlotDensity(kp, r0=at$r0, r1=at$r1, data=DetectableHet,window.size = 200000,col='grey',border=NA)
dev.off()

