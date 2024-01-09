setwd("~/../../../Volumes/antqueen/booster/KDL20240109_UploadCodeParthPaper/UnselfishMeioticDrive/TransgenicLineCharacterization/")# Load the data
Data <- read.table("./SourceData_TransgenicLineFitness.txt", header = T)

# Extract the column that contains the proportion of DsRed fluorescent offspring from the dataset.
DsRedProp <- Data$ProportionOfOffspringDsRed

# Set the expected proportion value for comparison purposes in the statistical test.
# The expected proportion is 1/20 because in each replicate colony 1 individual was dsRed and
# 19 individuals were wild-types
ExpectedProportion=0.05

# Perform a one-tailed one-sample Wilcoxon test to determine if the observed proportion of DsRed fluorescent offspring
# is significantly less than the expected proportion. This test is appropriate because the data are not normally distributed.
# 'mu' is set to the expected proportion, and 'alternative = "less"' specifies that we are testing for a decrease.
wilcox.test(DsRedProp, mu=ExpectedProportion, alternative = "less")

# Output the mean proportion of DsRed fluorescent offspring and its standard deviation
cat("Mean DsRed Offspring Proportion: ", mean(DsRedProp), "±", sd(DsRedProp), "SD\n")

# NOTE: Data was plotted in GraphPad Prism, which is not included in this script.
# But for a quick visualization of the data, I have provided R code that uses ggplot2 to plot these data below
library(ggplot2)
ggplot(data.frame(x = "DsRedProp", y = DsRedProp), aes(x=x, y=y)) +  
  # plotting the data
  geom_jitter(shape = 1, cex=5, size = 3, color = 'black') +
  # adding the mean
  geom_hline(yintercept = ExpectedProportion, linetype="dashed", color = "black") +
  # formatting the plot to my preferred aesthetics
  labs(y="Proportion of offspring dsRed") +
  theme_bw(base_size = 15) +
  theme(legend.position="none",    # Removes legend
        axis.title.y=element_blank(),  # Removes x-axis title
        axis.text.y=element_blank()) +  # Removes x-axis text
  scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.1)) +
  coord_flip()

