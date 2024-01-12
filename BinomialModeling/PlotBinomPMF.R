# This script uses binomial models to understand the number of chromosomes expected to incur segmental LOH
# and the proportion of offspring expected to bear segmental LOH not found in their mother

# In later sections we explore how developmental mortality of individuals bearing segmental LOH
# and co-segregation of recombined chromatids can influence rates of LOH.
# Based on these models, we identify the values of each of these parameters that would be required
# to explain the extent of heterozygosity maintenance observed in our study.

library(dplyr)
library(ggplot2)

########################################  
###                                  ###
##                                    ##
#       Simple Binomial PMF Plots      #
##                                    ##
###                                  ###  
########################################

# First, we plot Binomial Probability Mass Functions (PMFs) to illustrate the expected 
# number of chromosomes bearing CO-induced LOH assuming a single CO per chromosome, 
# independent assortment, and random segregation

# O. biroi has 14 chromosomes. Per generation, anywhere between zero and 14 chromosomes can
# incur segmental loss of heterozygosity
Chromosomes=0:14
# Set up the PDF file to hold the PMF plot
pdf("RplotBinomPMFnumchromsLOH.pdf",width = 5,height = 3)
# Plot the binomial PMF using "dbinom()", with 14 "trials", and Mendelian segregation (probability=0.5)
plot(Chromosomes,dbinom(Chromosomes, size=14, prob=.5),type='h',lwd=3,xlab = "Number of Chromosomes with Crossover-Distal LOH", ylab = "Probability", las=1)
# Close the pdf file
dev.off()
# Determine the probability that at least one segmental LOH occurs
sum(dbinom(1:14,size=14,prob=(.5*(5.5/14))))

# Later in the paper, we repeat these steps, assuming fewer crossovers per meiosis 
# (our empirically-observed number--5.5--which is likely to be an underestimate)

# Set up the PDF file to hold the PMF plot
pdf("RplotBinomPMFnumchromsLOH_5.5co_0.5LOH.pdf",width = 5,height = 3)
# Plot the binomial PMF using "dbinom()", with 14 "trials", and Mendelian segregation (probability=0.5)
plot(Chromosomes,dbinom(Chromosomes, size=14, prob=(.5*(5.5/14))),type='h',lwd=3,xlab = "Number of Chromosomes with LOH per Meiosis", ylab = "Probability", las=1)
# Close the pdf file
dev.off()
# Determine the probability that at least one segmental LOH occurs
sum(dbinom(1:14,size=14,prob=(.5*(5.5/14))))
sum(Chromosomes * dbinom(Chromosomes,size=14,prob=(.5*(5.5/14))))

# Then, using the same number of crossovers per meiosis (5.5), we substitute a 
# different rate at which segmental LOH occurs following a crossover. Under Mendelian segregation, this would be 0.5
# as used above. But under co-segregation of reciprocally recombined crossover products,
# this number could be lower. Here we use a 0.1 probability of segmental LOH following a crossover,
# which would be the case if the co-segregation probability was 0.9.

# Set up the PDF file to hold the PMF plot
pdf("RplotBinomPMFnumchromsLOH_5.5co_0.1LOH.pdf",width = 5,height = 3)
# Plot the binomial PMF using "dbinom()", with 14 "trials", and Mendelian segregation (probability=0.5)
plot(Chromosomes,dbinom(Chromosomes, size=14, prob=(.1*(5.5/14))),type='h',lwd=3,xlab = "Number of Chromosomes with LOH per Meiosis", ylab = "Probability", las=1)
# Close the pdf file
dev.off()
# Determine the probability that at least one segmental LOH occurs
sum(dbinom(1:14,size=14,prob=(.1*(5.5/14))))

#####################################################################################################  
###                                                                                               ###
##                                                                                                 ##
#       Exploring whether developmental mortality of homozygotes could maintain heterozygosity      #
##                                                                                                 ##
###                                                                                               ###  
#####################################################################################################

# Because we found that segmental LOH is incredibly rare, despite the fact that crossovers are common, 
# we wanted to determine the mechanism that maintains heterozygosity. 
# Thus, we extended the above binomial modelling logic to include additional aspects of the biological context.
# First, we consider the possibility that developmental mortality of individuals bearing segmental LOH maintains heterozygosity,
# We incorporated this into a model written as a function below ("CalculateLOHProportion".)
# For this analysis we are assuming Mendelian segregation.

# Define function to calculate proportion of offspring expected to bear segemental losses of heterozygosity (LOH)
CalculateLOHProportion <- function(CrossoversPerChromosome, MortalityRate) {
  # Number of chromosomes in the species' genome
  NumChromosomes <- 14
  # The probability of loss of heterozygosity per chromosome is one minus the co-segregation probability. 
  # For Mendelian segregation, the co-segregation probability is 0.5. 
  ProbabilityLOHPerChromosome <- CrossoversPerChromosome * 0.5
  # The probability that the egg bears a segmental loss of heterozygosity is taken from the binomial PMF, as used above
  # Here we determine the probability that at least one segmental LOH occurs in a meiosis
  ProbabilityEggBearsLOH <- sum(dbinom(1:NumChromosomes, NumChromosomes, ProbabilityLOHPerChromosome))
  # Adjusting mortality rate to not exceed the probability of LOH occurrence
  # This is because in our model, all mortality comes from segmental LOH.
  MortalityRate <- min(MortalityRate, ProbabilityEggBearsLOH)
  # Calculating the probability that pupae bear LOH, adjusting for mortality that may have occurred during development
  # from the egg stage through the larval stage until pupation
  ProbabilityPupaBearsLOH <- max(0, (ProbabilityEggBearsLOH - MortalityRate)) / ((1 - ProbabilityEggBearsLOH) + (ProbabilityEggBearsLOH - MortalityRate))
  return(ProbabilityPupaBearsLOH)
}

# Define the values for crossovers per chromosome and mortality rate, and Cosegregation probability
# 0.4 is our empirical value (5.5 per meiosis divided by 14 chromosomes)
# 0.3 and 0.49 are the lower and upper bounds of the 95% Confidence Interval
CrossoversPerChromosome_Values <- c(0.3, 0.4, 0.49)
# We will calculate the probability that a pupae bears segmental LOH across 
# essentially all possible values of the developmental lethality rate
DevelopmentalLethalityRate_Values <- seq(0, 0.99, by = 0.0125)

# Now we create a dataframe to calculate the probabilities across all of these values
df <- expand.grid(CrossoversPerChromosome_Values, DevelopmentalLethalityRate_Values)
# Setting the column names
colnames(df) <- c('Crossovers per Chromosome', 'Developmental Lethality Rate')
# Apply the CalculateLOHProportion function to each row in the data frame
df$`Probability Pupa Bears LOH` <- mapply(CalculateLOHProportion, df$`Crossovers per Chromosome`, df$`Developmental Lethality Rate`)
# Convert the 'Crossovers per Chromosome' column to a factor
df$`Crossovers per Chromosome` <- as.factor(df$`Crossovers per Chromosome`)
# Inspect the dataframe
head(df)

# Plot
# Setting up a pdf to save the plot
pdf("./results/CrossoversHomLethalProbNoLOHSamples_MortalityRates_0.4coPerChrom_Shorter.pdf", width = 6, height = 3)
# Plot the calculated values using geom_line() in ggplot2
ggplot(df, aes(x = `Developmental Lethality Rate`, y = `Probability Pupa Bears LOH`, color = `Crossovers per Chromosome`, shape = `Crossovers per Chromosome`)) +
  geom_line() +
  # Add the empirical mean and 95% confidence interval of the developmental mortality rate
  annotate(geom = "rect", xmin = 0.11, xmax = 0.24, ymin = -Inf, ymax = Inf, fill = "#DC2780", alpha = 0.2) +
  geom_vline(aes(xintercept = 0.17), color = "#DC2780") +  # Vertical line at x = 0.14
  # Add the empirical mean and 95% confidence interval of the proportion of pupae with at least one segmental LOH
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 0.2, fill = "#FF9F00", alpha = 0.2) +
  geom_hline(aes(yintercept = 0), color = "#FF9F00") +  # Horizontal line at x = 0.0
  # Plotting aesthetics
  theme_bw(base_size = 15) +
  labs(y = "Proportion of Pupae With\nat Least One Loss-of-Heterozygosity", x = "Developmental Lethality Rate") +
  scale_color_brewer(name = "Crossovers \nper Chromosome\nper Meiosis", palette = "Dark2") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), minor_breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), minor_breaks = seq(0, 1, by = 0.1), limits=c(0, 1))
dev.off()
# I shaded the area between the upper and lower bound of the 95% CI for proportion of pupae with at least one segmental LOH in adobe illustrator.

# Finding the minimal proportion of pupae bearing segmental LOH
# if the true developmental mortality rate fall withins the 95% confidence interval of our empirical estimate
CrossoversValue <- 0.3
MortalityValue <- 0.24
CalculateLOHProportion(CrossoversValue, MortalityValue)
# value is 0.864

# Finding the developmental mortality rate required to produce Proportions of pupae bearing segmental LOH
# that fall within the 95% confidence interval of our empirical estimate
CrossoversValue <- 0.3
ProbabilityPupaBearsLOHValue <- 0.2
# Defining a function (RootFunction) to find the developmental mortality rate given these conditions
# This function will be used in the root-finding method.
# The function returns the difference between the calculated LOH proportion (Output) and the target empirical estimate.
RootFunction <- function(MortalityValue) {
  Output <- CalculateLOHProportion(CrossoversValue, MortalityValue)
  return(Output - ProbabilityPupaBearsLOHValue)
}
# Using uniroot() to find the root of RootFunction.This method finds the value of MortalityValue for which RootFunction returns zero.
# We search only between 0 and 1, which are the possible values for the developmental mortality rate.
# The returned value will be the mortality rate that aligns the model's LOH proportion with the empirical estimate.
Solution <- uniroot(RootFunction, interval = c(0, 1))
# Display the result
print(Solution$root)
# result is 0.871


######################################################################################################  
###                                                                                                ###
##                                                                                                  ##
#       Exploring whether co-segregation of recombined chromatids could maintain heterozygosity      #
##                                                                                                  ##
###                                                                                                ###  
######################################################################################################


# Now, we will modify the function from above to include Cosegregation probability as a factor. 
# For Mendelian segregation, the co-segregation probability is 0.5. 
# Here we explore what happens if chromosomal segregation is non-Mendelian, and
# co-segregation probability is increased

# Define function to calculate proportion of offspring expected to bear segemental losses of heterozygosity (LOH)
CalculateLOHProportion <- function(CrossoversPerChromosome, MortalityRate, CosegregationProb) {
  # Number of chromosomes in the species' genome
  NumChromosomes <- 14
  # The probability of loss of heterozygosity per chromosome is one minus the co-segregation probability. 
  # For Mendelian segregation, the co-segregation probability is 0.5. 
  # With higher cosegregation probabilities, the probability of segmental LOH per crossover decreases
  ProbabilityLOHPerChromosome <- CrossoversPerChromosome * (1 - CosegregationProb)
  # The probability that the egg bears a segmental loss of heterozygosity is taken from the binomial PMF, as used above
  # Here we determine the probability that at least one segmental LOH occurs in a meiosis
  ProbabilityEggBearsLOH <- sum(dbinom(1:NumChromosomes, NumChromosomes, ProbabilityLOHPerChromosome))
  # Adjusting mortality rate to not exceed the probability of LOH occurrence
  # This is because in our model, all mortality comes from segmental LOH.
  MortalityRate <- min(MortalityRate, ProbabilityEggBearsLOH)
  # Calculating the probability that pupae bear LOH, adjusting for mortality that may have occurred during development
  # from the egg stage through the larval stage until pupation
  ProbabilityPupaBearsLOH <- max(0, (ProbabilityEggBearsLOH - MortalityRate)) / ((1 - ProbabilityEggBearsLOH) + (ProbabilityEggBearsLOH - MortalityRate))
  return(ProbabilityPupaBearsLOH)
}


# Create a large grid of values for a fine heatmap of "ProbabilityPupaBearsLOH" values
# across a range of values of "CosegregationProb" and "MortalityRate"
CosegregationVals <- seq(0.5, 1, length.out = 1000)
MortalityVals <- seq(0, 0.999, length.out = 1000)
# Set as a constant the empirical estimate of crossovers per chromosome per meiosis
CrossoversPerChromosomeEmpirical=5.5/14
df <- expand.grid(CosegregationProb = CosegregationVals, MortalityRate = MortalityVals)

# Apply the function across this grid of values to calculate "ProbabilityPupaBearsLOH"
df$`Probability Pupa Bears LOH` <- mapply(CalculateLOHProportion, CrossoversPerChromosomeEmpirical, results$MortalityRate, results$CosegregationProb)
head(df)

# NEXT WE WANT TO OBTAIN A CURVE THAT DEPICTS THE LOWEST VALUE OF CO-SEGREGATION PROBABILITY
# FOR WHICH THE PROBABILITY THAT A PUPA BEARS LOH FALLS WITHIN THE 95% CI OF THE EMPIRICAL ESTIMATE

# Filtering to include only rows where 'Probability Pupa Bears LOH' is <= 0.2
# Thus we are only including rows of the data frame for which that probability falls within the 
# 95% confidence interval of our empirical estimate.
dfProbLOHBelowThreshold <- df %>% 
  filter(`Probability Pupa Bears LOH` <= 0.2)
# For each MortalityRate value, find the value of CosegregationProb 
# corresponding to the highest 'Probability Pupa Bears LOH'
# Since we are working within the thresholded dataframe, this will finde the Cosegregation probability
# Value that produces the highest "Proportion" that is still within the 95% CI of our empirical estimate.
# This can be thought of as the minimal Cosegregation probability required to maintain heterozygosity
# to the degree observed in the study for each value of the homozygous lethality rate.
CosegregationValHighestProbLOH <- dfProbLOHBelowThreshold %>%
  group_by(MortalityRate) %>%
  filter(rank(-`Probability Pupa Bears LOH`, ties.method = "first") == 1) %>% # this ensures we select the highest 'value'
  select(MortalityRate, CosegregationProb, `Probability Pupa Bears LOH`) %>%
  arrange(MortalityRate)

# Cutting the "CosegregationValHighestProbLOH" curve once it reaches 0.5. 
# This ensures accurate plotting.
# We accomplish this by first re-arranging the HighestCosegregation data frame 
# to sort it by the CosegregationProb in descending order.
CosegregationValHighestProbLOH <- CosegregationValHighestProbLOH %>%
  arrange(desc(CosegregationProb))
# Then we find the first occurrence in the data where the CosegregationProb is exactly 0.5.
FirstOccurrence <- which(CosegregationValHighestProbLOH$CosegregationProb == 0.5)[1]
# Then we check if such an occurrence (CosegregationProb = 0.5) exists, and if it does, 
# we truncate the HighestCosegregation data frame to include only the data up to that first occurrence.
if (!is.na(FirstOccurrence)) {
  CosegregationValHighestProbLOH <- CosegregationValHighestProbLOH[1:FirstOccurrence, ]
}

# Plot the relationship
pdf("./HeatMapAndMinimalCosegregationProb_0.4COsPerChrom_PropPupLOH.pdf", width = 6, height = 4)
ggplot() +
  # First plot a heat map to show how "ProbabilityPupaBearsLOH" varies across Co-segregation Prob and Developmental Mortality Rate
  geom_tile(data = df, aes(x = MortalityRate, y = CosegregationProb, fill = `Probability Pupa Bears LOH`)) +
  scale_fill_distiller(palette = 'YlGn', direction = 1) +
  labs(fill = 'Prob.', y = "Co-Segregation Probability", x = "Homozygous Lethality Rate") +
  # Add the empirical estimate and 95% confidence interval of the developmental mortality rate
  annotate(geom = "rect", xmin = 0.11, xmax = 0.24, ymin = -Inf, ymax = Inf, fill = "#DC267F", alpha = 0.3) +
  geom_vline(aes(xintercept = 0.17), color = "#DC267F") +  # Vertical line at x = 0.14 
  # Add the line showing, for each value of the developmental mortality rate,
  # the lowest value of co-segregation prob. for which the "ProbabilityPupaBearsLOH" is 0.2 or less,
  # that is, it falls within the 95% confidence interval of the empirical estimate.
  geom_line(data = CosegregationValHighestProbLOH, aes(x = MortalityRate, y = CosegregationProb)) +
  # Adjust plotting aesthetics
  scale_x_continuous(breaks = seq(0, 1, by = 0.2), minor_breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0.5, 1, by = 0.1), minor_breaks = seq(0.5, 1, by = 0.05), limits=c(0.5, 1)) +
  theme_bw(base_size = 15)
dev.off()

# Defining a function ('root_function') to find the cosegregation probability where 
# the calculated proportion of "Probability Pupa Bears LOH" equals 0.2 (the upper bound of the
# 95% confidence interval of our empirical estimate) for a mortality rate of 0.24 (the upper 
# bound of the 95% confidence interval of our empirical mortality rate estimate)
# This allows us to identify the lowest co-segregation probability required to explain our data.
# This function calculates the difference between the calculated proportion and the target proportion (0.2).
RootFunction <- function(CosegregationProb) {
  CalculateLOHProportion(CrossoversPerChromosomeEmpirical, 0.24, CosegregationProb) - 0.2
}
# Using uniroot() to find the root of RootFunction.This method finds the value of CosegregationProb for which RootFunction returns zero.
# The search is conducted within the range from 0.5 to 1, the relevant values of co-segregation probability
# The returned value will be the Cosegregation Probability that aligns the model's LOH proportion with the empirical estimate.
Solution <- uniroot(RootFunction, c(0.5, 1))
# Extracting the root (cosegregation probability) from the result.
CosegregationProbThatYields0_2ProbPupaLOHatMortality_0_24 <- Solution$root
print(CosegregationProbThatYields0_2ProbPupaLOHatMortality_0_24)
# the result is 0.911

# Repeating the process with a mortality rate of 0 to see what level of
# cosegregation probability would be required to explain heterozygosity maintenance if 
# homozygous lethality played no role
RootFunction <- function(CosegregationProb) {
  CalculateLOHProportion(CrossoversPerChromosomeEmpirical, 0.0, CosegregationProb) - 0.2
}
Solution <- uniroot(RootFunction, c(0.5, 1))
CosegregationProbThatYields0_2ProbPupaLOHatMortality_0_0 <- Solution$root
print(CosegregationProbThatYields0_2ProbPupaLOHatMortality_0_0)
# the result is 0.9597
