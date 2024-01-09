# The purpose of this file is to calculate the mean and 95% confidence interval of 
# the developmental mortality rates from egg to pupation observed in replicate colonies

# Here are the raw developmental mortality rates
MortalityRates <- c(0.11, 0.11, 0.15, 0.15, 0.17, 0.22, 0.3)

# Use the t.test() function to calculate mean and 95% confidence interval
result <- t.test(MortalityRates)
MeanMortalityRate <- round(result$estimate, 3)
CIMortalityRate_Lower <- round(result$conf.int[1], 3)
CIMortalityRate_Upper <- round(result$conf.int[2], 3)

# Print results
cat("Mean:", MeanMortalityRate, "\n")
cat("95% CI:", CIMortalityRate_Lower, "to", CIMortalityRate_Upper, "\n")

# Plot the data
library(ggplot2)
library(ggbeeswarm)
pdf("./DevelopmentalMortalityRates.pdf", width = 3.5, height = 1)
ggplot(data.frame(x = "MortalityRates", y = MortalityRates), aes(x=x, y=y)) +
	# plotting the empirical data
	geom_beeswarm(shape = 1, cex=12, size = 3, color = 'black', corral = 'gutter', corral.width = 0.35) +
	# plotting the 95% CI
	geom_errorbar(aes(x = "MortalityRates", ymin = CIMortalityRate_Lower, ymax = CIMortalityRate_Upper), width = 0.3) +
	# adding the mean
	geom_hline(yintercept = MeanMortalityRate, linetype="dashed", color = "red") +
	# adding the mortality rate that would be required for complete heterozygosity maintenance if
	# segregation was Mendelian and one crossover occurs per chromosome
	geom_point(aes(y = 0.99), shape = 13, size = 3, color = "red") + 
	# formatting the plot to my preferred aesthetics
	labs(y="Developmental Mortality Rate") +
	theme_bw(base_size = 15) +
	theme(legend.position="none", # Removes legend
	      axis.title.y=element_blank(), # Removes x-axis title
	      axis.text.y=element_blank()) + # Removes x-axis text
	scale_y_continuous(limits = c(0, 1.0), breaks = seq(0, 1.0, 0.1)) +
	coord_flip()
dev.off()
