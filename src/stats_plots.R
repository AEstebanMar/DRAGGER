#! /usr/bin/env Rscript

message('\n\n============================================================================================\nStarting statistical analyses...\n\nLoading input data...')
library(ggplot2)
setwd('../data/processed')
Sign = readRDS("Sign.rds")
Rand = readRDS("Rand.rds")
TopQval = readRDS("TopQval.rds")
save.image('Debugging.RData')
setwd('../../output/plots')
message('Input loaded!\n')

### Extract data for statistican analysis. Some of it has already been loaded in the analysis script

BetaSign = Sign[,"beta"]
BetaRand = Rand[,"beta"]
pvalSign = Sign[,"P"]
pvalRand = Rand[,"P"]
ORSign = 10**BetaSign
ORRand = 10**BetaRand

genes = TopQval [,"gene_name"]
qvals = TopQval [,"qval"]
Slopes = TopQval [,"slope"]

# Build and plot data frames

message('Building Volcano plots...')

VolcanoSign = data.frame("Pvalue" = pvalSign, "Odds_Ratio" = ORSign)
VolcanoRand = data.frame("Pvalue" = pvalRand, "Odds_Ratio" = ORRand)
VolcanoGenes = data.frame("Qvalue" = qvals, "Slope" = Slopes)
### Construcción de los Volcano Plots. Código adaptado de https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
# RS Significativos y sus OR
png(file = "VolcanoSign.png")
ggplot(data=VolcanoSign, aes(x=log2(Odds_Ratio), y=-log10(Pvalue))) + geom_point() + theme_minimal() + geom_vline(xintercept=c(-0.6, 0.6), col="red") + geom_hline(yintercept=-log10(0.001), col="red") + ggtitle("SNP Significativos") + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
invisible(dev.off())
# RS Aleatorios y sus OR
png(file = "VolcanoRand.png")
ggplot(data=VolcanoRand, aes(x=log2(Odds_Ratio), y=-log10(Pvalue))) + geom_point() + theme_minimal() + geom_vline(xintercept=c(-0.6, 0.6), col="red") + geom_hline(yintercept=-log10(0.001), col="red") + ggtitle("SNP Aleatorios") + theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
invisible(dev.off())
### Chi2 para los eQTL
eQTLs = data.frame("eQTL" = c(439, 107),"no eQTL" = c(41480, 41812))
rownames(eQTLs) = c("SNP Significativos", "SNP aleatorios")

message('Volcano plots built!\n')

# Magnitud del efecto. Códigos tomados y adaptados de https://stats.stackexchange.com/questions/427864/how-to-calculate-an-effect-size-for-chi-square-in-r

Chi2 = function (Datos) {
	Chisq = chisq.test(Datos)
	n = sum(Datos)
	gl = min(dim(Datos) - 1)
	VCramer = sqrt((Chisq$statistic)/(n * gl))
	result = data.frame("X2" = Chisq$statistic, "p-value" = Chisq$p.value, "Cramer_V" = VCramer, row.names = NULL)
	return(result)
	}

message('Performing eQTL chi-squared test...')

save.image("Debugging.RData")

Chisqres = Chi2(eQTLs)

setwd("../")
cat("\n==========================================\nStatistics summary\n==========================================\n", file="Summary.txt", append = TRUE)
cat("Chi squared results for eQTL\n---------------------------------------------------------------\n", file="Summary.txt", append = TRUE)
cat("Chi squared value:\n", file="Summary.txt", append = TRUE)
cat(Chisqres$X2, file="Summary.txt", append = TRUE)
cat("\np-value:\n", file="Summary.txt", append = TRUE)
cat(Chisqres$p.value, file="Summary.txt", append = TRUE)
cat("\nMagnitude of effect (Cramer-V):\n", file="Summary.txt", append = TRUE)
cat(Chisqres$Cramer_V, file="Summary.txt", append = TRUE)
cat("\n---------------------------------------------------------------\n", file="Summary.txt", append = TRUE)

message('Chi-squared test done!\n')

message("Statistical analysis done!\n")
