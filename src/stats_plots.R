#! /usr/bin/env Rscript

message('============================================================================================\nStarting statistical analyses...\n\nLoading input data...')
library(ggplot2)
setwd('../data/processed')
Sign <- readRDS("Sign.rds")
Rand <- readRDS("Rand.rds")
SigMerged <- readRDS("SigMerged.rds")
RandMerged <- readRDS("RandMerged.rds")
TopQval <- readRDS("TopQval.rds")
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

### Construcción de los Volcano Plots. Código adaptado de https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
# RS Significativos y sus OR

VolcanoPlot <- function (x, y, ORcutoff=0.6, pvalcutoff=0.001, title){
	df <- data.frame("Odds_Ratio"=10**x, "Pvalue"=y)
	svg(file = paste(title,"_volcano_plot.svg"))
	plot <- ggplot(data=df, aes(x=log2(Odds_Ratio), y=-log10(Pvalue))) +
			geom_point() +
			theme_minimal() +
			geom_vline(xintercept=c(-ORcutoff, ORcutoff), col="red") +
			geom_hline(yintercept=-log10(pvalcutoff), col="red") +
			ggtitle(title) +
			theme(plot.title = element_text(size=20, hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
	print(plot)
	invisible(dev.off())
}

VolcanoPlot(x=BetaSign, y=pvalSign, ORcutoff=0.6, pvalcutoff=0.001, title="Significant_SNPs")
VolcanoPlot(x=BetaRand, y=pvalRand, ORcutoff=0.6, pvalcutoff=0.001, title="Random_SNPs")
message('Volcano plots built!\n')

### Chi2 para los eQTL

tableEQTL <- function(sigset,randset,fullset) {
	SigeQTL <- length(unique(sigset$rs_id))
	RandeQTL <- length(unique(randset$rs_id))
	totalSNP <- unique(fullset$RS)
	totalSNP <- length(totalSNP[!is.na(totalSNP)])
	df1 <- data.frame("eQTL"=c(SigeQTL,RandeQTL))
	df2 <- data.frame("non_eQTL"=rep(totalSNP,2))
	res = cbind(df1,df2)
	rownames(res) = c("Significant SNPs", "Random SNPs")
	return(res)
}

Chi2 = function (Data) {
	Chisq = chisq.test(Data)
	n = sum(Data)
	gl = min(dim(Data) - 1)
	VCramer = sqrt((Chisq$statistic)/(n * gl))
	result = data.frame("X2" = Chisq$statistic, "p-value" = Chisq$p.value, "Cramer_V" = VCramer, row.names = NULL)
	return(result)
	}

message('Performing eQTL chi-squared test...')

eQTL = tableEQTL(SigMerged, RandMerged, Sign)

Chisqres = Chi2(eQTL)

# Magnitud del efecto. Códigos tomados y adaptados de https://stats.stackexchange.com/questions/427864/how-to-calculate-an-effect-size-for-chi-square-in-r


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
setwd("./tables")
write.table(eQTL, file="eQTL_table.txt", sep ="\t", row.names = TRUE, col.names = TRUE)

message('Chi-squared test done!\n')

message("Statistical analysis done!\n")
