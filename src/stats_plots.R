#! /usr/bin/env Rscript

VolcanoPlot <- function (x, y, ORcutoff=0.6, pvalcutoff=0.001, title){ # Code adapted from https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
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

Chi2 <- function (Data) {
	Chisq = chisq.test(Data)
	n = sum(Data)
	gl = min(dim(Data) - 1)
	VCramer = sqrt((Chisq$statistic)/(n * gl)) # Effect size code adapted from https://stats.stackexchange.com/questions/427864/how-to-calculate-an-effect-size-for-chi-square-in-r
	result = data.frame("X2" = Chisq$statistic, "p-value" = Chisq$p.value, "Cramer_V" = VCramer, row.names = NULL)
	return(result)
}

groupSigVsRandBarplot <- function(set1, set2, attr, title) {
	SigData <- table(set1[,attr])
	SigData <- cbind(seq(1:length(SigData)),SigData,rep("Significant",length(SigData)))
	RandData <- table(set2[,attr])
	RandData <- cbind(seq(1:length(RandData)),RandData,rep("Random",length(RandData)))
	colnames(SigData) = colnames(RandData) <- c("Chr","Frequency","Group")
	totalData <- rbind(SigData, RandData)
	df1 <- data.frame(totalData)
	df1$Frequency = as.numeric(df1$Frequency)
	df1$Chr = as.numeric(df1$Chr)
	svg(file = paste(title,".svg"))
	plot <- ggplot(df1, aes(x=Chr, y=Frequency, fill=Group)) + 
			scale_x_continuous(breaks = seq(1, 23, by = 1)) +
   			geom_bar(position="dodge", stat="identity") +
   			scale_fill_manual(values = c("#000000","#7da1c4")) +
   			ggtitle(title) + xlab("Chromosome") +
   			theme(plot.title = element_text(hjust = 0.5)) +
   			scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
    print(plot)
	invisible(dev.off())
}

message('============================================================================================\nStarting statistical analyses...\n\nLoading input data...')
library(ggplot2)
setwd('../data/processed')
SignData <- list("SNP"=readRDS("Sign.rds"),"eQTL"=readRDS("SigMerged.rds"))
RandData <- list("SNP"=readRDS("Rand.rds"),"eQTL"=readRDS("RandMerged.rds"))
TopQval <- readRDS("TopQval.rds")
setwd('../../output/plots')

#### Build data frames and plots

message('Building bar plots\n')
groupSigVsRandBarplot(SignData$SNP, RandData$SNP, "CHR", "SNP distribution per chromosome")
groupSigVsRandBarplot(SignData$eQTL, RandData$eQTL, "CHR", "eQTL distribution per chromosome")
message('Building Volcano plots\n')
VolcanoPlot(x=SignData$SNP$beta, y=SignData$SNP$P, title="Significant_SNPs")
VolcanoPlot(x=RandData$SNP$beta, y=RandData$SNP$P, title="Random_SNPs")
# eQTL Chi squared analysis
message('Performing eQTL chi-squared test\n')
eQTL = tableEQTL(SignData$eQTL, RandData$eQTL, SignData$SNP)
Chisqres = Chi2(eQTL)

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

message("Statistical analysis done!\n")
