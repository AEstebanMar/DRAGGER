#! /usr/bin/env Rscript

VolcanoPlot <- function (x, y, ORcutoff=0.6, pvalcutoff=0.001, title, fontsize=32) { # Code adapted from https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
	
	df <- data.frame("Odds_Ratio"=10**x, "Pvalue"=y)

	tiff(file= paste0(title,".tiff"), units="cm",
		width=35, height=30, res=300)

	plot <- ggplot(data=df, aes(x=log2(Odds_Ratio), y=-log10(Pvalue))) +
			geom_point() +
			theme_minimal() +
			geom_vline(xintercept=c(-ORcutoff, ORcutoff), col="red") +
			geom_hline(yintercept=-log10(pvalcutoff), col="red") +
			ggtitle(title) +
			theme(
				plot.title=element_text(size=fontsize*1.2, hjust = 0.5),
				plot.subtitle=element_text(hjust = 0.5),
   				axis.title = element_text(size = fontsize),
				legend.text = element_text(size = fontsize),
				axis.text = element_text(size = fontsize*0.8),
				legend.title=element_blank()
				)

	print(plot)
	invisible(dev.off())
}

tableForChisq <- function(sigset,randset,fullset,testGroup,fullGroup,fixedTotal=FALSE) {

	sigData <- length(unique(sigset$RS))
	randData <- length(unique(randset$RS))

	if(!fixedTotal){
		totalData <- length(fullset$RS) # Some NAs exist in this list, but should not be eliminated as
								# they correspond to real, although currently unidentified, SNPs. 
	} else {
		totalData <- as.numeric(fullset)
	}

	df1 <- data.frame("first"=c(sigData,randData))
	df2 <- data.frame("second"=rep(totalData,2))
	df2 <- df2 - df1
	res <- cbind(df1,df2)

	colnames(res) <- c(testGroup,paste0("non_",testGroup))
	rownames(res) <- c(paste("Significant",fullGroup), paste("Random",fullGroup))
	return(res)
}

Chi2 <- function (Data) {

	Chisq <- chisq.test(Data)
	n <- sum(Data)
	gl <- min(dim(Data) - 1)

	VCramer <- sqrt((Chisq$statistic)/(n * gl)) # Effect size code adapted from https://stats.stackexchange.com/questions/427864/how-to-calculate-an-effect-size-for-chi-square-in-r
	result <- data.frame("X2" = Chisq$statistic, "p-value" = Chisq$p.value, "Cramer_V" = VCramer, row.names = NULL)
	return(result)
}

# This function builds a table and checks for missing categories. For example,
# if any chromosomes are missing in the dataset, empty rows will be created
# for them.

tableGroups <- function(set, attr) { 

	table <- table(set[,attr])
	df <- data.frame(rep(0, tail(names(table), 1)))

	for (n in 1:length(table)) {
		df[names(table[n]),] <- table[n]
	}

	df <- cbind(seq(1:nrow(df)), df)
	return(df)
}

groupSigVsRandBarplot <- function(set1, set2, attr, title, fontsize = 28) {

	SigData <- tableGroups(set=set1, attr=attr)
	SigData <- cbind(SigData,rep("Significant",length(SigData)))
	RandData <- tableGroups(set=set2, attr=attr)
	RandData <- cbind(RandData,rep("Random",length(RandData)))
	colnames(SigData) = colnames(RandData) <- c("Chr","Frequency","Group")

	df <- rbind(SigData, RandData)
	df$Frequency = as.numeric(df$Frequency)
	df$Chr = as.numeric(df$Chr)

	tiff(file= paste0(title,".tiff"), units="cm",
		width=40, height=30, res=300)

	plot <- ggplot(df, aes(x=Chr, y=Frequency, fill=Group)) + 
			scale_x_continuous(breaks = unique(df$Chr)) +
   			geom_bar(position="dodge", stat="identity") +
   			scale_fill_manual(values = c("#000000","#7da1c4")) +
   			ggtitle(title) + xlab("Chromosome") +
   			theme(plot.title = element_text(hjust = 0.5)) +
   			scale_y_continuous(labels = function(x) format(x, scientific = TRUE)) +
   			theme(
   				axis.title = element_text(size = fontsize),
				legend.text = element_text(size = fontsize),
				plot.title = element_text(size = fontsize*1.2),
				axis.text = element_text(size = fontsize*0.8),
				legend.title = element_blank()
				)

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

message('\nBuilding bar plots\n')
groupSigVsRandBarplot(SignData$SNP, RandData$SNP, "CHR", "SNP distribution per chromosome")
groupSigVsRandBarplot(SignData$eQTL, RandData$eQTL, "CHR", "eQTL distribution per chromosome")
message('Building Volcano plots\n')
VolcanoPlot(x=SignData$SNP$beta, y=SignData$SNP$P, title="Significant_SNPs")
VolcanoPlot(x=RandData$SNP$beta, y=RandData$SNP$P, title="Random_SNPs")
# eQTL Chi squared analysis
message('Performing eQTL chi-squared tests\n')
setwd("../tables")
eQTL <- tableForChisq(SignData$eQTL, RandData$eQTL, SignData$SNP,"eQTL","SNP")
eQTLChisqres <- Chi2(eQTL)

write.table(eQTL, file="eQTL_Chi2.txt", sep ="\t", row.names = TRUE, col.names = TRUE)
cat("Chi squared value:\n", file="eQTL_Chi2.txt", append = TRUE)
cat(eQTLChisqres$X2, file="eQTL_Chi2.txt", append = TRUE)
cat("\np-value:\n", file="eQTL_Chi2.txt", append = TRUE)
cat(eQTLChisqres$p.value, file="eQTL_Chi2.txt", append = TRUE)
cat("\nMagnitude of effect (Cramer-V):\n", file="eQTL_Chi2.txt", append = TRUE)
cat(eQTLChisqres$Cramer_V, file="eQTL_Chi2.txt", append = TRUE)

message("Statistical analysis done!")
