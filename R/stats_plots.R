
#' Volcano plot of variant p-value versus odds ratio
#' 
#' `plot_volcano` plots log2 of variant odds ratio versus log10 of p-value
#' for a given DAGGER-parsed dataframe containing variant data.
#' @param df A DAGGER-parsed data frame containing p-value and beta number
#' columns.
#' @param title Title of the resulting plot.
#' @param or_cutoff Value from which statistical significance line position for
#' odds ratio (X axis) will be calculated. Symmetrical line will be drawn in
#' opposite value in X-axis (log2 of odds-ratio distribution is symmetrical).
#' @param pval_cutoff Value from which statistical significance line position
#' for p-value (Y axis) will be calculated.
#' @returns A vulcano plot of the provided distributions.
#' @examples
#' GWAS_example <- GWAS_demo
#' colnames(GWAS_example) <- c("rs_id", "p_value", "beta_number")
#' plot_volcano(GWAS_example)
#' @export

plot_volcano <- function (df, title = "Odds Ratio vs Variant p-value",
							or_cutoff=1.05, pval_cutoff=0.001, fontsize=32) {

	df <- parse_column_names(df)
	plot <- ggplot2::ggplot(data=df, ggplot2::aes(x=log2(10**beta_number),
													y=-log10(p_value))) +
			ggplot2::geom_point() +
			ggplot2::theme_minimal() +
			ggplot2::geom_vline(xintercept=c(-log2(or_cutoff),
												log2(or_cutoff)),
								col="red") +
			ggplot2::geom_hline(yintercept=-log10(pval_cutoff), col="red") +
			ggplot2::ggtitle(title) +
			ggplot2::theme(
				plot.title = ggplot2::element_text(size=fontsize*1.2, hjust = 0.5),
				plot.subtitle = ggplot2::element_text(hjust = 0.5),
   				axis.title = ggplot2::element_text(size = fontsize),
				legend.text = ggplot2::element_text(size = fontsize),
				axis.text = ggplot2::element_text(size = fontsize*0.8),
				legend.title = ggplot2::element_blank()
				)
	return(plot)
}

.make_chisq_table <- function(df, col1, col1_val1, col1_val2,
								col2, col2_val1, col2_val2) {

	set1 <- .subset_df(df, col1, col1_val1)
	set1_1 <- .subset_df(set1, col2, col2_val1)
	set1_2 <- .subset_df(set1, col2, col2_val2)
	set2 <- .subset_df(df, col1, col1_val2)
	set2_1 <- .subset_df(set2, col2, col2_val1)
	set2_2 <- .subset_df(set2, col2, col2_val2)

	N_1 <- nrow(set1)
	N_2 <- nrow(set2)
	N_total <- nrow(fullset)

	res <- data.frame(N_1, N_2)
	res[2, ] <- res[2, ] - res[1, ]

	colnames(res) <- c(val1, val2)
	rownames(res) <- c(paste0("Not_",val1), paste0("Not_",val2))
	return(res)
}

.subset_df <- function (df, col, val) {
	expr <- paste0(".*", val, ".*")
	res <- df[grep(expr, df[, col], ignore.case = TRUE), ]
	return(res)
}

#' Perform chi squared test on dataframe
#' 
#' `test_chi2` plots log2 of variant odds ratio versus log10 of p-value
#' for a given DAGGER-parsed dataframe containing variant data.
#' @param df A DAGGER-parsed data frame containing p-value and beta number
#' columns.
#' @param title Title of the resulting plot.
#' @param or_cutoff Value from which statistical significance line position for
#' odds ratio (X axis) will be calculated. Symmetrical line will be drawn in
#' opposite value in X-axis (log2 of odds-ratio distribution is symmetrical).
#' @param pval_cutoff Value from which statistical significance line position
#' for p-value (Y axis) will be calculated.
#' @returns A vulcano plot of the provided distributions.
#' @examples
#' GWAS_example <- head(GWAS_demo)
#' colnames(GWAS_example)[1:2] <- c("rs_id", "p_value")
#' VolcanoPlot(GWAS_example)
#' @export

test_chi2 <- function (df, col, val1, val2) {

	chisq_table <- .make_chisq_table(df, col, val1, val2)

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

tableGroups <- function(set, col) { 

	table <- table(set[,col])
	df <- data.frame(rep(0, tail(names(table), 1)))

	for (n in 1:length(table)) {
		df[names(table[n]),] <- table[n]
	}

	df <- cbind(seq(1:nrow(df)), df)
	return(df)
}

groupSigVsRandBarplot <- function(set1, set2, col, title, fontsize = 28) {

	SigData <- tableGroups(set=set1, col=col)
	SigData <- cbind(SigData,rep("Significant",length(SigData)))
	RandData <- tableGroups(set=set2, col=col)
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

# message('============================================================================================\nStarting statistical analyses...\n\nLoading input data...')
# library(ggplot2)
# setwd('../data/processed')
# SignData <- list("SNP"=readRDS("Sign.rds"),"eQTL"=readRDS("SigMerged.rds"))
# RandData <- list("SNP"=readRDS("Rand.rds"),"eQTL"=readRDS("RandMerged.rds"))
# TopQval <- readRDS("TopQval.rds")
# setwd('../../output/plots')

# #### Build data frames and plots

# message('\nBuilding bar plots\n')
# groupSigVsRandBarplot(SignData$SNP, RandData$SNP, "CHR", "SNP distribution per chromosome")
# groupSigVsRandBarplot(SignData$eQTL, RandData$eQTL, "CHR", "eQTL distribution per chromosome")
# message('Building Volcano plots\n')
# VolcanoPlot(x=SignData$SNP$beta, y=SignData$SNP$P, title="Significant_SNPs")
# VolcanoPlot(x=RandData$SNP$beta, y=RandData$SNP$P, title="Random_SNPs")
# # eQTL Chi squared analysis
# message('Performing eQTL chi-squared tests\n')
# setwd("../tables")
# eQTL <- tableForChisq(SignData$eQTL, RandData$eQTL, SignData$SNP,"eQTL","SNP")
# eQTLChisqres <- Chi2(eQTL)

# write.table(eQTL, file="eQTL_Chi2.txt", sep ="\t", row.names = TRUE, col.names = TRUE)
# cat("Chi squared value:\n", file="eQTL_Chi2.txt", append = TRUE)
# cat(eQTLChisqres$X2, file="eQTL_Chi2.txt", append = TRUE)
# cat("\np-value:\n", file="eQTL_Chi2.txt", append = TRUE)
# cat(eQTLChisqres$p.value, file="eQTL_Chi2.txt", append = TRUE)
# cat("\nMagnitude of effect (Cramer-V):\n", file="eQTL_Chi2.txt", append = TRUE)
# cat(eQTLChisqres$Cramer_V, file="eQTL_Chi2.txt", append = TRUE)

# message("Statistical analysis done!")
