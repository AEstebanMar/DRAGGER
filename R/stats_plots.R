
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
			ggplot2::geom_vline(xintercept=c(-log2(or_cutoff), log2(or_cutoff)),
								col="red") +
			ggplot2::geom_hline(yintercept=-log10(pval_cutoff), col="red") +
			ggplot2::ggtitle(title) +
			ggplot2::theme(
				plot.title = ggplot2::element_text(size=fontsize*1.2,
													hjust = 0.5),
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

	subs_list <- .get_subsets(df, col1, col1_val1, col1_val2)
	subs_list <- lapply(subs_list,
						function(x) .get_subsets(x, col2, col2_val1, col2_val2))

	total_rows <- lapply(unlist(subs_list, recursive=FALSE), nrow)
	res <- data.frame(c(total_rows[[1]], total_rows[[2]]),
						c(total_rows[[3]], total_rows[[4]]))
	colnames(res) <- c(paste(col2, col2_val1, sep = "_"),
						paste(col2, col2_val2, sep = "_"))
	rownames(res) <- c(paste(col1, col1_val1, sep = "_"),
						paste(col1, col1_val2, sep = "_"))
	return(res)
}

.get_subsets <- function (df, col, val_1, val_2) {
	expr_1 <- paste0(".*", val_1, ".*")
	expr_2 <- paste0(".*", val_2, ".*")
	subset_1 <- df[grep(expr_1, df[, col], ignore.case = TRUE), ]
	subset_2 <- df[grep(expr_2, df[, col], ignore.case = TRUE), ]
	res <- list(subset_1, subset_2)
	names(res) <- c(val_1, val_2)
	return(res)
}

#' Perform chi squared test on dataframe
#' 
#' `test_chi2` takes a data frame, two column names, and two values for each
#' column and formats the data into a chi2-test compatible format. Then,
#' it performs the chi2 statistical tests, including effect size (Cramer-V)
#' calculation.
#' @param df A data frame with at least two columns.
#' @param col1,col2 Columns containing categories to test.
#' @param col1_val1,col1_val2,col2_val1,col2_val2 Values of the first and
#' second columns, respectively, to test.
#' @returns A list containing the chi2 table and test results, with chi2 value,
#' p-value and effect size (Cramer-V)
#' @examples
#' chisq_example <- data.frame(Tissue=c("Brain", "Brain", "Stomach"),
#' 								candidate=c(TRUE, TRUE, FALSE))
#' chisq_example <- chisq_example[rep(seq_len(nrow(chisq_example)),
#'									each = 20), ]
#' test_chi2(chisq_example, "Tissue", "brain", "stomach", "candidate", T, F)
#' @export

test_chi2 <- function (df, col1, col1_val1, col1_val2,
							col2, col2_val1, col2_val2) {
	chisq_table <- .make_chisq_table(df, col1, col1_val1, col1_val2,
										col2, col2_val1, col2_val2)
	n <- sum(chisq_table)
	degr_freed <- min(dim(chisq_table) - 1)
	
	chisq_results <- chisq.test(chisq_table)
	Cramer_V <- sqrt((chisq_results$statistic)/ (n * degr_freed))
	result <- data.frame("X2" = chisq_results$statistic,
						"p-value" = chisq_results$p.value,
						"Cramer_V" = Cramer_V, row.names = NULL)
	result <- list(table=chisq_table, chi2=result)
	return(result)
}

#' Draw a bar plot of counts by group
#' 
#' `barplot_by_groups` takes a data frame, a value to count, and a category
#' by which to group the counts. Default is counting number of variants across
#' groups.
#' @param df A data frame with at least two columns.
#' @param value Value to count.
#' @param group.by Column by which counts will be grouped.
#' @returns A bar plot representing value counts across groups.
#' @examples
#' barplot_example_df <- data.frame(rs_id=rep('rs7357', 10),
#'									relevance=c(rep('important', 2),
#'												rep('not_important', 8)))
#' barplot_by_groups(barplot_example_df, "rs_id", "relevance")
#' test_chi2(chisq_example, "Tissue", "brain", "stomach", "candidate", T, F)
#' @export

barplot_by_groups <- function(df, value = "rs_id", group.by) {
	subset <- df[, c(value, group.by)]
	colnames(subset) <- c("value", "group")
	subset_table <- plyr::count(subset, "group")
	plot <- ggplot2::ggplot(subset_table, ggplot2::aes(x = group, y = freq)) +
			ggplot2::geom_bar(stat = "identity") +
			ggplot2::scale_fill_grey() +
   			ggplot2::ggtitle(paste("RS count by", group.by)) +
   			ggplot2::xlab(group.by) +
   			ggplot2::ylab('count') +
   			ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
   							axis.text.x = ggplot2::element_text(angle = 90,
   														vjust = 0.5, hjust=1))
	return(plot)
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
