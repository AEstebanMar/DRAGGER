
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
#' @importFrom stats chisq.test
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
