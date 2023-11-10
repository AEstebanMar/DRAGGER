
#' Remove duplicate RS in GWAS dataframe
#' 
#' `remove_duplicate_rs` removes rows corresponding to duplicate rs, keeping
#' the first occurrence. Dataframe is first sorted by p-value if column exists.
#' @param df A DAGGER-parsed data frame of GWAS information. Must contain an
#' rs_id column. p_value column optional, but recommended.
#' @returns A data frame with no duplicate RS. Only the first occurrence of each
#' duplicate will remain, or the most statistically significant if p-value
#' column is found'
#' @examples
#' dupes_example <- data.frame(
#'						rs_id = c("rs2710888", "rs182532", "rs9660106"),
#' 						pvalue = c(2e-58, 1e-17, 2e-12),
#' 						beta_number = GWAS_demo$beta_number[c(1, 3, 5)]
#'														+ 0.002)
#' dupes_example <- rbind(head(GWAS_demo), dupes_example)
#' dupes_example <- remove_duplicate_rs(dupes_example)
#' @export

remove_duplicate_rs <- function(df) {

	if(!is.null(df$p_value)) {
		message("Sorting input by statistical significance")
		df <- df[order(df$p_value), ]
	}
	if(is.null(df$rs_id)) {
		stop('RS ID column not found or not properly parsed. Please run input
			through DAGGER::parse_column_names()')
	}

	if(!any(duplicated(df$rs_id))) {
		return(df)
	} else {
		message("Duplicates found in dataset, choosing first occurrence")
		res <- df[!duplicated(df$rs_id), ]
		return(res)
	}		
}

#' Filter dataframe by p-value
#' 
#' `filter_significance` filters out all rows in a data frame below input
#' p-value (default 0.05).
#' @param df A DAGGER-parsed data frame.
#' @param value P-value cutoff.
#' @returns A subset of the original data frame with all rows passing filter.
#' @examples
#' example_df <- head(GWAS_demo)
#' colnames(example_df)[2] <- "p_value"
#' example_df <- filter_significance(example_df, 1e-20)
#' @export

filter_significance <- function(df, value = 0.05) {
	if(is.null(df$p_value)) {
		message("No statistical significance column found in input or improperly
			parsed. You might want to run it through DAGGER::parse_column_names
			first. Returning it as-is")
		return(df)
	}
	res <- df[as.numeric(df$p_value) <= value, ]
	return(res)
}