
.rename_column <- function(df, expr, new_name, max_matches = Inf) {
	message(paste('Matching', new_name, "column..."))
	columns <- grep(expr, colnames(df), ignore.case=TRUE)
	if (length(columns) == 0) {
		message('No matches found')
	}

	if (length(columns) == 1) {
		message('Match found. Parsing...')
		colnames(df)[columns] <- new_name
	}

	if (length(columns) > 1) {
		warning('Multiple matches found. Choosing column with the most
			unique names', immediate. = TRUE)
		uniques <- lapply(X = df[columns], FUN = function(x) length(unique(x)))
		longest <- which.max(unlist(uniques))
		colnames(df)[columns[longest]] <- new_name
	}

	if (length(columns) > max_matches) {
		error('Matches exceed max allowed value for column type. Please
			choose one and remove or rename the others (e.g por expression data,
			choose p-value for strengh of variant-expression association')
	}

	return(df)
}

#' Parse column names to package standard
#' 
#' `parse_column_names` takes an input data frame and renames set columns in
#' order to make them recognisable by other DAGGER functions.
#' @param df A data frame. Must contain at least one of the columns expected
#' for analysis (RS ID, gene ID, drug name). If multiple columns exist for
#' any of this information, function will choose the one with most unique names.
#' Can also recognise, but does not require, a p-value column. Please ensure
#' only one of those columns exists in data frame, or rename those not wanted
#' for analysis to anything not matching "p-val" or "pval".
#' @returns A data frame of same dimensions as `df` with RS ID, gene ID,
#' drug name and p-value (if found) named to DAGGER standards "rs_id",
#' "gene_id", "drug_name", "pvalue". Functions further down the pipeline expect
#' columns with said names. A warning will be given if multiple matches are
#' found for any column other than p-value, which returns an error.
#' @examples
#' data(GWAS_demo)
#' print(GWAS_demo)
#' parse_column_names(GWAS_demo)
#' data(GTEx)
#' GTEx_demo <- head(GTEx)
#' print(GTEx_demo)
#' parse_column_names(GTEx_demo)

parse_column_names <- function(df) {
	df <- .rename_column(df, "pval|p-val", "p_value", 1)
	df <- .rename_column(df, "rs", "rs_id")
	df <- .rename_column(df, "gene.*id", "gene_id")
	df <- .rename_column(df, "drug.*name", "drug_name")
	return(df)
}	
