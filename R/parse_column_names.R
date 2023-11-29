
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
		stop('Matches exceed max allowed value for column type. Please
			choose one and remove or rename the others (e.g por expression data,
			choose p-value for strengh of variant-expression association')
	}

	return(df)
}

#' Parse column names to package standard
#' 
#' `parse_column_names` takes an input data frame and renames set columns in
#' order to simplify manipulation by other DRAGGER functions.
#' @param df A data frame. Must contain at least one of the columns expected
#' for analysis (RS ID, gene ID, drug name). If multiple columns exist for
#' any of this information, function will choose the one with most unique names.
#' Can also recognise, but does not require, a p-value column. Please ensure
#' only one of those columns exists in data frame, or rename those not wanted
#' for analysis to anything not matching "p-val" or "pval".
#' @param type Type of information contained in df. Renames p-value column
#' differently for variant (type "GWAS") and expression (type "GTEx") data.
#' @returns Input dataframe with RS ID, gene ID, drug name and p-value
#' (if found) renamed to DRAGGER standards "rs_id", "gene_id", "drug_name",
#' "pvalue". A warning will be given if multiple matches are found for any
#' column other than p-value, which returns an error.
#' @examples
#' parse_column_names(head(GWAS_demo), "GWAS")
#' parse_column_names(head(GTEx), "GTEx")
#' @export

parse_column_names <- function(df, type) {
	message(paste0("Parsing data frame as ", type, " data"))
	if (type == "GWAS") {
		df <- .rename_column(df, "p.*val|^p$", "p_val_variant", 1)
		df <- .rename_column(df, "beta$", "beta_value", 1)
	} else if (type == "GTEx") {
		df <- .rename_column(df, "p.*val|^p$", "p_val_nominal", 1)
	}
	df <- .rename_column(df, "rs|^snp.*$|^snp.*id|^variant.*id", "rs_id")
	df <- .rename_column(df, "gene.*name|gene.*symbol", "gene_symbol")
	df <- .rename_column(df, "drug.*name", "drug_name")
	return(df)
}	
