
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
		message("WARNING")
		warning('Multiple matches found. Choosing column with the most
			unique names for standardising')
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

parse_column_names <- function(df) {
	df <- .rename_column(df, "pval|p-val", "pvalue", 1)
	df <- .rename_column(df, "rs", "rs_id")
	df <- .rename_column(df, "gene.*id", "gene_id")
	df <- .rename_column(df, "drug.*name", "drug_name")
	return(df)
}	
