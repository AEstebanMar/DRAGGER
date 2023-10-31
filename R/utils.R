
normalise_columns <- function(df) {
	pval_column <- grep("pval|p-val", colnames(df), ignore.case=TRUE)
	rs_column <- grep("rs", colnames(df), ignore.case=TRUE)
	gene_id_column <- grep("rs", colnames(df), ignore.case=TRUE)
	drug_name_column <- grep("drug.*name", colnames(df), ignore.case=TRUE)
	if (length(pval_column) != 0) {
		message('P-value column found. Setting name to "pvalue"')
		colnames(df)[pval_column] <- "pvalue"
	}
	if (length(rs_column) != 0) {
		message('RS ID column found. Setting name to "rs_id"')
		colnames(df)[rs_column] <- "rs_id"
	} 
	if (length(gene_id_column) != 0) {
		message('Gene ID column found. Setting name to "gene_id"')
		colnames(df)[gene_id_column] <- "gene_id"
	}
	if (length(drug_name_column) > 1) {
		message('Interactions dataset has multiple drug name columns. Column
			with the most unique names has been chosen')
		drug_columns <- df[, drug_name_column]
		uniques <- lapply(X = drug_columns, FUN = function(x) length(unique(x)))
		longest <- which.max(unlist(uniques))
		colnames(df)[longest] <- "drug_name"
	}
	if (length(drug_name_column) == 1) {
		message('Drug name column found. Setting name to "drug_name"')
		colnames(df)[drug_name_column] <- "drug_name"
	}
	return(df)
}