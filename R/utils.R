
normalise_columns <- function(df) {
	pval_column <- grep("pval|p-val", colnames(df), ignore.case=TRUE)

	colnames(df)[pval_column] <- "pvalue"
	rs_column <- grep("rs", colnames(df), ignore.case=TRUE)
	colnames(df)[rs_column] <- "rs_id"
	gene_id_column <- grep("rs", colnames(df), ignore.case=TRUE)
	return(df)
}