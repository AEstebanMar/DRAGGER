#' Merge expression, GWAS and gene-drug interaction databases
#' 
#' `merge_gene_var_drug` takes three dataframes, one of GWAS, one of gene
#' expression and one of gene-drug interaction information and merges them.
#' It calls `remove_duplicate_rs`, `filter_significance` and
#' `parse_column_names` as needed. GWAS and GTEx dataframes are merged
#' by rs_id column, result is merged with DGIdb by gene_symbol column'
#' @param GWAS,GTEx,DGIdb GWAS, gene expression and drug-gene interaction
#' dataframes.
#' @returns Merged dataframe.
#' @examples
#' GWAS_example <- GWAS_demo[501:511, ]
#' GTEx_example <- GTEx[c(1:5, 1598, 146580, 219755), ]
#' DGIdb_example <- DGIdb[c((10:14), 5307, 5853, 8511), ]
#' merged <- merge_gene_var_drug(GWAS_example, GTEx_example, DGIdb_example)
#' @export

merge_gene_var_drug <- function(GWAS, GTEx, DGIdb) {

	message('\n\nParsing variant data\n\n')
	GWAS <- remove_duplicate_rs(
				filter_significance(parse_column_names(GWAS), 0.05)
				)
	GWAS_p_val_column <- grep("p_value", colnames(GWAS))
	if (length(GWAS_p_val_column) > 0) {
		colnames(GWAS)[GWAS_p_val_column] <- "p_val_variant"
	}
	message('\n\nParsing expression data\n\n')
	GTEx <- filter_significance(parse_column_names(GTEx), 0.05)
	GTEx_p_val_column <- grep("p_value", colnames(GTEx))
	if (length(GTEx_p_val_column) > 0) {
		colnames(GTEx)[GTEx_p_val_column] <- "p_val_nominal"
	}
	message('\n\nParsing drug-gene data\n\n')
	DGIdb <- parse_column_names(DGIdb)

	message('Merging genes and variants')
	gene_variants <- merge(GWAS, GTEx, by = "rs_id")

	message('Merging with drug database')
	res <- merge(gene_variants, DGIdb, by = "gene_symbol")
	return(res)
}