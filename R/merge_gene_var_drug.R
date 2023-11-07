merge_gene_var_drug <- function(GWAS, GTEx, DGIdb) {

	message('Parsing variant data')
	GWAS <- remove_duplicate_rs(
				filter_significance(parse_column_names(GWAS), 0.05)
				)
	message('Parsing expression data')
	GTEx <- filter_significance(parse_column_names(GTEx), 0.05)
	message('Parsing drug data')
	DGIdb <- parse_column_names(DGIdb)

	message('Merging genes and variants')
	gene_variants <- merge(GWAS, GTEx, by = "rs_id")
	message('Merging with drug database')
	res <- merge(gene_variants, DGIdb, by = "gene_symbol")
	return(res)
}