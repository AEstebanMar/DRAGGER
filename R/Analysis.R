
### FUNCTIONS


DAGGER <- function(GWAS, GTEx, DGIdb) {
	merged <- merge_gene_var_drug(GWAS = GWAS, GTEx = GTEx, DGIdb = DGIdb)
	DAGGER_df <- predict_effect(merged)
	res <- get_candidates(DAGGER_df)
	return(res)
}
