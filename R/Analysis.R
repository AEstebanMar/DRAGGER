
### FUNCTIONS

predict_effect <- function(gene_variant_df) {
	# Logic: negative betas mean protection, positive betas mean risk.
	# Negative slopes mean lower expression, positive slopes mean higher 
	# expression. If signs are opposite (prediction == TRUE), either the risk 
	# variant decreases expression or the protective variant increases it,
	# therefore an activator could be beneficial. If signs are equal
	# (prediction == FALSE), either the protective variant decreases expression
	# or the risk variant increases it. Either way, an inhibitor is desired.
	message('Predicting beneficial drug effect')
	betas <- gene_variant_df$beta < 0
	slopes <- gene_variant_df$slope > 0
	prediction <- betas == slopes
	prediction[prediction == TRUE] <- "activator"
	prediction[prediction == FALSE] <- "inhibitor"
	res <- cbind(gene_variant_df, prediction)
	return(res)
}

.is.candidate <- function(DAGGER_df, dict) {
	broad_type <- rep(NULL, nrow(DAGGER_df))
	broad_type[DAGGER_df$interaction_types %in% dict$activator] <- "activator"
	broad_type[DAGGER_df$interaction_types %in% dict$inhibitor] <- "inhibitor"
	candidates <- DAGGER_df$prediction == broad_type
	return(candidates)
}

get_candidates <- function(DAGGER_df) {
	message("Producing list of candidates for repositioning")
	dict <- list(
		activator = c("agonist", "activator", "positive modulator",
			"partial agonist", "inducer", "allosteric modulator"),
		inhibitor = c("inhibitor", "blocker", "antagonist", "inverse agonist",
			"negative modulator", "antisense oligonucleotide", "suppressor",
			"inhibitory allosteric modulator")
		)
	candidates <- .is.candidate(DAGGER_df, dict)
	DAGGER_df$candidate <- NULL
	DAGGER_df$candidate[candidates] <- TRUE
	return(DAGGER_df)
}

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

DAGGER <- function(GWAS, GTEx, DGIdb) {
	merged <- merge_gene_var_drug(GWAS = GWAS, GTEx = GTEx, DGIdb = DGIdb)
	DAGGER_df <- predict_effect(merged)
	res <- get_candidates(DAGGER_df)
	return(res)
}
