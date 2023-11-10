
#' Predict beneficial drug effect with variant risk and expression effect
#' 
#' `predict_effect` compares beta and slope columns in DAGGER dataframe
#' (containing variant and gene expression information) and predicts
#' beneficial drug effect.
#' @param gene_variant_df a dataframe resulted from merging a DAGGER-parsed
#' variant-risk dataframe and a DAGGER-parsed gene-variant dataframe.
#' @returns Input dataframe with an additional column with the predicted
#' beneficial drug effect according to beta and slope columns.
#' @examples
#' example_df <- merge_gene_var_drug(GWAS_demo, GTEx, DGIdb)[10:20, ]
#' example_df <- predict_effect(example_df)
#' @section Logic behind the function:
#' * Negative betas mean protection
#' * Positive betas mean risk.
#' * Negative slopes mean lower expression
#' * Positive slopes mean higher expression.
#' * If signs are opposite (prediction == TRUE), either the risk 
#' variant decreases expression or the protective variant increases it,
#' therefore an activator could be beneficial.
#' * If signs are equal (prediction == FALSE), either the protective variant
#' decreases expression or the risk variant increases it. Either way,
#' an inhibitor would be appropriate.
#' @export

predict_effect <- function(gene_variant_df) {

	message('Predicting beneficial drug effect')
	betas <- gene_variant_df$beta_number < 0
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

#' Filter potential candidates for drug repositioning
#' 
#' `get_candidates` compares predicted beneficial drug effect with described
#' gene-drug interaction. If they match, row is marked as potential candidate.
#' @param DAGGER_df A complete DAGGER dataset, containing risk-gene-drug
#' associations and predicted beneficial effect.
#' @returns Input dataset with an additional column recommending candidates
#' for drug repositioning.
#' @examples
#' example_df <- parse_column_names(head(GWAS_demo))
#' example_df <- filter_significance(example_df, 1e-20)
#' @export

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
	DAGGER_df$candidate <- FALSE
	DAGGER_df$candidate[candidates] <- TRUE
	return(DAGGER_df)
}
