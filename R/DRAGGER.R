
#' Wrapper for default workflow
#' 
#' `DRAGGER` calls `merge_gene_var_drug`, `predict_effect` and `get_candidates`,
#' thus running the default workflow. This simple structure is the most
#' bare-bones case, but is effective nonetheless. Resulting object may be
#' further processed for more specific results.
#' @inheritParams merge_gene_var_drug
#' @returns A data frame containing full repositioning analysis. Keeps all
#' colums from input data.
#' @examples
#' GWAS_example <- GWAS_demo[501:511, ]
#' GTEx_example <- GTEx[c(1:5, 1598, 146580, 219755), ]
#' DGIdb_example <- DGIdb[c((10:14), 5307, 5853, 8511), ]
#' DRAGGER_result <- DRAGGER(GWAS_example, GTEx_example, DGIdb_example)
#' @export


DRAGGER <- function(GWAS, GTEx, DGIdb) {
	merged <- merge_gene_var_drug(GWAS = GWAS, GTEx = GTEx, DGIdb = DGIdb)
	DRAGGER_df <- predict_effect(merged)
	res <- get_candidates(DRAGGER_df)
	return(res)
}
