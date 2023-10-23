# Structural analysis

enrich_variants <- function(rs_list) {
	rs_vec <- as.vector(rs_list)[[1]]
	res <- gprofiler2::gost(c(list))$result
	return(res)
}

locate_variant_genes <- function(rs_list) {
	rs_vec <- as.vector(rs_list)[[1]]
	res_list <- sapply(rs_vec, function(x) gprofiler2::gsnpense(query = x, filter_na = FALSE))
	res_df <- .multi_join(res_list)
	return(res_df)
}

.multi_join <- function(list) {
	joined_df <- Reduce(function(...) merge(..., all=T), list)
	return(joined_df)
}
