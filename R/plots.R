
#' Volcano plot of variant p-value versus odds ratio
#' 
#' `plot_volcano` plots log2 of variant odds ratio versus log10 of p-value
#' for a given DRAGGER-parsed data frame of variant data.
#' @param df A DRAGGER-parsed data frame containing p-value and beta number
#' columns.
#' @param pval_col Name of p-value column. Default is "p_val_variant".
#' @param title Title of the resulting plot.
#' @param or_cutoff Value from which statistical significance line position
#' X axis (odds ratio) will be calculated. Symmetrical line will be
#' drawn in opposite value in X-axis (log2 of distribution between 0 and 1,
#' intended use for function, is symmetrical).
#' @param pval_cutoff Value from which Y axis (statistical
#' significance) line position will be calculated.
#' @param fontsize Base font size for text. Titles will be bigger, axes
#' will be lower.
#' @returns A volcano plot of the provided distributions.
#' @examples
#' GWAS_example <- GWAS_demo
#' colnames(GWAS_example) <- c("rs_id", "p_val_variant", "beta_number")
#' plot_volcano(GWAS_example)
#' @export

plot_volcano <- function(df,
						 title = "Odds Ratio vs Variant p-value",
						 pval_col = "p_val_variant",
						 or_cutoff=1.05, pval_cutoff=0.001, fontsize=32) {
	plot_df <- data.frame(cbind(df$beta_number, df[, pval_col]))
	colnames(plot_df) <- c("beta_number", "p_value")
	plot <- ggplot2::ggplot(data=plot_df,
						ggplot2::aes(x=log2(10**beta_number),
									 y=-log10(p_value))) +
			ggplot2::geom_point() +
			ggplot2::theme_minimal() +
			ggplot2::geom_vline(xintercept=c(-log2(or_cutoff),
											 log2(or_cutoff)),
								col="red") +
			ggplot2::geom_hline(yintercept=-log10(pval_cutoff), col="red") +
			ggplot2::ggtitle(title) +
			ggplot2::theme(
				plot.title = ggplot2::element_text(size=fontsize*1.2,
													hjust = 0.5),
				plot.subtitle = ggplot2::element_text(hjust = 0.5),
   				axis.title = ggplot2::element_text(size = fontsize),
				legend.text = ggplot2::element_text(size = fontsize),
				axis.text = ggplot2::element_text(size = fontsize*0.8),
				legend.title = ggplot2::element_blank()
				)
	return(plot)
}

#' Draw a bar plot of counts by group
#' 
#' `barplot_by_groups` takes a data frame, a value to count, and a category
#' by which to group the counts. Default is counting number of variants across
#' groups.
#' @param df A data frame with at least two columns.
#' @param value Value to count.
#' @param group.by Column by which counts will be grouped.
#' @returns A bar plot representing value counts across groups.
#' @examples
#' barplot_example_df <- data.frame(rs_id=rep('rs7357', 10),
#'									relevance=c(rep('important', 2),
#'											rep('not_important', 8))
#'									)
#' barplot_by_groups(barplot_example_df, "rs_id", "relevance")
#' @export

barplot_by_groups <- function(df, value = "rs_id", group.by) {
	subset <- df[, c(value, group.by)]
	colnames(subset) <- c("value", "group")
	subset_table <- plyr::count(subset, "group")
	plot <- ggplot2::ggplot(subset_table, ggplot2::aes(x = group, y = freq)) +
			ggplot2::geom_bar(stat = "identity") +
			ggplot2::scale_fill_grey() +
   			ggplot2::ggtitle(paste("RS count by", group.by)) +
   			ggplot2::xlab(group.by) +
   			ggplot2::ylab('count') +
   			ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
   							axis.text.x = ggplot2::element_text(angle = 90,
   														vjust = 0.5, hjust=1))
	return(plot)
}
