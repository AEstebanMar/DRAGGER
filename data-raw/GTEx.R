## R code to prepare `GTEx` dataset

# Merged_eQTL.txt comes from bash processing of GWAS dataset downloaded from
# https://www.gtexportal.org/home/downloads/adult-gtex#qtl.
# File GTEx_Analysis_v8_eQTL.tar. Pre-processing detailed in GTEx.sh

GTEx <- read.table ("Merged_eQTL.txt", header=TRUE,
					sep="", dec=".", fill = TRUE)[, c("file_name.gene_id",
					"gene_name", "strand", "variant_id", "tss_distance",
					"rs_id_dbSNP151_GRCh38p7","pval_nominal","slope")]
GTEx$pval_nominal <- as.numeric(GTEx$pval_nominal)
colnames(GTEx)[4] <- "variant_information"


split_first_col <- do.call(rbind,
						strsplit(GTEx$file_name.gene_id, ".v8.egenes.txt,"))

GTEx <- cbind(split_first_col, GTEx[, -1])

colnames(GTEx)[1:2] <- c("Tissue","gene_id")

usethis::use_data(GTEx, overwrite = TRUE)
