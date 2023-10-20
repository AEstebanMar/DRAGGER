## R code to prepare `GTEx` dataset

# Merged_eQTL.txt comes from bash processing of GWAS dataset downloaded from
# https://www.gtexportal.org/home/downloads/adult-gtex#qtl. File GTEx_Analysis_v8_eQTL.tar.
# Pre-processing detailed in GTEx.sh

GTEx <- read.table ("Merged_eQTL.txt", header=TRUE, sep="", dec=".", fill = TRUE)
GTEx <- cbind(stringr::str_split_fixed(GTEx$file_name.gene_id, ".v8.egenes.txt,", 2), GTEx[, -1])
colnames(GTEx)[1:2] <- c("Tissue","gene_id")

usethis::use_data(GTEx, overwrite = TRUE)
