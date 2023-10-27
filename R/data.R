#' GWAS demo dataset
#'
#' A very small Alzheimer's Diseasevariant risk association taken from The GWAS
#' Catalog via the gwasrapidd package
#'
#' @format ## `GWAS_demo`
#' A data frame with 6 rows and 3 columns:
#' \describe{
#'	\item{rs_id}{Reference SNP ID}
#'	\item{pvalue}{Statistical significance of SNP association to AD risk}
#'	\item{beta_number}{Beta coefficient of SNP association to AD. Positive}
#' }
#' @source <https://www.ebi.ac.uk/gwas/>
"GWAS_demo"

#' GTEx dataset
#'
#' Merged contents GTEx_Analysis_v8_eQTL.tar, shows variant-gene expression
#' association. Processing detailed under data-raw/, files GTEx.sh and GTEx.R
#'
#' @format ## `GTEx`
#' A data frame with 1208024 rows and 9 columns:
#' \describe{
#'	\item{Tissue}{Tissue where gene expression is altered}
#'	\item{gene_id}{Ensembl ID of affected gene}
#'  \item{gene_name}{Official gene symbol}
#'	\item{variant_id}{Unique variant identifier. Contains chromosome, position,
#'	ref and alt alleles and genome build}
#'	\item{tss_distance}{Distance to gene transcription start site}
#'	\item{rs_id_dbSNP151_GRCh38p7}{Reference SNP ID}
#'	\item{slope}{Effect on gene expression}
#'	\item{qval}{statistical significance of gene expression change}
#' }
#' @source <https://www.gtexportal.org/home/downloads/adult-gtex#qtl>
"GTEx"