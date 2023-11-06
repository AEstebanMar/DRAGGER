#' GWAS demo dataset
#'
#' A height Genome-Wide Association Study taken from The GWAS
#' Catalog using the gwasrapidd package
#'
#' @format ## `GWAS_demo`
#' A data frame with 6 rows and 3 columns:
#' \describe{
#'	\item{rs_id}{Reference SNP ID}
#'	\item{pvalue}{Statistical significance of SNP association to AD risk}
#'	\item{beta_number}{Beta coefficient of SNP association to AD. Positive}
#' }
#' @source <https://www.ebi.ac.uk/gwas/studies/GCST90245848>
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
#'  \item{gene_name}{GENCODE gene name}
#'	\item{variant_id}{Unique variant identifier. Contains chromosome, position,
#'	ref and alt alleles and genome build}
#'	\item{tss_distance}{Distance to gene transcription start site}
#'	\item{rs_id_dbSNP151_GRCh38p7}{Reference SNP ID}
#'	\item{slope}{Effect on gene expression}
#'	\item{pval_nominal}{statistical significance of gene-variant association}
#' }
#' @source <https://www.gtexportal.org/home/downloads/adult-gtex#qtl>
"GTEx"

#' DGIdb dataset
#'
#' interactions.tsv file downloaded from dgidb.org
#'
#' @format ## `DGIdb`
#' A data frame with 85460 rows and 11 columns:
#' \describe{
#'	\item{gene_name}{GENCODE gene name}
#'	\item{gene_claim_name}{Alternate gene name, if any}
#'	\item{entrez_id}{Gene Entrez identifier}
#'	\item{interaction_claim_source}{Source of claimed interaction}
#'	\item{interaction_types}{Claimed interaction}
#'	\item{drug_claim_name}{Drug name}
#'	\item{drug_concept_id}{Unique CHEMBL drug identifier}
#'	\item{interaction_group_score}{Statistical significance of interaction.
#'	Detailed at https://www.dgidb.org/score}
#'	\item{PMIDs}{}
#' }
#' @source <https://www.gtexportal.org/home/downloads/adult-gtex#qtl>
"DGIdb"