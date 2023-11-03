## code to prepare `DGIdb` dataset goes here
## bash
## wget https://www.dgidb.org/data/monthly_tsvs/2022-Feb/interactions.tsv

DGIdb <- read.table('interactions.tsv', header = TRUE,sep = "\t",
					quote = "\"", dec = ".", fill = TRUE,
					comment.char = "")[, c("gene_name", "entrez_id",
								"interaction_claim_source", "interaction_types",
								"drug_claim_name", "drug_concept_id",
								"interaction_group_score", "PMIDs")]

DGIdb[DGIdb==""] <- "Missing"

usethis::use_data(DGIdb, overwrite = TRUE)
