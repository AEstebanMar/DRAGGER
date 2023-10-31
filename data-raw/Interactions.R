## code to prepare `Interactions` dataset goes here
## bash
## wget https://www.dgidb.org/data/monthly_tsvs/2022-Feb/interactions.tsv

Interactions <- read.table('interactions.tsv', header = TRUE,sep = "\t",
					quote = "\"", dec = ".", fill = TRUE,
					comment.char = "")[, c("gene_name", "entrez_id",
								"interaction_claim_source", "interaction_types",
								"drug_claim_name", "drug_concept_id",
								"interaction_group_score", "PMIDs")]

Interactions[Interactions==""] <- "Missing"

usethis::use_data(Interactions, overwrite = TRUE)
