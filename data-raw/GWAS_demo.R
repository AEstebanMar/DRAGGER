## code to prepare `GWAS_demo` dataset

GWAS_demo <- gwasrapidd::get_associations(study_id = "GCST90245848")

associations <- data.frame(GWAS_demo@associations)
dec_betas <- which(associations$beta_direction=="decrease")
associations$beta_number[dec_betas] <- associations$beta_number[dec_betas] * -1

risk_alleles <- data.frame(GWAS_demo@risk_alleles)

GWAS_demo <- merge(risk_alleles, associations,
	by = "association_id")[, c("variant_id", "pvalue", "beta_number")]

usethis::use_data(GWAS_demo, overwrite = TRUE)
