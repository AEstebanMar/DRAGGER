test_that("DRAGGER basic pipeline works", {
  test_GWAS <- data.frame(snp = c("rs01", "rs02", "rs03", "rs04"),
                          beta_number = c(1, -1, 1, -1),
                          p_value = rep(2e-7, 4))
  test_GTEx <- data.frame(rs = c("rs01", "rs02", "rs03", "rs04"),
                          gene_symbol = c("gene01", "gene02", "gene03", "gene04"),
                          p_value = rep(1e-08, 4),
                          slope = c(1, 1, -1, -1))
  test_DGIdb <- data.frame(drug_name = c("drug01", "drug02",
                                         "drug03", "drug04"),
                           interaction_types = c("activator", "agonist",
                                                 "antagonist", "inhibitor"),
                           gene_symbol = c("gene01", "gene02", "gene03", "gene04"))
  test_results <- suppressWarnings(DRAGGER(test_GWAS, test_GTEx, test_DGIdb))
  expected_results <- data.frame(gene_symbol = c("gene01", "gene02",
                                                 "gene03", "gene04"),
                                 rs_id = c("rs01", "rs02", "rs03", "rs04"),
                                 beta_number = c(1, -1, 1, -1),
                                 p_val_variant = rep(2e-7, 4), 
                                 p_val_nominal = rep(1e-08, 4),
                                 slope = c(1, 1, -1, -1),
                                 drug_name = c("drug01", "drug02",
                                               "drug03", "drug04"),
                                 interaction_types = c("activator", "agonist",
                                                    "antagonist", "inhibitor"),
                                 prediction = c("inhibitor", "activator",
                                                "activator", "inhibitor"),
                                 candidate = c(FALSE, TRUE, FALSE, TRUE))
  expect_equal(test_results, expected_results)
})
