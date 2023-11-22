test_that("merge_gene_var_drug works as intended", {
  dummy_GWAS <- data.frame(rs_id = c("rs01", "rs02", "rs03", "rs04"),
                           p_value = rep(0, 4))
  dummy_GWAS$p_value[3] <- 1
  dummy_GTEx <- data.frame(rs_id = c("rs01", "rs02", "rs03"),
                           gene_symbol = c("gene1", "gene2", "gene3"),
                           p_value = rep(0, 3))
  dummy_DGIdb <- data.frame(gene_symbol = c("gene2", "gene3"),
                            drug_name = c("drug2", "drug3"))
  expected_df <- data.frame(gene_symbol = c("gene2"),
                            rs_id = c("rs02"),
                            p_val_variant = 0,
                            p_val_nominal = 0,
                            drug_name = c("drug2"))
  expect_equal(merge_gene_var_drug(dummy_GWAS, dummy_GTEx, dummy_DGIdb),
               expected_df)
})

test_that("merge_gene_var_drug returns nothing if no matches are found", {
  lonely_GWAS <- data.frame(rs_id = c("rslonely"), p_value = 0)
  lonely_GTEx <- data.frame(rs_id = c("rsalsolonely"),
                            gene_symbol = "lonelygene", p_value = 0)
  lonely_DGIdb <- data.frame(gene_symbol = "noneyouknow", drug_name = "lembas")
  expected_df <- data.frame(gene_symbol = character(0), rs_id = character(0),
                            p_val_variant = numeric(0),
                            p_val_nominal = numeric(0),
                            drug_name = character(0))
  expect_equal(merge_gene_var_drug(lonely_GWAS, lonely_GTEx, lonely_DGIdb),
               expected_df)
  })

test_that("merge_gene_var_drug with missing GWAS p-value", {
  test_GWAS <- data.frame(rs_id = "rs01",
                          beta_value = 1)
  test_GTEx <- data.frame(rs = "rs01",
                          gene_symbol = "gene01",
                          p_value = 1e-08,
                          slope = 1)
  test_DGIdb <- data.frame(drug_name = "drug01",
                           interaction_types = "activator",
                           gene_symbol = "gene01")
  expect_no_error(suppressWarnings(merge_gene_var_drug(test_GWAS,
                                                      test_GTEx, test_DGIdb)))
  })

test_that("merge_gene_var_drug with missing GTEx p-value", {
  test_GWAS <- data.frame(rs_id = "rs01",
                          p_value = 1e-08,
                          beta_value = 1)
  test_GTEx <- data.frame(rs = "rs01",
                          gene_symbol = "gene01",
                          slope = 1)
  test_DGIdb <- data.frame(drug_name = "drug01",
                           interaction_types = "activator",
                           gene_symbol = "gene01")
  expect_no_error(suppressWarnings(merge_gene_var_drug(test_GWAS,
                                                      test_GTEx, test_DGIdb)))
  })

test_that("merge_gene_var_drug with no p-values", {
    test_GWAS <- data.frame(rs_id = "rs01",
                          beta_value = 1)
  test_GTEx <- data.frame(rs = "rs01",
                          gene_symbol = "gene01",
                          slope = 1)
  test_DGIdb <- data.frame(drug_name = "drug01",
                           interaction_types = "activator",
                           gene_symbol = "gene01")
  expect_no_error(suppressWarnings(merge_gene_var_drug(test_GWAS,
                                                      test_GTEx, test_DGIdb)))
  })
