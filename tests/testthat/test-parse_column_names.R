test_that("exact matches work", {
  test_df <- data.frame(p_value = integer(0), beta_value = integer(0),
                rs_id = integer(0), gene_symbol = integer(0),
                drug_name = integer(0))
  parsed_colnames <- colnames(parse_column_names(test_df))
  expected_colnames <- c("p_value", "beta_value", "rs_id", "gene_symbol",
                          "drug_name")
  expect_equal(parsed_colnames, expected_colnames)
})

test_that("weirder, but real-world cases also match", {
  test_df <- data.frame(P = integer(0), beta = integer(0), SNP = integer(0),
                      gene_names = integer(0),
                      drug_claim_primary_name = integer(0))
  parsed_colnames <- colnames(parse_column_names(test_df))
  expected_colnames <- c("p_value", "beta_value", "rs_id", "gene_symbol",
                          "drug_name")
  expect_equal(parsed_colnames, expected_colnames)
})

test_that("multiple matches to p_value column return error", {
  test_df <- data.frame(p = integer(0), p_value = integer(0))
  expect_error(parse_column_names(test_df, "exceed max allowed value"))
})

test_that("if multiple matches exist, throw warning and choose first", {
  test_df <- data.frame(gene_name = integer(0), gene_name_bad = integer(0))
  parsed_colnames <- colnames(suppressWarnings(parse_column_names(test_df)))
  expected_colnames <- c("gene_symbol", "gene_name_bad")
  expect_warning(parse_column_names(test_df), "Multiple matches")
  expect_equal(parsed_colnames, expected_colnames)
})

test_that("columns with no matches are left unchanged", {
    test_df <- data.frame(foo = integer(0), bar = integer(0))
    expect_equal(parse_column_names(test_df), test_df)
  })

test_that("contents of data frame are never changed", {
    test_df <- data.frame(p_value = 1, foo = "bar", gene_name = "TST1",
                          gene_name_bad = "BAD2", NotNumber = NaN,
                          NotApplicable = NA)
    parsed_df <- suppressWarnings(parse_column_names(test_df))
    parsed_values <- unlist(parsed_df, use.names = FALSE)
    expected_values <- c("1", "bar", "TST1", "BAD2", "NaN", NA)
    expect_equal(parsed_values, expected_values)
  })
