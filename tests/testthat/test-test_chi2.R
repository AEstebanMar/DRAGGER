test_that("chi2 tables are properly formed", {
  test_df <- data.frame(kingdom=c("Daggerfall", "Wayrest", "Orsinium"),
                              race=c("Breton", "Breton", "Orc"))
  res <- .make_chisq_table(test_df, "kingdom", "Daggerfall", "Orsinium",
                                    "race", "Breton", "Orc")
  expected_df <- data.frame(race_Breton = 1:0, race_Orc = 0:1)
  rownames(expected_df) <- c("kingdom_Daggerfall", "kingdom_Orsinium")
  expect_equal(res, expected_df)
})

test_that("chi2 test is properly done", {
  test_df <- data.frame(matrix(ncol = 2, nrow = 1000))
  colnames(test_df) <- c("kingdom", "arrow")
  test_df[1:500, 1] <- "Wayrest"
  test_df[501:1000, 1] <- "Sentinel"
  test_df[1:400, 2] <- "Hit"
  test_df[401:500, 2] <- "Miss"
  test_df[501:701, 2] <- "Hit"
  test_df[701:1000, 2] <- "Miss"
  res <- test_chi2(test_df, "kingdom", "Wayrest", "Sentinel",
                                  "arrow", "Hit", "Miss")
  expected_table <- data.frame(arrow_Hit = c(400, 100),
                               arrow_Miss = c(200, 300))
  rownames(expected_table) <- c("kingdom_Wayrest", "kingdom_Sentinel")
  expected_chi2 <- data.frame(X2 = 165.0041667, p.value = 9.1285076e-38,
                              Cramer_V = 0.4062070490)
  expect_equal(res$table, expected_table)
  expect_equal(res$chi2, expected_chi2)
})
