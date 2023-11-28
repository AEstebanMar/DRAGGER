test_that("remove_duplicate_rs fails if no rs_id column exists", {
  test_df <- data.frame(kingdoms=c("Daggerfall", "Sentinel",
                                    "Wayrest", "Orsinium"))
  expect_error(remove_duplicate_rs(test_df), "column not found")
})

test_that("remove_duplicate_rs sorts by statistical significance if possible", {
  test_df <- data.frame(rs_id = c("rs4", "rs3", "rs2", "rs1"),
                        p_val_variant = 4:1)
  test_df <- print(test_df, row.names=FALSE)
  expected_df <- data.frame(rs_id = c ("rs1", "rs2", "rs3", "rs4"),
                            p_val_variant = 1:4)
  rownames(expected_df) <- 4:1
  expect_equal(remove_duplicate_rs(test_df), expected_df)
  })

test_that("remove_duplicate_rs works as intended and chooses most significant one", {
  test_df <- data.frame(rs_id = c("rs1", "rs2", "rs1", "rs3"), p_val_variant = 4:1)
  expected_df <- data.frame(rs_id = c("rs3", "rs1", "rs2"), p_val_variant = 1:3)
  rownames(expected_df) <- 4:2
  expect_equal(remove_duplicate_rs(test_df), expected_df)
  })

test_that("filter_significance filters correctly by p-value", {
  test_df <- data.frame(foo = rep("bar", 50), p_val_variant = 1:50)
  expected_df <- data.frame(foo = rep("bar", 25), p_val_variant = 1:25)
  expect_equal(filter_significance(test_df, 25), expected_df)
  })

test_that("filter_significance fails if no p-value column exists", {
  test_df <- data.frame(sorceress = "Nulfaga")
  expect_warning(filter_significance(test_df), "No statistical significance")
  })

test_that("filter_significance fails if multiple p-value columns exist", {
  test_df <- data.frame(p_val_variant = 1, p_val_nominal = 2)
  expect_warning(filter_significance(test_df), "Multiple p-value")
  })