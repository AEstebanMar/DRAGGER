test_that("predict_effect works as intended", {
  test_df <- data.frame(beta_number = c(-1, -1, 1, 1), slope = c(1, -1, 1, -1))
  res <- predict_effect(test_df)
  expected_prediction <- c("activator", "inhibitor", "inhibitor", "activator")
  expect_equal(res$prediction, expected_prediction)
})

test_that("predict_effect returns NA for unexpected values", {
  beta_test_df <- data.frame(beta_number = c(0, NA, NaN, "Mantella"),
                        slope = c(1, -1, 1, -1))
  slope_test_df <- data.frame(beta_number = c(-1, -1, 1, 1),
                              slope = c(0, NA, NaN, "Totem"))
  expect_warning(predict_effect(beta_test_df), "beta number column")
  expect_warning(predict_effect(slope_test_df), "slope column")
  })
