test_that("predict_effect works as intended", {
  test_df <- data.frame(beta_value = c(-1, -1, 1, 1), slope = c(1, -1, 1, -1))
  res <- predict_effect(test_df)
  expected_prediction <- c("activator", "inhibitor", "inhibitor", "activator")
  expect_equal(res$prediction, expected_prediction)
})

test_that("predict_effect returns NA for unexpected values", {
  beta_test_df <- data.frame(beta_value = c(0, NA, NaN, "Mantella"),
                        slope = c(1, -1, 1, -1))
  slope_test_df <- data.frame(beta_value = c(-1, -1, 1, 1),
                              slope = c(0, NA, NaN, "Totem"))
  expect_warning(predict_effect(beta_test_df), "beta number column")
  expect_warning(predict_effect(slope_test_df), "slope column")
})

test_that("get_candidates works as expected", {
  test_activators <- data.frame(interaction_types = rep(c("agonist", "activator",
                        "positive modulator", "partial agonist", "inducer",
                        "allosteric modulator"), 2),
                        prediction = c(rep("activator", 6),
                                       rep("inhibitor", 6)))
  activator_res <- get_candidates(test_activators)
  expected_activator_candidates <- c(rep(TRUE, 6), rep(FALSE, 6))
  test_inhibitors <- data.frame(interaction_types = rep(c("inhibitor",
                                  "blocker", "antagonist", "inverse agonist",
                                  "negative modulator",
                                  "antisense oligonucleotide", "suppressor",
                                  "inhibitory allosteric modulator"), 2),
                                prediction = c(rep("activator", 8),
                                               rep("inhibitor", 8)))
  inhibitor_res <- get_candidates(test_inhibitors)
  expected_inhibitor_candidates <- c(rep(FALSE, 8), rep(TRUE, 8))
  expect_equal(activator_res$candidate, expected_activator_candidates)
  expect_equal(inhibitor_res$candidate, expected_inhibitor_candidates)
  })

test_that("get_candidates returns FALSE on interaction_types not in dict", {
  test_missing <- data.frame(interaction_types = c("Missing", "?", NA, NaN, 1,
                                                   0, -1, TRUE, FALSE),
                             prediction = c(rep("Activator", 5),
                                            rep("Inhibitor", 4)))
  expect_true(all(!get_candidates(test_missing)$candidate))
  })

test_that("get_candidates returns FALSE on missing prediction", {
  test_na_prediction <- data.frame(interaction_types = c("inhibitor", 
                                    "activator", NA, NaN, 1, 0, -1,
                                    TRUE, FALSE),
                                   prediction = rep(NA, 9))
  expect_true(all(!get_candidates(test_na_prediction)$candidate))
  })