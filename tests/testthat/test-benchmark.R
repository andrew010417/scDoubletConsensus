test_that("evaluate_against_truth scores consensus and intersection against a named truth vector", {
  consensus_df <- data.frame(
    cell = c("c1", "c2", "c3", "c4"),
    consensus_call = c("doublet", "doublet", "singlet", "singlet"),
    doublet_finder_call = c("Doublet", "Doublet", "Singlet", "Singlet"),
    scdblfinder_call = c("singlet", "doublet", "singlet", "singlet")
  )
  truth <- c(c1 = "doublet", c2 = "doublet", c3 = "singlet", c4 = "singlet")

  metrics <- evaluate_against_truth(
    consensus_df, truth,
    call_cols = c("doublet_finder_call", "scdblfinder_call")
  )

  consensus_row <- metrics[metrics$method == "consensus", ]
  intersection_row <- metrics[metrics$method == "intersection", ]

  # consensus calls c1 & c2 as doublet; both are true doublets -> perfect
  expect_equal(consensus_row$sensitivity, 1)
  expect_equal(consensus_row$precision, 1)
  expect_equal(consensus_row$f1, 1)

  # intersection only calls c2 -> misses c1, still perfect precision
  expect_equal(intersection_row$n_called, 1)
  expect_equal(intersection_row$sensitivity, 0.5)
  expect_equal(intersection_row$precision, 1)
})

test_that("evaluate_against_truth works with an unnamed, row-aligned truth vector", {
  consensus_df <- data.frame(
    cell = c("c1", "c2"),
    consensus_call = c("doublet", "singlet"),
    doublet_finder_call = c("Doublet", "Singlet"),
    scdblfinder_call = c("doublet", "singlet")
  )
  truth <- c("doublet", "singlet")

  metrics <- evaluate_against_truth(
    consensus_df, truth,
    call_cols = c("doublet_finder_call", "scdblfinder_call")
  )
  expect_equal(metrics$sensitivity[metrics$method == "consensus"], 1)
})

test_that("evaluate_against_truth errors on mismatched unnamed truth length", {
  consensus_df <- data.frame(
    cell = c("c1", "c2"),
    consensus_call = c("doublet", "singlet"),
    doublet_finder_call = c("Doublet", "Singlet"),
    scdblfinder_call = c("doublet", "singlet")
  )
  expect_error(evaluate_against_truth(
    consensus_df, c("doublet"),
    call_cols = c("doublet_finder_call", "scdblfinder_call")
  ))
})

test_that("optimize_weights finds the weight split favoring the more predictive tool", {
  # scDblFinder_score perfectly separates truth; pANN jitters with no relation to truth
  score_df <- data.frame(
    cell = paste0("cell", 1:10),
    pANN = c(0.50, 0.48, 0.52, 0.49, 0.51, 0.50, 0.52, 0.48, 0.51, 0.49),
    scDblFinder_score = c(0.9, 0.85, 0.8, 0.75, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35)
  )
  truth <- c(rep("doublet", 4), rep("singlet", 6))

  result <- optimize_weights(
    score_df, truth,
    weight_grid = c(0, 0.5, 1),
    method = "weighted",
    normalize_method = "minmax"
  )

  expect_true("grid" %in% names(result))
  expect_true("best_weights" %in% names(result))
  # weighting scDblFinder_score fully (1.0) should perform best since it's the
  # only informative signal
  expect_equal(unname(result$best_weights["scDblFinder_score"]), 1)
})
