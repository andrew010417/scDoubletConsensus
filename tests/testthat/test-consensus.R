test_that("normalize_scores zscore has mean 0 and sd 1", {
  x <- c(1, 2, 3, 4, 10)
  z <- normalize_scores(x, method = "zscore")
  expect_equal(mean(z), 0, tolerance = 1e-10)
  expect_equal(stats::sd(z), 1, tolerance = 1e-10)
})

test_that("normalize_scores minmax rescales to [0, 1]", {
  x <- c(5, 10, 15, 20)
  m <- normalize_scores(x, method = "minmax")
  expect_equal(min(m), 0)
  expect_equal(max(m), 1)
})

test_that("normalize_scores errors on constant input", {
  expect_error(normalize_scores(rep(1, 5), method = "zscore"))
  expect_error(normalize_scores(rep(1, 5), method = "minmax"))
})

test_that("consensus_rank ranks the cell with the best combined score first", {
  score_df <- data.frame(
    cell = paste0("cell", 1:5),
    toolA = c(0.1, 0.2, 0.9, 0.3, 0.95),
    toolB = c(0.05, 0.9, 0.85, 0.2, 0.99)
  )

  result <- consensus_rank(score_df, method = "weighted", normalize_method = "minmax")

  expect_true("consensus_score" %in% colnames(result))
  expect_true("consensus_rank" %in% colnames(result))
  # cell5 has the top score in both tools, should be rank 1
  expect_equal(result$cell[result$consensus_rank == 1], "cell5")
  # result must be sorted by consensus_rank ascending
  expect_true(!is.unsorted(result$consensus_rank))
})

test_that("consensus_rank respects custom weights", {
  score_df <- data.frame(
    cell = paste0("cell", 1:4),
    toolA = c(0.1, 0.9, 0.2, 0.3),
    toolB = c(0.9, 0.1, 0.2, 0.3)
  )

  # Weighting toolA heavily should push cell2 (high toolA) to the top
  result <- consensus_rank(
    score_df,
    weights = c(toolA = 0.9, toolB = 0.1),
    method = "weighted",
    normalize_method = "minmax"
  )
  expect_equal(result$cell[result$consensus_rank == 1], "cell2")
})

test_that("consensus_rank requires at least two score columns", {
  score_df <- data.frame(cell = c("a", "b"), toolA = c(0.1, 0.2))
  expect_error(consensus_rank(score_df))
})

test_that("classify_consensus_doublets top_n labels exactly top_n cells", {
  score_df <- data.frame(
    cell = paste0("cell", 1:10),
    toolA = seq(0.1, 1.0, length.out = 10),
    toolB = seq(0.1, 1.0, length.out = 10)
  )
  consensus_df <- consensus_rank(score_df, method = "weighted", normalize_method = "minmax")
  classified <- classify_consensus_doublets(consensus_df, top_n = 3)

  expect_equal(sum(classified$consensus_call == "doublet"), 3)
  expect_setequal(
    classified$cell[classified$consensus_call == "doublet"],
    c("cell8", "cell9", "cell10")
  )
})

test_that("classify_consensus_doublets intersection_median uses the intersection's median score as threshold", {
  consensus_df <- data.frame(
    cell = paste0("cell", 1:6),
    consensus_score = c(10, 8, 6, 4, 2, 0),
    consensus_rank = 1:6,
    doublet_finder_call = c("Doublet", "Doublet", "Doublet", "Singlet", "Singlet", "Singlet"),
    scdblfinder_call = c("doublet", "singlet", "doublet", "doublet", "singlet", "singlet")
  )
  # intersection (both tools agree) = cell1 (score 10), cell3 (score 6) -> median = 8
  result <- classify_consensus_doublets(
    consensus_df,
    cutoff_method = "intersection_median",
    call_cols = c("doublet_finder_call", "scdblfinder_call")
  )
  expect_setequal(result$cell[result$consensus_call == "doublet"], c("cell1", "cell2"))
})

test_that("classify_consensus_doublets intersection_median errors without call_cols", {
  consensus_df <- data.frame(
    cell = paste0("cell", 1:3),
    consensus_score = c(3, 2, 1),
    consensus_rank = 1:3
  )
  expect_error(classify_consensus_doublets(consensus_df, cutoff_method = "intersection_median"))
})

test_that("classify_consensus_doublets multiplet_rate calls round(rate * n) top cells", {
  consensus_df <- data.frame(
    cell = paste0("cell", 1:10),
    consensus_score = seq(10, 1, by = -1),
    consensus_rank = 1:10
  )
  result <- classify_consensus_doublets(
    consensus_df,
    cutoff_method = "multiplet_rate",
    multiplet_rate = 0.3
  )
  expect_equal(sum(result$consensus_call == "doublet"), 3)
  expect_setequal(result$cell[result$consensus_call == "doublet"], c("cell1", "cell2", "cell3"))
})

test_that("classify_consensus_doublets multiplet_rate defaults to estimate_multiplet_rate", {
  consensus_df <- data.frame(
    cell = paste0("cell", 1:6000),
    consensus_score = seq(6000, 1, by = -1),
    consensus_rank = 1:6000
  )
  result <- classify_consensus_doublets(consensus_df, cutoff_method = "multiplet_rate")
  expect_equal(sum(result$consensus_call == "doublet"), round(estimate_multiplet_rate(6000) * 6000))
})

test_that("compare_to_intersection recovers cells missed by strict intersection", {
  # This encodes the README's core motivating claim: a cell can score very
  # high on the tools' underlying scores yet be missed by binary-call
  # intersection, and the consensus should recover it.
  consensus_df <- data.frame(
    cell = c("c1", "c2", "c3", "c4"),
    consensus_score = c(3, 2, 1, 0),
    consensus_rank = c(1, 2, 3, 4),
    consensus_call = c("doublet", "doublet", "singlet", "singlet"),
    doublet_finder_call = c("Doublet", "Doublet", "Singlet", "Singlet"),
    scdblfinder_call = c("singlet", "doublet", "singlet", "singlet")
  )

  cmp <- compare_to_intersection(
    consensus_df,
    call_cols = c("doublet_finder_call", "scdblfinder_call")
  )

  expect_equal(cmp$summary$n_intersection_doublets, 1) # only c2
  expect_equal(cmp$summary$n_consensus_doublets, 2) # c1, c2
  expect_equal(cmp$recovered_cells, "c1") # high-confidence cell missed by intersection
  expect_equal(cmp$summary$n_recovered_by_consensus, 1)
  expect_equal(cmp$summary$n_lost_by_consensus, 0)
})
