test_that("estimate_multiplet_rate matches known 10x anchor points", {
  expect_equal(estimate_multiplet_rate(1000), 0.008)
  expect_equal(estimate_multiplet_rate(5000), 0.039)
  expect_equal(estimate_multiplet_rate(10000), 0.076)
})

test_that("estimate_multiplet_rate interpolates between anchor points", {
  rate <- estimate_multiplet_rate(1500)
  expect_true(rate > 0.008 && rate < 0.016)
})

test_that("estimate_multiplet_rate extrapolates beyond the table", {
  rate_at_10k <- estimate_multiplet_rate(10000)
  rate_beyond <- estimate_multiplet_rate(15000)
  expect_true(rate_beyond > rate_at_10k)
})

test_that("estimate_multiplet_rate errors on invalid input", {
  expect_error(estimate_multiplet_rate(-5))
  expect_error(estimate_multiplet_rate(c(1000, 2000)))
  expect_error(estimate_multiplet_rate("1000"))
})

test_that("choose_pcs picks fewer PCs when variance concentrates in early components", {
  # first few PCs carry almost all the variance -> should pick a small elbow
  stdev <- c(10, 9, 8, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
  pcs <- choose_pcs(stdev)
  expect_true(length(pcs) < length(stdev))
  expect_equal(pcs, seq_len(length(pcs)))
})

test_that("choose_pcs requires at least two PCs", {
  expect_error(choose_pcs(c(1)))
})

test_that("choose_pcs rejects non-numeric, non-Seurat input", {
  expect_error(choose_pcs("not a seurat object"))
})
