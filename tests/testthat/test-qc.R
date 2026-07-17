make_toy_seurat <- function() {
  set.seed(42)
  n_genes <- 50
  n_high_quality <- 10
  n_low_quality <- 10
  n_cells <- n_high_quality + n_low_quality

  # Low-quality cells get far fewer counts/genes detected, so a threshold
  # between the two groups produces a meaningful partial filter.
  high_quality <- matrix(
    rpois(n_genes * n_high_quality, lambda = 20),
    nrow = n_genes
  )
  low_quality <- matrix(
    rpois(n_genes * n_low_quality, lambda = 0.2),
    nrow = n_genes
  )
  counts <- cbind(high_quality, low_quality)
  dimnames(counts) <- list(
    c(paste0("MT-GENE", 1:5), paste0("GENE", 1:(n_genes - 5))),
    paste0("cell", 1:n_cells)
  )

  Seurat::CreateSeuratObject(counts = counts)
}

test_that("compute_qc_metrics adds percent.mt using the mt_pattern", {
  skip_if_not_installed("Seurat")
  obj <- make_toy_seurat()
  obj <- compute_qc_metrics(obj, mt_pattern = "^MT-")

  expect_true("percent.mt" %in% colnames(obj[[]]))
  expect_true(all(obj$percent.mt >= 0 & obj$percent.mt <= 100))
})

test_that("compute_qc_metrics rejects non-Seurat input", {
  expect_error(compute_qc_metrics(list(a = 1)))
})

test_that("filter_cells keeps all cells when thresholds are permissive", {
  skip_if_not_installed("Seurat")
  obj <- make_toy_seurat()
  obj <- compute_qc_metrics(obj, mt_pattern = "^MT-")

  filtered <- filter_cells(
    obj,
    min_genes = 0,
    max_genes = Inf,
    min_counts = 0,
    max_mt_percent = 100
  )
  expect_equal(ncol(filtered), ncol(obj))
})

test_that("filter_cells drops only the cells failing a threshold", {
  skip_if_not_installed("Seurat")
  obj <- make_toy_seurat()
  obj <- compute_qc_metrics(obj, mt_pattern = "^MT-")

  # Threshold sits between the toy high- and low-quality cell groups.
  partial <- filter_cells(
    obj,
    min_genes = 10,
    min_counts = 0,
    max_mt_percent = 100
  )
  expect_true(ncol(partial) > 0 && ncol(partial) < ncol(obj))
  expect_true(all(partial$nFeature_RNA >= 10))
})

test_that("filter_cells surfaces Seurat's error when every cell is excluded", {
  skip_if_not_installed("Seurat")
  obj <- make_toy_seurat()
  obj <- compute_qc_metrics(obj, mt_pattern = "^MT-")

  expect_error(
    filter_cells(obj, min_genes = 10^6, min_counts = 0, max_mt_percent = 100),
    "No cells found"
  )
})

test_that("filter_cells errors if QC metrics were not computed first", {
  skip_if_not_installed("Seurat")
  obj <- make_toy_seurat()
  expect_error(filter_cells(obj), "compute_qc_metrics")
})
