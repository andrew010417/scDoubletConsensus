make_toy_seurat <- function() {
  set.seed(42)
  n_genes <- 50
  n_cells <- 20
  counts <- matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes,
    dimnames = list(
      c(paste0("MT-GENE", 1:5), paste0("GENE", 1:(n_genes - 5))),
      paste0("cell", 1:n_cells)
    )
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

test_that("filter_cells keeps only cells passing all thresholds", {
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

  strict <- filter_cells(
    obj,
    min_genes = 10^6,
    min_counts = 0,
    max_mt_percent = 100
  )
  expect_equal(ncol(strict), 0)
})

test_that("filter_cells errors if QC metrics were not computed first", {
  skip_if_not_installed("Seurat")
  obj <- make_toy_seurat()
  expect_error(filter_cells(obj), "compute_qc_metrics")
})
