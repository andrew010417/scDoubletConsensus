#' Simulate a ground-truth doublet benchmark dataset
#'
#' No ground-truth doublet labels exist for real data without cell-hashing or
#' genotype demultiplexing, so accuracy can't be measured directly on a
#' typical dataset. This generates a synthetic dataset with known doublet
#' labels (via `scDblFinder::mockDoubletSCE()`, which builds doublets by
#' summing pairs of real single-cell profiles -- the standard approach used
#' in doublet-detection benchmark papers) and runs the full consensus
#' pipeline on it.
#'
#' @param n_cells Number of cells to simulate.
#' @param doublet_rate True doublet rate to simulate.
#' @param seed Random seed for reproducibility.
#' @param ... Additional arguments passed to [run_consensus_pipeline()].
#'
#' @return A list with `consensus_df` (as returned by
#'   [run_consensus_pipeline()], with an added `truth` column: `"doublet"`/
#'   `"singlet"`) and `truth` (the same labels as a standalone vector, in
#'   `consensus_df$cell` order).
#' @export
simulate_doublet_benchmark <- function(n_cells = 3000, doublet_rate = 0.1, seed = 1, ...) {
  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    stop("Package 'scDblFinder' is required for simulate_doublet_benchmark().")
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required for simulate_doublet_benchmark().")
  }
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for simulate_doublet_benchmark().")
  }

  set.seed(seed)
  # mockDoubletSCE()'s `ncells` is a per-cluster size vector (length >= 2),
  # not a total cell count, so split n_cells across two mock clusters.
  cluster_sizes <- c(ceiling(n_cells / 2), floor(n_cells / 2))
  sce <- scDblFinder::mockDoubletSCE(ncells = cluster_sizes, dbl.rate = doublet_rate)

  counts <- SingleCellExperiment::counts(sce)
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  truth <- ifelse(sce$type == "doublet", "doublet", "singlet")
  names(truth) <- colnames(counts)

  obj <- Seurat::CreateSeuratObject(counts = counts)
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, verbose = FALSE)
  obj <- Seurat::ScaleData(obj, verbose = FALSE)
  obj <- Seurat::RunPCA(obj, verbose = FALSE)
  obj <- Seurat::FindNeighbors(obj, dims = 1:10, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, verbose = FALSE)

  result <- run_consensus_pipeline(
    obj,
    run_qc = FALSE,
    PCs = 1:10,
    expected_doublet_rate = doublet_rate,
    seed = seed,
    ...
  )

  result$consensus_df$truth <- truth[result$consensus_df$cell]
  list(consensus_df = result$consensus_df, truth = truth)
}

#' Evaluate doublet calls against ground-truth labels
#'
#' Compares [run_consensus_pipeline()]'s `consensus_call` against the strict
#' intersection of `call_cols`, scoring both against known-correct labels.
#' Works with any `truth` vector -- synthetic (from
#' [simulate_doublet_benchmark()]) or real (e.g. cell-hashing/genotype
#' demultiplexing calls), so it doubles as the entry point for benchmarking
#' on a labeled dataset when one becomes available.
#'
#' @param consensus_df Output of [classify_consensus_doublets()], merged with
#'   the per-tool call columns named in `call_cols`.
#' @param truth Named character/logical vector of ground-truth labels
#'   (`"doublet"`/`"singlet"` or `TRUE`/`FALSE`), named by `cell_col` values,
#'   or an unnamed vector aligned row-for-row with `consensus_df`.
#' @param call_cols Character vector of per-tool call columns to intersect
#'   for the comparison baseline.
#' @param cell_col Name of the cell identifier column.
#'
#' @return A `data.frame` with one row per method (`"consensus"`,
#'   `"intersection"`) and columns `n_called`, `sensitivity`, `precision`,
#'   `specificity`, `f1`.
#' @export
evaluate_against_truth <- function(consensus_df, truth, call_cols, cell_col = "cell") {
  if (!"consensus_call" %in% colnames(consensus_df)) {
    stop("`consensus_df` must come from classify_consensus_doublets().")
  }
  if (length(call_cols) < 2 || !all(call_cols %in% colnames(consensus_df))) {
    stop("`call_cols` must name >= 2 existing per-tool call columns in `consensus_df`.")
  }

  is_doublet <- function(x) tolower(as.character(x)) %in% c("doublet", "true", "1")

  if (!is.null(names(truth))) {
    truth_vec <- is_doublet(truth[consensus_df[[cell_col]]])
  } else {
    if (length(truth) != nrow(consensus_df)) {
      stop("Unnamed `truth` must have one entry per row of `consensus_df`.")
    }
    truth_vec <- is_doublet(truth)
  }
  if (anyNA(truth_vec)) stop("`truth` is missing labels for some cells in `consensus_df`.")

  consensus_call <- is_doublet(consensus_df$consensus_call)
  intersection_call <- Reduce(`&`, lapply(call_cols, function(col) is_doublet(consensus_df[[col]])))

  score_one <- function(called) {
    tp <- sum(called & truth_vec)
    fp <- sum(called & !truth_vec)
    fn <- sum(!called & truth_vec)
    tn <- sum(!called & !truth_vec)
    sensitivity <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    precision <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
    specificity <- if (tn + fp > 0) tn / (tn + fp) else NA_real_
    f1 <- if (!is.na(sensitivity) && !is.na(precision) && (sensitivity + precision) > 0) {
      2 * sensitivity * precision / (sensitivity + precision)
    } else {
      NA_real_
    }
    data.frame(n_called = sum(called), sensitivity = sensitivity, precision = precision,
               specificity = specificity, f1 = f1)
  }

  result <- rbind(
    cbind(method = "consensus", score_one(consensus_call)),
    cbind(method = "intersection", score_one(intersection_call))
  )
  rownames(result) <- NULL
  result
}

#' Benchmark the consensus pipeline against a simulated ground truth
#'
#' Convenience wrapper combining [simulate_doublet_benchmark()] and
#' [evaluate_against_truth()] in one call.
#'
#' @inheritParams simulate_doublet_benchmark
#' @param call_cols Passed to [evaluate_against_truth()]. Defaults to the
#'   standard DoubletFinder/scDblFinder call columns produced by
#'   [run_consensus_pipeline()].
#'
#' @return A list with `metrics` (output of [evaluate_against_truth()]),
#'   `consensus_df`, and `truth`.
#' @export
benchmark_consensus <- function(n_cells = 3000,
                                 doublet_rate = 0.1,
                                 seed = 1,
                                 call_cols = c("doublet_finder_call", "scdblfinder_call"),
                                 ...) {
  sim <- simulate_doublet_benchmark(n_cells = n_cells, doublet_rate = doublet_rate, seed = seed, ...)
  metrics <- evaluate_against_truth(sim$consensus_df, sim$truth, call_cols = call_cols)
  list(metrics = metrics, consensus_df = sim$consensus_df, truth = sim$truth)
}

#' Grid-search the DoubletFinder/scDblFinder weight ratio against ground truth
#'
#' The package's default weighting between DoubletFinder's `pANN` and
#' scDblFinder's score is a literature-informed prior (see
#' `?run_consensus_pipeline`), not a dataset-specific optimum. When a labeled
#' benchmark is available (synthetic via [simulate_doublet_benchmark()], or
#' real ground truth), this searches a grid of weight splits and reports the
#' one that maximizes F1.
#'
#' @param score_df A `data.frame` with columns `cell`, `pANN`,
#'   `scDblFinder_score` (e.g. from merging [run_doubletfinder()] and
#'   [run_scdblfinder()] output).
#' @param truth Named or row-aligned ground-truth vector, as in
#'   [evaluate_against_truth()].
#' @param weight_grid Numeric vector in `[0, 1]`, each value the weight given
#'   to `scDblFinder_score` (with `1 - value` given to `pANN`).
#' @param method Passed to [consensus_rank()].
#' @param normalize_method Passed to [consensus_rank()].
#' @param multiplet_rate Doublet rate used for the top-N cutoff at each grid
#'   point (via `cutoff_method = "multiplet_rate"` in
#'   [classify_consensus_doublets()]). Defaults to the true rate implied by
#'   `truth`.
#' @param cell_col Name of the cell identifier column.
#'
#' @return A list with `grid` (a `data.frame` of every weight split and its
#'   resulting F1/sensitivity/precision) and `best_weights` (the named weight
#'   vector with the highest F1).
#' @export
optimize_weights <- function(score_df,
                              truth,
                              weight_grid = seq(0, 1, by = 0.05),
                              method = c("rank", "weighted"),
                              normalize_method = c("zscore", "minmax"),
                              multiplet_rate = NULL,
                              cell_col = "cell") {
  method <- match.arg(method)
  normalize_method <- match.arg(normalize_method)

  is_doublet <- function(x) tolower(as.character(x)) %in% c("doublet", "true", "1")
  if (!is.null(names(truth))) {
    truth_vec <- is_doublet(truth[score_df[[cell_col]]])
  } else {
    if (length(truth) != nrow(score_df)) {
      stop("Unnamed `truth` must have one entry per row of `score_df`.")
    }
    truth_vec <- is_doublet(truth)
  }
  if (anyNA(truth_vec)) stop("`truth` is missing labels for some cells in `score_df`.")

  if (is.null(multiplet_rate)) multiplet_rate <- mean(truth_vec)

  grid <- lapply(weight_grid, function(w_scdbl) {
    weights <- c(pANN = 1 - w_scdbl, scDblFinder_score = w_scdbl)
    consensus_df <- consensus_rank(
      score_df[, c(cell_col, "pANN", "scDblFinder_score")],
      cell_col = cell_col,
      weights = weights,
      method = method,
      normalize_method = normalize_method
    )
    n <- round(multiplet_rate * nrow(consensus_df))
    called <- consensus_df$consensus_rank <= n
    truth_ordered <- truth_vec[match(consensus_df[[cell_col]], score_df[[cell_col]])]

    tp <- sum(called & truth_ordered)
    fp <- sum(called & !truth_ordered)
    fn <- sum(!called & truth_ordered)
    sensitivity <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    precision <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
    f1 <- if (!is.na(sensitivity) && !is.na(precision) && (sensitivity + precision) > 0) {
      2 * sensitivity * precision / (sensitivity + precision)
    } else {
      NA_real_
    }
    data.frame(weight_scDblFinder = w_scdbl, weight_pANN = 1 - w_scdbl,
               sensitivity = sensitivity, precision = precision, f1 = f1)
  })
  grid <- do.call(rbind, grid)

  best_row <- grid[which.max(grid$f1), ]
  best_weights <- c(pANN = best_row$weight_pANN, scDblFinder_score = best_row$weight_scDblFinder)

  list(grid = grid, best_weights = best_weights)
}
