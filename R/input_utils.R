#' Estimate expected multiplet rate from cell count
#'
#' Interpolates 10x Genomics' published loading table (recovered cells vs.
#' expected multiplet rate) so a dataset-appropriate doublet rate can be
#' derived automatically instead of hard-coding `0.075` for every dataset.
#'
#' @param n_cells Number of cells (e.g. `ncol(seurat_obj)`).
#'
#' @return A single numeric multiplet rate (e.g. `0.046` for ~6,000 cells).
#' @export
estimate_multiplet_rate <- function(n_cells) {
  if (!is.numeric(n_cells) || length(n_cells) != 1 || n_cells <= 0) {
    stop("`n_cells` must be a single positive number.")
  }

  # 10x Genomics Chromium loading table: cells recovered -> multiplet rate.
  anchor_cells <- c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
  anchor_rate <- c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076)

  if (n_cells <= max(anchor_cells)) {
    stats::approx(anchor_cells, anchor_rate, xout = n_cells, rule = 2)$y
  } else {
    # extrapolate beyond 10,000 cells using the table's final slope
    slope <- (anchor_rate[length(anchor_rate)] - anchor_rate[length(anchor_rate) - 1]) /
      (anchor_cells[length(anchor_cells)] - anchor_cells[length(anchor_cells) - 1])
    anchor_rate[length(anchor_rate)] + slope * (n_cells - anchor_cells[length(anchor_cells)])
  }
}

#' Inspect a Seurat object's size and suggest a multiplet rate
#'
#' Run this on raw input before setting up a doublet-detection run. Detects
#' cell and gene counts and derives the expected multiplet rate from
#' [estimate_multiplet_rate()], so `expected_doublet_rate` does not need to be
#' guessed or hard-coded per dataset.
#'
#' @param seurat_obj A `Seurat` object.
#'
#' @return A list with `n_cells`, `n_genes`, and `estimated_multiplet_rate`.
#' @export
inspect_input <- function(seurat_obj) {
  if (!methods::is(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }

  n_cells <- ncol(seurat_obj)
  n_genes <- nrow(seurat_obj)
  rate <- estimate_multiplet_rate(n_cells)

  message(sprintf(
    "Input: %d cells x %d genes -> estimated multiplet rate %.3f (%.1f%%)",
    n_cells, n_genes, rate, rate * 100
  ))

  list(n_cells = n_cells, n_genes = n_genes, estimated_multiplet_rate = rate)
}

#' Choose principal components via an elbow heuristic
#'
#' Picks the smallest number of leading PCs such that cumulative variance
#' explained exceeds `cumulative_cutoff`% while the individual PC still
#' explains at least `variation_cutoff`% -- falling back to the point where
#' variance explained by consecutive PCs differs by less than 0.1 percentage
#' points. This mirrors the elbow-plot heuristic commonly used with Seurat's
#' `ElbowPlot()`, automated so `PCs = "auto"` can be passed to
#' [run_consensus_pipeline()] instead of eyeballing a plot.
#'
#' @param x Either a `Seurat` object with `RunPCA()` already run, or a numeric
#'   vector of per-PC standard deviations (`Seurat::Stdev(seurat_obj, "pca")`).
#' @param cumulative_cutoff Minimum cumulative percent variance explained.
#' @param variation_cutoff Maximum percent variance explained by a single PC
#'   at the candidate cutoff.
#' @param reduction Reduction name to read `Stdev()` from when `x` is a
#'   `Seurat` object.
#'
#' @return An integer vector `1:n` giving the chosen PCs.
#' @export
choose_pcs <- function(x, cumulative_cutoff = 90, variation_cutoff = 5, reduction = "pca") {
  if (methods::is(x, "Seurat")) {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
      stop("Package 'Seurat' is required for choose_pcs() on a Seurat object.")
    }
    if (!reduction %in% names(x@reductions)) {
      stop(sprintf(
        "`%s` reduction not found. Run Seurat::RunPCA() first, or pass PCs explicitly.",
        reduction
      ))
    }
    stdev <- Seurat::Stdev(x, reduction = reduction)
  } else if (is.numeric(x)) {
    stdev <- x
  } else {
    stop("`x` must be a Seurat object or a numeric vector of PC standard deviations.")
  }

  if (length(stdev) < 2) stop("Need at least two PCs to choose from.")

  pct_var <- 100 * stdev^2 / sum(stdev^2)
  cum_var <- cumsum(pct_var)

  co1 <- which(cum_var > cumulative_cutoff & pct_var < variation_cutoff)[1]
  co2 <- which(diff(pct_var) < 0.1)[1]

  n_pcs <- min(co1, co2, na.rm = TRUE)
  if (!is.finite(n_pcs)) n_pcs <- length(stdev)
  # A single PC breaks downstream matrix ops that expect >= 2 dimensions
  # (e.g. diag() on a length-1 numeric builds an identity matrix, not a
  # 1x1 diagonal one), so never return fewer than 2.
  n_pcs <- max(n_pcs, 2)

  seq_len(n_pcs)
}
