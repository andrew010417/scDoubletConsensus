#' Run the full QC + consensus doublet detection pipeline
#'
#' Convenience wrapper: QC filtering, DoubletFinder, scDblFinder, and
#' consensus scoring in one call. For finer control (e.g. inspecting QC
#' plots before choosing thresholds) call the individual functions directly.
#'
#' @param seurat_obj A `Seurat` object with raw counts, already normalized
#'   and dimensionality-reduced (`RunPCA()` already run) if `run_qc = FALSE`.
#' @param run_qc Whether to compute QC metrics and filter cells first.
#' @param qc_args Named list of arguments passed to [filter_cells()].
#' @param mt_pattern Passed to [compute_qc_metrics()].
#' @param PCs Principal components passed to [run_doubletfinder()], or
#'   `"auto"` to pick them via [choose_pcs()].
#' @param expected_doublet_rate Passed to [run_doubletfinder()], or `"auto"`
#'   to derive it from cell count via [estimate_multiplet_rate()] -- run
#'   automatically at the start of every call (see the reported estimate in
#'   the printed message from [inspect_input()]).
#' @param consensus_method `"rank"` or `"weighted"`, passed to
#'   [consensus_rank()].
#' @param weights Optional named weights for `pANN`/`scDblFinder_score`,
#'   passed to [consensus_rank()]. Defaults to `c(pANN = 0.4,
#'   scDblFinder_score = 0.6)`: scDblFinder is weighted higher based on
#'   published benchmarks (e.g. Xi & Li 2021) generally finding it more
#'   accurate and less parameter-sensitive than DoubletFinder. This is a
#'   literature-informed prior, not a dataset-specific optimum -- pass your
#'   own `weights`, or derive one from your data with [optimize_weights()].
#' @param top_n,quantile Passed to [classify_consensus_doublets()] when
#'   `cutoff_method = "legacy"`.
#' @param cutoff_method,multiplet_rate Passed to
#'   [classify_consensus_doublets()]. `call_cols` for
#'   `cutoff_method = "intersection_median"` is fixed to
#'   `c("doublet_finder_call", "scdblfinder_call")`.
#' @param seed Random seed for reproducibility.
#'
#' @return A list with `seurat_obj` (post-QC), `consensus_df` (per-cell
#'   scores and calls), and `comparison` (output of
#'   [compare_to_intersection()]).
#' @export
run_consensus_pipeline <- function(seurat_obj,
                                    run_qc = TRUE,
                                    qc_args = list(),
                                    mt_pattern = "^MT-",
                                    PCs = 1:10,
                                    expected_doublet_rate = 0.075,
                                    consensus_method = c("rank", "weighted"),
                                    weights = c(pANN = 0.4, scDblFinder_score = 0.6),
                                    top_n = NULL,
                                    quantile = 0.9,
                                    cutoff_method = c("legacy", "intersection_median", "multiplet_rate"),
                                    multiplet_rate = NULL,
                                    seed = 1) {
  consensus_method <- match.arg(consensus_method)
  cutoff_method <- match.arg(cutoff_method)

  inspect_input(seurat_obj)

  if (run_qc) {
    seurat_obj <- compute_qc_metrics(seurat_obj, mt_pattern = mt_pattern)
    seurat_obj <- do.call(filter_cells, c(list(seurat_obj = seurat_obj), qc_args))
  }

  if (identical(expected_doublet_rate, "auto")) {
    expected_doublet_rate <- estimate_multiplet_rate(ncol(seurat_obj))
  }

  df_scores <- run_doubletfinder(
    seurat_obj,
    PCs = PCs,
    expected_doublet_rate = expected_doublet_rate,
    seed = seed
  )
  scdbl_scores <- run_scdblfinder(seurat_obj, seed = seed)

  merged <- merge(df_scores, scdbl_scores, by = "cell")

  consensus_df <- consensus_rank(
    merged[, c("cell", "pANN", "scDblFinder_score")],
    weights = weights,
    method = consensus_method
  )
  consensus_df <- merge(
    consensus_df,
    merged[, c("cell", "doublet_finder_call", "scdblfinder_call")],
    by = "cell"
  )
  consensus_df <- classify_consensus_doublets(
    consensus_df,
    top_n = top_n,
    quantile = quantile,
    cutoff_method = cutoff_method,
    call_cols = c("doublet_finder_call", "scdblfinder_call"),
    multiplet_rate = multiplet_rate
  )
  consensus_df <- consensus_df[order(consensus_df$consensus_rank), ]
  rownames(consensus_df) <- NULL

  comparison <- compare_to_intersection(
    consensus_df,
    call_cols = c("doublet_finder_call", "scdblfinder_call")
  )

  list(
    seurat_obj = seurat_obj,
    consensus_df = consensus_df,
    comparison = comparison
  )
}
