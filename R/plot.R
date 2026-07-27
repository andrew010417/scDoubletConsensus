#' Plot where intersection-called cells fall in the consensus score distribution
#'
#' Answers the question "where do the strict-intersection doublets sit within
#' the full consensus score distribution?" -- every cell is plotted by
#' `consensus_rank` (x) vs. `consensus_score` (y), colored by whether it was
#' called a doublet by every tool in `call_cols` (the intersection), with
#' shape marking the final `consensus_call`. Optional cutoff lines from
#' [classify_consensus_doublets()]'s `"intersection_median"` and
#' `"multiplet_rate"` methods are overlaid for comparison.
#'
#' @param consensus_df Output of [classify_consensus_doublets()], merged with
#'   the per-tool call columns named in `call_cols`.
#' @param call_cols Character vector of per-tool call columns to intersect
#'   for the color aesthetic (e.g. `c("doublet_finder_call",
#'   "scdblfinder_call")`).
#' @param show_intersection_median_cutoff Draw a horizontal line at the
#'   intersection-median cutoff (see [classify_consensus_doublets()]).
#' @param show_multiplet_rate_cutoff Draw a vertical line at the multiplet-
#'   rate top-N cutoff (see [classify_consensus_doublets()]).
#' @param multiplet_rate Passed to [estimate_multiplet_rate()]'s override when
#'   drawing the multiplet-rate cutoff; defaults to the automatic estimate.
#'
#' @return A `ggplot` object.
#' @export
plot_consensus_overview <- function(consensus_df,
                                     call_cols,
                                     show_intersection_median_cutoff = TRUE,
                                     show_multiplet_rate_cutoff = TRUE,
                                     multiplet_rate = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_consensus_overview().")
  }
  if (!"consensus_score" %in% colnames(consensus_df) ||
        !"consensus_rank" %in% colnames(consensus_df)) {
    stop("`consensus_df` must come from consensus_rank()/classify_consensus_doublets().")
  }
  if (length(call_cols) < 2 || !all(call_cols %in% colnames(consensus_df))) {
    stop("`call_cols` must name >= 2 existing per-tool call columns in `consensus_df`.")
  }

  is_doublet <- function(x) tolower(as.character(x)) %in% c("doublet", "true", "1")
  consensus_df$intersection_call <- Reduce(`&`, lapply(call_cols, function(col) is_doublet(consensus_df[[col]])))
  consensus_df$intersection_call <- ifelse(consensus_df$intersection_call, "intersection doublet", "not in intersection")

  p <- ggplot2::ggplot(
    consensus_df,
    ggplot2::aes(
      x = consensus_rank,
      y = consensus_score,
      color = intersection_call,
      shape = consensus_call
    )
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::labs(
      x = "Consensus rank (1 = most doublet-like)",
      y = "Consensus score",
      color = "Per-tool intersection",
      shape = "Consensus call",
      title = "Consensus score distribution: where do intersection doublets fall?"
    ) +
    ggplot2::theme_minimal()

  if (show_intersection_median_cutoff && any(consensus_df$intersection_call == "intersection doublet")) {
    median_cutoff <- stats::median(
      consensus_df$consensus_score[consensus_df$intersection_call == "intersection doublet"],
      na.rm = TRUE
    )
    p <- p + ggplot2::geom_hline(
      yintercept = median_cutoff, linetype = "dashed", color = "black"
    ) +
      ggplot2::annotate(
        "text", x = -Inf, y = median_cutoff, label = "intersection median",
        hjust = -0.05, vjust = -0.5, size = 3
      )
  }

  if (show_multiplet_rate_cutoff) {
    if (is.null(multiplet_rate)) {
      multiplet_rate <- estimate_multiplet_rate(nrow(consensus_df))
    }
    n_cutoff <- round(multiplet_rate * nrow(consensus_df))
    p <- p + ggplot2::geom_vline(
      xintercept = n_cutoff, linetype = "dotted", color = "black"
    ) +
      ggplot2::annotate(
        "text", x = n_cutoff, y = Inf, label = "multiplet-rate cutoff",
        hjust = -0.05, vjust = 1.5, size = 3
      )
  }

  p
}
