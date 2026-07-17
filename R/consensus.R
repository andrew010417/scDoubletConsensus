#' Normalize a vector of doublet scores
#'
#' Puts raw scores from different detectors (e.g. DoubletFinder's `pANN`,
#' scDblFinder's probability) on a common scale so they can be combined.
#'
#' @param x Numeric vector of raw scores.
#' @param method `"zscore"` (mean 0, sd 1) or `"minmax"` (rescaled to
#'   `[0, 1]`).
#'
#' @return A numeric vector the same length as `x`.
#' @export
normalize_scores <- function(x, method = c("zscore", "minmax")) {
  method <- match.arg(method)
  if (!is.numeric(x)) stop("`x` must be numeric.")

  if (method == "zscore") {
    s <- stats::sd(x, na.rm = TRUE)
    if (is.na(s) || s == 0) stop("Cannot z-score a constant (zero-variance) vector.")
    (x - mean(x, na.rm = TRUE)) / s
  } else {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) stop("Cannot min-max scale a constant vector.")
    (x - rng[1]) / diff(rng)
  }
}

#' Build a rank- or weight-based doublet consensus
#'
#' This is the core innovation of the package: rather than intersecting
#' per-tool doublet calls, normalize each tool's score and combine them into
#' a single consensus score/rank. This recovers cells that individual tools
#' scored with high confidence but that fall outside the strict intersection
#' of binary calls.
#'
#' @param score_df A `data.frame` with one row per cell: an identifier column
#'   (`cell_col`) plus one numeric column of raw scores per detection tool
#'   (e.g. `pANN`, `scDblFinder_score`).
#' @param cell_col Name of the cell identifier column.
#' @param weights Optional named numeric vector giving a weight per score
#'   column (names must match `score_df` column names). Defaults to equal
#'   weighting. Only used when `method = "weighted"`.
#' @param method `"rank"` averages each tool's per-cell rank (weighted);
#'   `"weighted"` sums normalized scores directly (weighted).
#' @param normalize_method Passed to [normalize_scores()].
#'
#' @return `score_df` with added `<tool>_norm` columns (and `<tool>_rank`
#'   columns when `method = "rank"`), plus `consensus_score` and
#'   `consensus_rank` (1 = most doublet-like), sorted by `consensus_rank`.
#' @export
consensus_rank <- function(score_df,
                            cell_col = "cell",
                            weights = NULL,
                            method = c("rank", "weighted"),
                            normalize_method = c("zscore", "minmax")) {
  method <- match.arg(method)
  normalize_method <- match.arg(normalize_method)

  if (!cell_col %in% colnames(score_df)) {
    stop("`score_df` must contain the column named in `cell_col`.")
  }
  score_cols <- setdiff(colnames(score_df), cell_col)
  if (length(score_cols) < 2) {
    stop("Need at least two score columns to build a consensus.")
  }

  if (is.null(weights)) {
    weights <- stats::setNames(rep(1 / length(score_cols), length(score_cols)), score_cols)
  } else {
    if (is.null(names(weights)) || !all(score_cols %in% names(weights))) {
      stop("`weights` must be a named numeric vector covering every score column.")
    }
    weights <- weights[score_cols]
    weights <- weights / sum(weights)
  }

  result <- score_df
  norm_cols <- paste0(score_cols, "_norm")
  for (i in seq_along(score_cols)) {
    result[[norm_cols[i]]] <- normalize_scores(score_df[[score_cols[i]]], method = normalize_method)
  }

  if (method == "weighted") {
    consensus_score <- rep(0, nrow(result))
    for (i in seq_along(score_cols)) {
      consensus_score <- consensus_score + result[[norm_cols[i]]] * weights[[score_cols[i]]]
    }
  } else {
    rank_cols <- paste0(score_cols, "_rank")
    for (i in seq_along(score_cols)) {
      # higher normalized score = more doublet-like => rank 1 is most doublet-like
      result[[rank_cols[i]]] <- rank(-result[[norm_cols[i]]], ties.method = "average")
    }
    avg_rank <- rep(0, nrow(result))
    for (i in seq_along(score_cols)) {
      avg_rank <- avg_rank + result[[rank_cols[i]]] * weights[[score_cols[i]]]
    }
    consensus_score <- -avg_rank
  }

  result$consensus_score <- consensus_score
  result$consensus_rank <- rank(-result$consensus_score, ties.method = "min")
  result <- result[order(result$consensus_rank), ]
  rownames(result) <- NULL
  result
}

#' Classify cells as doublets from a consensus score
#'
#' @param consensus_df Output of [consensus_rank()].
#' @param top_n If set, the `top_n` cells by `consensus_rank` are called
#'   doublets. Takes precedence over `quantile`.
#' @param quantile Quantile of `consensus_score` (0-1) above which a cell is
#'   called a doublet. Ignored if `top_n` is set. Default `0.9`.
#'
#' @return `consensus_df` with an added `consensus_call` column
#'   (`"doublet"`/`"singlet"`).
#' @export
classify_consensus_doublets <- function(consensus_df, top_n = NULL, quantile = 0.9) {
  if (!"consensus_score" %in% colnames(consensus_df) ||
        !"consensus_rank" %in% colnames(consensus_df)) {
    stop("`consensus_df` must come from consensus_rank().")
  }

  if (!is.null(top_n)) {
    consensus_df$consensus_call <- ifelse(consensus_df$consensus_rank <= top_n, "doublet", "singlet")
  } else {
    threshold <- stats::quantile(consensus_df$consensus_score, probs = quantile, na.rm = TRUE)
    consensus_df$consensus_call <- ifelse(consensus_df$consensus_score >= threshold, "doublet", "singlet")
  }
  consensus_df
}

#' Compare consensus calls against simple intersection filtering
#'
#' Quantifies the README's motivating claim: intersection-based filtering
#' "can miss cells that receive very high confidence scores from individual
#' tools but don't appear in the intersection." Reports how many cells the
#' consensus call recovers relative to the intersection of raw per-tool
#' calls, and how many it drops.
#'
#' @param consensus_df Output of [classify_consensus_doublets()], optionally
#'   merged with each tool's raw call column (e.g.
#'   `doublet_finder_call`, `scdblfinder_call`).
#' @param call_cols Character vector naming the per-tool call columns to
#'   intersect. Values are matched case-insensitively against
#'   `"doublet"`/`"1"`/`"TRUE"`.
#' @param cell_col Name of the cell identifier column.
#'
#' @return A list with `summary` (one-row `data.frame` of counts) and
#'   `recovered_cells` (identifiers called doublets by consensus but not by
#'   the intersection of individual tools).
#' @export
compare_to_intersection <- function(consensus_df, call_cols, cell_col = "cell") {
  if (!"consensus_call" %in% colnames(consensus_df)) {
    stop("`consensus_df` must come from classify_consensus_doublets().")
  }
  if (length(call_cols) < 2) {
    stop("Need at least two per-tool call columns to form an intersection.")
  }

  is_doublet <- function(x) tolower(as.character(x)) %in% c("doublet", "true", "1")
  intersection_call <- Reduce(`&`, lapply(call_cols, function(col) is_doublet(consensus_df[[col]])))
  consensus_call <- tolower(consensus_df$consensus_call) == "doublet"

  recovered <- consensus_df[consensus_call & !intersection_call, cell_col, drop = FALSE]

  summary <- data.frame(
    n_cells = nrow(consensus_df),
    n_intersection_doublets = sum(intersection_call),
    n_consensus_doublets = sum(consensus_call),
    n_recovered_by_consensus = nrow(recovered),
    n_lost_by_consensus = sum(intersection_call & !consensus_call)
  )

  list(summary = summary, recovered_cells = recovered[[cell_col]])
}
