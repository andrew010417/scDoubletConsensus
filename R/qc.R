#' Compute per-cell QC metrics
#'
#' Adds gene count per cell, total RNA counts, and mitochondrial percentage
#' to the metadata of a Seurat object. No universal cutoff exists for any of
#' these metrics; use [plot_qc_metrics()] to pick dataset-specific thresholds
#' before calling [filter_cells()].
#'
#' @param seurat_obj A `Seurat` object.
#' @param mt_pattern Regex identifying mitochondrial genes. Defaults to the
#'   human convention `"^MT-"`; use `"^mt-"` for mouse.
#'
#' @return The `Seurat` object with `nFeature_RNA`, `nCount_RNA`, and
#'   `percent.mt` columns in `meta.data`.
#' @export
compute_qc_metrics <- function(seurat_obj, mt_pattern = "^MT-") {
  if (!methods::is(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }

  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(
    seurat_obj,
    pattern = mt_pattern
  )

  seurat_obj
}

#' Plot QC metrics to guide threshold selection
#'
#' Produces the violin and scatter plots used to eyeball reasonable QC
#' cutoffs, since gold-standard values do not exist across datasets.
#'
#' @param seurat_obj A `Seurat` object with QC metrics computed via
#'   [compute_qc_metrics()].
#'
#' @return A named list of `ggplot` objects: `violin`, `counts_vs_features`,
#'   and `counts_vs_mt`.
#' @export
plot_qc_metrics <- function(seurat_obj) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plot_qc_metrics().")
  }

  violin <- Seurat::VlnPlot(
    seurat_obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    ncol = 3,
    pt.size = 0
  )

  counts_vs_features <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
  )

  counts_vs_mt <- Seurat::FeatureScatter(
    seurat_obj,
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
  )

  list(
    violin = violin,
    counts_vs_features = counts_vs_features,
    counts_vs_mt = counts_vs_mt
  )
}

#' Filter low-quality cells
#'
#' Subsets a Seurat object on gene count per cell, total RNA counts, and
#' mitochondrial percentage. Thresholds must be chosen per-dataset (see
#' [plot_qc_metrics()]) -- there is no universally correct cutoff.
#'
#' @param seurat_obj A `Seurat` object with QC metrics computed via
#'   [compute_qc_metrics()].
#' @param min_genes Minimum `nFeature_RNA` to keep a cell.
#' @param max_genes Maximum `nFeature_RNA` to keep a cell (guards against
#'   multiplets that slipped past doublet detection).
#' @param min_counts Minimum `nCount_RNA` to keep a cell.
#' @param max_mt_percent Maximum `percent.mt` to keep a cell.
#'
#' @return The subset `Seurat` object.
#' @export
filter_cells <- function(seurat_obj,
                          min_genes = 200,
                          max_genes = Inf,
                          min_counts = 500,
                          max_mt_percent = 10) {
  if (!methods::is(seurat_obj, "Seurat")) {
    stop("`seurat_obj` must be a Seurat object.")
  }
  if (!"percent.mt" %in% colnames(seurat_obj[[]])) {
    stop("Run compute_qc_metrics() before filter_cells().")
  }

  keep <- seurat_obj$nFeature_RNA >= min_genes &
    seurat_obj$nFeature_RNA <= max_genes &
    seurat_obj$nCount_RNA >= min_counts &
    seurat_obj$percent.mt <= max_mt_percent

  subset(seurat_obj, cells = colnames(seurat_obj)[keep])
}
