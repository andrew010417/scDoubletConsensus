#' Run scDblFinder and extract a per-cell doublet score
#'
#' Thin wrapper around the `scDblFinder` package (deep-learning-style
#' probability estimation via simulated doublets and a k-NN classifier).
#'
#' @param seurat_obj A `Seurat` object with a raw counts assay.
#' @param assay Assay to pull raw counts from.
#' @param seed Random seed passed through for reproducibility.
#'
#' @return A `data.frame` with columns `cell`, `scDblFinder_score` (raw
#'   probability), and `scdblfinder_call` (`"doublet"`/`"singlet"`).
#' @export
run_scdblfinder <- function(seurat_obj, assay = "RNA", seed = 1) {
  if (!requireNamespace("scDblFinder", quietly = TRUE)) {
    stop(
      "Package 'scDblFinder' is required. Install via:\n",
      "  BiocManager::install('scDblFinder')"
    )
  }
  if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
    stop("Package 'SingleCellExperiment' is required.")
  }

  counts <- Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts)
  )

  set.seed(seed)
  sce <- scDblFinder::scDblFinder(sce)

  data.frame(
    cell = colnames(sce),
    scDblFinder_score = sce$scDblFinder.score,
    scdblfinder_call = sce$scDblFinder.class,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
