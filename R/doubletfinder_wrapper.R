#' Run DoubletFinder and extract a per-cell doublet score
#'
#' Thin wrapper around the `DoubletFinder` package (synthetic-doublet
#' nearest-neighbor comparison). Requires a normalized, dimensionality-reduced
#' Seurat object (`NormalizeData()`, `FindVariableFeatures()`, `ScaleData()`,
#' `RunPCA()` already run).
#'
#' @param seurat_obj A processed `Seurat` object.
#' @param PCs Integer vector of principal components to use.
#' @param expected_doublet_rate Expected doublet rate (e.g. `0.075` for the
#'   standard 10x ~7.5% at ~10k cells loaded).
#' @param sct Whether the object was normalized with `SCTransform()`.
#' @param seed Random seed passed through for reproducibility.
#'
#' @return A `data.frame` with columns `cell`, `pANN` (raw DoubletFinder
#'   score), and `doublet_finder_call` (`"Doublet"`/`"Singlet"`).
#' @export
run_doubletfinder <- function(seurat_obj,
                               PCs = 1:10,
                               expected_doublet_rate = 0.075,
                               sct = FALSE,
                               seed = 1) {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    stop(
      "Package 'DoubletFinder' is required. Install from GitHub via:\n",
      "  remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"
    )
  }

  set.seed(seed)

  sweep_res <- DoubletFinder::paramSweep(seurat_obj, PCs = PCs, sct = sct)
  sweep_stats <- DoubletFinder::summarizeSweep(sweep_res, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep_stats)
  optimal_pK <- as.numeric(as.character(
    bcmvn$pK[which.max(bcmvn$BCmetric)]
  ))

  homotypic_prop <- DoubletFinder::modelHomotypic(seurat_obj$seurat_clusters)
  n_exp <- round(expected_doublet_rate * ncol(seurat_obj))
  n_exp_adj <- round(n_exp * (1 - homotypic_prop))

  seurat_obj <- DoubletFinder::doubletFinder(
    seurat_obj,
    PCs = PCs,
    pN = 0.25,
    pK = optimal_pK,
    nExp = n_exp_adj,
    reuse.pANN = FALSE,
    sct = sct
  )

  meta <- seurat_obj[[]]
  pann_col <- grep("^pANN_", colnames(meta), value = TRUE)[1]
  call_col <- grep("^DF.classifications_", colnames(meta), value = TRUE)[1]

  data.frame(
    cell = colnames(seurat_obj),
    pANN = meta[[pann_col]],
    doublet_finder_call = meta[[call_col]],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
