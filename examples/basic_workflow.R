# Example end-to-end workflow for scDoubletConsensus.
#
# Requires Seurat, DoubletFinder, and scDblFinder installed:
#   install.packages(c("Seurat", "BiocManager"))
#   BiocManager::install("scDblFinder")
#   remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
#
# Replace `counts` below with your own raw count matrix (genes x cells).

devtools::load_all(".") # or: library(scDoubletConsensus)
library(Seurat)

counts <- Seurat::Read10X(data.dir = "path/to/10x/outs/filtered_feature_bc_matrix")
obj <- Seurat::CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200)

# 1. QC: compute metrics, inspect plots, then pick thresholds for this dataset.
obj <- compute_qc_metrics(obj, mt_pattern = "^MT-")
qc_plots <- plot_qc_metrics(obj)
print(qc_plots$violin)

obj <- filter_cells(
  obj,
  min_genes = 200,
  max_genes = 6000,
  min_counts = 500,
  max_mt_percent = 10
)

# 2. Standard preprocessing required by DoubletFinder.
obj <- Seurat::NormalizeData(obj)
obj <- Seurat::FindVariableFeatures(obj)
obj <- Seurat::ScaleData(obj)
obj <- Seurat::RunPCA(obj)
obj <- Seurat::FindNeighbors(obj, dims = 1:10)
obj <- Seurat::FindClusters(obj)

# 3. Run the full consensus pipeline (QC already applied above, so skip it here).
result <- run_consensus_pipeline(
  obj,
  run_qc = FALSE,
  PCs = 1:10,
  expected_doublet_rate = 0.075,
  consensus_method = "rank",
  quantile = 0.9
)

head(result$consensus_df)
result$comparison$summary
result$comparison$recovered_cells # cells the consensus catches that plain intersection misses
