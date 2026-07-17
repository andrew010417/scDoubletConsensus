# scDoubletConsensus

This project develops a quality control tool for scRNA-seq data that uses ensemble scoring from multiple doublet detection tools rather than simple intersection-based filtering.

## Core Innovation

Instead of taking the intersection of results from different doublet detectors, the tool normalizes individual confidence scores and creates a "consensus ranking" to identify problematic cells more sophisticated. As noted in the documentation, the existing approach "can miss cells that receive very high confidence scores from individual tools but don't appear in the intersection."

## Key Components

**The workflow combines two prominent detection methods:**
- DoubletFinder (synthetic doublet nearest-neighbor comparison)
- scDblFinder (deep learning probability estimation)

**The normalization and scoring pipeline:**
The approach involves standardizing scores across tools, then applying weighted combination or rank-based consensus to identify top-confidence doublets.

## Quality Control Foundation

Before doublet detection, the pipeline filters low-quality cells using three metrics: gene count per cell, RNA abundance, and mitochondrial gene percentage. The project acknowledges that "gold standard cutoff values don't exist," requiring manual threshold determination via visualization for each dataset.

## Installation

```r
install.packages(c("Seurat", "remotes", "BiocManager"))
BiocManager::install(c("SingleCellExperiment", "scDblFinder"))
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

remotes::install_github("andrew010417/scDoubletConsensus")
```

## Usage

```r
library(scDoubletConsensus)
library(Seurat)

obj <- CreateSeuratObject(counts = my_counts_matrix, min.cells = 3, min.features = 200)

# 1. QC: compute metrics, inspect plots, choose dataset-specific thresholds
obj <- compute_qc_metrics(obj, mt_pattern = "^MT-")
plot_qc_metrics(obj)
obj <- filter_cells(obj, min_genes = 200, max_genes = 6000, min_counts = 500, max_mt_percent = 10)

# 2. Standard Seurat preprocessing (required by DoubletFinder)
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA()
obj <- FindNeighbors(obj, dims = 1:10) |> FindClusters()

# 3. Run both detectors + consensus scoring in one step
result <- run_consensus_pipeline(obj, run_qc = FALSE, PCs = 1:10, consensus_method = "rank")

head(result$consensus_df)
result$comparison$summary          # consensus vs. intersection counts
result$comparison$recovered_cells  # high-confidence doublets missed by intersection
```

See `examples/basic_workflow.R` for the full script, and `R/` for the individual
building blocks (`compute_qc_metrics()`, `filter_cells()`, `run_doubletfinder()`,
`run_scdblfinder()`, `normalize_scores()`, `consensus_rank()`,
`classify_consensus_doublets()`, `compare_to_intersection()`) if you need finer
control than `run_consensus_pipeline()` gives you.

## Development Status

Algorithm design and an initial implementation of the QC + consensus-scoring
pipeline are in place (see `R/` and `tests/testthat/`). Outstanding roadmap
items: comparative validation against plain intersection-based filtering on
real datasets, and benchmarking against ground-truth doublet labels
(e.g. cell-hashing or genotype-multiplexed datasets).
