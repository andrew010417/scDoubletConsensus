# scDoubletConsensus

이 프로젝트는 여러 doublet(이중세포) 탐지 도구의 결과를 단순 교집합으로 필터링하는 대신, 앙상블 스코어링을 이용해 scRNA-seq 데이터의 품질을 관리하는 도구를 개발합니다.

## 핵심 아이디어

서로 다른 doublet 탐지기의 결과를 교집합으로 취하는 대신, 각 도구의 confidence score를 정규화한 뒤 "consensus ranking(합의 순위)"을 만들어 문제가 되는 세포를 더 정교하게 식별합니다. 기존 문서에서 언급하듯, 교집합 기반 접근은 "개별 도구에서는 매우 높은 confidence score를 받았지만 교집합에는 나타나지 않는 세포를 놓칠 수 있다"는 한계가 있습니다.

## 주요 구성 요소

**워크플로우는 대표적인 두 가지 탐지 방법을 결합합니다:**
- DoubletFinder (synthetic doublet를 이용한 nearest-neighbor 비교)
- scDblFinder (딥러닝 기반 확률 추정)

**정규화 및 스코어링 파이프라인:**
각 도구의 score를 표준화한 뒤, 가중 결합(weighted combination) 또는 순위 기반 합의(rank-based consensus)를 적용하여 confidence가 높은 doublet을 식별합니다.

## 품질 관리(QC) 기반 작업

Doublet 탐지에 앞서, 파이프라인은 세포당 유전자 수, RNA 총량, 미토콘드리아 유전자 비율이라는 세 가지 지표를 이용해 저품질 세포를 필터링합니다. "gold standard 컷오프 값은 존재하지 않는다"는 점을 인정하며, 데이터셋마다 시각화를 통해 직접 threshold를 정해야 합니다.

## 설치

```r
install.packages(c("Seurat", "remotes", "BiocManager"))
BiocManager::install(c("SingleCellExperiment", "scDblFinder"))
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")

remotes::install_github("andrew010417/scDoubletConsensus")
```

## 사용법

```r
library(scDoubletConsensus)
library(Seurat)

obj <- CreateSeuratObject(counts = my_counts_matrix, min.cells = 3, min.features = 200)

# 1. QC: 지표 계산 후 플롯을 보고 데이터셋에 맞는 threshold를 직접 결정
obj <- compute_qc_metrics(obj, mt_pattern = "^MT-")
plot_qc_metrics(obj)
obj <- filter_cells(obj, min_genes = 200, max_genes = 6000, min_counts = 500, max_mt_percent = 10)

# 2. DoubletFinder에 필요한 기본 Seurat 전처리
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA()
obj <- FindNeighbors(obj, dims = 1:10) |> FindClusters()

# 3. 두 탐지기 실행 + consensus scoring을 한 번에 수행
result <- run_consensus_pipeline(obj, run_qc = FALSE, PCs = 1:10, consensus_method = "rank")

head(result$consensus_df)
result$comparison$summary          # consensus vs. 교집합 결과 비교
result$comparison$recovered_cells  # 교집합에서는 놓쳤지만 consensus가 찾아낸 high-confidence doublet
```

전체 예시 스크립트는 `examples/basic_workflow.R`을 참고하세요. `run_consensus_pipeline()`보다
세밀한 제어가 필요하다면 `R/` 아래의 개별 함수들(`compute_qc_metrics()`, `filter_cells()`,
`run_doubletfinder()`, `run_scdblfinder()`, `normalize_scores()`, `consensus_rank()`,
`classify_consensus_doublets()`, `compare_to_intersection()`)을 직접 사용할 수 있습니다.

## 개발 현황

QC 및 consensus scoring 파이프라인의 알고리즘 설계와 초기 구현이 완료되었습니다 (`R/`,
`tests/testthat/` 참고). 남은 로드맵 항목: 실제 데이터셋에서 단순 교집합 필터링 대비
비교 검증(comparative validation), 그리고 ground-truth doublet 라벨(예: cell-hashing이나
genotype-multiplexing 기반 데이터셋)을 이용한 벤치마킹.
