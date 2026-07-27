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

# 0. 입력 확인: cell/gene 수를 자동 감지하고 10x loading table 기반 예상 multiplet rate를 출력
inspect_input(obj)

# 1. QC: 지표 계산 후 플롯을 보고 데이터셋에 맞는 threshold를 직접 결정
obj <- compute_qc_metrics(obj, mt_pattern = "^MT-")
plot_qc_metrics(obj)
obj <- filter_cells(obj, min_genes = 200, max_genes = 6000, min_counts = 500, max_mt_percent = 10)

# 2. DoubletFinder에 필요한 기본 Seurat 전처리
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA()
obj <- FindNeighbors(obj, dims = 1:10) |> FindClusters()

# 3. 두 탐지기 실행 + consensus scoring을 한 번에 수행
# PCs, expected_doublet_rate는 숫자를 직접 지정하거나 "auto"로 자동 산출 가능
result <- run_consensus_pipeline(
  obj,
  run_qc = FALSE,
  PCs = "auto",                     # choose_pcs()의 elbow 휴리스틱으로 자동 선택
  expected_doublet_rate = "auto",   # estimate_multiplet_rate(ncol(obj))로 자동 산출
  consensus_method = "rank",
  cutoff_method = "multiplet_rate"  # 아래 "Cutoff 방식" 참고
)

head(result$consensus_df)
result$comparison$summary          # consensus vs. 교집합 결과 비교
result$comparison$recovered_cells  # 교집합에서는 놓쳤지만 consensus가 찾아낸 high-confidence doublet

# 교집합 cell들이 consensus score 분포 중 어디 위치하는지 시각화
plot_consensus_overview(result$consensus_df, call_cols = c("doublet_finder_call", "scdblfinder_call"))
```

전체 예시 스크립트는 `examples/basic_workflow.R`을 참고하세요. `run_consensus_pipeline()`보다
세밀한 제어가 필요하다면 `R/` 아래의 개별 함수들(`compute_qc_metrics()`, `filter_cells()`,
`run_doubletfinder()`, `run_scdblfinder()`, `normalize_scores()`, `consensus_rank()`,
`classify_consensus_doublets()`, `compare_to_intersection()`)을 직접 사용할 수 있습니다.

## Consensus score cutoff 방식

`classify_consensus_doublets()`(및 `run_consensus_pipeline(cutoff_method = ...)`)는 세 가지
cutoff 방식을 지원합니다:

- `"legacy"` (기본값): 기존 방식대로 `top_n` 또는 `quantile`로 직접 자름.
- `"intersection_median"`: DoubletFinder와 scDblFinder 둘 다 doublet으로 부른 교집합
  cell들의 `consensus_score` **median** 값을 threshold로 사용. "교집합에 든 cell 정도의
  confidence면 doublet으로 인정"하는 기준.
- `"multiplet_rate"`: cell 수에 맞는 10x multiplet rate(`estimate_multiplet_rate()`)만큼
  상위 cell을 자름 — 예: 5,775개 cell이면 약 4.4%인 상위 ~254개.

두 방식 모두 `plot_consensus_overview()`에서 참조선으로 함께 표시되므로, cutoff를 바꿔가며
교집합 cell들이 분포상 어디 위치하는지 눈으로 비교할 수 있습니다.

## DoubletFinder / scDblFinder 가중치

기본 가중치는 5:5가 아니라 `c(pANN = 0.4, scDblFinder_score = 0.6)`입니다. scDblFinder는
여러 벤치마크 연구(예: Xi & Li, *Benchmarking Computational Doublet-Detection Methods for
Single-Cell RNA Sequencing Data*, 2021)에서 DoubletFinder 대비 전반적으로 더 우수하고
pK/pN 같은 파라미터에 덜 민감한 것으로 보고되어, 이를 문헌 기반 기본값(prior)으로 반영했습니다.
다만 이는 데이터셋마다 최적은 아닐 수 있으므로, ground-truth가 있는 경우
`optimize_weights()`로 직접 grid-search하여 데이터 기반 최적 비율을 찾을 수 있습니다.

## 벤치마킹: 실제로 얼마나 doublet을 잡아내는가

일반적인 scRNA-seq 데이터에는 cell-hashing/genotype-demultiplexing 같은 ground-truth
doublet 라벨이 없습니다. `simulate_doublet_benchmark()`는 `scDblFinder::mockDoubletSCE()`로
실제 세포쌍을 합성해 ground-truth 라벨이 있는 벤치마크 데이터를 만들고, 전체 파이프라인을
돌려 consensus call과 교집합 call을 truth와 비교합니다:

```r
bench <- benchmark_consensus(n_cells = 3000, doublet_rate = 0.1, seed = 1)
bench$metrics
#         method n_called sensitivity precision specificity   f1
#      consensus      ...         ...       ...          ... ...
#   intersection      ...         ...       ...          ... ...
```

`evaluate_against_truth()`는 truth 벡터를 인자로 받는 범용 함수이므로, 추후 실제
cell-hashing/genotype 라벨이 있는 데이터셋이 생기면 그대로 재사용해 검증할 수 있습니다.

## 개발 현황

QC 및 consensus scoring 파이프라인, cutoff 전략, 입력 자동 감지(`inspect_input()`,
`choose_pcs()`, `estimate_multiplet_rate()`), 합성 ground-truth 벤치마킹
(`simulate_doublet_benchmark()`, `evaluate_against_truth()`, `benchmark_consensus()`,
`optimize_weights()`), 분포 시각화(`plot_consensus_overview()`)까지 구현이 완료되었습니다
(`R/`, `tests/testthat/` 참고). 남은 로드맵 항목: 실제 cell-hashing/genotype-multiplexing
기반 ground-truth 데이터셋을 이용한 검증(합성 벤치마크는 진짜 생물학적 doublet의 근사치일 뿐).
