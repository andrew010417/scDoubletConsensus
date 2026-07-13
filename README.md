# scDoubletConsensus

여러 doublet detection tool의 **score를 앙상블(consensus)** 하여, 단순 교집합보다 더 정교하게 doublet을 판정하는 scRNA-seq QC 도구 개발 프로젝트.

## Background

Single-cell RNA-seq 데이터는 droplet 기반 실험 특성상 하나의 droplet에 세포 2개 이상이 함께 캡슐화되는 **doublet/triplet**이 섞여 들어온다. 이를 걸러내지 않으면 존재하지 않는 세포 타입이나 잘못된 발현 패턴이 downstream 분석(clustering, DEG 등)에 오염된다.

기존 워크플로우:

```
FASTQ → Cell Ranger → raw/filtered count matrix + QC report
                              ↓
                    QC trimming (기준: nFeature_RNA, nCount_RNA, percent.mt)
                              ↓
                    Doublet detection (DoubletFinder, scDblFinder 등)
```

Doublet detection에는 gold standard가 없어 데이터셋마다 FeatureScatter로 직접 cutoff를 잡아야 하고, 대표적인 tool들도 서로 다른 원리로 동작한다:

- **DoubletFinder**: 인공적으로 합성한 synthetic doublet을 각 세포의 최근접 이웃 비율(pANN)로 비교해 classification
- **scDblFinder**: deep learning 기반 모델로 doublet 확률을 추정

## Problem

현재 관행은 여러 tool의 결과를 **교집합(intersection)**으로 doublet을 확정하는 경우가 많다. 하지만 이 방식은:

1. 교집합에 들지 못해도 개별 tool에서 매우 높은 confidence score를 받은 세포를 놓칠 수 있음
2. 두 tool의 score 분포/스케일이 다른데 이를 반영하지 못함

## Goal

DoubletFinder와 scDblFinder(및 확장 가능한 다른 tool)의 **개별 confidence score를 정규화 후 종합 스코어링**하여, 교집합 여부와 무관하게 top-confidence 세포를 doublet으로 판정하는 consensus tool을 개발한다.

```
Tool A score ─┐
              ├─→ 정규화 → 가중합/랭크 기반 consensus score → top-N을 doublet으로 판정
Tool B score ─┘
```

## Repository Structure

```
scDoubletConsensus/
├── data/
│   ├── raw/            # Cell Ranger 산출물 (raw/filtered count matrix)
│   └── processed/       # QC trimming 이후 처리된 데이터
├── scripts/              # 분석 및 scoring 파이프라인 스크립트
├── notebooks/            # 탐색적 분석, 실험용 노트북
├── results/
│   ├── figures/          # QC plot, score 분포 등 시각화 결과
│   └── tables/           # consensus score 결과표
└── docs/                 # 참고자료, 발표자료 요약, 알고리즘 노트
```

## QC Trimming Criteria (선행 단계)

Doublet detection 이전에 아래 기준으로 저품질 세포를 우선 제거한다:

| 기준 | 설명 |
|---|---|
| `nFeature_RNA` | 세포당 검출된 유전자 수 |
| `nCount_RNA` | 세포당 총 RNA count |
| `percent.mt` | 미토콘드리아 유전자 비율 (세포 사멸 시 세포막·미토콘드리아막이 함께 파괴되어 mt-RNA 비율이 비정상적으로 높아짐) |

걸러내는 대상: low quality/dying cells, empty droplets, doublets/triplets, contaminated cells.

Gold standard cutoff가 없으므로, 각 샘플마다 FeatureScatter plot을 직접 확인하여 threshold를 결정한다.

## Roadmap

- [ ] DoubletFinder, scDblFinder GitHub 소스 분석 — score 산출 로직 파악
- [ ] 두 tool을 동일 데이터셋에 각각 적용, raw score 추출
- [ ] Score 정규화 방법 결정 (min-max, rank-based 등)
- [ ] Consensus scoring 알고리즘 설계 및 구현
- [ ] 교집합 기반 방식 vs consensus scoring 방식 성능 비교
- [ ] 벤치마크 데이터셋(ground-truth doublet label 존재하는 데이터)으로 검증

## References

- DoubletFinder: [github.com/chris-mcginnis-ucsf/DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
- scDblFinder: [github.com/plger/scDblFinder](https://github.com/plger/scDblFinder)
