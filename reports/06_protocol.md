# Stage 6 — Protocol Comparison (3x3 Grid): Does the Evaluation Split Inflate Performance?

Model: XGBoost (Morgan FP 2048-bit + RDKit descriptors), 5-fold CV. Three objectives (pooled regression / within-source pairwise ranking / within-source graded-relevance LambdaRank, ported from a sibling project — see scripts/lambdarank.py) x three splits (random, Butina-cluster scaffold-disjoint, leave-one-source-out). This is a from-scratch ablation built for this comparison, not a reproduction of any published architecture — the point is the *relative* change across cells.

## `LNPDB_vitro`

- Unit: distinct (source, lipid) pairs, n = 11,361, features = Morgan(2048) + RDKit descriptors(171) = 2219 dims

- LambdaRank relevance grades: hit_frac=0.1 -> 1,135 'hits' (rel=3) of 11,361 lipids; level counts [3404, 3400, 3422, 1135]

- Butina clusters (Tanimoto cutoff=0.4): 275

**pooled_regression / random_kfold**: pooled Spearman = **0.521**, within-experiment Spearman = **0.499** (n_sources=30), n_test_total=11361
  - top-k vs random (row-fragmented, partial admissibility): NDCG@k_e=0.635 (random=0.339); hit_rate@5=0.373 (random=0.109, EF@5=3.42x); hit-admissible fragments=137/145
**pooled_regression / butina_scaffold_disjoint_kfold**: pooled Spearman = **0.289**, within-experiment Spearman = **0.289** (n_sources=30), n_test_total=11361
  - top-k vs random (row-fragmented, partial admissibility): NDCG@k_e=0.529 (random=0.379); hit_rate@5=0.238 (random=0.097, EF@5=2.45x); hit-admissible fragments=88/101
**pooled_regression / leave_one_source_out**: pooled Spearman = **0.093**, within-experiment Spearman = **0.112** (n_sources=30), n_test_total=11361
  - LOSO folds at-or-below permutation-null chance level: **22/34**
  - top-k vs random (intact experiments): NDCG@k_e=0.352 (random=0.284); hit_rate@5=0.213 (random=0.101, EF@5=2.12x); hit-admissible fragments=30/30
**within_source_pairwise / random_kfold**: pooled Spearman = **0.289**, within-experiment Spearman = **0.320** (n_sources=30), n_test_total=11361
  - top-k vs random (row-fragmented, partial admissibility): NDCG@k_e=0.689 (random=0.339); hit_rate@5=0.435 (random=0.109, EF@5=3.98x); hit-admissible fragments=137/145
**within_source_pairwise / butina_scaffold_disjoint_kfold**: pooled Spearman = **0.227**, within-experiment Spearman = **0.247** (n_sources=30), n_test_total=11361
  - top-k vs random (row-fragmented, partial admissibility): NDCG@k_e=0.542 (random=0.379); hit_rate@5=0.269 (random=0.097, EF@5=2.76x); hit-admissible fragments=88/101
**within_source_pairwise / leave_one_source_out**: pooled Spearman = **0.158**, within-experiment Spearman = **0.190** (n_sources=30), n_test_total=11361
  - LOSO folds at-or-below permutation-null chance level: **17/34**
  - top-k vs random (intact experiments): NDCG@k_e=0.404 (random=0.284); hit_rate@5=0.253 (random=0.101, EF@5=2.52x); hit-admissible fragments=30/30
**within_source_lambdarank / random_kfold**: pooled Spearman = **0.311**, within-experiment Spearman = **0.316** (n_sources=30), n_test_total=11361
  - top-k vs random (row-fragmented, partial admissibility): NDCG@k_e=0.667 (random=0.339); hit_rate@5=0.386 (random=0.109, EF@5=3.54x); hit-admissible fragments=137/145
**within_source_lambdarank / butina_scaffold_disjoint_kfold**: pooled Spearman = **0.226**, within-experiment Spearman = **0.228** (n_sources=30), n_test_total=11361
  - top-k vs random (row-fragmented, partial admissibility): NDCG@k_e=0.515 (random=0.379); hit_rate@5=0.222 (random=0.097, EF@5=2.28x); hit-admissible fragments=88/101
**within_source_lambdarank / leave_one_source_out**: pooled Spearman = **0.130**, within-experiment Spearman = **0.147** (n_sources=30), n_test_total=11361
  - LOSO folds at-or-below permutation-null chance level: **21/34**
  - top-k vs random (intact experiments): NDCG@k_e=0.379 (random=0.284); hit_rate@5=0.220 (random=0.101, EF@5=2.19x); hit-admissible fragments=30/30


Figure: `figures/stage6_protocol_grid_LNPDB_vitro.png` (+ .svg)

## `LNPDB_vivo`

- Unit: distinct (source, lipid) pairs, n = 1,439, features = Morgan(2048) + RDKit descriptors(167) = 2215 dims

- LambdaRank relevance grades: hit_frac=0.1 -> 143 'hits' (rel=3) of 1,439 lipids; level counts [430, 427, 439, 143]

- Butina clusters (Tanimoto cutoff=0.4): 60

**pooled_regression / random_kfold**: pooled Spearman = **0.646**, within-experiment Spearman = **0.587** (n_sources=15), n_test_total=1439
  - top-k vs random (row-fragmented, partial admissibility): NDCG@k_e=0.779 (random=0.513); hit_rate@5=0.354 (random=0.152, EF@5=2.33x); hit-admissible fragments=49/71
**pooled_regression / butina_scaffold_disjoint_kfold**: pooled Spearman = **0.386**, within-experiment Spearman = **0.377** (n_sources=15), n_test_total=1439
  - top-k vs random (row-fragmented, partial admissibility): NDCG@k_e=0.604 (random=0.501); hit_rate@5=0.187 (random=0.143, EF@5=1.31x); hit-admissible fragments=24/33
**pooled_regression / leave_one_source_out**: pooled Spearman = **0.233**, within-experiment Spearman = **0.165** (n_sources=15), n_test_total=1439
  - LOSO folds at-or-below permutation-null chance level: **12/17**
  - top-k vs random (intact experiments): NDCG@k_e=0.367 (random=0.311); hit_rate@5=0.093 (random=0.098, EF@5=0.96x); hit-admissible fragments=15/15
**within_source_pairwise / random_kfold**: pooled Spearman = **0.506**, within-experiment Spearman = **0.553** (n_sources=15), n_test_total=1439
  - top-k vs random (row-fragmented, partial admissibility): NDCG@k_e=0.770 (random=0.513); hit_rate@5=0.333 (random=0.152, EF@5=2.19x); hit-admissible fragments=49/71
**within_source_pairwise / butina_scaffold_disjoint_kfold**: pooled Spearman = **0.319**, within-experiment Spearman = **0.356** (n_sources=15), n_test_total=1439
  - top-k vs random (row-fragmented, partial admissibility): NDCG@k_e=0.585 (random=0.501); hit_rate@5=0.244 (random=0.143, EF@5=1.71x); hit-admissible fragments=24/33
**within_source_pairwise / leave_one_source_out**: pooled Spearman = **0.090**, within-experiment Spearman = **0.061** (n_sources=15), n_test_total=1439
  - LOSO folds at-or-below permutation-null chance level: **11/17**
  - top-k vs random (intact experiments): NDCG@k_e=0.423 (random=0.311); hit_rate@5=0.107 (random=0.098, EF@5=1.09x); hit-admissible fragments=15/15
**within_source_lambdarank / random_kfold**: pooled Spearman = **0.584**, within-experiment Spearman = **0.564** (n_sources=15), n_test_total=1439
  - top-k vs random (row-fragmented, partial admissibility): NDCG@k_e=0.761 (random=0.513); hit_rate@5=0.354 (random=0.152, EF@5=2.33x); hit-admissible fragments=49/71
**within_source_lambdarank / butina_scaffold_disjoint_kfold**: pooled Spearman = **0.197**, within-experiment Spearman = **0.259** (n_sources=15), n_test_total=1439
  - top-k vs random (row-fragmented, partial admissibility): NDCG@k_e=0.650 (random=0.501); hit_rate@5=0.254 (random=0.143, EF@5=1.77x); hit-admissible fragments=24/33
**within_source_lambdarank / leave_one_source_out**: pooled Spearman = **0.021**, within-experiment Spearman = **-0.015** (n_sources=15), n_test_total=1439
  - LOSO folds at-or-below permutation-null chance level: **13/17**
  - top-k vs random (intact experiments): NDCG@k_e=0.337 (random=0.311); hit_rate@5=0.107 (random=0.098, EF@5=1.09x); hit-admissible fragments=15/15
  - LambdaRank dropped 1 distinct source(s) across all folds (IDCG=0 or single relevance level present in that training subset): [('SB_2024', 'idcg=0')]


Figure: `figures/stage6_protocol_grid_LNPDB_vivo.png` (+ .svg)
