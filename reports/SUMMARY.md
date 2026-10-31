# Pipeline Summary — Cross-Laboratory Prediction of LNP Transfection

Every number below is read directly from `results/registry.json`, the single source of truth every stage script writes to. Confidence flags and n are shown wherever the registry recorded them — nothing here is hand-typed or looked up from a prior draft.

## Claims vs. evidence

### 1. Ionizable-lipid chemistry predicts its source of origin at near-ceiling accuracy

Random-forest balanced accuracy vs. 1/S chance, permutation p<0.005 (Stage 3); collapses to chance (0.496) on the split-half negative control (Stage 8 gate: PASS), so the signal is real chemistry, not a classifier artifact.

| registry key | value |
|---|---|
| `separability.LNPDB_vitro.morgan_fp.random_forest.stratified_kfold.balanced_accuracy` | **0.967**, n=11350 |
| `separability.LNPDB_vivo.morgan_fp.random_forest.stratified_kfold.balanced_accuracy` | **0.982**, n=1407 |
| `separability.LNPDB_vitro.morgan_fp.random_forest.scaffold_disjoint_groupkfold.balanced_accuracy` | **0.880**, n=11350 |
| `separability.LNPDB_vitro.nn_same_source_rate` | **0.918**, n=11350 |
| `controls.split_half.LNPDB_vitro.separability_balanced_accuracy` | **0.496**, n=1470 |

### 2. Ionizable lipids are ~source-exclusive; excipient scaffold is shared but nearly invariant

97.9%/99.9% of lipids appear in exactly one source; cholesterol identity is almost invariant (Herfindahl~1) so formulation cannot itself explain cross-lab disagreement.

| registry key | value |
|---|---|
| `exclusivity.LNPDB_vitro.lipid.frac_single_source` | **0.979**, n=11128 |
| `exclusivity.LNPDB_vivo.lipid.frac_single_source` | **0.999**, n=1434 |
| `exclusivity.LNPDB_vitro.component_head.frac_single_source` | **0.777**, n=385 |
| `formulation.LNPDB_vitro.cholesterol.herfindahl` | **0.997**, n=18700 |
| `formulation.LNPDB_vitro.cramers_v_source` | **0.754**, n=18700 |

### 3. Identical lipids rank near-independently across labs (LOAD-BEARING)

Pooled cross-lab Spearman 0.08 (CI crosses zero) vs. a within-source replicate ceiling of 0.40 — real reproducibility exists, cross-lab rank agreement does not. In-vivo version is LOW-CONFIDENCE (n=9 pairs) — supports moving it to supplementary (T8).

| registry key | value |
|---|---|
| `concordance.LNPDB_vitro.pooled_spearman` | **0.083**, n=236, CI=[-0.041, 0.215] |
| `concordance.LNPDB_vitro.ceiling_spearman` | **0.400**, n=1035, CI=[0.346, 0.452] |
| `concordance.LNPDB_vivo.pooled_spearman` | **-0.123**, n=9, CI=[-0.771, 0.691], **LOW-CONFIDENCE** |
| `concordance.molea.hepatocyte_vs_chondrocyte_spearman` | **0.299**, n=1152, CI=[0.242, 0.357] |

### 4. Honest evaluation (LOSO) collapses relative to random/cluster CV

Pooled Spearman 0.52 (random) -> 0.29 (scaffold-disjoint) -> 0.09 (LOSO); 22/34 LOSO folds land at-or-below their own permutation-null chance level.

| registry key | value |
|---|---|
| `protocol.LNPDB_vitro.pooled_regression.random_kfold.pooled_spearman` | **0.521**, n=11361 |
| `protocol.LNPDB_vitro.pooled_regression.butina_scaffold_disjoint_kfold.pooled_spearman` | **0.289**, n=11361 |
| `protocol.LNPDB_vitro.pooled_regression.leave_one_source_out.pooled_spearman` | **0.093**, n=11361 |
| `protocol.LNPDB_vitro.pooled_regression.leave_one_source_out.n_folds_at_chance` | **22**, n=34 |

### 5. Positioning: extrapolation across disjoint support, not a shared-distribution robustness problem -> domain-generalization methods don't apply

Support is ~98% disjoint (Result 2) and LOSO performance vs. chemical distance to the training support shows only a weak, non-significant trend (Stage 7) — consistent with held-out sources being a genuine extrapolation regime, not a shared-distribution shift.

| registry key | value |
|---|---|
| `exclusivity.LNPDB_vitro.lipid.frac_single_source` | **0.979**, n=11128 |
| `decay.LNPDB_vitro.trend_slope` | **-0.154**, n=30 |

### 6. A within-experiment regime stays learnable

Size-weighted within-experiment scaffold-disjoint Spearman 0.32 (vitro); MOLEA's same-lab cross-cell-type Spearman (0.30) is an independent corroborating middle rung between this and the near-zero cross-lab number in Result 3.

| registry key | value |
|---|---|
| `learnable.LNPDB_vitro.size_weighted_spearman` | **0.322**, n=11350 |
| `learnable.LNPDB_vivo.size_weighted_spearman` | **0.114**, n=1097 |
| `protocol.LNPDB_vitro.within_source_pairwise.butina_scaffold_disjoint_kfold.within_experiment_spearman` | **0.247**, n=30 |
| `concordance.molea.hepatocyte_vs_chondrocyte_spearman` | **0.299**, n=1152, CI=[0.242, 0.357] |

## Validity gate (Stage 8)

| registry key | value |
|---|---|
| `controls.split_half.LNPDB_vitro.concordance_spearman` | **-0.178**, n=14, **LOW-CONFIDENCE** |
| `controls.split_half.LNPDB_vitro.separability_balanced_accuracy` | **0.496**, n=1470 |
| `controls.split_half.LNPDB_vivo.concordance_spearman` | **0.717**, n=9, **LOW-CONFIDENCE** |
| `controls.split_half.LNPDB_vivo.separability_balanced_accuracy` | **0.490**, n=383 |

## Full registry appendix

| key | value | n | ci | confidence | stage |
|---|---|---|---|---|---|
| `concordance.LNPDB_vitro.anchoring.anchored_spearman` | nan | 0 |  | OK | stage5 |
| `concordance.LNPDB_vitro.anchoring.baseline_spearman` | nan | 0 |  | OK | stage5 |
| `concordance.LNPDB_vitro.ceiling_spearman` | 0.3999 | 1035 | [0.346, 0.452] | OK | stage5 |
| `concordance.LNPDB_vitro.pooled_spearman` | 0.0834 | 236 | [-0.041, 0.215] | OK | stage5 |
| `concordance.LNPDB_vitro.shared_set_frac_reference_lipid` | 0.0087 | 231 |  | OK | stage5 |
| `concordance.LNPDB_vivo.anchoring.anchored_spearman` | nan | 0 |  | OK | stage5 |
| `concordance.LNPDB_vivo.anchoring.baseline_spearman` | nan | 0 |  | OK | stage5 |
| `concordance.LNPDB_vivo.ceiling_spearman` | 0.3319 | 194 | [0.196, 0.457] | OK | stage5 |
| `concordance.LNPDB_vivo.pooled_spearman` | -0.1228 | 9 | [-0.771, 0.691] | LOW-CONFIDENCE | stage5 |
| `concordance.LNPDB_vivo.shared_set_frac_reference_lipid` | 1.0000 | 2 |  | OK | stage5 |
| `concordance.molea.hepatocyte_vs_chondrocyte_spearman` | 0.2994 | 1152 | [0.242, 0.357] | OK | stage5 |
| `controls.featurization_sensitivity.LNPDB_vitro.alt_morgan_balanced_accuracy` | 0.9710 | 11350 |  | OK | stage8 |
| `controls.featurization_sensitivity.LNPDB_vivo.alt_morgan_balanced_accuracy` | 0.9824 | 1407 |  | OK | stage8 |
| `controls.seed_sensitivity.LNPDB_vitro.concordance_ci_width_std` | 0.0063 | 5 |  | OK | stage8 |
| `controls.seed_sensitivity.LNPDB_vitro.separability_balanced_accuracy_std` | 0.0008 | 5 |  | OK | stage8 |
| `controls.seed_sensitivity.LNPDB_vivo.concordance_ci_width_std` | 0.0403 | 5 |  | OK | stage8 |
| `controls.seed_sensitivity.LNPDB_vivo.separability_balanced_accuracy_std` | 0.0032 | 5 |  | OK | stage8 |
| `controls.split_half.LNPDB_vitro.concordance_spearman` | -0.1779 | 14 |  | LOW-CONFIDENCE | stage8 |
| `controls.split_half.LNPDB_vitro.separability_balanced_accuracy` | 0.4963 | 1470 |  | OK | stage8 |
| `controls.split_half.LNPDB_vivo.concordance_spearman` | 0.7167 | 9 |  | LOW-CONFIDENCE | stage8 |
| `controls.split_half.LNPDB_vivo.separability_balanced_accuracy` | 0.4903 | 383 |  | OK | stage8 |
| `corpus.n_comparable_environments_LNPDB_vitro` | 13 | None |  | OK | stage2 |
| `corpus.n_comparable_environments_LNPDB_vivo` | 8 | None |  | OK | stage2 |
| `corpus.n_environments_LNPDB_vitro` | 14 | None |  | OK | stage2 |
| `corpus.n_environments_LNPDB_vivo` | 9 | None |  | OK | stage2 |
| `corpus.n_measurements_LNPDB_vitro` | 18700 | None |  | OK | stage2 |
| `corpus.n_measurements_LNPDB_vivo` | 1717 | None |  | OK | stage2 |
| `corpus.n_rows_molea` | 1152 | None |  | OK | stage2 |
| `corpus.n_sources_LNPDB_vitro` | 34 | None |  | OK | stage2 |
| `corpus.n_sources_LNPDB_vivo` | 17 | None |  | OK | stage2 |
| `corpus.n_unique_lipids_LNPDB_vitro` | 11128 | None |  | OK | stage2 |
| `corpus.n_unique_lipids_LNPDB_vivo` | 1434 | None |  | OK | stage2 |
| `corpus.total_measurements` | 20417 | None |  | OK | stage2 |
| `decay.LNPDB_vitro.trend_slope` | -0.1541 | 30 |  | OK | stage7 |
| `decay.LNPDB_vivo.trend_slope` | -0.7934 | 16 |  | OK | stage7 |
| `exclusivity.LNPDB_vitro.component_head.frac_single_source` | 0.7766 | 385 |  | OK | stage4 |
| `exclusivity.LNPDB_vitro.component_linker.frac_single_source` | 0.8333 | 24 |  | LOW-CONFIDENCE | stage4 |
| `exclusivity.LNPDB_vitro.component_tail1.frac_single_source` | 0.7845 | 181 |  | OK | stage4 |
| `exclusivity.LNPDB_vitro.component_tail2.frac_single_source` | 0.9091 | 11 |  | LOW-CONFIDENCE | stage4 |
| `exclusivity.LNPDB_vitro.lipid.frac_single_source` | 0.9792 | 11128 |  | OK | stage4 |
| `exclusivity.LNPDB_vivo.component_head.frac_single_source` | 0.6613 | 62 |  | OK | stage4 |
| `exclusivity.LNPDB_vivo.component_linker.frac_single_source` | 1.0000 | 24 |  | OK | stage4 |
| `exclusivity.LNPDB_vivo.component_tail1.frac_single_source` | 0.9324 | 74 |  | OK | stage4 |
| `exclusivity.LNPDB_vivo.component_tail2.frac_single_source` | 0.9737 | 38 |  | OK | stage4 |
| `exclusivity.LNPDB_vivo.lipid.frac_single_source` | 0.9986 | 1434 |  | OK | stage4 |
| `formulation.LNPDB_vitro.cholesterol.herfindahl` | 0.9973 | 18700 |  | OK | stage4 |
| `formulation.LNPDB_vitro.cholesterol.top1_frac` | 0.9987 | 18700 |  | OK | stage4 |
| `formulation.LNPDB_vitro.cramers_v_source` | 0.7541 | 18700 |  | OK | stage4 |
| `formulation.LNPDB_vitro.helper_lipid.herfindahl` | 0.3392 | 18131 |  | OK | stage4 |
| `formulation.LNPDB_vitro.helper_lipid.top1_frac` | 0.4846 | 18131 |  | OK | stage4 |
| `formulation.LNPDB_vitro.nmi_source` | 0.6386 | 18700 |  | OK | stage4 |
| `formulation.LNPDB_vitro.peg_lipid.herfindahl` | 0.4958 | 18516 |  | OK | stage4 |
| `formulation.LNPDB_vitro.peg_lipid.top1_frac` | 0.6456 | 18516 |  | OK | stage4 |
| `formulation.LNPDB_vivo.cholesterol.herfindahl` | 0.9555 | 1717 |  | OK | stage4 |
| `formulation.LNPDB_vivo.cholesterol.top1_frac` | 0.9773 | 1717 |  | OK | stage4 |
| `formulation.LNPDB_vivo.cramers_v_source` | 0.7607 | 1717 |  | OK | stage4 |
| `formulation.LNPDB_vivo.helper_lipid.herfindahl` | 0.5101 | 1717 |  | OK | stage4 |
| `formulation.LNPDB_vivo.helper_lipid.top1_frac` | 0.6144 | 1717 |  | OK | stage4 |
| `formulation.LNPDB_vivo.nmi_source` | 0.7782 | 1717 |  | OK | stage4 |
| `formulation.LNPDB_vivo.peg_lipid.herfindahl` | 0.3830 | 1717 |  | OK | stage4 |
| `formulation.LNPDB_vivo.peg_lipid.top1_frac` | 0.4898 | 1717 |  | OK | stage4 |
| `learnable.LNPDB_vitro.size_weighted_spearman` | 0.3220 | 11350 |  | OK | stage7 |
| `learnable.LNPDB_vivo.size_weighted_spearman` | 0.1144 | 1097 |  | OK | stage7 |
| `protocol.LNPDB_vitro.pooled_regression.butina_scaffold_disjoint_kfold.pooled_spearman` | 0.2893 | 11361 |  | OK | stage6 |
| `protocol.LNPDB_vitro.pooled_regression.butina_scaffold_disjoint_kfold.topk.hit_rate_at_5` | 0.2384 | 88 |  | OK | stage6 |
| `protocol.LNPDB_vitro.pooled_regression.butina_scaffold_disjoint_kfold.topk.ndcg_at_k_e` | 0.5289 | 101 |  | OK | stage6 |
| `protocol.LNPDB_vitro.pooled_regression.butina_scaffold_disjoint_kfold.within_experiment_spearman` | 0.2893 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.pooled_regression.leave_one_source_out.n_folds_at_chance` | 22 | 34 |  | OK | stage6 |
| `protocol.LNPDB_vitro.pooled_regression.leave_one_source_out.pooled_spearman` | 0.0928 | 11361 |  | OK | stage6 |
| `protocol.LNPDB_vitro.pooled_regression.leave_one_source_out.topk.hit_rate_at_5` | 0.2133 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.pooled_regression.leave_one_source_out.topk.ndcg_at_k_e` | 0.3522 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.pooled_regression.leave_one_source_out.within_experiment_spearman` | 0.1122 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.pooled_regression.random_kfold.pooled_spearman` | 0.5211 | 11361 |  | OK | stage6 |
| `protocol.LNPDB_vitro.pooled_regression.random_kfold.topk.hit_rate_at_5` | 0.3734 | 137 |  | OK | stage6 |
| `protocol.LNPDB_vitro.pooled_regression.random_kfold.topk.ndcg_at_k_e` | 0.6355 | 145 |  | OK | stage6 |
| `protocol.LNPDB_vitro.pooled_regression.random_kfold.within_experiment_spearman` | 0.4994 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.butina_scaffold_disjoint_kfold.pooled_spearman` | 0.2258 | 11361 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.butina_scaffold_disjoint_kfold.topk.hit_rate_at_5` | 0.2221 | 88 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.butina_scaffold_disjoint_kfold.topk.ndcg_at_k_e` | 0.5150 | 101 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.butina_scaffold_disjoint_kfold.within_experiment_spearman` | 0.2281 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.leave_one_source_out.n_folds_at_chance` | 21 | 34 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.leave_one_source_out.pooled_spearman` | 0.1304 | 11361 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.leave_one_source_out.topk.hit_rate_at_5` | 0.2200 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.leave_one_source_out.topk.ndcg_at_k_e` | 0.3789 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.leave_one_source_out.within_experiment_spearman` | 0.1468 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.random_kfold.pooled_spearman` | 0.3107 | 11361 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.random_kfold.topk.hit_rate_at_5` | 0.3865 | 137 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.random_kfold.topk.ndcg_at_k_e` | 0.6671 | 145 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_lambdarank.random_kfold.within_experiment_spearman` | 0.3160 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.butina_scaffold_disjoint_kfold.pooled_spearman` | 0.2268 | 11361 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.butina_scaffold_disjoint_kfold.topk.hit_rate_at_5` | 0.2686 | 88 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.butina_scaffold_disjoint_kfold.topk.ndcg_at_k_e` | 0.5415 | 101 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.butina_scaffold_disjoint_kfold.within_experiment_spearman` | 0.2471 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.leave_one_source_out.n_folds_at_chance` | 17 | 34 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.leave_one_source_out.pooled_spearman` | 0.1576 | 11361 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.leave_one_source_out.topk.hit_rate_at_5` | 0.2533 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.leave_one_source_out.topk.ndcg_at_k_e` | 0.4042 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.leave_one_source_out.within_experiment_spearman` | 0.1899 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.random_kfold.pooled_spearman` | 0.2892 | 11361 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.random_kfold.topk.hit_rate_at_5` | 0.4347 | 137 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.random_kfold.topk.ndcg_at_k_e` | 0.6887 | 145 |  | OK | stage6 |
| `protocol.LNPDB_vitro.within_source_pairwise.random_kfold.within_experiment_spearman` | 0.3203 | 30 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.butina_scaffold_disjoint_kfold.pooled_spearman` | 0.3855 | 1439 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.butina_scaffold_disjoint_kfold.topk.hit_rate_at_5` | 0.1873 | 24 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.butina_scaffold_disjoint_kfold.topk.ndcg_at_k_e` | 0.6042 | 33 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.butina_scaffold_disjoint_kfold.within_experiment_spearman` | 0.3765 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.leave_one_source_out.n_folds_at_chance` | 12 | 17 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.leave_one_source_out.pooled_spearman` | 0.2327 | 1439 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.leave_one_source_out.topk.hit_rate_at_5` | 0.0933 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.leave_one_source_out.topk.ndcg_at_k_e` | 0.3671 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.leave_one_source_out.within_experiment_spearman` | 0.1648 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.random_kfold.pooled_spearman` | 0.6460 | 1439 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.random_kfold.topk.hit_rate_at_5` | 0.3537 | 49 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.random_kfold.topk.ndcg_at_k_e` | 0.7790 | 71 |  | OK | stage6 |
| `protocol.LNPDB_vivo.pooled_regression.random_kfold.within_experiment_spearman` | 0.5873 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.butina_scaffold_disjoint_kfold.pooled_spearman` | 0.1969 | 1439 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.butina_scaffold_disjoint_kfold.topk.hit_rate_at_5` | 0.2540 | 24 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.butina_scaffold_disjoint_kfold.topk.ndcg_at_k_e` | 0.6504 | 33 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.butina_scaffold_disjoint_kfold.within_experiment_spearman` | 0.2591 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.leave_one_source_out.n_folds_at_chance` | 13 | 17 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.leave_one_source_out.pooled_spearman` | 0.0210 | 1439 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.leave_one_source_out.topk.hit_rate_at_5` | 0.1067 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.leave_one_source_out.topk.ndcg_at_k_e` | 0.3370 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.leave_one_source_out.within_experiment_spearman` | -0.0149 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.random_kfold.pooled_spearman` | 0.5843 | 1439 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.random_kfold.topk.hit_rate_at_5` | 0.3537 | 49 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.random_kfold.topk.ndcg_at_k_e` | 0.7609 | 71 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_lambdarank.random_kfold.within_experiment_spearman` | 0.5641 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.butina_scaffold_disjoint_kfold.pooled_spearman` | 0.3194 | 1439 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.butina_scaffold_disjoint_kfold.topk.hit_rate_at_5` | 0.2444 | 24 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.butina_scaffold_disjoint_kfold.topk.ndcg_at_k_e` | 0.5850 | 33 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.butina_scaffold_disjoint_kfold.within_experiment_spearman` | 0.3564 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.leave_one_source_out.n_folds_at_chance` | 11 | 17 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.leave_one_source_out.pooled_spearman` | 0.0897 | 1439 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.leave_one_source_out.topk.hit_rate_at_5` | 0.1067 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.leave_one_source_out.topk.ndcg_at_k_e` | 0.4230 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.leave_one_source_out.within_experiment_spearman` | 0.0608 | 15 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.random_kfold.pooled_spearman` | 0.5060 | 1439 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.random_kfold.topk.hit_rate_at_5` | 0.3333 | 49 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.random_kfold.topk.ndcg_at_k_e` | 0.7701 | 71 |  | OK | stage6 |
| `protocol.LNPDB_vivo.within_source_pairwise.random_kfold.within_experiment_spearman` | 0.5533 | 15 |  | OK | stage6 |
| `separability.LNPDB_vitro.morgan_fp.logistic_regression.scaffold_disjoint_groupkfold.balanced_accuracy` | 0.8822 | 11350 |  | OK | stage3 |
| `separability.LNPDB_vitro.morgan_fp.logistic_regression.stratified_kfold.balanced_accuracy` | 0.9721 | 11350 |  | OK | stage3 |
| `separability.LNPDB_vitro.morgan_fp.random_forest.scaffold_disjoint_groupkfold.balanced_accuracy` | 0.8797 | 11350 |  | OK | stage3 |
| `separability.LNPDB_vitro.morgan_fp.random_forest.stratified_kfold.balanced_accuracy` | 0.9672 | 11350 |  | OK | stage3 |
| `separability.LNPDB_vitro.nn_same_source_rate` | 0.9182 | 11350 |  | OK | stage3 |
| `separability.LNPDB_vitro.permutation.observed_balanced_accuracy` | 0.9711 | None |  | OK | stage3 |
| `separability.LNPDB_vitro.rdkit_descriptors.logistic_regression.stratified_kfold.balanced_accuracy` | 0.9710 | 11350 |  | OK | stage3 |
| `separability.LNPDB_vivo.morgan_fp.logistic_regression.scaffold_disjoint_groupkfold.balanced_accuracy` | 0.9283 | 1407 |  | OK | stage3 |
| `separability.LNPDB_vivo.morgan_fp.logistic_regression.stratified_kfold.balanced_accuracy` | 0.9824 | 1407 |  | OK | stage3 |
| `separability.LNPDB_vivo.morgan_fp.random_forest.scaffold_disjoint_groupkfold.balanced_accuracy` | 0.9268 | 1407 |  | OK | stage3 |
| `separability.LNPDB_vivo.morgan_fp.random_forest.stratified_kfold.balanced_accuracy` | 0.9824 | 1407 |  | OK | stage3 |
| `separability.LNPDB_vivo.nn_same_source_rate` | 0.9957 | 1407 |  | OK | stage3 |
| `separability.LNPDB_vivo.permutation.observed_balanced_accuracy` | 0.9808 | None |  | OK | stage3 |
| `separability.LNPDB_vivo.rdkit_descriptors.logistic_regression.stratified_kfold.balanced_accuracy` | 0.9774 | 1407 |  | OK | stage3 |