# Stage 3 — Source Separability from Chemistry

Featurization: Morgan/ECFP fingerprints (radius=2, 2048 bits); RDKit full descriptor set as a robustness check. Classification unit = distinct (source, canonical ionizable-lipid SMILES) pairs (see module docstring for why raw measurement rows would leak). Minimum 20 unique lipids/source to be included.

## `LNPDB_vitro`

- Sources with >= 20 unique lipids (included): **29**
- Sources excluded (too few unique lipids): 5 — AC_2025 (n=7), LZ_2025 (n=1), SP_2020 (n=1), YZ_2022 (n=1), YZ_2024 (n=1)
- Classification unit: distinct (source, lipid) pairs = **11,350**
- Chance baseline (1/S): **0.0345**

- Morgan fingerprint valid: 11,350 / 11,350

**morgan_fp / logistic_regression / stratified_kfold**: balanced_acc=0.972 (+/-0.003), top1=0.952, macro_F1=0.952, chance=0.034  [8.7s]
**morgan_fp / logistic_regression / scaffold_disjoint_groupkfold**: balanced_acc=0.882 (+/-0.077), top1=0.892, macro_F1=0.821, chance=0.034  [8.6s]
**morgan_fp / random_forest / stratified_kfold**: balanced_acc=0.967 (+/-0.003), top1=0.953, macro_F1=0.953, chance=0.034  [33.5s]
**morgan_fp / random_forest / scaffold_disjoint_groupkfold**: balanced_acc=0.880 (+/-0.089), top1=0.893, macro_F1=0.834, chance=0.034  [33.0s]

**RDKit descriptor featurization (robustness check vs. Morgan FP):**

- RDKit descriptors (171 features) / logistic_regression / stratified_kfold: balanced_acc=0.971, chance=0.034

**Same-source nearest-neighbor rate**: observed = 0.918, expected under random source assignment = 0.068 (ratio = 13.60x)

**Permutation baseline** (200 shuffles, logistic regression, single 80/20 split): observed balanced_acc=0.971, null mean=0.034 (std=0.007), 0/200 null shuffles >= observed, empirical p < 0.0050 (resolution floor of 200 shuffles)

- UMAP embedding saved to `figures/stage3_embedding_LNPDB_vitro.png` (+ .svg)

## `LNPDB_vivo`

- Sources with >= 20 unique lipids (included): **13**
- Sources excluded (too few unique lipids): 4 — LL_2025 (n=1), SB_2024 (n=2), SS_2018 (n=10), SY_2025 (n=19)
- Classification unit: distinct (source, lipid) pairs = **1,407**
- Chance baseline (1/S): **0.0769**

- Morgan fingerprint valid: 1,407 / 1,407

**morgan_fp / logistic_regression / stratified_kfold**: balanced_acc=0.982 (+/-0.011), top1=0.995, macro_F1=0.982, chance=0.077  [0.6s]
**morgan_fp / logistic_regression / scaffold_disjoint_groupkfold**: balanced_acc=0.928 (+/-0.060), top1=0.830, macro_F1=0.858, chance=0.077  [0.5s]
**morgan_fp / random_forest / stratified_kfold**: balanced_acc=0.982 (+/-0.011), top1=0.995, macro_F1=0.982, chance=0.077  [1.8s]
**morgan_fp / random_forest / scaffold_disjoint_groupkfold**: balanced_acc=0.927 (+/-0.061), top1=0.830, macro_F1=0.872, chance=0.077  [1.7s]

**RDKit descriptor featurization (robustness check vs. Morgan FP):**

- RDKit descriptors (167 features) / logistic_regression / stratified_kfold: balanced_acc=0.977, chance=0.077

**Same-source nearest-neighbor rate**: observed = 0.996, expected under random source assignment = 0.145 (ratio = 6.88x)

**Permutation baseline** (200 shuffles, logistic regression, single 80/20 split): observed balanced_acc=0.981, null mean=0.075 (std=0.025), 0/200 null shuffles >= observed, empirical p < 0.0050 (resolution floor of 200 shuffles)

- UMAP embedding saved to `figures/stage3_embedding_LNPDB_vivo.png` (+ .svg)
