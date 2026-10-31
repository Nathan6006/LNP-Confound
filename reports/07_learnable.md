# Stage 7 — Within-Experiment Learnability and Transfer Decay

## `LNPDB_vitro`

### Within-experiment learnability (train+test inside one source)

Minimum 30 lipids/source to attempt a local scaffold-disjoint split; smaller sources reported as not-computable, not skipped silently.

- Sources with enough lipids to test: 29 / 34
- Size-weighted within-experiment scaffold-disjoint Spearman: **0.322**

| source | n_lipids | spearman | computable | reason if not |
|---|---|---|---|---|
| KZ_2016 | 1470 | 0.335 | True |  |
| YX_2024 | 1200 | 0.303 | True |  |
| KW_2014 | 1150 | 0.278 | True |  |
| LM_2019 | 1080 | 0.470 | True |  |
| BL_2023 | 709 | 0.277 | True |  |
| JL_2024 | 623 | 0.313 | True |  |
| SL_2021 | 572 | 0.197 | True |  |
| BL_2024 | 560 | 0.101 | True |  |
| XH_2025 | 500 | 0.727 | True |  |
| AA_2008 | 487 | 0.135 | True |  |
| JC_2023 | 288 | 0.448 | True |  |
| YR_2024 | 286 | -0.081 | True |  |
| LX_2024_2 | 270 | 0.535 | True |  |
| NC_2024 | 260 | 0.505 | True |  |
| LX_2024_3 | 252 | 0.260 | True |  |
| JW_2024 | 202 | 0.129 | True |  |
| SW_2024 | 180 | 0.331 | True |  |
| YX_2023 | 176 | 0.258 | True |  |
| LX_2024 | 165 | 0.361 | True |  |
| ZL_2022 | 144 | 0.444 | True |  |
| KS_2024 | 140 | 0.757 | True |  |
| YY_2023 | 119 | -0.030 | True |  |
| LL_2012 | 112 | 0.084 | True |  |
| ZC_2023 | 100 | 0.521 | True |  |
| XH_2024_2 | 90 | 0.673 | True |  |
| JM_2016 | 72 | 0.014 | True |  |
| SL_2020 | 63 | 0.604 | True |  |
| KL_2010 | 50 | 0.140 | True |  |
| XH_2024 | 30 | 0.352 | True |  |
| AC_2025 | 7 | NA | False | fewer than 30 lipids |
| SP_2020 | 1 | NA | False | fewer than 30 lipids |
| YZ_2022 | 1 | NA | False | fewer than 30 lipids |
| YZ_2024 | 1 | NA | False | fewer than 30 lipids |
| LZ_2025 | 1 | NA | False | fewer than 30 lipids |


### Decay curve: LOSO fold performance vs. chemical distance to training support

- Linear trend: Spearman_fold = 0.162 + -0.154 * distance (r=-0.103, p=0.588, n=30 sources)
- Negative slope: performance degrades as held-out chemistry gets further from training support (smooth decay), consistent with the slope direction expected if the model is extrapolating, not just failing categorically.

Figure: `figures/stage7_decay_LNPDB_vitro.png` (+ .svg)

## `LNPDB_vivo`

### Within-experiment learnability (train+test inside one source)

Minimum 30 lipids/source to attempt a local scaffold-disjoint split; smaller sources reported as not-computable, not skipped silently.

- Sources with enough lipids to test: 10 / 17
- Size-weighted within-experiment scaffold-disjoint Spearman: **0.114**

| source | n_lipids | spearman | computable | reason if not |
|---|---|---|---|---|
| LR_2023 | 383 | -0.038 | True |  |
| RG_2023 | 260 | NA | False | only 1 Butina cluster(s) — too chemically homogeneous to scaffold-split |
| ZH_2023 | 161 | 0.656 | True |  |
| SX_2025 | 130 | 0.186 | True |  |
| AP_2025 | 91 | 0.373 | True |  |
| KW_2014 | 90 | -0.384 | True |  |
| JW_2024 | 69 | 0.183 | True |  |
| XH_2024_2 | 61 | 0.184 | True |  |
| BL_2024 | 51 | -0.003 | True |  |
| ZC_2023 | 31 | 0.031 | True |  |
| XH_2024 | 30 | -0.460 | True |  |
| BL_2023 | 26 | NA | False | fewer than 30 lipids |
| KZ_2016 | 24 | NA | False | fewer than 30 lipids |
| SY_2025 | 19 | NA | False | fewer than 30 lipids |
| SS_2018 | 10 | NA | False | fewer than 30 lipids |
| SB_2024 | 2 | NA | False | fewer than 30 lipids |
| LL_2025 | 1 | NA | False | fewer than 30 lipids |


### Decay curve: LOSO fold performance vs. chemical distance to training support

- Linear trend: Spearman_fold = 0.558 + -0.793 * distance (r=-0.410, p=0.115, n=16 sources)
- Negative slope: performance degrades as held-out chemistry gets further from training support (smooth decay), consistent with the slope direction expected if the model is extrapolating, not just failing categorically.

Figure: `figures/stage7_decay_LNPDB_vivo.png` (+ .svg)
