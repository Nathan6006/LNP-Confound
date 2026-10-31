# Stage 5 — Cross-Lab Shared-Lipid Concordance (Centerpiece)

Within-source normalized rank is computed within (Experiment_ID, Experiment_method) groups; raw readouts are never compared across sources directly. See module docstring for the pooling method (all cross-source-pair points concatenated) and the assay-noise-ceiling construction.

## `LNPDB_vitro`

- Distinct lipids: 11,128; shared across >=2 sources: **231**
- Source-pairs sharing >=1 lipid: **13**
- Total pooled (source-pair, lipid) points: **236**

**Per-pair concordance** (top by n; full table in results/stage5_concordance.csv):

| source_a | source_b | n_shared | spearman | kendall |
|---|---|---|---|---|
| AA_2008 | KW_2014 | 135 | -0.040 | -0.033 |
| KW_2014 | NC_2024 | 55 | 0.488 | 0.366 |
| JC_2023 | YX_2024 | 30 | 0.262 | 0.181 |
| KZ_2016 | SL_2020 | 6 | 0.788 | 0.694 |
| AC_2025 | JW_2024 | 2 | -1.000 | -1.000 |
| AC_2025 | LZ_2025 | 1 | n<2 | n<2 |
| AC_2025 | SP_2020 | 1 | n<2 | n<2 |
| AC_2025 | YZ_2022 | 1 | n<2 | n<2 |
| AC_2025 | YZ_2024 | 1 | n<2 | n<2 |
| KZ_2016 | ZC_2023 | 1 | n<2 | n<2 |
| SP_2020 | YZ_2022 | 1 | n<2 | n<2 |
| SP_2020 | YZ_2024 | 1 | n<2 | n<2 |


**Pooled concordance (primary — all cross-source points concatenated)**: Spearman = **0.083** [-0.041, 0.215] (95% bootstrap CI), Kendall = 0.065, n = 236 — **OK**
**Robustness cross-check (n-weighted average of per-pair Spearman)**: 0.140

**Chemical-bias caveat**: of 231 shared lipids, 2 (0.9%) match a known reference-lipid name (MC3/SM-102/ALC-0315/etc.); 163 (70.6%) are fully acyclic (no Bemis-Murcko scaffold to compare — reported separately, not counted as 'congeneric' since an empty scaffold is not a real shared feature); of the 68 ring-containing shared lipids, 67 (98.5% of those, or 29.0% of all shared lipids) share a real Bemis-Murcko scaffold with another shared lipid (congeneric family member). Top single pair contributes 135 / 236 (57.2%) of pooled points.

**Assay-noise ceiling** (within-source replicate split-half): Spearman = **0.400** [0.346, 0.452], n = 1035 replicate-group pairs across 12 sources (dominant source `YZ_2022` contributes 34.3% of groups — noted, not hidden). **This is the ceiling the cross-lab concordance above should be compared against.**

**Reference-lipid anchoring test**: anchor = `CCCCC/C=C\C/C=C\CCCCCCCCC(CCCCCCCC/C=C\C...` shared by 4 sources (AC_2025, SP_2020, YZ_2022, YZ_2024). Baseline pooled Spearman (same restricted source-pairs, anchor excluded) = nan (n=0) -> after anchoring = **nan** (n=0). Comparison inconclusive (insufficient n).

Figure: `figures/stage5_concordance_LNPDB_vitro.png` (+ .svg)

## `LNPDB_vivo`

- Distinct lipids: 1,434; shared across >=2 sources: **2**
- Source-pairs sharing >=1 lipid: **9**
- Total pooled (source-pair, lipid) points: **9**

**Per-pair concordance** (top by n; full table in results/stage5_concordance.csv):

| source_a | source_b | n_shared | spearman | kendall |
|---|---|---|---|---|
| BL_2023 | SB_2024 | 1 | n<2 | n<2 |
| BL_2023 | ZC_2023 | 1 | n<2 | n<2 |
| BL_2023 | ZH_2023 | 1 | n<2 | n<2 |
| SB_2024 | SS_2018 | 1 | n<2 | n<2 |
| SB_2024 | SY_2025 | 1 | n<2 | n<2 |
| SB_2024 | ZC_2023 | 1 | n<2 | n<2 |
| SB_2024 | ZH_2023 | 1 | n<2 | n<2 |
| SS_2018 | SY_2025 | 1 | n<2 | n<2 |
| ZC_2023 | ZH_2023 | 1 | n<2 | n<2 |


**Pooled concordance (primary — all cross-source points concatenated)**: Spearman = **-0.123** [-0.771, 0.691] (95% bootstrap CI), Kendall = -0.065, n = 9 — **LOW-CONFIDENCE**
**Robustness cross-check (n-weighted average of per-pair Spearman)**: nan

**Chemical-bias caveat**: of 2 shared lipids, 2 (100.0%) match a known reference-lipid name (MC3/SM-102/ALC-0315/etc.); 2 (100.0%) are fully acyclic (no Bemis-Murcko scaffold to compare — reported separately, not counted as 'congeneric' since an empty scaffold is not a real shared feature); of the 0 ring-containing shared lipids, 0 (NA of those, or 0.0% of all shared lipids) share a real Bemis-Murcko scaffold with another shared lipid (congeneric family member). Top single pair contributes 1 / 9 (11.1%) of pooled points.

**Assay-noise ceiling** (within-source replicate split-half): Spearman = **0.332** [0.196, 0.457], n = 194 replicate-group pairs across 8 sources (dominant source `ZH_2023` contributes 82.5% of groups — noted, not hidden). **This is the ceiling the cross-lab concordance above should be compared against.**

**Reference-lipid anchoring test**: anchor = `CCCCC/C=C\C/C=C\CCCCCCCCC(CCCCCCCC/C=C\C...` shared by 4 sources (BL_2023, SB_2024, ZC_2023, ZH_2023). Baseline pooled Spearman (same restricted source-pairs, anchor excluded) = nan (n=0) -> after anchoring = **nan** (n=0). Comparison inconclusive (insufficient n).

Figure: `figures/stage5_concordance_LNPDB_vivo.png` (+ .svg)

## MOLEA cross-endpoint check (independent corroborating data point)

Within MOLEA (single lab, identical chemistry, only the assayed cell type changes): Hepatocyte-vs-Chondrocyte rank Spearman = **0.299** [0.242, 0.357], raw-value Pearson = 0.381, n = 1152. This is a same-lab, single-variable-changed data point — a middle rung between within-experiment learnable agreement and cross-lab shared-lipid agreement above.
