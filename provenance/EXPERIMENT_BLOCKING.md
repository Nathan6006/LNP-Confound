# How each paper's screening rounds are blocked for z-scoring

Policy **P4(c)**, 2026-11-10. `Experiment_value` is a per-experiment z-score, so
HOW a paper's rounds are grouped changes **every value in the folder**. A
pooling choice that lives only in a curator's head is unauditable, and this
project has already paid for both directions of getting it wrong. So every
paper that screens in more than one round records its blocking here.

## The rule

1. **Multiple screening rounds inside ONE paper share one `Experiment_ID`** and
   are z-scored together as one block.
2. **But pooling asserts a common scale**, and that assertion needs evidence.
   The evidence is a **tie**: a control, benchmark or lipid measured in both
   rounds, so the offset between them is measurable rather than assumed.
3. **No tie, no silent pooling.** Where rounds share nothing, say so on the row
   (`curation_note`) and treat the pooled z as carrying an unmeasured offset.

## Why rule 2 exists — both failure modes have already happened

**Splitting when the rounds ARE comparable destroys a real effect.** `QL_2026`'s
two libraries were standardised separately at first. That deleted the paper's
central finding: a **~37x library difference** (154,630 against 4,128 RLU
geometric mean) went to zero and **36% of cross-library pairs inverted**, while
every check still passed. Four shared benchmark LNPs tied the libraries together
and the pooled z restored it.

**Pooling when they are NOT comparable invents one.** This is the mirror, and it
is why a tie is required rather than assumed. `WT_2026` is the standing example
in the other direction: its authors' own `dataset_info.txt` forbids comparing
across their seven source files, so it is normalised per **`source_file`** and
declares that on every row.

**And the tie itself can be measuring our own error.** `QL_2026`'s first
cross-library offset was **1.1116**; three of its four benchmark anchors turned
out to be contaminated by a bar-top detector reading replicate dots (F-018), and
the corrected **1.0071 +/- 0.0247 is not distinguishable from 1.0**. A tie has
to be measured on clean anchors or it is worse than no tie.

## Settled, in the corpus today

| paper | `Experiment_ID` | rounds | blocking | tie |
|---|---|---|---|---|
| Liu 2026 thiosuccinimide | `QL_2026` | Fig 2b/2c 320 lipids + Fig 3c 24-member six-tail library | **one block**, 344 rows | 4 shared benchmark LNPs; offset measured, then dropped as indistinguishable from 1.0 (F-031) |
| Toh 2026 FALCON | `WT_2026` | 7 source files | **seven blocks**, declared `source_file` | the authors' own `nLnLE_Selectivity` ties cell types WITHIN a file; across files they forbid it |
| Lin 2026 aromatic | `LL_2026` + `LL_2026_2` | 840-lipid main screen + 147-lipid initial screen **at a different formulation** | **two `Experiment_ID`s** | none -- the formulation differs, so they are two screens, not two rounds |
| Tian 2026 | `ZT_2026` | 486-formulation DoE (vitro) + 51 lipids (vivo) | **one tag, two pools** | n/a: different `Model`, so `experiment` already means experiment-x-pool |
| He 2023 | `ZH_2023` | liver + spleen | **two blocks**, declared `experiment_x_model_target` | **none, and not recoverable** -- the per-organ z is upstream (F-080) |
| Li 2023 | `BL_2023` | in-vitro + lung_epithelium + muscle | **three blocks** | none; the block is the figure arm, five columns co-vary (F-080) |
| Gong 2025 | `NG_2025_2` | liver (IV) + muscle (IM) | **two blocks**, declared `experiment_x_model_target` | none; corrected from a false `experiment_x_model_type` on 2026-11-10 |

Note `LL_2026` / `LL_2026_2`: these are **two screens in one paper**, and
F-082 found they contribute **70 shared lipids** to what Result 3 reports as
cross-laboratory concordance. One paper, one lab, two tags. P4(a) keeps the two
tags (they are two screens) and the new `lab` column is what stops them being
read as two laboratories.

## Open — staged folders that cannot merge until this is answered

Each of these has an OPEN decision asking how many blocks it is. The answer
changes every value in the folder, so it is recorded here rather than in the
curator's reasoning.

| folder | the question | what a tie would be |
|---|---|---|
| `zhangz_40433897` | Fig S4b's 60 formulations and Fig S5's 18 are different rounds. One block or six? Its `01_blocks.csv` currently emits **six**. | a formulation measured in both rounds |
| `khanof_36652629` | six rounds, or one experiment? | the benchmark LNPs on each round's panel |
| `andersondg_41203968` | "one experiment or five?" -- the in-vivo panels are not one screen | the two benchmarks on Fig 2c |
| `lij_42130331` | Table S1 (50 values, median 380) and Table S2 (250, median 1,100) -- one experiment or two? The medians differ ~3x, so this is exactly the `QL_2026` situation | nothing identified yet |
| `namjh_39264000` | Fig S2c and S2e: two formulations in one experiment, or two experiments? Currently pooled into `SB_2025`, and "the roughly two-fold difference in their maxima survives" | the shared lipids between panels |
| `siegwart_40644303` | Fig 2b and Fig 3b: one in-vivo experiment or two? Both are "Avg. Liver Radiance" from the same assay, C57BL/6, 0.25 mg/kg, 6 h | same assay and dose is itself close to a tie |
| `whiteheadka_30637934` | does Fig S6B become a separate `Experiment_ID`? It shares four lipidoids and the identical nominal dose with S6A **and disagrees with it 0.50x to 1.28x** | the four shared lipidoids -- and they say the panels do NOT agree |
| `xuq_32520527` | — | — |
| `kooh_42690843` | `Experiment_ID`, and how many experiments this paper is | — |

`whiteheadka_30637934` is the one to watch: it has a tie AND the tie says the
two panels disagree by up to 2.6-fold. That is evidence for **separate** blocks,
not pooled ones, and it is the first staged case where the tie argues against
rule 1's default.
