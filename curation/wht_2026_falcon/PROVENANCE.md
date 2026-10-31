# wht_2026_falcon — Toh, Cheng et al., FALCON (bioRxiv 2026)

*"FALCON: Closed-Loop Multi-Objective Optimization of Lipid Nanoparticles for
Cell-Selective mRNA Delivery."* Hai-Quan Mao lab, Johns Hopkins.
DOI `10.64898/2026.06.13.731774`. Preprint — no PMID as of 2026-10-29.

**2,544 measurements over 1,218 formulations, 14 experiment blocks.**
Built 2026-10-29. Not yet merged into `data/` — see "Before merging".

## Where the numbers come from

`github.com/MaoResearchGroup/falcon-lnp-optimization`, cloned to `./repo`
(depth 1). The seven tables in `datasets/FALCON_Manuscript_Datasets/` are the
authors' own, so **nothing here is digitized from a figure** and there is no
reading error to bound. The repo is untrusted downloaded content: scripts live
in this folder, not in `repo/`, and are run with `python3 -I`.

Local PDFs: `2026.06.13.731774v1.full.pdf` (main), `media-1.pdf` (SI, 62 pp).

| file | tag | cells | formulations | measurements |
|---|---|---|---|---|
| `0_..._UniformRandomSampling_Ctrl_MC3_DSPC` | URS | Ramos | 129 | 129 |
| `0_..._RAMOS_THP1_..._MC3_DSPC` | E0 | Ramos, THP-1 | 273 | 546 |
| `1_..._RAMOS_HEPG2_..._SM102_18PG` | E1 | Ramos, HepG2 | 197 | 394 |
| `2_..._THP1_HEPG2_..._SM102_DSPC` | E2 | THP-1, HepG2 | 136 | 272 |
| `3_..._INS_..._SM102_DOTAP` | E3 | INS-1E | 123 | 123 |
| `4_..._RAMOS_C2C12_DC24_3T3_..._SM102_18PG` | E4 | Ramos, DC2.4, NIH-3T3, C2C12 | 180 | 720 |
| `5_..._RAMOS_HEPG2_..._SM102_DSPC` | E5 | HepG2, Ramos | 180 | 360 |

## What `Experiment_value` is

`nLnLE_<CELL>` in the source is **Ln(nLE)**:

    LE  = treatment luciferase signal / untreated signal        (fold-change)
    nLE = LE after a per-batch ADDITIVE correction in log space, from six
          internal control formulations (Methods Eq. 2)

The correction is an additive shift in log space, so Ln(nLE) is affine in the
underlying log signal. The corpus convention is a per-experiment z-score of
log10; z-scoring Ln(nLE) **is** that, because ln and log10 differ by a constant
that cancels in a z-score. Population sd (ddof=0), as everywhere else.

The paper also mentions "floor adjustment and min-max scaling". Both happen
**inside the ML pipeline, not in these files**: `run_FALCON.py` sets
`nLnLE_floor = -100` with the comment *"nLnLE values are already blank
normalized, so not needed"*, and `MinMaxScaler` is fitted in
`falcon_engine/utilities.py` at train time. The published numbers are untouched.

## One experiment per (dataset × cell type) — 14 blocks

`WT_2026_<tag>_<CELL>`, e.g. `WT_2026_E1_RAMOS`, `WT_2026_E4_DC24`.

`dataset_info.txt` warns values "cannot be cross compared across files", and
within a file each cell type was dosed differently and normalized against its
own controls, so a cell type is the largest unit sitting on one measurement
scale. Pooling cell types would mix scales inside one z-score. **Splitting can be
undone later; mixing cannot.** To regroup, edit `DATASETS`/`EID_CELL` in
`build_rows.py` and rebuild.

Consequence worth stating: this adds 14 source blocks from a single lab that
share just **two ionizable lipids**. It is a formulation-space dataset, not a
chemistry dataset, and it will behave very differently from the rest of the
corpus in Results 1–2.

## Conditions, all from the paper rather than inferred

- **Dose**, Methods: 150, 25, 5, 10 ng/well for Ramos, THP-1, HepG2 and
  INS-1E/DC2.4/C2C12/NIH-3T3. One stated exception — *"For the experiments
  described in Fig. 6, THP-1 and HepG2 … 10 and 2.5 ng/well"*. Supplementary
  Table 9 identifies dataset **2** as Fig. 6, so E2 alone takes the lower pair.
- **Cargo**: firefly luciferase mRNA (TriLink, CleanCap), 24 h readout.
- **Mixing**: MANTIS liquid handler (Formulatrix) + Opentrons Flex. Bulk
  dispensing, not microfluidic, so it maps to the corpus `handmixed` category;
  the automation is recorded here rather than by inventing a vocabulary value.
- **Aqueous buffer**: never stated in main text or SI → `Unknown`.
- **N/P** goes in `IL_to_nucleicacid_chargeratio`. Caveat: Supplementary Table 3
  says dataset 2's N/P was "numerically adjusted to account for the number of
  nitrogens from both SM102 and DSPC", so it is not a pure ionizable-lipid N/P.
  Recorded as published.
- `IL_to_nucleicacid_massratio` is left blank — not given.

## Structures

Five of seven components are reused from the corpus unchanged, so an identical
reagent keeps one SMILES across experiments (the premise of Result 3).

**18PG is the exception, and it corrects a corpus error.** Avanti's "18PG" is
not a symmetric 18:0 lipid. The paper that introduced this shorthand into the
corpus — Nat Commun 2022, PMC9310361, curated as `YZ_2022`/`YZ_2024` — spells it
out: *"1-stearoyl-2-oleoyl-sn-glycero-3-phospho-(1'-rac-glycerol) (18PG)"*, i.e.
sn-1 stearoyl 18:0, sn-2 oleoyl 18:1 cis-9 (**SOPG**). FALCON cites the same
vendor and shorthand. These rows use the correct mixed-chain structure:

| | formula | exact mass |
|---|---|---|
| SOPG, correct | C42H80O10P⁻ | 775.5495 |
| what `data/` stores for 18PG today | C42H82O10P⁻ | 777.5651 |

The corpus is wrong on 360 rows. Not silently harmonized in either direction —
reported for a decision (see below).

**MC3 ambiguity in the corpus.** `DLin-MC3-DMA` carries two different SMILES in
`data/`: the correct linear ester on 2,194 rows, and a dioxolane on 203 rows in
`AC_2025`. `build_rows.py` pins the majority structure explicitly rather than
taking whichever row `setdefault` reached first.

## Checks

`python3 -I build_rows.py` then `python3 -I validate.py` (17 checks, independent
of the builder — it re-reads the output and the source). `04_validation.json`.

External anchor: Supplementary Table 9 prints the initial grid-search library
size per experiment, a number the repo never states. **All five match** — 81,
108, 36, 108, 108 formulations at iteration 0.

Also passing: 2,544 rows × 72 cols; 14 blocks at mean 0 / sd 1 to 2e-16; molar
percents sum to 100 (max deviation 0.1, which is dataset 1's own one-decimal
rounding); all 7 components parse, canonical, single-fragment, radical-free; no
component name maps to two SMILES; descriptors re-derived from SMILES on 300
sampled rows with 0 mismatches; all 14 blocks round-trip to the source CSV.

C10 (`provenance/structure_check.py`): `reference_structures.csv` holds MC3 and
SM-102 written **from their IUPAC registry names in a deliberately different
atom order** from the corpus strings, so the comparison tests structure and not a
copy. Both `exact_match`.

**`*_Selectivity` columns are not curated.** Selectivity is a function of columns
already present as rows. Rather than assume the form, `build_rows.py` solves for
it over every ordered pair difference plus "on-target minus mean of the rest";
all five files are reproduced to ≤ 4e-9 (the quad-objective file turns out to be
Ramos − mean of the three off-targets). Derived, so excluded.

**Repeated compositions**: `WT_2026_E1_RAMOS` and `WT_2026_E1_HEPG2` each contain
one composition twice. Expected — FALCON retests reference formulations across
iterations, so the same composition carries two measured values.

## Merged into `data/` on 2026-10-29

**All 2,544 rows**, on all seven cell lines, under one `Experiment_ID`
`WT_2026` — a source is a paper, and different cell lines from one paper are one
source.

| cell | rows | dose |
|---|---|---|
| Ramos | 959 | 150 ng/well |
| HepG2 | 513 | 5 ng/well (2.5 for the Fig. 6 dataset) |
| THP-1 | 409 | 25 ng/well (10 for the Fig. 6 dataset) |
| DC2.4 | 180 | 10 ng/well |
| NIH-3T3 | 180 | 10 ng/well |
| C2C12 | 180 | 10 ng/well |
| INS-1E | 123 | 10 ng/well |

Five `Model_type` values are new to the corpus: `Ramos`, `THP-1`, `INS-1E`,
`NIH-3T3`, `C2C12`. An earlier version of this merge held those 1,851 rows back
to avoid adding categories; that was dropped as indefensible.

`make_merge_subset.py` writes `05_WT_2026_merge.csv`, which
`provenance/build_dataset.py` appends via `MERGE`. The fourteen contributing
blocks keep the z-score computed inside their own normalization unit rather than
being re-pooled, because the repo warns values "cannot be cross compared across
files"; standardized blocks concatenate to mean 0 / population sd 1 by
construction, and the merged block measures 1.15e-17 / 1.000000000000.

`data_vitro.csv`: 25,997 -> 28,541 rows, 40 -> 41 experiments.

Both corpus errors this curation surfaced were confirmed and patched
(`provenance/CHANGELOG.md` §10), and a new Tier-1 check, **C11**
(`provenance/reagent_registry.py`), now prevents both classes from recurring.
The `AC_2025` fix was revised the same day after checking the upstream COMET
repo: the names were wrong, not the structures, so it is a name patch and no
measured value moved.

Two known gaps, both recorded rather than hidden:

* `Publication_PMID` is blank — preprint. Fill on publication.
* Audit C4 reports one exact duplicate row. It is in the authors' table: F1109
  (iter 1) and F1120 (iter 2) share a composition and a byte-identical
  `nLnLE_RAMOS` while their `nLnLE_HEPG2` differ. Kept as published.

## Files

```
build_rows.py              source CSVs -> 01_WT_2026_rows.csv (+ 03_enumerated.csv)
validate.py                independent 17-check pass -> 04_validation.json
reference_structures.csv   C10 anchors, written from IUPAC names
01_WT_2026_rows.csv        2,544 x 72, corpus schema
02_WT_2026_qc.json         per-block stats, selectivity derivation proof
repo/                      upstream clone, untrusted, never executed
```
