# Wave-1 adversarial review — WT_2026 (Toh et al., FALCON, bioRxiv 2026)

Reviewer: Wave-1 agent, round 2026-10-30. Read-only pass.
Target: `curation/wht_2026_falcon/`, merged as 2,544 rows under `Experiment_ID = WT_2026`
in `data/data_vitro.csv`.

Everything below was re-derived with my own code, run as `python3 -I` from a
scratchpad directory, never from inside `repo/`. Source PDFs read with `pymupdf`;
structures built with `rdkit` from the IUPAC names by fragment concatenation, not
copied.

**Headline:** two findings I would act on. (1) The 14-block z-score deletes the
authors' own selectivity signal *inside* single source files — up to 56x mean and
6,225x per formulation — and inverts the sign of the real cell-type difference on
80 of 180 rows in one block. This is `ERROR_CLASSES` class 22 recurring, with the
tautological "pooled mean 0 / sd 1" assertion presented as the validation, exactly
as the class says never to accept. (2) One source file (`URS`, 129 rows) carries
its composition twice in two redundant parameterizations and **they disagree on
every row**; the curation read one of them and no existing check can see it. This
is a new error class.

---

## 1. CONFIRMED CORRECT

Each item names the method and gives the number I got independently.

### 1.1 Row and formulation counts — exact

My own count over the seven CSVs in `repo/datasets/FALCON_Manuscript_Datasets/`
(one pass, counting rows and `nLnLE_*` non-Selectivity columns per file):

| file | tag | formulations | measurement cols | measurements |
|---|---|---|---|---|
| `0_..._UniformRandomSampling_Ctrl_MC3_DSPC` | URS | 129 | 1 | 129 |
| `0_..._RAMOS_THP1_..._MC3_DSPC` | E0 | 273 | 2 | 546 |
| `1_..._RAMOS_HEPG2_..._SM102_18PG` | E1 | 197 | 2 | 394 |
| `2_..._THP1_HEPG2_..._SM102_DSPC` | E2 | 136 | 2 | 272 |
| `3_..._INS_..._SM102_DOTAP` | E3 | 123 | 1 | 123 |
| `4_..._RAMOS_C2C12_DC24_3T3_..._SM102_18PG` | E4 | 180 | 4 | 720 |
| `5_..._RAMOS_HEPG2_..._SM102_DSPC` | E5 | 180 | 2 | 360 |
| **total** | | **1,218** | | **2,544** |

Matches the claim exactly. Per-cell-line row counts measured by me from
`data/data_vitro.csv`: Ramos 959, HepG2 513, THP-1 409, DC2.4 180, NIH-3T3 180,
C2C12 180, INS-1E 123 — sum 2,544. Identical to the packet and to
`06_merge_subset_qc.json`.

### 1.2 Supplementary Table 9 initial-library sizes — all five re-derived

I extracted Suppl Table 9 from `media-1.pdf` p.11 myself and counted `Iter == "0"`
in each source file:

| experiment | Table 9 prints | I counted |
|---|---|---|
| 0 (Fig. 2) | 81 | **81** |
| 1 (Fig. 3) | 108 | **108** |
| 2 (Fig. 6) | 36 | **36** |
| 3 (Suppl.) | 108 | **108** |
| 4 (Suppl.) | 108 | **108** |

All five reproduce. Table 9 contains **no row for dataset 5** (its iteration-0
count is 108, unanchored) — so "all five" is the complete set, not a convenient
subset. `validate.py`'s `TABLE9_LIBRARY` is correctly scoped.

### 1.3 A printed statistic the curation did NOT use, which also reproduces

Table 9 also prints **"Batch Size per Iteration"**. I checked it against the
observed per-iteration counts, which the curation never did:

| dataset | Table 9 batch size | my observed Iter counts | verdict |
|---|---|---|---|
| 0 | 12, four strategy rows -> 48 | 1-4 all = 48 | exact |
| 1 | 15 (NSGA) + 10 (DA) = 25 | 24, 25, 25, 15 | 2 of 4 exact |
| 2 | 20 (10 per method) | 1-5 all = 20 | exact |
| 3 | 5 | 1-3 all = 5 | exact |
| 4 | 12 + 12 = 24 | 1-3 all = 24 | exact |

Four of five reproduce exactly, including the non-obvious dataset-0 case where
48 = 4 strategies x 12. Dataset 1 is short on iterations 1 and 4 (24 and 15 against
25); I flag it rather than excuse it, but a dropped formulation and a
single-method final iteration both explain it and nothing else in dataset 1
disagrees. Recorded here so it is not re-litigated.

### 1.4 The log is taken ONCE, not twice — the convention holds

This was the single biggest risk in the brief. It is fine, and here is why.

Methods, p.25 lines 810-827, read in full:

- Eq 1: `LE = LNP Treatment Luciferase Signal / No Treatment Luciferase Signal`
  — a fold-change, **not** a log.
- "These values were then natural log-transformed to improve distribution
  symmetry." — one log.
- Eq 2: `Ln(nLE)_j = Ln(LE)_j^batch + (1/6) * SUM_{i=1..6} ( Ln(LE)_i^ref - Ln(LE)_i^batch )`
  — the batch correction is an **additive shift in log space**, not a second log
  and not a rescaling.

So `nLnLE_<CELL>` = ln(fold-change) + a per-batch constant. A z-score of that
equals a z-score of log10 of the same quantity, because ln and log10 differ by the
constant factor ln(10) which cancels in `(v - mu)/sigma`. `Experiment_value` is
therefore a per-block z-score of a log10-equivalent quantity, as the corpus
convention requires. **No double log.** `PROVENANCE.md` section "What
`Experiment_value` is" is accurate.

The one thing that *could* have broken this — the same Methods paragraph says
"Both input features and output values were then scaled to a 0-1 range using
min-max normalization" — does not apply to the published tables. I verified in the
repo source (read, never executed):

- `falcon_engine/utilities.py` lines 105-121: `MinMaxScaler()` is fitted inside
  `extract_training_data()` on a `processed_data` copy at train time, per cell type.
- `run_FALCON.py` line 59: `nLnLE_floor = -100 #minimum cutoff for output values.
  nLnLE values are already blank normalized, so not needed.`
- The published values run **-2.1429 to 10.5429**, which is not a 0-1 scale.

Confirmed: the published numbers are untouched Ln(nLE).

### 1.5 The authors' cross-comparison warning — located and read in context

`repo/datasets/FALCON_Manuscript_Datasets/dataset_info.txt`, verbatim:

> Luciferase Expression (LE) values are reported as blank-subtracted,
> batch-normalized (across iterations) values.
>
> Please note that normalization was performed within each file/experiment. Thus,
> transfection values cannot be cross compared **across files**.

The warning is scoped to **across files only**. It says nothing about cell types
within a file, and the authors themselves do cross-compare cell types within a
file (see E1). Treating a file as the normalization unit is right; treating a
(file x cell type) as the unit is not licensed by this sentence.

There *is* a real between-file difference, and the per-file separation does delete
it — per-file raw means (cell types pooled) span **3.0992 (E0) to 6.0791 (E4)**, a
2.98 ln = **19.7x** spread, with per-file sds from 1.609 to 2.797. So the procedure
costs exactly that much. Unlike `ql_2026`, here the authors' own warning says that
difference is not interpretable, so deleting it is the correct call. The
`ql_2026` lesson applies to the *within*-file half of the same procedure, which is
Finding E1.

### 1.6 End-to-end re-derivation of all 2,544 rows — 0 mismatches

I rebuilt the expected corpus slice from the seven source CSVs with my own code —
my own molar-percent column resolution, my own population-sd z-score, no import of
`build_rows.py` — keyed each source measurement as
`(f"WT_2026_{FALCON_ID or tag_FormulaLabel}", cell)`, and joined to
`data/data_vitro.csv`.

**2,544 / 2,544 rows matched a unique source measurement. 0 mismatches** on
`IL_molratio`, `HL_molratio`, `CHL_molratio`, `PEG_molratio`,
`IL_to_nucleicacid_chargeratio`, `Experiment_value` (tolerance 1e-9), `IL_name`
and `HL_name`. The key is a bijection — no source measurement is used twice and
none is missing.

Fifteen rows drawn at random (seed 7) and printed end to end, source -> corpus
(molratios in IL/HL/CHL/PEG order):

```
WT_2026_F2122  HepG2  SM-102/DSPC  89.402/ 0.456/ 9.501/0.640 sum=100.000 NP=12.230 raw=+1.616781  z=-1.375685401
WT_2026_F0216  THP-1  MC3/DSPC     38.689/13.782/43.815/3.714 sum=100.000 NP= 4.000 raw=+4.087884  z=+0.004808049
WT_2026_E4_154 Ramos  SM-102/18PG  65.258/ 3.262/28.794/2.686 sum=100.000 NP=11.854 raw=+7.905670  z=+0.629584198
WT_2026_F0069  Ramos  MC3/DSPC     49.020/ 0.980/49.900/0.100 sum=100.000 NP=12.000 raw=-0.134671  z=-1.028947193
WT_2026_F0168  Ramos  MC3/DSPC     13.768/16.569/65.994/3.668 sum=100.000 NP= 4.006 raw=+4.235919  z=+0.964074357
WT_2026_E5_11  HepG2  SM-102/DSPC  25.000/25.000/47.500/2.500 sum=100.000 NP= 7.000 raw=+7.577282  z=+0.415911646
WT_2026_F0257  Ramos  MC3/DSPC     22.477/29.499/47.766/0.259 sum=100.000 NP= 6.818 raw=+4.595040  z=+1.127835983
WT_2026_E4_34  Ramos  SM-102/18PG  56.000/24.000/19.960/0.040 sum=100.000 NP=10.000 raw=+8.589522  z=+0.963744624
WT_2026_E5_24  Ramos  SM-102/DSPC  10.000/10.000/79.200/0.800 sum=100.000 NP= 4.000 raw=+4.031214  z=+1.191820961
WT_2026_F0109  Ramos  MC3/DSPC     17.251/34.511/47.004/1.234 sum=100.000 NP= 6.837 raw=+4.832785  z=+1.236249440
WT_2026_E4_75  C2C12  SM-102/18PG  18.000/ 2.000/76.000/4.000 sum=100.000 NP= 4.000 raw=+3.495112  z=-0.953291058
WT_2026_F1008  HepG2  SM-102/18PG  44.400/33.300/22.200/0.000 sum= 99.900 NP= 7.000 raw=+1.402998  z=-1.278125785
WT_2026_F0025  Ramos  MC3/DSPC     79.602/ 0.398/18.182/1.818 sum=100.000 NP= 4.000 raw=-1.193989  z=-1.512003783
WT_2026_F0224  Ramos  MC3/DSPC     38.002/23.878/37.156/0.964 sum=100.000 NP= 5.808 raw=+4.594432  z=+1.127558592
WT_2026_E4_133 DC2.4  SM-102/18PG  47.071/20.327/32.521/0.081 sum=100.000 NP= 8.088 raw=+10.461278 z=+1.727593407
```

Every `z` reproduced my independently computed value to the last printed digit.
Sums are 100.000 except file 1's documented one-decimal rounding (99.9 / 100.1,
max deviation 0.1000).

**Nothing is transposed between components.** I re-resolved the molar-percent
columns myself and checked IL/HL assignment in all seven files, which is
non-trivial because the headers differ (`MC3_PcntMolar`, `SM102_PcntMolar`,
`SM-102_PcntMolar`, `DMGPEG_PcntMolar`, `DMG-PEG2000_PcntMolar`). All seven
resolve correctly; no IL<->HL swap; helper lipid identity matches the filename and
the `Helper_lipid` column on every row.

### 1.7 Supplementary Table 5 — a direct external anchor on the molar percentages

Suppl Table 5 (`media-1.pdf` p.8) prints `IL% / HL% / Chol% / PEG%` for 16 named
formulations. **All 16 match the stored values to the printed precision.** This is
the strongest independent evidence that the composition columns were read from the
right place — and it is the anchor that settles Finding E2.

Two caveats I note rather than suppress, because both are errors in the *paper*:

- Table 5's **U0036** row has its last three cells (`(IL+HL)%`, `HL/(IL+HL)%`,
  `PEG/(PEG+Chol)%`) printed as `6.0, 50.8, 2.0`, which is just its own first three
  cells (`N/P, IL%, HL%`) repeated — a spreadsheet drag error.
- Table 5's last three columns for **F3017 / F3118** are offset by one row relative
  to the CSVs (F3017 is printed with F3118's ratio triple). The four percentage
  columns and `N/P` are correct in both rows; the CSVs are self-consistent there.

I therefore used only Table 5's four percentage columns and its `N/P`, and say so.

### 1.8 C10 anchors are real, independent, and the right compounds

I wrote both structures myself from the registry names by **explicit fragment
concatenation in Python** (built the arm string, the acid string, the head string,
then joined), not by typing a whole SMILES, and compared InChIKeys three ways:

| | formula | exact mass | InChIKey |
|---|---|---|---|
| **DLin-MC3-DMA**, folder `reference_structures.csv` | C43H79NO2 | 641.6111 | NRLNQCOGCKAESA-KWXKLSQISA-N |
| same, `03_enumerated.csv` (the enumeration) | C43H79NO2 | 641.6111 | NRLNQCOGCKAESA-KWXKLSQISA-N |
| same, **my fragment build from the IUPAC name** | C43H79NO2 | 641.6111 | NRLNQCOGCKAESA-KWXKLSQISA-N |
| **SM-102**, folder `reference_structures.csv` | C44H87NO5 | 709.6584 | BGNVBNJYBVCBJH-UHFFFAOYSA-N |
| same, `03_enumerated.csv` | C44H87NO5 | 709.6584 | BGNVBNJYBVCBJH-UHFFFAOYSA-N |
| same, **my fragment build from the IUPAC name** | C44H87NO5 | 709.6584 | BGNVBNJYBVCBJH-UHFFFAOYSA-N |

**Were they written by a different code path?** Yes, in the sense CLAUDE.md allows.
The two strings are *not* byte-identical before canonicalization: the folder
reference puts the ester in the parentheses
(`...CCCCCCCCC(OC(=O)CCCN(C)C)CCCCCCCC...`) while the enumeration puts the second
tail there (`...CCCCCCCCC(CCCCCCCC/C=C\C/C=C\CCCCC)OC(=O)CCCN(C)C`). So the
reference was not copied out of the enumeration. I record under UNVERIFIABLE that
this is a reordering by the same hand in the same session, which is weaker than the
Lin-style fragment assembly, and that for two catalogue reagents C11 already pins
them — this particular C10 anchor adds little beyond C11.

**Are the names the right compounds?** I parsed the names myself rather than
trusting the strings.

- MC3: *(6Z,9Z,28Z,31Z)-heptatriaconta-6,9,28,31-tetraen-19-yl
  4-(dimethylamino)butanoate*. C37 chain, ester at C19 (the exact midpoint), cis
  double bonds at 6/9 and 28/31 — i.e. two identical 18-carbon linoleyl arms
  (relative cis-9,12 from each terminus), plus `-OC(=O)CH2CH2CH2N(CH3)2`. Both
  stored strings encode exactly that. It is the **acyclic ester**, not the KC2
  dioxolane — the `AC_2025` trap. Correct.
- SM-102: *heptadecan-9-yl 8-[(2-hydroxyethyl)(6-oxo-6-undecoxyhexyl)amino]octanoate*.
  C17 ester alcohol attached at C9; octanoate with N on C8; N-(2-hydroxyethyl);
  N-hexyl terminating in `C(=O)O-C11H23`. Both stored strings encode exactly that.
  Correct.

Both anchors are recorded in `provenance/structure_check.json` as
`status: exact_match` with `skeleton_match` and `stereo_match` true.

### 1.9 The exact duplicate row — verified exactly as described

`1_Norm_RAMOS_HEPG2_Single_Dual_Objective_4ITER_SM102_18PG.csv`, source rows 109
and 120:

- `FALCON_ID` **F1109** (`Iter` 1) and **F1120** (`Iter` 2)
- identical composition parameter block: `NP_ratio 7.68567276`, `IL+HL 56.9`,
  `HL_IL+HL 19.3`, `PEG_PEG+Chol 4.6`, and identical molar percents
- `nLnLE_RAMOS` **byte-identical**: both the string `2.483377939`
- `nLnLE_HEPG2` **differs**: `3.237174629` vs `4.413790043`

In `data/data_vitro.csv` this produces **exactly one** duplicate group in
`WT_2026` (ignoring `Index`/`LNP_ID`/`Formulation_ID`, which is how
`audit.py c4_replicates` keys it): the two Ramos rows, both
`Experiment_value = -0.07815798875171255`. **Zero other duplicate groups** in the
2,544 rows. The reading "carried forward between iterations" is the only plausible
one and keeping the authors' published value is the right call.

I also searched for carry-forwards that would *not* trip C4 — see E6. The good news
first: **F1109/F1120 is the only cross-iteration repeat anywhere in the seven
files.** I checked every `nLnLE_*` column in every file for repeated value strings
and cross-referenced iteration; no other repeat spans two iterations. So the
curation's carry-forward note is complete in that respect.

### 1.10 The withdrawn restriction was executed cleanly

Checked all four artefacts against each other and against the data:

- `make_merge_subset.py` line 43: `keep, held = rows, []` — no filter applied;
  `held` is never used again (dead variable, cosmetic only). Docstring states the
  withdrawal and why.
- `05_WT_2026_merge.csv`: 2,544 rows, **one** `Experiment_ID` (`WT_2026`), all
  seven cell lines at 959/409/513/123/180/180/180.
- `provenance/build_dataset.py` lines 87-95: the `MERGE` comment says "all 2,544
  author-reported measurements under one WT_2026 block, spanning seven cell lines
  ... An earlier version merged only the subset on cell lines already in the corpus;
  that restriction was withdrawn as indefensible". Accurate.
- `06_merge_subset_qc.json` agrees on every count.
- I compared the *sets* of
  `(Formulation_ID, Model_type, Experiment_value, 4 molratios, N/P, dose)` in
  `05_WT_2026_merge.csv` and in `data/data_vitro.csv`'s `WT_2026` slice:
  **identical**.

No residue of the restriction anywhere. Clean.

### 1.11 Component structures and the 18PG propagation

- All seven component names map to **exactly one SMILES corpus-wide**: DSPC 1
  (6,401 rows), DOTAP 1 (2,526), 18PG 1 (1,474), Cholesterol 1 (28,374),
  DMG-PEG2000 1 (21,736), SM-102 1 (2,295), DLin-MC3-DMA 1 (3,081). This matters
  because `build_rows.py` harvests component SMILES out of the corpus with
  `setdefault` (first row wins) — with every name unique, that lookup was safe.
  MC3 was additionally pinned explicitly, correctly, against the `AC_2025`
  dioxolane.
- The 18PG correction propagated: one 18PG structure on all 1,474 rows across
  `WT_2026`, `YZ_2022`, `YZ_2024`. C42H80O10P-, 775.5495 — matching
  `PROVENANCE.md`'s table — with CIP **(R)** at the glycerol C2, the same
  configuration as the corpus's DSPC (InChIKey skeleton `NRJAVPSFFCBXDT`, the
  correct DSPC key). Internally consistent sn-glycero-3-phospho assignment.
  `provenance/reagent_registry.json` reports `rows_mismatched: 0` over 13
  registered reagents.

### 1.12 Metadata

- **Dose.** Methods p.25 confirms "150, 25, 5, and 10 ng/well for Ramos B, THP-1
  cells, HepG2, and INS-1E/DC2.4/C2C12/NIH-3T3" and the one exception "For the
  experiments described in Fig. 6, THP-1 and HepG2 cells were treated ... at 10 and
  2.5 ng/well". Suppl Table 9 labels experiment 2 as "**2 (Fig. 6)**", so the
  E2-only exception is correct. In `data/` I measured exactly 136 THP-1 rows at
  0.01 ug and 136 HepG2 rows at 0.0025 ug, all within the E2 block, and no other
  rows at those doses. Correct.
- **N/P caveat for dataset 2** confirmed at source: Suppl Table 3 footnote,
  "N/P ratios in the table above have been numerically adjusted to account for the
  number of nitrogens from both SM102 and DSPC, without changing the formulation
  itself." Recorded as published, which is right — but see E8, it has no home in
  `data/`.
- **Selectivity is derived and correctly excluded.** Re-confirmed; the solved
  forms in `02_WT_2026_qc.json` reproduce the printed columns to <=4e-9, including
  the quad-objective `Ramos - mean(rest)`. One source oddity worth recording: file
  4's column is named `nLnRLU_Selectivity` (RLU, not LE) while every other file
  uses `nLnLE_Selectivity`. The curation's `endswith("Selectivity")` match handles
  it; no error here.
- **Vocabulary values** all sit inside the corpus's existing sets:
  `Mixing_method = handmixed` (corpus: handmixed 28,375 / microfluidics 510),
  `Experiment_method = luminescence_normalized` (23,762 rows corpus-wide),
  `Aqueous_buffer = Unknown` (2,574), `Cargo_type = FLuc` (28,773),
  `Experiment_batching = individual` (20,614). No invented categories.
- **`Publication_PMID` blank.** Correct as a preprint; I note that CLAUDE.md's own
  instruction after the `QL_2026` lesson is to look the PMID up before recording a
  paper as unindexed. A DOI of `10.64898/2026.06.13.731774` with a June 2026
  posting date is consistent with it still being unindexed, but I had no network
  access and cannot confirm, so this stays unverified rather than confirmed.

---

## 2. ERRORS FOUND

### E1 — Per-cell-type z-scoring inside one source file deletes the authors' own selectivity signal

**Files.** `curation/wht_2026_falcon/build_rows.py` lines 24-28 (docstring),
79-95 (`DATASETS`), 207-213 (the per-`(file, cell)` z-score);
`make_merge_subset.py` lines 14-25 (docstring) and 49-61 (the assertions);
`PROVENANCE.md` lines 49-58 and 157-162. Affects the **2,292 rows** in the six
multi-cell-type blocks (2,544 minus URS 129 and E3 123).

**What is wrong.** The curation chose `(file x cell type)` as the normalization
unit — 14 blocks — on the stated ground that `dataset_info.txt` warns values
"cannot be cross compared across files" and that "within a file each cell type was
dosed differently and normalized against its own controls, so a cell type is the
largest unit sitting on one measurement scale."

The second half of that is contradicted by the paper itself, in three places:

1. **Eq 3, p.26 line 876**: `Selectivity Score (on-target cell / off-target cell)
   = on-target cell nLE / off-target cell nLE`. The authors divide one cell type's
   nLE by another's **within one row of one file**. In log space that is exactly
   the difference of the two `nLnLE_<CELL>` columns.
2. **p.26 lines 871-873**: "formulations were required to exhibit on-target
   Ln(nLE) > 2 and off-target Ln(nLE) values at least one unit lower than the
   reference formulation." An **absolute** threshold applied across cell types.
3. **The curation's own proof.** `02_WT_2026_qc.json`'s `selectivity_is_derived`
   shows the printed `Selectivity` column equals `nLnLE_RAMOS - nLnLE_HEPG2` (and
   friends) to a max residual of **4e-9**. That residual *is* the demonstration
   that the two columns share one scale. The curation established the common scale
   and then destroyed it.

The dose difference between cell types is a recorded experimental covariate —
`Dose_ug_nucleicacid` carries it faithfully — not a scale incomparability.
`dataset_info.txt` is scoped to *across files*, quoted in full in 1.5 above.

**Evidence — what the procedure costs.** All measured by me from the raw source.

Between-cell-type offsets in raw Ln(nLE), within a single file:

| block | pair | mean difference | fold | sd ratio |
|---|---|---|---|---|
| E5 | HepG2 - Ramos | **+4.0342 ln** | **56.5x** | 1.933 |
| E1 | HepG2 - Ramos | +2.2000 ln | 9.0x | 1.489 |
| E0 | THP-1 - Ramos | +1.9549 ln | 7.1x | 1.070 |
| E4 | DC2.4 - NIH-3T3 | +1.7185 ln | 5.6x | 1.044 |
| E4 | Ramos - NIH-3T3 | +1.3091 ln | 3.7x | — |
| E2 | HepG2 - THP-1 | +0.6712 ln | 2.0x | 1.002 |

Variance destroyed by using 14 blocks instead of 7 (between-cell-type component of
each file's own variance):

| block | between-cell variance | file total | fraction deleted |
|---|---|---|---|
| E0 | 0.9554 | 6.1106 | 15.6% |
| E1 | 1.2100 | 6.4086 | 18.9% |
| E2 | 0.1126 | 3.0276 | 3.7% |
| E4 | 0.5732 | 4.6131 | 12.4% |
| E5 | **4.0687** | 7.8243 | **52.0%** |
| weighted, six multi-cell blocks | | | **22.7%** |

And the concrete consequence in `data/`, block E5 (one file, one plate series,
180 formulations measured in both cell types):

```
Formulation_ID   raw HepG2  raw Ramos  raw diff     fold  |  z HepG2  z Ramos  z diff
WT_2026_E5_68       10.543      1.806    +8.736   6225.5  |   +1.634   -0.575  +2.209
WT_2026_E5_97       10.265      1.714    +8.551   5170.1  |   +1.520   -0.648  +2.168
WT_2026_E5_67         9.927      1.642    +8.286   3966.3  |   +1.381   -0.706  +2.087
...
WT_2026_E5_174       0.713      1.302    -0.589      0.6  |   -2.404   -0.976  -1.428
WT_2026_E5_165       0.577      1.414    -0.837      0.4  |   -2.460   -0.887  -1.573
```

- real HepG2 - Ramos: mean **+4.034 ln (56.5x)**, range -0.837 ... +8.736
- stored z difference: mean **+6.4e-17**, range -2.537 ... +2.209
- **80 of 180 formulations have the sign of the real cell-type difference inverted
  in the corpus.**

Because the sds also differ (E5: 2.4343 vs 1.2591, a 1.93x ratio), this is not a
pure offset and cannot be undone by adding a constant back.

**The tautology.** `make_merge_subset.py` lines 20-25 argue that concatenating
standardized blocks "still satisfies the corpus convention exactly, and not by
approximation: for blocks each having mean 0 and population sd 1, the pooled mean
is sum(n_i * 0)/N = 0 and the pooled second moment is sum(n_i * 1)/N = 1 ... Asserted
below rather than assumed." Lines 59-61 then assert it. That is algebraically
necessary and therefore evidence of nothing. `ERROR_CLASSES` class 22 says in so
many words: "assert the between-block gap SURVIVES; per-block mean 0/sd 1 is a
tautology." The merge QC presents the tautology as the validation.

The same reasoning is in `PROVENANCE.md` lines 157-162 and the `build_dataset.py`
`MERGE` comment, so the doc and the code agree with each other — they are just
both wrong about what the assertion proves.

**Correct handling.** One z-score per **source file** — seven blocks, not fourteen.
That keeps the across-file separation `dataset_info.txt` licenses, and preserves
the within-file between-cell-type structure the authors themselves compute and
threshold on. Concretely: change `build_rows.py` to compute `mu, sd` over the
pooled `nLnLE` values of all cell-type columns in a file, rather than per column.

**How I would verify the fix.** Add a check that inverts the stored value back to
the authors' own `Selectivity` column: for each multi-cell file, the stored
`Experiment_value` difference between the two cell types of one formulation must
equal the printed `Selectivity` divided by a **single per-file constant** (that
file's pooled sd), for every row, to ~1e-9. Under the current 14-block scheme this
check cannot pass for any constant; under 7 blocks it passes identically. This is
the "assert the gap survives" check class 22 asks for, and it is cheap and
permanent.

---

### E2 — The URS source file carries its composition twice and the two copies disagree on all 129 rows; the corpus inherits a wrong N/P

**Files.** Source
`repo/datasets/FALCON_Manuscript_Datasets/0_Norm_RAMOS_Single_Objective_UniformRandomSampling_Ctrl_MC3_DSPC.csv`,
all 129 rows. Corpus: the 129 `WT_2026_U####` rows in `data/data_vitro.csv`
(the `WT_2026_URS_RAMOS` block), specifically their
`IL_to_nucleicacid_chargeratio`.

**What is wrong.** Every FALCON dataset file states each formulation twice:

- as four **design parameters** — `NP_ratio`, `IL+HL`, `HL_IL+HL`,
  `PEG_PEG+Chol` (these are exactly the authors' model features:
  `run_FALCON.py` line 81, `input_param_names = ['NP_ratio','IL+HL','HL_IL+HL','PEG_PEG+Chol']`);
- as four **molar percents** — `*_PcntMolar`.

The second is an exact function of the first. I derived the conversion and
validated it on the other six files:

```
HL%   = (IL+HL) * (HL_IL+HL)/100
IL%   = (IL+HL) - HL%
PEG%  = (100 - (IL+HL)) * (PEG_PEG+Chol)/100
Chol% = (100 - (IL+HL)) - PEG%
```

Reconstructed vs printed `*_PcntMolar`, max absolute component difference:

| file | max discrepancy |
|---|---|
| `0_..._RAMOS_THP1_...` | **0.0000** |
| `1_..._RAMOS_HEPG2_...` | 0.1096 (its own one-decimal rounding) |
| `2_..._THP1_HEPG2_...` | **0.0000** |
| `3_..._INS_...` | **0.0000** |
| `4_..._Quad_Objective_...` | **0.0000** |
| `5_..._RAMOS_HEPG2_...` | **0.0000** |
| **`0_..._UniformRandomSampling_Ctrl_...` (URS)** | **63.2417** |

On the URS file **all 129 rows disagree**, median 30.72 percentage points, max
63.24. It is not noise — it is a **row permutation**: 126 of 129 rows' printed
`(IL+HL, HL_IL+HL, PEG_PEG+Chol)` triple matches a **distinct other row's**
`*_PcntMolar`-derived triple to <1e-6 (the remaining 3 have no exact partner, so
the permutation is not even a clean bijection — consistent with a botched
spreadsheet sort or paste). The `PEG_PEG+Chol` column is additionally in perfect
descending order down the file, which the derived values are not.

Worked example, source row 1 (`U0001`):

```
printed parameters : NP_ratio 24.8854  IL+HL 26.6090  HL_IL+HL  3.9985  PEG_PEG+Chol 8.9516
printed percents   : MC3 38.6871  DSPC 14.0164  Chol 46.7809  DMGPEG 0.5156  (sums to 100.0000)
percents implied by the parameters: 25.545 / 1.064 / 66.821 / 6.570
parameters implied by the percents: IL+HL 52.703  HL_IL+HL 26.595  PEG_PEG+Chol 1.0901
```

**Which half is right — settled, in the direction that matters.** Suppl Table 5
publishes **U0036**: `IL% 50.8, HL% 2.0, Chol% 43.0, PEG% 4.2`. The stored
`*_PcntMolar` for U0036 is `50.8371 / 1.9716 / 43.0412 / 4.1500` — an exact match.
The reconstruction from that row's parameters gives `12.7 / 20.7 / 62.3 / 4.3` — not
a match. **So the curation stored the right molar percents.** Good news, and it
closes the worse of the two possibilities.

**But `NP_ratio` is in the misaligned block, and the corpus stores it.**
`build_rows.py` line 224 puts `NP_ratio` into `IL_to_nucleicacid_chargeratio`, so
each URS corpus row pairs the (correct) molar percents with an N/P from the
(misaligned) parameter block. Two independent printed anchors say the stored N/P
is wrong:

1. Suppl Table 5's **"Ref."** row is a self-consistent benchmark —
   `N/P 6, IL% 50, HL% 10, Chol% 38.5, PEG% 1.5, (IL+HL) 60, HL/(IL+HL) 16.7,
   PEG/(PEG+Chol) 3.8` (50+10=60 ok, 10/60=16.7 ok, 1.5/40=3.75->3.8 ok). I
   searched **all seven files** for that composition at a generous +/-0.6 pp
   tolerance. Exactly one row matches, anywhere: URS row 3, `U0003`,
   `49.744 / 10.256 / 38.462 / 1.538`. Its stored `NP_ratio` is **16.766**, not 6.
2. Suppl Table 5's **U0036** N/P is `6.0`; the CSV row says **7.199**.

For contrast, **all 15 non-URS rows in Table 5 match their CSV `N/P` to the printed
precision** (F0162 25.1/25.056, F0178 21.4/21.412, F3118 24.3/24.259, ...). The only
two N/P values in the whole table that disagree are the two that bear on the URS
file. I weight anchor 1 more heavily than anchor 2 because the "Ref." row is
internally self-consistent whereas U0036's row is visibly garbled in its last
three cells.

Magnitude in the corpus: `IL_to_nucleicacid_chargeratio` on the 129 URS rows has a
median absolute error of **6.86** and a maximum of **21.54** against the
permutation-recovered value, on a stored range of 4.01-25.54. These are not small
perturbations; several rows are off by 4-5x.

**Why no existing check could see it.** Each `*_PcntMolar` row sums to exactly
100.0000, so C5 passes. The values are finite and in range, so C6 passes. No row is
a duplicate of another, so C4 passes. Nothing in the corpus schema holds the
parameter block, so no cross-field consistency is even expressible. `validate.py`
checks the sum and the N/P sign but never reconstructs one parameterization from
the other. **This is a new error class.**

**Correct value / repair.** The molar percents stay as they are. The correct N/P
for the 129 URS rows is **not recoverable** — the permutation is not a bijection
and `NP_ratio`'s true target row is unknowable from the published material. The
correct deliverable is therefore to **blank `IL_to_nucleicacid_chargeratio` on
those 129 rows** (the corpus already carries 19 blank columns on this slice, so a
blank is representable and honest) and to record the contradiction, rather than to
keep a number that is demonstrably attached to a different formulation. Ask the
authors for the MANTIS worksheet for that run.

**The permanent check I would add (new class).** Whenever a source table states the
same quantity in two redundant parameterizations, reconstruct one from the other on
**every row** and fail on any disagreement beyond the source's stated printing
precision. On FALCON this is four lines of arithmetic, it passes 1,089 rows at
0.0000 and fails 129 at 30-63 pp, and it would have caught this at curation time.
Generalized: *a source table that says the same thing twice is a free ground truth
— always cross-reconstruct it, and never read only one of the two.*

---

### E3 — The packet's sections 5 and 6 are wrong; `brief.py` is not finding this folder's anchors or reagents

**Files.** `provenance/review/round_2026-10-30/packet_WT_2026.md` sections 5 and 6;
generated by `provenance/review/brief.py`.

The packet says, section 5: "**No anchor records found for this folder.** Either
none was created — which is itself a finding — or `structure_check.py` does not look
at this folder. Determine which." And section 6 lists C11 reagents as `{}`.

Both are false, and the packet's own rule applies ("If the packet and the
repository disagree, the packet is the bug — say so"):

- `provenance/structure_check.json` contains **two** records with
  `"paper": "wht_2026_falcon"` — `DLin-MC3-DMA` and `SM-102`, both
  `"status": "exact_match"`, both `skeleton_match` and `stereo_match` true.
- `provenance/structure_check.py` line 91 globs
  `curation/*/reference_structures.csv`, so it **does** look at this folder. The
  third option the packet did not offer — the packet is looking in the wrong place —
  is the right one.
- `provenance/reagent_registry.json` lists `WT_2026` among the experiments for
  18PG (one SMILES, 1,474 rows, `status: ok`) and three other reagents;
  `registered: 13`, `rows_mismatched: 0`. All seven WT components are either in the
  registry or deliberately excluded (DMG-PEG2000, per the registry's own note on
  polymer-bearing reagents).

**Severity: metadata only** — but it matters operationally, because a reviewer told
"no anchor exists" may spend its budget rebuilding one, or report a missing anchor
as a finding. Fix `brief.py`'s lookup (likely a folder-name vs `Experiment_ID` key
mismatch: the JSONs key on `wht_2026_falcon`, the packet on `WT_2026`).

---

### E4 — One `Experiment_ID` for seven cell lines makes source and cell line non-identifiable for five of them

**Files.** `make_merge_subset.py` (the `EID = "WT_2026"` decision),
`provenance/build_dataset.py` `MERGE` comment, `PROVENANCE.md` lines 137-155.

The stated rule — a source is a paper — is defensible on its own terms. What is
not recorded anywhere is what it does to the analysis this corpus exists to
support. I measured it:

- `WT_2026` is **2,544 of 28,885** vitro rows = **8.81%**, and is now the
  **largest single block** in the corpus (next: `YX_2024` 2,400, `AC_2025` 1,978).
- **Five of the sixteen `Model_type` values in the corpus exist only inside
  `WT_2026`**: Ramos (959), THP-1 (409), INS-1E (123), NIH-3T3 (180), C2C12 (180)
  — **1,851 rows, 72.8% of the block**. I verified by counting those `Model_type`
  values among all non-`WT_2026` rows: zero in every case.
- Therefore those five `Model_type` one-hot columns are each a strict subset of the
  `WT_2026` source indicator. **A model with both source and cell-line effects
  cannot separate them for 1,851 rows.** For a manuscript whose subject is
  between-laboratory effects, that is a structural identifiability problem, not a
  stylistic one.
- Meanwhile `WT_2026` contributes **zero unique ionizable lipids**: 2 of the
  corpus's 18,410 distinct `IL_SMILES`, both already present from other
  experiments (MC3 on 3,081 rows total, SM-102 on 2,295).

**Effect on leave-one-source-out.** Dropping `WT_2026` removes 8.81% of the rows
and **five of sixteen cell-line categories entirely**, while removing **no
chemical diversity at all**. So LOSO on this source does not test what LOSO is for:
it tests whether the model can predict cell lines it has never seen, under the
label of testing whether it can predict a lab it has never seen. Conversely, the
fold that *holds out* any other source still contains a block where one source
supplies 9% of the rows across seven biological systems and five distinct dose
levels (0.0025 / 0.005 / 0.01 / 0.025 / 0.15 ug).

**What I would do.** The corpus already has the convention for this: CLAUDE.md
section 3, "A second screen from one paper gets a `_2` suffix (`LX_2024_2`,
`LL_2026_2`)." FALCON is seven independently normalized screens, which the authors
say are not cross-comparable. Mapping them to `WT_2026` ... `WT_2026_7` is
convention-compliant, fixes E1 at the same time (each file gets one honest
z-score), and makes the per-screen dose and cell-line composition visible to any
grouped cross-validation. At minimum, the identifiability fact above belongs in
`PROVENANCE.md` and in CLAUDE.md's "Things that will bite you".

**Severity: metadata**, in that no value or structure changes — but it changes what
every downstream source-level analysis means, so I would not treat it as cosmetic.

---

### E5 — `PROVENANCE.md` contradicts itself on merge status and points to a section that does not exist

`curation/wht_2026_falcon/PROVENANCE.md` line 8:

> Built 2026-10-29. **Not yet merged into `data/`** — see "Before merging".

Line 137 of the same file:

> ## Merged into `data/` on 2026-10-29

And there is no section titled "Before merging" anywhere in the file (I grepped).
This is `ERROR_CLASSES` class 21 in mild form — provenance text asserting one thing
while the code and the rest of the document say another. One line to fix.

Separately verified and **not** a finding: line 164's "`data_vitro.csv`: 25,997 ->
28,541 rows, 40 -> 41 experiments" is a correct delta record for this merge, not a
stale absolute; 28,541 + 344 (`QL_2026`) = 28,885 and 41 + 1 = 42, the live
figures.

**Severity: cosmetic.**

---

### E6 — Twelve undocumented byte-identical measurement values that C4 cannot see

The brief asked specifically for carried-forward values that are not exact full-row
duplicates. I searched every `nLnLE_*` column in all seven files for repeated
9-significant-figure value strings and cross-referenced composition and iteration.

**Cross-iteration:** only F1109/F1120, already documented. Good.

**Within-iteration, between different formulations — 12 groups, 25 rows, none
documented anywhere:**

| block | column | repeated value | FALCON_IDs / rows | Iter |
|---|---|---|---|---|
| E0 | `nLnLE_RAMOS` | -1.123268778 | F0027, F0072 | 0, 0 |
| E0 | `nLnLE_RAMOS` | -0.585175889 | F0044, F0045 | 0, 0 |
| E0 | `nLnLE_RAMOS` | -1.562404643 | F0052, F0081 | 0, 0 |
| E0 | `nLnLE_RAMOS` | -0.653002787 | F0096, F0116 | 1, 1 |
| E0 | `nLnLE_THP1` | 0.382577318 | F0117, F0127 | 1, 1 |
| E4 | `nLnLE_3T3` | 2.492402257 | rows 25, 87, 89 | 0, 0, 0 |
| E4 | `nLnLE_3T3` | 2.743061838 | rows 141, 144 | 2, 2 |
| E4 | `nLnLE_C2C12` | 3.537975508 | rows 64, 73 | 0, 0 |
| E4 | `nLnLE_C2C12` | 1.570863151 | rows 82, 90 | 0, 0 |
| E4 | `nLnLE_C2C12` | 1.666173331 | rows 85, 89 | 0, 0 |
| E4 | `nLnLE_C2C12` | 4.324019799 | rows 138, 139 | 2, 2 |
| E4 | `nLnLE_C2C12` | 5.617457203 | rows 163, 169 | 3, 3 |

Why this is worth recording: within one iteration the Eq-2 batch shift is a single
constant, so an identical `nLnLE` means an **identical raw luminescence reading**
for two different formulations. C4 cannot see any of it (the compositions differ,
so the rows are not duplicates), and nothing in the curation folder mentions it.

I do **not** claim these are errors — see UNVERIFIABLE. They are plausible if the
luminometer reports small integers for low-signal wells, and three of the twelve
pairs are adjacent rows (F0044/F0045, rows 138/139, rows 87/89) which cuts both
ways. The deliverable is a flag, not a repair.

**Severity: changes values if real, metadata if not.** 25 of 2,544 rows.

---

### E7 — File 1's molar percents taken at one-decimal precision when the exact values are recoverable

`build_rows.py` reads the `*_PcntMolar` columns for every file. In file 1 alone
those are printed to one decimal, which is why 394 rows sum to 99.9 or 100.1. The
curation documents this correctly as the authors' rounding (`validate.py` lines
64-70) — but does not note that file 1's parameter columns hold the full-precision
values, from which the exact percents follow (I reproduce the printed one-decimal
values from them to within 0.1096, i.e. the rounding itself). Taking the percents
from the parameters there would make all 2,544 rows sum to exactly 100 and remove
the only C5 tolerance in the slice.

**Severity: cosmetic.** Max effect 0.11 percentage points on one component.

---

### E8 — Quality flags that never reached the corpus (exhaustive audit)

The brief asked for this precisely, and it had not been audited. I started from the
schema itself.

**First, the structural fact.** I grepped all 72 column names for
`flag|confid|ambig|assum|note|comment|qual|warn|uncert`. **Zero hits.** The corpus
schema has **no flag, confidence, note or ambiguity column of any kind.**
CLAUDE.md section 5 says "`positional_ambiguity`, `confidence`, and explicit
'NOT extracted' notes are the deliverable when the source is genuinely silent",
and section 6 says the posture for the manuscript is "to flag those rows in the
data (`positional_ambiguity`, `*_assumed`, `confidence`)". There is nowhere in
`data/` for any of that to live. This is a corpus-wide gap, not a WT_2026 one, but
WT_2026 is where I can enumerate the cost.

**Second, the 19 columns blank on every `WT_2026` row** (measured, with how many
non-WT rows carry each — i.e. how visible the absence is):

`IL_protonated_SMILES` (16,253 elsewhere), `IL_head_name` (16,471),
`IL_head_SMILES` (12,603), `IL_linker_name` (2,114), `IL_linker_SMILES` (2,834),
`IL_tail1_name` (12,481), `IL_tail1_SMILES` (10,018), `IL_tail2_name` (3,248),
`IL_tail2_SMILES` (3,248), `IL_tail3_name`/`_SMILES` (0), `IL_tail4_name`/`_SMILES`
(0), `IL_to_nucleicacid_massratio` (25,511), `fifthcomponent_name` (486),
`fifthcomponent_SMILES` (486), `fifthcomponent_molratio` (16,739),
`Dialysis_buffer` (4,156), `Publication_PMID` (26,341).

Most are legitimately inapplicable (these are catalogue lipids, not a
combinatorial library, so there is no head/linker/tail decomposition). Two are
real recorded gaps: `IL_to_nucleicacid_massratio` (not given by the paper) and
`Publication_PMID` (preprint). Both are documented and blank, which is the right
outcome.

**Third — what the curation knew and `data/` does not.** Exhaustively:

| fact | where it lives now | reached `data/`? |
|---|---|---|
| the 14 (or 7) normalization units | `02_WT_2026_qc.json` `blocks`, `01_WT_2026_rows.csv` `Experiment_ID` | **implicitly, undocumented** — see below |
| dataset 2's N/P is not a pure ionizable-lipid N/P (Suppl Table 3 footnote: adjusted for DSPC nitrogens) — **272 rows** | `PROVENANCE.md` prose only | **no** |
| F1109/F1120 carried-forward Ramos value — 2 rows | `04_validation.json`, `PROVENANCE.md` | **no** (visible only as a C4 hit at audit time) |
| MANTIS + Opentrons Flex robotic bulk dispensing, recorded as `Mixing_method = handmixed` | `PROVENANCE.md` prose only | **no** |
| iterative closed-loop structure: 3-5 iterations per file, batch-normalized against six internal controls; the `Iter` column was not carried; stored as `Experiment_batching = individual` | nowhere in the curation either | **no** |
| `*_Selectivity` columns are derived, with the solved form and residual per file | `02_WT_2026_qc.json` | **no** (correctly excluded, but the exclusion is invisible) |
| the raw Ln(nLE) scale, per block min/max/mean/sd | `02_WT_2026_qc.json` `blocks` | **no** — once z-scored, the 56x E5 offset of E1 is unrecoverable from `data/` at any precision |
| 18PG means SOPG, not DSPG | `PROVENANCE.md`, and now C11 | **yes**, via `reagent_registry.py` — the one that survived |
| `Aqueous_buffer` never stated | `PROVENANCE.md` | **yes**, as `Unknown` in `data/` |
| `Publication_PMID` blank because preprint | `PROVENANCE.md` | **yes**, as a blank |

**One piece of good news I verified rather than assumed.** The 14 blocks *are*
recoverable from `data/` alone: `Formulation_ID` prefix (`U` -> URS, `F<d>` ->
dataset d, `E4_`/`E5_` for the two files with no `FALCON_ID`) crossed with
`Model_type` reproduces all 14 blocks at exactly the right counts — I ran it and got
129 / 273 / 273 / 197 / 197 / 136 / 136 / 123 / 180 / 180 / 180 / 180 / 180 / 180.
But this depends on a `FALCON_ID` convention explained only in the untrusted repo's
`dataset_info.txt`, it is documented nowhere in `provenance/`, and no code in
`provenance/` knows it. If E1 is fixed by regrouping to seven blocks, this recovery
rule should be written down in `PROVENANCE.md` and asserted in `validate.py`.

**Severity: metadata.** Nothing here changes a stored value; it changes what a
reader of `data/` can know. The dataset-2 N/P caveat is the one I would prioritize,
because it is a quantitative qualification on 272 rows of a numeric column.

---

## 3. UNVERIFIABLE

Distinguished from the findings above: these are genuinely unknowable from the
material, and the correct deliverable is a flag.

1. **Which of the two parameterizations the authors actually pipetted for the 129
   URS formulations.** Suppl Table 5 settles the *percentages* for one row
   (U0036) and I rely on that. It does not settle `NP_ratio`: its two
   URS-relevant N/P values are the only two in the table that disagree with the
   CSV, and one of those two rows is visibly garbled in the paper. Settled only by
   the authors, or by the MANTIS formulation worksheet for that run. Until then the
   honest state for those 129 rows is a blank N/P.
2. **Whether the 12 within-iteration byte-identical measurement values (E6) are
   genuine integer-RLU collisions or duplicated data entry.** Settled only by the
   raw luminometer output files, which are not published.
3. **Whether `brief.py`'s failure (E3) is a key mismatch or a path problem.** I can
   see that the packet disagrees with the JSONs and that `structure_check.py` does
   glob this folder; I did not read `brief.py` itself, as it is outside this
   target.
4. **`Publication_PMID`.** Correct to leave blank on the evidence I have, but I had
   no network access and could not run the lookup CLAUDE.md now requires after the
   `QL_2026` lesson. Someone should run it rather than inherit the assumption.
5. **The aqueous buffer.** Never stated in the main text or the 62-page SI; I
   searched both. `Unknown` is right.
6. **DMG-PEG2000's PEG block length.** The corpus stores a 45-mer
   (C122H242O50, MW 2509.2); the real reagent is a distribution. Deliberately
   outside C11 and correctly so. Not a WT_2026 issue.
7. **DOTAP's C2 stereocentre.** The registry SMILES asserts none and the paper
   states none; commercial DOTAP is supplied as a salt of unstated configuration.
   Not asserting it is the right call.
8. **Whether the C10 anchor's "deliberately different atom order" is a genuinely
   independent code path.** It satisfies CLAUDE.md's letter and is provably not a
   byte copy (1.8 above), but it is a reordering by the same hand in the same
   session rather than a fragment assembly, and the error-class registry says this
   is unfalsifiable after canonicalization by design. For two catalogue reagents
   C11 already pins the structures by InChIKey, so little rides on it here. Not a
   finding; recorded so it is not mistaken for a strong anchor.

---

## 4. SEVERITY

| # | finding | severity | rows affected |
|---|---|---|---|
| **E1** | per-cell-type z-score deletes the within-file selectivity signal | **changes values** — 22.7% of within-file variance removed, 52% in E5, sign of the real cell-type difference inverted on 80/180 E5 rows | **2,292** (the six multi-cell blocks) |
| **E2** | URS source table self-contradictory; stored N/P from the misaligned block | **changes values** (`IL_to_nucleicacid_chargeratio`, median error 6.86, max 21.54). Molar percents verified correct via Suppl Table 5, so structures and compositions are **not** affected | **129** |
| **E3** | packet sections 5/6 contradict `structure_check.json` and `reagent_registry.json` | **metadata only** (but misdirects reviewers) | 0 data rows |
| **E4** | one `Experiment_ID` makes source and cell line non-identifiable for 5 cell lines | **metadata only** in the data; changes the meaning of every source-level analysis | 1,851 of 2,544 |
| **E5** | `PROVENANCE.md` says "not yet merged" and cites a missing section | **cosmetic** | 0 |
| **E6** | 12 undocumented byte-identical within-iteration measurement values | **changes values if real; metadata if a genuine collision** | 25 |
| **E7** | file 1's percents taken at 1-dp when exact values are recoverable | **cosmetic** (<=0.11 pp) | 394 |
| **E8** | no flag/confidence column exists; 6 recorded facts never reached `data/` | **metadata only** | up to 2,544 (272 for the N/P caveat) |

**No structural errors.** Both ionizable lipids, all three helper lipids,
cholesterol and the PEG-lipid are correct; the C10 anchors reproduce
independently; the `AC_2025` MC3/KC2 trap was avoided explicitly; the 18PG
correction is right and propagated to all 1,474 rows.

---

## 5. CONFIDENCE

| # | confidence | basis, and what would change my mind |
|---|---|---|
| **E1** | **certain** that the within-file between-cell-type structure is deleted and that the authors cross-compare cell types within a file (Eq 3, the absolute Ln(nLE)>2 threshold, and the curation's own 4e-9 selectivity residual). **Certain** on every number in the tables, all measured from the raw source. The *judgement* that 7 blocks is the right unit rather than 14 is **probable** — it would be wrong if the authors confirmed the six internal controls were separately scaled per cell type in a way that destroys the cross-cell scale, but Eq 2 as written is a single additive shift and Eq 3 presumes the scale survives. |
| **E2** | the self-contradiction in the URS file is **certain** (126/129 exact permutation match; reconstruction exact to 0.0000 on five other files). That the **molar percents are the correct half** is **probable, strongly** (Suppl Table 5's U0036 matches them to 1 dp and does not match the alternative). That **`NP_ratio` is therefore wrong on 129 rows** is **probable**: two printed anchors disagree with it and 15/15 non-URS Table-5 N/P values agree with their CSVs, but I cannot rule out that the sort covered only the three columns to its right and left `NP_ratio` intact, in which case the URS N/P values are fine and Table 5's two URS N/P cells are both typesetting errors. The residual risk runs the other way too: if Table 5's U0036 percentages were themselves copied out of the same CSV, the percents could be the misaligned half, and then 129 rows' formulations are wrong by a median 30.7 pp. That is the worst case and it is why I would not leave this unflagged. |
| **E3** | **certain.** I read both JSONs and the glob in `structure_check.py`. |
| **E4** | **certain** on every count (5 `Model_type` values with zero non-WT rows; 1,851 rows; 8.81%; 2 of 18,410 `IL_SMILES`, 0 unique). The recommendation to split into seven IDs is a judgement, **probable**. |
| **E5** | **certain.** Two lines of the same file, and the cited section does not exist. |
| **E6** | **certain** that the 12 groups exist and that nothing records them. **Speculative** that any of them is a data-entry error rather than a genuine identical reading — hence listed under UNVERIFIABLE as well, and the deliverable is a flag. |
| **E7** | **certain**, and trivial. |
| **E8** | **certain** on the schema fact (zero flag columns in 72) and on each row of the table, all measured. The prioritization is a judgement. |

---

### Note on method

Per the attack brief, nothing here is a rerun of the curation's scripts. I built my
own molar-percent column resolver, my own population-sd z-score, my own
`(Formulation_ID, Model_type)` join, my own `(IL+HL, HL_IL+HL, PEG_PEG+Chol)` ->
percent conversion, and my own SMILES for both C10 anchors by fragment
concatenation from the IUPAC names. The source PDFs were read with `pymupdf` for
Methods Eq 1-2 and Suppl Tables 1-9. Nothing under `repo/` was executed; everything
was run as `python3 -I` from a scratchpad with paths passed as arguments.

The one place I agreed with the curation by reproducing its logic is 1.6, and there
the point is the agreement itself: 2,544 rows reproduced to 1e-9 by an independent
route means the extraction pipeline is sound. The two real findings are both about
**what the source tables mean**, which is exactly where the brief said the risk had
moved.
