# namjh_39264000 — provenance

**Status: EXTRACTED 2026-11-10. 108 rows emitted, NOT merged. 8 OPEN decisions
block the merge.**

| | |
|---|---|
| Paper | Bae, Choi, Lee, … Nam. *Rational Design of Lipid Nanoparticles for Enhanced mRNA Vaccine Delivery via Machine Learning.* Small 2025, 21, 2405618 |
| DOI | `10.1002/smll.202405618` |
| PMID | **39264000** — indexed. Carried in the PMC record's own self-citation string; no unindexed claim is made, so no `esummary_cache.json` entry is required (CLAUDE.md §3, F-027). |
| PMCID | PMC11855220 |
| Licence | **CC BY-NC-ND 4.0.** ND is *not* CC-BY-compatible. Release blocker — DECISIONS.md D-5. |
| Experiment_ID | `SB_2025` (72 rows) and `SB_2025_2` (36 rows) — provisional, DECISIONS.md D-2 |
| Rows | **108** in the 74-column curation schema (`06_rows.csv`) |
| Value source | see `00_value_sources.json` — **every emitted value is `5_pixels`** |

## Where the values come from, and why pixels

Tiers 1–4 were each executed and each is empty. The sweep is in
`00_value_sources.json` under `tier_1_to_4_sweep`; the short version:

| tier | result | evidence |
|---|---|---|
| 1 publisher Source Data | none exists | the PMC JATS lists exactly one supplementary element, `SMLL-21-2405618-s001.docx` |
| 2 SI tables | exhausted, no delivery value in any of them | `word/document.xml` holds exactly three `<w:tbl>`: the fingerprint dictionary, the LNP composition table, and a zeta-potential table |
| 3 author repository | none exists | six-query GitHub sweep, zero hits on the DOI or PMID; the paper's own Data Availability Statement says data are "available from the corresponding author upon reasonable request" |
| 4 PDF text table | none exists | the article PDF has no tables; every figure is an embedded raster (≤9 drawing paths per page) |

Table S2 **is** used, at tier 2, for every formulation field on every row.

**The SI is a `.docx` and its figures are pictures.** Figures S1–S4 are
`word/media/image9-12.tif`, not tables. This is the hazard the brief named and
it is real here: a triage that counted `<w:tbl>` elements would have found
three tables and concluded the per-lipid data were tabulated. They are not.

**The figure rasters used are the article PDF's own embedded streams**, not
re-renders. Four of the seven are PNG, i.e. lossless, which matters for a
flat-fill heatmap. PMC's `articles/instance/11855220/bin/…` route returned a
21 kB `text/html` error page (HTTP 200) for all eight figure names on
2026-11-10, so it supplied nothing.

## The paper's own cross-checks, and whether we reproduce them

All numbers are in `07_qc.json`. **Passing all of them is necessary and not
sufficient** — this curation contained a real connectivity error that passed
both HRMS anchors unchanged (see below).

| check | kind | result |
|---|---|---|
| HRMS, SI Fig S9, two printed FAB+ peaks | external mass | **−0.32 ppm** and **−0.46 ppm** |
| Tracer calibration on five SI-named reagents | external, by name | **8 / 8 exact** |
| Connectivity: FP3 and the acyl span | internal, against a printed statistic | **36 / 36**, with 6 members at FP3 = 3 |
| Each colourbar's own drawn zero | **drawn, never printed** | **0.75 / 0.50 / 1.00 px**, i.e. 0.22% / 0.15% / 0.22% of full scale |
| Grid measured twice (gridlines vs outer borders) | geometric | worst **0.039** of a cell |
| `heatmap.match` nearest-swatch distance | geometric | worst **1.73** against a refusal threshold of 5 |
| Fig 4c's printed `r = 0.7580` | redundant encoding | **0.768** from two independently digitised panels, Δ 0.010 |
| Fig 4c's four extreme markers | redundant encoding | three agree to **≤ 0.021 in ln**; the fourth is the floor cell, see below |
| Fig 2's nine printed group sizes | printed statistic | all sum to 213, and **eleven** derived identities over Table S2's six compositions all hold |
| Fig S2b/S2d pie percentages vs Table S2 | printed statistic | **8 / 8 exact** |

### The best check here is one the authors drew and never printed

Each colourbar's **bottom edge**. The calibration is fitted on four tick
notches alone, so where it puts value 0 is a prediction about a line nobody
wrote down. It lands within **1 px** on all three panels. Anchoring the scale
3 px off the notches — the F-019 shape, a colourbar read on its own border —
moves that to 3.5–4.0 px.

### And the most informative single row is a failure that is not one

Fig 4c's leftmost marker sits at ln(I.D.) = **2.494**, i.e. 12.1 fold. Our
lowest Fig S2e cell is stored at **53.09** fold. That is a 1.48 discrepancy in
ln and it is excluded from the agreement metric on purpose: that cell matches
the colourbar's **lightest swatch**, so it carries `value_is_bound = upper`
with a band of [0, 53.09]. Fig 4c draws the same LNP at a value our colour read
cannot express, **and it lands inside the band**. The bound is both necessary
and correct, and it was set by colour identity rather than by a value threshold
(F-040).

## The error this curation contained, and what found it

The first version of `stage03_enumerate.py` built each acyl fragment by string
surgery on the block's canonical SMILES — `O=C(O)X` → `[C:1](=O)X`. In SMILES
the mapped atom bonds to the **first** atom of the remainder, and
canonicalisation had already moved the carboxyl carbon's neighbour away from
the front. The result:

* **FH-1 became linear palmitoyl** instead of 2-hexyldecanoyl;
* **every EB–T chain put the glycerol ester on the far end of the tail**,
  leaving the ester bridge as a butyryl stub.

Not one atom changed. `C58H100N2O12` either way, exact mass identical to eight
decimals, **both HRMS anchors still matching at 0.3 and 0.5 ppm.** Mass cannot
see connectivity.

What found it was reading the emitted `IL_SMILES` and noticing the head was a
straight chain. What now *catches* it is
`checks.connectivity_fp3_and_acyl_span`, which rebuilds `FH-1-EB-5T-1` by the
broken route as a **historical replay** and records both halves: the formula is
identical and the tertiary-carbon count moves 3 → 2 while the acyl span moves
5 → 4. The required tertiary-carbon count is not chosen — Figure 2b prints
FP3 = 3 for exactly 36 of the 213 LNPs and 2 for the other 177, and FH-1 is the
only head with a branch point.

Fragments are now **built from the measured counts** in `stage01_blocks.py`,
with the attachment atom atom-mapped, and never derived by rewriting a
canonical string.

## How the structures were measured

Six functional heads, three ester bridges and two tails, all from Figure 3a,
all by one tracer: binarise → connected component → skeletonise → walk →
`curation/lib` `chains.raster_chain_length`. Ring substitution patterns come
from a hexagon's own vertex angles and the exocyclic bonds that point at them.
Nothing is counted by eye.

**The tracer is calibrated before it is trusted.** Five of the blocks have an
answer that does not come from the drawing, because the SI names the reagent:

| block | tracer measured | the named reagent requires |
|---|---|---|
| T-1 | 16 carbons, hexyl branch | 2-hexyl-1-decanol (1.00 g / 4.12 mmol → MW 242.7; C16H34O = 242.44) |
| T-2 | 20 carbons, octyl branch | 2-octyl-1-dodecanol (1.00 g / 3.35 mmol → MW 298.5; C20H42O = 298.55) |
| EB-5 | 5 backbone carbons | glutaric anhydride (942 mg / 8.24 mmol → MW 114.3; C5H6O3 = 114.10) |
| FH-5 | locants [3, 5] | 3,5-dihydroxybenzoic acid |
| FH-2 | zero ring substituents | drawn as plain phenyl |

Only then is it applied to **FH-1, FH-4, FH-6, EB-6 and EB-7**, which the paper
names in no procedure, no reagent list, no NMR assignment and no mass spectrum.
Those 24 members carry `confidence = medium` with the reason on the row.

## Contradictions

`08_contradictions.csv`, **11 entries**. The ones that bear on emitted values:

* "213 LNPs" is a count of **formulations, not molecules** — Figure 2's nine
  printed group sizes partition 213 across Table S2's six compositions and all
  reconcile. The 14 × 14 block set that would give the distinct lipids is drawn
  nowhere. **The tracker's "213 lipids" is wrong.**
* Figure 2b labels FP3 "Number of Tails (Ts)"; FP3 is *tertiary carbon*, and
  every member has exactly two tails. The count is right and the label is
  wrong — and the count is what made the connectivity check possible.
* The main text reads T-1/T-2 as a carbon length ("the T carbon length ranges
  from 1 to 2") when they are identifiers for 16- and 20-carbon groups. The
  project's "name tokens are not chain lengths" hazard, in the source itself.
* Three loci disagree on the I.D. injection volume (20 µL twice, 40 µL once)
  and two on the mouse age (6 vs 7 weeks). Neither changes an emitted value;
  the dose is 5 µg in every locus.

## C10 anchor

`make_reference.py` **builds** two anchors and writes
`reference_structures.csv`. The route shares as little as possible with the
enumerator, in both evidence and mechanism:

* **Evidence.** Every integer comes from a *prose name* in the SI —
  "3,5-dihydroxybenzoic acid", "glutaric anhydride", "2-hexyl-1-decanol", and
  the full IUPAC names of compounds 4–7, which state the whole connectivity.
  The enumerator's numbers come from traced pixels. Nothing in the generator
  reads a figure.
* **Mechanism.** Plain string substitution into a template traversed **head
  first**. No `RWMol`, no `AddBond`, no `RemoveAtom`, no `CombineMols`, no call
  into `assemble`, no import of `stage01` or `stage03`. The enumerator
  traverses tail first and bonds atom-mapped fragments with `assemble.attach`.

Verified: the two anchors are **not** byte-identical to the enumeration and
are **InChIKey-identical** to it (`RCLOSQNJIDJDJR-UHFFFAOYSA-N`,
`KHHPWFLXQHXBJR-UHFFFAOYSA-N`). A third row carries an **empty**
`reference_SMILES` on purpose, so the absence of an anchor for the other 34
members is visible rather than merely missing.

## What no check in this folder can see

* **Connectivity, from mass.** Stated above, and not hypothetically.
* **FH-1, FH-2, FH-4, FH-6, EB-6, EB-7.** No external anchor of any kind.
* **Which 177 of the 213 LNPs are not these 36**, and whether the 36 are a
  biased slice. They were chosen by the authors *because* the ML ranked their
  substructures highly, so they are enriched for actives — see DECISIONS.md
  D-1, which is the entry most likely to change what gets merged.
* **Whether `SB_2024` and this paper are the same laboratory.** `SB_2024`'s 60
  corpus rows share strain, tissue, route, dose, mixing method and buffer with
  this paper's I.D. screen. If they are one lab, they are not two environments,
  and that bears directly on Results 1 and 3 — DECISIONS.md D-2.

## Reproducing

```bash
cd curation/namjh_39264000
python3 -I stage01_blocks.py     # 01_blocks.csv, 01_blocks_qc.json
python3 -I stage03_enumerate.py  # 03_enumerated.csv
python3 -I stage04_values.py     # 05_values.csv, 04_values_qc.json
python3 -I stage06_verify.py     # 07_qc.json
python3 -I stage07_emit.py       # 06_rows.csv, 06_emit_qc.json
python3 -I make_reference.py     # reference_structures.csv
```

There is no `stage02` or `stage05`: block tracing and enumeration are one step
each, and the formulation metadata is a literal transcription of Table S2
inside `stage07_emit.py` rather than a stage of its own.
