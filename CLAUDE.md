# CLAUDE.md — Confound_analysis

Living project doc. Update it whenever a number is frozen, a decision is made, or a task changes status. Do not let facts drift between here, the code, and the manuscript — this file is the source of truth for **what is settled**; the frozen analysis script (see Task T3) is the source of truth for **every reported number**.

---

## 1. What this project is

**Paper (target: methods/ML audience):** *"Cross-laboratory prediction of lipid-nanoparticle transfection is structurally confound-limited."*
Authors: Nathan Liu, Qikui Zhu, Zheng-Rong Lu (CWRU BME). Corresponding: zxl125@case.edu.

**Central claim.** The standard LNP-ML workflow pools delivery measurements from many published screens and evaluates with random or cluster CV. We show that on the aggregated data now available, cross-laboratory prediction is **not learnable in a specific, measurable sense** — reported CV scores largely measure recognition of dataset *provenance*, not transferable structure–activity learning. The task is **extrapolation across near-disjoint chemical support**, not robustness to a shared distribution, so domain-generalization methods (IRM, group-DRO) that assume overlapping environments do not apply.

**The six manuscript results:**
1. Ionizable-lipid chemistry predicts its source of origin at near-ceiling accuracy.
2. Ionizable lipids are ~source-exclusive; excipient scaffold is shared but nearly invariant.
3. **Identical lipids rank near-independently across labs** (the load-bearing result — only place chemistry is held constant).
4. Honest evaluation (leave-one-source-out) collapses relative to random/cluster CV.
5. Positioning: this is extrapolation across disjoint support → DG methods can't help.
6. A within-experiment regime stays learnable (basis for a separate modeling+wet-lab paper, in prep).

A companion **modeling + wet-lab validation manuscript** (in preparation) is cross-referenced for the within-experiment regime and the model architecture used in Result 4.

---

## 2. Repository layout

```
confound/
├── CLAUDE.md                 ← this file
├── data/                     ← CANONICAL corpus. Only three files live here.
│   ├── data_vitro.csv            ★ 28,885 rows x 79 cols, 42 experiments, in vitro
│   ├── data_vivo.csv             ★  1,780 rows x 79 cols, 19 experiments, in vivo
│   └── Original_LNPDB.csv        untouched upstream reference (19,797 x 66)
├── provenance/               ← how data/ was built, and proof that it is right
│   ├── README.md                 lineage, precedence rules, every validation
│   ├── CHANGELOG.md              every change from original LNPDB to now
│   ├── build_dataset.py          rebuilds data/ from scratch, reproducibly
│   ├── audit.py                  INDEPENDENT audit (C1-C11), not a rerun of the builder
│   ├── dataset_validation.json / audit_findings.json
│   └── superseded/               pre-2026-10-28 corpus files + the official LNPDB download
├── curation/                 ← one folder per paper curated BY US
│   ├── PROTOCOL.md               stages 00-07, checks C1-C9, blind-validation rule
│   ├── siegwart_42634137/        Tian 2026 -> ZT_2026 (486 vitro + 51 vivo)
│   ├── lin_aromatic_2026/        Lin 2026  -> LL_2026 (840) + LL_2026_2 (147)
│   ├── wht_2026_falcon/          Toh 2026 FALCON -> WT_2026 (2,544 rows, 7 cell lines)
│   ├── xd_2026_lipidai/          Duan 2026 LipidAI -> XD_2026 (12 vivo) + standalone LipidDB
│   ├── ql_2026/                  Liu 2026 thiosuccinimide -> QL_2026 (344 rows, HeLa)
│   ├── liu_iphos/                blind validation vs SL_2021 (recall 97.7%, rho 1.000)
│   ├── amino_lipid_fg_2026/      blind validation vs FG_2026 (30/30 building blocks)
│   └── new_curation.csv          Siegwart rows in corpus schema, merged by build_dataset.py
├── lit_review/               ← candidate-source sweep + tracker UI
│   ├── index.html                sortable viewer; status = new | added | have | cite | dup
│   ├── build_data.py             ADDED registry records what we curated and merged
│   └── provenance/               the 20 verbatim PubMed queries + PMID lists
├── scripts/                  ← frozen analysis pipeline (Stages 0-9)
└── run_all.py                single entry point, seed 42
```

**⚠ The frozen pipeline currently points at moved files.** `scripts/stage00_inspect.py`
reads `data/LNPDB_vitro.csv` / `data/LNPDB_vivo.csv`, which moved to
`provenance/superseded/` on 2026-10-28. It will raise FileNotFoundError until
repointed. **Do not silently repoint it at `data_vitro.csv`** — the frozen
results in `results/` and `reports/` were computed on the 19,028-row corpus, and
the current corpus is 28,885 rows, so re-running changes every reported number.
That is a research decision (see §6), not a path fix.

Sibling dirs under `LNP/`: `LNP_ML/` (the companion modeling project; `new_data/`
holds the post-change-log processed corpus that `build_dataset.py` reads, and
`data_files/<lab>/` holds per-paper source tables that have repeatedly resolved
ambiguities), `LNPCD_web/`, `Papers/`.

## 3. Data dictionary (corpus rebuilt 2026-10-29 — see provenance/README.md)

**"Source" = `Experiment_ID`.** Human-readable `INITIALS_YEAR` (`AA_2008`, `KW_2014`, `ZT_2026`). Use it as the source/environment key everywhere. A second screen from one paper gets a `_2` suffix (`LX_2024_2`, `LL_2026_2`).

| | in vitro | in vivo |
|---|---|---|
| File | `data/data_vitro.csv` | `data/data_vivo.csv` |
| Rows | **28,885** | **1,780** |
| Columns | 79 | 79 (identical header) |
| Experiments | **42** | **19** |

79 columns = the 67 LNPDB curation columns + 5 RDKit descriptors from the ML
change log (`Num_carbon_in_tail`, `num_protonatable_nitrogens`,
`num_permanent_cationic_N`, `formal_net_charge`, `num_unsaturated_cc_bonds`)
+ **7 quality/provenance columns added 2026-10-30**, defined once in
`provenance/schema_quality.json` and read from there by the builder, by every
curation that emits a flag, and by `audit.py` C13. Do not restate that list in
code: these columns exist *because* per-curation flags were restated in six
places, drifted, and then died before the merge (review finding F-001).

| column | what it means |
|---|---|
| `Experiment_value_within_model` | `Experiment_value` re-standardised inside (experiment × `Model_type`). Exact, not approximate — a per-block z of an affine function of the log signal is the per-block z of the signal. Equals `Experiment_value` where an experiment has one `Model_type`. |
| `value_normalization_unit` | `experiment` (the §3 convention), `source_file` (`WT_2026`), or `experiment_x_model_type` (error class 22 — the between-assay offset has been deleted; only legitimate where it is unrecoverable, and then the reason is in `curation_note`). |
| `confidence` | `high`/`medium`/`low` — the **structure's** confidence, as assessed by the curation that built it. Must be that lipid's own assessment; a hard-coded literal once turned 24 disclosed `low` rows into `high` (F-011). |
| `positional_ambiguity` | 1 where which heteroatom bears which substituent is a convention, not a measurement. 987 rows. |
| `structure_assumed` | 1 where a specific feature was assumed: a quaternization site the paper draws nothing for, a commercial reagent sold as an isomer mixture, an ester connectivity inferred from one NMR shift. Pairs with `curation_note`, which must say *which* feature. |
| `value_is_bound` | `upper` (a heatmap cell on the colormap floor — the truth is at or below it) or `lower` (censored past the top of an axis). A model treating these as measurements is fitting the floor of a colour ramp. |
| `curation_note` | free text, one sentence, written for someone reading `data_vitro.csv` who has not opened the curation folder. |

**A blank is "not stated", not "clean."** `value_normalization_unit` is populated
on every z-scored row; the rest are populated only where a curation actually
assessed them, so a blank `confidence` on a curation that never assessed
confidence is an absence of information, not an assurance.

One-hot columns are deliberately **not** stored — the category sets change every
time a paper is added, so any baked-in encoding goes stale. Generate at model
time with `pd.get_dummies`.

### Things that will bite you

- **`ExactMolWt` was called `van.der.Waals.Molecular.Volume`.** It never held a
  van der Waals volume. Verified on 40,472/40,542 rows across all three upstream
  files; the 70 exceptions are stale rows in `JW_2024`/`LX_2024` whose
  `Molecular.Weight` also disagrees with their own SMILES. Renamed 2026-10-28.
- **`Experiment_value` is a per-experiment z-score of log10**, population sd
  (ddof=0). For a **discretized** source the quantity z-scored is the **bin
  ordinal**, not a midpoint — verified against `SL_2021`, whose four stored
  values are exactly equally spaced. Eleven experiments were off this convention
  and were rescaled; all 57 experiment blocks now sit at mean 0 / sd 1 to 1e-9.
- **`SB_2024` and `SY_2025` (vivo, 60 rows) are NOT z-scores.** They are
  `luminescence_relative_to_Spikevax` — benchmark ratios floored at 0, mean
  ≈0.44. Flagged, not fixed, not dropped: they are the vivo analogue of the
  `ZL_2022`/`YX_2023` problem that got those dropped from vitro. **Open decision.**
- **48 duplicate SMILES remain** (39 with conflicting values) across `KW_2014`,
  `BL_2023`, `AA_2008`, `JW_2024`. These are head-group collisions repeating
  across whole tail series. `LNP_ML/data_files/Whitehead_siRNA/`'s own note shows
  another curator hit this and deduplicated. Left flagged.
- **212 rows have an unassigned stereogenic C=C** (`FG_2026` 120,
  `YX_2026_SDL` 79, `NC_2024` 13). Matters only if you use stereo-aware
  fingerprints.
- **Several `BL_2023` `IL_name` values carry a trailing space.**
- **A reagent NAME is not a structure — check it, and check WHICH is wrong.**
  Two 2026-10-29 bugs: `18PG` is 18:0/18:1 (SOPG), not symmetric 18:0, and
  `AC_2025` had `DLin-MC3-DMA` and `KC2` carrying each other's structures — same
  formula, same exact mass (641.6111 Da), invisible to every composition check.
  For `AC_2025` the **names** were wrong, not the structures: the COMET source
  (`github.com/alvinchangw/COMET`) stores SMILES and labels with no names at
  all, so the values belong to the structures. Fixing it the other way would
  have moved every measured value onto the wrong molecule. Both fixed; `C11`
  (`provenance/reagent_registry.py`) now pins 13 named reagents to structures
  written from their chemical names and checks every row. Add an entry whenever
  a new catalogue reagent enters the corpus. PEG-lipids are deliberately absent
  — the block length is a distribution, so no single SMILES is right.
- **`WT_2026` is the only experiment with a blank `Publication_PMID`** — it is a
  preprint. Left blank rather than invented; fill it on publication. `QL_2026`
  was blank for the same stated reason and should not have been: PMID 41145968
  was already indexed and its DOI matches. **Look the PMID up before recording
  a paper as unindexed.**
- **288 IL names carry different structures in different experiments** (e.g.
  `A1_B1_C1` in both `JC_2023` and `YX_2024`). Generic combinatorial labels
  reused by unrelated papers, not an error, but never key on `IL_name` alone
  across experiments — use the SMILES.

### Rebuild / verify

```bash
python3 -I provenance/build_dataset.py   # rebuilds data/ from the archived inputs
python3 -I provenance/audit.py           # independent C1-C11 audit, writes audit_findings.json
```

## 4. Task tracker (advisor's tiered list)

Status keys: ☐ todo · ◐ in progress · ☑ done · ⊘ dropped/parked.

### Tier 1 — blocks everything else
- ☑ **T1. Fix shared-lipid concordance (Result 3).** Done in `scripts/stage05_concordance.py` (Stage 5, "centerpiece"). All 231 shared lipids across 13 in-vitro source-pairs; within-source-and-endpoint normalized ranks; **pooled** concordance (primary: all cross-source-pair points concatenated, Spearman **0.083** [CI −0.041, 0.215], n=236; robustness cross-check: n-weighted per-pair average = 0.140). Per-pair table in `results/stage5_concordance.csv` so no single pair dominates silently (top pair AA_2008↔KW_2014 is 57% of pooled points — flagged, not hidden). In-vivo version computed too but LOW-CONFIDENCE (n=9, CI [−0.771, 0.691]) — supports T8.
- ☑ **T2. Assay-noise control for T1.** Done in the same Stage 5 script: within-source replicate-row split-half ceiling = **0.400** [0.346, 0.452], n=1035 replicate-group pairs (vitro) vs. the cross-source 0.083 above — the airtight contrast the story needed. Also ran a reference-lipid anchoring test (MC3-anchored re-normalization); inconclusive for both pools since the sources sharing the anchor lipid share no *other* lipid to test anchoring against — reported honestly as inconclusive, not forced.
- ☑ **T3. Freeze ONE analysis script; regenerate every number.** Done. `scripts/stage00_inspect.py` → `stage09_summary.py`, single entry point `run_all.py`, global seed 42, pinned versions in `environment/requirements.txt`, `results/registry.json` is the single machine-readable source of truth, `reports/SUMMARY.md` is the human-readable rollup with a claims-vs-evidence table. Full end-to-end run verified byte-identical across two independent runs (~41 min on this M3/16GB machine). `environment/manifest.json` records versions/seed/timestamp/runtime per run (git commit hash is null until the repo's first commit).
- ☑ **T4. Decay analysis.** Done in `scripts/stage07_learnable.py` (Stage 7): LOSO fold Spearman vs. median NN-Tanimoto distance from the held-out source to the training support. Weak, non-significant negative trend for both pools (vitro: slope −0.154, r=−0.103, p=0.588, n=30; vivo: slope −0.793, r=−0.410, p=0.115, n=16) — a smooth-decay direction but not a strong or significant one; reported as such, not oversold.

### Tier 2 — strengthens results
- ☑ **T5. Protocol comparison under TWO metrics.** Done in `scripts/stage06_protocol.py` (Stage 6): 2×3 grid (pooled-regression / within-source-pairwise-ranking objective × random / Butina-scaffold-disjoint / LOSO split), XGBoost, both pooled Spearman and within-experiment Spearman reported per cell. Vitro pooled-regression: 0.521 (random) → 0.289 (scaffold-disjoint) → 0.093 (LOSO) — the inflation claim, cleanly monotone.
- ☑ **T6. Guard "falls toward chance."** Done in the same Stage 6 script: per-fold n and a per-fold permutation-null baseline (fast Ridge regressor on shuffled labels) for every LOSO fold. Vitro pooled-regression: **22/34** LOSO folds land at-or-below their own null; within-source-pairwise objective is more robust (**17/34**) — the gauge-invariant objective partially, not fully, rescues LOSO performance.
- ☑ **T7. Pre-empt "of course" on source classification.** Done in `scripts/stage03_separability.py` (Stage 3): scaffold-disjoint (Butina-cluster GroupKFold) classification alongside standard stratified k-fold. Vitro: 0.967 (stratified) → 0.880 (scaffold-disjoint); vivo: 0.982 → 0.927. Separability survives scaffold-disjointness — stated plainly as expected/the point, and validated as real (not a classifier artifact) by the Stage 8 split-half negative control, which *does* collapse to chance (0.496/0.490 vs. 0.500).
- ☑ **T8. In-vivo data fate → SUPPLEMENTARY.** Settled by data: only 2 shared lipids in vivo; "essentially perfect separability" rests on few comparable environments. Keep one sentence in main text, move detail to supp. Reconfirmed by the frozen script — Stage 5's in-vivo concordance is explicitly LOW-CONFIDENCE (n=9).

### Tier 3 — positioning
- ☐ **T11. Published-baseline reproduction (optional).** Show a real released-code architecture's reported number drop under LOSO. Most persuasive for methods audience, most work + politically delicate. Do it only with a clean released-code baseline, framed as protocol comparison (not criticism); else skip.

---

## 5. Curating a new paper into the corpus

One folder per paper under `curation/`, named `<lastname>_<pmid>` or
`<lastname>_<topic>_<year>`. `curation/PROTOCOL.md` has the full stage list.
What actually matters, learned the hard way:

**Figures are raster.** Most of these PDFs are flattened — no text layer inside
figures. Render with `pymupdf` at 600 dpi and *look at the image* with the Read
tool. Resolution matters more than cleverness: glyph clustering that was hopeless
at 300 dpi became exact at 600.

**Prefer symbolic over pixel.** Printed numerals > discrete-binned heatmap >
continuous heatmap > bar chart > scatter. A discrete legend (`<1e4 / 1e4-1e5 /
...`) means nearest-neighbour matching against swatches sampled from the same
raster — exact, not estimated. Report max colour distance and the count of cells
within 10 of a second swatch; **anything above ~5 means the geometry is wrong,
not the colours.** That metric has caught an off-by-one header row twice.

**OCR is not good enough.** Tesseract misread `1.2x10^5` as `109`. Segment
glyphs into connected components, binarize, preserve aspect ratio, cluster with
complete linkage, and label ~14 cluster representatives by eye.

**Name tokens are not chain lengths.** `KS_2024`'s `C24` means two C4 tails;
`KZ_2016`'s `SC1/SC2/SC3` are tail identifiers (isobutyl/isopentyl/neopentyl);
`YX_2024`'s `C1/C2/C3` are component indices. Never assume.

**Find the paper's own cross-checks and use them.** Marginal hit-rate bar charts,
stated library sizes, HRMS in the SI, NMR integrals, a named lead compound,
percentages printed next to ratios. Measure them by script so no human sits on
both sides of the comparison. Siegwart gave 7 such checks; Lin gave 48 hit rates
plus 4 HRMS anchors — though see the caveat below about how Lin's were
converted, which weakened all 48 into a shape comparison.

**But the scale that converts pixels to units must come from axis furniture
ONLY** — tick marks, labelled gridlines, printed anchors — and never from the
quantity being checked. `digitize_fig2fgh.py` fitted its pixels-per-percent by
least squares onto the very hit rates the panel was validating, so a uniform
error in the heatmap read was absorbed and invisible; the measured tick spacing
was computed and never used, and the file had no `assert` (error class 27). The
signature to grep for is a fit whose x-vector is the expected answer.

**And passing every printed statistic is necessary, not sufficient.** In
`QL_2026` every printed Fig 3d/3e hit rate and the "17 lipids" sentence still
reproduce on bar values that are up to 28% too high, because the error is
geometric and roughly multiplicative; the same held in `XD_2026`, where all 12
lipids stayed inside their stated Fig 5D band on tops that were 0-21 px wrong.
Geometry has to be checked geometrically.

**Better than a printed statistic: a round number the authors DREW but never
printed.** Fig 2c has a red dotted reference line; the corrected y-axis
constants put it at 50,020 RLU, 0.04% off a round 5e4, which decided between
three competing calibrations that the printed hit rates could not separate. The
`QL_2026` colorbar's zero landing on the centre of its own border is the same
kind of test. These are external to our pipeline in a way a recomputed marginal
is not, and they are cheap once you look for them: reference lines, thresholds
drawn as dashes, axis limits chosen at round values, a colorbar's own tick
notches.

**Mass cannot see connectivity.** HRMS confirms the molecular formula and
therefore the arm *counts* — it is blind to which nitrogen an arm sits on. For
regiochemistry you must read the described synthesis. The Lin SI's
protection/deprotection route is what revealed that arms go one-per-nitrogen,
not two-on-one; HRMS passed 4/4 either way.

**Flag, do not guess — ON THE ROW.** When the source is genuinely silent the
deliverable is a flag, not a structure. A combinatorial library made from crude
products really is a positional-isomer mixture; saying so is correct, not timid.

Write the flag into the **emitted row**, in the columns of
`provenance/schema_quality.json` (§3). A count in your folder's QC JSON does not
reach anyone: four of five curations computed their flags, wrote them to a QC
file, and emitted rows without them, so 1,152 corpus rows carried a convention
silently (F-001). Carry the per-item assessment through; never write a
`confidence` the upstream record already states, because a hard-coded literal
once relabelled 24 rows reading *"ASSUMED: the paper draws no structure"* as the
project's highest confidence (F-011).

**A log axis's decade gaps must be EQUAL — but check that on the ticks you
measured, not on numbers in a comment.** I marked a finding confirmed because a
comment listed anchors whose decade gaps were 406 and 408; they were the tick
positions rounded to integers, and the code was right while its comment was
wrong. Verifying a claim about code by reading a comment is the same error the
comment was an instance of.

**Create the C10 anchor before you call it done.** Find the structure the paper
draws in full — usually the lead compound (Siegwart Fig 2a) or the final product
of a "precise synthesis" scheme (Lin Fig S17, where a protection/deprotection
route also pins the regiochemistry). Assemble it by fragment, **not** string
concatenation, store it in `curation/<paper>/reference_structures.csv`, and run
`provenance/structure_check.py`. If the paper draws nothing, store an empty
`reference_SMILES` with a note so the gap is visible.

**Emit in the 72-column schema** and let `provenance/build_dataset.py` merge it
(add the path to `MERGE`). It reassigns Index/LNP_ID and recomputes descriptors.
Then update `lit_review/build_data.py`'s `ADDED` registry so the tracker shows
the paper as `added`, and rerun `python3 -I lit_review/build_data.py`.

---

## 6. Verification policy

Standing instruction from the user (2026-10-29). Two curated papers both
contained structural bugs that every check at the time was blind to, so
verification is a required stage, not a courtesy.

### The stopping rule is NOT "review until a reviewer finds nothing"

That rule is worse than it sounds. A clean report has two causes — there are no
bugs, or that reviewer could not see them — so stopping at the first clean
report systematically stops at the *weakest* reviewer in the sequence. It
selects for reviewer failure and reads it as correctness, while costing ~250k
tokens a pass.

Look at what actually found our bugs. Not one was found by looking harder:

| bug | why it survived | what found it |
|---|---|---|
iPhos apex counting | raster estimate, no ground truth | adding HRMS anchors |
Lin regiochemistry | molecular formula identical either way | reading the SI synthesis route |
Siegwart cholesterol attachment | O-count unchanged, only 1 H moved | a reviewer comparing to the *drawn* structure |
Siegwart radical carbon | parsed, canonicalised, descriptors self-consistent | same |

Each survived because **no check could see its class**. Coverage of error
classes is the lever; number of passes is not.

### Tier 1 — automated, every build, free

```bash
python3 -I provenance/build_dataset.py   # rebuild
python3 -I provenance/audit.py           # C1-C14, one command; EXITS NON-ZERO on failure
python3 -I provenance/test_checks.py .   # do the checks still have teeth? (slow, ~20 min)
```

**`audit.py` now gates.** Until 2026-10-30 it did not: `main()`'s return value
was discarded and the exit code came only from the C10/C11 delegates, so C1-C9
printed whatever they liked and the command exited 0. Eight simultaneous seeded
violations passed; of 43 mutations, 30 appeared in the report and 5 failed a
command (F-012). Every session, this one included, had been reading a clean
report that was hardcoded.

The gate is **`provenance/audit_baseline.json`**, not a blanket zero, because the
corpus has conditions we accept on purpose — 48 duplicate SMILES groups, 212
unassigned stereogenic C=C, 2,544 blank PMIDs on a preprint, 3 `SB_2024` rows
whose ratios do not sum to 100. Each metric is recorded there **with the reason
it is accepted**, and a departure in *either* direction fails. A decrease fails
deliberately: a known condition vanishing without anyone editing the baseline is
how a check that stopped looking announces itself, and that failure mode reads as
a clean report.

**When the gate goes red, do not edit the number to make it green.** Establish
what moved it, then change the baseline and its `why` in the same commit as the
change that moved it, and say in `CHANGELOG.md` which finding motivated it.

`provenance/test_checks.py` is the other half: it seeds one mutation into one row
of a *copy* and requires a non-zero exit, plus a control run that must still
pass. The two fail in opposite directions — a baseline with no teeth passes
everything, teeth with no baseline fail every clean run. Its first run found two
coverage gaps: gating only `blank_Publication_PMID` let a blanked `IL_name`
through, and **C8 was checking 7 of the 14 published descriptor columns**, so a
`LogP` of 42 passed. Add a mutation for every new check.

`provenance/ERROR_CLASSES.md` is the ledger: **every bug any reviewer finds gets
a row and, where possible, a permanent check. It only grows.** That is what
turns this from a treadmill into a ratchet — a class costs tokens to discover
once, then ~0 to prevent forever.

**C10 is the important one.** It compares each enumerated structure against the
structure the paper DRAWS (`provenance/structure_check.py` +
`curation/<paper>/reference_structures.csv`). It is the only Tier-1 check that
tests topology against external ground truth; C1-C9 all test internal
consistency, which is exactly why connectivity bugs sailed through them.

Per paper, read the drawn lead compound ONCE, assemble it **by a different code
path** than the enumerator, and store it. After that the check is automatic. If
a paper draws no full structure, store a row with an empty `reference_SMILES` and
say so, so the *absence* of an anchor is visible.

**"A different code path" has to be executed, not asserted.** The Siegwart
anchor's note claimed independent fragment assembly "through the C3 oxygen
(3-beta)" and the string it annotated was 3α — byte-identical to the
enumeration, same sha256, 223 bytes. C10 reported `exact_match` while comparing
a value to itself, on the project's flagship curation, for the one molecule it
exists to protect (F-013). So: write a `make_reference.py` in the curation folder
that *builds* the anchor, and prefer a route that shares as little as possible
with the enumerator — the replacement uses plain string substitution with no
`RWMol`/`AddBond`/`RemoveAtom` anywhere, and takes the sterol verbatim from
`reagent_registry.json` so the provenance claim is mechanically true. C10 now
reports `anchor_may_be_a_copy` when a stored reference is byte-identical to the
enumeration *and* the folder has no generator script. Seven anchors currently
trip it.

**Fragment assembly is not the safe option — it is a different hazard.** It was
mandated here after string concatenation fused piperazine to cholesterol via a
ring-closure collision; the replacement then inverted C3 on 537 rows, because
`RWMol.RemoveAtom` drops an atom out of a stereocentre's neighbour list without
correcting its chiral tag. Both routes need an **external, stereo-aware** check.
Where you must delete an atom, delete one that is not adjacent to a stereocentre
— for an esterification that means the acid loses its OH and the alcohol's
oxygen is retained, which is what the chemistry does anyway. And never bond to
"atom index 0 of the fragment": tag the attachment atom with an atom map.

### Tier 2 — one adversarial reviewer per paper

One well-briefed pass beats three generic ones. Brief it with:

1. **Attack, don't confirm.** "A clean report is only useful if you genuinely
   tried to break it."
2. **The live error-class list** from `ERROR_CLASSES.md`, especially the
   "not yet automatable" half — that is where it can add value Tier 1 cannot.
3. **What the automated checks cannot see** — say outright that mass and formula
   are blind to connectivity, or it will report "HRMS passes" and stop.
4. **Re-derive independently**: its own geometry, its own ≥10 spot-checks read
   off the rendered figure, not a rerun of our scripts.
5. **Read-only.**
6. Fixed report structure: CONFIRMED CORRECT (with numbers) / ERRORS FOUND
   (file, what, evidence, correct value) / UNVERIFIABLE (and why) / SEVERITY
   (changes values? changes structures? cosmetic?).

Tell it the tools: `pdfplumber`, `pymupdf`, `PIL`, `numpy`, `scipy`, `rdkit`,
`python3 -I`, and that it must render pages and *look at them*.

### Tier 3 — conditional, not routine

A second reviewer **only when Tier 2 found something structural**, and it
reviews **the fix plus whatever the first pass did not cover** — not the same
ground again. Two reasons: a structural find means that area was under-verified,
and the fix is new unreviewed code. Concretely, the Siegwart cholesterol fix
needed a second look because the same bug existed in a *second copy* inside
`curation/build_new_curation.py`, which re-reviewing the original file would
have missed.

### Tier 4 — a corpus-wide review every 3-5 merged papers

Standing instruction from the user (2026-10-30), after the fourth consecutive
curation turned out to contain a real error: *"since we keep on finding
mistakes… set up chains of subagents that go into every single paper we've added
so far and perform detailed quality control and validation."* Procedure and
round log live in **`provenance/review/REVIEW_PROGRAM.md`**.

**Trigger:** every 3-5 merged papers, and always before the frozen pipeline is
re-run or a number goes into the manuscript.

**Why an interval and not just per-paper.** Tier 2 runs at merge time and
catches what one paper's own material can reveal. What it structurally cannot
see is anything visible only *across* papers: a convention that drifted between
curations, a documented property the schema never had, a check that stopped
having teeth three merges ago. Those accumulate silently, so the unit of review
has to be the corpus, not the paper. Round 1 confirmed the point, expensively. Its four largest findings were all
invisible to a per-paper reviewer: the flagging posture this very section
promised was **not implemented in any column** (1,152 rows, F-001); **all 537
`ZT_2026` rows** stored epicholesterol because fragment assembly inverted a
chiral tag that no composition check can see (F-013); **5,040 rows across the
corpus's only three multi-cell-type experiments** — all three of them — were
standardised per cell type, deleting the between-assay offset (F-016); and
**the Tier-1 audit could not fail at all**, because `main()`'s return value was
discarded (F-012). Full ledger in `review/round_2026-10-30/FINDINGS.md`,
repairs in `CHANGELOG.md` §12.

**Shape:** three waves, documented in full in the program file.

| wave | agents | scope |
|---|---|---|
| 1 | one per curated paper | re-derive structures, values and metadata from the PDF/SI/repo, not from our scripts |
| 2 | three | `build_dataset.py`'s patches and precedence rules; the docs-vs-data contract; mutation-testing whether C1-C11 still fire |
| 3 | me, then a validation agent, then a post-fix agent | adjudicate, confirm the repair is the *right* repair, review the diff |

**Waves 1 and 2 are read-only.** The only file an agent may write is its own
report. A reviewer that edits destroys the independent second opinion on its own
change, and two reviewers editing one file in parallel corrupt it.

**Brief from a generated packet, never from prose.** `provenance/review/brief.py`
reads every figure out of `data/` and `curation/` at generation time. A prose
brief goes stale at the next merge, and a stale brief is worse than none — it
makes a reviewer confidently verify a number that is no longer in the data.

**Validation before any change** (Wave 3). Nothing is applied because a reviewer
said so: I reproduce the evidence myself, a fresh agent then judges whether the
proposed repair is correct *and whether it is the right half to repair*, and only
then is a file edited. That second step exists because of `AC_2025`, where the
first repair was confident, internally consistent, and assigned 419 measurements
to the wrong molecules.

### When to stop

Not "a reviewer found nothing", but:

> every known error class has an automated check **and** one adversarial pass
> found nothing outside those classes.

If a reviewer only re-finds classes Tier 1 already covers, coverage has
saturated — stop. If it finds a new class, add the check and go again.

### Treat reviewer findings as claims

Subagents are wrong too. Of the Siegwart reviewer's 12 findings, 3 were real
structural errors, 1 was half-right (a misleading column name in an intermediate
file, not a corpus error), and 1 was flatly wrong (it said the dose contradiction
was unrecorded; its grep scope had missed the file). Verify before acting, and
say plainly which findings survived.

### Residual risk, stated honestly

This does not reach zero. Some things are unverifiable from the paper alone —
which nitrogen reacts in Lin's C1-C4 cores, the exact branch of iso-decyl
acrylate. The posture for the manuscript is to state the verification procedure
and flag those rows in the data, so a reader can see which values are measured
and which are convention.

**This is now implemented** — 7 columns, §3 — having been documented here and
nowhere else for two days while 1,152 rows carried a convention silently
(F-001). The lesson is narrower than "write the code": the flags had been
*computed*, and were sitting in each curation's own QC JSON. They died at the
emit step, one stage before the merge, so no change to `build_dataset.py` could
have rescued them and nobody reading the builder would have found it. **A flag
that lives only in a per-folder QC file is not a flag, it is a memo.**

---

## 7. Working conventions

- **Reproducibility first.** One frozen script (T3) regenerates all numbers. When a number is finalized, update §3 / the manuscript AND note it here. Never hand-edit a reported number.
- **Seeds & versions pinned** and recorded (RDKit, scikit-learn, XGBoost/PyTorch). State them in Methods and in the script header.
- **Source key = `Experiment_ID`** everywhere. Chance baseline for S-way classification = 1/S (1/31 vitro, 1/17 vivo).
- **Featurization:** Morgan/ECFP (draft says radius 2, 2048 bits) for source-classification & Tanimoto NN. Confirm and record exact params in T3.
- **Primary metrics:** pooled Spearman vs within-experiment scaffold-disjoint Spearman — always be explicit which one a figure/number uses (the draft's main ambiguity).
- **Scratch/tmp files** go in the session scratchpad, not in `data/` or `scripts/`.
- **Never silently repoint a frozen script at new data.** Frozen numbers belong to the corpus they were computed on.
- **Verification is a stage, not a courtesy** (§6): Tier 1 automated every build, Tier 2 one briefed adversarial reviewer per paper, Tier 3 only if Tier 2 found something structural.
- **Every bug found becomes a permanent check** in `provenance/ERROR_CLASSES.md`. That file only grows.
- **A corpus-wide review runs every 3-5 merged papers** (§6 Tier 4, procedure in `provenance/review/REVIEW_PROGRAM.md`). Per-paper review cannot see cross-paper drift. Reviewers are read-only; a validation agent checks every repair before it is applied.
- Inclusion/exclusion criteria are meant to be **a priori, frozen, metric-blind** — document them in Methods when curation is finalized.

---

## 8. Open questions / decisions needed

- **`molea.csv`** = a **published combinatorial ionizable-lipid mRNA-transfection screen** (1,152 rows; `Markush code` = Linker×Head×tail enumeration, lipid-like SMILES). Targets: `mTP_(Log2)_Hepatocyte`, `mTP_(Log2)_Chondrocyte` (mTP = mRNA Transfection Potency, log2 scale, range ≈ −6…+14, per-source-uncentered), + `mTP_Difference` (a selectivity contrast — not a poolable readout). Hepatocyte = liver cell, chondrocyte = cartilage cell — just two cell types (the corpus already spans 12).
  - **Comparability:** in KIND it IS comparable — same domain (ionizable-lipid mRNA reporter transfection) as the corpus `luminescence_normalized` rows; poolable in principle after adding formulation features + per-cell-type z-scoring. Unknown to pin down for Methods: exactly how mTP was computed (single-dose log2 luminescence vs log2 fold-vs-control vs dose-response fitted potency).
  - **DECISION (2026-07-23, refined 2026-07-24): keep SEPARATE, do NOT pool.** Pooling buys ~nothing (zero canonical-SMILES overlap: 0/11,459 vitro, 0/1,434 vivo → can't feed T1; just one more disjoint silo). Its real value is STANDALONE as a **within-lab, cross-context concordance** point: identical lipids, same lab/library/experiment, only the cell type changes → **Hepatocyte vs Chondrocyte Spearman = 0.30 (Pearson 0.38, n=1,152)**. This is a middle rung between within-experiment learnable (~0.37, Result 6) and cross-lab shared-lipid (~0, Result 3): even one lab loses most rank agreement when only the target cell changes. Strengthens Result 3 and pre-empts the "cross-lab ≈0 is just assay noise" objection with a clean same-lab single-variable degradation.
  - Roles: (a) within-lab cross-cell-type concordance (companion to Result 3); (b) T2 reproducibility example; (c) Result 6 within-experiment regime. Single-source ⇒ NOT a T11/LOSO candidate.
- **Corpus provenance for Methods:** LNPDB-derived + which additional publications? Need the source list with citations (the `〘cite …〙` markers) — `Papers/` likely holds them.
- **Endpoint-comparability grouping** (§3 warning): define property type / readout / biological system / route thresholds → report honest comparable-environment counts for the in-vitro and in-vivo pools.
- **`Experiment_value` normalization:** confirm per-source scale/units and the within-source normalization used before any cross-source rank comparison.
- **Re-run the frozen pipeline on the 25,997-row corpus?** `scripts/` is pinned to the old 19,028-row files (now in `provenance/superseded/`). Re-running changes every number in `results/`, `reports/`, and the manuscript. Decide deliberately, then update §3/§4 together.
- **`SB_2024` / `SY_2025` (60 vivo rows) are not on the z-score scale** — drop them as we dropped `ZL_2022`/`YX_2023`, or keep and caveat?
- **48 duplicate SMILES with conflicting values** — resolve from source publications, or exclude those lipids from Result 3?
- **The in-vivo pool's only cross-source chemistry is now MC3 and SM-102.** `XD_2026` carries both, taking MC3 from 4 → 5 sources and SM-102 from 3 → 4, but it adds no *new* shared lipid (ALC-0315 is vitro-only) and shares at most 2 lipids with any single source, so its contribution lands in Stage 5's pooled set and not in any per-pair rho. Flagged because those two are the field's universal positive controls and they rank 2nd and 3rd of 12 here: a source whose only shared lipids are benchmarks that top every screen can push cross-lab concordance toward +1 for reasons that have nothing to do with transferable SAR. Decide before Stage 5 is re-run whether benchmark-only pairs belong in the pooled concordance at all — this bears directly on T8 and on Result 3.
- **git init + license check:** constituent-source licenses must be CC0/CC-BY-compatible before corpus release (Data & code availability section).

---

## 9. Changelog

- **2026-10-29** — **`XD_2026` merged: 12 in-vivo rows** from Duan et al., Acta Pharm Sin B 2026 (PMID 42592089), the LipidAI paper — E1–E9 plus ALC-0315/SM-102/MC3, BALB/c IV, whole-body IVIS luciferase at 6 h (Fig 5). `data_vivo.csv` 1,768 → **1,780 rows, 19 experiments**. Structures were traced **geometrically** from Fig 5A (`stage02_trace.py`: skeletonise → trace → Douglas-Peucker → count vertices), never counted by eye, with an annotated overlay per structure. **The benchmark control found a real error in the paper**: Fig 5A draws ALC-0315 with a 3-carbon hydroxyalkyl head where the compound has 4, while SM-102 and MC3 in the same panel reproduce the corpus exactly — and E1, drawn in the same style 390 px below, *does* have 4, so it is a localised authoring error, not a convention. The three benchmarks are therefore emitted with corpus SMILES and **E1–E9 carry `confidence = medium`**; they have no external anchor of any kind (supplied by a vendor, no name/formula/HRMS/SI), so `reference_structures.csv` holds an empty row for each and C10 records the ALC-0315 override as `expected_mismatch_paper_drawing_is_wrong` rather than as prose. Values came from two figures that agree: Fig 5C digitised against the 10^6–10^9 ticks (residuals <0.002 decades) lands inside its own Fig 5D A/B/C band **12/12**, and Fig 5D reproduces both printed numbers (75% overall, class B 1/3). `Model_target`, `Route` and `Dialysis_buffer` are recorded in `07_qc.json` as *readings*, not printed facts; dose and N:P are absent from the paper and left blank. **The repo's 485-row LipidDB was deliberately NOT merged** — it is a compilation across the 47 studies the paper cites, so one `Experiment_ID` would fabricate a lab; it is extracted standalone in the curation folder as a second, independent literature-pooled corpus. A Tier-2 adversarial reviewer re-derived every structure, value and metadata field at 300 dpi native and confirmed all of them; its nine findings were all verified and all fixed, none touching an emitted structure or value. Its best catch was that the LipidDB reference strings do not de-duplicate on a `doi:` cut — the authors append a per-row counter to a DOI *sometimes*, and to a PMCID or a page range otherwise. **Its replacement count of "299" was itself wrong** and was propagated into `CHANGELOG.md` before being checked: the curation folder's own `06_lipiddb_qc.json` records 56 / 50 against the paper's stated 47, and the round-1 review found nine defensible de-duplication rules giving 50-116. No single count is asserted now; see `CHANGELOG.md` §11. The lesson is the one §6 already states — **treat reviewer findings as claims** — and it was not applied to this one.


- **2026-10-29** — **Verification policy implemented (§6).** Rejected "review until a reviewer finds nothing" as the stopping rule: a clean report can mean no bugs or a weak reviewer, so that rule stops at the weakest reviewer in the sequence. Replaced with a three-tier policy whose stopping rule is *coverage*, not passes. Added `provenance/ERROR_CLASSES.md`, the ratchet ledger — 17 classes now carry an automated check, 5 are flagged as needing eyes. Added **check C10** (`provenance/structure_check.py`): each enumerated structure compared by InChIKey against the structure the paper DRAWS, reported as skeleton (connectivity, fails) vs stereo (flags). This is the only Tier-1 check testing topology against external ground truth; C1–C9 all test internal consistency, which is why three connectivity bugs passed them. Anchors created for both curated papers — Siegwart C6O2B2 (Fig 2a) and Lin 10A2-3B-EP10 (Fig S17, regiochemistry pinned by its protection/deprotection route) — both **exact matches**, each assembled by a different code path than the enumerator. Building the Siegwart anchor immediately caught a SMILES ring-closure collision that had fused the piperazine to cholesterol, which is why the playbook now says fragment assembly, not string concatenation. `audit.py` runs C10 as its final step so one command covers Tier 1.

- **2026-10-29** — **CLAUDE.md brought back in sync** (it had drifted since 2026-07-25: §2/§3 described files that no longer exist). Added §5 (curation playbook) and §6 (subagent verification as the default). Flagged that `scripts/` is pinned to the superseded corpus files.
- **2026-10-28** — **Corpus rebuilt and finalized: `data/data_vitro.csv` (25,997 x 72, 40 experiments) + `data/data_vivo.csv` (1,768 x 72)**, reproducible via `provenance/build_dataset.py`. Lineage = original LNPDB (vitro) + the ML change log + Bowen Li additions, with three deletions **reversed** at the user's instruction (`YX_2024` RAW264.7 1,200; `XH_2024` 30; `ZC_2023` 100 — all removed for *modelling* reasons, not data quality). Kept deleted: `ZL_2022`/`YX_2023` (z-scored without log) and 77 OOD cell-line rows. Precedence is explicit where the two upstream files disagree. Restoring `YX_2024` needed no rescaling — checked, not assumed: the two sources' z-scores are bit-identical (slope 1.000000, residual 7e-15). **676 rows had descriptors inconsistent with their own SMILES**; all recomputed. `van.der.Waals.Molecular.Volume` renamed to `ExactMolWt` after verifying on 40,472/40,542 rows. Eleven experiments were off the population-sd z-score convention and were rescaled. Two structure patches traced to source: `4A3-SC6` (held the SC8 octyl structure; series pattern pins MW 1354.99) and 42 `KZ_2016` `*-SC3` lipids (held the SC2 isopentyl structure; source metadata pins neopentyl). A fresh official LNPDB download is byte-identical to the archived reference (MD5 `32b5187...`), confirming the `NA` values and the SC6/SC8 bug are upstream.
- **2026-10-28** — **`provenance/audit.py` added**: an independent C1-C9 audit, deliberately not a rerun of the builder's own checks. Four of its first findings were the audit itself being wrong (name tokens that are not chain lengths; a partial replicate key; terpenoid alkenes that have no E/Z to lose) — the corrected scoping is what makes the clean result meaningful. Current state: 0 unparseable, 0 non-canonical, 0 duplicate rows, 0 descriptor disagreements across 18,066 unique lipids.
- **2026-10-28** — **Two papers curated and merged.** `ZT_2026` (Tian, Adv Mater 2026, PMID 42634137): 486-formulation DoE on one lipid + 51 lipids in vivo; 7 independent checks incl. 15/15 marginal bins exact. `LL_2026` + `LL_2026_2` (Lin, JACS 2026, PMID 41689530): 840 aromatic ionizable lipids + a 147-lipid initial screen at a different formulation; HRMS 4/4 exact to 0.0001 Da and 48/48 printed hit rates reproduced. **A regiochemistry bug was found and fixed in the Lin enumeration** — aromatic arms go one-per-nitrogen (pinned by the SI's protection/deprotection route), not two on one; the molecular formula is identical either way so HRMS could not detect it. This is the reason §6 exists.
- **2026-10-27** — **Literature sweep provenance rescued** into `lit_review/provenance/`: the 20 verbatim PubMed queries, raw/deduped PMID lists, the corpus map used for duplicate checking, and the XML parsers (including the buggy one, kept to document the bug). README states the sweep's four known blind spots.

- **2026-07-23** — Project set up. Profiled corpora; settled T8 (in-vivo → supplementary) from the 2-shared-lipid finding. Identified `molea.csv` as single-experiment hepatocyte-vs-chondrocyte training data (single-source; not a T11 candidate — a within-experiment-regime / T2 split-half example).
- **2026-07-23** — Switched canonical in-vitro corpus to `LNPDB_vitro.csv`, trimmed from the ORIGINAL `LNPDB.csv` (the `_processed` file was from a different modeling project). Trim = `Model=="in_vitro"` minus 90 `hemolysis_percent` rows → 17,050 rows, **S=33** (was 31), 11,455 lipids, **6** overlapping source-pairs (~3 substantive; was 12). All §3 in-vitro numbers updated and are the new pre-T3 baseline.
- **2026-07-23** — Re-derived in-vivo from `LNPDB.csv` → `LNPDB_vivo.csv` (same delivery-only definition, excl. diameter/zeta/uptake) = 1,717 rows, S=17, 1,436 lipids, 2 shared, 4 overlapping pairs (each 1 lipid). Row count exactly matches old `_del` file; supersedes `LNPDB_vivo_del.csv`.
- **2026-07-23** — Added `AC_2025` ("comet_lance", 1,978 rows) to `LNPDB_vitro.csv` from the processed file (it's absent from original `LNPDB.csv`). In-vitro now **19,028 rows, S=34**, 11,459 lipids, 230 shared, 9 overlap pairs. AC_2025 mapped into 66-col schema with blank decomposition metadata + RDKit-computed MolWt. This is the only source in the in-vitro corpus NOT sourced from `LNPDB.csv`.
- **2026-07-24** — Characterized `molea.csv`: it's a real ionizable-lipid mRNA-transfection screen (endpoint comparable in kind to the corpus), but confirmed decision to keep it SEPARATE. Key new finding: within-molea Hepatocyte-vs-Chondrocyte Spearman = **0.30** (identical chemistry, same lab, only cell type changes) — promotes molea from "reproducibility example" to a within-lab cross-context concordance point supporting Result 3. See §6.
- **2026-07-25** — **T3 done: full frozen pipeline built (Stages 0–9, `scripts/` + `run_all.py`).** `git init` run (no commits yet). Fixed two silent data-format bugs found in Stage 0 (not exclusions): `has_ester/has_carbonate/has_disulfide` serialize inconsistently across the raw CSVs (native bool vs. "TRUE"/"FALSE" string), and `Publication_PMID` reads as float64 in files with missing PMIDs vs. int64 in files without — both normalized before any content-hash/dedup comparison. Verified `LNPDB_vitro ∪ LNPDB_vivo` is an exact content-clean partition of `Original_LNPDB.csv` (1,030 excluded = 90 hemolysis + 269 blank-Model + 671 additional in-vivo assay exclusions; 1,978 added = AC_2025, absent from Original — reconciles exactly). Stage 1 curation excluded 328 in-vitro rows for missing `Experiment_value`; zero exact-duplicate rows found in either pool; final frozen corpus = 18,700 in-vitro / 1,717 in-vivo rows. All T1/T2/T4/T5/T6/T7 closed as a consequence (see §4); full headline numbers in `reports/SUMMARY.md`. Found and fixed one real analysis bug along the way: a naive Bemis-Murcko "congeneric family" check in Stage 5 was collapsing all acyclic (ring-less) shared lipids into one fake family via the shared empty-string scaffold — fixed to exclude acyclic lipids from that count and report them separately (163/231 shared vitro lipids are acyclic). Split-half negative control (Stage 8) **PASSES** (separability collapses to 0.496/0.490 vs. chance 0.500) — the pipeline's separability/concordance machinery is validated as measuring real chemistry, not producing signal from nothing. Full `run_all.py` run is byte-identical across two independent executions (~41 min end to end on this machine); one open item is `git_commit` in the manifest is null pending the repo's first commit.
- **2026-07-25** — **Added a third Stage 6 objective: graded-relevance LambdaRank (`scripts/lambdarank.py`).** Ported (loss math, pair sampler, gradient/hessian — all verbatim) from the companion modeling+wet-lab-validation project's `within_exp_lambdarank2.py` (found at `~/Downloads/LNP_ML/scripts/`, 4 identical copies across its subprojects). That project's `rel` labels come from wet-lab-curated hit flags (rel==3 = curated hit, ~3% of rows) via a `label_rel.py` script that does **not** exist in the download — our LNPDB corpus has no equivalent curated-hit column. Substituted: `rel` derived from the within-source-normalized rank already computed in Stages 5-7, with the **top 10%** within each source (not that project's ~3%) becoming "hits" — documented as a deliberate choice (our smaller per-source lipid counts would leave most sources with zero hits at 3%), easy to change via `LAMBDARANK_HIT_FRAC` in `stage06_protocol.py`. Ported self-test (finite-difference gradient check, grad-sum-to-zero, gauge invariance, non-negative Hessian, dZ sanity) passes exactly. Result: LambdaRank sits **between** the existing two objectives on vitro (LOSO pooled Spearman 0.130, vs. 0.093 pooled-regression / 0.158 within-source-pairwise) but is **worse** than both on vivo LOSO (0.021 vs. 0.233 / 0.090) — the added sophistication (graded relevance + hit-anchored sampling) does not reliably beat the simpler within-source pairwise ranker on this corpus, particularly in the smaller in-vivo pool. Not a dropped-experiment artifact (only 1 vivo source, `SB_2024`, excluded — already known to have just 2 unique lipids). Full grid now 3×3; `reports/06_protocol.md` and `results/stage6_protocol_grid.csv` updated.
- **2026-07-25** — **Added top-k ranking-quality metrics vs. a random baseline to every Stage 6 grid cell**, ported from the same sibling project's `analyze.py` evaluation pattern (graded NDCG@k_e/full/{3,5,10} and hit_rate@{3,5,10}, each paired with a 1000-permutation random-ordering floor, seed=0 fixed a priori; EF@k reported as a ratio of means, not mean of per-fragment ratios, per that project's own documented rationale for stability). Unit of analysis is an (experiment, fold) fragment — LOSO gives each source intact ("primary" for hit metrics); random/scaffold-disjoint k-fold fragments a source's rows across multiple folds' test sets ("secondary, partial admissibility" — reported honestly via an admissible-fragment count, never hidden). New output: `results/stage6_topk_metrics.csv`. **This surfaces a striking result Spearman alone missed:** vivo LOSO hit_rate@5 EF for pooled regression is **0.96x — literally sub-chance** at finding the best-performing lipids, while both within-source objectives (pairwise, LambdaRank) hold at ~1.09x. LambdaRank's overall rank correlation collapses hardest of the three under vivo LOSO (Spearman 0.021), yet its hit-finding EF@5 matches the pairwise ranker's — the rare-hit-focused loss preserves some ability to surface top performers even where general rank agreement breaks down. Vitro shows the same NDCG/EF monotone-collapse pattern as the Spearman numbers (e.g. pooled-regression EF@5: 3.42x random → 2.45x scaffold-disjoint → 2.12x LOSO).

- **2026-10-30** — **`QL_2026` merged: 344 in-vitro rows** (Liu et al., Adv Healthc Mater 2026, PMID 41145968) — 344 thiosuccinimide lipids in HeLa from a bispecific-linker library, 320 from Fig 2b/2c plus a 24-member six-tail library from Fig 3c. `data_vitro.csv` **28,541 → 28,885 rows, 41 → 42 experiments.** The four benchmark LNPs on the same panels are calibration, not library members, and are not emitted. One pooled z after a **1.1116×** scalar offset (sd 0.0234) tied library 2 to library 1 through their four shared benchmarks — standardising the two libraries separately had deleted the paper's central ~37× library difference (154,630 vs 4,128 RLU geometric mean) and inverted 36% of cross-library pairs while every check still passed. Its PMID had been recorded as blank "not indexed"; it was indexed all along, so **look the PMID up before recording a paper as unindexed.** Full detail in `provenance/CHANGELOG.md` §11, which also back-fills the `XD_2026` entry that section had been missing.
- **2026-10-30** — **Round-1 corpus review: repairs applied. Corpus 28,885 x 79 / 1,780 x 79; row counts unchanged.** Eight read-only reviewers over five curations, 21 findings, each repair sent to a separate validation agent first — and two validators changed the plan materially, which is the case for the step. Full detail in `provenance/CHANGELOG.md` §12 and `provenance/review/round_2026-10-30/FINDINGS.md`. The four that mattered: **(1)** All **537 `ZT_2026` rows** stored **epicholesterol** (3α). `RWMol.RemoveAtom` deleted cholesterol's hydroxyl out of C3's neighbour list without correcting C3's chiral tag. The epimers share molecular formula, exact mass to four decimals, the InChIKey *skeleton* block and all 17 descriptor columns, so C1-C9 were blind; C10 was blind because its anchor had the same sha256 as the enumeration it checked, while its note claimed independent assembly "through the C3 oxygen (3-beta)"; C11 was blind twice, comparing skeleton blocks only and inspecting `CHL_*` columns that are empty on every one of these rows. The validator settled 3β from Fig 1b's wedges, the 4.71-4.58 axial C3-H, and a ¹³C peak at 73.9 ppm — **and found a live second copy of the same bug in `build_new_curation.py` affecting 15 more structures**, so fixing one file would have left 15 wrong after the fix was announced. That copy is deleted; the file imports the one implementation. **(2)** **The Tier-1 audit could not fail** — `main()`'s return value was discarded, eight seeded violations exited 0, and `build_dataset.py` never imported `sys`. Now gated by `audit_baseline.json` (62 metrics, each with the reason it is accepted, failing on departure in *either* direction) plus `test_checks.py`, a seeded-mutation suite whose first run found that **C8 was checking 7 of 14 descriptor columns**. **(3)** **`WT_2026`'s 2,544 rows re-standardised**, 14 blocks to 7. Z-scoring each (file × cell type) column separately forced the authors' own published `nLnLE_Selectivity` to identically zero — raw HepG2 − Ramos is **+4.0342 ln (56.5×)** on file 5 with 173/180 signs one way, stored as 0.0000 with Spearman(raw, stored) = **−0.517**, anti-correlated. One z per file restores it and is now *checkable*: `z_a − z_b = Selectivity/σ`, asserted at max residual 1.8e-15. **(4)** **1,152 rows carried a convention and said nothing about it** (F-001) — the flagging posture §6 promised existed in no column anywhere, and the flags died one stage *before* the merge, so no builder change could have saved them. Schema 72 → **79**, defined once in `provenance/schema_quality.json`. F-011 had to land first: a hard-coded `"high"` in `build_new_curation.py` had relabelled the 24 lipids whose record reads *"ASSUMED: no ring N and the paper draws no structure"* as the project's highest confidence, so publishing `confidence` before fixing it would have been worse than publishing nothing. New checks **C12** (error class 22, with its specificity tested synthetically because the corpus had no correctly-pooled multi-cell-type experiment to serve as a negative control — all three were affected), **C13** (flag vocabulary), **C14** (cleaves every ester in `IL_SMILES` corpus-wide and pins the liberated sterol's *full* InChIKey against the registry). C14 earned itself on its first run by finding **F-021** in `ZH_2023`, an experiment no reviewer was assigned: a steroid that cleaves to C26H44O, missing the C19 angular methyl — not cholesterol, still open. Also fixed from two independent re-measurements: **F-019**, where `QL_2026`'s colorbar was anchored on the OUTER EDGE of its own black border (a second class-20 instance, on a colorbar) biasing every Fig 2b value high by ~1% at the hit threshold and **+38% at the floor** — the 48 white cells were stored at 1,310 RLU when their whole feasible band is [0, ~982], i.e. **above** the `value_is_bound = upper` flag this round had just given them; **F-018**, where a bar-top detector read replicate dots on 15 of 60 bars, which also means the stored cross-library offset of 1.1116 was largely measuring our own error — three of its four benchmark anchors were contaminated, and the corrected **1.0071 ± 0.0247 is not distinguishable from 1.0**; **F-020**, the same class in `XD_2026` (12/12 tops, mean +0.051 decades); and **F-010**, C2's coverage, 235 → **634** series across 11 experiments instead of 1, by testing per series whether the name token equals the measured chain length rather than keeping an allowlist. **F-017 I had confirmed WRONG**: its code was right and its comment was wrong — the anchors whose decade gaps looked unequal were tick positions rounded to integers. Error classes **29-36** added. Five findings this round were conclusions drawn from a secondary artifact rather than the primary one (a reviewer's count, an intermediate file, a comment, my own reimplementation of a mask); all five were wrong, and `REVIEW_PROGRAM.md` §5 now says to run the code's own objects and read the data itself. Corrections to my own earlier claims: F-016's "3 of N multi-cell experiments" is **3 of 3** — the corpus has only three; and F-011's mechanism was a hard-coded literal, not the join I had diagnosed.
- **2026-10-30** — **Review program established (§6 Tier 4).** A corpus-wide, three-wave subagent review now runs every 3–5 merged papers; procedure, attack brief and round log in `provenance/review/REVIEW_PROGRAM.md`, per-reviewer packets generated from live data by `provenance/review/brief.py`. Round 1 (5 papers, 8 agents) is under way in `provenance/review/round_2026-10-30/`. Two findings before the reviewers even reported: **(1)** the flagging posture §6 promises is not implemented — the 72-column schema has no `confidence`, `positional_ambiguity` or `*_assumed` column and neither does any curation output, so **≥266 rows across four of five curated papers** carry no flag where one is documented, including 113 of `QL_2026`'s 344 sitting on a colormap floor where identical RGB spans ~980 RLU. *Superseded the same day:* the real count is **1,152**, measured from the emitted corpus (my ≥266 counted only the flags enumerable from QC JSONs; a reviewer's 1,137, which I adopted next, still undercounted `ZT_2026` by 15), and the repair is **applied** — schema 72 → 79, see the round-1 entry above. **(2)** `build_manifest.py` asserted a frozen row-count literal and raised on a correct manifest the moment the merge landed; it now asserts against the live file length with a non-fatal baseline tripwire. Separately, `lit_review/build_data.py` could not run at all — its sweep inputs lived in an ephemeral session scratchpad — so they are recovered into `lit_review/provenance/`, verified by rebuilding and diffing to exactly the one expected row change.
