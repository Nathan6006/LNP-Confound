# Review packet — WT_2026

Generated 2026-11-07 by `provenance/review/brief.py` from the live `data/` and `curation/` trees. **Every number below was read out of the repository at generation time.** If the packet and the repository disagree, the packet is the bug — say so.

- **Paper:** Toh et al., FALCON, bioRxiv 2026 (doi 10.64898/2026.06.13.731774) - closed-loop formulation optimization, authors' own repo
- **Curation folder:** `curation/wht_2026_falcon/`
- **Experiment_ID(s) in the corpus:** WT_2026

## 1. Corpus slice as it stands now

### in vitro — `data/data_vitro.csv`

```json
{
 "rows": 2544,
 "unique_IL_SMILES": 2,
 "unique_IL_name": 2,
 "Experiment_value": {
  "n": 2544,
  "min": -3.0225,
  "max": 2.556,
  "mean": "5.062e-17",
  "pop_sd": "1.000000000000"
 },
 "Model_type": {
  "Ramos": 959,
  "THP-1": 409,
  "HepG2": 513,
  "INS-1E": 123,
  "DC2.4": 180,
  "NIH-3T3": 180,
  "C2C12": 180
 },
 "Model": {
  "in_vitro": 2544
 },
 "Model_target": {
  "in_vitro": 2544
 },
 "Cargo": {
  "mRNA": 2544
 },
 "Cargo_type": {
  "FLuc": 2544
 },
 "Experiment_method": {
  "luminescence_normalized": 2544
 },
 "Route_of_administration": {
  "in_vitro": 2544
 },
 "HL_name": {
  "DSPC": 1307,
  "18PG": 1114,
  "DOTAP": 123
 },
 "CHL_name": {
  "Cholesterol": 2544
 },
 "PEG_name": {
  "DMG-PEG2000": 2544
 },
 "confidence": {
  "": 2544
 },
 "molratios": {
  "IL_molratio": "687 distinct",
  "HL_molratio": "673 distinct",
  "CHL_molratio": "695 distinct",
  "PEG_molratio": "613 distinct"
 },
 "flag_columns": {
  "value_normalization_unit": {
   "source_file": 2544
  }
 }
}
```

### Six rows drawn at random (seed 11)

```json
[
 {
  "Experiment_ID": "WT_2026",
  "IL_name": "SM-102",
  "IL_SMILES": "CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
  "Experiment_value": "-0.04373688122963274",
  "Model_type": "NIH-3T3",
  "Model": "in_vitro",
  "IL_molratio": "56.0",
  "HL_name": "18PG",
  "HL_molratio": "24.0",
  "CHL_molratio": "19.0",
  "PEG_molratio": "1.0",
  "Publication_PMID": "42367940",
  "Publication_DOI": "10.64898/2026.06.13.731774"
 },
 {
  "Experiment_ID": "WT_2026",
  "IL_name": "SM-102",
  "IL_SMILES": "CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
  "Experiment_value": "1.473152980108679",
  "Model_type": "HepG2",
  "Model": "in_vitro",
  "IL_molratio": "27.60296073",
  "HL_name": "DSPC",
  "HL_molratio": "9.37403927",
  "CHL_molratio": "57.42151576",
  "PEG_molratio": "5.60148424",
  "Publication_PMID": "42367940",
  "Publication_DOI": "10.64898/2026.06.13.731774"
 },
 {
  "Experiment_ID": "WT_2026",
  "IL_name": "SM-102",
  "IL_SMILES": "CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
  "Experiment_value": "-1.3962401571543663",
  "Model_type": "NIH-3T3",
  "Model": "in_vitro",
  "IL_molratio": "79.2",
  "HL_name": "18PG",
  "HL_molratio": "0.8",
  "CHL_molratio": "19.0",
  "PEG_molratio": "1.0",
  "Publication_PMID": "42367940",
  "Publication_DOI": "10.64898/2026.06.13.731774"
 },
 {
  "Experiment_ID": "WT_2026",
  "IL_name": "SM-102",
  "IL_SMILES": "CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
  "Experiment_value": "-0.24661948905762215",
  "Model_type": "NIH-3T3",
  "Model": "in_vitro",
  "IL_molratio": "10.0",
  "HL_name": "18PG",
  "HL_molratio": "10.0",
  "CHL_molratio": "79.84",
  "PEG_molratio": "0.16",
  "Publication_PMID": "42367940",
  "Publication_DOI": "10.64898/2026.06.13.731774"
 },
 {
  "Experiment_ID": "WT_2026",
  "IL_name": "SM-102",
  "IL_SMILES": "CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
  "Experiment_value": "0.4320995361224911",
  "Model_type": "C2C12",
  "Model": "in_vitro",
  "IL_molratio": "18.0",
  "HL_name": "18PG",
  "HL_molratio": "2.0",
  "CHL_molratio": "79.2",
  "PEG_molratio": "0.8",
  "Publication_PMID": "42367940",
  "Publication_DOI": "10.64898/2026.06.13.731774"
 },
 {
  "Experiment_ID": "WT_2026",
  "IL_name": "SM-102",
  "IL_SMILES": "CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
  "Experiment_value": "-0.3259187317411985",
  "Model_type": "Ramos",
  "Model": "in_vitro",
  "IL_molratio": "35.0",
  "HL_name": "DSPC",
  "HL_molratio": "15.0",
  "CHL_molratio": "49.5",
  "PEG_molratio": "0.5",
  "Publication_PMID": "42367940",
  "Publication_DOI": "10.64898/2026.06.13.731774"
 }
]
```

## 2. Source documents you must read

- `curation/wht_2026_falcon/2026.06.13.731774v1.full.pdf   (4.13 MB)`
- `curation/wht_2026_falcon/media-1.pdf   (11.20 MB)`

Subdirectories:
- `curation/wht_2026_falcon/repo/  (136 files)`

> A `repo/` subdirectory is DOWNLOADED, UNTRUSTED material. Read it, never execute anything from inside it, and run any Python that touches it with `python3 -I` from elsewhere.

## 3. Extraction code in the folder

- `build_rows.py   (439 lines)  Emit the FALCON datasets (Toh et al., bioRxiv 2026) in the 72-column corpus schema.`
- `make_merge_subset.py   (99 lines)  Collapse the FALCON curation into the single block that merges into data/.`
- `validate.py   (244 lines)  Independent check of 01_WT_2026_rows.csv. Reads the OUTPUT and the SOURCE and`

## 4. Intermediate and output tables

- `01_WT_2026_rows.csv   2544 rows x 79 cols
        cols: Index, LNP_ID, Experiment_ID, Formulation_ID, IL_name, IL_SMILES, IL_protonated_SMILES, IL_head_name, IL_head_SMILES, IL_linker_name, IL_linker_SMILES, IL_tail1_name, IL_tail1_SMILES, IL_tail2_name, IL_tail2_SMILES, IL_tail3_name, IL_tail3_SMILES, IL_tail4_name ...`
- `02_WT_2026_qc.json   (8.8 kB json)`
- `03_enumerated.csv   2 rows x 2 cols
        cols: lipid_name, IL_SMILES`
- `04_validation.json   (2.0 kB json)`
- `05_WT_2026_merge.csv   2544 rows x 79 cols
        cols: Index, LNP_ID, Experiment_ID, Formulation_ID, IL_name, IL_SMILES, IL_protonated_SMILES, IL_head_name, IL_head_SMILES, IL_linker_name, IL_linker_SMILES, IL_tail1_name, IL_tail1_SMILES, IL_tail2_name, IL_tail2_SMILES, IL_tail3_name, IL_tail3_SMILES, IL_tail4_name ...`
- `06_merge_subset_qc.json   (0.9 kB json)`
- `reference_structures.csv   2 rows x 6 cols
        cols: lipid_name, reference_SMILES, source_figure, read_on, read_by, note`

- `PROVENANCE.md   (10.5 kB)`

## 5. C10 drawn-structure anchors for this folder

```json
[
 {
  "paper": "wht_2026_falcon",
  "lipid_name": "DLin-MC3-DMA",
  "source_figure": "IUPAC name",
  "read_on": "2026-10-29",
  "note": "Written from the registry name (6Z,9Z,28Z,31Z)-heptatriaconta-6,9,28,31-tetraen-19-yl 4-(dimethylamino)butanoate, in a deliberately different atom order from the corpus string, so the comparison tests the structure and not a copy. The paper draws no structures; MC3 is a catalogue reagent (MedKoo).",
  "reference_inchikey": "NRLNQCOGCKAESA-KWXKLSQISA-N",
  "enumerated_inchikey": "NRLNQCOGCKAESA-KWXKLSQISA-N",
  "skeleton_match": true,
  "stereo_match": true,
  "status": "exact_match"
 },
 {
  "paper": "wht_2026_falcon",
  "lipid_name": "SM-102",
  "source_figure": "IUPAC name",
  "read_on": "2026-10-29",
  "note": "Written from the registry name heptadecan-9-yl 8-[(2-hydroxyethyl)(6-oxo-6-undecoxyhexyl)amino]octanoate, atom order chosen independently of the corpus string. Catalogue reagent (BroadPharm BP-25499).",
  "reference_inchikey": "BGNVBNJYBVCBJH-UHFFFAOYSA-N",
  "enumerated_inchikey": "BGNVBNJYBVCBJH-UHFFFAOYSA-N",
  "skeleton_match": true,
  "stereo_match": true,
  "status": "exact_match"
 }
]
```

> An anchor can only be trusted if it was assembled by a DIFFERENT code path than the enumerator. Check that, do not assume it.

## 6. C11 named commercial reagents appearing in this slice

```json
{
 "DSPC (HL_name)": 1307,
 "Cholesterol (CHL_name)": 2544,
 "DLin-MC3-DMA (IL_name)": 675,
 "18PG (HL_name)": 1114,
 "SM-102 (IL_name)": 1869,
 "DOTAP (HL_name)": 123
}
```

> C11 checks these by InChIKey. It CANNOT check a name that is not in the registry, and it cannot tell you whether the paper's own shorthand means the lipid the registry assumes. `18PG` is the precedent: it was read as symmetric 18:0/18:0 across 360 rows and is in fact 18:0/18:1.

## 7. Folders that are NOT yours

- `curation/ng_2025/` — IN FLIGHT in another session as of 2026-10-30. Not merged. DO NOT READ, DO NOT REPORT ON, DO NOT TOUCH.
- `curation/sd_2026/` — queued, not started - PDFs only
- `curation/jl_2026/` — PAUSED pending the authors' reply about source data; the SI heatmaps are 192 dpi and the printed digits are 5.1x3.5 px, which we established is not extractable
- `curation/liu_iphos/` — blind-validation exercise against the existing SL_2021 rows, not a merge - reviewed under corpus target CORPUS-C
- `curation/amino_lipid_fg_2026/` — blind-validation exercise against existing FG_2026 rows, not a merge - reviewed under CORPUS-C

## 8. Error-class registry (verbatim)

The second table, *Known, not yet automatable*, is where your time is worth most: Tier-1 already covers the first table on every build.

# Error-class registry

The ratchet. Every bug any reviewer finds gets a row here and, wherever possible,
a permanent automated check. **This file only grows.** A class with a check can
never silently recur; a class without one is where review budget should go.

The point of keeping it: reviewers are expensive and find a class *once*; an
automated invariant is free and prevents it *forever*. Discovering a class is
worth paying for. Re-discovering one is not.

## Covered by an automated check

| # | Error class | Found | Check | Where |
|---|---|---|---|---|
| 1 | SMILES that does not parse | early | C1 | `audit.py` |
| 2 | Non-canonical SMILES (same molecule, two strings) | early | C1 | `audit.py` |
| 3 | Disconnected fragments / stray counter-ion | early | C1 | `audit.py` |
| 4 | **Radical electron from a deleted substituent** | 2026-10-29, Siegwart | C1 | `audit.py` + hard guard in both `quaternize()` implementations |
| 5 | Descriptor columns stale vs their own SMILES | 2026-10-28 | C8 | `audit.py` |
| 6 | Mass non-monotonic within a homologous series | 2026-10-28 | C2 | `audit.py` |
| 7 | Duplicate SMILES under two names in one experiment | 2026-10-28 | C3 | `audit.py` |
| 8 | Exact duplicate rows | 2026-10-28 | C4 | `audit.py` |
| 9 | Formulation ratios not summing to 100 | 2026-10-28 | C5 | `audit.py` |
| 10 | Non-finite or absurd `Experiment_value` | 2026-10-28 | C6 | `audit.py` |
| 11 | Field inconsistency; blank required fields | 2026-10-28 | C7 | `audit.py` |
| 12 | Lost E/Z on a stereogenic double bond | 2026-10-28 | C9 | `audit.py` |
| 13 | z-score convention drift (ddof, post-hoc row drops) | 2026-10-28 | rebuild | `build_dataset.py` |
| 14 | **Attachment at the wrong atom** (acyl on a side-chain methyl, not C3) | 2026-10-29, Siegwart | **C10** | `structure_check.py` |
| 15 | **Substituent on the wrong nitrogen** (regiochemistry) | 2026-10-28 Lin, 2026-10-29 Siegwart | **C10** | `structure_check.py` |
| 16 | Heatmap geometry off by a row/column (coloured header band) | 2026-10-28, twice | colour-distance gate | each `digitize_*.py`; >~5 means geometry, not colour |
| 17 | Wrong publication identifier (PMID ≠ DOI ≠ link) | 2026-10-29 | manifest build | `build_manifest.py` resolves every PMID/DOI against NCBI |
| 18 | **The PAPER draws a known benchmark wrong** (ALC-0315 head one carbon short) | 2026-10-29, LipidAI | **C10** `expected_mismatch` | `structure_check.py` — the override gets its own status so it stays visible instead of living in a prose note |
| 19 | **Chain length miscounted by eye from a figure** | standing risk, acted on 2026-10-29 | geometric vertex trace | `curation/xd_2026_lipidai/stage02_trace.py` — skeletonise, trace, RDP-simplify, count; never count by eye |
| 18 | **Two reagent names swapped** (MC3 ↔ KC2, identical formula and exact mass) | 2026-10-29, `AC_2025` | **C11** | `audit.py` — named-reagent registry |
| 19 | **Vendor shorthand expanded to the wrong lipid** (`18PG` read as symmetric 18:0/18:0; it is 18:0/18:1 SOPG) | 2026-10-29, `YZ_2022`/`YZ_2024` | **C11** | `audit.py` — named-reagent registry |
| 20 | **Log axis anchored on the panel box instead of its tick marks** (10³ tick sat 47 px above the bottom rule; every value high by up to 1.76×) | 2026-10-29, `ql_2026` Fig 3c | per-curation | reproduce every printed summary statistic, never a subset |
| 21 | **Provenance text asserts one thing, the code does another** (docs said the heatmap value was carried; the code took the conflicting bar-chart value) | 2026-10-29, `ql_2026` | per-curation | a check that inverts the stored value back to its source |
| 22 | **Blocks standardised separately when they should be pooled**, deleting a real between-block difference | 2026-10-29, `ql_2026` | per-curation | assert the between-block gap SURVIVES; per-block mean 0/sd 1 is a tautology |
| 23 | **Colormap floor/ceiling degeneracy left unflagged** (113 of 320 cells in a band where identical RGB spans ~980 RLU) | 2026-10-29, `ql_2026` Fig 2b | per-curation | flag floored and saturated cells symmetrically |
| 24 | **A C10 anchor documented as machine-derived but actually a typed literal** (the note and the README both said the B-a1 reference was read from vector bond coordinates; it was a string in the script, and nothing in any stage read the graph file) | 2026-10-30, `ng_2025` | per-curation | derive the anchor IN the script from the stored graph, so the provenance claim is executed rather than asserted; a claim no code runs is the same failure as copying the anchor from the enumeration |
| 25 | **Per-item drawings contradicting the paper's own reaction scheme** (about a third of the 104 dialkyl-ketone lipids are drawn as the open-chain mono-adduct, and every one of those also gives the amine nitrogen an extra methyl it cannot have) | 2026-10-30, `ng_2025` | per-curation | classify EVERY drawn member against the scheme and report the count; a spot-check gave 5 where the true figure was 15, understating the ambiguity six-fold |
| 26 | **A reconstructed rule validated on a subset but reported library-wide** (the >5x-MC3 hit threshold reproduces the 11 in-vivo lipids among the 64 phenols, but over all 310 it also catches a thiol lipid that was never taken into mice) | 2026-10-30, `ng_2025` | per-curation | run the reconstructed rule over the FULL set, and state the scope in the claim itself |
| 27 | **A validation calibrated on the quantity it validates** (`digitize_fig2fgh.py` fitted its pixels-per-percent scale by least squares onto the Fig 2e hit rates the panel was checking; the measured tick spacing was computed and never used, and the file had no `assert`) | 2026-10-30, `LL_2026` | **none possible in `audit.py`** — it is a property of a curation's validation code, not of the corpus | review rule: the pixel-to-unit scale must derive ONLY from axis furniture (ticks, labelled gridlines, printed anchors), never from the comparison target. Signature to grep for: a fit whose x-vector is the expected answer |
| 28 | **A measurement failure disguised as a measurement** (`else 0.0` fallback recorded a bar that exists as 0.0 because a `min_area` floor dropped its 450 px² component — and that fabricated zero became the set's single largest "disagreement") | 2026-10-30, `LL_2026` Fig 2h EP6 | per-curation | never let a detector miss fall through to a numeric literal: record `None` and fail, so an absent measurement cannot be averaged, deltaed or reported as agreement |

| 29 | **A hard-coded confidence literal discarding the per-item assessment** (`build_new_curation.py` read `lipid_name -> IL_SMILES` from the enumeration and then wrote `"high"` for every lipid it found, throwing away the `confidence` and `quaternization_site_assumed` columns beside it — so the 24 lipids whose record reads *"ASSUMED: no ring N and the paper draws no structure"* were labelled the project's highest confidence) | 2026-10-30, `ZT_2026` | **C13** | `audit.py` — plus: never write a confidence a per-item record already states; carry it. Inverting a disclosed assumption into an assurance is worse than publishing nothing, because the only surviving provenance then asserts the opposite of the truth |
| 30 | **A component looked up in the live corpus BY NAME, resolved by file order** (`comp.setdefault(name, smiles)` = first row wins, over a corpus written in alphabetical experiment order; 292 names map to more than one structure. `DLin-MC3-DMA` was *ambiguous* — 2,190 rows of the correct ester plus 203 dioxolane rows from `AC_2025` — so the lookup silently returned whichever experiment sorted first, and `MC3-CB1/2/3` were built on KC2) | 2026-10-30, `ZT_2026` | per-curation | the lookup must RAISE on an ambiguous name instead of picking one (`load_corpus_components` now returns the ambiguity set and `smi()` refuses). Regenerating after a `NAME_PATCHES` change fixes the symptom by luck of row order, not the cause; nothing tracks which artefacts were built from a corpus lookup |
| 31 | **Fragment assembly silently inverting a chiral tag** (`RWMol.RemoveAtom` deleted cholesterol's hydroxyl oxygen out of the middle of C3's neighbour list without correcting C3's tag, so all 537 `ZT_2026` rows stored epicholesterol, 3α. The epimers share molecular formula, exact mass to four decimals, the InChIKey SKELETON block and all 17 descriptor columns, so C1–C9 were blind; C10 was blind because its anchor was a byte copy of this code's own output; C11 was blind twice over, comparing skeleton blocks only and inspecting `CHL_*` columns that are empty on every one of these rows) | 2026-10-30, `ZT_2026` | **C14** | `audit.py` — cleave every ester inside `IL_SMILES` corpus-wide and pin the liberated sterol's FULL InChIKey against `reagent_registry.json`. Also: do the condensation the way the chemistry does (the ACID loses its OH, the alcohol's oxygen is retained), so no atom adjacent to a stereocentre is ever deleted |
| 32 | **A bar-top detector fooled by replicate markers** (replicate dots straddling the top edge occlude the MIDDLE of the bar width, so a "first row where >60% of the width is fill" rule fires several pixels INSIDE the bar; the bias is one-sided and non-constant, so it does not cancel in a z-score) | 2026-10-30, independently in `QL_2026` Fig 2c/3c and `XD_2026` Fig 5C **on the same day** | per-curation | measure from the FLANK columns, which a centred marker cannot reach, and require the left and right flank to agree to ≤1–2 px. Note the limit: every printed hit rate in the `QL_2026` panels still reproduces on the corrected values, so the paper's own statistics cannot discriminate this one — it has to be caught geometrically |
| 33 | **A check whose findings are printed but never gate** (`audit.py` discarded `main()`'s return value; C1–C9 printed and the exit code came only from the C10/C11 delegates, so eight simultaneous seeded violations exited 0. Of 43 mutations, 30 appeared in the report and 5 failed a command. `build_dataset.py` never imported `sys` at all) | 2026-10-30, corpus-wide | **the gate itself** | `audit_baseline.json` + `test_checks.py` — each metric has an accepted value WITH its reason, departure in EITHER direction fails (a vanished known condition is how a check that stopped looking announces itself), and a seeded-mutation suite proves the gate fires. The obvious fix is wrong on its own: failing on anything printed fails forever on conditions we accept on purpose, which is how it became print-only |
| 34 | **A quality flag produced, documented, and then dropped at the schema boundary** (1,137 rows carried a convention rather than a measurement; `CLAUDE.md` claimed the project flagged them "in the data (`positional_ambiguity`, `*_assumed`, `confidence`)" and none of those columns existed anywhere. The loss was one stage EARLIER than the merge — no file `build_dataset.py` reads had a flag column — so no change to the builder could have preserved them) | 2026-10-30, four of five curations | **C13** | `schema_quality.json` as the single source of the column names and their vocabulary; `audit.py` enforces the vocabulary and counts rows flagged without a note. A flag that lives only in a per-folder QC JSON is not a flag, it is a memo |

| 35 | **A mask change silently invalidating every window tuned against the old mask** (switching Fig 3c's bar detector from a blue-hue mask to a luminance mask made the panel box's own top rule visible at y = 2301-2314, and the search window started at 2311 — inside it — so every bar would have returned the same wrong top) | 2026-10-30, `QL_2026` Fig 3c | per-curation | assert the first scanned row is not already a full-width rule. Generally: a mask defines what counts as ink, so every threshold, window and offset calibrated under the old one has to be re-derived, not inherited |
| 36 | **A calibration measured mostly from our own detector error** (the `QL_2026` cross-library offset was computed from four benchmark LNPs plotted in both panels — and three of the four were dot-contaminated, so the 1.1116 "offset between the two screens" was largely the bar-top bug measuring itself. Correcting the detector did not sharpen the offset, it removed it: 1.0071 ± 0.0247, not distinguishable from 1.0) | 2026-10-30, `QL_2026` | per-curation | a quantity derived from a small set of shared anchors inherits every bias those anchors carry; check the anchors for the defect under investigation BEFORE trusting a constant derived from them |
| 37 | **A check made blind to a class by the normalization it performs before testing** (915 `IL_name` cells carried a trailing space from upstream. C7's required-field check called `.strip()` on each value *before* asking whether it was blank, so a padded value passed as present and the check could not see padding at all. `"A1I5R22C8 "` and `"A1I5R22C8"` are distinct dict keys carrying the *same* SMILES, so this is a silent join hazard, not a cosmetic one — and the condition was DISCLOSED in `CLAUDE.md` and worked around with two defensive `.strip()` calls in `build_dataset.py` for the corpus's whole life) | 2026-11-07, `BL_2023` + `ZH_2023` | **C7 `padded_cells`** | strip at the single emit choke point, and test for padding on EVERY cell rather than normalizing it away inside the test. Generally: a check that canonicalises its input before comparing has deliberately given up the ability to detect non-canonical input, so whatever it normalises away needs its own check. A known-and-worked-around defect is worse than an unknown one, because the workarounds are evidence everybody saw it and nobody owned it |
| 38 | **A tracked artifact whose input lives in a session scratchpad** (`build_manifest.py` read its NCBI cross-check cache from `.../claude-501/<old-session-uuid>/scratchpad/esum.json`. Scratchpads are session-scoped and get cleared, so the script raised `FileNotFoundError` and `data/experiment_manifest.csv` — committed to git — could not be regenerated from anything in the repository. The cache was also silently STALE: it predated the `QL_2026` and `XD_2026` merges, so the manifest recorded "PMID … not checked against NCBI" for both) | 2026-11-07, second instance — `lit_review/build_data.py` had it on 2026-10-30 | the rule, plus `fetch_esummary.py` | **no script in this repo may read an input from a session scratchpad**; anything a tracked output depends on is written into the repo together with the script that fetches it. Derive paths from `__file__`, not from an absolute path under `/Users/`. Two instances make this a missing rule rather than an oversight |
| 39 | **A number asserted in the project's own documentation that nothing measures** (`CLAUDE.md` stated the F-001 flag count as 1,152 and told the reader to re-measure it with `C13.rows_carrying_a_flag`, which returns 1,652 because it counts `confidence = high`; the figure matching the prose is 1,157. It also gave the chance baseline as `1/31` vitro — the S from before a July re-trim, matching neither the live corpus (42) nor the frozen run (34) — and described a trailing-space defect as affecting "several" rows when it was 915 cells. The count had moved four times, ≥266 → 1,137 → 1,152 → 1,157, each correction by hand) | 2026-11-07, corpus-wide (F-022/F-024/F-026) | **C15** | `provenance/check_docs.py` + `docs_claims.json` — every documented number paired with an expression that measures it, each tagged with the section that claims it. Deliberately does NOT parse the Markdown: the expected value is transcribed by hand, so changing a doc number means changing both and the pair either agrees or the build fails. **A quantity with a prose definition and no implementation gets re-implemented, slightly differently, by every pass that touches it.** Narrower cousin of class 21, which is a curation's provenance text contradicting its own code |
| 40 | **A source stating the same quantity in two parameterisations where only ONE is row-aligned** (`WT_2026`'s URS file gives each formulation as four design parameters AND four molar percents. On that one file of seven they describe different formulations: 0 of 129 design rows agree with their own percent row, while 126 derive exactly one OTHER row's percents, injectively. The design block's `PEG_PEG+Chol` is strictly monotone over all 129 rows with zero violations -- impossible under uniform random sampling, and the signature of a spreadsheet sort applied to a contiguous column selection. `IL_to_nucleicacid_chargeratio` was taken from `NP_ratio`, which travels with the sorted block, so it is wrong on 129 corpus rows: median error 6.86 against a file range of 4.01-25.54. C4, C5 and C6 are all blind because every row's percents still sum to exactly 100) | 2026-11-07, `WT_2026` | per-curation, automatable | where a source gives both parameterisations, DERIVE one from the other and require **row-wise** agreement. The positive signature of the defect is a cross-row injective match with ZERO self-matches. Note what does not work: a range check passes identically on misassigned values because a permutation preserves the multiset, which is why this curation's `validate.py` recorded `PASS N/P range 3.33-37.12` on 129 wrong rows |
| 41 | **A derived "within-model" column silently pooling normalisation units the source says are incomparable** (`Experiment_value_within_model` re-standardises inside (`Experiment_ID` x `Model_type`), but `WT_2026` is normalised per SOURCE FILE and one `Model_type` spans several files -- `Ramos` alone covers 959 rows across five. So the derived column straddles five scales the authors' own `dataset_info.txt` forbids comparing, and `value_normalization_unit` does not describe it: that column declares the unit of `Experiment_value` only. A consumer reaching for the "within model" column to remove an assay offset gets exactly the cross-file mixing the `source_file` flag exists to prevent) | 2026-11-07, `WT_2026` | **open** | a derived column needs its own declared unit, or must refuse to compute where the base unit is not the one it assumes. Generally: adding a re-normalised convenience column creates a SECOND normalisation contract, and a flag that describes only the first one is not a flag on the second |

## A note on the numbering above

Classes **18** and **19** each appear twice: the figure-drawing/vertex-trace pair
and the reagent-name/vendor-shorthand pair were added on the same day by
different passes. **Still duplicated, deliberately.** Round 1's reports are
written against these numbers, so renumbering now would make every reference in
`review/round_2026-10-30/` point at the wrong row. Renumber when round 2 opens,
in one pass, and retire the old numbers rather than reusing them.

**This file is edited concurrently by more than one session**, which is how the
duplication happened and how classes 27-28 below were first written as 24-25,
colliding with 24-26 added from another session minutes earlier. Before
appending, re-read the table and take the next free number; do not assume the
highest number you last saw is still the highest.

## Known, not yet automatable

These need eyes on a figure or judgement about a paper's text. This is where a
reviewer's time is actually worth spending.

| Error class | Why no check | Mitigation |
|---|---|---|
| Name token misread (`C24` = two C4 tails, not 24 carbons) | needs the library's own naming convention | C2 restricted to conventions verified against source metadata |
| Value misread from a figure | needs the figure | overlay PNG per extraction; spot-check ≥10 cells |
| Experimental metadata contradicting itself between caption and Methods | needs both read in full | record BOTH and flag UNRESOLVED (e.g. Siegwart 20 ng/n=3 vs 15 ng/n=5) |
| Commercial reagent that is an isomer mixture | no single true structure | flag `confidence = medium` and say so |
| Regiochemistry a paper never pins | genuinely unknowable | flag `positional_ambiguity` / `*_assumed`, confidence `low` |
| **A C10 anchor copied from the enumeration** (can only pass) | once both sides are canonicalised a copied reference is byte-identical to an honestly re-derived one, so a string comparison would fail the curations that did it right | assemble the reference by a different code path — fragment assembly, or just a different atom ordering — and say which in `read_by` / `note` |
| **De-duplicating reference strings on an assumed identifier format** | LipidAI appended a per-row counter to the DOI *sometimes*, and to a PMCID or a page range otherwise; cutting only at `doi:` over-counts badly. The class is real; the NUMBERS first recorded here were not. "299 against a true 47" was a reviewer claim propagated unchecked; the folder's own `06_lipiddb_qc.json` records 56/50 and the paper states 47, and nine defensible de-duplication rules give 50-116. | cut at the first identifier of ANY kind, cross-check against the count the paper states itself, and **report the rule with the number** - a de-duplicated count is meaningless without it |
| **A structure with no external anchor at all** | commercial lipids supplied without name, formula, HRMS or synthesis (LipidAI's E1–E9) | emit with `confidence = medium`, record an empty `reference_SMILES` row so the absence is visible, and never upgrade without an external structure |

## What the composition checks cannot see

Worth stating plainly, because it is the reason three structural bugs survived
every check that existed at the time:

- **Mass and molecular formula are blind to connectivity.** HRMS matched 4/4 to
  ±0.0001 Da on a molecule with the acyl group on the wrong end of the steroid.
- Deleting a hydroxyl left the **oxygen count unchanged**; only hydrogen moved,
  by one.
- Moving an arm between nitrogens changes **nothing** about the formula.

C10 exists specifically to close this gap, and it is the only Tier-1 check that
compares against an external ground truth rather than internal consistency.

**2026-10-29 added a second instance of exactly this.** `AC_2025` held the
DLin-KC2-DMA dioxolane under the name `DLin-MC3-DMA` and the DLin-MC3-DMA ester
under the name `KC2`. Both are C43H79NO2 at **641.6111 Da** — the same formula
and the same exact mass to four decimals. Every composition check passed over
419 rows for as long as they were in the corpus.

It also taught a second lesson, the harder one. Finding that a name and a
structure disagree does **not** tell you which of the two is wrong, and the two
repairs are not equivalent: rewriting the structure keeps the value with the
name, rewriting the name keeps the value with the structure, and they assign
every measurement to opposite molecules. The first fix here rewrote the
structures and was wrong. Only the upstream source settled it — COMET's master
table stores SMILES and labels with **no names at all**, so the values belong to
the structures and the names were a downstream annotation (built, per its own
comment, by matching molecular weight, which cannot tell these two apart).

**So: when a name and a structure disagree, go to the source and find out which
one the measurement was recorded against. Do not repair the cheaper half.** An
internal-consistency check cannot answer this — it is satisfied by either
repair, including the wrong one.

## What a figure's own printed summaries are for

The `ql_2026` Fig 3c axis error (class 20) is the cleanest case yet of a bug
that **no internal check could see** and that the paper itself gave away three
times over. The values were self-consistent, the bars were detected correctly,
the z-scores were perfect — and every number was wrong by up to 1.76×.

What caught it: the panel's own printed hit rates, and a sentence in the
Results. What let it survive a first pass: `validate.py` asserted two of the
three printed rates and **excused the third** as "a boundary case" — and the
excused one was the only one that discriminated.

The rule this earns: **when a paper prints a statistic derived from a figure you
are digitizing, reproduce ALL of it, and never excuse the part that disagrees.**
A disagreeing summary statistic is the single most informative thing available,
because it is computed from the data you are trying to recover and is
independent of how you recover it. Excusing one is how a calibration error gets
through a passing test suite.

C10 could not have caught it either: C10 only looks at lipids a *curation folder*
enumerated, and these are catalogue reagents that arrive named. That gap is what
C11 closes — a registry of named commercial reagents whose structure is fixed and
public, checked by InChIKey skeleton wherever that name appears in the corpus.
A name like `DSPC` or `SM-102` means one molecule; the corpus should be able to
prove it means that molecule everywhere, and a shorthand like `18PG` should be
expanded once, in one place, with the paper that defines it cited next to it.

