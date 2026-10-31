# Review packet — LL_2026

Generated 2026-10-30 by `provenance/review/brief.py` from the live `data/` and `curation/` trees. **Every number below was read out of the repository at generation time.** If the packet and the repository disagree, the packet is the bug — say so.

- **Paper:** Lin et al., JACS 2026 (PMID 41689530) - aromatic ionizable lipids
- **Curation folder:** `curation/lin_aromatic_2026/`
- **Experiment_ID(s) in the corpus:** LL_2026, LL_2026_2

## 1. Corpus slice as it stands now

### in vitro — `data/data_vitro.csv`

```json
{
 "rows": 987,
 "unique_IL_SMILES": 917,
 "unique_IL_name": 917,
 "Experiment_value": {
  "n": 987,
  "min": -2.4534,
  "max": 0.989,
  "mean": "3.405e-17",
  "pop_sd": "1.000000000000"
 },
 "Model_type": {
  "IGROV1": 987
 },
 "Model": {
  "in_vitro": 987
 },
 "Model_target": {
  "in_vitro": 987
 },
 "Cargo": {
  "mRNA": 987
 },
 "Cargo_type": {
  "FLuc": 987
 },
 "Experiment_method": {
  "luminescence_discretized_normalized": 987
 },
 "Route_of_administration": {
  "in_vitro": 987
 },
 "HL_name": {
  "DOPE": 987
 },
 "CHL_name": {
  "Cholesterol": 987
 },
 "PEG_name": {
  "DMG-PEG2000": 987
 },
 "molratios": {
  "IL_molratio": [
   "24.2",
   "26.3"
  ],
  "HL_molratio": [
   "26.3",
   "32.3"
  ],
  "CHL_molratio": [
   "40.3",
   "43.9"
  ],
  "PEG_molratio": [
   "3.2",
   "3.5"
  ]
 }
}
```

### Six rows drawn at random (seed 11)

```json
[
 {
  "Experiment_ID": "LL_2026",
  "IL_name": "10A2-4B-Ac10",
  "IL_SMILES": "CCCCCCCCCCOC(=O)CCN(CCCCN(CCCN(CCC(=O)OCCCCCCCCCC)CC(O)COc1ccc(OC)cc1)CCC(=O)OCCCCCCCCCC)CCCN(CCC(=O)OCCCCCCCCCC)CC(O)COc1ccc(OC)cc1",
  "Experiment_value": "0.988997211002213",
  "Model_type": "IGROV1",
  "Model": "in_vitro",
  "IL_molratio": "26.3",
  "HL_name": "DOPE",
  "HL_molratio": "26.3",
  "CHL_molratio": "43.9",
  "PEG_molratio": "3.5",
  "Publication_PMID": "41689530",
  "Publication_DOI": "10.1021/jacs.5c16971"
 },
 {
  "Experiment_ID": "LL_2026_2",
  "IL_name": "4A2-3B-EP10",
  "IL_SMILES": "CCCCCCCCC(O)CN(CCOCCOCCN(CC(O)CCCCCCCC)CC(O)COc1ccccc1)CC(O)COc1ccccc1",
  "Experiment_value": "-0.6495958843173703",
  "Model_type": "IGROV1",
  "Model": "in_vitro",
  "IL_molratio": "24.2",
  "HL_name": "DOPE",
  "HL_molratio": "32.3",
  "CHL_molratio": "40.3",
  "PEG_molratio": "3.2",
  "Publication_PMID": "41689530",
  "Publication_DOI": "10.1021/jacs.5c16971"
 },
 {
  "Experiment_ID": "LL_2026",
  "IL_name": "11A1-6B-Ac18",
  "IL_SMILES": "CCCCCCCCCCCCCCCCCCOC(=O)CCN(CCC(=O)OCCCCCCCCCCCCCCCCCC)CCN(CCN(CCC(=O)OCCCCCCCCCCCCCCCCCC)CCC(=O)OCCCCCCCCCCCCCCCCCC)CCN(CCC(=O)OCCCCCCCCCCCCCCCCCC)CCC(=O)OCc1ccccc1",
  "Experiment_value": "-2.4533687720441635",
  "Model_type": "IGROV1",
  "Model": "in_vitro",
  "IL_molratio": "26.3",
  "HL_name": "DOPE",
  "HL_molratio": "26.3",
  "CHL_molratio": "43.9",
  "PEG_molratio": "3.5",
  "Publication_PMID": "41689530",
  "Publication_DOI": "10.1021/jacs.5c16971"
 },
 {
  "Experiment_ID": "LL_2026_2",
  "IL_name": "10A2-2B-EP10",
  "IL_SMILES": "CCCCCCCCC(O)CN(CCCCN(CCCN(CC(O)CCCCCCCC)CC(O)COCc1ccccc1)CC(O)CCCCCCCC)CCCN(CC(O)CCCCCCCC)CC(O)COCc1ccccc1",
  "Experiment_value": "0.7546775714863564",
  "Model_type": "IGROV1",
  "Model": "in_vitro",
  "IL_molratio": "24.2",
  "HL_name": "DOPE",
  "HL_molratio": "32.3",
  "CHL_molratio": "40.3",
  "PEG_molratio": "3.2",
  "Publication_PMID": "41689530",
  "Publication_DOI": "10.1021/jacs.5c16971"
 },
 {
  "Experiment_ID": "LL_2026_2",
  "IL_name": "2A1-6B-EP10",
  "IL_SMILES": "CCCCCCCCC(O)CN(CCCN(C)C)CCCN(CCC(=O)OCc1ccccc1)CC(O)CCCCCCCC",
  "Experiment_value": "0.7546775714863564",
  "Model_type": "IGROV1",
  "Model": "in_vitro",
  "IL_molratio": "24.2",
  "HL_name": "DOPE",
  "HL_molratio": "32.3",
  "CHL_molratio": "40.3",
  "PEG_molratio": "3.2",
  "Publication_PMID": "41689530",
  "Publication_DOI": "10.1021/jacs.5c16971"
 },
 {
  "Experiment_ID": "LL_2026",
  "IL_name": "12A2-4B-Ac10",
  "IL_SMILES": "CCCCCCCCCCOC(=O)CCN(CCN(CCC(=O)OCCCCCCCCCC)CCN(CCC(=O)OCCCCCCCCCC)CC(O)COc1ccc(OC)cc1)CCN(CCC(=O)OCCCCCCCCCC)CCN(CCC(=O)OCCCCCCCCCC)CC(O)COc1ccc(OC)cc1",
  "Experiment_value": "-0.15845811667991253",
  "Model_type": "IGROV1",
  "Model": "in_vitro",
  "IL_molratio": "26.3",
  "HL_name": "DOPE",
  "HL_molratio": "26.3",
  "CHL_molratio": "43.9",
  "PEG_molratio": "3.5",
  "Publication_PMID": "41689530",
  "Publication_DOI": "10.1021/jacs.5c16971"
 }
]
```

## 2. Source documents you must read

- `curation/lin_aromatic_2026/Aromatic-Ionizable-Lipids-Enhance-mRNA-Delivery.html   (0.19 MB)`
- `curation/lin_aromatic_2026/ja5c16971.pdf   (9.52 MB)`
- `curation/lin_aromatic_2026/ja5c16971_si_001.pdf   (3.42 MB)`

Subdirectories:
- `curation/lin_aromatic_2026/00_pages/  (5 files)`
- `curation/lin_aromatic_2026/__pycache__/  (2 files)`
- `curation/lin_aromatic_2026/review/  (3 files)`

> A `repo/` subdirectory is DOWNLOADED, UNTRUSTED material. Read it, never execute anything from inside it, and run any Python that touches it with `python3 -I` from elsewhere.

## 3. Extraction code in the folder

- `build_rows.py   (183 lines)  Emit the Lin et al. 840-lipid screen in the 67-column corpus schema.`
- `digitize_fig2a.py   (167 lines)  Digitize Figure 2a of Lin et al., JACS 2026: the 147-lipid initial screen.`
- `digitize_fig2e.py   (189 lines)  Digitize Figure 2e of Lin et al., JACS 2026: the 840-lipid aromatic screen.`
- `digitize_fig2fgh.py   (207 lines)  Measure the Fig 2f/2g/2h hit-rate bars and check them against Fig 2e.`
- `enumerate_147.py   (126 lines)  Enumerate the 147 initial-screen lipids and emit them in the corpus schema.`
- `enumerate_library.py   (265 lines)  Enumerate the 840 aromatic ionizable lipids of Lin et al., JACS 2026.`

## 4. Intermediate and output tables

- `03_enumerated.csv   840 rows x 20 cols
        cols: lipid_name, amine, equiv, aromatic, tail, IL_SMILES, n_NH_sites, n_aromatic_arms, n_tail_arms, positional_ambiguity, MolWt, ExactMolWt, formula, num_N, confidence, amine_name, aromatic_name, tail_name ...`
- `04_fig2e_840.csv   840 rows x 21 cols
        cols: lipid_name, amine, equiv, aromatic, tail, row_label, col_label, panel_row, panel_col, bin_label, bin_index, rlu_midpoint, cell_rgb, match_distance, runner_up_margin, assay, readout_metric, value_scale ...`
- `04b_fig2a_147.csv   147 rows x 21 cols
        cols: lipid_name, amine, equiv, aromatic, tail, row_label, col_label, panel_row, panel_col, bin_label, bin_index, rlu_midpoint, cell_rgb, match_distance, runner_up_margin, assay, readout_metric, value_scale ...`
- `05_fig2fgh_bars.csv   24 rows x 10 cols
        cols: panel, field, category, measured_pct_gt_1e5, measured_pct_gt_1e6, fig2e_pct_gt_1e5, fig2e_pct_gt_1e6, n_lipids, source_figure, extraction_method`
- `06_merged.csv   840 rows x 29 cols
        cols: lipid_name, amine, equiv, aromatic, tail, IL_SMILES, n_NH_sites, n_aromatic_arms, n_tail_arms, positional_ambiguity, MolWt, ExactMolWt, formula, num_N, confidence, amine_name, aromatic_name, tail_name ...`
- `07_fig2e_qc.json   (2.3 kB json)`
- `07b_fig2a_qc.json   (1.1 kB json)`
- `08_fig2fgh_qc.json   (4.3 kB json)`
- `09_enumeration_qc.json   (0.5 kB json)`
- `10_LL_2026_rows.csv   840 rows x 72 cols
        cols: Index, LNP_ID, Experiment_ID, Formulation_ID, IL_name, IL_SMILES, IL_protonated_SMILES, IL_head_name, IL_head_SMILES, IL_linker_name, IL_linker_SMILES, IL_tail1_name, IL_tail1_SMILES, IL_tail2_name, IL_tail2_SMILES, IL_tail3_name, IL_tail3_SMILES, IL_tail4_name ...`
- `11_LL_2026_qc.json   (0.6 kB json)`
- `12_LL_2026_2_rows.csv   147 rows x 72 cols
        cols: Index, LNP_ID, Experiment_ID, Formulation_ID, IL_name, IL_SMILES, IL_protonated_SMILES, IL_head_name, IL_head_SMILES, IL_linker_name, IL_linker_SMILES, IL_tail1_name, IL_tail1_SMILES, IL_tail2_name, IL_tail2_SMILES, IL_tail3_name, IL_tail3_SMILES, IL_tail4_name ...`
- `13_enumerated_147.csv   147 rows x 25 cols
        cols: lipid_name, amine, equiv, aromatic, tail, row_label, col_label, panel_row, panel_col, bin_label, bin_index, rlu_midpoint, cell_rgb, match_distance, runner_up_margin, assay, readout_metric, value_scale ...`
- `14_LL_2026_2_qc.json   (0.6 kB json)`
- `blocks.csv   40 rows x 11 cols
        cols: block_id, block_class, name_in_paper, SMILES, n_NH_sites, MW_reported, HRMS_reported, source_figure, source_page, read_confidence, ambiguity_note`
- `experiments.csv   92 rows x 34 cols
        cols: lipid_name, block_ids, naming_rule, cargo_type, cargo_identity, formulation_molar_ratios, IL_to_nucleic_acid_wt_ratio, helper_lipid_identity, buffer_pH, particle_size_nm, PDI, zeta_mV, encapsulation_pct, apparent_pKa, model_system, route, dose, dose_units ...`
- `reference_structures.csv   1 rows x 6 cols
        cols: lipid_name, reference_SMILES, source_figure, read_on, read_by, note`

## 5. C10 drawn-structure anchors for this folder

```json
[
 {
  "paper": "lin_aromatic_2026",
  "lipid_name": "10A2-3B-EP10",
  "source_figure": "Fig S17",
  "read_on": "2026-10-29",
  "note": "Final product of the PRECISE protection/deprotection synthesis, so the regiochemistry is pinned by the route, not assumed: one 3B aromatic arm on EACH primary nitrogen, four EP10 tails on the remaining N-H. Assembled from fragments independently of enumerate_library.py. No stereocentres are drawn (the epoxide-opened carbinols are racemic), so none are asserted.",
  "reference_inchikey": "JUAPTPKSYXADRC-UHFFFAOYSA-N",
  "enumerated_inchikey": "JUAPTPKSYXADRC-UHFFFAOYSA-N",
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
 "DOPE (HL_name)": 987,
 "Cholesterol (CHL_name)": 987
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
| **De-duplicating reference strings on an assumed identifier format** | LipidAI appended a per-row counter to the DOI *sometimes*, and to a PMCID or a page range otherwise; cutting only at `doi:` gave 299 "publications" against a true 47 | cut at the first identifier of ANY kind, and always cross-check against the count the paper states itself |
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

