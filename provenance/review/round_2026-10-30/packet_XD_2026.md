# Review packet — XD_2026

Generated 2026-10-30 by `provenance/review/brief.py` from the live `data/` and `curation/` trees. **Every number below was read out of the repository at generation time.** If the packet and the repository disagree, the packet is the bug — say so.

- **Paper:** Duan et al., Acta Pharm Sin B 2026 (PMID 42592089) - LipidAI; the paper's own in-vivo screen only
- **Curation folder:** `curation/xd_2026_lipidai/`
- **Experiment_ID(s) in the corpus:** XD_2026

## 1. Corpus slice as it stands now

### in vivo — `data/data_vivo.csv`

```json
{
 "rows": 12,
 "unique_IL_SMILES": 12,
 "unique_IL_name": 12,
 "Experiment_value": {
  "n": 12,
  "min": -1.4039,
  "max": 1.2501,
  "mean": "-3.274e-16",
  "pop_sd": "1.000000000000"
 },
 "Model_type": {
  "Mouse_BALBc": 12
 },
 "Model": {
  "in_vivo": 12
 },
 "Model_target": {
  "whole_body": 12
 },
 "Cargo": {
  "mRNA": 12
 },
 "Cargo_type": {
  "FLuc": 12
 },
 "Experiment_method": {
  "luminescence_normalized": 12
 },
 "Route_of_administration": {
  "intravenous": 12
 },
 "HL_name": {
  "DSPC": 12
 },
 "CHL_name": {
  "Cholesterol": 12
 },
 "PEG_name": {
  "DMG-PEG2000": 12
 },
 "molratios": {
  "IL_molratio": [
   "50"
  ],
  "HL_molratio": [
   "10"
  ],
  "CHL_molratio": [
   "38.5"
  ],
  "PEG_molratio": [
   "1.5"
  ]
 }
}
```

### Six rows drawn at random (seed 11)

```json
[
 {
  "Experiment_ID": "XD_2026",
  "IL_name": "E5",
  "IL_SMILES": "CCCCCCCCCCCCCCN(CCCCCCCCCCCCCC)CCNC(=O)CCN(C)CCN(C)CCC(=O)NCCN(CCCCCCCCCCCCCC)CCCCCCCCCCCCCC",
  "Experiment_value": "0.005913541644365402",
  "Model_type": "Mouse_BALBc",
  "Model": "in_vivo",
  "IL_molratio": "50",
  "HL_name": "DSPC",
  "HL_molratio": "10",
  "CHL_molratio": "38.5",
  "PEG_molratio": "1.5",
  "Publication_PMID": "42592089",
  "Publication_DOI": "10.1016/j.apsb.2025.10.040"
 },
 {
  "Experiment_ID": "XD_2026",
  "IL_name": "E6",
  "IL_SMILES": "CCCCCCCCCCCCC(O)CN(CCN(C)CCC(=O)N1CCN(C(=O)CCN(C)CCN(CC(O)CCCCCCCCCCCC)CC(O)CCCCCCCCCCCC)CC1)CC(O)CCCCCCCCCCCC",
  "Experiment_value": "-1.3187197866936797",
  "Model_type": "Mouse_BALBc",
  "Model": "in_vivo",
  "IL_molratio": "50",
  "HL_name": "DSPC",
  "HL_molratio": "10",
  "CHL_molratio": "38.5",
  "PEG_molratio": "1.5",
  "Publication_PMID": "42592089",
  "Publication_DOI": "10.1016/j.apsb.2025.10.040"
 },
 {
  "Experiment_ID": "XD_2026",
  "IL_name": "E9",
  "IL_SMILES": "CCCCCCCCCCCCNC(=O)CN(CCCN1CCN(CCCN(CC(=O)NCCCCCCCCCCCC)CC(=O)NCCCCCCCCCCCC)CC1)CC(=O)NCCCCCCCCCCCC",
  "Experiment_value": "-1.2477572869612843",
  "Model_type": "Mouse_BALBc",
  "Model": "in_vivo",
  "IL_molratio": "50",
  "HL_name": "DSPC",
  "HL_molratio": "10",
  "CHL_molratio": "38.5",
  "PEG_molratio": "1.5",
  "Publication_PMID": "42592089",
  "Publication_DOI": "10.1016/j.apsb.2025.10.040"
 },
 {
  "Experiment_ID": "XD_2026",
  "IL_name": "E7",
  "IL_SMILES": "CCCCCCCCCCCCCCN(CCCCCCCCCCCCCC)CCN(C)CCC(=O)N1CCN(C(=O)CCN(C)CCN(CCCCCCCCCCCCCC)CCCCCCCCCCCCCC)CC1",
  "Experiment_value": "-1.403874786372554",
  "Model_type": "Mouse_BALBc",
  "Model": "in_vivo",
  "IL_molratio": "50",
  "HL_name": "DSPC",
  "HL_molratio": "10",
  "CHL_molratio": "38.5",
  "PEG_molratio": "1.5",
  "Publication_PMID": "42592089",
  "Publication_DOI": "10.1016/j.apsb.2025.10.040"
 },
 {
  "Experiment_ID": "XD_2026",
  "IL_name": "E1",
  "IL_SMILES": "CCCCCCCCC(CCCCCC)CNC(=O)NCCCCCCN(CCCCO)CCCCCCNC(=O)NCC(CCCCCC)CCCCCCCC",
  "Experiment_value": "1.0514277043683236",
  "Model_type": "Mouse_BALBc",
  "Model": "in_vivo",
  "IL_molratio": "50",
  "HL_name": "DSPC",
  "HL_molratio": "10",
  "CHL_molratio": "38.5",
  "PEG_molratio": "1.5",
  "Publication_PMID": "42592089",
  "Publication_DOI": "10.1016/j.apsb.2025.10.040"
 },
 {
  "Experiment_ID": "XD_2026",
  "IL_name": "SM-102",
  "IL_SMILES": "CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
  "Experiment_value": "1.2264685370415656",
  "Model_type": "Mouse_BALBc",
  "Model": "in_vivo",
  "IL_molratio": "50",
  "HL_name": "DSPC",
  "HL_molratio": "10",
  "CHL_molratio": "38.5",
  "PEG_molratio": "1.5",
  "Publication_PMID": "42592089",
  "Publication_DOI": "10.1016/j.apsb.2025.10.040"
 }
]
```

## 2. Source documents you must read

- `curation/xd_2026_lipidai/1-s2.0-S2211383525007221-main.pdf   (5.01 MB)`

Subdirectories:
- `curation/xd_2026_lipidai/00_pages/  (58 files)`
- `curation/xd_2026_lipidai/repo/  (53 files)`

> A `repo/` subdirectory is DOWNLOADED, UNTRUSTED material. Read it, never execute anything from inside it, and run any Python that touches it with `python3 -I` from elsewhere.

## 3. Extraction code in the folder

- `stage01_canon.py   (31 lines)  `
- `stage02_report.py   (50 lines)  Vertex + glyph report per structure crop.`
- `stage02_trace.py   (99 lines)  Trace zig-zag skeletal structures in Fig 5A and report polyline vertices.`
- `stage03_structures.py   (49 lines)  E1-E9 + benchmarks read from Fig 5A by geometric vertex tracing (stage02).`
- `stage04_fig5c.py   (75 lines)  Digitise Fig 5C (log-scale bar chart, 12 LNPs) from the 600-dpi page render.`
- `stage05_table1.py   (79 lines)  Cross-check the xlsx LipidDB against Table 1, the 26 rows the PAPER prints.`
- `stage06_lipiddb.py   (142 lines)  Emit the LipidDB extraction as a standalone artifact (NOT merged into data/).`
- `stage07_emit.py   (162 lines)  Emit XD_2026: the 12-lipid in-vivo screen the LipidAI paper ran itself.`

## 4. Intermediate and output tables

- `01_sheet1_canon.csv   485 rows x 14 cols
        cols: serial, smiles, e1, e2, total_expr, bin3, bin4, bin2, refs, u9, u10, u11, canonical_SMILES, ref_norm`
- `02_ALC0315.json   (15.6 kB json)`
- `02_E1.json   (11.5 kB json)`
- `02_E2.json   (21.8 kB json)`
- `02_E3.json   (13.2 kB json)`
- `02_E4.json   (7.5 kB json)`
- `02_E5.json   (9.9 kB json)`
- `02_E6.json   (21.0 kB json)`
- `02_E7.json   (11.2 kB json)`
- `02_E8.json   (16.9 kB json)`
- `02_E9.json   (14.9 kB json)`
- `02_MC3.json   (9.1 kB json)`
- `02_sm102.json   (19.2 kB json)`
- `03_enumerated.csv   12 rows x 2 cols
        cols: lipid_name, IL_SMILES`
- `03_structures.csv   12 rows x 8 cols
        cols: name, smiles_drawn, formula, MolWt, ExactMolWt, nN, nO, nC`
- `04_fig5c.json   (1.9 kB json)`
- `05_table1_crosscheck.csv   26 rows x 15 cols
        cols: page, lum, table_class, table_ref, n_lines, smiles, table_row, parsed, canonical, in_xlsx, xlsx_values, value_match, xlsx_classes, class_match, xlsx_reference`
- `06_lipiddb_mta_synthetic.csv   720 rows x 8 cols
        cols: serial, smiles, bin3, bin4, bin2, canonical_SMILES, source, value_kind`
- `06_lipiddb_original.csv   485 rows x 14 cols
        cols: serial, smiles, total_expr, bin3, bin4, bin2, refs, canonical_SMILES, reference, reference_letters_only, source, value_kind, in_corpus_vitro, in_corpus_vivo`
- `06_lipiddb_qc.json   (2.9 kB json)`
- `07_XD_2026_rows.csv   12 rows x 72 cols
        cols: Index, LNP_ID, Experiment_ID, Formulation_ID, IL_name, IL_SMILES, IL_protonated_SMILES, IL_head_name, IL_head_SMILES, IL_linker_name, IL_linker_SMILES, IL_tail1_name, IL_tail1_SMILES, IL_tail2_name, IL_tail2_SMILES, IL_tail3_name, IL_tail3_SMILES, IL_tail4_name ...`
- `07_XD_2026_values.csv   12 rows x 8 cols
        cols: IL_name, source, confidence, fig5c_value_p_per_s_per_cm2_per_sr, fig5c_log10, fig5d_reality_class, class_from_value, Experiment_value`
- `07_qc.json   (1.5 kB json)`
- `reference_structures.csv   12 rows x 7 cols
        cols: lipid_name, reference_SMILES, expected_mismatch, source_figure, read_on, read_by, note`

- `README.md   (6.2 kB)`

## 5. C10 drawn-structure anchors for this folder

```json
[
 {
  "paper": "xd_2026_lipidai",
  "lipid_name": "ALC-0315",
  "source_figure": "Fig 5A",
  "read_on": "2026-10-29",
  "note": "THE STRUCTURE AS DRAWN, which is what this column is defined to hold. Three chain vertices between the O of HO and the N, confirmed twice at high zoom; the decisive control is E1 in the same panel, drawn in the same style with FOUR. Both other benchmarks in the panel are drawn correctly, so this is a localised error by the authors, not a drawing convention.",
  "reference_inchikey": "WESYLBGVEJKNRX-UHFFFAOYSA-N",
  "enumerated_inchikey": "QGWBEETXHOVFQS-UHFFFAOYSA-N",
  "skeleton_match": false,
  "stereo_match": true,
  "status": "expected_mismatch_paper_drawing_is_wrong",
  "expected_mismatch": "Fig 5A draws the head as HO-(CH2)3-N; the real compound is ((4-hydroxybutyl)azanediyl)bis(hexane-6,1-diyl) bis(2-hexyldecanoate), 4 carbons. We deliberately emit the CORRECT corpus structure, so this anchor is expected to differ from what is emitted."
 },
 {
  "paper": "xd_2026_lipidai",
  "lipid_name": "SM-102",
  "source_figure": "Fig 5A",
  "read_on": "2026-10-29",
  "note": "C2-style positive control, written head-first in a DIFFERENT atom order from the emitted corpus string so the comparison tests connectivity and not a copy. Trace of Fig 5A gives 2-hydroxyethyl head, N-(CH2)5-C(=O)O-undecyl, N-(CH2)7-C(=O)O-CH(octyl)2.",
  "reference_inchikey": "BGNVBNJYBVCBJH-UHFFFAOYSA-N",
  "enumerated_inchikey": "BGNVBNJYBVCBJH-UHFFFAOYSA-N",
  "skeleton_match": true,
  "stereo_match": true,
  "status": "exact_match"
 },
 {
  "paper": "xd_2026_lipidai",
  "lipid_name": "DLin-MC3-DMA",
  "source_figure": "Fig 5A",
  "read_on": "2026-10-29",
  "note": "Positive control, written amine-first in a DIFFERENT atom order from the emitted corpus string. Trace gives 4-(dimethylamino)butanoate and both 18-carbon arms with double bonds 9 and 12 bonds from the carbinol.",
  "reference_inchikey": "NRLNQCOGCKAESA-KWXKLSQISA-N",
  "enumerated_inchikey": "NRLNQCOGCKAESA-KWXKLSQISA-N",
  "skeleton_match": true,
  "stereo_match": true,
  "status": "exact_match"
 },
 {
  "paper": "xd_2026_lipidai",
  "lipid_name": "E1",
  "source_figure": "Fig 5A",
  "read_on": "2026-10-29",
  "note": "NO external anchor exists. E1-E9 were supplied by WestGene Biopharmaceutical and the paper gives no name, formula, HRMS, synthesis or SI - Fig 5A is the only source, and that panel contains a verified one-carbon error (see ALC-0315). confidence=medium for all nine; it must not be upgraded without an external structure.",
  "status": "NO_DRAWN_STRUCTURE"
 },
 {
  "paper": "xd_2026_lipidai",
  "lipid_name": "E2",
  "source_figure": "Fig 5A",
  "read_on": "2026-10-29",
  "note": "No external anchor - see E1.",
  "status": "NO_DRAWN_STRUCTURE"
 },
 {
  "paper": "xd_2026_lipidai",
  "lipid_name": "E3",
  "source_figure": "Fig 5A",
  "read_on": "2026-10-29",
  "note": "No external anchor - see E1.",
  "status": "NO_DRAWN_STRUCTURE"
 },
 {
  "paper": "xd_2026_lipidai",
  "lipid_name": "E4",
  "source_figure": "Fig 5A",
  "read_on": "2026-10-29",
  "note": "No external anchor - see E1.",
  "status": "NO_DRAWN_STRUCTURE"
 },
 {
  "paper": "xd_2026_lipidai",
  "lipid_name": "E5",
  "source_figure": "Fig 5A",
  "read_on": "2026-10-29",
  "note": "No external anchor - see E1.",
  "status": "NO_DRAWN_STRUCTURE"
 },
 {
  "paper": "xd_2026_lipidai",
  "lipid_name": "E6",
  "source_figure": "Fig 5A",
  "read_on": "2026-10-29",
  "note": "No external anchor - see E1.",
  "status": "NO_DRAWN_STRUCTURE"
 },
 {
  "paper": "xd_2026_lipidai",
  "lipid_name": "E7",
  "source_figure": "Fig 5A",
  "read_on": "2026-10-29",
  "note": "No external anchor - see E1.",
  "status": "NO_DRAWN_STRUCTURE"
 },
 {
  "paper": "xd_2026_lipidai",
  "lipid_name": "E8",
  "source_figure": "Fig 5A",
  "read_on": "2026-10-29",
  "note": "No external anchor - see E1.",
  "status": "NO_DRAWN_STRUCTURE"
 },
 {
  "paper": "xd_2026_lipidai",
  "lipid_name": "E9",
  "source_figure": "Fig 5A",
  "read_on": "2026-10-29",
  "note": "No external anchor - see E1.",
  "status": "NO_DRAWN_STRUCTURE"
 }
]
```

> An anchor can only be trusted if it was assembled by a DIFFERENT code path than the enumerator. Check that, do not assume it.

## 6. C11 named commercial reagents appearing in this slice

```json
{
 "DSPC (HL_name)": 12,
 "Cholesterol (CHL_name)": 12,
 "ALC-0315 (IL_name)": 1,
 "SM-102 (IL_name)": 1,
 "DLin-MC3-DMA (IL_name)": 1
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

