# Review packet — QL_2026

Generated 2026-10-30 by `provenance/review/brief.py` from the live `data/` and `curation/` trees. **Every number below was read out of the repository at generation time.** If the packet and the repository disagree, the packet is the bug — say so.

- **Paper:** Liu et al., Adv Healthc Mater 2026 (PMID 41145968) - bispecific-linker thiosuccinimide lipids
- **Curation folder:** `curation/ql_2026/`
- **Experiment_ID(s) in the corpus:** QL_2026

## 1. Corpus slice as it stands now

### in vitro — `data/data_vitro.csv`

```json
{
 "rows": 344,
 "unique_IL_SMILES": 344,
 "unique_IL_name": 344,
 "Experiment_value": {
  "n": 344,
  "min": -1.2604,
  "max": 3.5243,
  "mean": "2.825e-16",
  "pop_sd": "1.000000000000"
 },
 "Model_type": {
  "HeLa": 344
 },
 "Model": {
  "in_vitro": 344
 },
 "Model_target": {
  "in_vitro": 344
 },
 "Cargo": {
  "mRNA": 344
 },
 "Cargo_type": {
  "FLuc": 344
 },
 "Experiment_method": {
  "luminescence_normalized": 344
 },
 "Route_of_administration": {
  "in_vitro": 344
 },
 "HL_name": {
  "DOPE": 344
 },
 "CHL_name": {
  "Cholesterol": 344
 },
 "PEG_name": {
  "DMG-PEG2000": 344
 },
 "molratios": {
  "IL_molratio": [
   "35.0"
  ],
  "HL_molratio": [
   "16.0"
  ],
  "CHL_molratio": [
   "46.5"
  ],
  "PEG_molratio": [
   "2.5"
  ]
 }
}
```

### Six rows drawn at random (seed 11)

```json
[
 {
  "Experiment_ID": "QL_2026",
  "IL_name": "114-L7-C8",
  "IL_SMILES": "CCCCCCCCSC1CC(=O)N(CCOCCOCCC(=O)NCCNCCNC(=O)CCOCCOCCN2C(=O)CC(SCCCCCCCC)C2=O)C1=O",
  "Experiment_value": "-0.5069080325449712",
  "Model_type": "HeLa",
  "Model": "in_vitro",
  "IL_molratio": "35.0",
  "HL_name": "DOPE",
  "HL_molratio": "16.0",
  "CHL_molratio": "46.5",
  "PEG_molratio": "2.5",
  "Publication_PMID": "41145968",
  "Publication_DOI": "10.1002/adhm.202502225"
 },
 {
  "Experiment_ID": "QL_2026",
  "IL_name": "T3A-L2-C8",
  "IL_SMILES": "CCCCCCCCSSCCC(=O)NCCCN(CCCNC(=O)CCSSCCCCCCCC)CCCNC(=O)CCSSCCCCCCCC",
  "Experiment_value": "1.3571804295481915",
  "Model_type": "HeLa",
  "Model": "in_vitro",
  "IL_molratio": "35.0",
  "HL_name": "DOPE",
  "HL_molratio": "16.0",
  "CHL_molratio": "46.5",
  "PEG_molratio": "2.5",
  "Publication_PMID": "41145968",
  "Publication_DOI": "10.1002/adhm.202502225"
 },
 {
  "Experiment_ID": "QL_2026",
  "IL_name": "114-L8-C12",
  "IL_SMILES": "CCCCCCCCCCCCSC1CC(=O)N(CCOCCOCCOCCC(=O)NCCNCCNC(=O)CCOCCOCCOCCN2C(=O)CC(SCCCCCCCCCCCC)C2=O)C1=O",
  "Experiment_value": "-0.3439928997388166",
  "Model_type": "HeLa",
  "Model": "in_vitro",
  "IL_molratio": "35.0",
  "HL_name": "DOPE",
  "HL_molratio": "16.0",
  "CHL_molratio": "46.5",
  "PEG_molratio": "2.5",
  "Publication_PMID": "41145968",
  "Publication_DOI": "10.1002/adhm.202502225"
 },
 {
  "Experiment_ID": "QL_2026",
  "IL_name": "110-L5-C6",
  "IL_SMILES": "CCCCCCSC1CC(=O)N(CCCCCC(=O)NCCN(CCNC(=O)CCCCCN2C(=O)CC(SCCCCCC)C2=O)CCNC(=O)CCCCCN2C(=O)CC(SCCCCCC)C2=O)C1=O",
  "Experiment_value": "-0.9029842420051316",
  "Model_type": "HeLa",
  "Model": "in_vitro",
  "IL_molratio": "35.0",
  "HL_name": "DOPE",
  "HL_molratio": "16.0",
  "CHL_molratio": "46.5",
  "PEG_molratio": "2.5",
  "Publication_PMID": "41145968",
  "Publication_DOI": "10.1002/adhm.202502225"
 },
 {
  "Experiment_ID": "QL_2026",
  "IL_name": "T3A-L5-C6",
  "IL_SMILES": "CCCCCCSC1CC(=O)N(CCCCCC(=O)NCCCN(CCCNC(=O)CCCCCN2C(=O)CC(SCCCCCC)C2=O)CCCNC(=O)CCCCCN2C(=O)CC(SCCCCCC)C2=O)C1=O",
  "Experiment_value": "-0.3439928997388166",
  "Model_type": "HeLa",
  "Model": "in_vitro",
  "IL_molratio": "35.0",
  "HL_name": "DOPE",
  "HL_molratio": "16.0",
  "CHL_molratio": "46.5",
  "PEG_molratio": "2.5",
  "Publication_PMID": "41145968",
  "Publication_DOI": "10.1002/adhm.202502225"
 },
 {
  "Experiment_ID": "QL_2026",
  "IL_name": "80-L4-C10",
  "IL_SMILES": "CCCCCCCCCCSC1CC(=O)N(CCCC(=O)NCCCN(C)C)C1=O",
  "Experiment_value": "-0.5069080325449712",
  "Model_type": "HeLa",
  "Model": "in_vitro",
  "IL_molratio": "35.0",
  "HL_name": "DOPE",
  "HL_molratio": "16.0",
  "CHL_molratio": "46.5",
  "PEG_molratio": "2.5",
  "Publication_PMID": "41145968",
  "Publication_DOI": "10.1002/adhm.202502225"
 }
]
```

## 2. Source documents you must read

- `curation/ql_2026/Adv Healthcare Materials - 2025 - Liu - Synthesis and Screening of Bispecific Linker‐based Combinatorial Libraries.pdf   (7.91 MB)`
- `curation/ql_2026/adhm70370-sup-0001-suppmat.docx   (6.55 MB)`

Subdirectories:
- `curation/ql_2026/00_pages/  (3 files)`
- `curation/ql_2026/__pycache__/  (1 files)`
- `curation/ql_2026/review/  (4 files)`

> A `repo/` subdirectory is DOWNLOADED, UNTRUSTED material. Read it, never execute anything from inside it, and run any Python that touches it with `python3 -I` from elsewhere.

## 3. Extraction code in the folder

- `build_rows.py   (296 lines)  Emit the 344 Liu et al. lipids in the 72-column corpus schema. NOT merged.`
- `digitize_fig2b.py   (185 lines)  Digitize Figure 2b of Liu et al., Adv Healthc Mater 2026: library 1, 320 lipids.`
- `digitize_fig2c.py   (146 lines)  Digitize Figure 2c: the three-tail lipids above 2e4, on a log axis.`
- `digitize_fig2def.py   (171 lines)  Digitize Figure 2d/2e/2f, the printed hit-rate bars, as an INDEPENDENT check`
- `digitize_fig3c.py   (177 lines)  Digitize Figure 3c: library 2 (24 lipids) plus four benchmark LNPs, log axis.`
- `enumerate_library.py   (217 lines)  Enumerate the 344 thiosuccinimide lipids of Liu et al., Adv Healthc Mater 2026.`
- `validate.py   (188 lines)  Independent check of the QL_2026 curation. Re-reads the OUTPUT files and the`

## 4. Intermediate and output tables

- `01_fig2b_320.csv   320 rows x 13 cols
        cols: lipid_name, amine, bl, tail, panel_row, panel_col, RLU, saturated, floored, cell_rgb, match_distance, source_figure, extraction_method`
- `02_fig2b_qc.json   (0.8 kB json)`
- `03_fig2def_bars.csv   16 rows x 7 cols
        cols: panel, label, value_pct, x0, x1, y_top, y_base`
- `04_fig2def_qc.json   (2.3 kB json)`
- `05_fig3c_24.csv   28 rows x 8 cols
        cols: lipid_name, is_benchmark, RLU, x0, x1, y_top, source_figure, extraction_method`
- `06_fig3c_qc.json   (1.6 kB json)`
- `07_fig2c_bars.csv   32 rows x 8 cols
        cols: lipid_name, is_benchmark, RLU, x0, x1, y_top, source_figure, extraction_method`
- `08_fig2c_qc.json   (0.7 kB json)`
- `09_enumerated_344.csv   344 rows x 11 cols
        cols: lipid_name, library, amine, bl, tail, n_arms, IL_SMILES, formula, MolWt, ExactMolWt, linkage`
- `10_enumeration_qc.json   (2.3 kB json)`
- `11_QL_2026_rows.csv   344 rows x 72 cols
        cols: Index, LNP_ID, Experiment_ID, Formulation_ID, IL_name, IL_SMILES, IL_protonated_SMILES, IL_head_name, IL_head_SMILES, IL_linker_name, IL_linker_SMILES, IL_tail1_name, IL_tail1_SMILES, IL_tail2_name, IL_tail2_SMILES, IL_tail3_name, IL_tail3_SMILES, IL_tail4_name ...`
- `12_QL_2026_qc.json   (1.1 kB json)`
- `13_validation.json   (1.7 kB json)`
- `liu2026_adhm202502225_blocks.csv   34 rows x 16 cols
        cols: block_id, block_class, name_in_paper, SMILES, n_reactive_sites, reactive_site_type, MW_reported, formula_reported, HRMS_reported, CAS_or_vendor, source_figure, source_table, source_page, read_method, read_confidence, ambiguity_note`
- `liu2026_adhm202502225_checks.csv   20 rows x 9 cols
        cols: check_id, check_description, what_it_would_catch, method, expected, observed, result, severity, note`
- `liu2026_adhm202502225_conditions.csv   20 rows x 45 cols
        cols: condition_id, block_ids, naming_rule, degree_of_substitution, reaction_description, library_size_stated, cargo_type, cargo_identity, formulation_molar_ratios_verbatim, IL_mol_pct, helper_lipid_identity, helper_mol_pct, sterol_identity, sterol_mol_pct, peg_lipid_identity, peg_mol_pct, IL_to_nucleic_acid_wt_ratio, N_to_P_ratio ...`
- `reference_structures.csv   1 rows x 6 cols
        cols: lipid_name, reference_SMILES, source_figure, read_on, read_by, note`

- `PROVENANCE.md   (12.0 kB)`

## 5. C10 drawn-structure anchors for this folder

```json
[
 {
  "paper": "ql_2026",
  "lipid_name": "T3A-L5-C10",
  "source_figure": "Fig 2a blocks + Methods",
  "read_on": "2026-10-29",
  "note": "Assembled by hand from the drawn blocks in Figure 2a (T3A tris(3-aminopropyl)amine; L5 = NHS ester of 6-maleimidohexanoate, five CH2; C10 = 1-decanethiol) and the stated one-pot chemistry: amide at each primary amine, thiosuccinimide at each maleimide. Written with a different ring traversal and atom order from enumerate_library.py's fragment assembly, so the comparison tests connectivity rather than echoing the builder.",
  "reference_inchikey": "DKFSQTNBZODLDT-UHFFFAOYSA-N",
  "enumerated_inchikey": "DKFSQTNBZODLDT-UHFFFAOYSA-N",
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
 "DOPE (HL_name)": 344,
 "Cholesterol (CHL_name)": 344
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

