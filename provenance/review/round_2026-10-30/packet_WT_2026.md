# Review packet — WT_2026

Generated 2026-10-30 by `provenance/review/brief.py` from the live `data/` and `curation/` trees. **Every number below was read out of the repository at generation time.** If the packet and the repository disagree, the packet is the bug — say so.

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
  "max": 3.4008,
  "mean": "1.154e-17",
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
 "molratios": {
  "IL_molratio": "687 distinct",
  "HL_molratio": "673 distinct",
  "CHL_molratio": "695 distinct",
  "PEG_molratio": "613 distinct"
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
  "Experiment_value": "0.3313630853772711",
  "Model_type": "NIH-3T3",
  "Model": "in_vitro",
  "IL_molratio": "56.0",
  "HL_name": "18PG",
  "HL_molratio": "24.0",
  "CHL_molratio": "19.0",
  "PEG_molratio": "1.0",
  "Publication_PMID": "",
  "Publication_DOI": "10.64898/2026.06.13.731774"
 },
 {
  "Experiment_ID": "WT_2026",
  "IL_name": "SM-102",
  "IL_SMILES": "CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
  "Experiment_value": "0.8641405543991477",
  "Model_type": "HepG2",
  "Model": "in_vitro",
  "IL_molratio": "27.60296073",
  "HL_name": "DSPC",
  "HL_molratio": "9.37403927",
  "CHL_molratio": "57.42151576",
  "PEG_molratio": "5.60148424",
  "Publication_PMID": "",
  "Publication_DOI": "10.64898/2026.06.13.731774"
 },
 {
  "Experiment_ID": "WT_2026",
  "IL_name": "SM-102",
  "IL_SMILES": "CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
  "Experiment_value": "-1.0904581120683585",
  "Model_type": "NIH-3T3",
  "Model": "in_vitro",
  "IL_molratio": "79.2",
  "HL_name": "18PG",
  "HL_molratio": "0.8",
  "CHL_molratio": "19.0",
  "PEG_molratio": "1.0",
  "Publication_PMID": "",
  "Publication_DOI": "10.64898/2026.06.13.731774"
 },
 {
  "Experiment_ID": "WT_2026",
  "IL_name": "SM-102",
  "IL_SMILES": "CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
  "Experiment_value": "0.11808242459243407",
  "Model_type": "NIH-3T3",
  "Model": "in_vitro",
  "IL_molratio": "10.0",
  "HL_name": "18PG",
  "HL_molratio": "10.0",
  "CHL_molratio": "79.84",
  "PEG_molratio": "0.16",
  "Publication_PMID": "",
  "Publication_DOI": "10.64898/2026.06.13.731774"
 },
 {
  "Experiment_ID": "WT_2026",
  "IL_name": "SM-102",
  "IL_SMILES": "CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
  "Experiment_value": "0.8377925304707449",
  "Model_type": "C2C12",
  "Model": "in_vitro",
  "IL_molratio": "18.0",
  "HL_name": "18PG",
  "HL_molratio": "2.0",
  "CHL_molratio": "79.2",
  "PEG_molratio": "0.8",
  "Publication_PMID": "",
  "Publication_DOI": "10.64898/2026.06.13.731774"
 },
 {
  "Experiment_ID": "WT_2026",
  "IL_name": "SM-102",
  "IL_SMILES": "CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC",
  "Experiment_value": "0.8779746603626277",
  "Model_type": "Ramos",
  "Model": "in_vitro",
  "IL_molratio": "35.0",
  "HL_name": "DSPC",
  "HL_molratio": "15.0",
  "CHL_molratio": "49.5",
  "PEG_molratio": "0.5",
  "Publication_PMID": "",
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

- `build_rows.py   (322 lines)  Emit the FALCON datasets (Toh et al., bioRxiv 2026) in the 72-column corpus schema.`
- `make_merge_subset.py   (89 lines)  Collapse the FALCON curation into the single block that merges into data/.`
- `validate.py   (170 lines)  Independent check of 01_WT_2026_rows.csv. Reads the OUTPUT and the SOURCE and`

## 4. Intermediate and output tables

- `01_WT_2026_rows.csv   2544 rows x 72 cols
        cols: Index, LNP_ID, Experiment_ID, Formulation_ID, IL_name, IL_SMILES, IL_protonated_SMILES, IL_head_name, IL_head_SMILES, IL_linker_name, IL_linker_SMILES, IL_tail1_name, IL_tail1_SMILES, IL_tail2_name, IL_tail2_SMILES, IL_tail3_name, IL_tail3_SMILES, IL_tail4_name ...`
- `02_WT_2026_qc.json   (6.2 kB json)`
- `03_enumerated.csv   2 rows x 2 cols
        cols: lipid_name, IL_SMILES`
- `04_validation.json   (1.8 kB json)`
- `05_WT_2026_merge.csv   2544 rows x 72 cols
        cols: Index, LNP_ID, Experiment_ID, Formulation_ID, IL_name, IL_SMILES, IL_protonated_SMILES, IL_head_name, IL_head_SMILES, IL_linker_name, IL_linker_SMILES, IL_tail1_name, IL_tail1_SMILES, IL_tail2_name, IL_tail2_SMILES, IL_tail3_name, IL_tail3_SMILES, IL_tail4_name ...`
- `06_merge_subset_qc.json   (1.2 kB json)`
- `reference_structures.csv   2 rows x 6 cols
        cols: lipid_name, reference_SMILES, source_figure, read_on, read_by, note`

- `PROVENANCE.md   (9.5 kB)`

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

