# Review packet — `levkinpa_29558113`

Generated from live files. Queue says **5** emitted rows; status `extracted_below_threshold`; proposed `Experiment_ID` **none**; PMID 29558113.

## Rows read from `curation/levkinpa_29558113/06_rows.csv` — 45 rows  **(queue says 5 — DISCREPANCY)**

**19 columns.** Corpus schema is 79; this file has 13 of them.

## Structures
- distinct full InChIKeys: **5**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **0** of 5
- **NEW to the corpus: 5**

## Values
- `Experiment_ID`: {'PL_2018': 45}
- `Model_type`: {'HEK293T': 45}
- `Model_target`: {'eGFP plasmid DNA': 45}
- `value_normalization_unit`: {'': 45}

## Quality columns (schema_quality.json)
- `confidence`: 5/45 populated — {'': 40, 'high': 3, 'medium': 2}
- `positional_ambiguity`: 5/45 populated — {'': 40, '0': 5}
- `structure_assumed`: 5/45 populated — {'': 40, '0': 3, '1': 2}
- `value_is_bound`: 0/45 populated — {'': 45}
- `curation_note`: 5/45 populated — {'': 40, 'the alpha carbon is racemic -- the core is DL-homocysteine thiolactone (SI p.2), so the stored SMILES carries no stereo descriptor there by design': 3, "A5's three-carbon linker is measured from the Figure S1b drawing but its two N-ethyl branches are read, not measured (the raster is at its noise floor); the paper gives A5 no name, formula or mass, so there is no anchor orthogonal to the drawing; the alpha carbon is racemic -- the core is DL-homocysteine thiolactone (SI p.2), so the stored SMILES carries no stereo descriptor there by design": 2}

## Folder inventory
  - `curation/levkinpa_29558113/00_acquisition.json` (3291 B)
  - `curation/levkinpa_29558113/00_value_sources.json` (7765 B)
  - `curation/levkinpa_29558113/01_amines.json` (5542 B)
  - `curation/levkinpa_29558113/01_blocks.csv` (5147 B)
  - `curation/levkinpa_29558113/02_chain_counts.json` (5239 B)
  - `curation/levkinpa_29558113/03_block_hrms.json` (1311 B)
  - `curation/levkinpa_29558113/03_enumerated.csv` (24826 B)
  - `curation/levkinpa_29558113/03_reconciliation.json` (614 B)
  - `curation/levkinpa_29558113/04_figS3_verification.json` (7255 B)
  - `curation/levkinpa_29558113/05_values.csv` (16371 B)
  - `curation/levkinpa_29558113/05_values_summary.json` (771 B)
  - `curation/levkinpa_29558113/06_emit_summary.json` (951 B)
  - `curation/levkinpa_29558113/06_rows.csv` (10032 B)
  - `curation/levkinpa_29558113/07_qc.json` (12697 B)
  - `curation/levkinpa_29558113/08_contradictions.csv` (6139 B)
  - `curation/levkinpa_29558113/09_anchor_build.json` (643 B)
  - `curation/levkinpa_29558113/10_seed_checks.py` (8853 B)
  - `curation/levkinpa_29558113/10_seeded.json` (6358 B)
  - `curation/levkinpa_29558113/DECISIONS.md` (8775 B)
  - `curation/levkinpa_29558113/make_reference.py` (8589 B)
  - `curation/levkinpa_29558113/reference_structures.csv` (1442 B)
  - `curation/levkinpa_29558113/stage01_chains.py` (10115 B)
  - `curation/levkinpa_29558113/stage02_amines.py` (15531 B)
  - `curation/levkinpa_29558113/stage03_enumerate.py` (13434 B)
  - `curation/levkinpa_29558113/stage04_figS3.py` (14045 B)
  - `curation/levkinpa_29558113/stage05_values.py` (6993 B)
  - `curation/levkinpa_29558113/stage06_emit.py` (5617 B)
  - `curation/levkinpa_29558113/stage07_qc.py` (13472 B)

## DECISIONS.md — **9 OPEN**
  - D-1
  - D-2
  - D-3
  - D-4
  - D-5
  - D-6
  - D-7
  - D-8
  - D-9

## 08_contradictions.csv — 9 rows
