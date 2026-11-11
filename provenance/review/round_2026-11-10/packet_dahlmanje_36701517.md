# Review packet — `dahlmanje_36701517`

Generated from live files. Queue says **37** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 36701517.

## Rows read from `curation/dahlmanje_36701517/06_rows.csv` — 37 rows

**79 columns.** Corpus schema is 79; this file has 79 of them.

## Structures
- distinct full InChIKeys: **1**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **1** of 1
  | InChIKey | rows here | in corpus as |
  |---|---|---|
  | `ZISVTYVLWSZJAL-UHFFFAOYSA-N` | 37 | AC_2025:CKK-E12, BL_2023:MD-1, JW_2024:other |
- **NEW to the corpus: 0**

## Values
- `Experiment_value`: n=37, mean **-1.08108e-07**, sd(pop) **1**, min -1.354, max 1.723
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'AR_2023': 37}
- `Model_type`: {'Mouse_Ai14': 37}
- `Model_target`: {'liver': 13, 'lung_epithelium': 9, 'kidney': 6, 'heart': 5, 'spleen': 4}
- `Cargo_type`: {'Cre': 37}
- `Experiment_method`: {'editing_efficiency_normalized': 37}
- `Experiment_batching`: {'individual': 37}
- `value_normalization_unit`: {'experiment': 37}

## Quality columns (schema_quality.json)
- `confidence`: 37/37 populated — {'medium': 37}
- `positional_ambiguity`: 0/37 populated — {'': 37}
- `structure_assumed`: 37/37 populated — {'1': 37}
- `value_is_bound`: 0/37 populated — {'': 37}
- `curation_note`: 37/37 populated — {'cell type ECs; Cat-LNP (helper lipid 18:0 DDAB). cKK-E12 is the only ionizable lipid in this paper and the paper draws no structure for it; the stored SMILES is the corpus string, which carries no stereocentres although cKK-E12 has two. Value digitised from Fig S13 at a measured 220 dpi.': 4, 'cell type ECs; Neu-LNP (helper lipid DSPC). cKK-E12 is the only ionizable lipid in this paper and the paper draws no structure for it; the stored SMILES is the corpus string, which carries no stereocentres although cKK-E12 has two. Value digitised from Fig S13 at a measured 220 dpi.': 4, 'cell type DCs; Cat-LNP (helper lipid 18:0 DDAB). cKK-E12 is the only ionizable lipid in this paper and the paper draws no structure for it; the stored SMILES is the corpus string, which carries no stereocentres although cKK-E12 has two. Value digitised from Fig S13 at a measured 220 dpi.': 3, 'cell type DCs; Neu-LNP (helper lipid DSPC). cKK-E12 is the only ionizable lipid in this paper and the paper draws no structure for it; the stored SMILES is the corpus string, which carries no stereocentres although cKK-E12 has two. Value digitised from Fig S13 at a measured 220 dpi.': 3, 'cell type Monocytes; Cat-LNP (helper lipid 18:0 DDAB). cKK-E12 is the only ionizable lipid in this paper and the paper draws no structure for it; the stored SMILES is the corpus string, which carries no stereocentres although cKK-E12 has two. Value digitised from Fig S13 at a measured 220 dpi.': 2}

## Folder inventory
  - `curation/dahlmanje_36701517/00_acquisition.json` (5109 B)
  - `curation/dahlmanje_36701517/00_value_sources.json` (8682 B)
  - `curation/dahlmanje_36701517/01_blocks.csv` (2587 B)
  - `curation/dahlmanje_36701517/02_corpus_overlap.json` (7615 B)
  - `curation/dahlmanje_36701517/03_enumerated.csv` (1138 B)
  - `curation/dahlmanje_36701517/04_digitized.json` (61548 B)
  - `curation/dahlmanje_36701517/04b_figS6c_cationic.json` (4278 B)
  - `curation/dahlmanje_36701517/05_values.csv` (42922 B)
  - `curation/dahlmanje_36701517/06_rows.csv` (31640 B)
  - `curation/dahlmanje_36701517/07_qc.json` (35315 B)
  - `curation/dahlmanje_36701517/08_contradictions.csv` (7168 B)
  - `curation/dahlmanje_36701517/09_seedgate.json` (6301 B)
  - `curation/dahlmanje_36701517/DECISIONS.md` (6471 B)
  - `curation/dahlmanje_36701517/PROVENANCE.md` (6513 B)
  - `curation/dahlmanje_36701517/_qc_narrative_sections.json` (6906 B)
  - `curation/dahlmanje_36701517/make_reference.py` (14615 B)
  - `curation/dahlmanje_36701517/reference_structures.csv` (1432 B)
  - `curation/dahlmanje_36701517/seed_runner.py` (5051 B)
  - `curation/dahlmanje_36701517/stage01_digitize.py` (19475 B)
  - `curation/dahlmanje_36701517/stage02_overlap.py` (7635 B)
  - `curation/dahlmanje_36701517/stage03_emit.py` (29138 B)
  - `curation/dahlmanje_36701517/stage04_selfcheck.py` (13441 B)
  - `curation/dahlmanje_36701517/stage05_seedgate.py` (7498 B)
  - `curation/dahlmanje_36701517/stage06_finalqc.py` (5012 B)

## DECISIONS.md — **7 OPEN**
  - D-1 Does a paper with ONE ionizable lipid belong in a corpus keyed on ionizable-lipid structure?
  - D-2 The pooled barcode screen: admissible in principle, not transcribed here
  - D-3 The aVHH arm (Fig S17) is a second cargo-expression readout and is not extracted
  - D-4 The schema has no cell-type column, so these rows duplicate on their key
  - D-5 `Cargo_type = "Cre"` is a new vocabulary value
  - D-6 cKK-E12's stereochemistry is not pinned by anything in this paper
  - D-7 Six An-LNP liver values are unreadable and are not reconstructed

## 08_contradictions.csv — 8 rows
