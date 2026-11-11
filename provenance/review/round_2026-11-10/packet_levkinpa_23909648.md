# Review packet — `levkinpa_23909648`

Generated from live files. Queue says **29** emitted rows; status `extracted`; proposed `Experiment_ID` **LL_2013**; PMID 23909648.

## Rows read from `curation/levkinpa_23909648/06_rows.csv` — 29 rows

**28 columns.** Corpus schema is 79; this file has 17 of them.

## Structures
- distinct full InChIKeys: **3**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **0** of 3
- **NEW to the corpus: 3**

## Values
- `Experiment_ID`: {'LL_2013': 29}
- `Model_type`: {'NOT STATED IN THE SI': 29}
- `Model_target`: {'NOT STATED IN THE SI': 29}
- `value_normalization_unit`: {'': 29}

## Quality columns (schema_quality.json)
- `confidence`: 29/29 populated — {'high': 29}
- `positional_ambiguity`: 29/29 populated — {'0': 29}
- `structure_assumed`: 29/29 populated — {'0': 29}
- `value_is_bound`: 1/29 populated — {'': 28, 'upper': 1}
- `curation_note`: 29/29 populated — {"Structure from the full IUPAC name printed on SI p.13, independently confirmed by HRMS (-0.006 ppm, our computed [M+H]+ vs the authors' calculated m/z); the screened material for the unpurified library members of the same name is a one-pot mixture; THIS row is the column-purified single compound.": 10, "Structure from the full IUPAC name printed on SI p.10, independently confirmed by HRMS (-0.072 ppm, our computed [M+H]+ vs the authors' calculated m/z); the SI's 1H NMR TEXT lists -1 H against the formula C22H46N2, which is a printing error in the paper, not a structural doubt: the spectrum's own integral sum is 44.45 H (see 08_contradictions.csv); the screened material for the unpurified library members of the same name is a one-pot mixture; THIS row is the column-purified single compound.": 9, "Structure from the full IUPAC name printed on SI p.7, independently confirmed by HRMS (-0.053 ppm, our computed [M+H]+ vs the authors' calculated m/z); the SI's 1H NMR TEXT lists -3 H against the formula C30H62N2, which is a printing error in the paper, not a structural doubt: the spectrum's own integral sum is 63.2 H (see 08_contradictions.csv); the screened material for the unpurified library members of the same name is a one-pot mixture; THIS row is the column-purified single compound.": 9, "Structure from the full IUPAC name printed on SI p.10, independently confirmed by HRMS (-0.072 ppm, our computed [M+H]+ vs the authors' calculated m/z); the SI's 1H NMR TEXT lists -1 H against the formula C22H46N2, which is a printing error in the paper, not a structural doubt: the spectrum's own integral sum is 44.45 H (see 08_contradictions.csv); the screened material for the unpurified library members of the same name is a one-pot mixture; THIS row is the column-purified single compound. This row's value is printed as 0 with an SD of 0, i.e. a non-detect: the true transfection efficiency is at or below the assay's detection limit, not exactly zero. It also cannot be log-transformed, which the corpus convention requires -- DECISIONS.md D-5.": 1}

## Folder inventory
  - `curation/levkinpa_23909648/00_acquisition.json` (5538 B)
  - `curation/levkinpa_23909648/00_value_sources.json` (14360 B)
  - `curation/levkinpa_23909648/01_blocks.csv` (5155 B)
  - `curation/levkinpa_23909648/01_blocks.json` (25778 B)
  - `curation/levkinpa_23909648/02_drawn_geometry.json` (6390 B)
  - `curation/levkinpa_23909648/03_enumerated.csv` (1694 B)
  - `curation/levkinpa_23909648/03_enumerated.json` (9104 B)
  - `curation/levkinpa_23909648/03_not_determinable.csv` (3532 B)
  - `curation/levkinpa_23909648/04_corpus_overlap.json` (8505 B)
  - `curation/levkinpa_23909648/05_values.csv` (13901 B)
  - `curation/levkinpa_23909648/05_values_geometry.json` (3248 B)
  - `curation/levkinpa_23909648/06_emit_summary.json` (939 B)
  - `curation/levkinpa_23909648/06_rows.csv` (21376 B)
  - `curation/levkinpa_23909648/06_skipped.csv` (7070 B)
  - `curation/levkinpa_23909648/07_qc.json` (8425 B)
  - `curation/levkinpa_23909648/08_contradictions.csv` (8285 B)
  - `curation/levkinpa_23909648/10_seed_checks.py` (14919 B)
  - `curation/levkinpa_23909648/10_seeded.json` (12737 B)
  - `curation/levkinpa_23909648/DECISIONS.md` (9152 B)
  - `curation/levkinpa_23909648/PROVENANCE.md` (13188 B)
  - `curation/levkinpa_23909648/make_reference.py` (8048 B)
  - `curation/levkinpa_23909648/reference_structures.csv` (2982 B)
  - `curation/levkinpa_23909648/stage01_blocks.py` (17360 B)
  - `curation/levkinpa_23909648/stage02_drawn.py` (11132 B)
  - `curation/levkinpa_23909648/stage03_enumerate.py` (13491 B)
  - `curation/levkinpa_23909648/stage04_corpus_overlap.py` (7297 B)
  - `curation/levkinpa_23909648/stage05_values.py` (16595 B)
  - `curation/levkinpa_23909648/stage06_emit.py` (11292 B)
  - `curation/levkinpa_23909648/stage07_qc.py` (25213 B)

## DECISIONS.md — **8 OPEN**
  - D-1
  - D-2
  - D-3
  - D-4
  - D-5
  - D-6
  - D-7
  - D-8

## 08_contradictions.csv — 10 rows
