# Review packet — `sd_2026`

Generated from live files. Queue says **18** emitted rows; status `extracted_below_threshold`; proposed `Experiment_ID` **none**; PMID 42806111.

## Rows read from `curation/sd_2026/03_enumerated.csv` — 384 rows  **(queue says 18 — DISCREPANCY)**

**18 columns.** Corpus schema is 79; this file has 8 of them.

## Structures
- distinct full InChIKeys: **384**
- structure column: **`canonical_SMILES`**
- unparseable structures: **0**
- already in the corpus: **0** of 384
- **NEW to the corpus: 384**

## Values
- `value_normalization_unit`: {'': 384}

## Quality columns (schema_quality.json)
- `confidence`: 384/384 populated — {'low': 383, 'medium': 1}
- `positional_ambiguity`: 384/384 populated — {'0': 384}
- `structure_assumed`: 384/384 populated — {'1': 383, '0': 1}
- `value_is_bound`: 0/384 populated — {'': 384}
- `curation_note`: 384/384 populated — {'block assignment is the combinatorial enumeration the paper describes; this specific member is not named, drawn or measured anywhere in the paper; the Ugi creates a stereocentre at the linker-derived carbon and the reaction is not stereoselective, so this is stored without stereochemistry.': 239, 'block assignment is the combinatorial enumeration the paper describes; this specific member is not named, drawn or measured anywhere in the paper.': 80, 'block assignment is the combinatorial enumeration the paper describes; this specific member is not named, drawn or measured anywhere in the paper; the Ugi creates a stereocentre at the linker-derived carbon and the reaction is not stereoselective, so this is stored without stereochemistry; H8 is drawn with no wedge at its ring stereocentre.': 24, 'block assignment is the combinatorial enumeration the paper describes; this specific member is not named, drawn or measured anywhere in the paper; the Ugi creates a stereocentre at the linker-derived carbon and the reaction is not stereoselective, so this is stored without stereochemistry; H11 is drawn with no wedge at its ring stereocentre.': 24, 'block assignment is the combinatorial enumeration the paper describes; this specific member is not named, drawn or measured anywhere in the paper; H8 is drawn with no wedge at its ring stereocentre.': 8}

## Folder inventory
  - `curation/sd_2026/00_acquisition.json` (8149 B)
  - `curation/sd_2026/00_value_sources.json` (12156 B)
  - `curation/sd_2026/01_blocks.csv` (3247 B)
  - `curation/sd_2026/01_blocks_qc.json` (6680 B)
  - `curation/sd_2026/03_enumerated.csv` (169355 B)
  - `curation/sd_2026/03_enumerated_qc.json` (5497 B)
  - `curation/sd_2026/05_values.csv` (12530 B)
  - `curation/sd_2026/05_values_qc.json` (4122 B)
  - `curation/sd_2026/07_qc.json` (25618 B)
  - `curation/sd_2026/08_contradictions.csv` (5658 B)
  - `curation/sd_2026/41587_2026_3298_MOESM1_ESM.pdf` (47319200 B)
  - `curation/sd_2026/DECISIONS.md` (9080 B)
  - `curation/sd_2026/make_reference.py` (7776 B)
  - `curation/sd_2026/reference_structures.csv` (5135 B)
  - `curation/sd_2026/s41587-026-03298-8.pdf` (7813793 B)
  - `curation/sd_2026/stage01_blocks.py` (23652 B)
  - `curation/sd_2026/stage03_enumerate.py` (14628 B)
  - `curation/sd_2026/stage05_values.py` (16608 B)
  - `curation/sd_2026/stage07_qc.py` (43100 B)

## DECISIONS.md — **8 OPEN**
  - D-1 The 384-lipid screen has values that cannot be read, and structures that can
  - D-2 Which of the 384 members are SC-1…SC-5 and LC-2…LC-5?
  - D-3 Supplementary Fig. 53: four more organs, by raster digitisation
  - D-4 A second, separate lipid library in Supplementary Fig. 68
  - D-5 The licence is not Creative Commons
  - D-6 Three-lipid experiments and the z-score convention
  - D-7 What Fig. 1n's radiance is a radiance *of*

## 08_contradictions.csv — 9 rows
