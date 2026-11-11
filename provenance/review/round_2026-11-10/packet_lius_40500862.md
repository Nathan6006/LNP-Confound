# Review packet — `lius_40500862`

Generated from live files. Queue says **4** emitted rows; status `extracted_below_threshold`; proposed `Experiment_ID` **none**; PMID 40500862.

## Rows read from `curation/lius_40500862/05_values.csv` — 23 rows  **(queue says 4 — DISCREPANCY)**

**30 columns.** Corpus schema is 79; this file has 11 of them.

## Structures
- distinct full InChIKeys: **13**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **0** of 13
- **NEW to the corpus: 13**

## Values
- `Model_type`: {'C57BL/6 mouse': 19, 'IGROV1': 4}
- `Model_target`: {'spleen_over_liver_ratio': 15, 'luminescence_normalized': 4, 'UNSTATED': 4}
- `value_normalization_unit`: {'': 23}

## Quality columns (schema_quality.json)
- `confidence`: 23/23 populated — {'high': 18, 'medium': 5}
- `positional_ambiguity`: 23/23 populated — {'0': 23}
- `structure_assumed`: 23/23 populated — {'0': 18, '1': 5}
- `value_is_bound`: 0/23 populated — {'': 23}
- `curation_note`: 23/23 populated — {"Structure confirmed by the SI's own HRMS on p.20: our independently assembled SMILES gives [M+H]+ 789.6199 against the printed calculated 789.6199 (-0.03 ppm).": 4, "Structure confirmed by the SI's own HRMS on p.22: our independently assembled SMILES gives [M+H]+ 689.5674 against the printed calculated 689.5674 (+0.06 ppm). VALUE TIE: in Figure S22C the top rules of A5-CE-C7-6 and A5-C7-6 are detected on the same pixel row (230), so their stored values are identical by construction; at 69.0 px per decade the panel resolves no better than 3.4% and cannot order these two.": 4, "Structure confirmed by the SI's own HRMS on p.31: our independently assembled SMILES gives [M+H]+ 741.6351 against the printed calculated 741.6351 (+0.04 ppm).": 3, "Structure confirmed by the SI's own HRMS on p.32: our independently assembled SMILES gives [M+H]+ 657.5776 against the printed calculated 657.5776 (+0.02 ppm). VALUE TIE: in Figure S22C the top rules of A5-CE-C7-6 and A5-C7-6 are detected on the same pixel row (230), so their stored values are identical by construction; at 69.0 px per decade the panel resolves no better than 3.4% and cannot order these two.": 3, 'ASSUMED: the number of carbonate tails. A5-CE-C6 is named only in Figure S5 and has no HRMS, NMR or drawn structure anywhere in the SI; 3 tails is inferred from the linear members of its own series that do carry HRMS. The headgroup, the carbonate linkage and the alcohol identity are not assumed.': 1}

## Folder inventory
  - `curation/lius_40500862/00_value_sources.json` (12035 B)
  - `curation/lius_40500862/01_blocks.csv` (1875 B)
  - `curation/lius_40500862/01_digitised.json` (13964 B)
  - `curation/lius_40500862/02_structure_qc.json` (1410 B)
  - `curation/lius_40500862/03_enumerated.csv` (2966 B)
  - `curation/lius_40500862/04_corpus_overlap.json` (727 B)
  - `curation/lius_40500862/05_values.csv` (16344 B)
  - `curation/lius_40500862/07_qc.json` (9241 B)
  - `curation/lius_40500862/08_contradictions.csv` (7374 B)
  - `curation/lius_40500862/DECISIONS.md` (7071 B)
  - `curation/lius_40500862/PROVENANCE.md` (5115 B)
  - `curation/lius_40500862/make_reference.py` (10143 B)
  - `curation/lius_40500862/reference_structures.csv` (2031 B)
  - `curation/lius_40500862/stage01_values.py` (15240 B)
  - `curation/lius_40500862/stage02_structures.py` (12797 B)
  - `curation/lius_40500862/stage03_emit.py` (8707 B)
  - `curation/lius_40500862/stage07_checks.py` (13939 B)
  - `curation/lius_40500862/stage07_seed_test.py` (5827 B)

## DECISIONS.md — **8 OPEN**
  - D-1
  - D-2
  - D-3
  - D-4
  - D-5
  - D-6
  - D-7
  - D-8

## 08_contradictions.csv — 8 rows
