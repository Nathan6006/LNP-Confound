# Review packet — `andersondg_41203968`

Generated from live files. Queue says **132** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 41203968.

## Rows read from `curation/andersondg_41203968/06_emitted_rows.csv` — 132 rows

**26 columns.** Corpus schema is 79; this file has 9 of them.

## Structures
- distinct full InChIKeys: **120**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **0** of 120
- **NEW to the corpus: 120**

## Values
- `value_normalization_unit`: {'': 132}

## Quality columns (schema_quality.json)
- `confidence`: 132/132 populated — {'high': 114, 'medium': 18}
- `positional_ambiguity`: 132/132 populated — {'0': 132}
- `structure_assumed`: 132/132 populated — {'0': 119, '1': 13}
- `value_is_bound`: 2/132 populated — {'': 130, 'upper': 2}
- `curation_note`: 132/132 populated — {'the SI draws no structure for this member; it is built from the headgroup and tail blocks, each of which was derived independently from other lipids in the same series.': 12, "structure read from the vector coordinates of the SI's own drawing; its printed HRMS m/z reproduced to 0.6 ppm; and re-assembling the headgroup and tail blocks reproduces it exactly.": 6, "structure read from the vector coordinates of the SI's own drawing; its printed HRMS m/z reproduced to 0.1 ppm; and re-assembling the headgroup and tail blocks reproduces it exactly.": 6, "structure read from the vector coordinates of the SI's own drawing; its printed HRMS m/z reproduced to 3.2 ppm; and re-assembling the headgroup and tail blocks reproduces it exactly.": 5, "structure read from the vector coordinates of the SI's own drawing; the SI prints no HRMS beside it, so the only external mass check is absent; and re-assembling the headgroup and tail blocks reproduces it exactly.": 5}

## Folder inventory
  - `curation/andersondg_41203968/00_acquisition_raw.json` (5380 B)
  - `curation/andersondg_41203968/00_value_sources.json` (7031 B)
  - `curation/andersondg_41203968/01_blocks.csv` (5231 B)
  - `curation/andersondg_41203968/02_blocks_consensus.json` (5435 B)
  - `curation/andersondg_41203968/02_hrms.json` (16653 B)
  - `curation/andersondg_41203968/02_structures.csv` (29158 B)
  - `curation/andersondg_41203968/03_enumerated.csv` (23999 B)
  - `curation/andersondg_41203968/05_values.csv` (41756 B)
  - `curation/andersondg_41203968/06_emitted_rows.csv` (77380 B)
  - `curation/andersondg_41203968/06_not_emitted.json` (360 B)
  - `curation/andersondg_41203968/07_qc.json` (6692 B)
  - `curation/andersondg_41203968/08_contradictions.csv` (6626 B)
  - `curation/andersondg_41203968/DECISIONS.md` (8698 B)
  - `curation/andersondg_41203968/chemdraw_vector.py` (15755 B)
  - `curation/andersondg_41203968/make_reference.py` (5117 B)
  - `curation/andersondg_41203968/reference_structures.csv` (740 B)
  - `curation/andersondg_41203968/stage01_structures.py` (8423 B)
  - `curation/andersondg_41203968/stage02_blocks.py` (19552 B)
  - `curation/andersondg_41203968/stage03_values.py` (7767 B)
  - `curation/andersondg_41203968/stage04_qc.py` (17409 B)
  - `curation/andersondg_41203968/stage05_emit.py` (6791 B)

## DECISIONS.md — **9 OPEN**
  - D-1 — One experiment or five? The in vivo panels are not one screen.
  - D-2 — Four lipids whose structures this curation could not establish.
  - D-3 — The two benchmarks on Figure 2c.
  - D-4 — "the next best lipid AMG1010" names a lipid that is not in the panel.
  - D-5 — Two cells are exactly 0 RFU and the corpus z-scores log10.
  - D-6 — AMG106 fails its own HRMS, and the paper's formula is impossible.
  - D-7 — Four library members the SI never draws.
  - D-8 — The in vitro screen publishes a mean, not wells.
  - D-9 — Licence, for the §8 release check.

## 08_contradictions.csv — 9 rows
