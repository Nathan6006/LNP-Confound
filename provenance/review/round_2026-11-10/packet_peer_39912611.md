# Review packet — `peer_39912611`

Generated from live files. Queue says **144** emitted rows; status `extracted`; proposed `Experiment_ID` **GN_2025**; PMID 39912611.

## Rows read from `curation/peer_39912611/06_rows.csv` — 144 rows

**40 columns.** Corpus schema is 79; this file has 37 of them.

## Structures
- distinct full InChIKeys: **12**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **0** of 12
- **NEW to the corpus: 12**

## Values
- `Experiment_value`: n=144, mean **2.33706e-18**, sd(pop) **1**, min -2.697, max 1.979
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'GN_2025': 144}
- `Model_type`: {'HeLa': 36, 'RAW264.7': 36, 'HepG2': 36, 'Mouse_C57BL6J_female_8to10wk': 36}
- `Model_target`: {'in_vitro': 108, 'lungs': 12, 'liver': 12, 'spleen': 12}
- `Cargo_type`: {'mRNA': 144}
- `Experiment_method`: {'luminescence_normalized': 108, 'radiance_normalized': 36}
- `Experiment_batching`: {'individual': 144}
- `value_normalization_unit`: {'experiment': 144}

## Quality columns (schema_quality.json)
- `confidence`: 144/144 populated — {'medium': 138, 'low': 6}
- `positional_ambiguity`: 0/144 populated — {'': 144}
- `structure_assumed`: 0/144 populated — {'': 144}
- `value_is_bound`: 6/144 populated — {'': 138, 'upper': 6}
- `curation_note`: 144/144 populated — {'value digitised from a 700 px publisher raster (the only one reachable without a browser); 1 px is 0.02 decades in Fig 2A and 0.04 in Fig 2B/2C': 126, "structure has no IUPAC name in this SI (lipid 14 is from the authors' ref 41); it is read off the drawn scheme Supp Fig 15 and corroborated against ZT_2026's Lipid 14-CB1/2/3 in the corpus; value digitised from a 700 px publisher raster (the only one reachable without a browser); 1 px is 0.02 decades in Fig 2A and 0.04 in Fig 2B/2C": 12, 'no bar is drawn for this lipid in its panel, so the value is at or below the axis floor; stored at the floor as an UPPER BOUND, not a measurement; value digitised from a 700 px publisher raster (the only one reachable without a browser); 1 px is 0.02 decades in Fig 2A and 0.04 in Fig 2B/2C': 6}

## Folder inventory
  - `curation/peer_39912611/00_acquisition.json` (7588 B)
  - `curation/peer_39912611/00_value_sources.json` (9131 B)
  - `curation/peer_39912611/01_blocks.csv` (552 B)
  - `curation/peer_39912611/02_corpus_overlap.json` (2175 B)
  - `curation/peer_39912611/03_enumerated.csv` (4091 B)
  - `curation/peer_39912611/04_fig2_qc.json` (2852 B)
  - `curation/peer_39912611/04_fig2_values.csv` (7914 B)
  - `curation/peer_39912611/05_fig3_qc.json` (5767 B)
  - `curation/peer_39912611/05_fig3_values.csv` (3270 B)
  - `curation/peer_39912611/05_values.csv` (66428 B)
  - `curation/peer_39912611/06_rows.csv` (95815 B)
  - `curation/peer_39912611/07_qc.json` (19248 B)
  - `curation/peer_39912611/08_contradictions.csv` (4854 B)
  - `curation/peer_39912611/DECISIONS.md` (5531 B)
  - `curation/peer_39912611/PROVENANCE.md` (11225 B)
  - `curation/peer_39912611/make_reference.py` (10439 B)
  - `curation/peer_39912611/reference_structures.csv` (8799 B)
  - `curation/peer_39912611/stage01_structures.py` (12460 B)
  - `curation/peer_39912611/stage04_fig2.py` (13478 B)
  - `curation/peer_39912611/stage05_fig3.py` (14586 B)
  - `curation/peer_39912611/stage06_emit.py` (9452 B)
  - `curation/peer_39912611/stage07_checks.py` (28963 B)

## DECISIONS.md — **6 OPEN**
  - D-1 The in-vitro dose series: 0.5/0.25/0.125 (Methods) or 0.25/0.15/0.025 (Figure 2)?
  - D-2 Which `Experiment_ID`, and should the two pools share one?
  - D-3 The six upper-bound rows
  - D-4 The Untreated control bar is not measured
  - D-5 Stereochemistry of the 2-hexyldecanoate centres
  - D-6 `value_normalization_unit` and the two multi-model pools

## 08_contradictions.csv — 8 rows
