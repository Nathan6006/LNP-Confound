# Review packet — `whiteheadka_38888972`

Generated from live files. Queue says **49** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 38888972.

## Rows read from `curation/whiteheadka_38888972/05_values.csv` — 49 rows

**74 columns.** Corpus schema is 79; this file has 74 of them.

## Structures
- distinct full InChIKeys: **49**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **41** of 49
  | InChIKey | rows here | in corpus as |
  |---|---|---|
  | `WFJQIMKJKNLLBO-UHFFFAOYSA-N` | 1 | KW_2014:A123_C10, NC_2024:A123_O10 |
  | `ZLXRYQCCHKMLEW-UHFFFAOYSA-N` | 1 | KW_2014:A304_C10, NC_2024:A304_O10 |
  | `QDOPKTMAJAXNAB-UHFFFAOYSA-N` | 1 | KW_2014:A205_C14, NC_2024:A205_O14 |
  | `AZDRBKOPHTYJAM-UHFFFAOYSA-N` | 1 | KW_2014:A111_C10 |
  | `MVJNZLFMDDDBLQ-UHFFFAOYSA-N` | 1 | NC_2024:A211_O10 |
  | `ATUZDWRFUUIPHR-UHFFFAOYSA-N` | 1 | NC_2024:A601_O10 |
  | `IMIGQTSDHBDPSH-UHFFFAOYSA-N` | 1 | KW_2014:A503_C10, NC_2024:A503_O10 |
  | `WANMECLPAVOWMI-UHFFFAOYSA-N` | 1 | KW_2014:A503_C11, NC_2024:A503_O11 |
- **NEW to the corpus: 8**

## Values
- `Experiment_value`: n=49, mean **4.10669e-17**, sd(pop) **1**, min -2.118, max 1.983
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'KW_2024': 49}
- `Model_type`: {'Mouse_B6': 49}
- `Model_target`: {'whole_body': 49}
- `Cargo_type`: {'FLuc': 49}
- `Experiment_method`: {'luminescence_normalized': 49}
- `Experiment_batching`: {'individual': 49}
- `value_normalization_unit`: {'experiment': 49}

## Quality columns (schema_quality.json)
- `confidence`: 49/49 populated — {'high': 33, 'low': 10, 'medium': 6}
- `positional_ambiguity`: 49/49 populated — {'0': 49}
- `structure_assumed`: 49/49 populated — {'0': 33, '1': 16}
- `value_is_bound`: 0/49 populated — {'': 49}
- `curation_note`: 49/49 populated — {'Value: in-vivo total luminescent flux digitised from Figure 2B (log axis, 6 decade ticks, tick-fit residual 0.12 px; the drawn 10^8 reference line lands at 0.997x nominal).': 31, 'The but-2-ene bridge of head 516 is drawn with a WAVY bond in Figure 1B and its Table S1 name gives no E/Z, so the geometry is stored unassigned, as declared. Value: in-vivo total luminescent flux digitised from Figure 2B (log axis, 6 decade ticks, tick-fit residual 0.12 px; the drawn 10^8 reference line lands at 0.997x nominal).': 2, 'Whitehead 2024 lipidoid 303X2: Figure 1B draws head 303 differently from its Table S1 name; the catalogue name is stored (08_contradictions.csv). Value: in-vivo total luminescent flux digitised from Figure 2B (log axis, 6 decade ticks, tick-fit residual 0.12 px; the drawn 10^8 reference line lands at 0.997x nominal).': 1, 'Whitehead 2024 lipidoid 110X7: head 110 has 6 N-H but the synthesis used 4 equivalents of acrylate, so full substitution is a limiting case, not a product. Value: in-vivo total luminescent flux digitised from Figure 2B (log axis, 6 decade ticks, tick-fit residual 0.12 px; the drawn 10^8 reference line lands at 0.997x nominal).': 1, 'Whitehead 2024 lipidoid 111O10: head 111 has 7 N-H but the synthesis used 4 equivalents of acrylate, so full substitution is a limiting case, not a product. Value: in-vivo total luminescent flux digitised from Figure 2B (log axis, 6 decade ticks, tick-fit residual 0.12 px; the drawn 10^8 reference line lands at 0.997x nominal).': 1}

## Folder inventory
  - `curation/whiteheadka_38888972/00_acquisition.json` (6902 B)
  - `curation/whiteheadka_38888972/00_value_sources.json` (9243 B)
  - `curation/whiteheadka_38888972/01_blocks.csv` (12333 B)
  - `curation/whiteheadka_38888972/01_blocks_qc.json` (4696 B)
  - `curation/whiteheadka_38888972/03_enumerated.csv` (241805 B)
  - `curation/whiteheadka_38888972/03_enumerated_qc.json` (620 B)
  - `curation/whiteheadka_38888972/04_fig2b_bars.csv` (3810 B)
  - `curation/whiteheadka_38888972/04_fig2b_qc.json` (1309 B)
  - `curation/whiteheadka_38888972/05_values.csv` (55652 B)
  - `curation/whiteheadka_38888972/05_values_qc.json` (1546 B)
  - `curation/whiteheadka_38888972/06_tns.csv` (37963 B)
  - `curation/whiteheadka_38888972/06_tns_qc.json` (2290 B)
  - `curation/whiteheadka_38888972/07_fig2c_points.csv` (1105 B)
  - `curation/whiteheadka_38888972/07_fig2c_qc.json` (1599 B)
  - `curation/whiteheadka_38888972/07_qc.json` (14816 B)
  - `curation/whiteheadka_38888972/08_contradictions.csv` (9832 B)
  - `curation/whiteheadka_38888972/08_pies.csv` (393 B)
  - `curation/whiteheadka_38888972/08_pies_qc.json` (2038 B)
  - `curation/whiteheadka_38888972/09_seeded_checks.json` (2527 B)
  - `curation/whiteheadka_38888972/09_seeded_checks.py` (8622 B)
  - `curation/whiteheadka_38888972/DECISIONS.md` (12540 B)
  - `curation/whiteheadka_38888972/PROVENANCE.md` (9661 B)
  - `curation/whiteheadka_38888972/make_reference.py` (8089 B)
  - `curation/whiteheadka_38888972/reference_structures.csv` (1300 B)
  - `curation/whiteheadka_38888972/run_all.sh` (969 B)
  - `curation/whiteheadka_38888972/stage01_blocks.py` (16588 B)
  - `curation/whiteheadka_38888972/stage02_fig2a.py` (10399 B)
  - `curation/whiteheadka_38888972/stage03_enumerate.py` (15977 B)
  - `curation/whiteheadka_38888972/stage04_fig2b.py` (11831 B)
  - `curation/whiteheadka_38888972/stage05_fig2c.py` (8287 B)
  - `curation/whiteheadka_38888972/stage06_rows.py` (12157 B)
  - `curation/whiteheadka_38888972/stage07_fig3c_pie.py` (10419 B)
  - `curation/whiteheadka_38888972/stage08_crosschecks.py` (22023 B)

## DECISIONS.md — **0 OPEN**
  - D-1  Three amine heads: the Table S1 NAME and the Figure 1B DRAWING disagree
  - D-2  What to store for the five heads that 4 equivalents cannot fully substitute
  - D-3  Is TNS fluorescence an admissible corpus endpoint?
  - D-4  The in-vitro HeLa values exist but carry no identity
  - D-5  `Dose_ug_nucleicacid` is not derivable
  - D-6  Licence: this paper is not CC
  - D-7  The `Experiment_ID`
  - D-8  Per-organ rows: liver / lungs / spleen are drawn but not recoverable
  - D-9  `Model_target` for the emitted rows

## 08_contradictions.csv — 15 rows
