# Review packet — `linyx_41851281`

Generated from live files. Queue says **140** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 41851281.

## Rows read from `curation/linyx_41851281/05_values.csv` — 140 rows

**33 columns.** Corpus schema is 79; this file has 7 of them.

## Structures
- distinct full InChIKeys: **119**
- structure column: **`smiles`**
- unparseable structures: **0**
- already in the corpus: **18** of 119
  | InChIKey | rows here | in corpus as |
  |---|---|---|
  | `NRLNQCOGCKAESA-KWXKLSQISA-N` | 1 | AC_2025:DLin-MC3-DMA, BL_2023:DLin-MC3-DMA, NG_2025_2:DLin-MC3-DMA |
  | `QGWBEETXHOVFQS-UHFFFAOYSA-N` | 2 | AC_2025:ALC-0315, LZ_2025:ALC-0315, XD_2026:ALC-0315 |
  | `BGNVBNJYBVCBJH-UHFFFAOYSA-N` | 1 | AC_2025:SM-102, SB_2024:SM-102, SS_2018:SM-102 |
  | `DEGZKZQHROPBLO-UHFFFAOYSA-N` | 1 | AA_2008:80_O10, KW_2014:A80_C10 |
  | `GPQCDUXAUCIAHP-UHFFFAOYSA-N` | 1 | KW_2014:A113_C14 |
  | `ZKURLUHVSGHTAE-UHFFFAOYSA-N` | 1 | KW_2014:A113_C12 |
  | `YDGXKCWVUGEOBQ-UHFFFAOYSA-N` | 1 | KW_2014:A113_C10 |
  | `KSLFMYIEIFRPOW-UHFFFAOYSA-N` | 1 | AA_2008:109_O18 |
- **NEW to the corpus: 101**

## Values
- `value_normalization_unit`: {'': 140}

## Quality columns (schema_quality.json)
- `confidence`: 140/140 populated — {'high': 134, 'medium': 4, 'low': 2}
- `positional_ambiguity`: 140/140 populated — {'0': 136, '1': 4}
- `structure_assumed`: 140/140 populated — {'0': 138, '1': 2}
- `value_is_bound`: 0/140 populated — {'': 140}
- `curation_note`: 140/140 populated — {"Enumerated structure reproduces the paper's printed exact mass to +0.0000 Da.": 67, "Enumerated structure reproduces the paper's printed exact mass to -0.0000 Da.": 62, '2 arms on a head with 4 N-H sites: 2 distinct regiochemistries fit the formula and the paper neither draws nor states which. The maximal-dispersion convention was used. DECISIONS.md D-1.': 3, "Commercial benchmark; structure taken from the corpus reagent registry, not from this paper. Formulated DIFFERENTLY from the library lipids: ALC-0315:DSPC:Chol:ALC-0159 = 46.3:9.4:42.7:1.6, against the library's ionisable lipid:DOPE:Chol:DMG-PEG2000 = 45:15:38.5:1.5.": 2, 'Structure NOT derived: this lipid is outside the H-T combinatorial library and is drawn only in Supplementary Fig. 12. Its arm formula is pinned by HRMS (see 07_qc.json) but the connectivity was not re-derived here. DECISIONS.md D-5.': 2}

## Folder inventory
  - `curation/linyx_41851281/00_acquisition_raw.json` (7505 B)
  - `curation/linyx_41851281/00_names.json` (2458 B)
  - `curation/linyx_41851281/00_value_sources.json` (10974 B)
  - `curation/linyx_41851281/01_blocks.csv` (7920 B)
  - `curation/linyx_41851281/02_sfig1a.png` (1309275 B)
  - `curation/linyx_41851281/02_trace.json` (65465 B)
  - `curation/linyx_41851281/02_trace_overlay.png` (1497482 B)
  - `curation/linyx_41851281/03_enumerated.csv` (23358 B)
  - `curation/linyx_41851281/05_values.csv` (77991 B)
  - `curation/linyx_41851281/07_qc.json` (16031 B)
  - `curation/linyx_41851281/08_contradictions.csv` (6839 B)
  - `curation/linyx_41851281/DECISIONS.md` (6510 B)
  - `curation/linyx_41851281/make_reference.py` (6254 B)
  - `curation/linyx_41851281/reference_structures.csv` (1595 B)
  - `curation/linyx_41851281/stage01_blocks.py` (13206 B)
  - `curation/linyx_41851281/stage02_trace.py` (9408 B)
  - `curation/linyx_41851281/stage03_enumerate.py` (9928 B)
  - `curation/linyx_41851281/stage05_values.py` (13361 B)
  - `curation/linyx_41851281/stage07_qc.py` (19391 B)

## DECISIONS.md — **6 OPEN**
  - D-1
  - D-2
  - D-3
  - D-4
  - D-5
  - D-6
  - D-7

## 08_contradictions.csv — 11 rows
