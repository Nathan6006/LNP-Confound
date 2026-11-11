# Review packet — `khanof_36652629`

Generated from live files. Queue says **18** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 36652629.

## Rows read from `curation/khanof_36652629/05_values.csv` — 18 rows

**38 columns.** Corpus schema is 79; this file has 9 of them.

## Structures
- distinct full InChIKeys: **10**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **2** of 10
  | InChIKey | rows here | in corpus as |
  |---|---|---|
  | `BGNVBNJYBVCBJH-UHFFFAOYSA-N` | 4 | AC_2025:SM-102, SB_2024:SM-102, SS_2018:SM-102 |
  | `NRLNQCOGCKAESA-KWXKLSQISA-N` | 1 | AC_2025:DLin-MC3-DMA, BL_2023:DLin-MC3-DMA, NG_2025_2:DLin-MC3-DMA |
- **NEW to the corpus: 8**

## Values
- `value_normalization_unit`: {'': 18}

## Quality columns (schema_quality.json)
- `confidence`: 18/18 populated — {'high': 18}
- `positional_ambiguity`: 18/18 populated — {'0': 18}
- `structure_assumed`: 18/18 populated — {'0': 18}
- `value_is_bound`: 0/18 populated — {'': 18}
- `curation_note`: 18/18 populated — {'Benchmark lipid; structure from the project reagent registry (C11), not from this paper. Value digitised from SI p2, Figure S1B; it is RAW total flux in p/s and is NOT comparable across panels -- Table S2 gives this panel formulation 1 and other panels a different one.': 1, "Structure enumerated from the SI synthesis; our formula C69H141N5O2 equals the SI's printed formula. Value digitised from SI p2, Figure S1B; it is RAW total flux in p/s and is NOT comparable across panels -- Table S2 gives this panel formulation 2 and other panels a different one.": 1, "Structure enumerated from the SI synthesis; our formula C69H141N5O6 equals the SI's printed formula. Stored without stereochemistry: the four carbinol centres come from ring-opening a RACEMIC 1,2-epoxyalkane, so the compound is a mixture of stereoisomers and no single configuration is correct. Value digitised from SI p2, Figure S1B; it is RAW total flux in p/s and is NOT comparable across panels -- Table S2 gives this panel formulation 3 and other panels a different one.": 1, "Structure enumerated from the SI synthesis; our formula C69H141N5O6 equals the SI's printed formula. Stored without stereochemistry: the four carbinol centres come from ring-opening a RACEMIC 1,2-epoxyalkane, so the compound is a mixture of stereoisomers and no single configuration is correct. Value digitised from SI p3, Figure S2B; it is RAW total flux in p/s and is NOT comparable across panels -- Table S2 gives this panel formulation 4 and other panels a different one.": 1, "Structure enumerated from the SI synthesis; our formula C53H109N5O6 equals the SI's printed formula. Stored without stereochemistry: the four carbinol centres come from ring-opening a RACEMIC 1,2-epoxyalkane, so the compound is a mixture of stereoisomers and no single configuration is correct. Value digitised from SI p3, Figure S2B; it is RAW total flux in p/s and is NOT comparable across panels -- Table S2 gives this panel formulation 6 and other panels a different one.": 1}

## Folder inventory
  - `curation/khanof_36652629/00_scan_for_values.json` (2552 B)
  - `curation/khanof_36652629/00_value_sources.json` (10006 B)
  - `curation/khanof_36652629/01_blocks.csv` (1035 B)
  - `curation/khanof_36652629/03_enumerated.csv` (2804 B)
  - `curation/khanof_36652629/04_corpus_overlap.json` (2204 B)
  - `curation/khanof_36652629/05_digitize.json` (19182 B)
  - `curation/khanof_36652629/05_values.csv` (15692 B)
  - `curation/khanof_36652629/07_qc.json` (25549 B)
  - `curation/khanof_36652629/08_contradictions.csv` (9425 B)
  - `curation/khanof_36652629/DECISIONS.md` (6258 B)
  - `curation/khanof_36652629/PROVENANCE.md` (14605 B)
  - `curation/khanof_36652629/make_reference.py` (10251 B)
  - `curation/khanof_36652629/reference_structures.csv` (4886 B)
  - `curation/khanof_36652629/stage00_scan_for_values.py` (4523 B)
  - `curation/khanof_36652629/stage01_enumerate.py` (11815 B)
  - `curation/khanof_36652629/stage04_corpus_overlap.py` (3732 B)
  - `curation/khanof_36652629/stage05_digitize.py` (21924 B)
  - `curation/khanof_36652629/stage06_values.py` (8580 B)
  - `curation/khanof_36652629/stage07_qc.py` (29294 B)

## DECISIONS.md — **7 OPEN**
  - D-1 Six rounds, or one experiment?
  - D-2 The only shared chemistry with the corpus is benchmarks
  - D-3 The mRNA dose is nowhere in the material we hold
  - D-4 Fig S9A disagrees with the panels it is derived from
  - D-5 Fig S4D's OC2-K3-E10 characterisation matches no OC2-K3-E10 formulation
  - D-6 Fig S4C uses a different cargo
  - D-7 Six lipids are stereoisomer mixtures stored as flat SMILES

## 08_contradictions.csv — 11 rows
