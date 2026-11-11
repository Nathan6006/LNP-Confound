# Review packet — `xuq_32662132`

Generated from live files. Queue says **13** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 32662132.

## Rows read from `curation/xuq_32662132/06_rows.csv` — 13 rows

**79 columns.** Corpus schema is 79; this file has 79 of them.

## Structures
- distinct full InChIKeys: **1**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **0** of 1
- **NEW to the corpus: 1**

## Values
- `Experiment_value`: n=13, mean **7.69231e-08**, sd(pop) **1**, min -2.77, max 0.879
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'XZ_2020': 13}
- `Model_type`: {'primary_human_CD8_T_lymphocyte': 13}
- `Model_target`: {'in_vitro': 13}
- `Cargo_type`: {'FLuc': 13}
- `Experiment_method`: {'luminescence_normalized': 13}
- `Experiment_batching`: {'individual': 13}
- `value_normalization_unit`: {'experiment': 13}

## Quality columns (schema_quality.json)
- `confidence`: 13/13 populated — {'medium': 13}
- `positional_ambiguity`: 13/13 populated — {'0': 13}
- `structure_assumed`: 13/13 populated — {'0': 13}
- `value_is_bound`: 0/13 populated — {'': 13}
- `curation_note`: 13/13 populated — {"Head 93 is not drawn anywhere in this paper's Supporting Information and its one mass spectrum (Figure S1) is illegible at the published raster resolution; the head is 1-(3-aminopropyl)imidazole on the agreement of three external records (AA_2008 and KW_2014 corpus rows, and curation/xuq_32520527's read of a drawn figure). The tail is built from the named reagents of its own synthesis in SI section 1.2. The excipient ratio is the WEIGHT ratio 16:4:1:1 (lipidoid:cholesterol:DOPE:DSPE-PEG); it is NOT converted to mole percent because DSPE-PEG's block length is a distribution, so HL/CHL/PEG_molratio are left blank.": 3, "Head 93 is not drawn anywhere in this paper's Supporting Information and its one mass spectrum (Figure S1) is illegible at the published raster resolution; the head is 1-(3-aminopropyl)imidazole on the agreement of three external records (AA_2008 and KW_2014 corpus rows, and curation/xuq_32520527's read of a drawn figure). The tail is built from the named reagents of its own synthesis in SI section 1.2. The excipient ratio is the WEIGHT ratio 16:8:1:1 (lipidoid:cholesterol:DOPE:DSPE-PEG); it is NOT converted to mole percent because DSPE-PEG's block length is a distribution, so HL/CHL/PEG_molratio are left blank.": 1, "Head 93 is not drawn anywhere in this paper's Supporting Information and its one mass spectrum (Figure S1) is illegible at the published raster resolution; the head is 1-(3-aminopropyl)imidazole on the agreement of three external records (AA_2008 and KW_2014 corpus rows, and curation/xuq_32520527's read of a drawn figure). The tail is built from the named reagents of its own synthesis in SI section 1.2. The excipient ratio is the WEIGHT ratio 16:1:1:1 (lipidoid:cholesterol:DOPE:DSPE-PEG); it is NOT converted to mole percent because DSPE-PEG's block length is a distribution, so HL/CHL/PEG_molratio are left blank.": 1, "Head 93 is not drawn anywhere in this paper's Supporting Information and its one mass spectrum (Figure S1) is illegible at the published raster resolution; the head is 1-(3-aminopropyl)imidazole on the agreement of three external records (AA_2008 and KW_2014 corpus rows, and curation/xuq_32520527's read of a drawn figure). The tail is built from the named reagents of its own synthesis in SI section 1.2. The excipient ratio is the WEIGHT ratio 16:0:1:1 (lipidoid:cholesterol:DOPE:DSPE-PEG); it is NOT converted to mole percent because DSPE-PEG's block length is a distribution, so HL/CHL/PEG_molratio are left blank.": 1, "Head 93 is not drawn anywhere in this paper's Supporting Information and its one mass spectrum (Figure S1) is illegible at the published raster resolution; the head is 1-(3-aminopropyl)imidazole on the agreement of three external records (AA_2008 and KW_2014 corpus rows, and curation/xuq_32520527's read of a drawn figure). The tail is built from the named reagents of its own synthesis in SI section 1.2. The excipient ratio is the WEIGHT ratio 16:4:8:1 (lipidoid:cholesterol:DOPE:DSPE-PEG); it is NOT converted to mole percent because DSPE-PEG's block length is a distribution, so HL/CHL/PEG_molratio are left blank.": 1}

## Folder inventory
  - `curation/xuq_32662132/00_acquisition.json` (3543 B)
  - `curation/xuq_32662132/00_value_sources.json` (11638 B)
  - `curation/xuq_32662132/01_blocks.csv` (2911 B)
  - `curation/xuq_32662132/01_blocks_qc.json` (4283 B)
  - `curation/xuq_32662132/02_corpus_overlap.json` (9177 B)
  - `curation/xuq_32662132/03_enumerated.csv` (5996 B)
  - `curation/xuq_32662132/04_library_membership.csv` (481 B)
  - `curation/xuq_32662132/04_library_membership.json` (4095 B)
  - `curation/xuq_32662132/05_figS2_geometry.json` (13815 B)
  - `curation/xuq_32662132/05_values.csv` (18607 B)
  - `curation/xuq_32662132/06_rows.csv` (16761 B)
  - `curation/xuq_32662132/06_skipped.json` (6922 B)
  - `curation/xuq_32662132/07_qc.json` (13215 B)
  - `curation/xuq_32662132/08_contradictions.csv` (6480 B)
  - `curation/xuq_32662132/08_retracted_suspicions.csv` (5827 B)
  - `curation/xuq_32662132/09_figS1_mass_test.json` (8900 B)
  - `curation/xuq_32662132/DECISIONS.md` (6318 B)
  - `curation/xuq_32662132/PROVENANCE.md` (14368 B)
  - `curation/xuq_32662132/digitize_figS2.py` (9460 B)
  - `curation/xuq_32662132/make_reference.py` (8197 B)
  - `curation/xuq_32662132/measure_figS1.py` (12211 B)
  - `curation/xuq_32662132/reference_structures.csv` (1550 B)
  - `curation/xuq_32662132/stage01_blocks.py` (12785 B)
  - `curation/xuq_32662132/stage02_overlap.py` (5828 B)
  - `curation/xuq_32662132/stage03_enumerate.py` (4220 B)
  - `curation/xuq_32662132/stage04_library.py` (6375 B)
  - `curation/xuq_32662132/stage05_values.py` (7864 B)
  - `curation/xuq_32662132/stage06_rows.py` (8509 B)
  - `curation/xuq_32662132/stage07_qc.py` (22708 B)

## DECISIONS.md — **9 OPEN**
  - D-1 Does a single-lipid formulation DoE belong in the corpus at all?
  - D-2 Should the four Figure S3 %GFP values for 93-O17S be emitted too?
  - D-3 Is `O16B` the same tail as `O16S-S`, and `N16B` the same as `N16S-S`?
  - D-4 What does the "C" in `O17C` designate?
  - D-5 Should this paper be curated at all without the main text?
  - D-6 Is a catalogue number from *other* papers an acceptable structure source?
  - D-7 The four 9322-O17S measurements — now with a molecule.
  - D-8 Weight ratios into a mole-percent schema.
  - D-9 Licence not captured.

## 08_contradictions.csv — 10 rows
