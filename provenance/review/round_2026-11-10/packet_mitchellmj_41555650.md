# Review packet — `mitchellmj_41555650`

Generated from live files. Queue says **48** emitted rows; status `extracted`; proposed `Experiment_ID` **EH_2026**; PMID 41555650.

## Rows read from `curation/mitchellmj_41555650/06_rows.csv` — 48 rows

**46 columns.** Corpus schema is 79; this file has 29 of them.

## Structures
- distinct full InChIKeys: **5**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **2** of 5
  | InChIKey | rows here | in corpus as |
  |---|---|---|
  | `RVHYPUORVDKRTM-UHFFFAOYSA-N` | 6 | AC_2025:C12-200, JW_2024:other |
  | `BGNVBNJYBVCBJH-UHFFFAOYSA-N` | 6 | AC_2025:SM-102, SB_2024:SM-102, SS_2018:SM-102 |
- **NEW to the corpus: 3**

## Values
- `Experiment_value`: n=48, mean **4.16667e-08**, sd(pop) **1**, min -1.305, max 2.912
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'EH_2026': 48}
- `Model_type`: {'in_vivo': 48}
- `Model_target`: {'brain': 8, 'heart': 8, 'kidney': 8, 'liver': 8, 'lung': 8, 'spleen': 8}
- `Experiment_method`: {'luminescence_normalized': 48}
- `Experiment_batching`: {'individual': 48}
- `value_normalization_unit`: {'experiment': 48}

## Quality columns (schema_quality.json)
- `confidence`: 48/48 populated — {'high': 48}
- `positional_ambiguity`: 48/48 populated — {'0': 48}
- `structure_assumed`: 48/48 populated — {'0': 48}
- `value_is_bound`: 0/48 populated — {'': 48}
- `curation_note`: 48/48 populated — {"4 stereocentres are unassigned because the 2-alkyloxirane reagent is racemic, which is the paper's chemistry and not an omission; value digitised from SI Fig S2A (tier 5) and cross-validated against SI Fig S1A at median 6.2% (check P11).": 6, "Table S1's four molar-ratio figures for b-LNP 14 sum to 110.0, not 100, so they are parts rather than mol% and the true molar composition is not recoverable (see DECISIONS.md D-5); stored verbatim as printed; 4 stereocentres are unassigned because the 2-alkyloxirane reagent is racemic, which is the paper's chemistry and not an omission; value digitised from SI Fig S2A (tier 5) and cross-validated against SI Fig S1A at median 6.2% (check P11).": 6, "5 stereocentres are unassigned because the 2-alkyloxirane reagent is racemic, which is the paper's chemistry and not an omission; value digitised from SI Fig S2A (tier 5) and cross-validated against SI Fig S1A at median 6.2% (check P11).": 6, "Table S1's four molar-ratio figures for b-LNP 70 sum to 102.5, not 100, so they are parts rather than mol% and the true molar composition is not recoverable (see DECISIONS.md D-5); stored verbatim as printed; 5 stereocentres are unassigned because the 2-alkyloxirane reagent is racemic, which is the paper's chemistry and not an omission; value digitised from SI Fig S2A (tier 5) and cross-validated against SI Fig S1A at median 6.2% (check P11).": 6, "Table S1's four molar-ratio figures for b-LNP 96 sum to 102.5, not 100, so they are parts rather than mol% and the true molar composition is not recoverable (see DECISIONS.md D-5); stored verbatim as printed; 7 stereocentres are unassigned because the 2-alkyloxirane reagent is racemic, which is the paper's chemistry and not an omission; value digitised from SI Fig S2A (tier 5) and cross-validated against SI Fig S1A at median 6.2% (check P11).": 6}

## Folder inventory
  - `curation/mitchellmj_41555650/00_value_sources.json` (10544 B)
  - `curation/mitchellmj_41555650/01_blocks.csv` (1679 B)
  - `curation/mitchellmj_41555650/02_corpus_overlap.json` (52257 B)
  - `curation/mitchellmj_41555650/02_overlap.py` (12388 B)
  - `curation/mitchellmj_41555650/03_enumerate.py` (10301 B)
  - `curation/mitchellmj_41555650/03_enumerated.csv` (4880 B)
  - `curation/mitchellmj_41555650/04_digitize_figS2A.py` (19466 B)
  - `curation/mitchellmj_41555650/04_figS1S2_measurements.json` (38761 B)
  - `curation/mitchellmj_41555650/05_tables.json` (18600 B)
  - `curation/mitchellmj_41555650/05_values.csv` (305260 B)
  - `curation/mitchellmj_41555650/05_values.py` (15674 B)
  - `curation/mitchellmj_41555650/05_values_qc.json` (2908 B)
  - `curation/mitchellmj_41555650/06_rows.csv` (49412 B)
  - `curation/mitchellmj_41555650/06_rows.py` (9827 B)
  - `curation/mitchellmj_41555650/07_build_qc.py` (13059 B)
  - `curation/mitchellmj_41555650/07_qc.json` (33928 B)
  - `curation/mitchellmj_41555650/07_seed_checks.py` (14015 B)
  - `curation/mitchellmj_41555650/07_seed_results.json` (4197 B)
  - `curation/mitchellmj_41555650/08_contradictions.csv` (5433 B)
  - `curation/mitchellmj_41555650/09_pool_membership.json` (9676 B)
  - `curation/mitchellmj_41555650/09_pool_membership.py` (5577 B)
  - `curation/mitchellmj_41555650/DECISIONS.md` (9314 B)
  - `curation/mitchellmj_41555650/PROVENANCE.md` (6466 B)
  - `curation/mitchellmj_41555650/make_reference.py` (9448 B)
  - `curation/mitchellmj_41555650/reference_report.json` (2345 B)
  - `curation/mitchellmj_41555650/reference_structures.csv` (4731 B)
  - `curation/mitchellmj_41555650/work/p13_figS1.png` (123770 B)
  - `curation/mitchellmj_41555650/work/p14_figS2.png` (139622 B)

## DECISIONS.md — **8 OPEN**
  - D-1. The headline screen is a POOLED mRNA-BARCODE experiment. Is any of it corpus-admissible?
  - D-2. The pooled screen's normalisation residual was NOT computed.
  - D-3. The main article was never acquired, and that is the binding constraint.
  - D-4. `Experiment_ID` is proposed, not reserved.
  - D-5. Table S1's four molar figures do not sum to 100 on 75 of 93 rows.
  - D-6. Licence is not release-compatible.
  - D-7. `data/` carries a conflicting structure for `C12-200`, and we do not edit `data/`.
  - D-8. Three `Model_target` values would be new to the corpus.

## 08_contradictions.csv — 9 rows
