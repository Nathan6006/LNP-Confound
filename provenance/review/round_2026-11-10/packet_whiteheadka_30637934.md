# Review packet — `whiteheadka_30637934`

Generated from live files. Queue says **11** emitted rows; status `extracted`; proposed `Experiment_ID` **KW_2019**; PMID 30637934.

## Rows read from `curation/whiteheadka_30637934/06_rows.csv` — 11 rows

**74 columns.** Corpus schema is 79; this file has 74 of them.

## Structures
- distinct full InChIKeys: **11**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **6** of 11
  | InChIKey | rows here | in corpus as |
  |---|---|---|
  | `BPPZRZOKGADTQE-UHFFFAOYSA-N` | 1 | NC_2024:A306_Oi10 |
  | `PGKBCAOUJZLRIS-UHFFFAOYSA-N` | 1 | KW_2014:A306_C10, NC_2024:A306_O10 |
  | `WYYPKZMMNDPYQD-UHFFFAOYSA-N` | 1 | KW_2014:A306_C11, NC_2024:A306_O11 |
  | `ZBKHKDUXOVRTRE-UHFFFAOYSA-N` | 1 | KW_2014:A306_C12, NC_2024:A306_O12 |
  | `MSCOHYYLRYFQHL-UHFFFAOYSA-N` | 1 | KW_2014:A306_C13, NC_2024:A306_O13 |
  | `LRQXEZVLSSESNB-UHFFFAOYSA-N` | 1 | KW_2014:A306_C14, NC_2024:A306_O14 |
- **NEW to the corpus: 5**

## Values
- `Experiment_value`: n=11, mean **6.8632e-16**, sd(pop) **1**, min -2.104, max 1.067
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'KW_2019': 11}
- `Model_type`: {'HeLa': 11}
- `Model_target`: {'in_vitro': 11}
- `Cargo_type`: {'FLuc': 11}
- `Experiment_method`: {'luminescence_normalized': 11}
- `Experiment_batching`: {'individual': 11}
- `value_normalization_unit`: {'experiment': 11}

## Quality columns (schema_quality.json)
- `confidence`: 11/11 populated — {'high': 10, 'medium': 1}
- `positional_ambiguity`: 11/11 populated — {'0': 11}
- `structure_assumed`: 11/11 populated — {'1': 11}
- `value_is_bound`: 0/11 populated — {'': 11}
- `curation_note`: 11/11 populated — {"In-vitro values were measured on CRUDE 1:4 reaction product (SI Methods p.2: 'In vitro experiments were conducted with crude lipidoids'); the stored structure is the nominal tetra-substituted species, so structure_assumed=1. The 306 core's four N-H sites are two symmetric pairs, so the tetra-product has no regiochemical choice and positional_ambiguity is 0 -- the uncertainty is substitution DEGREE, not position. Excipient mol% read as helper 16 / cholesterol 46.5 against the SI sentence's positional order (cholesterol 16 / DOPE 46.5); see DECISIONS.md D-1. Value: Fig S6A bar top, tier 4b PDF vector coordinates (4,751 RLU, s.d. 1,186); the fill/stroke frame offset of 0.057 pt is worth 1,166 RLU, i.e. 24.5% of this value": 1, "In-vitro values were measured on CRUDE 1:4 reaction product (SI Methods p.2: 'In vitro experiments were conducted with crude lipidoids'); the stored structure is the nominal tetra-substituted species, so structure_assumed=1. The 306 core's four N-H sites are two symmetric pairs, so the tetra-product has no regiochemical choice and positional_ambiguity is 0 -- the uncertainty is substitution DEGREE, not position. Excipient mol% read as helper 16 / cholesterol 46.5 against the SI sentence's positional order (cholesterol 16 / DOPE 46.5); see DECISIONS.md D-1. Value: Fig S6A bar top, tier 4b PDF vector coordinates (741,646 RLU, s.d. 139,302); the fill/stroke frame offset of 0.057 pt is worth 1,166 RLU, i.e. 0.2% of this value": 1, "In-vitro values were measured on CRUDE 1:4 reaction product (SI Methods p.2: 'In vitro experiments were conducted with crude lipidoids'); the stored structure is the nominal tetra-substituted species, so structure_assumed=1. The 306 core's four N-H sites are two symmetric pairs, so the tetra-product has no regiochemical choice and positional_ambiguity is 0 -- the uncertainty is substitution DEGREE, not position. Excipient mol% read as helper 16 / cholesterol 46.5 against the SI sentence's positional order (cholesterol 16 / DOPE 46.5); see DECISIONS.md D-1. Value: Fig S6A bar top, tier 4b PDF vector coordinates (95,220 RLU, s.d. 19,051); the fill/stroke frame offset of 0.057 pt is worth 1,166 RLU, i.e. 1.2% of this value": 1, "Tail is ISODECYL acrylate, a commercial branched-isomer mixture; stored as 8-methylnonyl, matching the corpus's NC_2024 A306_Oi10 (full-InChIKey identical). The specific branch is not pinned by the paper. In-vitro values were measured on CRUDE 1:4 reaction product (SI Methods p.2: 'In vitro experiments were conducted with crude lipidoids'); the stored structure is the nominal tetra-substituted species, so structure_assumed=1. The 306 core's four N-H sites are two symmetric pairs, so the tetra-product has no regiochemical choice and positional_ambiguity is 0 -- the uncertainty is substitution DEGREE, not position. Excipient mol% read as helper 16 / cholesterol 46.5 against the SI sentence's positional order (cholesterol 16 / DOPE 46.5); see DECISIONS.md D-1. Value: Fig S6A bar top, tier 4b PDF vector coordinates (1,176,187 RLU, s.d. 138,096); the fill/stroke frame offset of 0.057 pt is worth 1,166 RLU, i.e. 0.1% of this value": 1, "In-vitro values were measured on CRUDE 1:4 reaction product (SI Methods p.2: 'In vitro experiments were conducted with crude lipidoids'); the stored structure is the nominal tetra-substituted species, so structure_assumed=1. The 306 core's four N-H sites are two symmetric pairs, so the tetra-product has no regiochemical choice and positional_ambiguity is 0 -- the uncertainty is substitution DEGREE, not position. Excipient mol% read as helper 16 / cholesterol 46.5 against the SI sentence's positional order (cholesterol 16 / DOPE 46.5); see DECISIONS.md D-1. Value: Fig S6A bar top, tier 4b PDF vector coordinates (1,741,649 RLU, s.d. 253,583); the fill/stroke frame offset of 0.057 pt is worth 1,166 RLU, i.e. 0.1% of this value": 1}

## Folder inventory
  - `curation/whiteheadka_30637934/00_acquisition.json` (3149 B)
  - `curation/whiteheadka_30637934/00_value_sources.json` (8592 B)
  - `curation/whiteheadka_30637934/01_blocks.csv` (1173 B)
  - `curation/whiteheadka_30637934/01_mass_identity.json` (391 B)
  - `curation/whiteheadka_30637934/02_corpus_overlap.json` (62362 B)
  - `curation/whiteheadka_30637934/02_fig_s6_geometry.json` (3721 B)
  - `curation/whiteheadka_30637934/02_fig_s6_panelA.csv` (930 B)
  - `curation/whiteheadka_30637934/02_fig_s6_panelB.csv` (454 B)
  - `curation/whiteheadka_30637934/03_enumerated.csv` (3411 B)
  - `curation/whiteheadka_30637934/04_table_s1.csv` (398 B)
  - `curation/whiteheadka_30637934/04_table_s2.csv` (414 B)
  - `curation/whiteheadka_30637934/04_tables_qc.json` (3338 B)
  - `curation/whiteheadka_30637934/05_values.csv` (22451 B)
  - `curation/whiteheadka_30637934/05_values_qc.json` (2605 B)
  - `curation/whiteheadka_30637934/06_rows.csv` (18461 B)
  - `curation/whiteheadka_30637934/06_rows_qc.json` (756 B)
  - `curation/whiteheadka_30637934/07_qc.json` (11595 B)
  - `curation/whiteheadka_30637934/07_seeded_checks.json` (5013 B)
  - `curation/whiteheadka_30637934/08_contradictions.csv` (10219 B)
  - `curation/whiteheadka_30637934/DECISIONS.md` (8706 B)
  - `curation/whiteheadka_30637934/PROVENANCE.md` (8221 B)
  - `curation/whiteheadka_30637934/make_reference.py` (10297 B)
  - `curation/whiteheadka_30637934/reference_structures.csv` (1821 B)
  - `curation/whiteheadka_30637934/stage01_enumerate.py` (10384 B)
  - `curation/whiteheadka_30637934/stage02_fig_s6.py` (15635 B)
  - `curation/whiteheadka_30637934/stage03_corpus_overlap.py` (9169 B)
  - `curation/whiteheadka_30637934/stage04_tables.py` (8511 B)
  - `curation/whiteheadka_30637934/stage05_values.py` (8355 B)
  - `curation/whiteheadka_30637934/stage06_rows.py` (13002 B)
  - `curation/whiteheadka_30637934/stage07_seeded_checks.py` (15003 B)
  - `curation/whiteheadka_30637934/stage08_qc.py` (14393 B)

## DECISIONS.md — **5 OPEN**
  - D-1  Which excipient takes 16 mol% and which takes 46.5 mol%
  - D-2  Whether the in-vivo panels are digitised and emitted at all
  - D-3  Whether Figure S6B becomes a separate Experiment_ID
  - D-4  The isodecyl branch in 306Oi10
  - D-5  Whether all 11 rows should carry structure_assumed = 1 for crude product

## 08_contradictions.csv — 10 rows
