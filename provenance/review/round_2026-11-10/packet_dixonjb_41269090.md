# Review packet — `dixonjb_41269090`

Generated from live files. Queue says **28** emitted rows; status `extracted`; proposed `Experiment_ID` **EM_2025**; PMID 41269090.

## Rows read from `curation/dixonjb_41269090/06_rows.csv` — 28 rows

**79 columns.** Corpus schema is 79; this file has 79 of them.

## Structures
- distinct full InChIKeys: **2**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **2** of 2
  | InChIKey | rows here | in corpus as |
  |---|---|---|
  | `ZISVTYVLWSZJAL-UHFFFAOYSA-N` | 23 | AC_2025:CKK-E12, BL_2023:MD-1, JW_2024:other |
  | `NRLNQCOGCKAESA-KWXKLSQISA-N` | 5 | AC_2025:DLin-MC3-DMA, BL_2023:DLin-MC3-DMA, NG_2025_2:DLin-MC3-DMA |
- **NEW to the corpus: 0**

## Values
- `Experiment_value`: n=28, mean **-3.57143e-08**, sd(pop) **1**, min -2.073, max 1.269
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'EM_2025': 28}
- `Model_type`: {'Mouse_B6': 28}
- `Model_target`: {'lymph_node_popliteal': 8, 'lymph_node_axillary': 6, 'lymph_node_brachial': 6, 'lymph_node_axillary_nondraining': 2, 'lymph_node_brachial_nondraining': 2, 'lymph_node_popliteal_contralateral_nondraining': 2}
- `Cargo_type`: {'aVHH': 28}
- `Experiment_method`: {'protein_adundance_normalized': 28}
- `Experiment_batching`: {'individual': 28}
- `value_normalization_unit`: {'experiment': 28}

## Quality columns (schema_quality.json)
- `confidence`: 28/28 populated — {'medium': 23, 'high': 5}
- `positional_ambiguity`: 28/28 populated — {'0': 28}
- `structure_assumed`: 28/28 populated — {'0': 28}
- `value_is_bound`: 0/28 populated — {'': 28}
- `curation_note`: 28/28 populated — {"Value is a RATIO OF MEANS derived from SI Table S2: the average aVHH+ cell count divided by the average LEC count. The paper's own figure plots a MEAN OF RATIOS (one point per independent experiment); the two statistics differ by up to 2.5x on one cell of this paper and only the ratio of means is derivable from published material. Structure is the corpus string already carried by AC_2025, JW_2024 and BL_2023 for this same molecule; it declares no stereochemistry, while cKK-E12 is lysine-derived and the commercial reagent (BroadPharm BP-29590) is a single diastereomer, so the two diketopiperazine centres are undeclared. confidence=medium records that, not a doubt about the constitution. Cohort: study_A (both paws injected, Figure 2 / Table S2).": 18, "Value is a RATIO OF MEANS derived from SI Table S3/S4: the average aVHH+ cell count divided by the average LEC count. The paper's own figure plots a MEAN OF RATIOS (one point per independent experiment); the two statistics differ by up to 2.5x on one cell of this paper and only the ratio of means is derivable from published material. Structure is the corpus string already carried by AC_2025, JW_2024 and BL_2023 for this same molecule; it declares no stereochemistry, while cKK-E12 is lysine-derived and the commercial reagent (BroadPharm BP-29590) is a single diastereomer, so the two diketopiperazine centres are undeclared. confidence=medium records that, not a doubt about the constitution. Cohort: study_B (one hindlimb injected, Figure 4 / Tables S3-S4). This node does NOT drain the injection site in that cohort, which Model_target records as _nondraining; it is the paper's own off-target control.": 3, "Value is a RATIO OF MEANS derived from SI Table S3/S4: the average aVHH+ cell count divided by the average LEC count. The paper's own figure plots a MEAN OF RATIOS (one point per independent experiment); the two statistics differ by up to 2.5x on one cell of this paper and only the ratio of means is derivable from published material. Benchmark structure taken from provenance/reagent_registry.py. Cohort: study_B (one hindlimb injected, Figure 4 / Tables S3-S4). This node does NOT drain the injection site in that cohort, which Model_target records as _nondraining; it is the paper's own off-target control.": 3, "Value is a RATIO OF MEANS derived from SI Table S3/S4: the average aVHH+ cell count divided by the average LEC count. The paper's own figure plots a MEAN OF RATIOS (one point per independent experiment); the two statistics differ by up to 2.5x on one cell of this paper and only the ratio of means is derivable from published material. Structure is the corpus string already carried by AC_2025, JW_2024 and BL_2023 for this same molecule; it declares no stereochemistry, while cKK-E12 is lysine-derived and the commercial reagent (BroadPharm BP-29590) is a single diastereomer, so the two diketopiperazine centres are undeclared. confidence=medium records that, not a doubt about the constitution. Cohort: study_B (one hindlimb injected, Figure 4 / Tables S3-S4).": 2, "Value is a RATIO OF MEANS derived from SI Table S3/S4: the average aVHH+ cell count divided by the average LEC count. The paper's own figure plots a MEAN OF RATIOS (one point per independent experiment); the two statistics differ by up to 2.5x on one cell of this paper and only the ratio of means is derivable from published material. Benchmark structure taken from provenance/reagent_registry.py. Cohort: study_B (one hindlimb injected, Figure 4 / Tables S3-S4).": 2}

## Folder inventory
  - `curation/dixonjb_41269090/00_acquisition.json` (6557 B)
  - `curation/dixonjb_41269090/00_value_sources.json` (10094 B)
  - `curation/dixonjb_41269090/01_blocks.csv` (2576 B)
  - `curation/dixonjb_41269090/01_screening1.csv` (1832 B)
  - `curation/dixonjb_41269090/02_corpus_overlap.json` (10108 B)
  - `curation/dixonjb_41269090/03_enumerated.csv` (5031 B)
  - `curation/dixonjb_41269090/03_enumerated_qc.json` (2751 B)
  - `curation/dixonjb_41269090/04_si_tables.json` (6706 B)
  - `curation/dixonjb_41269090/05_values.csv` (142655 B)
  - `curation/dixonjb_41269090/06_rows.csv` (46320 B)
  - `curation/dixonjb_41269090/06_rows_z.json` (128 B)
  - `curation/dixonjb_41269090/07_qc.json` (20921 B)
  - `curation/dixonjb_41269090/08_contradictions.csv` (5505 B)
  - `curation/dixonjb_41269090/DECISIONS.md` (6700 B)
  - `curation/dixonjb_41269090/PROVENANCE.md` (17418 B)
  - `curation/dixonjb_41269090/make_reference.py` (10903 B)
  - `curation/dixonjb_41269090/reference_structures.csv` (1713 B)
  - `curation/dixonjb_41269090/stage00_scan_for_values.py` (2796 B)
  - `curation/dixonjb_41269090/stage01_blocks.py` (6588 B)
  - `curation/dixonjb_41269090/stage02_corpus_overlap.py` (11377 B)
  - `curation/dixonjb_41269090/stage03_enumerate.py` (6710 B)
  - `curation/dixonjb_41269090/stage05_values.py` (11897 B)
  - `curation/dixonjb_41269090/stage07_emit.py` (13866 B)
  - `curation/dixonjb_41269090/stage08_qc.py` (31640 B)

## DECISIONS.md — **6 OPEN**
  - D-1 Does the pooled SANDS barcode screen belong in the corpus at all?
  - D-2 Does a screen with ONE ionizable lipid belong in a cross-laboratory chemistry corpus?
  - D-3 Should the eight non-LEC cell populations be emitted?
  - D-4 Is a ratio of means an acceptable `Experiment_value` where the paper's own statistic is a mean of ratios?
  - D-5 `Model_target` needs five new terms, and `Experiment_method` is a spelling the corpus got wrong
  - D-6 `Dose_ug_nucleicacid` is blank on all 28 rows

## 08_contradictions.csv — 6 rows
