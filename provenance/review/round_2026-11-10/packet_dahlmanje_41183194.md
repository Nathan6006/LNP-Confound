# Review packet — `dahlmanje_41183194`

Generated from live files. Queue says **10** emitted rows; status `extracted`; proposed `Experiment_ID` **KG_2025**; PMID 41183194.

## Rows read from `curation/dahlmanje_41183194/06_rows.csv` — 10 rows

**34 columns.** Corpus schema is 79; this file has 34 of them.

## Structures
- distinct full InChIKeys: **2**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **2** of 2
  | InChIKey | rows here | in corpus as |
  |---|---|---|
  | `BGNVBNJYBVCBJH-UHFFFAOYSA-N` | 6 | AC_2025:SM-102, SB_2024:SM-102, SS_2018:SM-102 |
  | `ZISVTYVLWSZJAL-UHFFFAOYSA-N` | 4 | AC_2025:CKK-E12, BL_2023:MD-1, JW_2024:other |
- **NEW to the corpus: 0**

## Values
- `Experiment_value`: n=10, mean **2e-07**, sd(pop) **1**, min -2.401, max 0.9829
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'KG_2025': 10}
- `Model_type`: {'Mouse_B6': 10}
- `Model_target`: {'spleen': 5, 'liver': 5}
- `Cargo_type`: {'Thy1.1': 10}
- `Experiment_method`: {'protein_adundance_normalized': 10}
- `Experiment_batching`: {'individual': 10}
- `value_normalization_unit`: {'experiment': 10}

## Quality columns (schema_quality.json)
- `confidence`: 10/10 populated — {'medium': 9, 'high': 1}
- `positional_ambiguity`: 10/10 populated — {'0': 10}
- `structure_assumed`: 10/10 populated — {'1': 8, '0': 2}
- `value_is_bound`: 0/10 populated — {'': 10}
- `curation_note`: 10/10 populated — {"value digitised from Fig. 2B (no numeral is printed); the same panel's GlycoLNP61 bars reproduce Fig. 3C's printed 54 and 78 to +0.67 and -0.16 points; the helper lipid's sphingoid base length is taken from the reagent name, not measured: the drawing gives 14 bonds by mid-line crossing and 13 by apex counting, though its N-acyl chain measures 16 carbons by all three methods.": 4, "value is a numeral printed inside Fig. 3C; the helper lipid's sphingoid base length is taken from the reagent name, not measured: the drawing gives 14 bonds by mid-line crossing and 13 by apex counting, though its N-acyl chain measures 16 carbons by all three methods.": 2, "value digitised from Fig. 2B (no numeral is printed); the same panel's GlycoLNP61 bars reproduce Fig. 3C's printed 54 and 78 to +0.67 and -0.16 points; the helper lipid's sphingoid base length is taken from the reagent name, not measured: the drawing gives 14 bonds by mid-line crossing and 13 by apex counting, though its N-acyl chain measures 16 carbons by all three methods; the paper prints 'DMG-PEG' with no PEG block length; the 2000 Da form is assumed.": 2, 'value is a numeral printed inside Fig. 3C.': 1, 'value is a numeral printed inside Fig. 3C; Fig. 3C prints 81 for this value and the Results text says 82; the figure numeral is emitted and the conflict is recorded in 08_contradictions.csv.': 1}

## Folder inventory
  - `curation/dahlmanje_41183194/00_acquisition.json` (6276 B)
  - `curation/dahlmanje_41183194/00_value_sources.json` (11347 B)
  - `curation/dahlmanje_41183194/01_blocks.csv` (2126 B)
  - `curation/dahlmanje_41183194/02_corpus_overlap.json` (8992 B)
  - `curation/dahlmanje_41183194/03_enumerated.csv` (798 B)
  - `curation/dahlmanje_41183194/03_seedgate.json` (4852 B)
  - `curation/dahlmanje_41183194/04_chains.json` (7550 B)
  - `curation/dahlmanje_41183194/04_fig2b.json` (2857 B)
  - `curation/dahlmanje_41183194/05_emit_checks.json` (939 B)
  - `curation/dahlmanje_41183194/05_values.csv` (7342 B)
  - `curation/dahlmanje_41183194/06_rows.csv` (10602 B)
  - `curation/dahlmanje_41183194/07_fig3d.json` (1629 B)
  - `curation/dahlmanje_41183194/07_qc.json` (15939 B)
  - `curation/dahlmanje_41183194/08_contradictions.csv` (5036 B)
  - `curation/dahlmanje_41183194/DECISIONS.md` (13621 B)
  - `curation/dahlmanje_41183194/PROVENANCE.md` (10108 B)
  - `curation/dahlmanje_41183194/make_reference.py` (18937 B)
  - `curation/dahlmanje_41183194/reference_structures.csv` (2741 B)
  - `curation/dahlmanje_41183194/stage01_corpus_overlap.py` (11467 B)
  - `curation/dahlmanje_41183194/stage02_fig2b.py` (15103 B)
  - `curation/dahlmanje_41183194/stage03_seedgate.py` (9163 B)
  - `curation/dahlmanje_41183194/stage04_chains.py` (14306 B)
  - `curation/dahlmanje_41183194/stage05_emit.py` (23350 B)
  - `curation/dahlmanje_41183194/stage06_checks.py` (20295 B)
  - `curation/dahlmanje_41183194/stage07_fig3d.py` (6135 B)

## DECISIONS.md — **10 OPEN**
  - D-1. Does a source whose only chemical variable is the helper lipid belong in this corpus?
  - D-2. The pooled 109-LNP screen publishes ratios and unlabelled points, and no attributable per-LNP value.
  - D-3. Six of the ten emitted rows rest on a pixel read of Fig. 2B.
  - D-4. Fig. 3A/3C and Fig. 3D disagree, in opposite directions, on the same condition.
  - D-5. `Control LNP61`'s liver value has two printed candidates, 81 and 82.
  - D-6. The glycolipid's sphingoid base length is not pinned by the raster.
  - D-7. What `KG_2025` would actually add to the in-vivo shared-lipid set — measured, not assumed.
  - D-8. The lung, kidney and heart values are printed but not emitted.
  - D-9. The article licence fails the §8 release requirement.
  - D-10. `GlycoLNP60` has a complete formulation and no delivery value.
  - D-11. Two papers from one laboratory are two `Experiment_ID`s but one laboratory.
  - D-12. `KG_2025` is the sixth instance of F-067's one-lipid formulation screen, and it is already curated the way F-067 recommends.

## 08_contradictions.csv — 6 rows
