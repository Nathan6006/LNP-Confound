# Review packet — `maohq_39375194`

Generated from live files. Queue says **5112** emitted rows; status `extracted`; proposed `Experiment_ID` **LC_2024**; PMID 39375194.

## Rows read from `curation/maohq_39375194/06_rows.csv` — 5112 rows

**79 columns.** Corpus schema is 79; this file has 79 of them.

## Structures
- distinct full InChIKeys: **1**
- structure column: **`IL_SMILES`**
- unparseable structures: **0**
- already in the corpus: **1** of 1
  | InChIKey | rows here | in corpus as |
  |---|---|---|
  | `NRLNQCOGCKAESA-KWXKLSQISA-N` | 5112 | AC_2025:DLin-MC3-DMA, BL_2023:DLin-MC3-DMA, NG_2025_2:DLin-MC3-DMA  **… and 8 more of 11 total** |
- **NEW to the corpus: 0**

## Values
- `Experiment_value`: n=5112, mean **-1.56534e-17**, sd(pop) **1**, min -2.019, max 3.065
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'LC_2024': 4680, 'LC_2024_2': 432}
- `Model_type`: {'HEK293T': 1152, 'B16F10': 1152, 'PC3': 1152, 'Neuro2a': 792, 'ARPE-19': 792, 'HepG2': 72}
- `Model_target`: {'in_vitro': 5112}
- `Cargo_type`: {'FLuc': 5112}
- `Experiment_method`: {'luminescence_normalized': 5112}
- `Experiment_batching`: {'individual': 5112}
- `value_normalization_unit`: {'experiment': 5112}

## Quality columns (schema_quality.json)
- `confidence`: 5112/5112 populated — {'high': 4140, 'medium': 972}
- `positional_ambiguity`: 5112/5112 populated — {'0': 5112}
- `structure_assumed`: 5112/5112 populated — {'0': 4140, '1': 972}
- `value_is_bound`: 0/5112 populated — {'': 5112}
- `curation_note`: 5112/5112 populated — {"Values are the authors' own ln(mean RLU) from github.com/MaoResearchGroup/ml_lnp_design_rules@82cee58 Raw_Data/, not digitised. Dose 0.1 ug pDNA/well and 48 h come from this paper's Methods, which describe the same screen that Nat Commun 2022 (PMID 35879315) describes at 0.2 ug/well and 24 h -- see 08_contradictions.csv C-1.": 3168, "Helper lipid '18PG' is emitted as DIOLEOYL PG (18:1/18:1), which is what SI Table 2 draws and what SI Table 3's Double_bonds = 2 states; the corpus registry pins the name '18PG' to SOPG (18:0/18:1) on the authority of PMID 35879315, which describes THIS SAME library. Unresolved: DECISIONS.md D-2. Values are the authors' own ln(mean RLU) from github.com/MaoResearchGroup/ml_lnp_design_rules@82cee58 Raw_Data/, not digitised. Dose 0.1 ug pDNA/well and 48 h come from this paper's Methods, which describe the same screen that Nat Commun 2022 (PMID 35879315) describes at 0.2 ug/well and 24 h -- see 08_contradictions.csv C-1.": 972, "DOTAP is emitted without a stereocentre, matching the reagent registry and 2,564 existing corpus rows; SI Table 2 draws a defined wedge at the propane C2, so the emitted stereochemistry is a convention, not a reading. See 08_contradictions.csv C-3. Values are the authors' own ln(mean RLU) from github.com/MaoResearchGroup/ml_lnp_design_rules@82cee58 Raw_Data/, not digitised. Dose 0.1 ug pDNA/well and 48 h come from this paper's Methods, which describe the same screen that Nat Commun 2022 (PMID 35879315) describes at 0.2 ug/well and 24 h -- see 08_contradictions.csv C-1.": 972}

## Folder inventory
  - `curation/maohq_39375194/00_value_sources.json` (11288 B)
  - `curation/maohq_39375194/01_blocks.csv` (3250 B)
  - `curation/maohq_39375194/01_checks.json` (12970 B)
  - `curation/maohq_39375194/01_seed_results.json` (5148 B)
  - `curation/maohq_39375194/02_corpus_overlap.json` (13311 B)
  - `curation/maohq_39375194/02b_duplicate_quantified.json` (989 B)
  - `curation/maohq_39375194/03_enumerated.csv` (1149 B)
  - `curation/maohq_39375194/03_structures_qc.json` (1502 B)
  - `curation/maohq_39375194/04_c10_anchor_check.json` (2637 B)
  - `curation/maohq_39375194/04_emit_seed_results.json` (1833 B)
  - `curation/maohq_39375194/05_value_counts.json` (1185 B)
  - `curation/maohq_39375194/05_values.csv` (4094214 B)
  - `curation/maohq_39375194/06_rows.csv` (6612001 B)
  - `curation/maohq_39375194/07_qc.json` (43979 B)
  - `curation/maohq_39375194/08_contradictions.csv` (8582 B)
  - `curation/maohq_39375194/DECISIONS.md` (9126 B)
  - `curation/maohq_39375194/PROVENANCE.md` (10899 B)
  - `curation/maohq_39375194/make_reference.py` (18665 B)
  - `curation/maohq_39375194/reference_structures.csv` (4627 B)
  - `curation/maohq_39375194/scripts/_extract_imgs.py` (714 B)
  - `curation/maohq_39375194/scripts/_jats.py` (443 B)
  - `curation/maohq_39375194/scripts/_peek.py` (370 B)
  - `curation/maohq_39375194/scripts/_render.py` (742 B)
  - `curation/maohq_39375194/scripts/_sitext.py` (245 B)
  - `curation/maohq_39375194/scripts/_xlsx.py` (379 B)
  - `curation/maohq_39375194/scripts/seed_checks.py` (9636 B)
  - `curation/maohq_39375194/scripts/seed_emit.py` (1811 B)
  - `curation/maohq_39375194/scripts/seed_emit_one.py` (1359 B)
  - `curation/maohq_39375194/scripts/stage01_invert.py` (26427 B)
  - `curation/maohq_39375194/scripts/stage02_overlap.py` (13911 B)
  - `curation/maohq_39375194/scripts/stage02b_duplicate.py` (3720 B)
  - `curation/maohq_39375194/scripts/stage03_structures.py` (14530 B)
  - `curation/maohq_39375194/scripts/stage03b_c10.py` (2729 B)
  - `curation/maohq_39375194/scripts/stage04_values.py` (18834 B)
  - `curation/maohq_39375194/scripts/stage05_qc.py` (12349 B)

## DECISIONS.md — **6 OPEN**
  - D-1 — Licence: neither the article nor the repository is CC0/CC-BY compatible
  - D-2 — `18PG` is SOPG in one paper from this lab and DOPG in the other, for the same library
  - D-3 — The dose and timepoint for this screen are stated twice, differently
  - D-4 — DOTAP: the paper draws a stereocentre the corpus does not store
  - D-5 — `RLU_HepG2` is already in the corpus; confirm the exclusion
  - D-6 — Is a 1,080-formulation screen of ONE ionizable lipid wanted in the corpus at all?

## 08_contradictions.csv — 10 rows
