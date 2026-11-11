# Review packet — `goulart_41138890`

Generated from live files. Queue says **3** emitted rows; status `extracted`; proposed `Experiment_ID` **MR_2025**; PMID 41138890.

## Rows read from `curation/goulart_41138890/06_rows.csv` — 3 rows

**79 columns.** Corpus schema is 79; this file has 79 of them.

## Structures
- distinct full InChIKeys: **3**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **3** of 3
  | InChIKey | rows here | in corpus as |
  |---|---|---|
  | `RVHYPUORVDKRTM-UHFFFAOYSA-N` | 1 | AC_2025:C12-200, JW_2024:other |
  | `BGNVBNJYBVCBJH-UHFFFAOYSA-N` | 1 | AC_2025:SM-102, SB_2024:SM-102, SS_2018:SM-102 |
  | `QGWBEETXHOVFQS-UHFFFAOYSA-N` | 1 | AC_2025:ALC-0315, LZ_2025:ALC-0315, XD_2026:ALC-0315 |
- **NEW to the corpus: 0**

## Values
- `Experiment_value`: n=3, mean **3.33333e-11**, sd(pop) **1**, min -1.054, max 1.343
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'MR_2025': 3}
- `Model_type`: {'primary_rat_cortical_neuron': 3}
- `Model_target`: {'in_vitro': 3}
- `Cargo_type`: {'GFP': 3}
- `Experiment_method`: {'GFP_fluorescence_normalized': 3}
- `Experiment_batching`: {'': 3}
- `value_normalization_unit`: {'experiment': 3}

## Quality columns (schema_quality.json)
- `confidence`: 3/3 populated — {'high': 3}
- `positional_ambiguity`: 3/3 populated — {'0': 3}
- `structure_assumed`: 3/3 populated — {'0': 3}
- `value_is_bound`: 0/3 populated — {'': 3}
- `curation_note`: 3/3 populated — {'Structure NOT drawn in the paper -- the name alone is given in Fig 4G; pinned from provenance/reagent_registry.py. Fig 4G bar digitised from the publisher raster (no SI table, no Source Data, no repo, article paywalled); PEG_SMILES left blank because the paper DRAWS a phosphoethanolamine-carbamate-PEG45 while the corpus name C14-PEG2000 stores a glycerol ether; formulation is the one optimised for C12-200, so the two benchmark arms are not formulation-matched to their own optimum.': 2, 'Structure drawn in full in Fig 1A and taken from the corpus (AC_2025, 498 rows, one structure); an anchor assembled independently from that drawing matches it at an exact full InChIKey. Fig 4G bar digitised from the publisher raster (no SI table, no Source Data, no repo, article paywalled); PEG_SMILES left blank because the paper DRAWS a phosphoethanolamine-carbamate-PEG45 while the corpus name C14-PEG2000 stores a glycerol ether; formulation is the one optimised for C12-200, so the two benchmark arms are not formulation-matched to their own optimum.': 1}

## Folder inventory
  - `curation/goulart_41138890/00_acquisition_raw.json` (263 B)
  - `curation/goulart_41138890/00_value_sources.json` (9903 B)
  - `curation/goulart_41138890/01_blocks.csv` (1057 B)
  - `curation/goulart_41138890/02_corpus_overlap.json` (9347 B)
  - `curation/goulart_41138890/02_si_tables.json` (8166 B)
  - `curation/goulart_41138890/02_tableS1.csv` (954 B)
  - `curation/goulart_41138890/02_tableS3.csv` (682 B)
  - `curation/goulart_41138890/03_crosscheck.json` (924 B)
  - `curation/goulart_41138890/03_enumerated.csv` (854 B)
  - `curation/goulart_41138890/04_fig4g.json` (4510 B)
  - `curation/goulart_41138890/05_anchor_check.json` (1960 B)
  - `curation/goulart_41138890/05_values.csv` (5155 B)
  - `curation/goulart_41138890/06_emit_qc.json` (707 B)
  - `curation/goulart_41138890/06_rows.csv` (4568 B)
  - `curation/goulart_41138890/07_qc.json` (11124 B)
  - `curation/goulart_41138890/08_contradictions.csv` (4139 B)
  - `curation/goulart_41138890/09_seeds.json` (5730 B)
  - `curation/goulart_41138890/DECISIONS.md` (7656 B)
  - `curation/goulart_41138890/PROVENANCE.md` (9488 B)
  - `curation/goulart_41138890/make_reference.py` (12124 B)
  - `curation/goulart_41138890/reference_structures.csv` (2010 B)
  - `curation/goulart_41138890/stage02_si_tables.py` (11886 B)
  - `curation/goulart_41138890/stage03_crosscheck.py` (4088 B)
  - `curation/goulart_41138890/stage04_fig4g.py` (14360 B)
  - `curation/goulart_41138890/stage05_anchor_check.py` (8217 B)
  - `curation/goulart_41138890/stage06_emit.py` (14927 B)
  - `curation/goulart_41138890/stage07_qc.py` (22296 B)
  - `curation/goulart_41138890/stage09_seeds.py` (11395 B)

## DECISIONS.md — **5 OPEN**
  - D-1 — Should `MR_2025` be merged at all? All three of its ionizable lipids are universal controls already in the corpus.
  - D-2 — `Experiment_method = GFP_fluorescence_normalized` is a new token in a three-value vocabulary.
  - D-3 — `Model_type = primary_rat_cortical_neuron`: a new cell type, and the second rat cell in the corpus.
  - D-4 — `PEG_SMILES` is blank because the corpus name resolves to the wrong molecule.
  - D-5 — `Dose_ug_nucleicacid` left blank for Figure 4G.
  - D-6 — The three-lipid ranking is confounded by formulation.
  - D-7 — Only the 72 h timepoint is emitted, of three measured.
  - D-8 — The in-vivo and in-vitro arms are NOT pooled, and the in-vivo arm is not emitted.

## 08_contradictions.csv — 8 rows
