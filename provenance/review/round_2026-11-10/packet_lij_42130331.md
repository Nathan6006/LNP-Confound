# Review packet — `lij_42130331`

Generated from live files. Queue says **300** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 42130331.

## Rows read from `curation/lij_42130331/07_QG_2026_rows.csv` — 300 rows

**79 columns.** Corpus schema is 79; this file has 79 of them.

## Structures
- distinct full InChIKeys: **300**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **0** of 300
- **NEW to the corpus: 300**

## Values
- `Experiment_value`: n=300, mean **-2e-08**, sd(pop) **1**, min -2.267, max 2.425
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'QG_2026': 300}
- `Model_type`: {'': 300}
- `Model_target`: {'': 300}
- `Cargo_type`: {'': 300}
- `Experiment_method`: {'': 300}
- `Experiment_batching`: {'': 300}
- `value_normalization_unit`: {'experiment': 300}

## Quality columns (schema_quality.json)
- `confidence`: 300/300 populated — {'medium': 300}
- `positional_ambiguity`: 300/300 populated — {'0': 300}
- `structure_assumed`: 300/300 populated — {'0': 300}
- `value_is_bound`: 0/300 populated — {'': 300}
- `curation_note`: 300/300 populated — {"beta-hydroxy thioether lipid enumerated from the Fig. S1A scheme and the Fig. S1B / Fig. S2 block drawings; the paper publishes no HRMS, NMR table, SMILES or structure file, so the thiolactone ring-opening regiochemistry and the thiol-epoxide regiochemistry are read from the drawing and have no independent check. Value is the authors' own printed number from SI Table S2 (Fig. 1E); the assay behind that figure is not stated anywhere in the SI and the main text is paywalled -- see DECISIONS.md D-1.": 250, "beta-hydroxy thioether lipid enumerated from the Fig. S1A scheme and the Fig. S1B / Fig. S2 block drawings; the paper publishes no HRMS, NMR table, SMILES or structure file, so the thiolactone ring-opening regiochemistry and the thiol-epoxide regiochemistry are read from the drawing and have no independent check. Value is the authors' own printed number from SI Table S1 (Fig. 1C); the assay behind that figure is not stated anywhere in the SI and the main text is paywalled -- see DECISIONS.md D-1.": 50}

## Folder inventory
  - `curation/lij_42130331/00_acquisition.json` (3832 B)
  - `curation/lij_42130331/00_pages/glyph_montage.png` (9985 B)
  - `curation/lij_42130331/00_pages/si_text.txt` (20674 B)
  - `curation/lij_42130331/00_pages/trace_figS1B_row1.png` (31017 B)
  - `curation/lij_42130331/00_pages/trace_figS1B_row2.png` (20509 B)
  - `curation/lij_42130331/00_pages/trace_figS1B_row3.png` (23938 B)
  - `curation/lij_42130331/00_pages/trace_figS2_legend.png` (31548 B)
  - `curation/lij_42130331/00_pages/zoom_R_C8.png` (5775 B)
  - `curation/lij_42130331/00_pages/zoom_R_OG8.png` (8907 B)
  - `curation/lij_42130331/00_pages/zoom_acylchloride.png` (23457 B)
  - `curation/lij_42130331/00_pages/zoom_amines_70_75_80.png` (20146 B)
  - `curation/lij_42130331/00_pages/zoom_amines_96_400.png` (21968 B)
  - `curation/lij_42130331/00_value_sources.json` (7213 B)
  - `curation/lij_42130331/01_blocks.csv` (2623 B)
  - `curation/lij_42130331/02_traces.json` (25636 B)
  - `curation/lij_42130331/03_enumerated.csv` (128428 B)
  - `curation/lij_42130331/03_enumeration_summary.json` (115 B)
  - `curation/lij_42130331/04_glyph_clusters.json` (6071 B)
  - `curation/lij_42130331/04_glyph_labels.json` (7774 B)
  - `curation/lij_42130331/04_glyphs.json` (9498 B)
  - `curation/lij_42130331/04_merged_glyph_verification.json` (5379 B)
  - `curation/lij_42130331/05_value_checks.json` (541 B)
  - `curation/lij_42130331/05_values.csv` (24827 B)
  - `curation/lij_42130331/07_QG_2026_rows.csv` (240934 B)
  - `curation/lij_42130331/07_emit_summary.json` (543 B)
  - `curation/lij_42130331/07_qc.json` (15849 B)
  - `curation/lij_42130331/08_contradictions.csv` (4122 B)
  - `curation/lij_42130331/DECISIONS.md` (9027 B)
  - `curation/lij_42130331/PROVENANCE.md` (9111 B)
  - `curation/lij_42130331/make_reference.py` (5878 B)
  - `curation/lij_42130331/reference_structures.csv` (1014 B)
  - `curation/lij_42130331/stage02_blocks.py` (18849 B)
  - `curation/lij_42130331/stage03_enumerate.py` (12218 B)
  - `curation/lij_42130331/stage04_tables.py` (16741 B)
  - `curation/lij_42130331/stage05_values.py` (6277 B)
  - `curation/lij_42130331/stage07_emit.py` (6719 B)
  - `curation/lij_42130331/stage09_qc.py` (17022 B)

## DECISIONS.md — **0 OPEN**
  - D-n  <one line>
  - D-1  What assay produced Fig. 1C and Fig. 1E -- i.e. what the 300 values ARE
  - D-2  Are Library 1 and Library 2 one experiment or two?
  - D-3  What formulation was the screen run at?
  - D-4  Does a library whose connectivity has no independent check belong in the corpus?
  - D-5  The source is not licensed for redistribution

## 08_contradictions.csv — 6 rows
