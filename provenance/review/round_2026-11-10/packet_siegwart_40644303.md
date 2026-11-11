# Review packet — `siegwart_40644303`

Generated from live files. Queue says **505** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 40644303.

## Rows read from `curation/siegwart_40644303/09_rows_vitro.csv` — 465 rows  **(queue says 505 — DISCREPANCY)**

Other schema-valid files in this folder (a reviewer should know which one the merge would take):
  - `curation/siegwart_40644303/09_rows_vitro.csv` — 465 rows
  - `curation/siegwart_40644303/09_rows_vivo.csv` — 40 rows

**79 columns.** Corpus schema is 79; this file has 79 of them.

## Structures
- distinct full InChIKeys: **465**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **360** of 465
  | InChIKey | rows here | in corpus as |
  |---|---|---|
  | `NBRVOPVWEPRENW-UHFFFAOYSA-N` | 1 | KZ_2016:3A1-SC1 |
  | `WKKVESBHBQRSJW-UHFFFAOYSA-N` | 1 | KZ_2016:3A1-SC2 |
  | `BPKMHWVKTMJOPH-UHFFFAOYSA-N` | 1 | KZ_2016:3A1-SC4 |
  | `JLNUXINXEONXJX-UHFFFAOYSA-N` | 1 | KZ_2016:3A1-SC5 |
  | `BRVFIBJORLJFFM-UHFFFAOYSA-N` | 1 | KZ_2016:3A1-SC6 |
  | `VGJUPJYFPBQKGH-UHFFFAOYSA-N` | 1 | KZ_2016:3A1-SC7 |
  | `UWCGSIHPQSDOTN-UHFFFAOYSA-N` | 1 | KZ_2016:3A1-SC8 |
  | `AUYJBZCTIRQTGK-UHFFFAOYSA-N` | 1 | KZ_2016:3A1-SC9 |
- **NEW to the corpus: 105**

## Values
- `Experiment_value`: n=465, mean **-2.43011e-10**, sd(pop) **1**, min -0.6089, max 2.114
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'SL_2025': 465}
- `Model_type`: {'IGROV1': 465}
- `Model_target`: {'in_vitro': 465}
- `Cargo_type`: {'FLuc': 465}
- `Experiment_method`: {'luminescence_discretized_normalized': 465}
- `Experiment_batching`: {'individual': 465}
- `value_normalization_unit`: {'experiment': 465}

## Quality columns (schema_quality.json)
- `confidence`: 465/465 populated — {'high': 341, 'medium': 124}
- `positional_ambiguity`: 465/465 populated — {'0': 341, '1': 124}
- `structure_assumed`: 465/465 populated — {'0': 465}
- `value_is_bound`: 68/465 populated — {'': 397, 'lower': 68}
- `curation_note`: 465/465 populated — {'value is the ordinal of a three-level band read from Figure 1a, not a continuous measurement': 296, 'the core carries 5 arms on 6 N-H and the paper does not say which nitrogen keeps its hydrogen; the stored structure follows the stated convention (primary N-H first, then secondary from the terminus inward); value is the ordinal of a three-level band read from Figure 1a, not a continuous measurement': 54, "Figure 1a is a three-level heat map; this cell's band is the open top band '> log2', so the value is a LOWER bound.": 45, 'the core carries 6 arms on 7 N-H and the paper does not say which nitrogen keeps its hydrogen; the stored structure follows the stated convention (primary N-H first, then secondary from the terminus inward); value is the ordinal of a three-level band read from Figure 1a, not a continuous measurement': 29, 'the core carries 5 arms on 7 N-H and the paper does not say which nitrogen keeps its hydrogen; the stored structure follows the stated convention (primary N-H first, then secondary from the terminus inward); value is the ordinal of a three-level band read from Figure 1a, not a continuous measurement': 18}

## Folder inventory
  - `curation/siegwart_40644303/00_acquisition.json` (7439 B)
  - `curation/siegwart_40644303/00_value_sources.json` (9403 B)
  - `curation/siegwart_40644303/01_blocks.csv` (5419 B)
  - `curation/siegwart_40644303/01_blocks_qc.json` (5435 B)
  - `curation/siegwart_40644303/02_si_table_S3.csv` (1341 B)
  - `curation/siegwart_40644303/03_enumerate_qc.json` (774 B)
  - `curation/siegwart_40644303/03_enumerated.csv` (146176 B)
  - `curation/siegwart_40644303/04_fig1_heatmap.json` (194294 B)
  - `curation/siegwart_40644303/04_invivo_bars.json` (14410 B)
  - `curation/siegwart_40644303/05_values.csv` (295845 B)
  - `curation/siegwart_40644303/06_scatter_check.json` (17475 B)
  - `curation/siegwart_40644303/07_qc.json` (20919 B)
  - `curation/siegwart_40644303/08_contradictions.csv` (6526 B)
  - `curation/siegwart_40644303/09_emit_qc.json` (651 B)
  - `curation/siegwart_40644303/09_rows_vitro.csv` (318445 B)
  - `curation/siegwart_40644303/09_rows_vivo.csv` (34266 B)
  - `curation/siegwart_40644303/DECISIONS.md` (11887 B)
  - `curation/siegwart_40644303/PROVENANCE.md` (4401 B)
  - `curation/siegwart_40644303/make_reference.py` (7980 B)
  - `curation/siegwart_40644303/reference_structures.csv` (1468 B)
  - `curation/siegwart_40644303/stage01_blocks.py` (14100 B)
  - `curation/siegwart_40644303/stage02_si_table.py` (4299 B)
  - `curation/siegwart_40644303/stage03_enumerate.py` (8655 B)
  - `curation/siegwart_40644303/stage04_fig1_heatmap.py` (21423 B)
  - `curation/siegwart_40644303/stage04_invivo_bars.py` (11431 B)
  - `curation/siegwart_40644303/stage05_values.py` (5085 B)
  - `curation/siegwart_40644303/stage06_scatter_check.py` (12628 B)
  - `curation/siegwart_40644303/stage07_emit.py` (11528 B)
  - `curation/siegwart_40644303/stage08_qc.py` (34160 B)

## DECISIONS.md — **10 OPEN**
  - D-n  <one line>
  - D-1  Are Figure 2b and Figure 3b one in-vivo experiment or two?
  - D-2  Should the siRNA half of Figure 1a be emitted at all?
  - D-3  Two in-vivo lipids have no readable value. Drop, or record a bound?
  - D-4  SC3: n-propyl (this paper's Figure 1b) or neopentyl (corpus KZ_2016)?
  - D-5  Is Figure 1a's lowest mRNA band closed at zero, or open below?
  - D-6  3A2: piperazine (Figure 1b) or the corpus's two-nitrogen core?
  - D-7  The paper is not CC0/CC-BY-compatible. Can it enter a released corpus?
  - D-8  The partial-substitution convention on 5A2 and 6A1
  - D-9  Is a three-level ordinal an acceptable `Experiment_value`?
  - D-10  This experiment shares all 465 of its lipids with KZ_2016. Is that a shared-lipid pair, or a duplicate?

## 08_contradictions.csv — 10 rows
