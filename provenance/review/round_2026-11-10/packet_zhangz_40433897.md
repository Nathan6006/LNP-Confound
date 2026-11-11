# Review packet — `zhangz_40433897`

Generated from live files. Queue says **75** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 40433897.

## Rows read from `curation/zhangz_40433897/06_rows.csv` — 75 rows

**18 columns.** Corpus schema is 79; this file has 13 of them.

## Structures
- distinct full InChIKeys: **3**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **0** of 3
- **NEW to the corpus: 3**

## Values
- `Experiment_value_raw`: n=75, mean **53.5305**, sd(pop) **22.8503**, min 1.207, max 92.57
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Model_type`: {'4T1-Luc': 66, 'MC38-Luc': 9}
- `Cargo_type`: {'FLuc': 75}
- `Experiment_method`: {'siLuc_knockdown_normalized': 75}
- `value_normalization_unit`: {'experiment': 75}

## Quality columns (schema_quality.json)
- `confidence`: 75/75 populated — {'low': 75}
- `positional_ambiguity`: 0/75 populated — {'': 75}
- `structure_assumed`: 75/75 populated — {'1': 75}
- `value_is_bound`: 0/75 populated — {'': 75}
- `curation_note`: 75/75 populated — {'Structure PROPOSED, not pinned: the paper draws no VNxCyOz structure anywhere, so this is inferred from the Fig S2 HRMS formula (C70H121N2O6, core C6H16N2) plus the CAS list drawn in Fig S1a - Exactly one drawn amine has formula C6H16N2 (CAS 100-36-7, N,N-diethylethylenediamine) and exactly one of its nitrogens carries two hydrogens, so the two arms have only one site. Its +CH2 partner 104-78-9 is also drawn, which is what the measured CH2 ladder predicts. Our build reproduces the printed [M+H]+ calc to -0.08 mDa. Value digitised from the FigS4b embedded raster at 219.9 dpi (tier 5: no Source Data, no numeric value table, and zero vector paths on any figure page).': 48, "Structure PROPOSED, not pinned: the paper draws no VNxCyOz structure anywhere, so this is inferred from the Fig S2 HRMS formula (C70H121N2O6, core C6H16N2) plus the CAS list drawn in Fig S1a - Exactly one drawn amine has formula C6H16N2 (CAS 100-36-7, N,N-diethylethylenediamine) and exactly one of its nitrogens carries two hydrogens, so the two arms have only one site. Its +CH2 partner 104-78-9 is also drawn, which is what the measured CH2 ladder predicts. Our build reproduces the printed [M+H]+ calc to -0.08 mDa. Value digitised from the FigS4b embedded raster at 219.9 dpi (tier 5: no Source Data, no numeric value table, and zero vector paths on any figure page). Bar top accepted at a 2-row rule, that bar's own measured rule thickness.": 10, "Structure PROPOSED, not pinned: the paper draws no VNxCyOz structure anywhere, so this is inferred from the Fig S2 HRMS formula (C71H123N2O6, core C7H18N2) plus the CAS list drawn in Fig S1a - TWO drawn amines have formula C7H18N2: 104-78-9 (N,N-diethyl-1,3-propanediamine) and 10061-68-4 (N,N'-diethyl-1,3-propanediamine). The CH2 ladder discriminates: VN3C2 = VN2C2 + CH2 to 0.05 mDa, and the linker homologue of VN2C2's core (100-36-7) is 104-78-9, not 10061-68-4. Our build reproduces the printed [M+H]+ calc to -0.03 mDa. Value digitised from the FigS5a embedded raster at 219.5 dpi (tier 5: no Source Data, no numeric value table, and zero vector paths on any figure page).": 3, 'Structure PROPOSED, not pinned: the paper draws no VNxCyOz structure anywhere, so this is inferred from the Fig S2 HRMS formula (C69H119N2O6, core C5H14N2) plus the CAS list drawn in Fig S1a - Exactly one drawn amine has formula C5H14N2 (CAS 109-55-7, N,N-dimethyl-1,3-propanediamine) and exactly one of its nitrogens carries two hydrogens. Our build reproduces the printed [M+H]+ calc to -0.03 mDa. Value digitised from the FigS5a embedded raster at 219.5 dpi (tier 5: no Source Data, no numeric value table, and zero vector paths on any figure page).': 3, 'Structure PROPOSED, not pinned: the paper draws no VNxCyOz structure anywhere, so this is inferred from the Fig S2 HRMS formula (C70H121N2O6, core C6H16N2) plus the CAS list drawn in Fig S1a - Exactly one drawn amine has formula C6H16N2 (CAS 100-36-7, N,N-diethylethylenediamine) and exactly one of its nitrogens carries two hydrogens, so the two arms have only one site. Its +CH2 partner 104-78-9 is also drawn, which is what the measured CH2 ladder predicts. Our build reproduces the printed [M+H]+ calc to -0.08 mDa. Value digitised from the FigS5b embedded raster at 219.5 dpi (tier 5: no Source Data, no numeric value table, and zero vector paths on any figure page).': 3}

## Folder inventory
  - `curation/zhangz_40433897/00_acquisition.json` (1563 B)
  - `curation/zhangz_40433897/00_value_sources.json` (13435 B)
  - `curation/zhangz_40433897/01_amines.json` (16039 B)
  - `curation/zhangz_40433897/01_blocks.csv` (1351 B)
  - `curation/zhangz_40433897/02_bars.json` (49166 B)
  - `curation/zhangz_40433897/02_corpus_overlap.json` (4915 B)
  - `curation/zhangz_40433897/03_enumerated.csv` (2863 B)
  - `curation/zhangz_40433897/03_proposed_structures.json` (5164 B)
  - `curation/zhangz_40433897/05_counts.json` (281 B)
  - `curation/zhangz_40433897/05_values.csv` (48727 B)
  - `curation/zhangz_40433897/06_rows.csv` (72155 B)
  - `curation/zhangz_40433897/07_qc.json` (21572 B)
  - `curation/zhangz_40433897/08_contradictions.csv` (5053 B)
  - `curation/zhangz_40433897/DECISIONS.md` (8299 B)
  - `curation/zhangz_40433897/PROVENANCE.md` (11587 B)
  - `curation/zhangz_40433897/make_reference.py` (9842 B)
  - `curation/zhangz_40433897/reference_structures.csv` (4679 B)
  - `curation/zhangz_40433897/stage01_amines.py` (12923 B)
  - `curation/zhangz_40433897/stage02_bars.py` (23959 B)
  - `curation/zhangz_40433897/stage03_checks.py` (16575 B)
  - `curation/zhangz_40433897/stage03b_structures.py` (10284 B)
  - `curation/zhangz_40433897/stage04_overlap.py` (7237 B)
  - `curation/zhangz_40433897/stage05_emit.py` (16243 B)
  - `curation/zhangz_40433897/stage06_seeds.py` (14569 B)

## DECISIONS.md — **5 OPEN**
  - D-1 — No VNxCyOz structure in this paper is PINNED. Three are PROPOSED. Do they go into the corpus?
  - D-2 — Figure S4b's printed DMG-PEG correlation has the wrong sign. Report it, or treat our digitisation as suspect?
  - D-3 — Ten of 118 bars are refused. Accept the gaps, or find another route?
  - D-4 — Figure S4b's 60 formulations and Figure S5's 18 are different rounds. One block or six?
  - D-5 — `Dose`, `timepoint` and all formulation mol% for Figures S5, S7, S8 and S13 are not recoverable from the SI.

## 08_contradictions.csv — 7 rows
