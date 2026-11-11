# Review packet — `mckinlaycj_39105751`

Generated from live files. Queue says **5** emitted rows; status `extracted_below_threshold`; proposed `Experiment_ID` **none**; PMID 39105751.

## Rows read from `curation/mckinlaycj_39105751/05_values.csv` — 5 rows

**38 columns.** Corpus schema is 79; this file has 9 of them.

## Structures
- distinct full InChIKeys: **5**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **0** of 5
- **NEW to the corpus: 5**

## Values
- `Experiment_value_z`: n=5, mean **-1.11022e-17**, sd(pop) **1**, min -1.95, max 0.8455
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `value_normalization_unit`: {'experiment': 5}

## Quality columns (schema_quality.json)
- `confidence`: 5/5 populated — {'high': 5}
- `positional_ambiguity`: 5/5 populated — {'0': 5}
- `structure_assumed`: 5/5 populated — {'0': 5}
- `value_is_bound`: 0/5 populated — {'': 5}
- `curation_note`: 5/5 populated — {"Ehx (2-ethylhexyl) is a racemic centre the paper draws without a wedge, stored flat; value is Table 1's luciferase expression (p/s); the organ is not stated in Table 1 (liver -- the text cites 'Table 1 and Figure 5D' together and Figure 5D is captioned 'Total Fluc bioluminescence in balb/c liver').": 3, "value is Table 1's luciferase expression (p/s); the organ is not stated in Table 1 (liver -- the text cites 'Table 1 and Figure 5D' together and Figure 5D is captioned 'Total Fluc bioluminescence in balb/c liver').": 1, "structure drawn in full in SI Figure S7A and matched monomer-for-monomer; Ehx (2-ethylhexyl) is a racemic centre the paper draws without a wedge, stored flat; value is Table 1's luciferase expression (p/s); the organ is not stated in Table 1 (liver -- the text cites 'Table 1 and Figure 5D' together and Figure 5D is captioned 'Total Fluc bioluminescence in balb/c liver').": 1}

## Folder inventory
  - `curation/mckinlaycj_39105751/00_acquisition.json` (2968 B)
  - `curation/mckinlaycj_39105751/00_value_sources.json` (13456 B)
  - `curation/mckinlaycj_39105751/01_blocks.csv` (6516 B)
  - `curation/mckinlaycj_39105751/01_blocks_qc.json` (8502 B)
  - `curation/mckinlaycj_39105751/03_enumerated.csv` (8956 B)
  - `curation/mckinlaycj_39105751/05_values.csv` (8298 B)
  - `curation/mckinlaycj_39105751/07_qc.json` (7829 B)
  - `curation/mckinlaycj_39105751/08_contradictions.csv` (6903 B)
  - `curation/mckinlaycj_39105751/DECISIONS.md` (9456 B)
  - `curation/mckinlaycj_39105751/PROVENANCE.md` (6453 B)
  - `curation/mckinlaycj_39105751/blocks.py` (10747 B)
  - `curation/mckinlaycj_39105751/make_reference.py` (10613 B)
  - `curation/mckinlaycj_39105751/reference_structures.csv` (2640 B)
  - `curation/mckinlaycj_39105751/seqtables.py` (8166 B)
  - `curation/mckinlaycj_39105751/stage00_scan_for_values.py` (3259 B)
  - `curation/mckinlaycj_39105751/stage01_blocks.py` (7271 B)
  - `curation/mckinlaycj_39105751/stage03_enumerate.py` (7695 B)
  - `curation/mckinlaycj_39105751/stage05_values.py` (12132 B)
  - `curation/mckinlaycj_39105751/stage07_qc.py` (21114 B)

## DECISIONS.md — **5 OPEN**
  - D-1. Does a peptoid belong in a corpus of small-molecule ionizable lipids?
  - D-2. Is a 5-row experiment admissible, and against what is it z-scored?
  - D-3. Should SI Figure S10 (236 vs SM-102 vs MC3) be digitised?
  - D-4. Source licences were not captured.
  - D-5. `Model_target` for the Table 1 values: liver, or whole animal?

## 08_contradictions.csv — 11 rows
