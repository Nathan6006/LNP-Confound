# Review packet — `andersondg_26469188`

Generated from live files. Queue says **36** emitted rows; status `extracted`; proposed `Experiment_ID` **KK_2015**; PMID 26469188.

## Rows read from `curation/andersondg_26469188/06_rows.csv` — 36 rows

**79 columns.** Corpus schema is 79; this file has 79 of them.

## Structures
- distinct full InChIKeys: **1**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **1** of 1
  | InChIKey | rows here | in corpus as |
  |---|---|---|
  | `RVHYPUORVDKRTM-UHFFFAOYSA-N` | 36 | AC_2025:C12-200, JW_2024:other |
- **NEW to the corpus: 0**

## Values
- `Experiment_value`: n=36, mean **2.77778e-08**, sd(pop) **1**, min -1.857, max 1.207
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'KK_2015': 36}
- `Model_type`: {'Mouse_unspecified': 36}
- `Model_target`: {'liver': 36}
- `Cargo_type`: {'EPO': 36}
- `Experiment_method`: {'protein_adundance_normalized': 36}
- `Experiment_batching`: {'individual': 36}
- `value_normalization_unit`: {'experiment': 36}

## Quality columns (schema_quality.json)
- `confidence`: 36/36 populated — {'medium': 36}
- `positional_ambiguity`: 36/36 populated — {'0': 36}
- `structure_assumed`: 36/36 populated — {'1': 36}
- `value_is_bound`: 0/36 populated — {'': 36}
- `curation_note`: 36/36 populated — {"C12-200's structure is NOT attested by this paper: no drawing, no formula, no HRMS, no MW, no synthesis anywhere in the SI, and the main text is closed-access. It is carried over from the corpus (AC_2025 498 rows + JW_2024 20 rows, InChIKey RVHYPUORVDKRTM-UHFFFAOYSA-N). The PEG-lipid is never named in the SI either and is emitted as 'Unknown' with no structure.": 35, "C12-200's structure is NOT attested by this paper: no drawing, no formula, no HRMS, no MW, no synthesis anywhere in the SI, and the main text is closed-access. It is carried over from the corpus (AC_2025 498 rows + JW_2024 20 rows, InChIKey RVHYPUORVDKRTM-UHFFFAOYSA-N). The PEG-lipid is never named in the SI either and is emitted as 'Unknown' with no structure. Cholesterol mol% is printed as 31.5 (four components sum to 110); emitted as the over-determined 21.5 -- see 08_contradictions.csv C-01.": 1}

## Folder inventory
  - `curation/andersondg_26469188/00_value_sources.json` (11780 B)
  - `curation/andersondg_26469188/01_blocks.csv` (1306 B)
  - `curation/andersondg_26469188/02_corpus_overlap.json` (14538 B)
  - `curation/andersondg_26469188/03_enumerated.csv` (3309 B)
  - `curation/andersondg_26469188/05_values.csv` (24368 B)
  - `curation/andersondg_26469188/06_rows.csv` (35638 B)
  - `curation/andersondg_26469188/07_qc.json` (38797 B)
  - `curation/andersondg_26469188/08_contradictions.csv` (8269 B)
  - `curation/andersondg_26469188/DECISIONS.md` (9186 B)
  - `curation/andersondg_26469188/PROVENANCE.md` (18361 B)
  - `curation/andersondg_26469188/_checks.json` (10252 B)
  - `curation/andersondg_26469188/_emit.json` (434 B)
  - `curation/andersondg_26469188/_parsed.json` (21036 B)
  - `curation/andersondg_26469188/_seeds.json` (4734 B)
  - `curation/andersondg_26469188/_structures.json` (4878 B)
  - `curation/andersondg_26469188/make_reference.py` (20432 B)
  - `curation/andersondg_26469188/reference_structures.csv` (3352 B)
  - `curation/andersondg_26469188/stage01_extract.py` (16477 B)
  - `curation/andersondg_26469188/stage02_checks.py` (18115 B)
  - `curation/andersondg_26469188/stage03_corpus.py` (22175 B)
  - `curation/andersondg_26469188/stage04_structures.py` (10361 B)
  - `curation/andersondg_26469188/stage05_emit.py` (16033 B)
  - `curation/andersondg_26469188/stage06_qc.py` (21593 B)
  - `curation/andersondg_26469188/test_seeds.py` (15126 B)

## DECISIONS.md — **10 OPEN**
  - D-1 — Admit a one-ionizable-lipid formulation screen at all? (F-067)
  - D-2 — A-03's cholesterol: emit the printed 31.5, or the over-determined 21.5?
  - D-3 — The PEG-lipid is never identified
  - D-4 — Licence: this source fails §8
  - D-5 — `confidence = medium` or `low` for a structure this paper never attests?
  - D-6 — `Cargo_type`: `EPO` or the corpus's `hEPO`?
  - D-7 — `Model_type`: the mouse strain is not stated
  - D-8 — Tables S3 and S4: 8 admissible values, not emitted
  - D-9 — F-061: this paper cannot adjudicate it, and `data/` is untouched
  - D-10 — `Model_target = liver`, but the measurement is serum

## 08_contradictions.csv — 11 rows
