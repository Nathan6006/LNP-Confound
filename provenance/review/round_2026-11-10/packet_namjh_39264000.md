# Review packet — `namjh_39264000`

Generated from live files. Queue says **108** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 39264000.

## Rows read from `curation/namjh_39264000/06_rows.csv` — 108 rows

**74 columns.** Corpus schema is 79; this file has 74 of them.

## Structures
- distinct full InChIKeys: **36**
- unparseable `IL_SMILES`: **0**
- already in the corpus: **0** of 36
- **NEW to the corpus: 36**

## Values
- `Experiment_value`: n=108, mean **5.84666e-18**, sd(pop) **1**, min -1.428, max 1.792
  - *mean 0 / sd 1 is a TAUTOLOGY after z-scoring — it is not evidence of anything.* What matters is whether the z was taken over the right block (error class 22).
- `Experiment_ID`: {'SB_2025': 72, 'SB_2025_2': 36}
- `Model_type`: {'Mouse_ICR': 108}
- `Model_target`: {'ear': 72, 'serum': 36}
- `Cargo_type`: {'RLuc': 72, 'hEPO': 36}
- `Experiment_method`: {'luminescence_normalized': 72, 'protein_adundance_normalized': 36}
- `Experiment_batching`: {'individual': 108}
- `value_normalization_unit`: {'experiment': 108}

## Quality columns (schema_quality.json)
- `confidence`: 108/108 populated — {'medium': 96, 'high': 12}
- `positional_ambiguity`: 108/108 populated — {'0': 108}
- `structure_assumed`: 108/108 populated — {'0': 108}
- `value_is_bound`: 12/108 populated — {'': 96, 'upper': 12}
- `curation_note`: 108/108 populated — {"Block(s) FH-2 are drawn only in Fig 3a and are named in no synthetic procedure, no HRMS and no NMR assignment; they were measured from that drawing by a tracer calibrated on the blocks the SI does name (2-hexyl-1-decanol, 2-octyl-1-dodecanol, glutaric anhydride, 3,5-dihydroxybenzoic acid). 3 sp3 stereocentres are stored unassigned because the paper assigns none: the glycerol centre comes from racemic 3-iodopropane-1,2-diol and the tail branch carbons from racemic branched alcohols. the value is an UPPER BOUND: this cell matches the colourbar's lightest swatch, so any value in [0.0, 56.72] fold renders identically (Fig 4c draws this curation's lowest cell at 12.1 fold against a stored 53.1, inside the band). value digitised from Fig S2c's colour ramp; the raw reading is 56.72 fold change vs the Nil group, band [0.0, 56.72].": 2, "Block(s) FH-2/EB-6 are drawn only in Fig 3a and are named in no synthetic procedure, no HRMS and no NMR assignment; they were measured from that drawing by a tracer calibrated on the blocks the SI does name (2-hexyl-1-decanol, 2-octyl-1-dodecanol, glutaric anhydride, 3,5-dihydroxybenzoic acid). 3 sp3 stereocentres are stored unassigned because the paper assigns none: the glycerol centre comes from racemic 3-iodopropane-1,2-diol and the tail branch carbons from racemic branched alcohols. the value is an UPPER BOUND: this cell matches the colourbar's lightest swatch, so any value in [0.0, 56.72] fold renders identically (Fig 4c draws this curation's lowest cell at 12.1 fold against a stored 53.1, inside the band). value digitised from Fig S2c's colour ramp; the raw reading is 56.72 fold change vs the Nil group, band [0.0, 56.72].": 2, "Block(s) EB-6 are drawn only in Fig 3a and are named in no synthetic procedure, no HRMS and no NMR assignment; they were measured from that drawing by a tracer calibrated on the blocks the SI does name (2-hexyl-1-decanol, 2-octyl-1-dodecanol, glutaric anhydride, 3,5-dihydroxybenzoic acid). 3 sp3 stereocentres are stored unassigned because the paper assigns none: the glycerol centre comes from racemic 3-iodopropane-1,2-diol and the tail branch carbons from racemic branched alcohols. value digitised from Fig S2c's colour ramp; the raw reading is 1111.11 fold change vs the Nil group, band [1111.11, 1111.11].": 2, "Block(s) FH-2/EB-7 are drawn only in Fig 3a and are named in no synthetic procedure, no HRMS and no NMR assignment; they were measured from that drawing by a tracer calibrated on the blocks the SI does name (2-hexyl-1-decanol, 2-octyl-1-dodecanol, glutaric anhydride, 3,5-dihydroxybenzoic acid). 3 sp3 stereocentres are stored unassigned because the paper assigns none: the glycerol centre comes from racemic 3-iodopropane-1,2-diol and the tail branch carbons from racemic branched alcohols. the value is an UPPER BOUND: this cell matches the colourbar's lightest swatch, so any value in [0.0, 56.72] fold renders identically (Fig 4c draws this curation's lowest cell at 12.1 fold against a stored 53.1, inside the band). value digitised from Fig S2c's colour ramp; the raw reading is 56.72 fold change vs the Nil group, band [0.0, 56.72].": 2, "Block(s) FH-1/EB-6 are drawn only in Fig 3a and are named in no synthetic procedure, no HRMS and no NMR assignment; they were measured from that drawing by a tracer calibrated on the blocks the SI does name (2-hexyl-1-decanol, 2-octyl-1-dodecanol, glutaric anhydride, 3,5-dihydroxybenzoic acid). 4 sp3 stereocentres are stored unassigned because the paper assigns none: the glycerol centre comes from racemic 3-iodopropane-1,2-diol and the tail branch carbons from racemic branched alcohols. value digitised from Fig S2e's colour ramp; the raw reading is 72.69 fold change vs the Nil group, band [66.36, 79.63].": 2}

## Folder inventory
  - `curation/namjh_39264000/00_acquisition.json` (14545 B)
  - `curation/namjh_39264000/00_value_sources.json` (12379 B)
  - `curation/namjh_39264000/01_blocks.csv` (4105 B)
  - `curation/namjh_39264000/01_blocks_qc.json` (7089 B)
  - `curation/namjh_39264000/03_enumerated.csv` (7418 B)
  - `curation/namjh_39264000/04_values_qc.json` (5229 B)
  - `curation/namjh_39264000/05_values.csv` (16028 B)
  - `curation/namjh_39264000/06_emit_qc.json` (1270 B)
  - `curation/namjh_39264000/06_rows.csv` (153177 B)
  - `curation/namjh_39264000/07_qc.json` (22085 B)
  - `curation/namjh_39264000/08_contradictions.csv` (6865 B)
  - `curation/namjh_39264000/DECISIONS.md` (11378 B)
  - `curation/namjh_39264000/PROVENANCE.md` (11322 B)
  - `curation/namjh_39264000/make_reference.py` (7423 B)
  - `curation/namjh_39264000/reference_structures.csv` (2065 B)
  - `curation/namjh_39264000/stage01_blocks.py` (26003 B)
  - `curation/namjh_39264000/stage03_enumerate.py` (12406 B)
  - `curation/namjh_39264000/stage04_values.py` (23896 B)
  - `curation/namjh_39264000/stage06_verify.py` (43425 B)
  - `curation/namjh_39264000/stage07_emit.py` (16395 B)

## DECISIONS.md — **8 OPEN**
  - D-n  <one line>
  - D-1  The paper's headline is 213 LNPs; only 36 ionizable lipids are drawn anywhere.
  - D-2  The `Experiment_ID`, and its collision risk with `SB_2024`.
  - D-3  The vaccine readouts are not emitted. Should they be?
  - D-4  Table S2's "Mass ratio" cannot be mapped onto `IL_to_nucleicacid_massratio`.
  - D-5  The source licence is CC BY-NC-**ND**, which is not CC-BY-compatible.
  - D-6  `Model_target` for the I.M. hEPO experiment: `serum` or `muscle`?
  - D-7  Are Fig S2c and Fig S2e two formulations within one experiment, or two experiments?
  - D-8  36 of 36 structures are emitted; only 2 of 36 have an external anchor.

## 08_contradictions.csv — 11 rows
