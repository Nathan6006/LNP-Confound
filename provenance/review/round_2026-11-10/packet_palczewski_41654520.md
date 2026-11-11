# Review packet — `palczewski_41654520`

Generated from live files. Queue says **64** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 41654520.

## Rows read from `curation/palczewski_41654520/05_values.csv` — 64 rows

**33 columns.** Corpus schema is 79; this file has 7 of them.

## Structures
- distinct full InChIKeys: **30**
- structure column: **`canonical_SMILES`**
- unparseable structures: **0**
- already in the corpus: **0** of 30
- **NEW to the corpus: 30**

## Values
- `value_normalization_unit`: {'': 64}

## Quality columns (schema_quality.json)
- `confidence`: 64/64 populated — {'high': 59, 'low': 5}
- `positional_ambiguity`: 1/64 populated — {'': 63, '1': 1}
- `structure_assumed`: 5/64 populated — {'': 59, '1': 5}
- `value_is_bound`: 0/64 populated — {'': 64}
- `curation_note`: 64/64 populated — {'CBB11: CBB21 bis-acylated on both terminal primary amines with linoleoyl (linoleic acid, (9Z,12Z)-octadeca-9,12-dienoic); connectivity from the SI IUPAC name, arm count from its HRMS.': 3, "CBB14: CBB13 exhaustively methylated with methyl iodide; the SI name says 'dimethylammonio' and 'N,N-dimethylpropan-1-aminium', so both linker nitrogens are quaternary. Emitted as the dication.": 3, "CBB1: both terminal amines carbamoylated with cholesteryl chloroformate. Cholesterol is taken verbatim from provenance/reagent_registry.py and esterified by assemble.esterify, so the acid's OH leaves and C3-beta is untouched; the SI name's (3S,...)-cyclopenta[a]phenanthren-3-yl is cholesteryl, 3-beta.": 2, "CBB2: both terminal amines carbamoylated with cholesteryl chloroformate. Cholesterol is taken verbatim from provenance/reagent_registry.py and esterified by assemble.esterify, so the acid's OH leaves and C3-beta is untouched; the SI name's (3S,...)-cyclopenta[a]phenanthren-3-yl is cholesteryl, 3-beta.": 2, 'CBB3: CBB21 bis-acylated on both terminal primary amines with dha (docosahexaenoic acid, (4Z,7Z,10Z,13Z,16Z,19Z)-docosa-4,7,10,13,16,19-hexaenoic); connectivity from the SI IUPAC name, arm count from its HRMS.': 2}

## Folder inventory
  - `curation/palczewski_41654520/00_acquisition.json` (4736 B)
  - `curation/palczewski_41654520/00_value_sources.json` (8095 B)
  - `curation/palczewski_41654520/01_blocks.csv` (5212 B)
  - `curation/palczewski_41654520/03_enumerated.csv` (23364 B)
  - `curation/palczewski_41654520/05_controls.json` (1284 B)
  - `curation/palczewski_41654520/05_values.csv` (73587 B)
  - `curation/palczewski_41654520/07_qc.json` (35658 B)
  - `curation/palczewski_41654520/08_contradictions.csv` (8086 B)
  - `curation/palczewski_41654520/DECISIONS.md` (8907 B)
  - `curation/palczewski_41654520/make_reference.py` (7589 B)
  - `curation/palczewski_41654520/paper_library.py` (16499 B)
  - `curation/palczewski_41654520/reference_structures.csv` (1403 B)
  - `curation/palczewski_41654520/stage01_blocks.py` (4445 B)
  - `curation/palczewski_41654520/stage03_enumerate.py` (17188 B)
  - `curation/palczewski_41654520/stage05_values.py` (10917 B)
  - `curation/palczewski_41654520/stage07_qc.py` (40580 B)

## DECISIONS.md — **6 OPEN**
  - D-1
  - D-2
  - D-3
  - D-4
  - D-5
  - D-6

## 08_contradictions.csv — 12 rows
