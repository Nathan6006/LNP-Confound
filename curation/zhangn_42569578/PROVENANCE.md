# PROVENANCE — `zhangn_42569578` → proposed `MY_2026`

PMID **42569578** · PMC13450591 · doi 10.1016/j.mtbio.2026.103483 ·
*Materials Today Bio* 2026 · **CC BY 4.0**

> "High-throughput microfluidics and machine learning-assisted screening of
> lipid nanoparticle formulations for siRNA delivery"
> Yu M, Wang L, Liu D, Guo J, Liu J, Zhou J, Mathew A, He Z, Wu J, Zhang N.

**First author: Mingzhi Yu.** The folder is named for **Nan Zhang**, the last
author and PI. `Experiment_ID` derives from the first author, so the key is
**`MY_2026`**, not `ZN_*`. Checked against all 54 corpus `Experiment_ID`s **by
PMID**: 42569578 is absent, and no `MY_*` prefix exists (`MO_2026` is PMID
42050330, a different paper).

---

## 1. The anchor posture — recorded by hand

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run** (out of scope for this task), so
the C10/C16/C18 posture is stated here instead of being measured by them.

**There is no C10 anchor for this paper, and there cannot be one.** The paper
draws **no chemical structure anywhere**. Fig 1 and the graphical abstract
render the four lipid classes as coloured cartoon blobs; Figs 2–7 are DLS
scatter, Spearman/SHAP panels, bar charts and IVIS organ photographs. All four
ionizable lipids are catalogue reagents — SM-102, ALC-0315 and DLin-MC3-DMA
from Avanti, CKK-E12 from MedChemExpress — so there is no "structure the paper
draws" to compare against.

`reference_structures.csv` therefore carries **four rows with an empty
`reference_SMILES`**, which is the posture `CLAUDE.md` asks for: store the empty
row so the *absence* of an anchor is visible rather than invisible.

### What this means for the guards, stated plainly

`make_reference.py` implements both required guards and `test_seeds.py`
demonstrates that each one fires when seeded:

| guard | mechanism | seeded result |
|---|---|---|
| A — independence | walks the module's **AST at run time**, inspecting only `Import` / `ImportFrom` / `Call` / `Attribute` nodes, so the ban list (a plain list of strings in the same file) cannot match itself | `import corpus_components` → 0 → **1**, exit 1 |
| B — template provenance | no string literal ≥ 25 chars may be a substring of any enumerated SMILES | a 40-char slice of an emitted SMILES → 0 → **1**, exit 1 |

**Against an empty anchor both guards are VACUOUS, not passing.** There is no
stored structure for a copy to corrupt. Writing "C10 clean" here would be the
F-013 / F-037 mistake in its purest form — comparing a value to itself and
reporting a match.

Byte-identity of an emitted anchor was deliberately **not** used as the signal
(F-037): a correct re-derivation is byte-identical to a correct enumeration
once canonicalised, which is exactly why C10's own `anchor_may_be_a_copy`
missed the defect it existed to catch.

### What *is* non-vacuous

A connectivity probe, run on every invocation: a constitutional isomer of
SM-102 built by moving one methyl along a chain has **identical molecular
formula (C44H87NO5)** and an **exact-mass difference of 0.0000000000 Da**, and
is separated by full InChIKey and by skeleton block. So the gate tests
connectivity, not arithmetic.

### Structure provenance, and one honest non-independence

| lipid | source | independent of the corpus? |
|---|---|---|
| ALC-0315 | `provenance/reagent_registry.py`, written from the IUPAC name | **yes** |
| SM-102 | `provenance/reagent_registry.py` | **yes** |
| MC-3 (DLin-MC3-DMA) | `provenance/reagent_registry.py` | **yes** |
| CKK-E12 | the live corpus, via `curation/corpus_components.smi` | **NO** |

cKK-E12 has no registry entry. Its structure was taken *from* the corpus, so
its agreement *with* the corpus is circular and **is not counted as a check**.
It is declared non-independent rather than tallied, and the honest remedy is a
`reagent_registry.py` entry written from the chemical name — which is outside
this task's write scope.

No structure was hand-written. The registry exists precisely because a sibling
hand-wrote SM-102, MC3 and ALC-0315 and all three were wrong — one CH₂ short,
one CH₂ short on a tail, one ester reversed — all parsing and canonicalising
cleanly, with only the InChIKey catching it.

---

## 2. Value provenance — tier 1, with no pixel read

Every candidate value comes from the **publisher Source Data workbook**
(`mmc2.xlsx`) and the **SI tables** (`mmc1.docx`). **No digitiser was written
or run for this paper**: no axis calibration, no bar-top detection, no heatmap
LUT. The digitisation defect surface is **zero**, which `PIPELINE_V2` §3 records
as the single highest-value property a curation can have (0 defects for the two
tier-1 curations against 8 for the four that digitised).

The figures were used **only to verify**, which is strictly stronger than using
them to extract:

- Fig 5c's y-axis — "S100P/GAPDH mRNA expression (fold change relative to
  control)", with an Untreated bar at 1.0 — independently establishes what
  `mmc2.xlsx` column V actually holds, against a header that calls it
  "Fluorescence intensity".
- Figs 4c, 5a and 5b confirm that all three measurement conditions share one
  25-formulation set.

The intake gate (`00_value_sources.json`) is satisfied in substance: no
digitiser was ever written, so there was nothing for the gate to gate.

---

## 3. Acquisition — verified, and completed

The handover's route was **re-verified, not trusted**. Both supplements carry
PK-ZIP magic bytes and unzip to valid OOXML; md5 and sha256 for all 13 source
files are in `00_value_sources.json`; the workbook was re-dumped from its raw
sheet XML rather than read off the handover's summary.

Two things the handover got wrong or left incomplete:

1. **The table is 25 formulations, not 51.** `mmc2.xlsx` has 51 `<row>`
   elements of which **26 carry cells** (1 header + 25 data); the other 25 are
   empty. 25 matches the paper's "25 candidates for siRNA loading" exactly.
   The 22 columns are real, but the header occupies **18** of them — A, H, M
   and R are blank spacers.
2. **Figures 6 and 7 were missing.** The handover supplied `ga1` and `gr1–gr5`;
   the JATS references `gr1–gr7`. Both missing figures are inadmissible
   endpoints, so no emitted value depended on them, but the set was incomplete.
   Fetched `gr4/5/6/7_lrg.jpg` from the open Elsevier CDN route — Fig 5 at
   2371 × 2857 px against 790 × 952 for the PMC copy, 3× linear and 9× pixels.

No paywall and no bot check was bypassed.

### A contamination incident worth recording (F-065)

The session scratchpad **already contained** a directory named `repo/` holding a
**different project** — `ml_lnp_design_rules`, an ACS Nano paper, complete with
its own `.git`, timestamped 07:41 the same day. My first repository fetch landed
*inside* it. It was detected by comparing against this paper's own git-tree
listing, which contains only four `.py` files plus `.idea`; the foreign
directory had `Raw_Data/`, `utilities.py`, a `LICENSE` and a PDF.

Had I grepped `repo/Raw_Data/` believing it was this paper's repository, I would
have manufactured a finding out of another paper's data. Everything was
re-fetched into a uniquely named empty directory and nothing from the
contaminated one was read. **A shared scratchpad directory name is not a safe
container.**

`00_acquisition_raw.json` was left in the folder by the handover at **0 bytes**;
it is superseded by `00_value_sources.json` and carries no information.

---

## 4. Accession sweep — explicit NULL result

Searched **1,895,389 characters** across all three machine-readable sources
(`PMC13450591.jats.xml` 89,772; `mmc1.docx/word/document.xml` 1,783,844;
`mmc2.xlsx` sheet + sharedStrings 21,773) for
`GSE\d+`, `SRP\d+`, `PRJNA\d+`, `PRJEB\d+`, `SR[RX]\d+`, `E-MTAB-\d+`,
`accession`, `deposit`, `github`, `zenodo`, `dryad`, `figshare`, `osf.io`,
`gitlab`.

**NULL for every data-repository pattern.** The single non-null hit is one
`github` occurrence in the main text:
`github.com/luhan611/siRNA-LNP-ML-Optimization` (HTTP 200, public, owned by the
second author Luhan Wang, who is credited with "Software"). It holds **four
`.py` files and no data** — the input workbooks are referenced by absolute local
paths under `/home/luhan/ucd/` and are not committed. The repository carries
**no licence**.

A repository search was run for both the first author and the PI; this is the
only one associated with the paper.

The Data availability statement reads *"Data will be made available on
request"*, which **understates** what is published: `mmc2.xlsx` is a complete
Source Data table and is the reason this curation is tier 1.

---

## 5. Files

| file | what it is |
|---|---|
| `00_value_sources.json` | intake gate, acquisition, md5/sha256, accession sweep |
| `01_blocks.csv` | the design factors and levels (Table 1) — **not** a chemical enumeration |
| `02_corpus_overlap.json` | InChIKey joins, three positive controls, token-pattern test |
| `03_enumerated.csv` | the four ionizable lipids, with structure source and independence |
| `05_values.csv` | **300 values**, each with source cell, readout as printed, and `corpus_admissible` + reason |
| `06_rows.csv` | **header-only** — 79 columns read live from `data/data_vitro.csv` |
| `07_qc.json` | every check with numbers; regenerated by `stage06_qc.py` |
| `08_contradictions.csv` | 15 entries, including one retracted suspicion with the measurement that retired it |
| `DECISIONS.md` | 5 entries, all `OPEN` |
| `reference_structures.csv` | 4 rows, empty `reference_SMILES` |
| `stage01_stage.py` · `stage03_corpus_overlap.py` · `stage06_qc.py` · `make_reference.py` · `test_seeds.py` | the pipeline |

Nothing outside this folder was written. `data/` was not touched.
