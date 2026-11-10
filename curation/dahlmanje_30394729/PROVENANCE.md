# PROVENANCE — dahlmanje_30394729

**Paper.** Sago CD, Lokugamage MP, Islam FZ, Krupczak BR, Sato M, Dahlman JE.
*Nanoparticles That Deliver RNA to Bone Marrow Identified by in Vivo Directed
Evolution.* J Am Chem Soc 140(49), December 2018. PMID **30394729**, PMCID
PMC6556374, DOI 10.1021/jacs.8b08976.

**Lab and first author, verified from the Europe PMC `core` record, not the
folder name.** `authorString` = "Sago CD, Lokugamage MP, Islam FZ, Krupczak BR,
Sato M, Dahlman JE." The folder is named after the senior author; the proposed
`Experiment_ID` derives from the **first** author and is **CS_2018**. The SI
title page lists the same six authors in the same order, all at the Wallace H.
Coulter Dept. of BME, Georgia Tech. `CS` does not collide with any of the 43
vitro / 20 vivo Experiment_IDs in the live corpus. `SS_2018` in `data_vivo.csv`
is PMID 29653760 (Sabnis et al., Moderna) — checked by PMID, not by name.

**Endpoint.** Two of them, and they are not the same kind of thing.
*Screen*: `Normalized DNA delivery (%)` — DNA barcode counts recovered from
FACS-sorted **ICAM-2-low bone marrow endothelial cells**, normalised within the
sample and then to the input pool, i.e. a **compositional, pooled** readout over
115 and 31 LNPs. *Individual*: `ICAM-2 Protein Silencing (MFI, Normalized)` and
`ICAM-2 Indel (%)` by TIDE — both genuine delivery endpoints in `data/`'s
vocabulary.

**Rows emitted: ZERO.** `06_rows.csv` is header-only so the zero is visible.

---

## What the paper actually contributes

**One ionizable lipid.** 7C1, in every formulation of both libraries and in every
individually tested LNP. The libraries vary cholesterol, DSPC and the PEG-lipid;
the ionizable lipid never changes. A corpus keyed on ionizable-lipid structure
gains **one** lipid from 160 formulations — and gains **zero**, because 7C1 is a
polydisperse branched-PEI(Mn ≈ 600)/C15-epoxide polymer for which no single
SMILES is correct, while `data/` has **0 rows with a blank `IL_SMILES`** in
either pool. See `DECISIONS.md` D-1.

| | count |
|---|---|
| distinct ionizable lipids | **1** |
| formulations formulated | 160 (120 + 40) |
| formulations pooled and measured | 146 (115 + 31) |
| formulations individually measured | 5 |
| per-LNP values the design supports | 146 |
| per-LNP values actually attributable | **5** |

The five individually measured formulations are the 'original' 80:20 7C1
(80/0/20 C14PEG2000), the screen-1 winner (60/37.5/2.5 C16PEG2000 = library 1
LNP 25), LNP1 (70/10/20 C16PEG2000 = library 2 LNP 6), **BM1** (80/2.5/17.5
C18PEG2000 = library 2 LNP 18) and LNP3 (80/0/20 C18PEG2000 = library 2 LNP 19).

---

## Where the values are, and the tier

Detail in `00_value_sources.json`. Summary: **no pixel was read in this
curation.** The composition tables, the barcode list and the normalisation
worked example all come back as native PDF text tables (tier 4). Every
value-bearing panel is a raster, and digitisation was ruled out before any
geometry code was written — for the pooled panels because nothing is
attributable, and for the individual panels because the structure blocker is
upstream of the value and six of the seven headline numbers are printed as
numerals in the Results text anyway.

`triage.py`'s text heuristic reported `table_pages = [1]`, which is wrong in both
directions: page 1 is the title page and it missed all four real tables.
`native_tables` is the field that was read, by hand, as PIPELINE_V2 warns.

Tier 3 is an explicit null: `find_repo.py --pmid 30394729 --author Sago` found
nothing, and `GSE[0-9]+` / `SRP` / `PRJNA` / `SRA` / `E-MTAB` / accession /
deposit / github / zenodo / dryad / figshare return **zero hits** in the full
JATS body (34,392 chars) and in all 23,227 characters of SI text. There is no
Data Availability statement.

---

## Checks, with numbers

Full output in `07_qc.json`.

- **The normalisation, reconstructed from the paper's own worked example**
  (SI Fig S1E): `nd_i = [(c_i/Σc)/(in_i/Σin)] / Σ_k[same] × 100%`. Raw column
  totals 4/4 exact; **stage 3: 9/9 cells exact, max residual 0.0**. Two printed
  cells of the *intermediate* table are wrong (34 where 34.78→35, and **3** where
  19.57→20), which the exactness of stage 3 identifies as typos rather than a
  different formula. Confirmed by eye on a 400 dpi render, not only from the text
  layer.
- **A round number drawn and never computed**: `Total (%) = 100` sits under all
  three tables. The stage-2 columns actually sum to **100 / 83 / 101 / 99** and
  every stage-3 column sums to **101**.
- **Factorial proved, not assumed**: library 1 = 5 PEG types × 24 ratios = 120;
  library 2 = 2 × 20 = 40. Every PEG variant carries the **identical** ratio set
  exactly once, zero duplicates. Fig 1D draws "2×; 160 LNPs tested in vivo" and
  never computes it; 120 + 40 = 160.
- **mol% sums**: **158 of 160** formulations sum to exactly 100. The two that do
  not are library 2 LNP 36 and LNP 40, both 60/45/5/0 = **110**, the same ratio in
  both PEG blocks at stride 4. Not a tautology of our own derivation — the four
  columns are transcribed independently and nothing computes one from the others.
- **The 'Formulated' pies against the composition tables (check E1)**: 5 of 6
  reconcile **exactly**. Fig 4A prints 12/25/25/21/17 and the table gives
  12.5/25/25/20.83/16.67 → 12/25/25/21/17, 5/5, which also pins the band
  convention as `[lo,hi)` with ">20" meaning ≥20; Fig 4B 20×5, Fig 5I 50/50,
  Fig 5J 20/40/20/20, Fig 5K 60/40 all exact. **Fig 5H does not reconcile** under
  either convention (printed 40/40/20 vs 45/20/35 or 45/40/15).
- **156 barcodes**: 156 listed, all 8 nt, all distinct, numbered 1..156 — matches
  the claim exactly.
- **The five individually tested formulations** are located in the library
  tables: BM1 = library 2 LNP 18, LNP1 = LNP 6, LNP3 = LNP 19, screen-1 winner =
  library 1 LNP 25. The screen-1 winner's ratio also appears as library 2 LNP 34,
  which is by design (library 2 was evolved from library 1).
- **PEG repeat count**: the drawn `(OCH2CH2)45OCH3` gives 45 × 44.05 = 1982.4 Da,
  0.9 % from the printed "2000"; never multiplied out in the paper.

### Seeded-mutation gate — `09_seeds.py`, 6/6 credited

One child process per seed; the verdict is the child's return code, and a seed
whose `METRIC` line is absent is **not** credited however it exited.

| seed | metric | control → seeded |
|---|---|---|
| N1 raw count +10 | stage-3 exact cells | 9 → 6 |
| M1 one cholesterol value | library-1 mol% violations | 0 → 1 |
| F1 one PEG type | factorial holds | 1 → 0 |
| B1 drop a barcode | n barcodes | 156 → 155 |
| W1 shift every IL mol% by 0.5 | winners uniquely located | 3 → 0 |
| E1 move one ratio in all 5 PEG blocks | pies reconciling | 5 → 4 |

**A seed that was refused.** The first E1 seed moved one library-1 PEG value from
20 to 22.5 and was a **no-op**, because the half-open convention already puts 20
in the ">20" band. It exited non-zero — on the factorial, not on E1 — while E1's
metric never moved, and the harness refused it. That is why `credited` requires
the metric to move and not merely a non-zero exit. It was replaced by a seed that
moves the same ratio in all five PEG blocks, which leaves the factorial, the row
count and every mol% sum intact and can therefore only be caught by E1.

### Emit guards — exercised, not assumed

`07_emit.py` runs the five formulations that *would* have been emitted through
`emit.emit_rows`, and both guards fire on purpose: a literal `confidence` raises
(`emit: confidence='low' is a LITERAL`, error class 29 / F-011 / F-025 / F-046)
and a caveat with an empty note raises (F-024). The control passes with 5 rows,
`rows_with_a_caveat` = 5, each note naming the polymer problem on the row.

---

## The C10 anchor — posture, recorded by hand

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **NOT run** from this folder: they write tracked
files other agents read. The comparison was done in-process by
`make_reference.py` and is recorded here and in `10_anchor_report.json`.

**The ionizable lipid is never drawn**, so `reference_structures.csv` carries a
row for `7C1` with an **empty** `reference_SMILES` and a note saying why — the
gap is visible rather than absent. The anchors are the things the paper does
draw.

| anchor | route | result |
|---|---|---|
| `DSPC` | fragment assembly: sn-glycero-3-phosphocholine + 2 stearoyl, through `assemble.attach` with atom-mapped attachment atoms; no atom deleted | **full InChIKey matches** the registry and the corpus (`NRJAVPSFFCBXDT-HUESYALOSA-N`; 6,401 vitro + 637 vivo corpus rows, 1 distinct structure each) |
| `7C1_representative_arm` | fragment assembly: diethylamine + the epoxide-opened C15 chain, N→CH2 pinned by atom map | `SSPNMOBMZCGVOW-UHFFFAOYSA-N`, flagged `expected_mismatch` — it is a cartoon of the polymer, not the lipid |
| `Cholesterol` | hand re-traversal from the drawn structure, C3-OH first | **connectivity only**, `expected_mismatch = stereo_not_independently_derived` |
| `7C1` | — | empty, on purpose |

**Independence is executed, not asserted.** `_ast_guard()` walks
`make_reference.py`'s own AST at run time and refuses if the module imports a
stage module of this folder, names a file this folder emits, or contains a string
literal byte-identical to a `reagent_registry` SMILES. A substring scan would
flag the module's own docstring; an AST walk does not. The guard's first run
**failed on the module itself**, because its ban lists are string literals naming
the files they ban; the fix exempts exactly those two assignments by node
identity. `--selftest` injects each violation into a temp copy in a child process
and requires the refusal: **3/3 fire with the right message, the unmutated
control passes.** (The first version of the import seed was `import 06_checks`,
which is a `SyntaxError` — the child died before reaching the guard and was
refused for having no metric line, exactly as designed; it was replaced with a
parseable module name.)

**Byte-identity is not the signal (F-037).** DSPC's assembled anchor
canonicalises byte-identical to the registry string. That is neither evidence of
copying nor of independence: a correct re-derivation is byte-identical to a
correct original once canonicalised, which is precisely why C10's
`anchor_may_be_a_copy` missed F-037.

**The gate is shown to test connectivity.** The seeded defect is the Markovnikov
regiochemistry of the epoxide opening — amine on the substituted carbon instead
of the terminal one. Identical formula `C19H41NO`, identical exact mass
`299.318815` (difference **0.0 Da**), different InChIKey and different skeleton
block. Mass and formula cannot see it.

**And the anchor caught a real error — ours.** The first cholesterol anchor wrote
all eight centres by hand from the drawn wedges and produced
`HVYWMOMLDIMFJA-OXMOYTRBSA-N`: a **diastereomer** with four of eight centres
inverted, identical in formula, in exact mass to six decimals, and in the
InChIKey skeleton block — the F-013 class reproduced from scratch. It is not
ent-cholesterol either (`-HFRBDLMJSA-N`). The honest repair is **not** to flip
tokens until the key matches the registry, which is fitting the anchor to the
answer it exists to test, so the stereo was **discarded** and the gap declared on
the row. The discarded string is kept verbatim in `make_reference.py` and in
`10_anchor_report.json`. The same build also refused itself once before that: a
bracket oxygen written `[O:3]` is a radical, and `assemble.check_no_radicals`
raised.

No excipient SMILES in this folder is typed from memory: `01_blocks.csv` takes
cholesterol and DSPC from `provenance/reagent_registry.py`.

---

## Corpus overlap

- **Token-pattern test.** A substring join on the paper's own label `7C1`
  captures **36 corpus rows spanning 35 distinct structures across 4
  experiments** — `KL_2010` `117C14`, `LL_2012` `A7C10`, `LM_2019` `A1Iso7C11`,
  `ZH_2023` `A4I18R7C18-2` — and **not one** is Dahlman's 7C1 lipomer. The
  paper's `7C1` is a compound index from Dahlman 2014: the `7` is the amine and
  the `C1` the epoxide index, while the drawn tail is `C13` and the epoxide is
  C15. Three numbers in one name, none of them a chain length.
- **Substructure sweep, labelled as a SUPERSET.** The loose amino-alcohol motif
  (`[NX3]CC(O)C`) matches **1,247 of 18,720** vitro and **276 of 1,506** vivo
  structures — a chemotype, not an identification. The strict pattern (the full
  drawn C13 arm) matches **0 and 0**.
- **Excipients, compared by InChIKey rather than name** (324 of 19,919 corpus
  structures carry more than one `IL_name`): cholesterol
  `HVYWMOMLDIMFJA-DPAQBDIFSA-N` and DSPC `NRJAVPSFFCBXDT-HUESYALOSA-N`, one
  distinct structure each across 28,684 + 1,713 and 6,401 + 637 corpus rows.
- **In-vivo shared chemistry.** The pool has exactly **two** structures shared
  across experiments — DLin-MC3-DMA (6 sources) and SM-102 (4). **This paper adds
  none**, and uses no benchmark ionizable lipid at all: no MC3, no SM-102, no
  ALC-0315 anywhere. Its only comparator is its own 80:20 7C1. The "do benchmarks
  top the screen" question therefore does not arise; and BM1 topping screen 2 is
  not an independent finding, since BM1 *is* the top performer by definition.
- **PEG names.** `C14/C16/C18-PEG2000/3000` are not in the corpus `PEG_name`
  vocabulary. The drawn skeleton implies DMPE-/DPPE-/DSPE-PEG2000, two of which
  are; that is an inference and is recorded, not applied (`DECISIONS.md` D-4).

---

## Licence

**Article: NOT CC-BY.** Europe PMC `core`: `isOpenAccess = N`, `license = null`;
the PMC copy is an NIH author manuscript whose `<permissions>` block reads
"available for text mining ... fair use", and the `supplementaryFiles` endpoint
refuses with "not open access". Text-mining permission is not a redistribution
licence, so under `CLAUDE.md` §8 this source fails the CC0/CC-BY requirement for
corpus release. **Repository: none exists.**

---

## Files

| file | what |
|---|---|
| `00_acquisition.json` | sources, hashes, routes, licence |
| `00_value_sources.json` | the intake gate, written before any extractor |
| `01_si_tables.py` → `04_si_tables.json`, `01_blocks.csv`, `03_enumerated.csv` | native SI tables and the 160-formulation enumeration |
| `02_corpus_overlap.py` → `02_corpus_overlap.json` | token test, substructure superset, excipients, shared chemistry |
| `05_values.py` → `05_values.csv`, `05_values_summary.json` | every value-bearing panel, `corpus_admissible` with a reason, counts summing (662 catalogued, 0 admissible) |
| `06_checks.py` → `07_qc_checks.json` | N1/F1/M1/B1/W1/P1/E1 |
| `07_emit.py` → `06_rows.csv` (header-only), `06_rows_report.json` | the zero, and the emit guards exercised |
| `09_seeds.py` → `09_seedgate.json` | the seeded-mutation gate, 6/6 |
| `make_reference.py` → `reference_structures.csv`, `10_anchor_report.json` | the anchors and the AST independence guard |
| `10_qc.py` → `07_qc.json` | the assembled QC |
| `08_contradictions.csv` | 11 entries, 3 open |
| `DECISIONS.md` | D-1..D-5, all OPEN |

Run order: `01_si_tables.py`, `02_corpus_overlap.py`, `05_values.py`,
`06_checks.py`, `07_emit.py`, `09_seeds.py`,
`make_reference.py --selftest`, `10_qc.py`.
