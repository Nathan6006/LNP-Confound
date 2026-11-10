# PROVENANCE — `dahlmanje_30216729`

**Verdict: a well-argued ZERO. 0 rows emittable, 1 distinct ionizable lipid, and
that lipid has no writable structure.**

---

## 1. Identity, resolved from the Europe PMC `core` record

Single call: `search?query=EXT_ID:30216729 AND SRC:MED&resultType=core`.

| field | value |
|---|---|
| full title | **Modifying a Commonly Expressed Endocytic Receptor Retargets Nanoparticles in Vivo.** |
| first author | **Sago CD** — Cory D Sago |
| last author (PI) | Dahlman JE — James E Dahlman |
| journal / year | Nano Letters, 2018 |
| PMID / PMCID / DOI | 30216729 / PMC6426696 / 10.1021/acs.nanolett.8b03149 |
| `abstractText` length | **1497** chars (so not a no-original-data record) |
| Methods / Results present | yes, both |
| `isOpenAccess` | `N` |
| `license` | `null` |

### `Experiment_ID` would be `CS_2018`, not anything with "Dahlman" in it

The folder is named for the **PI**. `CLAUDE.md` §3 keys on the **first author**,
so the ID derives from Sago CD → `CS_2018`.

Checked against **all 54 live corpus IDs by PMID, not by the year token**:

- `CS_2018` is **absent** from the corpus and from every sibling curation folder
  (grepped repo-wide).
- The corpus's existing `SS_2018` is **PMID 29653760**, a different paper — the
  trap a sibling hit, confirmed not to apply here.
- No ID collision, so no `_2` suffix would be needed.

Moot in practice: zero rows are emitted, so no ID is claimed.

---

## 2. The endpoint, and why it is not a delivery endpoint

**What is measured:** the abundance of the administered DNA barcode itself —
ddPCR absolute counts, and Illumina deep-sequencing "normalized delivery" (%).

The cargo is a **91-nt chemically modified single-stranded DNA barcode**
(Methods, "DNA Barcoding"; sequence in SI Fig S1a). It encodes no gene, is never
transcribed and is never translated. The readout is PCR amplification of that
same DNA molecule.

**`readout_requires_expression`: NO.** Nothing downstream of translation is
measured anywhere in the paper.

The paper's own vocabulary agrees: the title says *retargets nanoparticles*, the
abstract says *"we designed a platform named QUANT to study nanoparticle
biodistribution"*, the Discussion says *"macrophage uptake changes with tissue
type"*, and the axis labels are `Normalized Nanoparticle Biodistribution (%)`,
`Nanoparticle Accumulation (Au)` and `% Normalized Nanoparticle **Uptake**`.
Biodistribution and uptake are both on the brief's explicit not-delivery list.

Two panels are labelled `Absolute Delivery (AU)` (Fig 5f, Fig 6b). That is a
label on the product of a compositional share and a ddPCR count; the molecule
measured is still the barcode. The word on the axis was not treated as changing
what was measured.

**The `LR_2023` precedent does not reach this paper.** Pooling is not the
obstacle — `LR_2023` contributes 407 of 1,803 in-vivo rows as `barcoded` — but
that precedent requires a **translated** readout, and a screen sequencing cargo
abundance before translation matches none of `data/`'s in-vivo
`Experiment_method` values.

### The one admissible readout class in the paper has no data

Methods, "Endothelial RNA interference", describes an **ICAM2 protein MFI**
readout, which *is* admissible. But only the non-targeting control siLuciferase
is ever injected and its ICAM2 MFI is *defined* as 100%, so no knockdown exists.
Measured: `ICAM` occurs **2** times in the 165,357-char article XML, both inside
that one Methods paragraph, and **0** times in the 25,876-char SI text. Methods
also mentions "in vivo gene editing" timing for an experiment the paper does not
contain. Both recorded in `08_contradictions.csv`.

---

## 3. How many ionizable lipids: ONE, and it cannot be written

**`7C1`, on 248 of 248 formulations.** Measured from the two SI formulation
tables as native text (`stage01_si_tables.py`):

| | library 1 (SI p14) | library 2 (SI p15) |
|---|---|---|
| rows | **128** (paper states 128) | **120** (paper states 120) |
| distinct `Lipomer Type` | `['7C1']` | `['7C1']` |
| helper lipid | DSPC on all 128 | none on all 120 |
| PEG-lipid combinations | 4, **32 each** | 5, **24 each** |

The paper says so twice in its own words: Fig 4c, *"LNP library 1 was
synthesized with the amine 7C1, Cholesterol, DSPC, and PEG compounds at variable
molar ratios"*, and Results, *"all 226 LNPs included the 7C1 compound"*.

**The variable here is the biology.** The screens are run in wild-type and
**Cav1⁻/⁻** mice (JAX 007083); the experimental factors are the genotype and the
excipient molar ratios. The ionizable chemistry never changes. This is the
seventh one-lipid screen of the run and the cleanest case of it — the title's
"modifying a receptor" really does mean the variable is not chemistry.

**Positive control on that finding** (F-070): a constant column and an unread
column look identical in a report, so the sweep also checks a column known to
**vary** — `Diameter (nm)` returns **66 distinct values** across library 1 — and
a column known constant by a different mechanism (`Chol Type`). Had the
extractor been returning blanks, every column would have read as one value and
the `7C1` result would have been vacuous.

**7C1 has no atom-resolved structure.** SI p14 draws it as a generic scheme: an
epoxide labelled `R=C13` (a C15 epoxide) + `PEI600`, 90 °C / 72 h / EtOH.
PEI600 is branched polydisperse polyethylenimine, nominal MW 600, and the degree
of substitution is never stated; no formula, no HRMS, no NMR. It is the
ionizable-component analogue of a PEG-lipid, which `CLAUDE.md` §3 deliberately
keeps out of the reagent registry because the block length is a distribution.
**No representative monomer was written** — see `DECISIONS.md` D-3.

---

## 4. Where the values are, and the tier

`00_value_sources.json` is the intake gate and was written **first**. The gate
closes the paper at **tier 2**:

- **tier 1** (publisher Source Data) — does not exist; ACS published no Source
  Data file for this 2018 paper.
- **tier 2 / 4** (SI tables, native text) — **exist and were fully extracted**:
  both formulation tables (248 rows) and the normalisation worked example.
  They contain formulation metadata and diameters, and **not one delivery
  endpoint**.
- **tier 3** (author repo) — **explicit null result**, see §7.
- **tier 5** (pixels) — **not used, and nothing was rendered**. The per-LNP
  values exist only as unlabelled scatter points, so digitising them fails on
  *attributability*, not resolution, and no dpi fixes that. No
  `axis.native_tol`, no `calibrate`, no `bars.bar_top` call was made anywhere in
  this folder; `07_qc.json:C_vector_geometry_posture` records that no native dpi
  was measured **because nothing was rasterised**, rather than guessing one.

`05_values.csv` inventories **44 panels**, every one carrying
`corpus_admissible` + a reason. Counts sum: **0 yes / 44 no / 44 total**.

---

## 5. Compositional, with the residual measured

**Is the readout compositional? YES, twice over.** `stage01_si_tables.py`
reconstructs the paper's own worked example (SI Fig S4e, SI p13 — three stacked
tables for 3 LNPs × 3 mice × 2 organs):

| stage | formula recovered | cells | exact | max abs residual |
|---|---|---|---|---|
| A | `pct_i = raw_i / Σ_j raw_j × 100` | **21** | **21/21** | **0.0** |
| B | `(pct_i / pct_i,input)` re-closed to sum 100 | **18** | **18/18** | **0.0** |

Stage B is computed from the **unrounded** raw counts, so what is reconstructed
is the pipeline, not the printed integers.

**State the precision honestly.** The residual is 0.0 *against integers the
authors rounded*. The pre-rounding mean |discrepancy| is **0.299** and **0.291**
percentage points, which is the signature of rounding to whole percent (expected
mean 0.25 for uniform rounding), not of a wrong formula. **Measured sensitivity
limit: perturbations smaller than ±0.5 pp are not detectable by this check.**
A seeded +5 pp change to one printed cell moves stage A from 21/21 → 20/21; a
seeded removal of the input-normalisation step moves stage B from 18/18 → 1/18.

The SI prints `Total (%) = 100` on every column of both stages. A value is a
**share of a pool**: it moves when any other pool member moves, and it is not
comparable to a value from a pool of different membership.

**The Methods formula alone would not have been enough.** The entire
specification is one sentence — *"Counts for each particle, per cell type, were
normalized to the barcoded LNP mixture applied to cells or injected into the
mouse."* — which names only the input normalisation and states neither the
per-sample compositional step nor the re-closure. Both were recovered from
Fig S4e. Recorded in `08_contradictions.csv`.

---

## 6. Attributability — the independent blocker

| quantity | count |
|---|---|
| distinct ionizable lipids | **1** |
| formulations | **248** (128 + 120) |
| formulations individually measured (diameter) | **248** |
| pooled and injected, per the paper | **226** (111 + 115) |
| per-LNP delivery values the design supports | **≈1,356 + 666 + 690** across Fig 5f, Fig 6b, SI S5f, SI S6c–e |
| per-LNP delivery values **attributable** | **0** |

Every per-LNP panel (Fig 5f: 111 LNPs × 3 cell types × 2 genotypes; SI S5f: 115;
SI S6c–e: 226) is an **unlabelled scatter** — no per-point label, no ordering
key, no legend mapping a point to an LNP number. And **which members were pooled
is never published**: the SI formulation tables do not mark the 111 of 128 or the
115 of 120 that were injected. Library 2 makes that unrecoverability concrete —
only 110 of its 120 printed diameters are even inside the paper's own 20–200 nm
band, against the 115 it says were pooled (`DECISIONS.md` D-5).

So even granting the endpoint, a digitised point could not be attached to a
formulation. The 398 attributable values in `05_values.csv` are diameters and
formulation metadata; **0 of them are a delivery endpoint**.

*(Note on the 4,242 design-supported total in `05_values_summary.json`: it
double-counts panels that present the same underlying measurements — Fig 4d's
128 diameters and the SI p14 table's 128 are the same numbers. It is an upper
bound on panel cells, not a count of distinct measurements.)*

---

## 7. Deposit / repository sweep — EXPLICIT NULL RESULT

| searched | chars | hits |
|---|---|---|
| full JATS (`PMC6426696.jats.xml`) | **165,357** | **0** |
| full SI text (19 pages) | **25,876** | **0** |

Patterns: `GSE[0-9]+`, `SRP[0-9]+`, `PRJNA[0-9]+`, `SRX[0-9]+`, `SRR[0-9]+`,
`\bSRA\b`, `accession`, `deposit`, `github`, `zenodo`, `dryad`, `figshare`,
`osf.io`.

`find_repo.py` was run for **both** authors, as instructed, because a wrong
author term returns a clean-looking zero:

- `--author Sago` → 0 candidates (2 unrelated GitHub repos, both filtered out)
- `--author Dahlman` → 0 candidates, 0 repo hits

Methods, "Data Access", verbatim: *"The data, analyses, and scripts used to
generate all figures in the paper are available upon request to
james.dahlman@bme.gatech.edu."*

**Tier 3 of the intake gate does not exist for this paper.**

---

## 8. Corpus overlap

All queries in `stage02_corpus_overlap.py`; **every one carries a positive
control**, because a sweep that returns 0 for a broken reason is
indistinguishable in a report from a sweep that returns 0 because nothing is
there.

### Positive controls (all pass)

| control | result | passes |
|---|---|---|
| name scan, token `DSPC` | **7,038** rows | ✔ |
| InChIKey join, `DLin-MC3-DMA` | **3,098** rows, under names `DLin-MC3-DMA` + `Dlin-MC3-DMA` | ✔ |
| SMARTS sweep, tertiary aliphatic amine | **18,472** structures | ✔ |
| SMARTS sweep, sterol skeleton | **63** structures | ✔ |

### The sterol control was wrong first, and that is the finding

The first version hand-wrote a steroid SMARTS. It returned **6** corpus
structures and reported `passes: True` — while **not matching the registry's
cholesterol at all** (`HasSubstructMatch → False`). A non-zero count made a
mis-specified control look sound. This is F-070 one step further on: **a positive
control needs a positive control.**

Fixed by deriving the pattern from the registry cholesterol's own Murcko
scaffold with bond orders generified (`MakeScaffoldGeneric` writes explicit
single bonds, which cannot match cholesterol's own Δ5 alkene — the guard caught
exactly that), and by making every control **assert its own premise**: the sweep
now raises rather than reporting any count if a control's pattern does not match
the molecule it names. Result: **63** structures, premise verified.

### Token-pattern test (mandatory)

| token | substring rows | distinct structures | what it actually matches |
|---|---|---|---|
| `7C1` exact name | **0** | — | nothing |
| `PEI600` / `PEI` exact | **0** | — | nothing |
| `7C1` substring | **36** | **35** | `117C14`, `A7C10`, `B7C10` — Anderson-lab lipidoids; **none is 7C1** |
| `7` substring | **4,005** | — | anything with a 7 in its name |
| `C14` substring | **1,048** | — | `112C14`, `117C14`, `200C14` — lipidoid **tail lengths** |
| `C18` substring | **581** | — | `103C18`, `117C18`, `63C18` — same |

`C14` and `C18` are **real tokens printed in this paper's own tables**, where
they name the **PEG alkyl tail**. In the corpus they name a lipidoid tail. Token,
match and row are all real and the join is still meaningless — the
wrong-namespace failure, not a numbering failure. The substring result for `7C1`
reproduces the brief's measured 36 rows / 35 structures exactly.

**Attestation is structural, not nominal:** 7C1 cannot be joined on InChIKey
because it has no InChIKey. That is a stronger statement than "absent".

### Substructure sweep — labelled a SUPERSET

**19,962** distinct component structures scanned, **0** unparseable.

| query | structures | |
|---|---|---|
| amino-alcohol from epoxide ring-opening (7C1's chemotype) | **1,535** | superset |
| same, tightened to ≥5 further tail carbons | **988** | superset |
| PEI-like run of three ethylenimine nitrogens | **3,128** | superset |

**Claimed as nothing more than a chemotype count.** A chemotype match is not the
same molecule, and 7C1 has no structure, so it cannot be shown identical to
anything.

### In-vivo shared chemistry — measured this run, not remembered

The in-vivo pool has **exactly 2** structures shared across experiments:

| InChIKey | names | sources |
|---|---|---|
| `NRLNQCOGCKAESA-KWXKLSQISA-N` | `DLin-MC3-DMA`, `Dlin-MC3-DMA` | 6 |
| `BGNVBNJYBVCBJH-UHFFFAOYSA-N` | `SM-102` | 4 |

**This paper adds zero.** It cannot join that set even in principle, and it emits
no in-vivo row. It neither adds a third shared in-vivo structure nor raises the
source count on either existing one.

Both excipients the paper draws **are** already in the corpus in bulk —
cholesterol **30,397** rows, DSPC **7,038** rows — which is why the drawn
structures offered an anchor target at all (§10).

---

## 9. Benchmarks — the question is moot

**This paper uses no benchmark ionizable lipid.** There is no MC3, no SM-102, no
ALC-0315, no C12-200 anywhere: all 248 formulations are 7C1. The negative control
is **naked DNA barcode**, not a benchmark LNP, and it ranks below the LNPs as
expected (Fig 4e–f, SI Fig S4c–d).

So "do benchmarks top the screen?" has no referent here, and there is no
selected-winner tautology to worry about either — no second screen selects a
winner from a first.

---

## 10. Anchor posture — recorded BY HAND, as instructed

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run** (out of scope for this task).

**Net C10 result for this folder: ZERO independent anchors — honestly zero.**

`reference_structures.csv`, header exactly
`lipid_name,reference_SMILES,source_figure,read_on,read_by,expected_mismatch,note`:

| lipid | `reference_SMILES` | `expected_mismatch` |
|---|---|---|
| `7C1` | **empty** | `no_structure_is_drawn` |
| `Cholesterol` | registry string | `anchor_is_a_registry_pin_NOT_independent` |
| `DSPC` | registry string | `anchor_is_a_registry_pin_NOT_independent` |
| `C14PEG2000` | **empty** | `block_length_is_a_distribution` |
| `C16PEG2000` | **empty** | `block_length_is_a_distribution` |
| `C18PEG2000` | **empty** | `block_length_is_a_distribution` |
| `C14PEG3000` | **empty** | `block_length_is_a_distribution` |

**Why the two registry pins are NOT counted.** They come from
`provenance/reagent_registry.py`, the same source the corpus uses, so each would
compare a value to itself — the F-013 shape. Declaring that is the point;
`CLAUDE.md` §6 records C10's independence being asserted twice and false twice.

- **Cholesterol** was deliberately **not hand-written**: a hand-built
  cholesterol has already come out a diastereomer with 4 of 8 centres inverted —
  identical formula, identical exact mass to six decimals, identical skeleton
  InChIKey block. The correct response to that is to discard the anchor, not to
  flip tokens until it matches.
- **DSPC** could be built by `assemble.esterify` from stearic acid +
  sn-glycero-3-phosphocholine, a genuinely different code path. But the
  headgroup fragment written from the chemical name,
  `COP(=O)([O-])OCC[N+](C)(C)C`, is a **28-character substring of the registry's
  DSPC string**, so the template-provenance guard correctly refuses to certify
  the result as independent — and writing sn-glycerol's one stereocentre by hand
  is the sterol gamble with one centre instead of eight.
- **No bond trace** of either drawn structure was attempted. Stated as a
  decision (`DECISIONS.md` D-4), not as an impossibility.

### The guards, and making them fail on purpose

**Byte-identity is not the signal and testing for it is a trap** (F-037): a
correct re-derivation is byte-identical to a correct enumeration once
canonicalised. So:

1. **AST-walking independence guard** — walks this module's own AST at run time,
   inspecting **`Import` / `ImportFrom` / `Call` / `Attribute` nodes only**,
   never string literals, so the ban list cannot match itself (one sibling's
   guard failed on its own module for exactly that reason, and another's
   substring scan flagged its own docstring). It reports
   `make_reference.py` as `is_independent = false`, detecting
   `reagent_registry` and the `REGISTRY` attribute load — an **honest
   self-report**, which is the whole purpose.
2. **Template-provenance guard** — no string literal of ≥25 chars may be a
   substring of any target structure. **0 violations** over 32 long literals.
3. **Made to fail on purpose**, in child processes:
   `python3 -I make_reference.py --selftest` → **3/3 credited**: a seeded
   `reagent_registry` import is reported non-independent; a seeded 40-char
   literal copied out of the cholesterol string is flagged; and a **clean
   control comes back clean on both**, without which the guards would be
   always-on and certify nothing.

### What the drawing DID let us check: a number drawn and never computed

The four PEG-lipids are drawn with explicit polyether subscripts — located
**geometrically** (a smaller-font numeric span within 6 pt right of a `)` body
glyph), not by grepping the page for plausible integers.

| drawn label | drawn repeat | nominal MW | predicted `MW/44.0526` | residual |
|---|---|---|---|---|
| `C14PEG2000` | **45** | 2000 | 45.40 | −0.40 |
| `C16PEG2000` | **45** | 2000 | 45.40 | −0.40 |
| `C18PEG2000` | **45** | 2000 | 45.40 | −0.40 |
| `C14PEG3000` | **67** | 3000 | 68.10 | −1.10 |

Perfectly consistent within each MW class. **Ratio test, the sharper form**
(it cancels the end-group mass and the vendor's rounding of 1982 → "2000"):
drawn **67/45 = 1.48889** vs nominal **1.5** → **0.74% relative residual**.

**Not a tautology:** the subscripts live in the vector drawing, the PEG MW column
lives in the table, neither was derived from the other here, and the paper never
prints a PEG molecular weight computed from a repeat count.

And a second drawn-vs-table tie: library 1's **4 drawn PEG-lipids exactly equal
its 4 table combinations**, at **32 formulations each = 128**, matching the
drawn-but-uncited caption *"Library 1: 128 LNPs (4 PEG-lipids)"*. Library 2 uses
a **fifth** combination, `C18PEG3000`, that the paper never draws.

---

## 11. Checks, with numbers

| check | result |
|---|---|
| library row counts vs paper | **128/128** and **120/120** ✔ |
| distinct `Lipomer Type` over 248 rows | **1** (`7C1`) |
| positive control: `Diameter` varies | **100** distinct values ✔ |
| **mol% sums to 100** | **128/128** + **120/120** = **248/248**, 0 offenders, 0 unparsed |
| mol% tautology check | not a tautology — all percents printed independently, none derived here |
| inclusion criterion, library 1 | **119** in 20–200 nm vs **111** stated pooled → superset margin **+8** ✔ |
| inclusion criterion, library 2 | **110** in band vs **115** stated → margin **−5** ✘ (D-5, OPEN) |
| normalisation stage A | **21/21** exact, max residual **0.0** |
| normalisation stage B | **18/18** exact, max residual **0.0** |
| measured sensitivity limit | ±0.5 pp (the authors' own rounding) — stated, not hidden |
| drawn PEG repeats, ratio test | **0.74%** relative residual |
| library 1 drawn = table | **4/4** ✔, caption count matches ✔ |
| accession / repo sweep | **0** hits over **165,357** + **25,876** chars |
| `find_repo` both authors | **0** / **0** candidates |
| corpus exact name `7C1` | **0** rows |
| in-vivo shared structures | **2** (MC3, SM-102); this paper adds **0** |
| anchors | **7** rows, **2** with a structure, **5** empty on purpose, **0 independent** |
| anchor guard selftest | **3/3** credited (2 seeded + 1 clean control) |
| **seeded mutations** | **10/10 credited**, folder artifacts untouched ✔ |
| panels inventoried | **44**, all with `corpus_admissible` + reason; **0 yes / 44 no** |
| contradictions | **7** recorded (**5** open, **2** retired with the measurement that retired them) |
| **rows emitted** | **0** — `06_rows.csv` is header-only, **79** columns read live |

### Seeded mutations — 10/10, and one seed was wrong first

`seed_checks.py`. **Writes only to `tempfile.mkdtemp()`** (F-069), one child
process per seed, verdicts from `returncode` and stdout — never `python3 … |
tail`, which returns tail's exit code. Tracked-artifact sha256 hashes are
compared across the whole suite, which is the only check that sees a write from a
run that exited non-zero.

**A seed that only raises is not credited.** The harness requires a
`METRIC <name> <value>` line in **both** the clean and the seeded run and
refuses to credit a seed whose metric line is absent. That caught a real defect
in my own harness: the first stage-A seed inflated a **raw** count, which trips
`check_normalisation`'s printed-total guard and **raises before the metric is
reached** — reported `NOT CREDITED: the seeded run printed no metric line`.
**I left the guard alone and retargeted the seed** to a printed percentage, then
credited the guard with a seed of its own (`printed_total_inconsistent`,
0 → 1). Widening the gate instead is the documented wrong move.

| seed | metric | clean → seeded |
|---|---|---|
| `molpct_perturbed` | `molpct_lib1_ok` | 128 → 127 |
| `molpct_unparseable` | `molpct_lib1_unparsed` | 0 → 1 |
| `lipomer_type_overwritten` | `n_distinct_ionizable` | 1 → 2 |
| `diameter_column_blanked` | `diam_positive_control` | 1 → 0 |
| `printed_pct_perturbed` | `normA_exact` | 21 → 20 |
| `printed_total_inconsistent` | `normA_total_guard_raises` | 0 → 1 |
| `stage_B_skips_input_normalisation` | `normB_exact` | 18 → 1 |
| `peg_repeat_shifted` | `peg_all_within_one` | 1 → 0 |
| `peg_repeat_mw_classes_merged` | `peg_ratio_residual_pct` | 0.7407 → 33.3333 |
| `drawn_peg_label_dropped` | `peg_n_drawn` | 4 → 3 |

---

## 12. Source integrity (F-065)

`source/` holds exactly one pre-staged file:

```
nl8b03149_si_001.pdf
sha256 b08cbc1705ef96417b062eb63aa82f902f52d6cdcca903698f0afaa7e13e3f99
13,069,079 bytes, 19 pages
```

**Nothing was copied into `source/`.** No glob copy out of the shared
scratchpad was performed at any point; the two figure JPEGs that
`pmc_fetch.py` downloaded stayed in the scratchpad and were never used for a
value.

**Attestation is by content, not by file count or size** — those would have
passed on 17 figures from other papers. Three in-file ties, each re-checkable:

1. every page carries the watermark
   `Downloaded from pubs.acs.org/nalefd/article-supplement/739037/pdf/nl8b03149_si_001/`,
   and `nl8b03149` is the ACS manuscript id inside this paper's own DOI
   `10.1021/acs.nanolett.8b03149`;
2. SI page 1, line 2 is the paper's title verbatim;
3. the SI's FACS panel (CD31, CD45, CD68, CD11b) is the exact four-clone set in
   this paper's Methods, and SI Fig S1a's 91-nt barcode contains
   `GCTCTCATACGAACTCGTCC`, which is primer `D` in the SI's own primer table on
   p5.

**Stated gap, not papered over:** the sibling trick of sha256-matching the staged
file against the file inside Europe PMC's own `supplementaryFiles` ZIP could
**not** be done — that endpoint returned **HTTP 500** for `MED/30216729`, as did
`fullTextXML` for `PMC6426696`. The JATS was obtained via `efetch` instead
(`provenance/tools/pmc_fetch.py`), after confirming a real `pubmed_pmc`
**self-link** (→ `6426696`), not merely `pubmed_pmc_refs`.

All downloads were validated by **magic bytes**: both figure JPEGs begin
`FF D8 FF E0 ... JF`. Two PMC fetch attempts returned HTML error pages and were
correctly rejected by the tool rather than saved under a `.jpg` name.

---

## 13. Licence — NOT CC

Off the Europe PMC record, not the landing page:

- `isOpenAccess`: **`N`**; `license`: **`null`**
- PMC permissions block, verbatim: *"This file is available for text mining. It
  may also be used consistent with the principles of fair use under the
  copyright law."*

This is an **NIH author manuscript** (`NIHMS991269`) deposited under the
public-access policy; the publisher of record is ACS (*Nano Letters*),
subscription. **"Text mining / fair use" is not a CC licence**, so it fails the
`CLAUDE.md` §8 CC0/CC-BY-compatible requirement for corpus release. Moot here,
since zero rows are emitted.

**Repo licence: none, because no repository exists** (§7).

---

## 14. What I could NOT verify

- **The 7C1 structure.** Unverifiable in principle from this paper: a
  polydisperse branched polymer with an unstated degree of substitution. No
  formula, no HRMS, no NMR.
- **Byte-provenance of the staged SI against Europe PMC's own ZIP** — both
  relevant Europe PMC endpoints returned HTTP 500. Three content ties were used
  instead (§12).
- **Which LNPs were pooled.** Never published for either library. This is why
  attributable per-LNP delivery values are **0** and why D-5's library-2
  discrepancy cannot be resolved.
- **Connectivity of the drawn cholesterol and DSPC.** Not traced, by decision
  (D-4), so this folder's anchors are registry pins and are declared as such.
- **Whether the 13th cell type in Fig 2 is a spleen population.** SI Fig S1f
  plots 12 labels where the text says 13; recorded in `08_contradictions.csv`.
- **The `N:P` ratio / nucleic-acid charge ratio.** Not stated anywhere in the
  paper. Not reconstructed, not guessed.

---

## 15. Files

| file | what |
|---|---|
| `00_value_sources.json` | **the intake gate, written first** |
| `01_blocks.csv` | the 248-formulation combinatorial design |
| `01_si_tables.json` | library tables, mol% sums, inclusion criterion, normalisation residual |
| `02_corpus_overlap.json` | token test, structural attestation, superset sweep, in-vivo shared |
| `03_enumerated.csv` | one row: `7C1`, empty SMILES, `confidence = low` |
| `03_drawn_checks.json` | drawn PEG repeats, drawn-vs-table, vector-geometry posture |
| `04_anchor_qc.json` | anchor posture + both guard reports |
| `05_values.csv` | **44 panels**, `corpus_admissible` + reason on every row |
| `05_values_summary.json` | the counts, made to sum |
| `06_rows.csv` | **header-only**, 79 columns read live from `data_vivo.csv` |
| `07_qc.json` | every check, with numbers |
| `07_seed_results.json` | 10/10 seeded mutations + hash comparison |
| `08_contradictions.csv` | 7 contradictions, retired ones keep the measurement that retired them |
| `DECISIONS.md` | D-1 … D-6, all `STATUS: OPEN` |
| `make_reference.py` | anchor builder + AST and template guards + `--selftest` |
| `reference_structures.csv` | 7 rows, 5 empty on purpose, **0 independent** |
| `stage01_si_tables.py`, `stage02_corpus_overlap.py`, `stage03_drawn_checks.py`, `stage05_values.py`, `stage07_qc.py`, `seed_checks.py` | the pipeline |

**Nothing outside this folder was modified.** `curation/lib/` was not edited, the
shared Playwright browser was not used, `provenance/structure_check.py`,
`check_anchor_independence.py` and `check_anchor_coverage.py` were not run,
nothing was merged into `data/`, and no script in this folder reads an input from
a session scratchpad.
