# PROVENANCE — maohq_39375194

Cheng L, Zhu Y, Ma J, Aggarwal A, Toh WH, Shin C, Sangpachatanaruk W, Weng G,
Kumar R, Mao H-Q. *Machine Learning Elucidates Design Features of Plasmid
Deoxyribonucleic Acid Lipid Nanoparticles for Cell Type-Preferential
Transfection.* ACS Nano 2024;18(42). PMID **39375194**, PMCID **PMC11512640**,
DOI **10.1021/acsnano.4c07615**.

Full published title confirmed from Europe PMC's REST `core` record in one call.
That record also settles that this is an original-data paper by the stated
discriminator: `abstractText` is **1,770 characters** (not 0) and the JATS full
text carries METHODS and RESULTS AND DISCUSSION with named subsections —
`pubType` was not used.

**Endpoint.** In vitro firefly-luciferase reporter expression from **plasmid
DNA** LNPs, read as RLU on a luminometer **48 h** after dosing, in six monolayer
cell lines; stored by the authors as **ln(mean RLU)**.

**Cargo, exactly as published:** *firefly luciferase pDNA (Aldevron)*. `Cargo =
pDNA`. Not pooled with mRNA rows.

---

## What this paper is

| | |
|---|---|
| distinct ionizable-lipid **structures** | **1** (DLin-MC3-DMA) |
| distinct helper-lipid structures | 6 |
| **formulations** | **1,152** = 1,080 training + 72 de-novo validation |
| **individually measured** transfection values | **6,192** = 5,760 + 432 |
| of those, already in the corpus | 1,080 (`RLU_HepG2` ≡ `YZ_2022`) |
| emitted | **5,112** |

The 1,080 rows are **formulations of one lipid**, not 1,080 lipids. The one-hot
group is the *helper* lipid.

---

## Value sources

Every emitted value is **tier 3 — the authors' own repository**. No figure was
digitised; no page was rendered for a value. `00_value_sources.json` holds the
full gate, including why tiers 1, 2 and 4 could not supply values and why tier 5
was not permitted.

```
github.com/MaoResearchGroup/ml_lnp_design_rules
commit 82cee58057a6168b61739c5863c6ac8a4ff19a7e   2024-10-31 14:22:46 -0400

Raw_Data/1080_LNP_OHE.csv          md5 98baa22d8d8e9a92648e1c9b9904198f
Raw_Data/Compiled_validation.xlsx  md5 5e1ea0cd1e8f7c114ff355d58c4206b4
source/nn4c07615_si_001.pdf        sha256 27ef09af…c5f39e8a
```

**Untrusted-data posture.** The repository was cloned `--depth 1` into its own
new, empty directory in the session scratchpad — outside this folder. Nothing in
it was executed; `utilities.py` was *read* to establish how the deposited values
relate to the paper's text. The two data files were copied into `source/github/`
and verified byte-identical by md5. Every script in `scripts/` lives in a
different directory from the clone, takes its input paths as arguments, and is
run with `python3 -I`.

The SI is rendered only in `renders/`, and only to read **structures** — each
image extracted at its native embedded resolution (220–336 dpi, measured, not
assumed) and paired with its row label by y-overlap geometry rather than by
assuming table order.

---

## C10 anchor posture — recorded by hand

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` write tracked files other agents read, so **none of
them was run**. The local equivalent is `scripts/stage03b_c10.py`, whose output
is `04_c10_anchor_check.json`.

| anchor | chain lengths | posture | result |
|---|---|---|---|
| DLin-MC3-DMA | C37, 4×C=C | **DERIVED** | exact match (full InChIKey) |
| DSPC | 18:0/18:0 | **DERIVED** | exact match |
| DOPE | 18:1/18:1 | **DERIVED** | exact match |
| 18PG (as drawn) | 18:1/18:1 | **DERIVED** | exact match |
| 14PA | 14:0/14:0 | **DERIVED** | exact match |
| DOTAP | 18:1/18:1 | **DERIVED** | exact match |
| DDAB | 2 × C18 alkyl | **DERIVED** | exact match |
| Cholesterol | sterol | **NO ANCHOR**, empty `reference_SMILES` | gap recorded, see below |
| DMG-PEG2000 | 45 EO | **NO ANCHOR**, empty `reference_SMILES` | deliberate, see below |

**7 of 7 anchors are exact matches on the FULL InChIKey**, 0 skeleton
mismatches, 0 stereo mismatches. Chain lengths deliberately span C14, C18
saturated, C18 mono-unsaturated and C37 tetra-unsaturated.

### How independence is established — and how it is not

All 7 anchors are **byte-identical** to the emitted structures once
canonicalised. **That is expected and is not the signal.** F-037's whole lesson
is that a correct re-derivation is byte-identical to a correct original, which
is exactly why C10's own `anchor_may_be_a_copy` missed a re-copied template.

Independence is established **by provenance, at run time**, by
`make_reference.py::ast_independence()`, which walks *this module's own parsed
AST* — not a text scan; a sibling's substring version flagged its own
docstring — and confirms:

- **no import** of `reagent_registry`, `reagents`, `corpus_components`,
  `assemble` or `structure_check` (the sources `stage03_structures.py` uses);
- **no call** to `RWMol`, `EditableMol`, `AddBond`, `RemoveAtom`,
  `CombineMols`, `ReplaceSubstructs`, `esterify`, `attach` or `SetChiralTag`.

Measured: `forbidden_imports_found: []`, `forbidden_calls_found: []`.

The **traversal is opposite**. The registry writes these lipids *tail-first* —
the sn-1 acyl chain first, head group in a branch. `make_reference.py` composes
*head-first*: from the choline nitrogen / ethanolamine nitrogen / head glycerol
/ phosphate outward, acyl chains appended last. The backbone stereocentre is
**not copied**: both tags are generated and the one whose CIP code is **R** is
selected, R being a name-derived constant (in every head group here the C2
priorities are O-acyl > CH₂-O-P > CH₂-O-acyl > H).

### The gate was made to fire, with identical formula and identical exact mass

`python3 -I make_reference.py --selftest` — **3 of 3 seeds credited**:

| seed | class | formula | exact mass | clean key | seeded key |
|---|---|---|---|---|---|
| `mc3_ester_at_C18` | regiochemistry | C43H79NO2 *(identical)* | 641.611081 *(identical)* | CLYUJYGYPALXTO | RBLCKFZMIVYNAG |
| `dspc_1_3_diacyl` | connectivity | C44H88NO8P *(identical)* | 789.624755 *(identical)* | NRJAVPSFFCBXDT | QVOSSHZMHWUKKG |
| `pg_oleoyl_delta10` | fragment direction | C42H78O10P⁻ *(identical)* | 773.533809 *(identical)* | DSNRWDQKZIEDDB | KZKARNHNCHWPRE |

**The anchor earned itself before it was trusted.** `make_reference.py`'s first
version had the MC3 branch attached by the wrong end — the fragment was written
C1→C18 and appended as a branch off C19, which bonds **C1** to C19. The result
had the correct molecular formula and the correct exact mass to six decimals and
a different InChIKey. It was caught by the anchor comparison and fixed by
re-deriving from the chemical name, not by tuning until the key matched. Logged
as `08_contradictions.csv` C-10.

### The two gaps, stated so they are visible

- **Cholesterol** — `reference_SMILES` empty on purpose. Hand-writing an
  8-stereocentre steroid in this file would be a worse anchor than none, and
  cholesterol is already covered corpus-wide by **C14**, which cleaves every
  ester in `IL_SMILES` and pins the liberated sterol's *full* InChIKey.
- **DMG-PEG2000** — `reference_SMILES` empty on purpose: the block length is a
  distribution, so no single SMILES is the reagent (`CLAUDE.md` §3). The drawn
  member is reproducible as `dmg_peg2000(45)` and *is* what `06_rows.csv`
  emits, because the corpus column carries a PEG SMILES and C17 gates one
  structure per reagent name.

---

## Checks, and the proof they can fail

Nothing here is credited on a clean run alone.

| suite | result |
|---|---|
| `scripts/stage01_invert.py` | all file checks pass, exit 0 |
| `scripts/seed_checks.py` | **13 / 13 seeds credited**, control passes |
| `make_reference.py --selftest` | **3 / 3 seeds credited** + AST independence |
| `scripts/seed_emit.py` | **5 / 5 seeds credited** (literal confidence, caveat without note, bad vocabulary, missing record, control) |
| `scripts/stage05_qc.py` | 0 failures |

A seed is credited **only when a named numeric metric moves** *and* the gate
exits non-zero, each seed in its own process with `returncode` read directly. On
the first run **two seeds moved a metric while the gate stayed silent**
(`measurement_blanked`, `values_already_minmaxed`) — a real coverage gap. The
**gate was widened, not the seed**: a structural "a helper lipid is wholly
measured or wholly absent in a cell line" invariant, and a refusal of any target
column with no value above 1.0.

One check is reported but **explicitly not credited**: the four molar percents
summing to 100. HL+IL = a and CHL+PEG = 100−a by construction, so it is 100
identically for any input. A check that cannot fail is not a check.

---

## Licences

| | |
|---|---|
| **Article** | **Closed / subscription.** Europe PMC `isOpenAccess = N`, `license = null`. The PMC deposit is an NIH author manuscript; `<permissions>` reads *"available for text mining … consistent with the principles of fair use"*. **Not CC0/CC-BY.** |
| **GitHub repository** | **GPL-3.0** (verbatim text at the repo root). A copyleft *software* licence. **Not CC0/CC-BY**, and not a permissive data licence. |

`CLAUDE.md` §8 requires CC0/CC-BY compatibility before corpus release; neither
side qualifies. See `DECISIONS.md` **D-1**. Note that `YZ_2022` (Nat Commun,
PMID 35879315), which covers this library's HepG2 arm and is already merged,
**is** open access — so the licence-clean and licence-unclear versions of one
experiment sit in different papers.

---

## Files

| file | what |
|---|---|
| `00_value_sources.json` | the intake gate, written first |
| `01_blocks.csv` | the 9-component inventory with structures and how each resolved |
| `01_checks.json` / `01_seed_results.json` | file-level checks and their 13 seeds |
| `02_corpus_overlap.json` | full-key, skeleton, name, formulation overlap; the 18PG conflict |
| `02b_duplicate_quantified.json` | the `RLU_HepG2` ≡ `YZ_2022` duplication, in ln units |
| `03_enumerated.csv` | the single ionizable lipid, with corpus-derived descriptors |
| `03_structures_qc.json` | structure resolution, incl. the 18PG head-stereo derivation |
| `04_c10_anchor_check.json` | the local C10 comparison |
| `04_emit_seed_results.json` | emit's refusals, seeded |
| `05_values.csv` | **10,368** rows — every value, admissible or not, with a reason |
| `05_value_counts.json` | the counts, asserted to sum |
| `06_rows.csv` | **5,112** emitted rows in the 79-column schema |
| `07_qc.json` | every check with numbers |
| `08_contradictions.csv` | 10 contradictions, 3 open, 7 retired with the measurement that retired them |
| `DECISIONS.md` | D-1 … D-6, all OPEN |
| `make_reference.py` / `reference_structures.csv` | the anchors |
| `renders/` | SI structure images at native resolution |
| `scripts/` | the stages; none reads from a scratchpad, none writes outside this folder |

**Not merged into `data/`.** Nothing outside this folder was modified.
