# PROVENANCE — XQ_2020 / proposed `YL_2020`

Li Y, Li R, Chakraborty A, Ogurlu R, Zhao X, Chen J, Xu Q.
*Combinatorial Library of Cyclic Benzylidene Acetal-Containing pH-Responsive
Lipidoid Nanoparticles for Intracellular mRNA Delivery.*
**Bioconjugate Chem.** 2020. PMID **32520527**, PMCID **PMC9069324**,
DOI **10.1021/acs.bioconjchem.0c00295**.

**Endpoint: intracellular GFP mRNA delivery to HeLa cells, read as the
percentage of GFP-positive cells by flow cytometry at 24 h.** The run brief
warned that a sibling paper turned out to be about antimicrobial activity; the
title and abstract were re-read from the PMC record before anything else and
this one is in scope.

---

## 1. What is here

| file | what |
|---|---|
| `00_acquisition.json` | URLs, sha256, UTC dates, the licence string, what could not be fetched and why |
| `00_value_sources.json` | **the intake gate** — why tiers 1-4 fail, panel by panel |
| `01_blocks.csv` / `_qc.json` | 11 amine heads + the O16CBA tail; the named-reagent MW ladder |
| `03_enumerated.csv` / `_qc.json` | the 11 lipidoids, by exhaustive regiochemistry enumeration |
| `04_corpus_overlap.json` | full-InChIKey / skeleton / `IL_name` comparison against `data/` |
| `05_values.csv` / `_qc.json` | 33 values, each with its pixel window and row |
| `07_qc.json` | 30 checks with numbers, 29 gating |
| `08_contradictions.csv` | 4 places the paper contradicts itself |
| `08_retracted_suspicions.csv` | 5 suspicions, each with the measurement that retired it |
| `DECISIONS.md` | 6 OPEN decisions — this folder is not merge-ready, deliberately |
| `make_reference.py` → `reference_structures.csv`, `anchor_independence.json` | the C10 anchors and their independence evidence |

Rebuild, in order:

```bash
python3 -I stage01_blocks.py
python3 -I stage03_enumerate.py
python3 -I stage05_values.py
python3 -I make_reference.py
python3 -I stage07_qc.py ; echo $?     # the gate. NOT through a pipe.
python3 -I stage07_qc.py --self-test   # every check, seeded
```

---

## 2. Where the numbers came from

**Values: Figure 5 only.** There is no Source Data file (ACS published one SI
asset for this article and the PMC record declares exactly one), no value table
in the SI (all 8 pages read as text — they are a synthetic scheme, two 1H NMR
spectra, three ESI-MS panels and a TEM plate), and no repository (`find_repo.py`
over the DOI, PMID, title and author returned only third-party scrapers). The
paper has no data-availability statement at all.

The RESULTS prose prints **seven** of the 33 numbers in words, plus two counts,
two maxima, two ranges and three pairwise-sign statements. Those are used as
**checks**, never as the extraction — covering 7 of 33 is not a value source.

**Structures: Figure 1C (all 11 amine heads drawn in full), Figure 1B and
Figure S1 (the O16CBA tail drawn in full), and Figure 2F (the molecular formula
and calculated/observed m/z of all 11 lipidoids).**

The tail is additionally pinned **without any drawing**: the Experimental
Section names every reagent with a gram and a millimole figure, so each
reagent's molecular weight is recoverable by division. Nine reagents agree with
the structures we assign to within **1.63%** (worst case: triethylamine, printed
as "2.2 g" — two significant figures). This is what fixes the hexyl chain
length; no pixel is counted anywhere in this curation for a chain.

---

## 3. Anchor posture — verified BY HAND

This run may not execute `provenance/structure_check.py`,
`check_anchor_independence.py` or `check_anchor_coverage.py`, because they write
tracked files other agents are using. The equivalent evidence was produced
inside this folder instead and is recorded here.

**Coverage.** 11 of 11 lipidoids carry an anchor. None is empty.

**What the paper actually draws, stated so the anchor is not oversold.**
Figure 1C draws all eleven amine heads in full and Figure 1B draws the O16CBA
tail in full, but **the lipidoid itself is drawn only as `R-N[tail]₂`** — R is a
letter and the arm count is a bracket subscript. So no complete lipidoid is
drawn anywhere in the paper, and every anchor here is a **composition of two
separately drawn fragments**. The anchor therefore tests that the composition is
right; it cannot test a drawn whole, because there is none. That sentence is in
the `note` column of every row of `reference_structures.csv`.

**Independence — executed, not asserted.** `make_reference.py` builds a
**string**; `stage03_enumerate.py` builds a **molecule graph**:

| | `stage03_enumerate.py` | `make_reference.py` |
|---|---|---|
| arm | derived by reducing O16CBA's alkene under a unique SMARTS match, cross-checked against an `assemble.esterify` route | hand-written from Figure 1B, ester-first |
| join | `assemble.attach` → `Chem.RWMol` + `AddBond` on atom-mapped atoms, once per arm, breadth-first over every N–H site | `str.replace` |
| heads | `01_blocks.csv`, written substituent-first | re-written amine-first, in this file; `01_blocks.csv` is never opened |
| RDKit used to build? | yes | **no** — the first RDKit call is `MolFromSmiles` on the finished string |

`make_reference.py` contains no `RWMol`, no `AddBond`, no `RemoveAtom`, no atom
map, no SMARTS, and does not import `curation/lib/assemble`.

**The independence test is on the GENERATOR'S COMPOSED TEMPLATE, not on the
emitted anchor.** Once canonicalised, a correct re-derivation is byte-identical
to a correct enumeration, so testing the emitted anchor for byte-difference
would flag correct work and miss the real defect — which is exactly how C10's
own `anchor_may_be_a_copy` missed F-037. Measured (`anchor_independence.json`):

- templates byte-identical to the enumeration: **0 of 11**
- worst shared literal run between a template and the enumeration: **18
  characters** of templates 102+ characters long (it is the common
  `c(OC)cc(OCCCCCC)` aryl fragment, which both files necessarily contain)
- exact **full**-InChIKey match, anchor vs enumeration: **11 of 11**

**Both guards were made to fail on purpose**, in a copy of the folder under the
scratchpad, before being trusted:

| seeded defect | result |
|---|---|
| head 78's template replaced by the enumeration's own canonical SMILES | **exit 1**, "the composed TEMPLATE is byte-identical to the enumeration for ['78-O16CBA']" |
| head 80's linker written with 2 carbons instead of 3; head 306's central methyl dropped | **exit 1**, "anchor and enumeration DISAGREE for ['80-O16CBA', '306-O16CBA']" |
| head 400's three arms written 1+2 instead of 2+1 | **exit 1** at parse — a 5-valent nitrogen |
| **positive control**: 113's and 306's 2+1 split written on the *opposite* nitrogen | **exit 0**, 11/11 exact — confirming by a second route that the two placements are the same molecule |

That last row matters: it is independent confirmation of the symmetry argument
that makes the three-armed members' regiochemistry unique, rather than my
assertion of it.

---

## 4. The two chemistry traps the brief named, and what settled them

**Acetal regiochemistry — which diol oxygens close the ring.** 1,1,1-tris(hydroxymethyl)ethane
has three CH₂OH groups and only two close the 1,3-dioxane. Molecular formula is
blind to the choice. **Settled by computation, not by claim**: the three CH₂OH
groups have identical canonical rank (`breakTies=False` → `[1, 1, 1]`), so all
three pairings give one constitution. `07_qc.json` P10b, seeded as `acetal_sym`.

**Arm placement on the multi-nitrogen heads.** Heads 113, 306 and 400 carry three
nitrogens each and take three arms. Rather than choose, `stage03_enumerate.py`
**enumerates every placement**: breadth-first over all N–H sites at every step,
with the arm count solved exactly from the formula Figure 2F prints (no
remainder tolerated). The branches reconverge — 113 and 306 expand 1 → 2 → 1,
400 expands 2 → 2 → 1 — and **all eleven lipidoids end with exactly one distinct
constitution**. There is nothing to flag because there is nothing to choose.

**The fragment was not built by string surgery on a canonical SMILES.** The run
brief names a sibling that did exactly that and mislocated every ester while both
HRMS anchors still passed at 0.3-0.5 ppm. Here the Michael arm is built twice, by
routes that share no code — reducing O16CBA's alkene under a unique SMARTS match,
and `assemble.esterify`-ing propanoic acid onto HexDMBAH — and the two must agree
or stage 03 raises. The attachment atom is atom-mapped in both.

---

## 5. What is NOT verified

- **Connectivity, by any mass check.** All 33 mass comparisons in `07_qc.json`
  pass on any isomer. 77-O16CBA and 82-O16CBA are an in-paper demonstration of
  this: same formula, same printed calculated *and* observed m/z, different
  molecules. Check P10e asserts our two SMILES differ while their masses agree.
- **Stereochemistry, twice.** See `DECISIONS.md` D-2 and D-3. The second half is
  the one worth repeating: **`curation/lib`'s stereo guard cannot see the acetal
  pair at all**, so `check_stereo_preserved` passing is not evidence here.
- **That the x-axis group labels were read correctly.** The order
  75/76/77/78/80/81/82/93/113/306/400 is read off a 1050 px raster by eye. It is
  constrained hard — five groups are identified by name in the prose and the rest
  is pinned by eleven sign statements and two counts, and the `swap_groups` seed
  shows a 113↔306 transposition fires P4 — but no check reads the glyphs.
- **The in-paper ESI-MS spectra themselves.** Only the m/z *labels* were read,
  as PDF text spans. No peak was located in a spectrum.
- **Anything about the nanoparticles.** Size, zeta, encapsulation and the
  degradation kinetics were not extracted; they are not delivery endpoints.

---

## 6. Corpus overlap (report only — a merge-level fact, not acted on)

**Zero overlap.** None of the eleven R-O16CBA lipidoids matches a corpus lipid by
full InChIKey, by skeleton InChIKey, or by `IL_name`, across 19,919 unique
structures in `data_vitro.csv` + `data_vivo.csv` (read 2026-11-10).

**But the head groups are shared, and that is invisible to the above.** These are
the Xu-lab/Tufts amine catalogue numbers, and four of them — **75, 80, 113, 400**
— are already in the corpus as the heads of `QL_2026`'s bispecific-linker library
(160 rows). Our structures for those four, read off Figure 1C, are byte-identical
to `QL_2026`'s, which were taken from a vendor catalogue by chemical name. That
is an **external cross-check on four of eleven heads by an independent curator on
a different paper** (`07_qc.json` P11) — and it is also the only chemistry this
paper shares with any other source.

**The name-token trap.** 163 corpus `IL_name`s begin with one of our eleven amine
numbers and they do **not** all mean the same thing: `QL_2026`'s `113-L1-C6` uses
113 as the Xu amine number and does contain that head, while `BL_2024`'s `113-9`,
`113-22` and `113-23` use 113 as an unrelated series index and do not. Keying on
the leading token across experiments would merge them. `04_corpus_overlap.json`
has the per-number breakdown.

---

## 7. Results in one table

| | |
|---|---|
| distinct ionizable structures | **11** |
| formulations | **33** (11 lipidoids × {no helper lipid, chol/DOPE 4/1/1, chol/DOPE 4/1/4}) |
| emittable rows | **33** |
| rows with a caveat (`emit.caveat_count`) | **3** — amine 77's three rows |
| checks | **30**, of which **29 gating**; **29 pass**, 1 non-gating failure (P6b = contradiction C-1) |
| gate exit code | 0 clean, 1 on every seeded defect — verified directly, not through a pipe |
| seeds | 22, each moving the agreement count; all 29 gating checks covered |
| C10 anchors | 11/11, exact full-InChIKey match, 0 copied templates |
