# PROVENANCE — proposed `XZ_2020`

Zhao X, Chen J, Qiu M, Li Y, Glass Z, Xu Q.
*Imidazole-Based Synthetic Lipidoids for In Vivo mRNA Delivery into Primary
T Lymphocytes.*
**Angew. Chem. Int. Ed.** 2020. PMID **32662132**, DOI **10.1002/anie.202008082**.
PMCID: **none** — no author manuscript was ever deposited.

**Endpoint: mRNA delivery.** FLuc luminescence (RLU) and percent EGFP-positive
cells in primary human CD8+ T lymphocytes ex vivo; FLuc IVIS and Cre/tdTomato
recombination in mouse spleen in vivo. Verified from the Europe PMC
`resultType=core` record **before any chemistry** (953-character abstract,
`pubType` = Journal Article + Research Support N.I.H. Extramural, an 11-page SI
with its own Experimental Procedures). It is original research, in scope.

---

## 1. The headline

**The values this paper exists for are not reachable.** The screen — 18 amine
heads against up to 8 tails — is plotted only in the main-text figures, and the
article is subscription-only with no deposited author manuscript, no PMC record,
no Europe PMC full text, and no open copy in OpenAlex or Unpaywall. All three
Wiley routes return HTTP 403 to a plain client. No bot check was met and none
was bypassed; no login was attempted; the shared browser was not used.

What we hold is the Supporting Information, which contains delivery values for
**two** lipidoids and a full synthesis for **five** tails.

| | count |
|---|---|
| amine heads NAMED (Fig S4/S5 axes) | 18 |
| tail labels NAMED anywhere | 11 |
| tails STRUCTURALLY SPECIFIED (SI sections 1.1–1.4) | 5 |
| amine heads ENUMERABLE (93, 9322 — both from corpus attestation) | **2** |
| amine heads NOT enumerable (no attestation, grammar unverified) | **16** |
| ionizable structures BUILT | **10** |
| structures with a reachable DELIVERY value | 2 |
| structures with **both** a structure and a value | **2** |
| delivery datapoints individually measured, in reach | 38 |
| corpus rows emitted | **13** |

## 2. Files

| file | what |
|---|---|
| `00_acquisition.json` | what was fetched, sha256, every route that failed and its status code |
| `00_value_sources.json` | **the intake gate**, written first; why tiers 1–3 fail panel by panel |
| `01_blocks.csv` / `_qc.json` | head 93 + five tails, each built from the named reagents of its own synthesis; the name-token rule |
| `02_corpus_overlap.json` | full / skeleton InChIKey, `IL_name`, substructure and Tanimoto NN against `data/` |
| `03_enumerated.csv` | the five lipidoids |
| `04_library_membership.csv` / `.json` | the library's head and tail labels; per-head membership **refused**, with the numbers |
| `05_values.csv` | 44 values with exact source location, readout, units, cargo, cell, dose, timepoint, assay |
| `09_figS1_mass_test.json` | Figure S1 base-peak m/z by axis geometry; the head-9322 confirmation and the four refusals |
| `05_figS2_geometry.json` | every tick position, calibration, bar window, top row and flank delta |
| `06_rows.csv` / `06_skipped.json` | 13 corpus-schema rows; 31 skipped, each with a reason, counts summing to 44 |
| `07_qc.json` | **14** checks with numbers, all gating |
| `08_contradictions.csv` | 10 places the paper contradicts itself |
| `08_retracted_suspicions.csv` | 6 entries — 5 suspicions retired by measurement, and **R-06, a retraction of one of my own reported findings** |
| `DECISIONS.md` | **9 OPEN decisions** — this folder is not merge-ready, deliberately |
| `make_reference.py` → `reference_structures.csv` | the C10 anchor and its own gate |

Rebuild, in order:

```bash
python3 -I stage01_blocks.py
python3 -I stage03_enumerate.py
python3 -I stage02_overlap.py
python3 -I stage04_library.py
python3 -I digitize_figS2.py
python3 -I measure_figS1.py
python3 -I stage05_values.py
python3 -I stage06_rows.py
python3 -I make_reference.py
python3 -I stage07_qc.py ; echo $?          # the gate. NOT through a pipe.
python3 -I stage07_qc.py --self-test        # every check, seeded
```

## 3. Where the numbers came from

**Structures: the SI's own Experimental Procedures, tier 4 (PDF text).** Each of
O17O, O17S, O17Se, O16S-S and N16S-S is built reagent by named reagent from its
own synthesis. Nothing structural is read off a figure, because the SI draws no
structure anywhere.

**Neither head is from this paper.** Both are resolved from the corpus through
`curation/corpus_components.py`'s guarded lookup, which **raises** on a name
carrying more than one structure instead of picking by row order (error class
30, F-014/F-044):

* **93 = 1-(3-aminopropyl)imidazole** `NCCCn1ccnc1` — `IL_head_name` `93`
  (`AA_2008`, 11 rows) and `A93` (`KW_2014`, 5 rows) resolve to the *same single*
  structure, two spellings in two experiments; `curation/xuq_32520527` read it
  off a drawn figure in a third Xu-lab paper and agrees.
* **9322 = 1-(3-aminopropyl)-2-methylimidazole** `Cc1nccn1CCCN` — `IL_head_name`
  `A9322` (`YR_2024`, 22 rows), one structure, **and confirmed against this
  paper's own Figure S1** (below).

Every emitted row carries `confidence = medium` and a note saying the head comes
from outside this paper. D-6 is OPEN.

**The other sixteen heads are NOT enumerable.** `9310`–`9316`, `9321`, `9323`,
`9324`, `9331`–`9334`, `9341`, `9351`, `9352` have **zero** corpus attestation —
no row in either pool carries any of those tokens in `IL_name` or
`IL_head_name`. One data point (`9322` = the 2-methyl analogue of `93`) does not
define a digit grammar, and the paper states none. **The grammar is unverified
and no structure is invented.**

**Values: Figures S2 and S3 only, tier 5.** Figure S3's ten percentages are
*printed numerals* at ~22 px glyph height — transcribed, not estimated. Figure
S2's 33 bars are measured through `curation/lib` (`axis.calibrate` +
`bars.bar_top`); no geometry was written in this folder.

**Figure S1: the glyph read is still refused; the panel is not.** Its base-peak
m/z labels are 5–7 px tall and remain illegible at 28× LANCZOS — that stands.
What the first pass got wrong was treating an unreadable label as an unusable
panel, and concluding that because the axis cannot *pin* a formula it cannot
*test* one. Those are different questions.

With head 9322 recovered, `measure_figS1.py` measures the base-peak **position**
against the m/z axis (ticks fixed by a 3-digit/4-digit label-width signature, not
by reading them) at **3.93 Da/px**:

* **9322-O17O measured 765.5 against a predicted 764.651 — +0.8 Da = +0.20 px.**
* Head discrimination: the plain-imidazole head 93 predicts 750.635, **+3.78 px**
  away. The panel picks the 2-methyl head by **18.9×** on a one-CH₂ difference.
* The three candidate labels the first pass could not separate: **764.73 at
  0.20 px**, 784.73 at 4.89 px, 794.75 at 7.44 px. Decided.
* 9322-O17C, which has no candidate tail, yields a **derived per-arm mass of
  310.19 ± 1.99 Da** — excluding a carbonate (−74.1 Da) and consistent with a
  plain C17 alkyl acrylate (−0.1 Da) without pinning it. D-4 stays OPEN.

**4 of 6 panels refuse** the tick-label calibration and were **not** rescued by
loosening a threshold (O17S misses at 4.65 % spacing nonuniformity against 4 %,
O16B at 18.3 %, EC18 at 37.8 %, N16B never reaches 10 label boxes). What remains
true: this axis cannot separate O from S (16 Da ≈ 4 px) and cannot separate
O16B from N16B (1.97 Da = 0.47 px), so **D-3 is still OPEN**.

## 4. The checks that carried weight

All 14 gate; `--self-test` seeds each one and every seed **moves its agreement
count** (none is a no-op that merely raises).

| check | number |
|---|---|
| tail name-token rule | **5/5** tails: the number in the name equals the heavy-atom count of the chain off the ester O / amide N, and O↔ester, N↔amide, with zero free parameters |
| SI reagent molar-mass ladder | 16/18 within 5 %; the **two intermediates** at **0.25 %** and **0.17 %**, which is what pins C14 and C12 |
| repeated-formulation identity | **3/3** — 16:4:1:1 is drawn three times and three independent window detections and three independent `bar_top` calls return top row **274.0**, spread **0.000 px** |
| untreated viability = 100 | read **99.65 %**, error **0.35 pp = 1.52 px**, against a round number the authors drew and never printed |
| printed optima reproduced | 2/2 — the SI's "30 µg/mL" is Fig S2d's argmax; its "~6 h" is consistent with Fig S2c (no 6 h bar; 7 h is the max) |
| axis calibration | 4/4 panels; worst fit residual 0.305 px, worst **leave-one-out** residual 0.763 px against 1.5 px |
| bar flank agreement | **33/33** bars at flank delta **0.0 px** — no s.d. cap was read as a top (error class 32) |
| corpus overlap | **0** full-InChIKey, **0** skeleton hits; near-miss 277 vitro + 25 vivo rows on the same head |
| Figure S1 confirms head 9322 | **2/2** independent tests: head 9322 vs head 93 at **0.20 px vs 3.78 px** (18.9× margin), and the glyph-candidate test picking 764.73 at **0.20 px** over 784.73 (4.89 px) and 794.75 (7.44 px) |
| C10 anchor gate | `--check` rc 0, `--seed` rc **1**, `--provenance` rc 0 |
| emit refuses a literal `confidence` | raises (F-011) |
| row counts | 13 + 31 = 44 |
| z-normalisation | mean 7.7e-08, sd 1.000000105 at 6 dp |
| flags on the row | 13/13 carry a caveat **and** a note |

## 5. ANCHOR POSTURE — recorded by hand

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run**: they write tracked files other
agents read, and this run is forbidden to touch anything outside this folder.
The posture is therefore stated here and gated locally by
`make_reference.py --check`.

* **Anchors stored:** 3 rows in `reference_structures.csv`, header exactly
  `lipid_name,reference_SMILES,source_figure,read_on,read_by,expected_mismatch,note`.
  **Two** carry a structure (93-O17S and, since the 2026-11-10 recovery,
  9322-O17S); **one carries an empty `reference_SMILES` on purpose** (the 93xx
  analogue series) so that the *absence* of an anchor is visible rather than
  silent. Both stored anchors match their enumeration on the full InChIKey.
* **What the anchor is derived from:** no figure. The paper draws nothing we can
  see. The anchor is typed out by hand from chemical names — SI 1.2 for the
  tail, SI 1.5 for the 1 : 2.4 amine-to-acrylate stoichiometry, and the
  externally established head.
* **Independence is by PROVENANCE, not by bytes.** F-037 is explicit that a
  correct re-derivation is byte-identical to a correct enumeration once
  canonicalised, so byte-identity is not the signal and testing for it is a
  trap. `make_reference.py --provenance` audits the code body and reports
  **0 occurrences** of `assemble`, `RWMol`, `AddBond`, `RemoveAtom`,
  `ReactionFromSmarts`, `AllChem`, `stage03` and `01_blocks`: the anchor SMILES
  are source string literals and the module imports no fragment joiner and no
  reaction machinery, while the enumerator builds by running two Michael-addition
  SMARTS reactions. The audit is **executed**, not asserted.
* **The gate was made to fail on purpose.** `make_reference.py --seed` moves one
  Michael arm from the amine nitrogen onto imidazole C2. Molecular formula
  `C44H83N3O4S2` and exact mass **781.582500** are identical to six decimal
  places, so no mass check can see it; the skeleton InChIKey goes
  `OLOIOPYNTGMAAP` → `GNLMIRXUIQYCQW` and the gate exits **1**. Clean run exits
  **0**. Both return codes are read from `returncode`, never through a pipe, and
  both are recorded in `07_qc.json` under `c10_anchor_gate`.
* **Anchor coverage, by hand:** 2 of 10 enumerated lipidoids carry an anchor —
  one per head. The other eight share a head and the same two-arm construction
  and differ only in a tail typed from its own SI synthesis; none of them is
  emitted as a corpus row, so the anchors cover **13 of 13 emitted rows**.

## 6. What could NOT be verified

* Everything in the main text: the library screen, the lead lipidoid, the in
  vivo figures. Subscription wall, no licensed route available to this run.
* The **sixteen** 93xx analogues other than 9322 — no corpus attestation, no
  stated digit grammar, no structure. (Head **9322 itself is now resolved**; the
  first pass's claim that it was unobtainable is retracted in
  `08_retracted_suspicions.csv` R-06.)
* The identity of `O17C` (D-4), `O16B`/`N16B` (D-3), `O18S-S`, `EC18`, `EC16`.
  For `O17C` the discriminating test — ¹³C for a carbonate near 155 ppm versus
  an ester near 173, or an O-CH₂ proton count — was **searched for by name and
  is absent from this SI**. It is not manufactured.
* The dose behind Figure S2a: the SI states 17 µg/mL and 30 µg/mL in adjacent
  sentences (C-06), so `Dose_ug_nucleicacid` is left blank on all 13 rows.
* The SI's licence string (D-9).


---

## 7. Correction log

**2026-11-10 — head 9322 recovered; one reported finding retracted.**
The first pass reported that the 22 corpus rows whose `IL_name` contains `9322`
were "BL_2023 family, an unrelated index with no imidazole", and on that basis
recorded head 9322 as unobtainable — calling it "the most damaging fact in this
curation". **Both halves of that were false.** The rows are `YR_2024`, and 22 of
22 carry a 3-aminopropylimidazole head; their `IL_head_SMILES` is
`Cc1nccn1CCCN`, the 2-methyl analogue of amine 93.

The claim was never measured. It came from eyeballing a name in an exploratory
grep run before any substructure search, and my own leading-token regex required
the name to *start* with `93`, so a name beginning with `A` never reached any of
my three keys. That is a conclusion drawn from a secondary artifact instead of
the primary one — the pattern `REVIEW_PROGRAM.md` §5 exists for, landing on me.

What changed as a result: `01_blocks.csv` 6 → 7 blocks, `03_enumerated.csv`
5 → 10 lipidoids, a second C10 anchor, a 14th gating check, and the four
Figure S3 values for 9322-O17S gained a molecule. **No emitted row changed** —
`06_rows.csv` is still the same 13 Figure S2a formulations of 93-O17S, byte-for-
byte, because the emission rule never depended on head 9322. One *skip reason*
changed, from "no structure" to the dose-response reason, and that is recorded
rather than silently rewritten. Full detail in `08_retracted_suspicions.csv`
R-06; R-01's wrong cell is corrected in place and marked, not deleted.
