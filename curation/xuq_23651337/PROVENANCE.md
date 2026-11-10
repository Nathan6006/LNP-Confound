# PROVENANCE — `xuq_23651337`

Wang M, Sun S, Alberti KA, Xu Q. **"A combinatorial library of unsaturated
lipidoids for efficient intracellular gene delivery."**
ACS Synth Biol 2012;1(9). PMID 23651337 · DOI 10.1021/sb300023h.
Curated 2026-11-10.

---

## 1. The endpoint, and where its values are

**The paper's endpoint is EGFP expression in HeLa cells after transfection
with EGFP-encoding plasmid DNA and with EGFP mRNA**, extended to fibroblasts
and to several cancer cell lines. That is a reporter-protein-abundance
delivery endpoint and it is **admissible in kind**.

**Every one of those values is in the main article, which is not in hand.**

Title and endpoint verified from the Europe PMC REST `core` record
(`EXT_ID:23651337 AND SRC:MED`, fetched 2026-11-10). `abstractText` is **797
characters** and describes the paper's own synthesis, screen and mechanistic
attribution, so this is ordinary research despite `pubType = Letter` — the
no-original-data discriminator is a zero-length abstract plus no
Methods/Results, and this paper fails that test in the right direction.

## 2. Acquisition — the hard part, and it was not solved

Only `source/sb300023h_si_001.pdf` (3 pages, sha256
`92256854c364c0b0e32e3b4ff956177e7dbcb9c0defc2572b49f627fbf1ae353`) was
available. It was already in the folder; this curation fetched nothing and used
no browser. Its own watermark shows it came from an institutional subscription
session.

The main article was **not obtained**, with these checkable negatives:

| probe | result |
|---|---|
| Europe PMC `core` | `isOpenAccess=N`, `inEPMC=N`, `inPMC=N`, `hasSuppl=N`; one `fullTextUrl`, `availability='Subscription required'` |
| `elink dbfrom=pubmed db=pmc id=23651337` | **only** `pubmed_pmc_refs`, 17 links (papers citing it). **No `pubmed_pmc` self-link → no PMC deposit** |
| `esearch db=pmc term=10.1021/sb300023h` | **count=0**, "No items found" |
| `esearch db=pmc term=sb300023h` | **count=0**, `quotedphrasesnotfound=['sb300023h']` |
| `find_repo.py --pmid 23651337 --title "unsaturated lipidoid library" --author Xu` | no candidate repository; all hits unrelated string collisions |
| `find ~/LNP -iname '*sb300023*'` | only the SI already in hand; `~/LNP/Papers/` holds no Wang 2012/2013 |

ACS needs a browser and no paywall or bot-check bypass is permitted, so this
was not pushed further. See `DECISIONS.md` **D-1**.

## 3. Tier, and why

`00_value_sources.json` is the intake gate and was written **before** any
digitiser.

* Tiers **1, 2 and 3 are all unavailable**: ACS Synth Biol 2012 published no
  Source Data (`hasSuppl=N`, and the only supplementary asset is this PDF);
  `pdfplumber.find_tables()` returns **0 tables on all 3 pages** and the SI has
  no numeric table, no synthesis section and no experimental section; and there
  is no author repository.
* Tier **4b — `4b_pdf_vector_geometry`** is what was used, for the one panel
  that carries per-item numbers. **Nothing was rendered and no pixel was
  read.**

### The vector case, stated as asked

`page.get_drawings()` on SI p. S3 returns **239 path groups**:

| objects | stroke width | what they are |
|---|---|---|
| 16 filled rectangles (black) | — (fill) | the bars |
| 16 stroked rectangles | **0.4596 pt** | their outlines, geometrically coincident with the fills |
| 48 line segments | **0.4596 pt** | bar edge strokes |
| 64 line segments | **0.2757 pt** | the error bars — 2 stems + 2 caps per bar |
| 94 line segments | **1.1949 pt** | axis rules, major and minor ticks |

**Bar body, outline, error bar and axis furniture separate by EXACT STROKE
WIDTH.** They are selected by that width, not thresholded apart, so **error
class 32 cannot arise here** — F-018 and F-020 were raster bar-top detectors
reading an overlaid marker, and that failure mode has no analogue in a vector
read. The error-bar stems' shared endpoint is independently available as a
*third* attestation of each bar top.

### Native dpi and tolerances

**Not applicable — the panel is vector.** There is no embedded raster to
measure, `axis.native_tol()` was deliberately **not** called (it refuses a
guessed dpi, and a guess is all we could give it), and **no tolerance was
widened anywhere in this folder.**

Both tolerance figures, as asked:

* `curation/lib`'s own, applied where it belongs: **`axis.LOO_TOL_PX = 1.5`
  px**, and `axis.calibrate`'s hard **2.5 px** residual limit. Every coordinate
  is converted from PDF points to 600-dpi pixels (×600/72) before being handed
  to `axis`, which is exact and lossless on vector coordinates and is what
  makes those tolerances mean what they were written to mean. Measured:
  calibrate residual **0.4427 px**, worst LOO **0.5405 px** (= 0.0649 pt).
* This folder's own, for drawn objects meant to be coincident: **0.01 pt**
  (= 0.083 px at 600 dpi), stated in `stage05_values.py`. A vector read has no
  detector noise — its only slack is PDF coordinate rounding — so 1.5 px would
  be ~75× too loose to catch anything, and using it would have been a silent
  widening in the other direction.

## 4. Library size: distinct / formulations / individually measured

| quantity | value |
|---|---|
| distinct structures **determinable** | **0** |
| distinct structures that **exist** | **NOT DETERMINABLE** |
| formulations | **NOT RECOVERABLE** |
| individually measured, in the file we have | **16** |
| individually measured **on a delivery endpoint** | **0** |

* **0 determinable.** A lipidoid here is (amine head + *n* × oleyl
  acrylamide), and the SI names, draws and formulates **no amine**.
* **Not determinable how many exist.** 16 unsaturated lipidoids are
  individually measured in Figure S3, and at least one saturated counterpart
  exists (lipidoid 16's other form, Figure S2). Whether the saturated series
  has 16 members or only the one shown is stated nowhere in the SI.
* **Formulations not recoverable at all.** No lipidoid:DNA weight ratio, no
  N:P, no helper lipid, no dose, no particle concentration, no timepoint. The
  SI has no experimental section.

The two **tail reagents** are the only structures this paper yields:

| | SMILES | formula | exact mass | confidence |
|---|---|---|---|---|
| oleyl acrylamide | `C=CC(=O)NCCCCCCCC/C=C\CCCCCCCC` | C21H39NO | 321.3032 | high |
| n-octadecyl acrylamide | `C=CC(=O)NCCCCCCCCCCCCCCCCCC` | C21H41NO | 323.3188 | medium (linkage assumed — D-3) |

## 5. Corpus overlap, and the free external check

* **Lipidoid-level overlap: none possible** — nothing enumerated, so the key
  that matters (a lipidoid's own full InChIKey) has nothing to compare.
* **Tail reagents: 0 full-InChIKey matches, 0 skeleton-block matches.** Zero is
  the expected answer; a reagent is a monomer, and a hit would have meant the
  corpus stored a reagent as a library member.
* **Name tokens: 0 exact matches.** **5,075** distinct corpus `IL_name`s merely
  *start* with one of this paper's tokens — reported to size the trap, never as
  candidates. A leading number token is not an identity.
* **Near misses, as a superset then strict.** The loose SMARTS for a Michael
  adduct of oleyl acrylamide returned **19** unique corpus SMILES; walking each
  chain atom by atom (unbranched, ring-free, all-carbon, CH₃-terminated) leaves
  **0 of 19**. All 19 were a C18:1 chain read from the methyl end into a
  substituted interior carbon. **The corpus contains no Michael adduct of oleyl
  acrylamide.** Recorded as `08_retracted_suspicions.csv` **R-2**.

### The free external check this paper did give

`AA_2008` (Akinc/Anderson 2008), curated into this corpus independently and
from a different paper, names its tails `<amine>_<N|O><n>` — `N` = acrylamide,
`O` = acrylate, `n` = alkyl carbons. So **`N18` is N-octadecylacrylamide, this
paper's own saturated comparator.** Two results:

1. **43 of 43** `AA_2008` `N18` lipids, retro-Michaelised, liberate **exactly**
   the structure this folder built from the chemical name alone (full
   canonical-SMILES identity; 0 armless, 0 liberating anything else). This
   attests *the structure for the name*. It does **not** attest that this paper
   used the amide — see D-3.
2. **The catalogue-number hypothesis is measurably rejected.** The
   Anderson/Akinc amine numbers actually in use are **47 values and they are
   not contiguous**: 1, 6, 7, 10, 11, 13, 15, 17, 20, 21, 22, 24, 25, 26, 28,
   31, 32, 33, 34, 36, 38, 60–64, 70, 75–77, 79–82, 86, 87, 90, 91, 93–96,
   98–100, 103, 109. **Nine of this paper's sixteen tokens — 2, 3, 4, 5, 8, 9,
   12, 14, 16 — are not catalogue numbers at all**, while Figure S3 plots a
   contiguous run 1..16. Under the hypothesis all 16 would have to be members.
   **Rejected.**

   This is the most valuable thing the sweep produced. `AA_2008` *does* store a
   structure for 7 of our 16 tokens, so reading the labels as catalogue numbers
   was available and would have put up to 16 measured values on wrong
   molecules. It is closed by measurement, not by preference, and it is seeded
   (`seed_tokens_are_catalogue`).

## 6. Anchor posture — recorded BY HAND

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run**, by instruction: they write tracked
files other agents read. This section is the hand-written substitute.

**The paper draws nothing.** The 3-page SI contains two TEM micrographs, a line
plot and a bar chart, and no chemical structure of any kind.
`reference_structures.csv` therefore holds **three rows**:

| `lipid_name` | `reference_SMILES` | `expected_mismatch` |
|---|---|---|
| oleyl acrylamide | `C=CC(=O)NCCCCCCCC/C=C\CCCCCCCC` | *(empty — matches)* |
| n-octadecyl acrylamide | `C=CC(=O)NCCCCCCCCCCCCCCCCCC` | *(empty — matches)* |
| lipidoid 1-16 (the library itself) | **empty, on purpose** | `no_anchor_possible_paper_draws_nothing` |

Header is exactly
`lipid_name,reference_SMILES,source_figure,read_on,read_by,expected_mismatch,note`.

**C10-equivalent result: 2/2 anchors InChIKey-identical to the enumeration.**
The third is an empty anchor, so the *absence* is visible rather than silent.

### Independence, established by provenance and not by byte-identity

Byte-identity of the emitted anchor is **not** the signal and testing for it is
a trap (F-037): a correct re-derivation is byte-identical to a correct
enumeration once canonicalised, which is exactly why C10's own
`anchor_may_be_a_copy` missed the defect it was written for. So:

* **`make_reference.py::assert_no_assembler()` walks its own module's AST at
  run time** (1,500 nodes) and confirms no import of `assemble`/`axis`/`bars`/
  `heatmap`/`chains`/`reagents`/`emit`/`lib` and no call to `RWMol`,
  `AddBond`, `RemoveAtom`, `CombineMols`, `esterify`, `attach`, `ReplaceCore`,
  `ReplaceSidechains`, `emit_rows` or `calibrate`. **Executed, not asserted in
  a comment** — the F-013 note claimed independent assembly "through the C3
  oxygen (3-beta)" about a string that was 3α and a byte copy.
* **The AST guard was made to fail on purpose.** `stage07_qc.py` copies
  `make_reference.py` into a temp directory, inserts `import assemble` *and* a
  `Chem.RWMol()` call, runs it, and requires a non-zero exit. It fires.
  (`make_reference.py` itself is untouched; the copy is deleted.) My own
  instinct was to trust the clean line — C16's first two versions reported "all
  derived" while examining nothing.
* **The traversal is opposite.** `stage03_enumerate.py` starts at the acryloyl
  nitrogen and `assemble.attach`-es a fragment whose first atom is C1.
  `make_reference.py` writes the chain **methyl-end-first** (C18 → C1) by pure
  string construction and appends the amide last, computing every index from
  the opposite end.
* **The gate is proved to test connectivity, not mass.**
  `make_reference.py --seed reverse` shifts the alkene one carbon further from
  nitrogen. Formula **C21H39NO** and exact mass **321.3032** are *identical*
  either way and both SMILES parse and canonicalise; only the InChIKey
  separates them (`QWACECOCWJDWKP-KHPPLWFESA-N` vs
  `JROSKIFCGJOXDY-QXMHVHEDSA-N`). The gate fails and the agreement count moves
  2 → 1.
* **Chain-length coverage, honestly.** "Pick anchors of different chain
  lengths" **cannot be satisfied from this paper** — both tails are C18 and
  that is the whole point of the comparison it makes. Instead the *generator*
  is exercised at n = 12, 14, 16, 18 (formula checked against a hand-written
  C(n+3)H(2n+5)NO / C(n+3)H(2n+3)NO expectation) and stage03 carries
  `seed_chain_length`. That covers the off-by-one-that-scales-with-n hazard;
  no shorter product is claimed to be in the paper.

### Stereo posture — VACUOUS, not passing

**Neither product has a single stereogenic carbon.**
`assemble.check_stereo_preserved` compares CIP code *counts*, so on these
molecules it compares 0 to 0: **it is vacuous, and a clean line from it is
evidence of nothing.** `assemble.parity_preserved` was not called — it compares
spatial parity of assigned stereocentres, of which there are none.

The stereochemistry that *does* exist is the oleyl **C9=C10 Z** double bond,
which CIP-centre counting cannot see at all. It is therefore asserted
explicitly by `stage03_enumerate.py::double_bond_stereo()` and seeded
(`seed_alkene_stereo`, which also leaves formula and exact mass untouched).

## 7. Every check, with numbers

Eight panel checks on Figure S3, all passing (`05_values_qc.json`):

| check | result | measured |
|---|---|---|
| CA y label → tick, absolute nearest-tick | **6/6** | worst margin **29.81 pt** vs 5.0 pt required |
| CB y calibration from ticks only | **6/6** | residual **0.4427 px**, worst LOO **0.5405 px** vs 1.5 px |
| CC minor-tick **hold-out** (unprinted) | **6/6** | worst **0.0964 pt** = 0.803 px |
| CD bar top, **three independent objects** | **16/16** | worst spread **0.00269 pt** |
| CE baseline on the 0 % tick | **16/16** | worst **0.00161 pt** |
| CF bar → label assignment, absolute | **16/16** | worst margin **12.3591 pt** vs 5.0 pt |
| CG frame top = **110 %** (unprinted) | **1/1** | **0.00191 pt** = 0.016 px |
| CH one x grid, 4 independent series | **5/5** | see below |

**Numbers the authors drew and never printed** (external to our read in a way
a recomputed marginal is not):

* **The plot box is exactly 16 bins wide.** Measured width **199.6795 pt**;
  16 × the measured 12.48 pt pitch is **199.6799 pt**. Agreement **0.00048 pt**.
* **The frame's top rule lands on the unlabelled 110 % level**, to 0.00191 pt.
  110 % is not a value any measurement produces and nothing in our read could
  have been fitted to it.
* **Six minor y ticks, all unlabelled**, predicted at 110/90/70/50/30/10 % by a
  fit that *excluded them* — they were passed to `axis.calibrate` as
  `forbid_target`, making error class 27 a type error. Worst error 0.0964 pt.
* **One x grid under four independently drawn series.** Bar right edges
  12.479997 pt (residual **0.00111 pt** — an exact grid); bar left edges
  12.481350; minor x ticks 12.480023; major x ticks 12.475629.

Structure and overlap checks:

| check | result |
|---|---|
| tail formula / chain length / alkene position / alkene geometry | **3/3** both tails |
| anchor vs enumeration, full InChIKey | **2/2** |
| generator formula scaling at n = 12, 14, 16, 18 | **4/4** |
| `AA_2008` N18 retro-Michael → our octadecyl acrylamide | **43/43** |
| catalogue-number disproof | **9 of 16** tokens non-catalogue → rejected |
| `emit` guards (literal `confidence`, caveat w/o note, correct call) | **3/3** |

### Seeding — 16 of 16 seeds credited

**A check is credited only when its AGREEMENT COUNT MOVES.** Each seed runs in
**its own process** and the verdict is read from `subprocess.run(...).returncode`.

**The harness caught one of my own no-op seeds.** `seed_major_tick` at +3 pt
exited non-zero but moved nothing, because `axis.calibrate` *refuses* (residual
20.92 px vs its 2.5 px limit) before any count can be computed. It was split:
`seed_major_tick` is now a **subtle +0.25 pt** displacement that passes
calibrate's leverage-masked residual test and fails leave-one-out (CB 6→5),
while `seed_major_tick_gross` keeps +3 pt and is credited by the **named
library refusal**, recorded as a distinct route (`credited_by`) so the two are
never blurred.

Two other things worth recording because they nearly went wrong:

* **`python3 … | tail` gives `tail`'s exit code.** Two stage02 seeds printed
  `FAILED` next to `EXIT=0` for exactly that reason during development. That is
  F-012's shape, re-encountered. `stage07_qc.py` never infers a verdict from
  stdout or from a pipeline's status.
* **A refusal was respected, not loosened.** CF originally required the bar
  centre to coincide with its x tick within 0.02 pt and **refused on 8 of 16
  bars at 0.0607 pt**. The threshold was *not* widened. The source was found:
  the panel's x geometry is device-quantised (bar right edges exact to
  0.00111 pt, everything else alternating ±0.06 pt about the same grid, the
  alternation being half the 0.1191 pt difference between the file's two bar
  widths). An exact-coincidence criterion was testing a claim the file does not
  make, so for a *categorical* axis it was replaced by an absolute nearest-tick
  criterion and CH now **measures** the quantisation. Recorded as
  `08_retracted_suspicions.csv` **R-1**.

## 8. What is emitted

**Zero rows.** `06_rows.csv` is written with its header and no data rows — an
empty-but-headed file, the same posture as the empty `reference_SMILES`, so the
zero is visible rather than indistinguishable from a stage never run.

`05_values.csv` holds **18 rows**: 16 Figure S3 viability values plus one row
each for Figure S2 and Figure S1, so the accounting is complete rather than
silently short. **`corpus_admissible = no` on 18 of 18**, each with its own
reason — viability is cytotoxicity, EB exclusion is nucleic-acid binding in a
cuvette, and micrographs carry no number. `data/`'s `Experiment_method`
vocabulary is delivery-only. Counts sum: 16 + 2 = 18, and 0 + 18 = 18.

The `emit` guards are still exercised, against **synthetic** items, because
"we would have flagged correctly" is not a claim this folder gets for free:
`emit.emit_rows` refuses a hard-coded literal `confidence` (F-011/F-025/F-046,
error class 29) and refuses a caveat with an empty note (F-024), and a control
call succeeds so the suite can distinguish "refuses the right things" from
"refuses everything".

## 9. What could NOT be verified

* **Any delivery value.** The main article is unobtainable. This is the whole
  of the paper's corpus relevance.
* **All 16 amine heads**, therefore all 16 lipidoids. Not named, drawn or
  formulated anywhere in the SI.
* **Which of `16I`/`16II` is unsaturated** (D-2).
* **Whether the saturated comparator is the amide or the ester** (D-3). The
  corpus attests the *structure for the name*, not this paper's *choice*.
* **Every formulation parameter**: ratio, N:P, dose, timepoint, helper lipid,
  particle concentration. The SI has no experimental section.
* **The number of distinct structures in the library.**
* **The licence** of the one file we have (D-4), and `isOpenAccess = N`.
* **The oleyl reagent against anything external.** The corpus contains zero
  Michael adducts of oleyl acrylamide, so the unsaturated tail — the paper's
  entire point — has no external attestation. Its structure rests on the
  chemical name alone.

## 10. Files

```
00_value_sources.json        the intake gate, written FIRST
00_acquisition.json          hashes, licence gap, checkable negatives
01_blocks.csv / _qc.json     2 tail reagents; the 16-amine gap, explicit
02_corpus_overlap.json       InChIKey + skeleton + name + superset sweep + (d)
03_enumerated.csv / _qc.json 0 lipidoids, 2 tail reagents, stereo posture
05_values.csv / _qc.json     18 located values, 0 admissible, 8 panel checks
06_rows.csv                  header + 0 rows
06_emit_reconciliation.json  why zero, and the emit guards
07_qc.json                   every check and every seed, counts and exit codes
08_contradictions.csv        4 internal contradictions
08_retracted_suspicions.csv  4, each with the measurement that retired it
DECISIONS.md                 D-1 .. D-5, all OPEN
make_reference.py            the anchor generator + AST independence proof
reference_structures.csv     2 anchors + 1 deliberately empty
stage01/02/03/05/06/07*.py   the pipeline; all import curation/lib, none copies it
```

Reproduce: `python3 -I stage01_blocks.py && python3 -I stage03_enumerate.py &&
python3 -I stage02_corpus_overlap.py && python3 -I stage05_values.py &&
python3 -I stage06_emit.py && python3 -I make_reference.py &&
python3 -I stage07_qc.py`
