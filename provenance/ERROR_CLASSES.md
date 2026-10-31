# Error-class registry

The ratchet. Every bug any reviewer finds gets a row here and, wherever possible,
a permanent automated check. **This file only grows.** A class with a check can
never silently recur; a class without one is where review budget should go.

The point of keeping it: reviewers are expensive and find a class *once*; an
automated invariant is free and prevents it *forever*. Discovering a class is
worth paying for. Re-discovering one is not.

## Covered by an automated check

| # | Error class | Found | Check | Where |
|---|---|---|---|---|
| 1 | SMILES that does not parse | early | C1 | `audit.py` |
| 2 | Non-canonical SMILES (same molecule, two strings) | early | C1 | `audit.py` |
| 3 | Disconnected fragments / stray counter-ion | early | C1 | `audit.py` |
| 4 | **Radical electron from a deleted substituent** | 2026-10-29, Siegwart | C1 | `audit.py` + hard guard in both `quaternize()` implementations |
| 5 | Descriptor columns stale vs their own SMILES | 2026-10-28 | C8 | `audit.py` |
| 6 | Mass non-monotonic within a homologous series | 2026-10-28 | C2 | `audit.py` |
| 7 | Duplicate SMILES under two names in one experiment | 2026-10-28 | C3 | `audit.py` |
| 8 | Exact duplicate rows | 2026-10-28 | C4 | `audit.py` |
| 9 | Formulation ratios not summing to 100 | 2026-10-28 | C5 | `audit.py` |
| 10 | Non-finite or absurd `Experiment_value` | 2026-10-28 | C6 | `audit.py` |
| 11 | Field inconsistency; blank required fields | 2026-10-28 | C7 | `audit.py` |
| 12 | Lost E/Z on a stereogenic double bond | 2026-10-28 | C9 | `audit.py` |
| 13 | z-score convention drift (ddof, post-hoc row drops) | 2026-10-28 | rebuild | `build_dataset.py` |
| 14 | **Attachment at the wrong atom** (acyl on a side-chain methyl, not C3) | 2026-10-29, Siegwart | **C10** | `structure_check.py` |
| 15 | **Substituent on the wrong nitrogen** (regiochemistry) | 2026-10-28 Lin, 2026-10-29 Siegwart | **C10** | `structure_check.py` |
| 16 | Heatmap geometry off by a row/column (coloured header band) | 2026-10-28, twice | colour-distance gate | each `digitize_*.py`; >~5 means geometry, not colour |
| 17 | Wrong publication identifier (PMID ≠ DOI ≠ link) | 2026-10-29 | manifest build | `build_manifest.py` resolves every PMID/DOI against NCBI |
| 18 | **The PAPER draws a known benchmark wrong** (ALC-0315 head one carbon short) | 2026-10-29, LipidAI | **C10** `expected_mismatch` | `structure_check.py` — the override gets its own status so it stays visible instead of living in a prose note |
| 19 | **Chain length miscounted by eye from a figure** | standing risk, acted on 2026-10-29 | geometric vertex trace | `curation/xd_2026_lipidai/stage02_trace.py` — skeletonise, trace, RDP-simplify, count; never count by eye |
| 18 | **Two reagent names swapped** (MC3 ↔ KC2, identical formula and exact mass) | 2026-10-29, `AC_2025` | **C11** | `audit.py` — named-reagent registry |
| 19 | **Vendor shorthand expanded to the wrong lipid** (`18PG` read as symmetric 18:0/18:0; it is 18:0/18:1 SOPG) | 2026-10-29, `YZ_2022`/`YZ_2024` | **C11** | `audit.py` — named-reagent registry |
| 20 | **Log axis anchored on the panel box instead of its tick marks** (10³ tick sat 47 px above the bottom rule; every value high by up to 1.76×) | 2026-10-29, `ql_2026` Fig 3c | per-curation | reproduce every printed summary statistic, never a subset |
| 21 | **Provenance text asserts one thing, the code does another** (docs said the heatmap value was carried; the code took the conflicting bar-chart value) | 2026-10-29, `ql_2026` | per-curation | a check that inverts the stored value back to its source |
| 22 | **Blocks standardised separately when they should be pooled**, deleting a real between-block difference | 2026-10-29, `ql_2026` | per-curation | assert the between-block gap SURVIVES; per-block mean 0/sd 1 is a tautology |
| 23 | **Colormap floor/ceiling degeneracy left unflagged** (113 of 320 cells in a band where identical RGB spans ~980 RLU) | 2026-10-29, `ql_2026` Fig 2b | per-curation | flag floored and saturated cells symmetrically |
| 24 | **A C10 anchor documented as machine-derived but actually a typed literal** (the note and the README both said the B-a1 reference was read from vector bond coordinates; it was a string in the script, and nothing in any stage read the graph file) | 2026-10-30, `ng_2025` | per-curation | derive the anchor IN the script from the stored graph, so the provenance claim is executed rather than asserted; a claim no code runs is the same failure as copying the anchor from the enumeration |
| 25 | **Per-item drawings contradicting the paper's own reaction scheme** (about a third of the 104 dialkyl-ketone lipids are drawn as the open-chain mono-adduct, and every one of those also gives the amine nitrogen an extra methyl it cannot have) | 2026-10-30, `ng_2025` | per-curation | classify EVERY drawn member against the scheme and report the count; a spot-check gave 5 where the true figure was 15, understating the ambiguity six-fold |
| 26 | **A reconstructed rule validated on a subset but reported library-wide** (the >5x-MC3 hit threshold reproduces the 11 in-vivo lipids among the 64 phenols, but over all 310 it also catches a thiol lipid that was never taken into mice) | 2026-10-30, `ng_2025` | per-curation | run the reconstructed rule over the FULL set, and state the scope in the claim itself |
| 27 | **A validation calibrated on the quantity it validates** (`digitize_fig2fgh.py` fitted its pixels-per-percent scale by least squares onto the Fig 2e hit rates the panel was checking; the measured tick spacing was computed and never used, and the file had no `assert`) | 2026-10-30, `LL_2026` | **none possible in `audit.py`** — it is a property of a curation's validation code, not of the corpus | review rule: the pixel-to-unit scale must derive ONLY from axis furniture (ticks, labelled gridlines, printed anchors), never from the comparison target. Signature to grep for: a fit whose x-vector is the expected answer |
| 28 | **A measurement failure disguised as a measurement** (`else 0.0` fallback recorded a bar that exists as 0.0 because a `min_area` floor dropped its 450 px² component — and that fabricated zero became the set's single largest "disagreement") | 2026-10-30, `LL_2026` Fig 2h EP6 | per-curation | never let a detector miss fall through to a numeric literal: record `None` and fail, so an absent measurement cannot be averaged, deltaed or reported as agreement |

| 29 | **A hard-coded confidence literal discarding the per-item assessment** (`build_new_curation.py` read `lipid_name -> IL_SMILES` from the enumeration and then wrote `"high"` for every lipid it found, throwing away the `confidence` and `quaternization_site_assumed` columns beside it — so the 24 lipids whose record reads *"ASSUMED: no ring N and the paper draws no structure"* were labelled the project's highest confidence) | 2026-10-30, `ZT_2026` | **C13** | `audit.py` — plus: never write a confidence a per-item record already states; carry it. Inverting a disclosed assumption into an assurance is worse than publishing nothing, because the only surviving provenance then asserts the opposite of the truth |
| 30 | **A component looked up in the live corpus BY NAME, resolved by file order** (`comp.setdefault(name, smiles)` = first row wins, over a corpus written in alphabetical experiment order; 292 names map to more than one structure. `DLin-MC3-DMA` was *ambiguous* — 2,190 rows of the correct ester plus 203 dioxolane rows from `AC_2025` — so the lookup silently returned whichever experiment sorted first, and `MC3-CB1/2/3` were built on KC2) | 2026-10-30, `ZT_2026` | per-curation | the lookup must RAISE on an ambiguous name instead of picking one (`load_corpus_components` now returns the ambiguity set and `smi()` refuses). Regenerating after a `NAME_PATCHES` change fixes the symptom by luck of row order, not the cause; nothing tracks which artefacts were built from a corpus lookup |
| 31 | **Fragment assembly silently inverting a chiral tag** (`RWMol.RemoveAtom` deleted cholesterol's hydroxyl oxygen out of the middle of C3's neighbour list without correcting C3's tag, so all 537 `ZT_2026` rows stored epicholesterol, 3α. The epimers share molecular formula, exact mass to four decimals, the InChIKey SKELETON block and all 17 descriptor columns, so C1–C9 were blind; C10 was blind because its anchor was a byte copy of this code's own output; C11 was blind twice over, comparing skeleton blocks only and inspecting `CHL_*` columns that are empty on every one of these rows) | 2026-10-30, `ZT_2026` | **C14** | `audit.py` — cleave every ester inside `IL_SMILES` corpus-wide and pin the liberated sterol's FULL InChIKey against `reagent_registry.json`. Also: do the condensation the way the chemistry does (the ACID loses its OH, the alcohol's oxygen is retained), so no atom adjacent to a stereocentre is ever deleted |
| 32 | **A bar-top detector fooled by replicate markers** (replicate dots straddling the top edge occlude the MIDDLE of the bar width, so a "first row where >60% of the width is fill" rule fires several pixels INSIDE the bar; the bias is one-sided and non-constant, so it does not cancel in a z-score) | 2026-10-30, independently in `QL_2026` Fig 2c/3c and `XD_2026` Fig 5C **on the same day** | per-curation | measure from the FLANK columns, which a centred marker cannot reach, and require the left and right flank to agree to ≤1–2 px. Note the limit: every printed hit rate in the `QL_2026` panels still reproduces on the corrected values, so the paper's own statistics cannot discriminate this one — it has to be caught geometrically |
| 33 | **A check whose findings are printed but never gate** (`audit.py` discarded `main()`'s return value; C1–C9 printed and the exit code came only from the C10/C11 delegates, so eight simultaneous seeded violations exited 0. Of 43 mutations, 30 appeared in the report and 5 failed a command. `build_dataset.py` never imported `sys` at all) | 2026-10-30, corpus-wide | **the gate itself** | `audit_baseline.json` + `test_checks.py` — each metric has an accepted value WITH its reason, departure in EITHER direction fails (a vanished known condition is how a check that stopped looking announces itself), and a seeded-mutation suite proves the gate fires. The obvious fix is wrong on its own: failing on anything printed fails forever on conditions we accept on purpose, which is how it became print-only |
| 34 | **A quality flag produced, documented, and then dropped at the schema boundary** (1,137 rows carried a convention rather than a measurement; `CLAUDE.md` claimed the project flagged them "in the data (`positional_ambiguity`, `*_assumed`, `confidence`)" and none of those columns existed anywhere. The loss was one stage EARLIER than the merge — no file `build_dataset.py` reads had a flag column — so no change to the builder could have preserved them) | 2026-10-30, four of five curations | **C13** | `schema_quality.json` as the single source of the column names and their vocabulary; `audit.py` enforces the vocabulary and counts rows flagged without a note. A flag that lives only in a per-folder QC JSON is not a flag, it is a memo |

| 35 | **A mask change silently invalidating every window tuned against the old mask** (switching Fig 3c's bar detector from a blue-hue mask to a luminance mask made the panel box's own top rule visible at y = 2301-2314, and the search window started at 2311 — inside it — so every bar would have returned the same wrong top) | 2026-10-30, `QL_2026` Fig 3c | per-curation | assert the first scanned row is not already a full-width rule. Generally: a mask defines what counts as ink, so every threshold, window and offset calibrated under the old one has to be re-derived, not inherited |
| 36 | **A calibration measured mostly from our own detector error** (the `QL_2026` cross-library offset was computed from four benchmark LNPs plotted in both panels — and three of the four were dot-contaminated, so the 1.1116 "offset between the two screens" was largely the bar-top bug measuring itself. Correcting the detector did not sharpen the offset, it removed it: 1.0071 ± 0.0247, not distinguishable from 1.0) | 2026-10-30, `QL_2026` | per-curation | a quantity derived from a small set of shared anchors inherits every bias those anchors carry; check the anchors for the defect under investigation BEFORE trusting a constant derived from them |

## A note on the numbering above

Classes **18** and **19** each appear twice: the figure-drawing/vertex-trace pair
and the reagent-name/vendor-shorthand pair were added on the same day by
different passes. **Still duplicated, deliberately.** Round 1's reports are
written against these numbers, so renumbering now would make every reference in
`review/round_2026-10-30/` point at the wrong row. Renumber when round 2 opens,
in one pass, and retire the old numbers rather than reusing them.

**This file is edited concurrently by more than one session**, which is how the
duplication happened and how classes 27-28 below were first written as 24-25,
colliding with 24-26 added from another session minutes earlier. Before
appending, re-read the table and take the next free number; do not assume the
highest number you last saw is still the highest.

## Known, not yet automatable

These need eyes on a figure or judgement about a paper's text. This is where a
reviewer's time is actually worth spending.

| Error class | Why no check | Mitigation |
|---|---|---|
| Name token misread (`C24` = two C4 tails, not 24 carbons) | needs the library's own naming convention | C2 restricted to conventions verified against source metadata |
| Value misread from a figure | needs the figure | overlay PNG per extraction; spot-check ≥10 cells |
| Experimental metadata contradicting itself between caption and Methods | needs both read in full | record BOTH and flag UNRESOLVED (e.g. Siegwart 20 ng/n=3 vs 15 ng/n=5) |
| Commercial reagent that is an isomer mixture | no single true structure | flag `confidence = medium` and say so |
| Regiochemistry a paper never pins | genuinely unknowable | flag `positional_ambiguity` / `*_assumed`, confidence `low` |
| **A C10 anchor copied from the enumeration** (can only pass) | once both sides are canonicalised a copied reference is byte-identical to an honestly re-derived one, so a string comparison would fail the curations that did it right | assemble the reference by a different code path — fragment assembly, or just a different atom ordering — and say which in `read_by` / `note` |
| **De-duplicating reference strings on an assumed identifier format** | LipidAI appended a per-row counter to the DOI *sometimes*, and to a PMCID or a page range otherwise; cutting only at `doi:` over-counts badly. The class is real; the NUMBERS first recorded here were not. "299 against a true 47" was a reviewer claim propagated unchecked; the folder's own `06_lipiddb_qc.json` records 56/50 and the paper states 47, and nine defensible de-duplication rules give 50-116. | cut at the first identifier of ANY kind, cross-check against the count the paper states itself, and **report the rule with the number** - a de-duplicated count is meaningless without it |
| **A structure with no external anchor at all** | commercial lipids supplied without name, formula, HRMS or synthesis (LipidAI's E1–E9) | emit with `confidence = medium`, record an empty `reference_SMILES` row so the absence is visible, and never upgrade without an external structure |

## What the composition checks cannot see

Worth stating plainly, because it is the reason three structural bugs survived
every check that existed at the time:

- **Mass and molecular formula are blind to connectivity.** HRMS matched 4/4 to
  ±0.0001 Da on a molecule with the acyl group on the wrong end of the steroid.
- Deleting a hydroxyl left the **oxygen count unchanged**; only hydrogen moved,
  by one.
- Moving an arm between nitrogens changes **nothing** about the formula.

C10 exists specifically to close this gap, and it is the only Tier-1 check that
compares against an external ground truth rather than internal consistency.

**2026-10-29 added a second instance of exactly this.** `AC_2025` held the
DLin-KC2-DMA dioxolane under the name `DLin-MC3-DMA` and the DLin-MC3-DMA ester
under the name `KC2`. Both are C43H79NO2 at **641.6111 Da** — the same formula
and the same exact mass to four decimals. Every composition check passed over
419 rows for as long as they were in the corpus.

It also taught a second lesson, the harder one. Finding that a name and a
structure disagree does **not** tell you which of the two is wrong, and the two
repairs are not equivalent: rewriting the structure keeps the value with the
name, rewriting the name keeps the value with the structure, and they assign
every measurement to opposite molecules. The first fix here rewrote the
structures and was wrong. Only the upstream source settled it — COMET's master
table stores SMILES and labels with **no names at all**, so the values belong to
the structures and the names were a downstream annotation (built, per its own
comment, by matching molecular weight, which cannot tell these two apart).

**So: when a name and a structure disagree, go to the source and find out which
one the measurement was recorded against. Do not repair the cheaper half.** An
internal-consistency check cannot answer this — it is satisfied by either
repair, including the wrong one.

## What a figure's own printed summaries are for

The `ql_2026` Fig 3c axis error (class 20) is the cleanest case yet of a bug
that **no internal check could see** and that the paper itself gave away three
times over. The values were self-consistent, the bars were detected correctly,
the z-scores were perfect — and every number was wrong by up to 1.76×.

What caught it: the panel's own printed hit rates, and a sentence in the
Results. What let it survive a first pass: `validate.py` asserted two of the
three printed rates and **excused the third** as "a boundary case" — and the
excused one was the only one that discriminated.

The rule this earns: **when a paper prints a statistic derived from a figure you
are digitizing, reproduce ALL of it, and never excuse the part that disagrees.**
A disagreeing summary statistic is the single most informative thing available,
because it is computed from the data you are trying to recover and is
independent of how you recover it. Excusing one is how a calibration error gets
through a passing test suite.

C10 could not have caught it either: C10 only looks at lipids a *curation folder*
enumerated, and these are catalogue reagents that arrive named. That gap is what
C11 closes — a registry of named commercial reagents whose structure is fixed and
public, checked by InChIKey skeleton wherever that name appears in the corpus.
A name like `DSPC` or `SM-102` means one molecule; the corpus should be able to
prove it means that molecule everywhere, and a shorthand like `18PG` should be
expanded once, in one place, with the paper that defines it cited next to it.
