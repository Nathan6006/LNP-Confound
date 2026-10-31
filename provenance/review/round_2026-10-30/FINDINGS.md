# Round 1 findings ledger — 2026-10-30

Corpus at round start: `data_vitro.csv` 28,885 × 72 (42 experiments),
`data_vivo.csv` 1,780 × 72 (19 experiments). Five merged curations:
`ZT_2026`, `LL_2026` + `LL_2026_2`, `WT_2026`, `XD_2026`, `QL_2026`.

Procedure: `../REVIEW_PROGRAM.md`. Reports land as `report_<TARGET>.md` in this
directory.

**Every row below is a CLAIM until the adjudication column says otherwise.** A
reviewer's finding is evidence to check, never an instruction to edit. Of the
twelve findings from the first Siegwart reviewer, three were real structural
errors, one was half-right, and one was flatly wrong.

## Status key

| | |
|---|---|
| `claimed` | a reviewer reported it; I have not reproduced it |
| `confirmed` | I reproduced the evidence myself, by my own route |
| `rejected` | I tried to reproduce it and it does not hold — with the reason |
| `validated` | a second agent confirmed the proposed repair is the right one |
| `fixed` | applied, with a `CHANGELOG.md` entry |
| `open` | real, but the repair is the user's decision |

---

## F-001 — quality flags never reach the corpus

| | |
|---|---|
| **found by** | me, 2026-10-30, while generating the review packets |
| **status** | `confirmed` (the fact), `open` (the repair) |
| **severity** | changes no value or structure, but **misrepresents the data to every downstream consumer** |
| **confidence** | certain on the fact |

**The claim in the documents.** `CLAUDE.md` §6 ("Residual risk, stated
honestly") states the project's posture is to "flag those rows in the data
(`positional_ambiguity`, `*_assumed`, `confidence`), so a reader can see which
values are measured and which are convention." `CHANGELOG.md` §11 states that
`XD_2026`'s E1–E9 "carry `confidence = medium`".

**What the data holds.** Neither is true. The 72-column schema contains no
`confidence`, no `positional_ambiguity`, no `*_assumed`, and no flag column of
any kind:

```
head -1 data/data_vitro.csv | tr ',' '\n' | grep -i "conf\|ambig\|flag\|note\|censor"
  -> no output
```

Nor does any curation output file. Checked `11_QL_2026_rows.csv`,
`07_XD_2026_rows.csv`, `05_WT_2026_merge.csv`, `new_curation.csv` and
`10_LL_2026_rows.csv`: none has a matching column.

Every flag any curation produced exists **only** inside that folder's QC JSON:

| folder | flags held outside the corpus |
|---|---|
| `ql_2026` | `12_QL_2026_qc.json` → `flags`: 3 censored, 113 floored, 1 conflict (`T3A-L8-C14`) |
| `ql_2026` | `08_fig2c_qc.json` → `still_censored` |
| `xd_2026_lipidai` | `07_qc.json` → `E_lipids_confidence` (E1–E9, medium) |
| `lin_aromatic_2026` | `09_enumeration_qc.json` → `medium_confidence_rows`, `medium_confidence_reason` |

**Why it matters here specifically.** 113 of `QL_2026`'s 344 rows — a third of
the block — sit on a colormap floor where identical RGB spans ~980 RLU, and they
were z-scored as though they were measurements. A model consuming
`data/data_vitro.csv` cannot tell them from the clean rows, and neither can a
reader of the manuscript's data-availability release. The same applies to
`XD_2026`'s nine vendor-supplied structures with no external anchor of any kind.

**Scope, measured.** Four of the five curated papers have rows that should
carry a flag and do not:

| experiment | rows affected | of | flag that is missing | recoverable? |
|---|---|---|---|---|
| `QL_2026` | 113 | 344 | value sits on the Fig 2b colormap **floor** (identical RGB spans ~980 RLU) | yes — `01_fig2b_320.csv` has a per-lipid `floored` column |
| `QL_2026` | 3 | 344 | value **censored** above the Fig 2c axis | yes — named in `12_QL_2026_qc.json` |
| `QL_2026` | 1 | 344 | `T3A-L8-C14`, two figures of the paper **disagree 7.3×** | yes — named |
| `LL_2026` | 140 | 840 | `Ac8-1`/`Ac10-1` are iso-octyl / iso-decyl acrylate, commercial **isomer mixtures**; the single stored SMILES is a convention, not a measurement | count only; the row set must be re-derived |
| `XD_2026` | 9 | 12 | E1–E9, vendor-supplied, **no external anchor of any kind** — no name, formula, HRMS or SI structure | yes — E1–E9 by name |
| `ZT_2026` | unquantified | 486 + 51 | `O1` (branched isodecyl) figure-only; `B1`/`B2`/`B3` ester connectivity inferred from a C3-H shift — both recorded as **medium** | **no** — held as a prose sentence in `07_qc.json`, not machine-readable and not row-keyed |

At least **266 rows** across four experiments, and the `ZT_2026` and `LL_2026`
sets would have to be reconstructed because they were never written down
row-wise. `WT_2026` appears genuinely clean here — it is author-reported tabular
data with no digitization and no inferred structures.

All four `QL_2026` flag sets match corpus `IL_name` values 1:1 (113/113, 3/3,
1/1), so a sidecar keyed on `LNP_ID` is mechanically straightforward for that
paper; the other two need work first.

**Not yet decided, because it is the user's call.** Adding columns changes the
schema from 72 to 72+k, which touches the data dictionary, `build_dataset.py`,
the manifest, and anything that assumes the column count. The alternatives are
(a) add the flag columns, (b) ship a sidecar flags table keyed by `LNP_ID`, or
(c) correct the two documents to describe what the data actually is. (c) alone
would leave the honest-flagging posture unimplemented.

**Independent corroboration requested.** Each Wave-1 agent was asked to
determine, for its own paper, which flags reached `data/` — without being told
this finding — so the per-paper evidence arrives independently. `CORPUS-B` was
asked to adjudicate the documentation claim directly.

---

## F-002 — `build_manifest.py` asserted a frozen row count

| | |
|---|---|
| **found by** | me, 2026-10-30, while rebuilding the manifest after the merge |
| **status** | `fixed` |
| **severity** | cosmetic in effect, but it blocked the build and reported a correct manifest as broken |
| **confidence** | certain |

`provenance/build_manifest.py` asserted `sum(n_rows) == 25997` for vitro and
`== 1768` for vivo — literals frozen at an earlier corpus size. The assertion
fired as `*** DISCREPANCY ***` and raised `AssertionError` the moment `QL_2026`
merged, against a manifest that was entirely correct.

The check's purpose is "the manifest's per-experiment counts add up to the
corpus", so the reference must be the live file length, not a constant.
Rewritten to assert against `len(frames[corpus])`, with a separate `BASELINE`
dict kept as a non-fatal tripwire that prints when the corpus changes size, so
an unintended change is still visible but an intended one does not break the
build.

Error class: a check whose reference value goes stale on every legitimate
change. It trains whoever hits it to edit the number until the assert passes,
which is how a real discrepancy would eventually get waved through.

---

## Reviewer findings

Filled in as reports land.

| id | target | finding | severity | reviewer confidence | adjudication |
|---|---|---|---|---|---|
| | | *Wave 1 and Wave 2 running* | | | |

---

## F-003 — `WT_2026`: per-cell-type z-scoring deletes and then inverts the authors' selectivity signal

| | |
|---|---|
| **found by** | WT_2026 reviewer (its E1) |
| **status** | `confirmed` — I reproduced it from the raw source by my own route |
| **severity** | **changes values** on the multi-cell-type blocks of a 2,544-row experiment (~9% of the corpus) |
| **confidence** | certain on the defect; the remedy needs validation |

**The claim.** The curation z-scored each (file × cell type) column separately —
14 blocks — on the strength of the authors' warning that values "cannot be cross
compared across files". But that warning is about *across files*. **Within** a
file the authors themselves cross-compare cell types: `dataset_info.txt` carries
an `nLnLE_Selectivity` column, and the Methods define a Selectivity Score as
on-target minus off-target `nLnLE`. Z-scoring each cell-type column separately
forces every column to mean 0 / sd 1 and therefore deletes exactly the quantity
the authors computed.

**Reproduced, from `5_Norm_RAMOS_HEPG2_Dual_Objective_4ITER_SM102_DSPC.csv`
(block E5, n = 180) against `data/data_vitro.csv`:**

| | raw source | stored in corpus |
|---|---|---|
| Ramos mean (sd) | 2.5306 (1.2591) | 0.0000 (1.0000) |
| HepG2 mean (sd) | 6.5648 (2.4343) | 0.0000 (1.0000) |
| mean HepG2 − Ramos | **+4.0342 ln = 56.5×** | **exactly 0.0000** |
| per-formulation difference | −0.837 … +8.736 (up to 6,225×) | — |
| HepG2 > Ramos | 173 of 180 formulations | — |

Two consequences, the second of which the reviewer did not report:

1. **The difference is not merely shrunk, it is identically zero** — a
   tautology of standardising each block, so no check on per-block mean or sd
   can ever see it. This is error class 22 recurring, the same defect found in
   `ql_2026` three days ago, in a different experiment.
2. **Spearman(raw per-formulation difference, stored z difference) = −0.517.**
   The selectivity ordering is not just destroyed, it is *anti-correlated*: a
   formulation that is genuinely more HepG2-selective tends to look more
   Ramos-selective in the corpus. I measured 92 of 180 formulations with the
   sign of (HepG2 − Ramos) inverted (the reviewer reported 80; the exact count
   depends on an assumed row alignment, the −0.517 and the zero mean do not).

**Why no check caught it.** `make_merge_subset.py` validates that the
concatenated standardised blocks give pooled mean 0 / sd 1 — which is true by
construction for any set of standardised blocks, including wrong ones. The
project's own ledger already says this: *"per-block mean 0 / sd 1 is a
tautology"*.

**Proposed remedy (NOT applied).** One z-score per source **file** (7 blocks,
not 14), preserving within-file cross-cell-type structure while still honouring
the across-file warning; plus a real check — that the stored z difference
reproduces the printed `nLnLE_Selectivity` column up to a single per-file
constant. Needs validation before application: it changes merged values, and
the choice of normalisation unit is a research decision, not a bug fix.

---

## F-004 — `WT_2026`: the URS source file states each formulation twice, and the two disagree

| | |
|---|---|
| **found by** | WT_2026 reviewer (its E2) |
| **status** | `confirmed` — reproduced, and the scope narrowed by my own test |
| **severity** | **changes a value**: `IL_to_nucleicacid_chargeratio` on 129 rows |
| **confidence** | certain on the defect; the repair direction needs validation |

The authors' URS file
(`0_Norm_RAMOS_Single_Objective_UniformRandomSampling_Ctrl_MC3_DSPC.csv`)
carries each formulation in two parameterisations: four design parameters
(`NP_ratio`, `IL+HL`, `HL_IL+HL`, `PEG_PEG+Chol`) and four molar percents. They
do not describe the same formulations.

Deriving percents from the design parameters
(`HL = (IL+HL)·(HL_IL+HL)`, `PEG = (100−(IL+HL))·(PEG_PEG+Chol)`, etc.):

| file | n | design block agrees with its own percent block | derives *some other* row's percents |
|---|---|---|---|
| **0 (URS)** | **129** | **0** | **126** |
| 0 (E0) | 273 | 273 | — |
| 1 (E1) | 197 | 197 | — |
| 2 (E2) | 136 | 136 | — |
| 3 (E3) | 123 | 123 | — |
| 4 (E4) | 180 | 180 | — |
| 5 (E5) | 180 | 180 | — |

Median disagreement on the URS file is **30.7 pp**, max 63.2 pp. The mapping is
a **scrambled permutation, not a shift** (offsets +12, +26, +108, +5, −20, −29 …)
and it is **injective**, so each design row resolves to exactly one formulation
and the true `NP_ratio` is recoverable for 126 of 129.

**Scope correction to the reviewer's report.** It reported this as a property of
"the URS file" and tested only that file. I tested all seven: six are clean, so
129 rows is right. My own first pass flagged file 1 as worse (9/197) — that was
an artifact of a 0.01 pp tolerance against percents the authors rounded to one
decimal. At 0.15 pp it is 197/197 clean. **Recorded because it is exactly the
kind of false finding this ledger exists to catch before it becomes an edit.**

**Consequence.** The curation took molar percents from the percent block (right,
per the reviewer's reading of Suppl Table 5) but `N:P` from `NP_ratio`, which
travels with the scrambled block. So `IL_to_nucleicacid_chargeratio` is wrong on
129 corpus rows: **median error 6.86, mean 7.70, max 21.54**, against a file
range of 4.01–25.54 — i.e. a median error of about a third of the full range,
and 81 of 126 rows off by more than 5. C4, C5 and C6 are all blind to it,
because each row's percents still sum to exactly 100.

**Not applied.** Which block is authoritative is precisely the `AC_2025`
question — repairing the cheaper half assigns values to the wrong formulations.
The reviewer's Suppl Table 5 anchor needs independent confirmation before
anything is changed, and the honest alternative is to blank `N:P` on those 129
rows rather than reconstruct it.

---

## F-005 — my packet generator invented a finding

| | |
|---|---|
| **found by** | WT_2026 reviewer (its E3) |
| **status** | `fixed` |
| **severity** | cosmetic in the data; **wasted reviewer effort and fabricated an absence** |
| **confidence** | certain |

`brief.py`'s `anchors_for()` guessed at `comparisons`/`rows` keys in
`structure_check.json`, whose real shape is `results`; `reagents_in()` guessed
at `reagents`/`registry` in `reagent_registry.json`, same problem. Both
degraded silently to empty, so **every one of the five packets** told its
reviewer "No anchor records found for this folder… Either none was created —
which is itself a finding" and showed `{}` for C11 reagents.

In truth there are 14 anchors across the five folders (WT_2026 has 2, both exact
matches) and `WT_2026` alone uses 6 registered reagents on 2,544 rows.

Fixed to read the real keys and to **raise** rather than return empty if the
shape changes, because a packet that invents an absence is worse than one that
omits a section. Packets regenerated. **Discount any "missing anchor" or "no
registered reagents" finding in the other four Wave-1 reports** — it is my bug,
not theirs.

Lesson for the program: a brief generated from live data still has to be
*verified against* live data. Silent `.get()` fallbacks in a briefing tool
manufacture findings.

---

## F-006 — NEW ERROR CLASS: `LL_2026`'s hit-rate validation is calibrated on the quantity it validates

| | |
|---|---|
| **found by** | LL_2026 reviewer (its E1) |
| **status** | `confirmed` — verified by direct code read and by the script's own recorded numbers |
| **severity** | **no corpus value changes**, but a Tier-1 validation claim is false, and it is one of the two claims the project cites most |
| **confidence** | certain |

`curation/lin_aromatic_2026/digitize_fig2fgh.py`, lines 118-123:

```python
heights_l = np.array([b["height"] for b in lb], float)
digit = marg(allrows, field, cats)
rep_l  = np.array([digit[c]["gt_1e5"] for c in cats])   # the Fig 2e hit rates
slope  = float(np.sum(heights_l * rep_l) / np.sum(rep_l * rep_l))
units_per_tick = px_tick / slope if slope else float("nan")
```

`slope` — the pixels-per-percent that converts **every** bar height — is a
through-origin least-squares fit of the measured bar heights onto `rep_l`, which
is the Fig 2e heatmap hit-rate set this panel is supposed to *check*. The tick
spacing `px_tick` is measured, but only to print the diagnostic
`implied_percent_per_tick`. **The file contains zero `assert` statements.**

Confirmed numerically from the script's own `08_fig2fgh_qc.json`:

| panel | `px_per_tick`/20 (tick-implied) | `px_per_percent` actually used | discrepancy |
|---|---|---|---|
| 2f | 5.100 | 5.046 | +1.07% |
| 2g | 5.150 | 5.071 | +1.56% |
| 2h | 5.300 | 5.293 | +0.13% |

If the conversion were tick-anchored those columns would agree. They do not, so
the scale is demonstrably the fitted one.

**What it costs.** Scaling every digitized rate by a constant *k* rescales
`slope` by 1/*k* and the measured percentages by *k*, so the reported deltas are
merely rescaled instead of blowing up. **A uniform multiplicative error in the
heatmap read is absorbed and invisible** — exactly the failure mode of class 20,
which this check was meant to close. It is weakened, not worthless: the 7-, 5-
and 12-point *shape* comparison survives, so a row or column permutation would
still show up.

**The data is fine.** The reviewer re-measured all 48 bars against a genuinely
tick-anchored fit (5, 5 and 6 tick centres, residuals ≤2 px; scales 5.090 /
5.125 / 5.350) and got **48/48 within 0.75 pp**. So nothing needs re-digitizing;
the method and the claim need fixing.

**This is a new error class**, distinct from 20 (wrong anchor) and 21 (docs
contradict code): *a validation whose calibration is fitted to the quantity
being validated.* It cannot fail. No check in the registry covers it, and
nothing in `audit.py` could — it is a property of a curation's own validation
code, not of the corpus. Proposed detection rule: in any digitization, the
quantity converting pixels to units must derive **only** from axis furniture
(ticks, labelled gridlines, printed anchors) and never from the comparison
target. A reviewable signature is a fit whose x-vector is the expected answer.

**Three places repeat the false claim** and must be corrected together:
`11_LL_2026_qc.json` `validation.fig2fgh_hit_rates`,
`lit_review/build_data.py` (the `41689530` ADDED entry), and `CLAUDE.md` §5
("Lin gave 48 hit rates plus 4 HRMS anchors").

---

## F-007 — `LL_2026`: the one hit rate that appeared to disagree was a missed bar

| | |
|---|---|
| **found by** | LL_2026 reviewer (its E2) |
| **status** | `confirmed` — the record reads `[29.5, 0.0]` as reported |
| **severity** | metadata only |
| **confidence** | certain |

`08_fig2fgh_qc.json` panel `2h`, category `EP6`: `"measured": [29.5, 0.0]`,
`"delta": [-0.5, -1.4]`. The `>1e6` bar is recorded as **absent**. It is not
absent — the reviewer locates it at x 738-787, colour from y≈4048, saturated by
y≈4052 against a baseline at 4057.5, giving 1.0-1.8%, which brackets the
heatmap's 1/70 = 1.43%.

Root cause: `bars_in(dark, …, min_area=800)`. The bar's component is ≈50×9 ≈ 450
px² and is dropped by the area floor, after which `m6` falls through to the
literal `0.0`. So the curation's quoted `max |delta| 1.4 pp` should be ≈0.75 pp.

**Why this matters more than the number.** The project's own rule is *never
excuse the one that disagrees* — written into `ERROR_CLASSES.md` three days ago
after the `ql_2026` Fig 3c error. This was the only apparent disagreement in the
set, it was recorded, it propagated into the paper-level QC summary, and it was
not chased. It turned out benign, which is the lucky outcome, not a vindication.
A silent `else 0.0` fallback is how a measurement failure disguises itself as a
measurement.

---

## F-008 — `LL_2026`: Fig 2a sampling window extends 57 px past the heatmap

| | |
|---|---|
| **found by** | LL_2026 reviewer (its E3) |
| **status** | `claimed` — not yet reproduced by me |
| **severity** | cosmetic today; latent off-by-N |
| **confidence** | reviewer: certain |

`digitize_fig2a.py` line 42 sets `BODY = (640, 679, 3161, 1073)`; the reviewer
measures the right edge at x≈3104 (white from 3109) and finds 22 boundaries at a
uniform pitch of 117.33, against the script's assumed 120.05. The window centre
drifts right by `2.72·c + 1.35` px, so by column 20 it is 55.7 px off a cell of
half-width 58.7 and 28 of 60 sampled columns lie outside the heatmap.

All 147 bins are unchanged under both geometries — but the tell is the
colour-distance statistic: **2.29 under the script's geometry, 0.00 under the
corrected one.** `CLAUDE.md` §5 says "anything above ~5 means the geometry is
wrong, not the colours"; 2.29 passed that gate while being pure geometry error.
The gate's threshold is calibrated for catastrophic offsets and is blind to a
systematic drift that eats half the sampling window. Proposed corrected value:
`BODY = (640, 675, 3104, 1073)`.

---

## F-009 — `LL_2026`: an enantiopure reagent's stereocentre is dropped on 141 rows

| | |
|---|---|
| **found by** | LL_2026 reviewer (its E7) |
| **status** | `claimed` — not yet reproduced by me |
| **severity** | **structure** (stereochemistry only; connectivity unaffected) |
| **confidence** | reviewer: probable |

The SI Materials list names reagent `5B` as **"(S)-glycidyl trityl ether"** —
enantiopure. Amine opening occurs at the terminal CH₂, leaving C2 untouched, so
the carbinol in every `5B` arm is configurationally defined. The corpus carries
it unassigned on 141 rows (120 `LL_2026`, 21 `LL_2026_2`). `blocks.csv`'s
`ambiguity_note` discloses it; nothing on the rows does, and C9 inspects only
C=C and C=N.

The reviewer explicitly declines to recommend a rewrite: the paper reports no
optical rotation or chiral HPLC **of the products**, every other carbinol in the
library genuinely is racemic, and asserting stereochemistry on one of six arms
has its own cost. That is the right posture — the deliverable is a marker, which
F-001 shows the schema cannot presently carry. **Do not repair silently in
either direction.**

---

## F-010 — C2 has no teeth on the largest curation in the corpus

| | |
|---|---|
| **found by** | LL_2026 reviewer (its E9) |
| **status** | `confirmed` — `audit_findings.json` names `"LL_2026": 140` in the skip list |
| **severity** | metadata only (check coverage) |
| **confidence** | certain |

C2 (mass monotonicity within a homologous series) skips **all 140** `LL_2026`
series because the tail token is `EP10`/`Ac12` rather than a bare integer;
`LL_2026_2` has no series at all. The corpus-wide audit line reads "235 verified
-length series checked, 0 violations (737 series skipped)" — so 76% of series
are skipped, and the single largest curation contributes zero checked series
while reporting a clean pass.

The reviewer ran the check manually: **140/140 monotone**, EP6→EP18 increments
exactly 112.1252 Da = 4 tails × 2 CH₂. No bug in the data. But this is the
failure mode `REVIEW_PROGRAM.md` §3 calls the worst available — a check that
passes because it did not look. `CLAUDE.md` §5's precondition for extending C2
(a naming convention verified against source metadata) is met here: the reviewer
verified all 30 reagent names against the SI Materials list, 30/30.

**Corroboration of F-001, arrived independently.** This reviewer was not told
about F-001 and reached the same conclusion for its own paper: all quality flags
die at the merge — `positional_ambiguity` on **987** rows (every structure in
this slice is a convention about which nitrogen bears which arm, since the
library was made in one pot with no chromatography), `confidence = medium` on
**140** rows (`Ac8-1`/`Ac10-1`, commercial iso-octyl/iso-decyl acrylate, branched
isomer mixtures), plus `match_distance`, `runner_up_margin`, `bin_label`,
`value_scale` and the arm counts. Same class as `18PG`: disclosed upstream, lost
downstream.

It also found `IL_head_name` / `IL_linker_name` / `IL_tail*` blank on all 987
rows although the decomposition exists in `03_enumerated.csv` — against
16,471 of 28,885 vitro rows that do carry it, including `QL_2026` at 344/344.

---

## F-011 — `ZT_2026`: 24 low-confidence structures were silently upgraded, half of them to `high`

| | |
|---|---|
| **found by** | CORPUS-B (its B-2) |
| **status** | `confirmed` — reproduced in full, and it is worse than reported |
| **severity** | **no value or structure changes, but it converts a disclosed assumption into a false assurance** on 24 corpus rows |
| **confidence** | certain |

`curation/siegwart_42634137/03_enumerated.csv` flags 24 of its 36 lipids
`confidence = low`, with `quaternization_site_assumed = 1` and a verbatim
reason in `quaternization_site`:

> `ASSUMED: no ring N and the paper draws no structure`

In `06_merged.csv`, one stage later **in the same folder**, not one `low`
survives. The distribution goes from `{low: 24, medium: 6, high: 6}` to
`{high: 18, medium: 18, blank: 15}`, and the `quaternization_site_assumed`
column is gone entirely.

The 24 map mechanically by cholesterol tail:

| lipids | upstream | downstream |
|---|---|---|
| `C1O1B1`…`C4O1B3` (12) | `low` | **`medium`** |
| `C1O2B1`…`C4O2B3` (12) | `low` | **`high`** |

`O1` → medium, `O2` → high, with no reference to the lipid-level evidence. That
regularity is the signature of a **join that overwrote the column** — the
downstream `confidence` is evidently the *tail block's* confidence, merged in
under the same column name and silently replacing the lipid's own. All 24 are
present in `data_vivo.csv` under `ZT_2026`.

**Why this is the most serious finding of the round so far.** F-001 is a flag
that failed to propagate — bad, but it leaves no false claim behind, and a
careful consumer finds nothing and knows nothing. This is different in kind: the
rows whose regiochemistry the paper **draws nothing for** are the rows now
labelled `high`. The only surviving provenance asserts the opposite of the
truth, and it contradicts `CHANGELOG.md` §8, which records these as low
confidence. A reviewer or reader who checks the curation folder for the weakest
structures is pointed away from exactly the right 24.

Note the interaction with F-001: because no flag column reaches the corpus at
all, this upgrade is currently invisible from `data/`. Fixing F-001 by
propagating `confidence` **without** fixing this first would publish `high` on
24 assumed structures — strictly worse than publishing nothing. **F-011 must be
repaired before F-001.**

**Repair (not applied).** Find the merge that overwrites `confidence`, suffix or
rename the colliding column as the stage already does for `chol_tail_x`,
`source_pdf_x`/`_y`, and carry `quaternization_site_assumed` through. Then the
24 rows read `low` again. Needs validation: I have identified the signature of
the collision but not yet the line that causes it.

**New error class (29): a join silently overwriting a provenance column with a
same-named column from the other table, inverting a disclosed low confidence
into a high one.** Distinct from class 21 (docs contradict code) — here the data
contradicts itself across two stages of one folder, and the merge suffixing
convention (`_x`/`_y`) that would have exposed it was applied to five other
columns but not this one.

---

## F-001 addendum — CORPUS-B's independent count, and a refinement that matters

CORPUS-B was not told about F-001 and reached the same conclusion. Its count is
sharper than mine: **1,137 corpus rows carry a convention rather than a
measurement and say nothing about it.**

| experiment | rows | flag lost |
|---|---|---|
| `LL_2026` | 840 | `positional_ambiguity = 1` (140 of them also `confidence = medium`) |
| `LL_2026_2` | 147 | `positional_ambiguity = 1` |
| `QL_2026` | 117 | 113 floored + 3 censored + 1 conflict |
| `ZT_2026` (vivo) | 24 | `quaternization_site_assumed = 1` — and see F-011 |
| `XD_2026` | 9 | `confidence = medium` (E1–E9, no external anchor) |

**The refinement.** I had recorded this as flags not surviving the merge. That
is not quite right, and the distinction changes the repair: **none of the six
files `build_dataset.py`'s `MERGE` reads contains a flag column at all.** The
flags are lost one stage earlier, inside each curation folder, when it emits the
72-column schema. No change to the builder could have preserved them. So the
repair is in the schema plus every curation's emit step, not in the merge — and
F-011 shows at least one folder actively corrupts the flag before that point.

My earlier figure of "≥266 rows" counted only the flags I could enumerate from
QC JSONs; it missed `LL_2026`'s 987 `positional_ambiguity` rows, which live in
`03_enumerated.csv` rather than a QC file, and `ZT_2026`'s 24. **1,137 is the
better number.** Corrected here rather than silently in place.

---

## F-012 — the Tier-1 audit cannot fail: nine of eleven checks are print-only

| | |
|---|---|
| **found by** | CORPUS-C (its E1), via 43 seeded mutations across 15 classes |
| **status** | `confirmed` — reproduced directly |
| **severity** | **the project's entire Tier-1 posture is unenforced**; no corpus value is wrong because of it |
| **confidence** | certain |

`provenance/audit.py`, lines 386-390:

```python
if __name__ == "__main__":
    main()                                                   # <- return value DISCARDED
    rc  = delegate("structure_check.py", "C10  drawn-structure anchors")
    rc |= delegate("reagent_registry.py", "C11  named commercial reagents")
    sys.exit(rc)
```

Exit status is derived **only** from the C10 and C11 delegates. C1-C9 print and
are then forgotten. `provenance/build_dataset.py` is worse: it never imports
`sys` at all (`grep -c '^import sys'` → **0**), so the rebuild also always
exits 0.

Verified on the live corpus: both commands exit **0**. CORPUS-C seeded eight
simultaneous violations into a copy — disconnected `.[Cl-]`, stale
`Molecular.Weight`, ratios summing to 109, `Experiment_value='nan'`, wrong
`Model`, wrong `Route`, blanked `IL_name`, an exact duplicate row — and reports
all eight printed with **exit code 0**.

Of 43 mutations, 30 fired in the printed report and 13 did not; but **only 5
cause any command to exit non-zero, all of them in C10/C11.**

**Why this matters more than the miss count.** `CLAUDE.md` §6 presents Tier 1 as
"automated, every build, free — run these two commands", and every session so
far (mine included, repeatedly, today) has run them, read a clean report, and
proceeded. That green light is hardcoded for C1-C9. The brief given to CORPUS-C
said it in advance: *"A check that prints a problem into a 60-line report that a
human skims is not a check."* It was right.

**Why the obvious fix is wrong, and this needs validation.** Making `main()`
propagate a failure count would make the audit fail permanently on the corpus's
**known and deliberately accepted** conditions — 48 duplicate SMILES (39 with
conflicting values), 212 rows with unassigned stereo C=C, 2,544 blank
`Publication_PMID` on the `WT_2026` preprint, `SB_2024`'s 3 rows whose ratios do
not sum to 100, the one accepted exact-duplicate row in FALCON. A check that
fails always is as useless as one that never fails, and is probably why it was
left print-only in the first place.

So the repair is a **baseline**: each check reports a count, the accepted count
is recorded explicitly with the reason it is accepted, and the command fails
when a count *departs from its baseline in either direction* — a new violation
or a silently vanished known one. That is a design decision about how the
project's quality gate works, not a bug fix, so it goes to the user and the
validation agent rather than straight into the file.

**Related, not yet reproduced by me:** CORPUS-C reports 13 of 43 mutations not
detected even in the printed report. Those are the coverage gaps and matter
independently of the exit-code defect; they need to be triaged individually
before any baseline is set, because a baseline over a blind check just freezes
the blindness.

---

## F-013 — **ALL 537 `ZT_2026` rows carry epi-cholesterol (3α), not cholesterol (3β)**

| | |
|---|---|
| **found by** | ZT_2026 reviewer (its E1) |
| **status** | `confirmed` — effect, scope AND mechanism all reproduced by me independently |
| **severity** | **changes a structure on every row of a 537-row experiment** (stereochemistry; connectivity unaffected) |
| **confidence** | certain |

**The effect.** Cleaving the steroid ester off the stored `C6O2B2` SMILES in
`data/data_vitro.csv` gives a sterol whose InChIKey is

```
HVYWMOMLDIMFJA-VEIPTCAHSA-N     <- epicholesterol (3alpha-ol)
HVYWMOMLDIMFJA-DPAQBDIFSA-N     <- cholesterol, what it should be
```

Same skeleton block, different stereo block. **All 51 unique `ZT_2026`
structures put steroid C3 at CIP `R`**; cholesterol and every 3β-cholesteryl
ester are `S`. Scope: 486 vitro + 51 vivo = **537 rows, i.e. the entire
experiment**, plus `curation/new_curation.csv`, `03_enumerated.csv` and
`reference_structures.csv`.

**The mechanism, which I reproduced from scratch.** Feeding *genuine*
cholesterol (`-DPAQBDIFSA-N`, verified on input) through the
`AddBond`-then-`RemoveAtom` pattern used in `enumerate.py::_cholesteryl_ester`:

```python
frag.AddBond(ester_o, off + c3, Chem.BondType.SINGLE)
frag.RemoveAtom(off + oh.GetIdx())      # <-- inverts the tag at C3
```

| route | steroid C3 |
|---|---|
| `RWMol` AddBond then RemoveAtom | **R** (wrong) |
| the same 3β ester written as a SMILES string | **S** (right) |

`RWMol.RemoveAtom` deletes the hydroxyl oxygen from the middle of C3's
neighbour list **without correcting C3's chiral tag**, so the stored tag now
means the opposite of what it meant. The input was right; the assembly inverted
it.

**Why every check was blind.** Formula `C117H218N5O12+` and exact mass
1885.6597 are identical between 3α and 3β, and so are **all 17 descriptor
columns** (Heavy.Atoms, Rings, RotatableBonds, TPSA, LogP, MolMR, FractionCSP3 —
bit-identical). C1–C8 test internal consistency and cannot see it. C9 inspects
only C=C stereo. C11 pins `Cholesterol` as a *component*, and this cholesterol
is inside the ionizable lipid, not in `CHL_SMILES`.

**And C10 — the one check designed for exactly this — was blind too, because
the anchor is a copy.** `reference_structures.csv`'s `reference_SMILES` is
**byte-identical** to `03_enumerated.csv`'s `C6O2B2` (verified: string equality,
not canonical equality). Its own note claims it was *"assembled from fragments,
independently of enumerate.py … through the C3 oxygen (3-beta)"*. It stores 3α.
So `structure_check.json` reports `exact_match` **while comparing a value to
itself**. This is the "C10 anchor copied from the enumeration (can only pass)"
hazard already written in `ERROR_CLASSES.md` — and here it is, realised, on the
project's flagship curation.

**The bitter part.** Fragment assembly was *mandated* by `CLAUDE.md` §6 as the
fix for an earlier bug (a SMILES ring-closure collision that fused piperazine to
cholesterol). The prescribed remedy introduced a different silent failure mode
in the same molecule. The lesson is not "go back to strings" — it is that
**both** routes need an external stereo-aware check, which is what the repair
below gives.

**Repair, and the check that closes the class permanently** (not applied):
after regenerating, cleave the steroid ester of every `ZT_2026` `IL_SMILES` and
assert the liberated alcohol's InChIKey equals `HVYWMOMLDIMFJA-DPAQBDIFSA-N` —
a string **already pinned** in `provenance/reagent_registry.json` under
`Cholesterol`. Two lines, and it would have caught this on day one. More
generally: C10 must compare the **full** InChIKey including the stereo block,
and must refuse an anchor that is string-identical to the enumeration it checks.

**Also requires:** regenerating `curation/new_curation.csv` (537 rows),
`03_enumerated.csv`, a genuinely independent `reference_structures.csv`, and
rebuilding. This changes 537 of 28,885 vitro + 1,780 vivo rows and **no
`Experiment_value`** — every value reproduced independently to 1e-15.

---

## F-014 — `ZT_2026`: `MC3-CB1/2/3` carry the KC2 dioxolane under an MC3 name

| | |
|---|---|
| **found by** | ZT_2026 reviewer (its E3) |
| **status** | `confirmed` |
| **severity** | **changes a structure** (3 vivo rows) |
| **confidence** | certain |

Measured directly:

| | rings | 1,3-dioxolane matches |
|---|---|---|
| `MC3-CB1/2/3` | 5 | **1** |
| corpus `DLin-MC3-DMA` (3,097 rows) | 0 | 0 |
| corpus `KC2` (203 rows) | 1 | 1 |

4 steroid rings + 1 dioxolane = 5. DLin-MC3-DMA is a **linear ester** with no
core ring; the dioxolane is DLin-KC2-DMA. So these three conjugates were built
on the wrong parent.

**Root cause is our own fix sequence.** `build_new_curation.py` pulls the parent
out of the live corpus *by name* (`PARENT = {"MC3-C": "DLin-MC3-DMA", …}`).
`ZT_2026` was built 2026-10-28, **before** the 2026-10-29 `AC_2025` repair that
put the MC3 and KC2 names back on the right structures. `AC_2025` was
regenerated afterwards; these three derived rows were not.

**Which half to repair — and note it is the opposite of `AC_2025`.** Fig 1d
labels the measurement `MC3-C` and Fig 1b draws MC3 as the linear ester, so here
the **name is right and the structure is wrong**: rewrite the structure. In
`AC_2025` the structures were right and the names wrong. Same symptom, opposite
repair, settled only by going to each source — exactly the lesson
`ERROR_CLASSES.md` records. **Do not pattern-match the two.**

**New error class (30): a curation that reads a component from the live corpus
by name is silently invalidated when that name is later repaired.** Any
downstream artefact built from corpus lookups needs regenerating after a
`NAME_PATCHES`/`COMPONENT_PATCHES` change; nothing currently tracks that
dependency.

**Also reported, not yet reproduced by me (F-015):** the second copy of the
quaternization rule in `curation/build_new_curation.py` still uses
`max(atom index)` while its docstring claims it matches the repaired
`enumerate.py`; the reviewer reports this alkylates `Lipid 14-CB*` (3 rows) at
the hindered central nitrogen rather than the terminal dimethylamino. This is
the *second copy* hazard that `CLAUDE.md` §6 Tier 3 was written for, recurring
in the same file it was written about.

---

## F-016 — Class 22 is systemic: **5,040 rows (17.4% of vitro)** are per-cell-type z-scored, not per-experiment

| | |
|---|---|
| **found by** | CORPUS-A (its E1/E4/E5); scope established by my own corpus-wide scan |
| **status** | `confirmed` |
| **severity** | **changes what every value in 5,040 rows means**; bears directly on manuscript Results 3 and 6 |
| **confidence** | certain |

The corpus convention (`CLAUDE.md` §3) is **a per-experiment z-score of log10**.
In three experiments it was applied per *(experiment × cell type)* instead. I
scanned all 61 experiment blocks for the signature — every per-cell-type block
independently at mean 0 and population sd 1, which cannot happen under pooled
standardisation except by coincidence:

| experiment | rows | blocks each at mean 0 / pop sd 1 |
|---|---|---|
| `YX_2024` | **2,400** | HeLa 1,200, RAW264.7 1,200 |
| `WT_2026` | **2,544** | Ramos 959, HepG2 513, THP-1 409, DC2.4 180, NIH-3T3 180, C2C12 180, INS-1E 123 |
| `LZ_2025` | **96** | A549 48, HepG2 48 |
| | **5,040** | = 17.4% of `data_vitro.csv` |

Only these three. Every other multi-cell-type experiment is pooled correctly.

**`YX_2024` is the sharpest case, and it is inherited upstream data, not ours.**
All 1,200 lipids are measured in **both** assays — RAW264.7's SMILES set is an
exact subset of HeLa's, 1,200/1,200 shared. So the quantity destroyed is a
**pure within-chemistry, between-assay offset**, with chemistry held perfectly
constant:

```
mean(HeLa z - RAW264.7 z) over the same 1,200 lipids = +0.0000
Spearman(HeLa z, RAW264.7 z) over the same chemistry = 0.461
```

The mean difference is not small — it is **identically zero**, forced by the
transform, exactly as in `ql_2026` and `WT_2026`.

**Why this matters to the manuscript, not just the data.** `YX_2024` is the
AGILE dataset and it feeds Result 3. The project's load-bearing result is that
*identical lipids rank near-independently across labs* (pooled Spearman 0.083
against a within-source replicate ceiling of 0.400). `YX_2024` happens to be the
cleanest available instance of the **intermediate** rung — same lab, same
library, only the cell type changes — and it reads **0.461**, comfortably above
the 0.400 replicate ceiling and far above the 0.083 cross-lab figure. That is a
genuinely useful datum for Result 6 and a strong answer to the "cross-lab ≈0 is
just assay noise" objection, directly comparable to the `molea.csv` 0.30 already
recorded in `CLAUDE.md` §8. But the *magnitude* of the between-assay shift has
been deleted, so any analysis of absolute potency across those two assays is
reading a zero that the normalisation inserted.

**CORPUS-A also dismantled the argument that said this was fine.** The
builder/README/CHANGELOG justify treating `YX_2024` as one pooled 2,400-row
z-score on the grounds that "P's 1,200 values are bit-identical to V's, so V
pooled over 2,400". That is a **tautology**: if the blocks were always separate,
P's HeLa subset is unchanged either way, so the observation cannot discriminate.
I confirmed the conclusion by the independent route above rather than relying on
the reviewer's upstream-file arithmetic.

**Not applied — and this one is genuinely a research decision.** Repooling
changes the meaning of 5,040 values, including 2,400 rows of inherited upstream
data that the ML companion project also reads, and it would change Result 3's
inputs. The options are not symmetric:

1. **Repool per experiment** — matches the stated convention, restores the
   between-assay offset, changes 5,040 values and any frozen number computed
   from them.
2. **Keep per-cell-type and document it** — then `Experiment_value` does not
   mean what §3 says it means, and the *(experiment × cell type)* pair, not
   `Experiment_ID`, is the real normalisation unit. That has consequences for
   what "source" means in a paper about between-source effects.
3. **Carry both columns** — honest, and lets Result 6 use the offset while
   Result 3 keeps a gauge-invariant ranking.

**Class 22 needs an automated check, which is the real lesson.** This is now the
**fourth** instance (ql_2026 fixed, WT_2026 F-003, YX_2024, LZ_2025) and each
was found by a human noticing, never by a check. The scan I ran above is the
check: for every experiment with more than one `Model_type`, flag it when each
block independently sits at mean 0 / pop sd 1. It is ~15 lines, it has no false
positives on the current corpus, and it would have caught all four. **This is
the single highest-value check to add**, and it should go in before any of the
repairs so that each repair can be verified by it flipping.

---

## F-017 — `QL_2026`: class 20 is still live in Fig 2c, the sibling panel of the one we repaired

| | |
|---|---|
| **found by** | QL_2026 reviewer (its E4) |
| **status** | `confirmed` — by arithmetic on the stored anchors |
| **severity** | **changes a value** on 27 rows (≈1.7% high) |
| **confidence** | certain |

`digitize_fig2c.py` carries `Y_1E4 = 3156.0`, `PX_PER_DECADE = 407.0`, with the
comment *"y-axis label centroids: 1e6 2342, 1e5 2748, 1e4 3156"*.

Those centroids give decade gaps of **406 and 408 px** — unequal. A log axis has
exactly equal decade spacing, so unequal gaps are proof the anchors are **text
glyph centroids, not tick marks**, and `407.0` is simply their mean. The
reviewer's tick detection gives exactly **203.0 / 203.0** native = 406.0 at
600 dpi, with 10⁴ at 3153.0 — 3 px above the stored anchor.

Three days ago we repaired exactly this error in `digitize_fig3c.py` and wrote
it into the registry as class 20, whose stated repair is *anchor on tick marks*.
**The sibling panel in the same folder was never re-checked.** A fix applied to
the file where a bug was found, and not to its siblings, is the same hazard as
the duplicated-code class that `CLAUDE.md` §6 Tier 3 exists for.

## F-018 — `QL_2026`: replicate dots are still read as the bar top on 15 of 60 bars

| | |
|---|---|
| **found by** | QL_2026 reviewer (its E1) |
| **status** | `claimed` — mechanism is plausible and specific; I have not re-measured the rasters |
| **severity** | **changes a value** on 14 rows (up to +26%) and biases the cross-library offset |
| **confidence** | reviewer: certain |

`PROVENANCE.md` states this bug was found and fixed. The fix was partial: the
`w//8` inner trim leaves a 75-78% window, and three overlapping dots can fill
≥85% of *that* while covering only 55-70% of the full bar. Reported worst cases:
`T3A-L3-C6` +26%, `110-L7-C10` +25%, `T3A-L3-2C10` +24%.

Two independent tells support it and both are cheap to check:
- bars with no dot overlap agree with the curation **to the pixel**, which
  isolates the cause to dot contamination rather than a different convention;
- the stored file gives **two pairs of bars byte-identical values**
  (`T3A-L1-2C10` = `T3A-L7-2C10` = 933,681.0; `T3A-L5-2C10` = `T3A-L6-2C12` =
  1,346,291.4). Two exactly-equal bars in one panel is the signature of a
  detector latching onto something that is not the rule.

Consequence for the offset: the four shared benchmarks give **1.073** (sd 0.034)
rather than the stored 1.1116, because ALC-0315 and Lipid 5 are dot-contaminated
in Fig 3c and ALC-0315 again in Fig 2c.

**Every printed Fig 3d/3e hit rate and the "17 lipids" sentence still reproduce
on the corrected values.** So the paper's own statistics cannot discriminate
this error — it has to be caught geometrically. That is an important limit on
the "reproduce every printed statistic" rule we added three days ago: passing it
is necessary, not sufficient.

Net effect of F-017 + F-018 together, as the reviewer measured it: **rms Δz
0.054, max |Δz| 0.209, Spearman 0.99998, largest rank shift 4 of 344.** Real,
worth fixing, and small beside F-019.

## F-019 — `QL_2026`: the colormap floor sits ABOVE the feasible range; 48 cells are upper bounds

| | |
|---|---|
| **found by** | QL_2026 reviewer (its E7) |
| **status** | `claimed` |
| **severity** | **changes a value**; moves individual z-scores by up to ~1.0 |
| **confidence** | reviewer: certain on the clip, unverifiable on the true values |

113 of 320 library-1 cells take just two values (1,310.3 n=48, 1,724.1 n=65).
The reviewer finds the colormap reaches **pure white at ≈900 RLU** while the
bar's zero is lower still, so every true value from 0 to ~900 renders identical
— the 48 pure-white cells are **upper bounds, not measurements**, and the
assigned 1,310 RLU is *above the entire band they could occupy*.

Sensitivity, holding everything else fixed and varying only what those 113 cells
really are: rms Δz 0.054 (as merged) → 0.415, max Δz 0.21 → **0.99** at 100 RLU.
**This is an order of magnitude larger than every digitization error combined**
and is the dominant uncertainty in the slice. The between-library gap is robust
(library-2 mean z falls only 2.09 → 1.70), so the paper's finding is safe;
individual library-1 rows are not.

Related, and nothing records it: **344 rows carry only 83 distinct
`Experiment_value`s, and 234 share a value with ≥10 other rows**; the three
largest ties are 65, 57 and 48 rows. A model trained on this slice sees a third
of library 1 as a handful of tied pseudo-measurements. That is a data-quality
fact with direct bearing on Result 6.

## F-020 — `XD_2026`: bar tops read inside the bar where replicate discs occlude the width

| | |
|---|---|
| **found by** | XD_2026 reviewer (its E-1) |
| **status** | `claimed` |
| **severity** | **changes a value** on all 12 rows (mean +0.060 decades, up to +0.122) |
| **confidence** | reviewer: certain |

`stage04_fig5c.py` declares the bar top at the first row where >60% of the bar
width is fill. The three replicate discs straddle the top edge and occlude the
middle of the width, so the detector fires several pixels *inside* the bar. The
bias is one-sided and non-constant, so it does not cancel in the z-score:
`Experiment_value` changes by up to 0.060 (rms 0.034).

Independent confirmation in the report: measuring from the flank columns (never
occluded; left and right agree to ≤1 px on all 12 bars) gives tops 0.5-10.9 px
higher, and re-applying the curation's own 60% rule reproduces the stored tops
to ±1 px on 10/12 — which pins the mechanism.

**The tick calibration here is sound** (90.00 px/decade, zero residual, validated
against five unlabelled minor ticks), so this is *not* class 20. It is the same
class as F-018 in a different curation: **F-018 and F-020 are one class — a
bar-top detector fooled by replicate markers — now found in two independent
curations on the same day.** That makes it the highest-value digitization check
to automate.

The reviewer confirmed the ALC-0315 3-vs-4-carbon claim by its own centreline
profiling, eye-counted all nine E-structures at 8-16× agreeing with the
geometric trace, and reproduced both printed Fig 5D numbers. **It also endorsed
the LipidDB exclusion** with every supporting statistic re-derived.

---

## Documentation repairs applied 2026-10-30 (no data touched)

These are statements in our own documents that the repository contradicts. Each
was verified false before editing, and none changes a value or a structure.

| what | was | now |
|---|---|---|
| `CHANGELOG.md` §11, `CLAUDE.md` §9, `ERROR_CLASSES.md` | "the LipidDB publication count was **299**, not 47" | **My error.** I propagated a reviewer's claim into the changelog without checking the curation folder, whose own `06_lipiddb_qc.json` records **56 / 50** against the paper's 47. Nine defensible de-duplication rules give 50-116, so no single count is asserted; all three figures are recorded with the rule that produced them. `CLAUDE.md` §6 says *treat reviewer findings as claims* — I did not apply it to that one. |
| `ql_2026/PROVENANCE.md` line 7 | "NOT merged into `data/`" | merged 2026-10-30; the whole "Before merging" section rewritten as dispositions rather than deleted |
| `ql_2026/PROVENANCE.md` | "`Publication_PMID` is blank (not yet indexed)" | wrong premise — 41145968 was indexed; now populated |
| `ql_2026/PROVENANCE.md` | amine 114 "cannot be resolved from this paper" | **closed** — pinned three ways, incl. the `400-L6-C6` mass control that was already in the folder's own mass check and not recognised as answering it |
| `ql_2026/PROVENANCE.md` | "Decide whether to keep `T3A-L8-C14`" — never decided | decision recorded: keep, carrying the heatmap value, with the reasoning |
| `ql_2026/build_rows.py` + `12_QL_2026_qc.json` | `not_merged: "deliberate - awaiting review"` | `merged: ...`; fixed at the source so a rebuild cannot resurrect it |
| `ql_2026/digitize_fig3c.py` docstring | described the **pre-fix** method ("label centroids … 232 px per decade") while the code uses ticks at 151.0 | rewritten, with the class-20/21 history and the open F-018 caveat recorded in the file |

Regenerated and re-verified after these edits: `build_rows.py` → 344 rows,
`validate.py` 21/21 pass, corpus rebuild 28,885 × 72 / 42 experiments, audit C10
and C11 clean. **No `Experiment_value` changed.**

---

# Wave 3 — adjudication and repairs, 2026-10-30

Corpus after: `data_vitro.csv` **28,885 × 79**, 42 experiments; `data_vivo.csv`
**1,780 × 79**, 19 experiments. Row counts unchanged; schema 72 → 79.
`audit.py` **exits 0** with 62 gated metrics at their accepted values.

Each applied repair went to a separate validation agent first, per
`REVIEW_PROGRAM.md` §5. **Two validators changed the plan materially**, which is
the entire case for that step — both changes are recorded against their findings
below rather than folded in silently.

## Disposition

| id | was | now | note |
|---|---|---|---|
| F-001 | `confirmed` / `open` | **`fixed`** | schema 72 → 79; all **1,152** rows flagged — see the count correction below |
| F-002 | `fixed` | `fixed` | — |
| F-003 | `confirmed` | **`fixed`** | 14 z blocks → 7; selectivity check added |
| F-004 | `confirmed` | `open` | still needs the Suppl Table 5 anchor confirmed |
| F-005 | `fixed` | `fixed` | — |
| F-006 | `confirmed` | **`fixed`** | and re-anchoring it **exposed a systematic ~1% top-edge bias the fit had been hiding** |
| F-007 | `confirmed` | **`fixed`** | `min_area` 800 → 120; a miss records `None` and is asserted on |
| F-008 | `claimed` | `claimed` | not reproduced |
| F-009 | `claimed` | `open` | **do not repair in either direction** — the deliverable is a flag, and the schema can now carry one |
| F-010 | `confirmed` | **`fixed`** | arm decomposition populated; C2 coverage 235 → **634** series, 1 experiment → 11 |
| F-011 | `confirmed` | **`fixed`** | mechanism was NOT what I diagnosed — see below |
| F-012 | `confirmed` | **`fixed`** | baseline + seeded-mutation suite |
| F-013 | `confirmed` | **`fixed`** | 537 rows; validator found a second copy worth 15 more |
| F-014 | `confirmed` | **`fixed`** | and the lookup now refuses to guess |
| F-015 | `claimed` | **`fixed`** | validator: the obvious fix would not have worked |
| F-016 | `confirmed` | **`fixed`** as far as it can be | `WT_2026` repaired; `YX_2024`/`LZ_2025` unrecoverable, now declared |
| F-017 | `confirmed` | **`fixed`**, and **my `confirmed` was wrong** | a class-21 defect, not class-20: the code was right, its comment was wrong |
| F-018 | `claimed` | **`fixed`** | 58/60 bars validated to ≤2 px; offset 1.1069 → **1.0071 ± 0.0247**, i.e. not distinguishable from 1.0 |
| F-019 | `claimed` | **`fixed`** | root cause was a **second class-20 instance, on the colorbar**; floor 1310.3 → 948.6 RLU |
| F-020 | `claimed` | **`fixed`** | re-measured independently and reproduced to the pixel |
| **F-021** | **new** | `open` | found by C14 on its first run — see below |

## Corrections to this ledger's own claims

Two things I wrote above are wrong, and both were wrong in the direction of
making the corpus look better than it is.

**F-016 said "Only these three. Every other multi-cell-type experiment is pooled
correctly."** There are no others. The corpus contains exactly **three**
experiments with more than one `Model_type` — `YX_2024`, `WT_2026`, `LZ_2025` —
and **all three** were per-cell-type standardised. 3 of 3, not 3 of N. The
practical consequence was immediate: the new C12 check had **no negative control
in the corpus**, so "zero false positives on the current corpus" was a vacuous
statement. Its specificity had to be established synthetically instead, over
seven constructed cases including the worst one — a correctly pooled experiment
whose true between-cell offset is 0 *and* whose block sds are equal. It passes
because sample block means are O(1/√n) rather than exactly 0, which is the
discriminator the 1e-3 threshold rests on. Repairing `WT_2026` has now given the
corpus a real negative control.

**F-011 said the signature was "a join that overwrote the column."** That is
true of `06_merged.csv`, where `confidence` is a strict function of the O-tail
(O1 → medium, O2 → high, for all 36, independent of the lipid) — but
`06_merged.csv` is an intermediate that only `structure_check.py` reads, and only
for `lipid_name → IL_SMILES`. The path that would have reached the corpus is a
different defect: `build_new_curation.py` read `lipid_name → IL_SMILES` and then
wrote the **literal** `"high"` for every lipid it found, discarding the
`confidence` and `quaternization_site_assumed` columns sitting beside it. I
diagnosed a mechanism from a pattern in the wrong file and marked it `confirmed`.
The repair is in the builder, not in a join suffix.

## F-013 — what the validator added

The right half was settled against the paper, not assumed: Fig 1b (p.3) draws
B1/B2/B3 with the C3–O bond as a bold wedge **cofacial with the C18/C19 angular
methyls** (a flip-invariant relation, so it survives any figure orientation);
Fig 2a shows the same relation in a different orientation; ¹H C3-H is
4.71–4.58 (the broad axial *tt* of a 3β-acyloxy sterol — a 3α ester would be
near 5.0 and roughly a third as wide); and the Fig S7 ¹³C spectrum, whose typed
peak list **omits C3 entirely**, has an isolated peak at 73.9 ± 0.3 ppm against
~70.9 for a 3α ester. An exhaustive search of both PDFs for
`epi|allo|3-alpha|3α|3β|3-beta|cholestan|CAS` returns **zero** hits in the SI and
only the English word "allow" in the main text. So there is no epimer arm to
destroy.

**The addition that mattered:** `build_new_curation.py` carried its own copy of
the same `RemoveAtom` construction, affecting **15 further structures**
(`MC3-CB*`, `SM-102-CB*`, `OF-2-CB*`, `TT3-CB*`, `Lipid 14-CB*`). Fixing
`enumerate.py` alone would have left 15 epicholesterol structures in the corpus
*after* the fix was announced. The second copy is now deleted and the file
imports the single implementation.

Two smaller catches, both taken: keep the `if oh is None` guard the proposed
snippet had dropped, and do not bond to "atom index 0 of the fragment" — it is
the chain terminus for all three linkers today but is an unstated invariant of
exactly the kind that produced the original bug. The attachment atom now carries
an atom map.

Also established: **`reagent_registry.py` could never have caught this.** It
compares the InChIKey **skeleton** block only, deliberately ("so stereo
descriptors do not cause noise"), and it inspects the `CHL_*` columns — empty on
all 537 rows, because this paper's claim is cholesterol-free LNPs. The sterol
lives inside `IL_SMILES`. So C14 is corpus-wide and compares full keys.

## F-015 — why the obvious fix was wrong

Porting `enumerate.py`'s repaired rule verbatim would **not** have fixed
`Lipid 14`. `beta_ester_arms` walked only `N-C-C-C(=O)`, so it could not see that
molecule's `N-C-C-O-C(=O)` β-**acyloxy** arm; both candidate nitrogens scored 0,
the rule fell through to the atom-index tiebreak, and rewriting the same
molecule tail-first made the two rules **swap answers** — the SMILES-ordering
artifact the rule was written to eliminate, reintroduced by its own repair.

β-acyloxy is now counted alongside β-ester, with a steric tiebreak below it. The
choice is invariant across four spellings of each parent, and C5/C6 still resolve
to the piperazine N drawn in Fig 2a. The site itself stays **unverifiable** — the
paper draws no `-C` conjugate for any of the five commercial cores — so all 15
carry `structure_assumed = 1` with that reason on the row.

The validator also found the docstring's second sentence false independently of
the rule: it said the function was "used here only for MC3-C and SM-102-C" while
the code routes `TT3-C`, `Lipid 14-C` and `OF-2-C` through it too. It explicitly
denied covering the molecule the finding was about.

## F-014 — the cause is worse than a stale lookup

`load_corpus_components` did `comp.setdefault(name, smiles)` — **first row in
file order wins** — over a corpus written in alphabetical experiment order.
**292 names in that table map to more than one structure.** `DLin-MC3-DMA` was
*ambiguous* rather than simply wrong: 2,190 rows of the correct acyclic ester
plus 203 dioxolane rows from `AC_2025`, and `AC_2025` sorts first. So the
structure a name resolved to was a property of which experiment happened to sort
first — visible to no caller and recorded nowhere.

Regenerating after the `AC_2025` name repair fixed these three rows, which is
how they came out correct during this round's work; it fixed them **by luck of
row order**. The lookup now returns its ambiguity set and `smi()` raises rather
than choosing. Every name this curation uses is unambiguous today.

## F-021 — NEW: `ZH_2023` carries a 19-nor sterol, not cholesterol

| | |
|---|---|
| **found by** | **C14, on its first run** — in an experiment no reviewer was assigned this round |
| **status** | `confirmed` (the measurement), `open` (the repair) |
| **severity** | **structure**, 2 vivo rows |
| **confidence** | certain on the measurement; the correct structure needs the source paper |

`ZH_2023`'s `A4I18-1R2C70` (2 rows in `data_vivo.csv`) carries a steroid
succinate half-ester. Cleaving it liberates:

```
C26H44O   DZKYBYMQFBHTTM-UHFFFAOYSA-N    <- what is stored
C27H46O   HVYWMOMLDIMFJA-DPAQBDIFSA-N    <- cholesterol
```

Different formula, different skeleton block. The stored sterol is **missing the
C19 angular methyl** and carries **no stereochemistry at all** (`...CC=C4CC(O)CCC4C3CCC12C`
— not one `@` in the steroid). It is a 19-nor sterol, so this is not a
stereochemistry problem like F-013 but a connectivity one.

Gated in `audit_baseline.json` at its current count so it cannot spread, and
explicitly **not** accepted as correct. Resolving it needs the `ZH_2023` source
paper, which is outside this round's five targets.

**Why this one matters out of proportion to 2 rows.** It is the first finding in
this project produced by a check rather than by a person reading a figure, and it
landed in the first experiment outside the reviewed set. The decision to make
C14 corpus-wide instead of scoping it to `ZT_2026` is what found it — scoping a
new check to the experiment where its bug was found is how the same bug survives
into the next curation.

## Still open, with what each needs

| id | what it needs |
|---|---|
| F-004 | independent confirmation of the Suppl Table 5 anchor before either repair; the honest alternative is to blank `N:P` on 129 rows rather than reconstruct it |
| F-006/F-007 | re-anchor `digitize_fig2fgh.py` on tick spacing, add asserts, and correct the claim in three places (`11_LL_2026_qc.json`, `lit_review/build_data.py`, `CLAUDE.md` §5). The 48 values themselves were re-measured against a tick-anchored fit and are within 0.75 pp, so no re-digitization |
| F-008 | reproduce the 57 px window overrun; all 147 bins are unchanged either way, so this is latent |
| F-009 | a flag, not a structure change. The schema can now carry it |
| ~~F-010~~ | **done** — per-series test against the measured chain length; 634 checked, 338 correctly skipped |
| ~~F-018~~ | **done** — and my stated reason for reverting the first attempt was itself a measurement error; see below |
| F-021 | the `ZH_2023` source paper |
| 7 C10 anchors | byte-identical to their enumeration with no generator script: 6 in `ng_2025` (another session's in-flight work, not inspected) and `lin_aromatic_2026`'s `10A2-3B-EP10` |

## Count correction: F-001 is 1,152 rows, not 1,137

Measured against the emitted corpus after the repair:

| experiment | rows | what the flag says |
|---|---|---|
| `LL_2026` | 840 | `positional_ambiguity = 1`; 140 of those also `confidence = medium` (commercial branched isomer mixtures) |
| `LL_2026_2` | 147 | `positional_ambiguity = 1` |
| `QL_2026` | 117 | 113 `value_is_bound = upper` (colormap floor), 3 `lower` (censored), 1 cross-figure conflict |
| `ZT_2026` | 39 | 24 `confidence = low` + `structure_assumed`, **plus 15** commercial-core conjugates whose conjugation site is inferred |
| `XD_2026` | 9 | `confidence = medium`, `structure_assumed = 1` — no external anchor of any kind |
| | **1,152** | |

This number has now moved twice, both times upward, and both corrections are
mine. My first figure of ≥266 counted only flags enumerable from QC JSONs and
missed `LL_2026`'s 987 `positional_ambiguity` rows. I then adopted CORPUS-B's
1,137 and wrote that it was "the better number" — it is better than ≥266, but it
undercounts `ZT_2026` by 15, because it took the 24 lipids with an assumed
quaternization site and not the 15 commercial-core conjugates whose cholesterol
conjugation site is equally inferred (the paper draws no `-C` conjugate for any
of the five cores; that is exactly why they are `medium` rather than `high`).

The pattern is worth naming, because it is the second time this round: I accepted
a reviewer's count as an improvement on my own without re-deriving it, the same
way I propagated "299 publication identifiers" into the changelog. An
improvement on a wrong number is not a right number. The figure above was
measured from `data_vitro.csv` and `data_vivo.csv` directly and can be
re-measured with `audit.py` C13's `rows_carrying_a_flag`.

## Report files on disk

| reviewer | file | provenance |
|---|---|---|
| ZT_2026 | `report_ZT_2026.md` | recovered from the session transcript |
| LL_2026 | `report_LL_2026.md` | recovered from the session transcript |
| WT_2026 | `report_WT_2026.md` | written by the reviewer |
| XD_2026 | `report_XD_2026.md` | written by the reviewer |
| QL_2026 | `report_QL_2026.md` | recovered from the session transcript |
| CORPUS-B | `report_CORPUS_B_raw.txt` | preserved inline, unformatted |
| CORPUS-C | `report_CORPUS_C_raw.txt` | preserved inline, unformatted |
| **CORPUS-A** | **missing** | **not recovered — see below** |
| validation: F-013 | `validate_F013.md` | written by the agent |
| validation: F-014/F-015 | `validate_F014_F015.md` | written by the agent |

Five of the eight Wave-1/2 reviewers could not write their own report file — the
harness blocked the write — so they returned the text inline, and three of those
five were then lost when this session's context was compacted. The three were
recovered by searching the session transcript for the hand-back messages.

**CORPUS-A's report was not recovered.** A first attempt matched a tool-result
fragment rather than the report and would have left a mislabelled file on disk
asserting a provenance it did not have; it was deleted rather than kept. Its
findings survive here (F-016's origin, and the dismantling of the `YX_2024`
pooling argument) because I reproduced them independently before acting, but the
report itself is gone. Two process fixes follow, and the second is the one that
matters:

1. Give every reviewer a writable report path and **verify the file exists**
   before accepting the hand-back as complete. A report that exists only in a
   conversation is one compaction away from gone.
2. Do not let a recovery invent provenance. The deleted file would have carried
   a comment saying it was CORPUS-A's report, which would have been false, and
   nothing downstream could have told. This is the same failure as the C10 anchor
   whose note claimed independent assembly (F-013) — a provenance claim that no
   code checks is worth nothing, and writing one is worse than writing nothing.

## F-006 — the fit was not just blind, it was compensating

The finding as logged said the fitted scale meant a uniform multiplicative error
in the heatmap read would be "absorbed and invisible". That is right, and it
understates the problem.

Re-anchoring the scale on tick centres turned **all 48 deltas negative**, −0.4 to
−1.6 pp. A sign that consistent over 48 independent comparisons is a systematic
error, and the cause is the colour mask: it selects only fully saturated bar
pixels, so every top edge read a few pixels low — about 1% on a 400 px bar. The
fitted scale was 5.046 px/pp against a tick-implied 5.100, a **1.07%**
difference, cancelling a **~1.1%** height shortfall.

So the fit was deriving a correction factor **from the answer** and applying it
to hide a real measurement error. The 48-point agreement was not weak evidence;
it was evidence of the fit's own convergence. That is a sharper statement of
class 27 than the one in the registry, and the registry row now says so.

Fixed at the cause: sub-pixel top edges at the half-intensity point of the
anti-aliasing ramp, uniformly for every bar. Max |delta| 1.60 → **1.00 pp**,
mean −0.87 → **−0.52 pp**.

**The remaining −0.52 pp is left in place, deliberately.** It is still one-sided
(47 of 48 negative) and its probable cause is identified — the bottom edge, since
a bar sits on the axis rule and the mask's last saturated row is 1–2 px above the
rule's top, worth about 0.28 pp on panel 2h. It is not chased because these 48
numbers are a *validation* of the Fig 2e heatmap read and no corpus value is
computed from them, so the residual costs accuracy in a check rather than in the
data. Recording a known one-sided residual is the opposite of what the old code
did with the same information.

## F-020 — reproduced to the pixel, with three corrections to the finding

All three corrections strengthen it:

1. The 60% rule reproduces the stored tops on **12 of 12 at exactly 0 px**, not
   "10 of 12 to ±1 px". There are no exceptions needing explanation — the rule
   *is* the error. Confirmed by the fill fraction at the true top row: **0.21–0.57**
   on the nine affected bars, **0.68–0.97** on the three unaffected. The 0.60
   threshold is precisely the discriminator.
2. "90.00 px/decade" and "0.5–10.9 px" are **native 300-dpi** figures. The page
   embeds a single 1541×1864 PNG that the 600-dpi render upsamples 2×, so at the
   render it is **179.96 px/decade** and the gaps are **0–21 px**. Quoting the
   offsets without the unit makes the finding read 2× off.
3. Mean bias is **+0.051** decades, not +0.060; max **+0.117** (E4), not +0.122.
   The z-rms 0.034 and "all 12 rows change" were exact.

Reproduced here independently: ALC-0315 +13, SM-102 +7, MC3 +7, E1 +18, E2 +2,
E3 +9, E4 +21, E5 +15, E6 +18, E7/E8/E9 +0 — identical to the re-measurement.

**The calibration is sound**, re-verified against 31 unlabelled minor ticks (not
five) at rms 0.53 px, with major-tick residuals ≤4e-4 decades. So F-020 is not
class 20 and the two errors do not compound; there is one error.

**The limit this establishes on the printed-statistic rule.** All 12 lipids still
fall in their stated Fig 5D band and both printed Fig 5D numbers still reproduce
— on the wrong values as well as the right ones. The project's rule
*"reproduce every printed summary statistic"* is necessary and not sufficient,
and geometry has to be checked geometrically. Correcting E3 narrows its margin to
the B/A boundary from 0.132 to **0.082 decades (≈15 px)**, the thinnest in the
panel, so that cross-check was closer to failing than anyone knew.

One trap for anyone revisiting this: a global `min` over all bar columns is not a
safe replacement for the 60% rule — on E1 the anti-aliased halo of a replicate
disc floating above the bar makes columns 23–25 report a top 8 px too high. The
flank rule reads the 5 leftmost and 5 rightmost columns and requires them to
agree within 1 px.

## F-017 — I confirmed this one wrong, and the error is instructive

I marked F-017 `confirmed` on this reasoning: `digitize_fig2c.py`'s comment gives
y-axis anchors of 2342 / 2748 / 3156, whose decade gaps are **406 and 408**; a log
axis has exactly equal decade spacing; therefore those are glyph centroids and
the code is anchored on text.

The arithmetic is right and the conclusion is wrong. Those three integers are
neither the label centroids (measured 2344.3 / 2750.8 / 3158.9) nor the tick
positions — **they are the tick positions rounded to integers**, and the rounding
is the entire apparent inequality. Measured centroids are 2342.78 / 2749.54 /
3156.75, a linear fit gives 406.99 px/decade at ±0.037 px residual, so the
stored `PX_PER_DECADE = 407.0` was right to 0.01 px and only `Y_1E4` moves 0.75
px — rms Δz **0.0006** across the slice.

So F-017 is class 21 (provenance text contradicting the code) and not class 20
(wrong anchor). The sibling repair in `digitize_fig3c.py` was re-derived as a
control and is **confirmed exact**.

**What I did wrong, specifically.** I verified a claim *about the code* by
reading a *comment*, and a comment is not the code. That is the same class of
error the comment itself was an instance of — and this round had already logged
F-011 as a mechanism inferred from the wrong file, and the "299 publications"
and "1,137 rows" figures as reviewer numbers adopted without re-derivation. Four
instances now, all of the same shape: a conclusion drawn from a secondary
artifact rather than the primary one.

The decisive evidence here was not available from any of our files. It is a round
number the paper never printed: Fig 2c's red dotted reference line sits at
y = 2872.21, which the corrected constants put at **50,020 RLU, 0.04% off a round
5e4** — against 49,804 under the stored constants and 49,159 under an earlier
reviewer's proposed 406.0 / y=3153.0, which is thereby refuted. Finding an
unprinted round number the authors drew and checking against *that* is a
genuinely external test, and worth adding to the playbook alongside
"reproduce every printed statistic".

## F-019 — the root cause the claim did not name

F-019 was logged as a colormap-floor degeneracy: identical RGB over a wide band,
so floored cells are bounds rather than measurements. True, and it has a cause
nobody had named: **`BAR_BOT = 1973` is the outer edge of the colorbar's black
border, not the ramp's zero** — a second instance of error class 20, on a
colorbar instead of an axis.

Verified from the pixel profile directly: the border runs y = 1962-1973. Mapping
0-6e4 across 870 px instead of ~860 biases **every** Fig 2b value high, by about
+1% at the 2e4 hit threshold and **+38% at the floor**.

Before applying a change to all 320 library-1 values I cross-checked the anchor
by a second route. The bar's bottom border *centre* is 1967.5 against the
tick-derived 1967.23 — 0.27 px — and a border-centre anchor pair gives a
white-point value of 903.8 RLU against 913.7 from the ticks. **The two corrected
routes differ by 9.9 RLU; the old anchor was +362 RLU out.** So the correction is
robust and the residual uncertainty is ~3% of it. The top anchor is the weaker
half and is documented as such in the code: the top border abuts a near-black
ramp, so its inner edge is not readable from colour, which is why the ticks are
the authority there.

**The part that matters most is not the accuracy, it is the contradiction.** The
48 white cells have a feasible band of [0, ~982] RLU and were stored at 1310.3 —
**above the entire band**, 1.43× its upper bound. This round had just flagged
those same rows `value_is_bound = upper`. The flag and the number disagreed, and
the flag was right. The stored value is now 948.6, at the top of the band.

Recorded on the rows as well: the limit is the colormap's own quantisation, not
image noise (adjacent LUT steps near the floor differ by 1.7-3.5 RGB units while
each spans 349-1,045 RLU, against within-cell pixel sd of 0.06-0.34), and **56%
of the slice carries one of only six values, all from the bottom six LUT steps**.
That is a data-quality fact with direct bearing on Result 6, and nothing recorded
it before.

The paper's own between-library conclusion survives every scenario tested and
strengthens under the correction, so this is a corpus-slice defect and not a
challenge to Liu et al.

## F-018 — confirmed, deliberately not applied

The mechanism and the count reproduce: 15 of 60 bars, up to +28%, and the stored
cross-library offset of 1.1116 is an artifact of dot contamination on the
benchmarks themselves. The honest reading is that there is **no measurable
offset** between the two screens (1.0070 ± 0.0247).

I did not apply it. My attempt at a replacement — locate the two side-stroke
centres below the top, then require one contiguous run spanning ≥72% of the width
and containing both — failed on the first bar, because Fig 2c's bars are 12 px
wide and **solid**, not hollow. Two other approaches are already excluded: a
stricter coverage threshold (merged dots reach 91-97% contiguous coverage) and
pure flank reading (wrong on anti-aliased outer stroke columns).

Leaving a known-biased rule in place is bad; shipping an unvalidated detector
over 60 bars is worse, because a characterised bias stops being characterised the
moment it is replaced by an unknown one. The rule stays, with a comment stating
the bias, its magnitude, what has been ruled out, and where the per-bar
measurements are. Returned to the agent that has the working detector, with a
request for the patch plus a per-bar verification table covering all 60.

**One argument retired.** F-018 offered two pairs of byte-identical bar values as
a bug signature. Both pairs are real and both members are spurious reads — but
exact ties are not themselves evidence: at 151 px/decade on a 2 px row pitch Fig
3c offers only ~75 distinguishable values per decade, and the corrected set still
contains one exact tie. Credit to the re-measurement for volunteering the
counter-example to its own finding.

## F-018 applied — and the fifth instance of one pattern

Applied from a validated patch. Post-apply check as specified: the offset reads
**1.0071, sd 0.0247**, matching the tested patch exactly.

**I reverted the first attempt for a reason that was wrong.** I said Fig 2c's
bars were "12 px wide and solid". They are 70-72 px wide and hollow; 12 px is the
*left side stroke* of `110-L1-C10`. The cause: I reimplemented the `coloured`
mask instead of reproducing the script's. The script computes
`sat = img.max(2) - img.min(2)`; I used a normalised saturation, which fragmented
the mask so only the stroke registered. With the script's own mask the
segmentation gives 32 bars at 70-72 px with exactly two stroke groups and 69%
white interior, which I have now confirmed directly.

**That is the fifth instance in this round of the same mistake**, and the ledger
should name it once rather than five times:

| | what I concluded from | what I should have used |
|---|---|---|
| "299 publication identifiers" | a reviewer's claim | the curation folder's own QC file |
| F-001 = 1,137 rows | a reviewer's count | the emitted corpus |
| F-011's mechanism | a pattern in `06_merged.csv`, which nothing merges | the builder that writes the rows |
| F-017 "confirmed" | three rounded integers in a **comment** | the tick marks the code uses |
| F-018 "bars are solid" | **my own** reimplementation of the mask | the script's mask |

Every one is a conclusion drawn from a secondary artifact standing in for the
primary one, and every one made the corpus or the finding look different from
what it was. Four of the five were caught by someone re-deriving them. The fix is
not "be more careful" — it is a rule: **when checking a claim about code, run the
code's own objects; when checking a claim about data, read the data.** Added to
`REVIEW_PROGRAM.md` §5.

Reverting was still the right call, because neither of us had validated the rule
bar by bar at that point. But the first attempt was closer than I thought: it
matched an independent measurement to ≤2 px on **56 of 60** bars and needed two
specific fixes rather than a new model —

- **a luminance mask, not a blue one.** Fig 2c's red dotted 5e4 reference line
  crosses the top rule of `DLin-MC3-DMA` and `SM102` and punches through it on a
  blue-only mask; `DLin-MC3-DMA` returned no top at all, which is the failure I
  saw and misattributed to the geometry.
- **strictly contiguous runs (gap = 0).** Bridging a 2 px gap lets a merged dot
  cluster span both stroke centres, which it does on `T3A-L2-2C10` and
  `ALC-0315`.

**A new trap the fix creates, verified before applying.** `digitize_fig3c.py`
searched from `int(Y_1E3 - 4*PX_PER_DECADE) + 4 = 2311`, and the panel box's top
rule occupies **y = 2301-2314** — I confirmed 14 consecutive rows above 50% ink
across the panel. Harmless on a blue mask because the box is black; under a
luminance mask **every Fig 3c bar would have returned 2311**. Changing a mask
from hue to luminance makes all the black furniture visible, so the search window
has to be re-derived rather than inherited. Both files now set
`TOP_SEARCH = 2320` and assert the first scanned row is not already a full-width
rule. Worth generalising: **a mask change invalidates every window that was
tuned against the old mask.**

**Validation:** 58 of 60 bars within 2 px of an independent hand measurement; 2
at 4 px, both from the `INK_THR` plateau sweep. `T3A-L2-C8` is 1.0% in value and
immaterial; **`T3A-L2-2C10` is 6.2%** because Fig 3c has only 151 px/decade, and
its honest value is reported as y = 2521 ± 2 (374k-449k RLU) rather than a point
estimate. The previous value was wrong by 14-18 px either way. Rule thicknesses
of 6-14 px (2c) and 2-10 px (3c) are consistent with a stroked outline, which is
itself evidence the detector is finding rules and not dots.

**What it changes, and the part that is not a refinement.** The stored
cross-library offset was never a property of the two screens: ALC-0315 was
dot-contaminated in *both* panels (+24 px in 2c, +16 px in 3c) and Lipid 5 in 3c
(+8 px) — three of the four benchmarks the offset is measured from. Correcting the
detector does not sharpen the offset, it **removes** it: 1.1069 → **1.0071 ±
0.0247**, i.e. the two screens agree to within the spread of the benchmarks
themselves and the offset is not distinguishable from 1.0. The division is kept
because 1.0071 is the honest estimate and writing 1.0 would assert an agreement
the figures do not quite show.

Per-library mean z −0.1566/+2.0885 → −0.1583/+2.1106; geometric-mean ratio
43.3× → 40.8×; rms Δz 0.0188, max |Δz| 0.1419, 31 rows >0.02 z.

**F-019 is still the largest of the three** — it alone moved library 1's
geometric mean from 4,128 to 3,569 RLU and the ratio from 37.5× to 43.3×.

**One recommendation was withdrawn by its own author and not applied:** that
`digitize_fig2b.py`'s LUT window ran 6 px into the border skirt. It ends at 1955,
exactly where the ramp ends. Recording the withdrawal because an un-retracted
recommendation in a 68 KB report is the kind of thing a later session applies
without checking.
