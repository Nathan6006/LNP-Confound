# Round 2 — findings ledger

Opened 2026-11-07. Procedure in `../REVIEW_PROGRAM.md`. Round 1's ledger is
`../round_2026-10-30/FINDINGS.md`; numbering continues from it, so this round
starts at **F-022**.

**This round started from round 1's open list, not a blank sheet** — that was
the explicit lesson written into `REVIEW_PROGRAM.md` §7 after round 1, because
an open finding nobody re-reads becomes closed by attrition. Carried in:

| carried from round 1 | assigned to |
|---|---|
| F-004 — `WT_2026` `N:P` wrong on 129 rows, repair direction unsettled | the `WT_2026` reviewer, as its priority |
| F-008 — `LL_2026` Fig 2a sampling window 57 px overrun, unreproduced | the `LL_2026` reviewer |
| F-009 — enantiopure `5B` stereocentre unassigned on 141 rows | the `LL_2026` reviewer |
| F-021 — `ZH_2023` carries a 19-nor sterol | **blocked**: needs the Angew paper (He et al. 2023), which is not held |
| 7 C10 anchors byte-identical to their enumeration | 1 (`10A2-3B-EP10`) to the `LL_2026` reviewer; 6 in `ng_2025`, not inspected |

**Wave 1 agents:** five, one per curated paper, read-only, briefed from packets
generated out of live data by `brief.py` on the day. Each was told what round 1
already fixed in its curation and instructed **not** to spend budget
re-confirming it, which is the difference between a round and a re-run.

---

## F-022 — 915 `IL_name` cells carry a trailing space, and no check could see it

| | |
|---|---|
| **found by** | me, incidentally, while pulling `ZH_2023` rows for F-021 |
| **status** | **applied** 2026-11-07 |
| **severity** | silent **join hazard**; no value, structure or row count changes |
| **confidence** | certain — measured on the emitted corpus both before and after |

`"A1I5R22C8 "` and `"A1I5R22C8"` are different dict keys carrying the **same**
`IL_SMILES`. Any `groupby`, join or set operation keyed on `IL_name` splits one
lipid into two.

**Extent**, measured from `data_vitro.csv` and `data_vivo.csv` directly:

| | cells | where |
|---|---|---|
| vitro | 720 | all `BL_2023` |
| vivo | 195 | 150 `ZH_2023`, 45 `BL_2023` |
| **total** | **915** | `IL_name` only — **no other column in the corpus carries padding** |

All of it is **trailing**, never leading. In `ZH_2023` the padding tracks the
organ *exactly* — every liver row padded, every spleen row not — which is the
fingerprint of two source tables concatenated upstream, and is why the pairs
look like duplicates without being any.

**Why it survived.** It was *known*. `CLAUDE.md` §3 listed it under "Things that
will bite you" as "Several `BL_2023` `IL_name` values carry a trailing space" —
understating it by two orders of magnitude — and `build_dataset.py` called
`.strip()` defensively at two comparison sites rather than fixing it. Meanwhile
C7's required-field check **also** called `.strip()` before asking whether a
value was blank, so a padded value passed as present and the check was
structurally incapable of seeing the class. That is **error class 37**: a check
that canonicalises its input before comparing has given up the ability to detect
non-canonical input.

A known-and-worked-around defect is worse than an unknown one. The workarounds
are evidence that several passes saw it, each did the locally cheap thing, and
nobody owned it.

**Verified safe before applying**, because stripping identifiers can merge rows:

| test | result |
|---|---|
| distinct names collapsing into one inside any experiment | **0** |
| stripped name mapping onto a second `IL_SMILES` | **0** |
| duplicate `(Experiment_ID, IL_name, Model_target, Model_type)` keys | **40 → 40**, unchanged |
| the 12 `ZH_2023` names existing in both forms | one lipid per **two organs** (liver + spleen) — two legitimate rows |

**Applied.** `strip_cells()` in `build_dataset.py`, called from `main()` before
`write()` *and* before `validate()`. The first attempt stripped inside `write()`
and was wrong in an instructive way: the emitted file was clean while
`validate()` went on reading the padded rows still in memory, so
`dataset_validation.json` named `A55-T3 ` as a duplicate-SMILES partner for a
lipid the corpus no longer contained. **A fix applied downstream of the
validator makes the validator lie** — the same docs-disagree-with-data split
this project keeps hitting, this time created by the repair itself.

**Post-apply diff against the pre-fix corpus:** 720 vitro + 195 vivo rows
changed, exactly 915 cells, `IL_name` the only column touched, every change
exactly a strip, row order and `LNP_ID` alignment preserved, row counts
unchanged at 28,885 / 1,780.

**Now gated.** `C7_fields.padded_cells` checks **every** cell in both corpora,
baselined at 0 with its reason; the baseline is 62 → **64** metrics. Two seeded
mutations added to `test_checks.py` (a trailing space on `IL_name`, a leading
space on `Model_type`), so the gate is proven to fire rather than assumed to.

---

## F-023 — `build_manifest.py` could not run: it read its input from a deleted session scratchpad

| | |
|---|---|
| **found by** | me, running it as part of this round's regeneration |
| **status** | **applied** 2026-11-07 |
| **severity** | **reproducibility** — a tracked output with no recoverable lineage |
| **confidence** | certain — it raised `FileNotFoundError` |

`build_manifest.py` hard-coded

```
ESUM = ".../claude-501/-Users-nathanliu-Documents-LNP-confound/<uuid>/scratchpad/esum.json"
```

an NCBI esummary cache sitting in a **Claude session scratchpad**. Scratchpads
are session-scoped and get cleared. By today the file was gone, the script
raised, and `data/experiment_manifest.csv` — **tracked in git** — could not be
regenerated from anything in the repository. The manifest had quietly become an
artifact with no lineage, which is the one thing `provenance/` exists to prevent.

**Second instance of the class.** `lit_review/build_data.py` was unrunnable for
exactly this reason on 2026-10-30 and its sweep inputs were recovered into
`lit_review/provenance/`. One occurrence is an oversight; two is a missing rule,
so the rule is now written down: **a script in this repository may not read an
input from a session scratchpad.** Scratch is for intermediates nobody needs
twice; anything a tracked output depends on belongs in the repo.

**Applied.** New `provenance/fetch_esummary.py` derives the PMID list from the
emitted corpus, fetches the records from NCBI E-utilities in batches (generic
`tool` name, no email — the user's address is not ours to hand to a third-party
service), and writes `provenance/esummary_cache.json` (49 records, 163 kB).
`build_manifest.py` now reads that, and its paths are derived from `__file__`
instead of hard-coded under `/Users/nathanliu/`.

**The recovery is faithful, and it closed two gaps.** The rebuilt manifest is
byte-identical to the 2026-10-30 version on 59 of 61 rows. The two that changed
both *improved*: `QL_2026` and `XD_2026` previously read **"PMID … not checked
against NCBI"**, because the stale cache predated those two merges, and now
carry a verified citation with `PMID → DOI` confirmed against esummary. So the
cache had been silently stale as well as ephemeral, and the manifest was
recording two unverified identifiers as a note nobody read.

Cached rather than fetched per run on purpose: the manifest build stays offline
and deterministic, and a rebuild years from now reproduces today's answer
instead of picking up a publisher's later metadata correction. A PMID missing
from the cache degrades to a note, not a crash — a PMID the corpus holds that
PubMed does not know is a finding worth seeing.

---

## F-024 — "the number of flagged rows" has no definition in code, and the pointer we gave returns a different number

| | |
|---|---|
| **found by** | me, Wave-2 style docs-vs-data sweep: every number `CLAUDE.md` asserts, checked against the live corpus |
| **status** | **applied** 2026-11-07 |
| **severity** | **a number quoted in the manuscript's supporting material is wrong**; no data changes |
| **confidence** | certain — both counts now come out of `audit.py`'s own code |

Round 1's ledger states the F-001 figure as **1,152** and tells the reader it
"can be re-measured with `audit.py` C13's `rows_carrying_a_flag`." That metric
returns **1,652**.

Neither is a typo. They answer different questions, and nothing wrote the
question down:

| metric | count | what it means |
|---|---|---|
| `rows_carrying_a_flag` | **1,652** | any quality column populated at all — `confidence = high` included |
| `rows_with_a_caveat` | **1,157** | the structure or value is genuinely qualified |
| round 1's figure | 1,152 | a prose definition, re-implemented by hand |

The 495-row gap between the first two is rows carrying **`confidence = high`
and nothing else**: 486 `ZT_2026` vitro, 6 `ZT_2026` vivo, 3 `XD_2026` vivo
benchmarks. `high` is an *assurance*. Counting it as a flag inverts the meaning
of the number.

**And 1,152 was wrong twice over, in ways that nearly cancelled:**

| experiment | round 1 said | actually | why |
|---|---|---|---|
| `LL_2026` | 840 | 840 | — |
| `LL_2026_2` | 147 | 147 | — |
| `QL_2026` | 117 | **116** | 113 `upper` + 3 `lower` = 116. The 117th was a cross-figure conflict that **carries no flag column at all**, so it was counted in a total of flagged rows while not being one. |
| `ZT_2026` | 39 | **45** | 39 is the `structure_assumed = 1` count. Six more rows carry `confidence = medium` with no other marker and were missed. |
| `XD_2026` | 9 | 9 | — |
| | 1,152 | **1,157** | −1 + 6 |

**This count has now moved four times** — ≥266 → 1,137 → 1,152 → 1,157 — and
every correction came from someone re-deriving it by hand. The ledger already
named the pattern ("an improvement on a wrong number is not a right number")
and it kept happening, because naming a pattern does not stop it. **A quantity
with a prose definition and no implementation gets re-implemented, slightly
differently, by every pass that touches it.**

**Applied.** C13 now emits both counts with the distinction stated in code, and
both are gated — 64 → **70** metrics. `CLAUDE.md` §6 carries the two-row table
and says which to quote. Round 1's ledger is left as written, with the
correction recorded here; rewriting history would lose the fact that the number
was wrong.

**It also earned itself immediately.** A third new metric,
`confidence_without_a_note`, found **6 `ZT_2026` vivo rows** —
`C5O1B1`, `C6O1B1`, `C5O1B2`, `C6O1B2`, `C5O1B3`, `C6O1B3` — carrying
`confidence = medium` with an **empty** `curation_note`. The existing
`flagged_without_a_note` check returned 0 because it only ever inspected
`structure_assumed` and `value_is_bound`, never `confidence`. So the corpus says
six structures are less certain and does not say why, which is error class 34 in
miniature. Gated at 6 **with the instruction to re-baseline to 0 once the
curation supplies the reason** — not to clear the flag instead.

---

## F-025 — 15 rows share one assumption and one note, and two different confidences

| | |
|---|---|
| **found by** | me, following F-024's `ZT_2026` discrepancy into the curation |
| **status** | **applied** 2026-11-07, with the validator's modification — see the verdict below |
| **severity** | **metadata**, 12 vivo rows touched; no structure, value or row-count changes |
| **confidence** | certain on the mechanism; the repair direction was settled by validation |

All 15 `ZT_2026` commercial-core cholesterol conjugates carry `structure_assumed = 1`
and the **identical** `curation_note`:

> "Cholesterol conjugation site inferred: the paper draws no -C conjugate for any commercial core."

Their `confidence` is split:

| rows | `confidence` |
|---|---|
| `Lipid 14-CB1/2/3`, `TT3-CB1/2/3` | **`high`** |
| `MC3-CB1/2/3`, `SM-102-CB1/2/3`, `OF-2-CB1/2/3` | `medium` |

Same evidence, same assumption, same sentence, two answers.

**The cause, traced in `curation/build_new_curation.py`.** `COMMERCIAL_CORES`
stores `(core_SMILES, confidence)` where that confidence describes **how well
the parent core's structure is known** — `TT3` and `Lipid 14` are `"high"`,
`OF-2` is `"medium"`. The resolution branch

```python
elif core.endswith("-C") and core[:-2] in COMMERCIAL_CORES:
    parent_smiles, conf = COMMERCIAL_CORES[core[:-2]]
```

then attaches that to the **conjugate** and sets `structure_assumed = "1"`
unconditionally. The sibling branch (`MC3-C` / `SM-102-C`) hard-codes
`conf = "medium"`, which is why the split falls where it does.

So `confidence` here is answering *"how well do we know the core?"* while
`schema_quality.json` defines it as the confidence of **the structure that was
emitted**. The emitted structure is a conjugate with an inferred attachment
site. Not the same question.

**Not applied.** The obvious repair — cap the branch at `medium`, moving 6 rows
`high → medium` — is out with a validation agent, because the alternative
reading is live: if the paper *does* draw the conjugation site for `TT3` and
`Lipid 14` and not for the others, then the 6 `high` rows are correct and the
note is what is wrong. The validator was also asked whether
`COMMERCIAL_CORES` should stop carrying a confidence at all, since it is being
read as something it is not, and whether the same value is consumed anywhere
else. Round 1's lesson there was that fixing one copy left a second live copy
affecting 15 more structures.

This is the F-011 class — a `confidence` that reports something other than the
emitted structure's — surviving the round that was supposed to have cleared it.

### The validator changed the repair, and was right to

Verdict: **APPLY WITH MODIFICATION.** Full write-up in `validation_F025.md`. I
re-derived every mechanical claim it made before applying any of it; all of them
reproduce exactly.

**What it added that I had missed, and it is the sharper defect.** The local
`quaternize()` wrapper calls the shared `_EN.quaternize()` and assigns its second
return value to `_note` — i.e. **discards the site classification the shared
function already computed**. I ran `enumerate.py::_quaternization_site` on all
five cores myself: every one returns the literal string `"ASSUMED: no ring N and
the paper draws no structure"`, and that exact string is paired with
`confidence = low` on **24 of 24** rows in `03_enumerated.csv`. So the curation's
own implementation already classifies these five into the `low` condition, and
the branch threw that away and wrote a softer label from a table that means
something else.

**Why a flat `medium` was wrong, and why flat `low` would also have been wrong.**
I reproduced its symmetry analysis — tertiary amines distinct under
`CanonicalRankAtoms(breakTies=False)`:

| core | tertiary N | distinct classes | the site is |
|---|---|---|---|
| MC3 | 1 | 1 | forced |
| SM-102 | 1 | 1 | forced |
| TT3 | 3 | **1** | a choice among equivalents — not a choice |
| OF-2 | 2 | **1** | likewise |
| **Lipid 14** | 2 | **2** | **a genuine two-way choice** |

So the labelling was **inverted**: the three rows with the most structural
uncertainty (`Lipid 14-CB*`) carried the project's *highest* confidence, while
six rows whose site is forced by having one amine carried `medium`. And
`Lipid 14` is not hypothetical — `PROVENANCE.md` records it was previously
alkylated at the *hindered central* nitrogen, the answer falling out of a
SMILES-order tiebreak (F-015).

**Three project documents already said `medium`; the code was the outlier.**
`quaternize()`'s own docstring 140 lines above the branch, `PROVENANCE.md` §4,
and `new_curation_qc.json`'s `confidence_basis` — while the *same* qc file's
`by_confidence` lists `Lipid 14-CB*` and `TT3-CB*` under `high`. A file
contradicting itself.

**Applied, as the validator recommended rather than as I proposed.** A cap at
`medium` was rejected: it would leave `"TT3": (..., "high")` in the table for the
next reader to misuse. `COMMERCIAL_CORES` is now **SMILES-only**, and a new
`_site_confidence()` derives both the confidence and its note from the core's
own structure — `low` where more than one chemically distinct amine could carry
the conjugate, `medium` where the site is forced. The `PARENT` branch's
hard-coded `"medium"` now comes from the same function. This matches the pattern
already established elsewhere in the repo (`enumerate.py:309` forces `low` over
the component confidence when the site is ASSUMED; `lin_aromatic_2026` starts at
`high` and downgrades) — the commercial-core branch had been the only place
passing a component confidence through with no downgrade while setting
`structure_assumed = 1`.

**Measured effect**, position-aligned diff against the pre-fix curation:
**21 of 537 rows changed; columns touched `curation_note` (21), `confidence` (6),
`structure_assumed` (6). No SMILES changed.**

| rows | change |
|---|---|
| `Lipid 14-CB1/2/3` | `high` → **`low`** |
| `TT3-CB1/2/3` | `high` → **`medium`** |
| `MC3-CB*`, `SM-102-CB*`, `OF-2-CB*` | unchanged at `medium`, new more specific note |
| `C5O1B*`, `C6O1B*` | `structure_assumed` `''` → **`1`** + a note (see below) |

**Scope check, run myself:** `grep -rn COMMERCIAL_CORES` gives the definition and
**exactly one** read. No second copy — round 1's duplication hazard is closed
because the file now imports `enumerate.py`. The validator also swept every other
curation's confidence literals and found `siegwart/enumerate.py` and
`lin_aromatic/enumerate_library.py` both implement the *correct* weakest-link
downgrade. It flagged one thing it did **not** verify —
`lin_aromatic_2026/enumerate_147.py:111` injecting `{**r, "confidence": "high"}`
at a call site — which is passed to the `LL_2026` reviewer rather than taken on
its word.

### A second finding fell out of the same trace

The 6 rows `C5O1B1`/`C6O1B1`/`C5O1B2`/`C6O1B2`/`C5O1B3`/`C6O1B3` are the ones
F-024's new `confidence_without_a_note` metric had just flagged: `medium`
confidence, empty note. **The reason was in `03_enumerated.csv` all along** —
their `O1` arm is `isodecyl acrylate (branched terminus)`, a commercial branched
isomer mixture sold under one name, so the stored SMILES is one isomer of a
mixture. The `O2` rows, which are `high`, use dodecyl acrylate, a single
compound.

And the qc file's own `confidence_basis` gets this wrong too: it attributes their
`medium` to "the core is C5/C6, whose piperazine N is drawn in Fig 2a" — which
cannot be the discriminator, because the `high` O2 rows have **the same cores**.
The note is now written onto the row and `structure_assumed = 1` set, matching
how `LL_2026` flags the identical condition for iso-octyl / iso-decyl acrylate.

`C13.confidence_without_a_note` went **6 → 0**, which made the audit **FAIL** on
"a known condition VANISHED" — the baseline's decrease-also-fails design catching
a real improvement for the first time. Re-baselined to 0 in the same change, with
the reason, per `CLAUDE.md` §6.

---

## F-027 — `WT_2026`'s PMID was blank on a reason that was false when it was written

| | |
|---|---|
| **found by** | the round-2 `WT_2026` reviewer (its E2) |
| **status** | **applied** 2026-11-07 |
| **severity** | **metadata**, 2,544 rows — plus a false `why` in the audit baseline, which is worse |
| **confidence** | certain — verified independently, two ways |

`curation/wht_2026_falcon/build_rows.py` read `PMID = ""` with the comment
*"preprint, not in PubMed as of 2026-10-29"*. **The record had been live since
2026-06-29**, four months earlier.

Verified independently before applying, not taken from the reviewer:

```
esearch db=pubmed term="10.64898/2026.06.13.731774[doi]"   ->  42367940
esearch db=pubmed term="2026.06.13.731774[All Fields]"     ->  42367940
esummary 42367940 -> "FALCON: Closed-Loop Multi-Objective Optimization of Lipid
   Nanoparticles for Cell-Selective mRNA Delivery.", bioRxiv,
   doi 10.64898/2026.06.13.731774, entrez 2026/06/29, first author Toh WH
```

The second query is the one that matters: it searches PubMed for **the DOI this
curation already stored**, so the check needed no new information at any point.
PubMed has indexed bioRxiv deposits under the NIH Preprint Pilot since 2020, so
"it is a preprint" was never a reason for a blank PMID.

**Why this one stings.** `CLAUDE.md` §3 carried, on a single bullet, both the
previous instance (`QL_2026`, PMID 41145968, also already indexed) **and** the
rule it earned — *"Look the PMID up before recording a paper as unindexed"* —
and then asserted `WT_2026` was unindexed on the very next line. `CHANGELOG.md`
§11 did the same: the sentence recording the `QL_2026` correction ends *"`WT_2026`
(a preprint) is now the only experiment with a blank PMID."* The rule was written
down, in the same paragraph, and not applied to the one case in front of it.

**And `audit_baseline.json` was accepting 2,544 blanks with that false reason as
its justification** — so the gate was not merely blind, it was holding a
statement of fact that was wrong. A baseline entry's `why` is load-bearing: it is
the thing a future session reads instead of re-deriving.

**Applied.** `PMID = "42367940"` in `build_rows.py` with the query recorded
beside it; WT regenerated (validator **18/18 pass**, all 7 blocks round-trip,
selectivity residuals ≤ 4.0e-9); corpus rebuilt; `blank_Publication_PMID`
re-baselined **2,544 → 0** with a new `why`; the stale assertions in
`CLAUDE.md` §3, `CHANGELOG.md` §11 and §12, and the curation's `PROVENANCE.md`
all corrected. **Every paper in the corpus now has a PMID**, and C15 guards all
three facts (`blank_pmid_vitro`, `blank_pmid_vivo`, `WT_2026_pmid`) so the claim
cannot go stale silently again.

**The rule is strengthened, because "look it up" demonstrably did not work.**
Recording a paper as unindexed now requires a dated query result stored in
`provenance/esummary_cache.json`, not a statement in a comment.

---

## F-004 — settled by the round-2 reviewer: blank `N:P`, do not reconstruct it

| | |
|---|---|
| **carried from** | round 1, where it was left open pending independent confirmation |
| **status** | **evidence settled; the repair is the user's call and is NOT applied** |
| **severity** | changes a **value** — `IL_to_nucleicacid_chargeratio` on 129 rows; no endpoint moves |
| **confidence** | reviewer: certain on direction, certain that reconstruction is unsafe |

Round 1 found that `WT_2026`'s URS source file states each formulation twice —
four design parameters and four molar percents — and that on that one file of
seven they describe different formulations (0 of 129 self-agree, 126 derive
exactly one *other* row's percents, injectively). The open question was **which
block is authoritative**, because repairing the wrong half assigns real
measurements to the wrong formulations. That is the `AC_2025` question, and
`AC_2025` is why this project does not guess.

**The reviewer settled it, and the strongest argument needs no external source
at all.** `PEG_PEG+Chol` in the URS design block is **strictly monotone
decreasing across all 129 rows, zero violations**, while nothing else in the file
is sorted (every other column has ~60 ascending and ~60 descending violations).
A perfectly sorted column cannot arise from uniform random sampling, and no other
file has one. That is the signature of a spreadsheet sort applied to a
**contiguous column selection** with the sort key last — the design block moved,
not the percents. Three further ties agree: the paper's own text pins
`FALCON_ID` → response (`U0036` is named the top URS composition and carries the
file-maximum value), Supplementary Table 5 pins `FALCON_ID` → the **percent**
block to the printed decimal, and URS row 1's percents are byte-identical to
`F0162`'s in the clean E0 file.

**So the stored molar percents are right** — which is the half that matters, and
it means no structure or endpoint is affected.

**But `N:P` is not recoverable, and this is the part that changes the repair.**
`NP_ratio` sits inside the same column block, and nothing observable says whether
it was inside the sorted selection. For the one URS formulation the paper
tabulates there are **three mutually inconsistent candidates**: 7.199 (stored),
16.766 (reconstructed through the permutation), **6.0** (Supplementary Table 5) —
and no row in any of the seven files has an `NP_ratio` within 0.02 of 6.0.
Neither alignment carries signal against the response (ρ = +0.042 stored,
−0.041 reconstructed, against a random-pairing 95th percentile of 0.167). The
reviewer also found Suppl Table 5 is **itself** block-pasted in two places, so it
cannot rescue the reconstruction either.

**Reviewer's recommendation, which matches mine from before this round:** blank
`IL_to_nucleicacid_chargeratio` on the 129 `WT_2026_U*` rows, leave the percents
untouched, and put the reason in `curation_note`. **Not applied** — I flagged
this as the user's decision and it still is. Three options remain on the table:
blank (recommended), reconstruct, or email the Mao lab.

**Error class 40 added** for the general case, with the detection rule: derive
one parameterisation from the other and require **row-wise** agreement, because a
range check cannot see this — a permutation preserves the multiset, which is
exactly why this curation's `validate.py` recorded `PASS N/P range 3.33–37.12`
on 129 misassigned rows.

---

## F-028 — `Experiment_value_within_model` pools five normalisation units the source forbids comparing

| | |
|---|---|
| **found by** | the round-2 `WT_2026` reviewer (its E4) — **new, not a round-1 carry-over** |
| **status** | **open** — the repair is a schema decision, not a bug fix |
| **severity** | changes a **value** for any consumer of that column, 2,544 rows |
| **confidence** | reviewer: certain; the mechanism is arithmetic and I have not yet re-derived the table |

`Experiment_value_within_model` re-standardises inside (`Experiment_ID` ×
`Model_type`). `WT_2026` is the corpus's only experiment normalised per **source
file**, because the authors' own `dataset_info.txt` forbids comparing across the
seven files — and one `Model_type` spans several of them. `Ramos` alone covers
959 rows across URS + E0 + E1 + E4 + E5, and the derived column's per-file means
run from **−0.488 to +0.613** instead of being centred.

So the column manufactures a normalisation unit straddling five scales the source
says are incomparable, and **nothing declares it**: `value_normalization_unit`
documents the unit of `Experiment_value` only. A consumer who reaches for the
"within model" column precisely to remove an assay offset gets exactly the
cross-file mixing the `source_file` flag exists to prevent.

This is the error-class-22 shape one level up, created by a convenience column
rather than by a curation — **error class 41**. Options, none of them free:
blank the column wherever `value_normalization_unit != 'experiment'`; compute it
per (`Experiment_ID` × `Model_type` × source file) and say so; or add a second
column declaring the derived unit. Left open because it is a schema change
affecting how a documented column is defined, and `Experiment_value_within_model`
exists to answer a question the user asked for.

---

## Also reported by the `WT_2026` reviewer, dispositions

| its # | what | disposition |
|---|---|---|
| E3 | URS percent rows 1–3 are **not** uniform-random samples: `U0001` and `U0002` are byte-identical to `F0162`/`F0192`, and `U0003` is a round designed composition with `(IL+HL)% = 60.0000` exactly. The 3-and-3 unmatched rows landing at positions 1–3 is 1 in 357,760. | **open.** Certain on the observation, speculative on the cause, and the reviewer correctly declined to recommend exclusion. Needs the authors. Their values sit *above* the file mean, so they are not inert in a URS-vs-FALCON comparison. |
| E5 | The five quality columns are blank on all 2,544 rows although at least three conventions are in force: `18PG`'s expansion to SOPG is inherited from `YZ_2022` and never stated by this paper (1,114 rows); experiment 2's N/P is "numerically adjusted to account for the number of nitrogens from both SM102 and DSPC" per Suppl Table 3 (272 rows), recorded in `PROVENANCE.md` but not on the rows; `Mixing_method = handmixed` is a vocabulary mapping for robotic dispensing. | **open** — the F-001 failure mode recurring in the one curation that had no flags at all. Each needs a `structure_assumed`/`curation_note` pair. Deferred only because it is additive and the round is still running. |
| E6 | `PROVENANCE.md` still describes the pre-F-003 curation — "14 experiment blocks", "72 cols", stale corpus totals. | **already closed.** I checked by content rather than by the reviewer's line numbers: the file says 7 blocks and 79 columns throughout, and the two surviving strings are *deliberate* historical annotations (`*(was "per (dataset x cell type) -- 14 blocks"; corrected 2026-10-30, F-003)*` and a total explicitly scoped "at the time of this merge"). The reviewer's line numbers had shifted under the PMID edit. Left as written; the wording already pre-empts a future reviewer re-flagging it. |
| E7 | `validate.py`'s `PASS N/P range 3.33-37.12` cannot see F-004, because a permutation preserves the multiset. | **accepted and folded into error class 40** as the worked example of a check that passes on misassigned data. |
| §1 | F-003's repair re-derived independently and holds: all 2,544 rows round-trip at ≤2.7e-15, selectivity identity at 1.78e-15, Spearman(raw, stored) = **1.000000**, **0 of 180** sign inversions where round 1 measured 92. Both C10 anchors confirmed **byte-different** from the corpus strings and InChIKey-identical, so not the F-013 failure mode. All seven component names map to exactly one SMILES corpus-wide. Doses re-read from Methods verbatim, including the ng→µg conversion. | **F-003 closed.** |

---

## Wave 1 reports

## F-026 — the chance baseline in `CLAUDE.md` §7 matched no corpus this project has had

| | |
|---|---|
| **found by** | me, same docs-vs-data sweep |
| **status** | **applied** 2026-11-07 |
| **severity** | **metadata / manuscript risk** — it is a number that would go into Methods |
| **confidence** | certain |

§7 read: *"Chance baseline for S-way classification = 1/S (1/31 vitro, 1/17 vivo)."*

The live corpus is **S = 42 vitro / 19 vivo**. The corpus the frozen pipeline
actually ran on was **S = 34 vitro / 17 vivo** (§9, 2026-07-23). So `1/31`
described neither — it is the S from *before* the in-vitro corpus was
re-trimmed in July, left in place through every corpus change since.

Corrected to name the live S and to say explicitly that `results/` and
`reports/` are still the S = 34 / 17 run. **Quote the S of the corpus a number
was computed on, never a remembered one.**

### What the rest of the sweep found

Checked against the live corpus: vitro/vivo row, column and experiment counts;
`Original_LNPDB` dimensions; the 7 quality columns; the 987
`positional_ambiguity` rows; all six per-curation row counts; the 60
`SB_2024`/`SY_2025` rows; the 13 registered reagents. **All correct.**

Two of my own apparent findings were **my** error, not the docs', which is worth
recording because it is the fifth-instance pattern from round 1 pointing at me
again:

- *"288 IL names with >1 structure is stale, it is 292."* It is not. §3 says
  288 names carry different structures **in different experiments** — exactly
  288, measured both ways. 292 is the count of names with >1 structure
  regardless of experiment, which is what error class 30 cites. Two different
  quantities, both documents right, my check conflating them.
- *"the reagent registry has 4 entries, not 13."* I counted the top-level keys
  of `reagent_registry.json`, which is the audit's **output** (`what`,
  `registered`, `rows_mismatched`, `results`). The registry is in
  `reagent_registry.py`, and the audit prints `reagents clean: 13`.

Both are the same mistake as round 1's five: reading a secondary artifact
instead of the primary one.

---

## Wave 1 reports

Reports land in this directory as `report_<TARGET>.md`. Adjudication follows the
Wave 3 rule in `../REVIEW_PROGRAM.md` §5: nothing is applied because a reviewer
said so — I reproduce the evidence myself, and a reviewer's finding is a claim
until then.
