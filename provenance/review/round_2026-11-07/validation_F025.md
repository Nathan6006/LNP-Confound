# Validation — F-025: `confidence` on the 15 `ZT_2026` commercial-core cholesterol conjugates

**Validator:** read-only validation agent, round 2026-11-07
**Verdict: APPLY WITH MODIFICATION.**

The cap to `medium` is correct in direction and is the value three independent
project documents already say these rows should carry — apply it. Two
modifications: **(1)** remove the confidence from `COMMERCIAL_CORES` entirely
rather than cap its output, because nothing else reads it and leaving it there
leaves the next reader the same trap; **(2)** `Lipid 14-CB1/2/3` should go to
**`low`**, not `medium` — it is the only one of the five cores where the
conjugation site is a genuine selection between two chemically distinct
nitrogens, it is the one that was already resolved the wrong way once (F-015),
and this curation's own rule (`enumerate.py:309`) assigns `low` to exactly that
condition. A flat `medium` across all 15 is an improvement on the status quo but
still erases the one distinction that matters.

No baseline metric moves under either version. No structure, value or row count
changes.

---

## 1. The finding reproduces, and the description is accurate

### 1.1 The 6/9 split and the identical note — confirmed in both primary files

`data/data_vivo.csv`, `Experiment_ID == ZT_2026`: 51 rows, of which 15 match
`-CB`. All 15 carry `structure_assumed = 1.0`, blank `positional_ambiguity`,
blank `value_is_bound`, and **one** distinct `curation_note`:

> `Cholesterol conjugation site inferred: the paper draws no -C conjugate for any commercial core.`

| `confidence` | rows |
|---|---|
| `high` | `Lipid 14-CB1/2/3`, `TT3-CB1/2/3` (6) |
| `medium` | `MC3-CB1/2/3`, `SM-102-CB1/2/3`, `OF-2-CB1/2/3` (9) |

Identical in `curation/new_curation.csv` (rows 486–524), the file
`build_dataset.py` merges. So the split originates upstream of the merge, as
described.

### 1.2 The mechanism — confirmed by reading and by running the code

`curation/build_new_curation.py`:

- **line 92** `COMMERCIAL_CORES = {name: (core_SMILES, confidence)}` — `TT3`
  `"high"`, `Lipid 14` `"high"`, `OF-2` `"medium"`. The comment immediately
  above (lines 86–91) states plainly what that confidence is *about*: TT3 and
  Lipid 14 are "machine-counted" because their chains are saturated, OF-2's
  unsaturated tails "were read by eye from the alkene pattern and carry a lower
  confidence flag." It is a **parent-core tracing** confidence.
- **line 373–374** `elif core.endswith("-C") and core[:-2] in COMMERCIAL_CORES:
  parent_smiles, conf = COMMERCIAL_CORES[core[:-2]]` — that core confidence is
  taken verbatim.
- **line 378–385** the same branch then sets `"structure_assumed": "1"`
  unconditionally with the note above.
- **line 368–372** the `PARENT` branch (`MC3-C`, `SM-102-C`) hard-codes
  `conf = "medium"`.
- **line 422–426** `row["confidence"] = res[2]` writes it to the row.

Every part of the finding's description of the mechanism is correct. I found
nothing in it to correct.

### 1.3 One addition the finding does not mention, and it is the sharper defect

`quaternize()` (line 346) is `smi, _note = _EN.quaternize(core_smiles, n_link)`.
The shared implementation returns a **site classification** alongside the SMILES,
and this branch throws it away into `_note` and substitutes hand-written text.

Running the shared function on all five commercial cores
(`curation/siegwart_42634137/enumerate.py::_quaternization_site`, invoked
directly on the live `COMMERCIAL_CORES` strings and the corpus MC3/SM-102
structures):

| core | chosen N idx | N in ring | returned note |
|---|---|---|---|
| TT3 | 6 | False | `ASSUMED: no ring N and the paper draws no structure` |
| Lipid 14 | 1 | False | `ASSUMED: no ring N and the paper draws no structure` |
| OF-2 | 8 | False | `ASSUMED: no ring N and the paper draws no structure` |
| DLin-MC3-DMA | 43 | False | `ASSUMED: no ring N and the paper draws no structure` |
| SM-102 | 19 | False | `ASSUMED: no ring N and the paper draws no structure` |

That string is not incidental. In `curation/siegwart_42634137/03_enumerated.csv`
it is the *verbatim* reason attached to **24/24** rows carrying
`confidence = low`, and the folder's confidence rule is literally

```python
"confidence": ("high" if (hk in ("C5", "C6") and ak == "O2")
               else ("low" if site_note.startswith("ASSUMED") else aconf)),
```

So the project's single implementation already classifies all five commercial
cores into the `low` condition, and the commercial-core branch discards that
return value and writes a different, softer label from a table that means
something else. The conflation the finding identifies is real; it sits on top of
a discarded classification.

---

## 2. Is `medium` the right direction? Yes. Alternatives (a) and (c) fail.

### 2.1 The paper draws no `-C` conjugate for ANY commercial core — all five, not some

This is the question the task asked me to settle from the PDFs, because if the
conjugate were drawn for TT3 and Lipid 14 the 6 `high` rows would be correct and
the repair wrong. It is not.

**Main text, Fig 1b (p. 3, rendered at 600 dpi and read):** the panel is a
component table in three bands. "Commercially available cores" draws **SM-102,
OF-2, TT3, Dlin-MC3-DMA (MC3), Lipid 14 — all five unconjugated**, as bare
cores. "Cholesterol tails" draws **B1, B2, B3 separately, as the free
ω-bromoalkyl cholesteryl esters** (each `Br`-terminated). The conjugate of a
commercial core with a cholesterol tail is drawn nowhere in the panel. Fig 1d's
heatmap and Fig 1e's bars label them `MC3-C`, `SM-102-C`, `Lipid 14-C`,
`OF-2-C`, `TT3-C` as names only, with no structures.

**Main text, Fig 2a (p. 5, rendered and read):** the only drawn conjugate in the
paper is **C6O2B2**, a *synthesized*-core lipid, showing its N⁺ on a piperazine
ring nitrogen. That is the anchor the `high` tier in `03_enumerated.csv` is built
on, and it says nothing about any commercial core.

**Supporting Information (all 23 pages' text layer enumerated):** NMR is
supplied for **C6O2, C6O1, C5O2, B1, B2, B3 and C6O2B2** (Figs S1–S7) — the
synthesized cores, the three free cholesterol tails, and the one lead conjugate.
**There is no NMR, no MS, and no structure for any of `MC3-C`, `SM-102-C`,
`TT3-C`, `Lipid 14-C`, `OF-2-C`.** `TT3`, `Lipid 14` and `OF-2` do not occur in
the main PDF's text layer at all — they exist only inside the Fig 1b raster.

**Consistent with this, `curation/siegwart_42634137/reference_structures.csv`
holds exactly one anchor, `C6O2B2`.** There is no commercial-core conjugate
anchor because there is nothing in the paper to anchor against. C10 is therefore
unaffected by this repair.

→ **Alternative (a) is dead.** No commercial-core conjugate is drawn or
characterized. Under this folder's own meaning of the word — `high` = "drawn in
Fig 2a" / parent machine-counted *and* site not assumed — no commercial-core
conjugate can earn `high`. The 6 `high` values assert the opposite of the only
evidence there is, which is the exact shape of F-011.

### 2.2 Three project documents already say `medium`. The code is the outlier.

1. **`build_new_curation.py`'s own `quaternize()` docstring** (lines 228–230):
   *"The paper draws NO -C conjugate for any of the five, so the conjugation site
   is inferred for every one of them. These rows carry `confidence = medium` and
   `structure_assumed = 1` with the reason on the row."* The function's
   docstring contradicts the branch 140 lines below it.
2. **`curation/siegwart_42634137/PROVENANCE.md` §4** (the F-011 write-up):
   *"the 15 commercial-core conjugates read `medium` with `structure_assumed = 1`
   — the paper draws no `-C` conjugate for any of the five, so every one of those
   sites is inferred."*
3. **`curation/new_curation_qc.json` → `IL_SMILES_resolution.confidence_basis`**:
   *"medium = the core is C5/C6 … **or the cholesterol CONJUGATION SITE is
   inferred for a commercial core (the paper draws no -C conjugate at all)**, or
   the parent has unsaturation the chain counter cannot validate (OF-2)."*
   The **same JSON file** then lists `Lipid 14-CB1/2/3` and `TT3-CB1/2/3` under
   `by_confidence.high`. It is internally self-contradictory: its stated rule
   says medium, its enumeration says high.

The `high` is an implementation slip, not a considered assessment. Nothing
anywhere argues for it.

### 2.3 Alternative (c) — "leave `confidence` alone, `structure_assumed` carries it" — fails

`provenance/schema_quality.json` defines `confidence` as *the structure's*
confidence as assessed by the curation that built it, and CLAUDE.md §3 adds
"Must be that lipid's own assessment." The emitted structure on these rows is the
**conjugate**, not the core. A consumer filtering `confidence == 'high'` to find
structures the curation vouches for would pull in six molecules whose
attachment chemistry appears in no figure and no spectrum. `structure_assumed`
cannot rescue that: the two columns are read independently, and F-011 is the
precedent — a wrong `confidence` beside a correct caveat was judged *worse than
no label*, because it asserts the opposite of the evidence.

### 2.4 Alternative (b) — "all 15 should be `low`" — over-corrects for 12 of 15

This is where the folder's `ASSUMED → low` rule and the documents' `medium`
genuinely disagree, and the rule's proxy is what is wrong. `_quaternization_site`
decides `ASSUMED` purely from *"is the chosen N in a ring"*, which is a proxy for
*"is it the piperazine N drawn in Fig 2a"*. For the commercial cores that proxy
produces false positives, because for four of five there is no choice to make.

Counting distinct tertiary-amine alkylation sites with RDKit symmetry classes
(`CanonicalRankAtoms(breakTies=False)` over non-amide, non-aryl N with no H):

| core | tertiary-amine N | **distinct symmetry classes** | is the site genuinely a choice? |
|---|---|---|---|
| DLin-MC3-DMA | 1 | 1 | **No** — one amine, forced |
| SM-102 | 1 | 1 | **No** — one amine, forced |
| TT3 | 3 | **1** | No — the three N(C₁₂)₂ arms are equivalent; mono-alkylation at any gives the same molecule |
| OF-2 | 2 | **1** | No — the two arms off the diketopiperazine are equivalent |
| **Lipid 14** | 2 | **2** | **Yes** — terminal NMe₂ vs. the central trialkyl N, chemically distinct |

For MC3 and SM-102 the row's note — "conjugation site inferred" — is in fact an
**over-claim**: there is exactly one tertiary amine, and the paper states in text
that the cores were "functionalized with cholesterol-derived alkyl bromides to
generate a cholesterol-based quaternary amino lipid structure" (p. 2). Given that
sentence, the product is determined. `low` would be wrong for those.

So blanket `low` is not supported. `medium` is the right landing for 12 of the 15.

### 2.5 …but `Lipid 14-CB1/2/3` is the real thing, and a flat `medium` hides it

Lipid 14 is the single commercial core with two genuinely distinct candidate
nitrogens and no drawing. The current rule picks atom idx 1 — the terminal
dimethylamino — and this is **not a case where the heuristic has never been
wrong**: `PROVENANCE.md` §3 records that `Lipid 14-CB1/2/3` were previously
*"alkylated at the hindered central nitrogen"*, because the β-arm counter walked
only `N-C-C-C(=O)` and could not see Lipid 14's `N-C-C-O-C(=O)` β-acyloxy arm, so
both nitrogens scored 0 and the answer fell through to an atom-index tiebreak —
i.e. it was a function of SMILES spelling. The rule is better now and the current
answer is the chemically sensible one (NMe₂ is both more nucleophilic and far
less hindered than an N flanked by two C₈ ester chains), but it remains an
*inference the paper does not confirm*, on the one molecule where the inference
has a live alternative that was once selected.

Under a flat `medium` with an identical note across all 15, nothing on the row
distinguishes "site forced by a single amine" (MC3, SM-102) or "forced by
symmetry" (TT3, OF-2) from "chosen between two distinct nitrogens by heuristic,
previously chosen wrongly" (Lipid 14). That is precisely the information the flag
columns exist to carry, and `low` is what this curation's own rule assigns to it.

**Note the direction this cuts.** The current labelling is not merely
inconsistent, it is *inverted*: the three rows with the most structural
uncertainty (`Lipid 14-CB*`) carry the project's highest confidence, while the
six with the least (`MC3-CB*`, `SM-102-CB*`, single forced amine) carry `medium`.

---

## 3. Is a cap the right shape? No — remove the confidence from the table.

The finding's own instinct here is right and I endorse it over the proposal as
written. Reasons:

- **Nothing else reads it.** See §4 — line 374 is the only consumer in the repo.
- **A cap leaves the mislabelled value in place.** `COMMERCIAL_CORES["TT3"][1] ==
  "high"` would still read, to the next person, as "the TT3 conjugate is
  high-confidence", with the capping logic a branch away. That is the half-fix
  the task warns about.
- **No information is lost.** The tracing provenance that the confidence encodes
  — TT3/Lipid 14 machine-counted at bond-gap CV < 0.02, OF-2's unsaturation read
  by eye — is already written out in full in the comment at lines 86–91, in
  `confidence_basis`, and in PROVENANCE.md. The tuple is not where that lives.
- **It matches the project's established shape.** Both sibling enumerators
  implement weakest-link semantics rather than pass-through:
  `enumerate.py:309` forces `low` when the site note starts with `ASSUMED`,
  overriding the component confidence `aconf`; `lin_aromatic_2026/enumerate_library.py:198`
  starts at `"high"` and downgrades to `"medium"` for a `SUPPLIED` tail. The
  commercial-core branch is the **only** place in the repo that takes a component
  confidence and applies no downgrade at all while simultaneously setting
  `structure_assumed = 1`.

**Recommended shape:**

```python
COMMERCIAL_CORES = {"TT3": "...", "Lipid 14": "...", "OF-2": "..."}   # SMILES only
```

and in the branch, derive the confidence from the conjugate's own situation
rather than the core's — the number of *distinct* candidate tertiary amines,
which is computable from the parent and is the thing the label is about:

- more than one distinct candidate N → `low` (Lipid 14)
- exactly one distinct candidate N → `medium` (TT3, OF-2, MC3, SM-102)

If a derived rule is judged too clever for a 15-row branch, a flat literal
`"medium"` with `"low"` spelled out for `Lipid 14` is acceptable — but it should
be written **in the branch beside `structure_assumed = "1"`**, where the reason
is visible, not carried in a table that means something else. In either case
`PARENT`'s hard-coded `"medium"` (line 372) then becomes the same literal from
the same place, and the two branches stop disagreeing by accident.

**One thing to preserve if the tuple goes:** the OF-2 eye-read caveat currently
has no other route onto the row. Today it is invisible anyway (OF-2-C and MC3-C
both read `medium` with the same note), so nothing regresses — but if per-core
notes are differentiated as §5 suggests, OF-2's should mention it.

---

## 4. Scope check — the conflation is confined to this one branch

`grep -rn "COMMERCIAL_CORES"` across the repo (excluding `curation/ng_2025/`,
which I did not read):

| hit | what |
|---|---|
| `curation/build_new_curation.py:92` | the definition |
| `curation/build_new_curation.py:373–374` | **the only read** |
| `provenance/review/round_2026-11-07/FINDINGS.md:233,239,240,257` | this finding's own write-up |
| `provenance/review/round_2026-10-30/validate_F014_F015.md:194,233,327,328` | round-1 prose |
| `provenance/review/round_2026-10-30/validate_F013.md:309` | round-1 prose |

**One live consumer.** No second copy of this branch exists — the duplication
hazard that burned round 1 (`build_new_curation.py` carrying its own copies of
the ester assembly and the site rule, F-014/F-015) was already closed: the file
now imports `enumerate.py`'s implementation, and `COMMERCIAL_CORES` was never
duplicated.

I also swept every `curation/**/*.py` for hard-coded confidence literals to see
whether the same component→derived conflation recurs elsewhere:

- `curation/siegwart_42634137/enumerate.py:47–49, 309` — component confidences
  for O1/O2 alkyls, **correctly** combined with weakest-link semantics (O1's
  unconfirmed branch position is why `C5O1*`/`C6O1*` are `medium` while
  `C5O2*`/`C6O2*` are `high`). Not the same defect.
- `curation/lin_aromatic_2026/enumerate_library.py:56–57, 198–200` — same
  pattern, same correct downgrade for `SUPPLIED` tails. Not the same defect.
- `curation/xd_2026_lipidai/stage07_emit.py:122, 132, 158` — `high` for the three
  benchmarks, `medium` for E1–E9, which is the documented `XD_2026` posture
  (CLAUDE.md §9, 2026-10-29: E1–E9 have no external anchor of any kind). Not a
  conflation.
- `curation/lin_aromatic_2026/enumerate_147.py:111` — `{**r, "confidence": "high"}`
  injected at the call site. This is a *hard-coded literal* of the F-011 family,
  though not the COMMERCIAL_CORES conflation. **Outside the scope of this
  validation and I did not verify whether it is justified for the 147-lipid
  screen — flagging it as worth a look in this round rather than asserting it is
  a bug.**

---

## 5. Downstream — nothing depends on these 6 being `high`, and no baseline moves

I recomputed `audit.py`'s C13 metrics directly from `data/data_vivo.csv` with the
live `schema_quality.json`, replicating the function body, under three states:

| state | `bad_values` | `rows_carrying_a_flag` | `rows_with_a_caveat` | `flagged_without_a_note` | `confidence_without_a_note` |
|---|---|---|---|---|---|
| current | 0 | 63 | 54 | 0 | 6 |
| proposed (6 → `medium`) | 0 | 63 | 54 | 0 | 6 |
| modified (TT3 → `medium`, Lipid 14 → `low`) | 0 | 63 | 54 | 0 | 6 |

All four gated `vivo.C13_quality_vocabulary.*` values in
`provenance/audit_baseline.json` (`bad_values` 0, `rows_carrying_a_flag` 63,
`rows_with_a_caveat` 54, `confidence_without_a_note` 6) are reproduced exactly
and **none of them moves.** Why:

- `rows_carrying_a_flag` / `rows_with_a_caveat`: all 15 rows already carry
  `structure_assumed = 1`, so they are already counted as a caveat; changing the
  confidence from one populated value to another populated value cannot move
  either count.
- `confidence_without_a_note`: all 15 have a non-empty `curation_note`, so
  demoting them adds nothing to this counter.
- `bad_values`: `medium` and `low` are both in the schema's allowed vocabulary.

**The baseline `why` strings do NOT need editing, but one should be re-read.**
`audit_baseline.json` line 339 says *"9 of the 63 carry only confidence = high (6
ZT_2026, 3 XD_2026 benchmarks)"*. I checked which rows those are: the 6 ZT_2026
are `C5O2B1/B2/B3` and `C6O2B1/B2/B3` — the **enumerated** Fig-2a-anchored
lipids, not the CB rows (which carry `structure_assumed` and so are not
"only confidence = high"). That `why` stays true after the repair. Likewise line
349's F-024 set (`C5O1*`/`C6O1*` medium with empty note) is untouched.

Other downstream:

- **`provenance/test_checks.py:90`** — the C13 mutation is
  `confidence="med"` on a **vitro** row, a vocabulary test. Unaffected; still has
  teeth.
- **`provenance/structure_check.py` / C10** — `reference_structures.csv` holds
  only the `C6O2B2` anchor. No commercial-core conjugate is anchored, because
  none is drawn. Unaffected.
- **`curation/new_curation_qc.json`** — `IL_SMILES_resolution.by_confidence` will
  change (`Lipid 14-CB*` and `TT3-CB*` move out of `high`). This file is
  **regenerated by `build_new_curation.py`'s `qc()`**, so it updates itself; it
  needs no hand-edit, and its `confidence_basis` prose is already correct for the
  `medium` outcome. If `Lipid 14` goes to `low`, `confidence_basis` **does** need
  one clause added, since its current text routes all commercial-core conjugates
  to `medium`.
- **`curation/siegwart_42634137/PROVENANCE.md` §4** — its sentence *"the 15
  commercial-core conjugates read `medium`"* becomes true under the flat cap. If
  `Lipid 14` goes to `low`, that sentence must be amended in the same commit, or
  the docs-vs-data contract breaks in the other direction.
- **Rebuild chain.** The fix is in `curation/build_new_curation.py`, so it needs
  `build_new_curation.py` → `build_dataset.py` → `audit.py`, not an edit to
  `data/data_vivo.csv` or `curation/new_curation.csv`. No `Experiment_value`,
  `IL_SMILES`, descriptor or row count is touched; `ZT_2026` stays 51 vivo / 486
  vitro.

---

## 6. Summary of what I confirmed, and what I am asserting

**Confirmed against the primary objects:**
- 6/9 split and the single identical note, in `data/data_vivo.csv` **and**
  `curation/new_curation.csv` (read the data).
- The mechanism at lines 92 / 373–374 / 378–385 / 422 (read the code), plus the
  discarded `site_note` at line 346 that the finding does not mention.
- All five commercial cores return `ASSUMED: no ring N and the paper draws no
  structure` from the shared `_quaternization_site` (**ran the function**, not
  read its docstring).
- No `-C` conjugate for any commercial core is drawn in Fig 1b or Fig 2a, or
  characterized anywhere in the 23-page SI (**rendered pp. 3 and 5 at 600/400 dpi
  and looked at them**; enumerated the SI text layer).
- Symmetry-class counts per core by RDKit (**ran it**).
- All four gated vivo C13 metrics, recomputed and matched, unmoved under both
  candidate repairs (**ran it**).

**I am asserting, as judgment rather than measurement:**
- That `low` is the right label for `Lipid 14-CB1/2/3`. This follows from this
  curation's own rule applied to a genuine 2-way site choice, but it is a
  judgment about where a three-level vocabulary should cut, and it is the one
  place where the folder's code rule and the folder's prose disagree. If the
  decision is to keep all 15 at a flat `medium`, that is defensible — but then
  `Lipid 14-CB*`'s `curation_note` must say that the site was selected between
  two distinct nitrogens by heuristic and was previously selected wrongly
  (F-015), so the distinction survives somewhere on the row.
- That `curation/lin_aromatic_2026/enumerate_147.py:111`'s hard-coded
  `"confidence": "high"` is worth a look. I did not verify it and make no claim
  about whether it is wrong.
