# Blind model eval — Sonnet 5.5 vs Opus on adversarial re-derivation

Run 2026-11-09, to answer the standing question: can a cheaper model take any
review role without losing quality?

## Read this first: the eval was flawed, and it was flawed in my favour

Two design defects, both mine, both found by the reviewers I was scoring:

1. **The "clean control" was not clean.** `emit_quality.py` was supposed to be
   the correct form, so that a defect reported in it would score as a false
   positive. Both arms found real defects in it, and I verified them: its
   `ARMS` dict declares `O1` as "isodecyl acrylate, branched isomer mixture"
   and the loop binds `arm_desc` and **never reads it**, so the caveat dies at
   emit — which is F-001's mechanism and F-025's content, recurring on the arm
   instead of the core. `schema_quality.json` names a branched-isomer-mixture
   reagent as exactly the `structure_assumed = 1` case. So **the false-positive
   rate is unmeasured**: this eval has no clean control and cannot produce one
   retrospectively.

2. **The answer key's mechanism was wrong on one case.** For the genuine
   pre-F-033 file I wrote the key from `CHANGELOG`'s account — "left the 0%
   line as three hard-coded literals" — when the code computes
   `height = bar["bottom"] - top_subpixel` from the **per-bar colour-mask
   bottom**, with the literal `3335` serving only as a ±6 px acceptance gate.
   Both arms read the code and got it right; the key was wrong. Logged as
   **error class 50**: ground truth transcribed from our own prose instead of
   from the code it describes. It is the 2026-10-30 secondary-artifact lesson,
   one level worse, because the secondary artifact was being used to score
   somebody.

So the only question this eval actually answers is the narrow one: **given a
file with a known real defect, did the arm find it and get the mechanism
right?** It says nothing about false positives.

## The cases

Both arms received the same three files under neutral names, read-only, with the
same brief, and no access to this repo's findings ledger. No two files were
near-identical, so neither arm could find an answer by diffing.

| given as | actually | real defect |
|---|---|---|
| `digitize_marginals.py` | `lin_aromatic_2026/digitize_fig2fgh.py` @ git `81763b5` — **genuine pre-F-033 code** | origin not tick-derived; uniform additive bias, mean −0.52 pp, 47/48 deltas negative |
| `build_conjugates.py` | hand-written replica of the F-013 mechanism | `RemoveAtom` on an atom adjacent to C3 ⇒ epicholesterol; identical formula, exact mass to 4 dp, and InChIKey **skeleton** block |
| `emit_quality.py` | intended as the clean control | **not clean** — see above |

## Result

| | Sonnet 5.5 | Opus |
|---|---|---|
| defects found, right mechanism | **3 of 3** | **3 of 3** |
| false positives | **unmeasurable** (no clean control) | **unmeasurable** |
| ran the code | 2 of 3 (the digitizer's inputs are absent) | 2 of 3, and **built a synthetic figure with known ground truth** for the third |
| independent re-derivation of the structure bug | yes — string substitution, no `RWMol` | yes — **two** routes, plus a per-step neighbour-list parity trace |
| additional real findings beyond the planted defect | 1 | **7** |
| tokens | 93k | 141k |

Both arms identified the F-013 mechanism correctly and proved it the way the
project's own policy demands — by an independent route, against an external
reference, with the stereo block as the discriminator. Both noted unprompted
that formula, mass and the skeleton block cannot see it. Sonnet also hit the
ring-closure collision on its first reference attempt, which is the bug that
motivated the fragment-assembly rule in the first place.

**Where they differed was evidence, not conclusion.** Sonnet reasoned from the
code and quoted the file's own reported numbers. Opus constructed a 2900×4200
synthetic page at the file's own hard-coded geometry with 48 bars of exactly
known height, recovered `px_per_percent` = 5.3500 against a ground truth of
5.3500, and then measured the origin error directly: **+0.079 pp** with a clean
bottom edge, **−0.288 pp** when the axis rule overdraws by 2 px — which
independently reproduces the file's own 0.28 pp estimate and the real page's
one-sided −0.52 pp. That is a different class of evidence from reading the
code, and it is what let it separate "the scale is fine" from "the origin is
not" with numbers.

## What the Opus arm found that nobody planted — and that is LIVE at HEAD

Four findings against the **current** `digitize_fig2fgh.py`, verified by me
after its report (see `review/FINDINGS_2026-11-09.md`): a panel that fails
detection is excluded from the gate and "48/48" can be asserted on 34; the
F-007 `min_area` repair was applied to one of the two masks; `len(db)` is never
checked so a spurious dark component displaces a real bar and reports as a
measurement; and the docstring claims the tick interval "is asserted to round
to `TICK_PCT`" when no assert mentions `TICK_PCT` at all.

## Conclusion, stated at the strength the evidence supports

- **For detecting a known defect and naming its mechanism, Sonnet 5.5 was not
  distinguishable from Opus on these three cases.** Three cases is a small
  sample and two of them it could not execute.
- **For depth of evidence and for finding what nobody was looking for, Opus was
  clearly ahead** — 7 additional real findings against 1, and it manufactured
  its own ground truth rather than accepting the file's account of itself.
- That maps onto the role split the planning document proposed (now
  `PIPELINE_V2.md` §5), and gives
  it some evidence rather than only a prior: **Sonnet for the mechanical roles
  and for checking a specific suspicion; Opus for Wave-1 adversarial
  re-derivation and Wave-3 adjudication**, which is where "find what nobody
  asked about" is the whole job.
- **Do not read this as a licence to run Wave 1 on Sonnet.** The thing Wave 1
  exists for is the 7-versus-1 column, not the 3-of-3 column.

## If this is re-run, fix the eval first

1. Build the clean control **by seeding nothing**: take a file at a commit where
   it was later confirmed correct by an independent reviewer, and verify it
   yourself before using it as a control.
2. Derive every answer-key mechanism **from the code**, citing file and line,
   and re-derive it before it scores anyone (error class 50).
3. Use more cases. Three is enough to show a model can do the task and not
   enough to rank two models.
