# Round 3 — reviewing the 2026-11-10 run's ~7,200 staged rows

Opened at the user's instruction: *"set up thorough subagent review cycles to
try to check literally every single thing added in this latest round (the 7k)."*

Follows `REVIEW_PROGRAM.md`: three waves, **Waves 1 and 2 read-only**, every
agent briefed from a **generated** packet rather than prose, and **I adjudicate
every finding myself before anything is edited.**

## What is in scope

**29 staged folders, 7,207 rows, ~1,148 distinct lipids (740 new).** None is
merged; `data/` is untouched at 29,195 x 79 / 1,803 x 79 and every folder
carries OPEN decisions, so `build_dataset.py` refuses them all.

This differs from rounds 1 and 2, which reviewed **merged** papers. Nothing here
is in `data/` yet, so `brief.py` — which reads the corpus — cannot generate the
packets. **New: `provenance/review/brief_staged.py`**, which reads each folder's
own emitted-rows file, identified by **schema** (`IL_SMILES` + a column starting
`Experiment_value`) rather than by filename. That matters: the folders used six
different names for that file, and a scan keyed on `06_rows.csv` found 14 of 27
and put the run's new-lipid count at **79 when it is ~740**.

## What the packets already surfaced, before any reviewer read them

| condition | folders |
|---|---|
| **no schema-valid row file at all** — values and structures never joined | `linyx_41851281`, `palczewski_41654520`, `sd_2026`, `xuq_32520527` (**255 rows**) |
| **row count disagrees with the queue** | `levkinpa_29558113` (45 vs 5), `lius_40500862` (23 vs 4), `xuq_32438474` (120 vs 90), `siegwart_40644303` (465 vs 505) |
| **72 `HepG2` rows present** where that whole arm was reported excluded as already-merged | `maohq_39375194` |

That is 8 of 29 folders flagged by a generator that only reads files. It is the
argument for generated packets in one line.

## Wave 1 + 2 — launched

| agent | scope | why first |
|---|---|---|
| **`maohq_39375194`** | 5,112 rows = **71% of everything staged**, on **one** structure (0 new lipids) | the 72 HepG2 rows; the `LC_2024`/`LC_2024_2` split under a single pooled z; F-057's `18PG` on 972 rows |
| **`siegwart_40644303`** | 505 rows, **489 structures** — the run's largest chemistry contribution | claims a **360-lipid overlap with `KZ_2016`** when the corpus's largest is 135 and Result 3's pooled n is 236; same lab as the F-013 epicholesterol incident |
| **the unmergeable four** | 255 rows across 4 folders | is the gap clerical or substantive? `sd_2026` is **384 structures against 18 values**; `palczewski_41654520` is **CC BY**, one of only four §8-clean sources |
| **CORPUS-A** | all 29 folders, cross-folder only | z-block honesty (error class 22), quality-column semantics (F-001/F-011/F-024), vocabulary drift, **`Experiment_ID` collisions** (`CS_2018` claimed twice; `YZ_2024` already taken), duplicate structures across folders |

## Still to launch

Wave 1 for the remaining 23 folders, in descending row order, and **CORPUS-B**
(docs-vs-data) and **CORPUS-C** (do the checks still have teeth — seed a bug of
each known class into a *copy* and confirm C1-C18 catch it).

## The stopping rule

Not "a reviewer found nothing". Per §6: **every known error class has an
automated check, and one adversarial pass found nothing outside those classes.**
If a reviewer only re-finds classes Tier 1 already covers, coverage has
saturated. If it finds a new class, add the check and go again.

## Standing rule for me

**Treat every reviewer finding as a claim.** Of one past reviewer's twelve, three
were real, one half-right, one flatly wrong. Five findings in round 1 were
conclusions drawn from a secondary artifact rather than the primary one and all
five were wrong — a pattern that has already landed on four agents and twice on
me today. Nothing is applied until I reproduce the evidence myself and a fresh
agent judges the repair.
