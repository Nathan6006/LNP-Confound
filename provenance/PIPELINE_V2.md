# Pipeline v2 — what changed on 2026-11-09, and how to curate from here

Replaces `OPTIMIZATION_PLAN.md`, which was the plan. This is the plan plus what
happened when it was built, including the three places the plan turned out to be
wrong. Read **§3 (how to curate a paper now)** before starting a new paper; the
rest is why it is shaped that way.

Nothing here changed `data/`. The corpus is byte-identical: **29,195 × 79 vitro,
1,803 × 79 vivo**.

---

## 1. The diagnosis

The ratchet applied to **checks** and not to **extraction code**.
`ERROR_CLASSES.md` and `audit.py` only grow, and that works. But every paper's
extraction code was written from scratch in its own folder — **50 scripts, 9,606
lines, of which 20 files and 3,853 lines read pixels or build structures** — so
a defect fixed in one folder stayed live in the others:

- **class 20** (scale or origin anchored on a panel box, border or glyph instead
  of detected ticks): four independent instances in four separately written
  files, the fourth of them (F-033) in a *sibling file of a folder where the
  class had already been found and fixed*;
- **class 32** (bar-top detector reading an overlaid replicate marker): `ql_2026`
  and `xd_2026_lipidai`, independently, days apart;
- **class 30** (unguarded corpus name → structure lookup): "repaired" in round 1,
  still live in three emitters when F-044 found it.

That is now **error class 49**. The fix is a shared, repeatedly-reviewed library,
and `curation/corpus_components.py` was the model for it.

---

## 2. What was built

### `curation/lib/` — one audited implementation per hazard

Every function **refuses** rather than returning a plausible number. Error class
28 is a detector miss falling through to a literal, and F-033 is what it costs.

| module | what it guarantees | retires |
|---|---|---|
| `axis.py` | scale **and** origin both fitted from detected tick positions; **there is no API to supply an origin separately**, because that is the shape F-033 took; asserts uniform spacing on the positions it was given, not on a comment; `forbid_target=` makes class 27 a type error | 20, 27, 33 |
| `bars.py` | full-width coverage over `nrow` rows **plus left/right flank agreement** — a centred marker is in neither flank; refuses on a miss | 28, 32, 35 |
| `heatmap.py` | grid from axis labels; raises above ~5 colour distance; **raises on an exact LUT tie** rather than letting `argmin` pick by row order; band **midpoints** with bounds carried; floor and ceiling flagged by **colour identity, symmetrically** | 16, 23, 43 |
| `chains.py` | vector-first (the PDF's own bond coordinates); **two independent methods must agree or it raises** | 48 |
| `assemble.py` | the acid loses its OH and the alcohol keeps its oxygen; atom-map attachment, never "index 0"; **stereocentre-count and radical assertions on every build** | 14, 15, 31 |
| `reagents.py` | refuses ambiguous names; cross-experiment skeleton-InChIKey comparison | 30, 46 |
| `emit.py` | the one writer of the quality columns; **refuses a literal `confidence`** | 29, 34 |

Two existing assets were promoted rather than rewritten: `ng_2025/stage01_vector.py`
(reads the PDF's bond coordinates, 20/20 exact) is now `chains`' primary path, and
`xd_2026_lipidai/stage02_trace.py` is its raster fallback.

### `provenance/test_lib.py` — 47 tests

24 seeded mutations, 21 controls, **2 historical replays**. A seeded mutation
proves a check *can* fire; a replay proves it fires on the thing that actually
happened. The F-013 replay rebuilds a sterol ester through `RemoveAtom` on an
atom adjacent to C3 and confirms the product is invisible downstream — identical
formula, identical exact mass to 4 dp, **identical InChIKey skeleton block
`RWTQCZGAMKTBRV`** — and that `check_stereo_preserved` catches it.

### New gates

| | what it does | state |
|---|---|---|
| **C17** `check_reagent_crosscheck.py` | one reagent name, one structure, across experiments. Gates excipients and registered reagents; **reports** the 288 generic IL labels, which legitimately collide | found **F-049** on its first run |
| **C18** `check_anchor_coverage.py` | C10 anchor coverage and copy posture per folder | baselines 3 known gaps, each naming work not done |
| **C15 extended** | 33 → **50** claims, across `CLAUDE.md`, `provenance/README.md`, `curation/PROTOCOL.md`; each claim carries a `file` and the file's existence is checked | caught two of my own mis-derived entries immediately |
| **`DECISIONS.md` merge gate** | `build_dataset.py` refuses a folder with an `OPEN` entry | verified: exit 1, no partial write, rebuild byte-identical |
| **`new_curation.py`** | scaffolds a folder, stage stubs that import `curation/lib`, the intake gate, `make_reference.py`, `DECISIONS.md` | stubs refuse to be skipped |

`PROTOCOL.md`'s per-paper checks are renumbered **P1–P9**; **C1–C18** are
corpus-level. The two schemes previously collided (P1 was HRMS, C1 is SMILES
parse), with per-folder QC using one and review reports the other.

---

## 3. How to curate a paper now

### Stage 00 — acquire, and record it

```bash
python3 -I provenance/new_curation.py <folder> <paper.pdf> [si.pdf ...]
```

Fill `00_acquisition.json`: URL, sha256, bytes, UTC timestamp, **and the licence
string off the landing page**. §8 of `CLAUDE.md` needs every source
CC0/CC-BY-compatible before release and that is far cheaper to capture now.

**Getting the files** — measured behaviour is in `BROWSER_ACQUISITION.md`:

- A **403 on a deep link is evidence about the request, not the site.** Visit the
  publisher's homepage first, then the article. This is how ScienceDirect went
  from "hard block" to working.
- Springer Nature (71 of the corpus) needs **no login**: `static-content.springer.com`
  serves Source Data to a plain request. The browser is only a URL enumerator.
- Wiley, PNAS, ACS, RSC need a **browser** but no login. Cloudflare's challenge
  self-resolves in ~8 s; do not treat the first 403 as final.
- Elsevier needs the **institutional login to enumerate** supplementary assets
  and to fetch the article PDF — but `ars.els-cdn.com` then serves those assets
  **with no session**.
- **Never scrape PubMed.** It serves a reCAPTCHA. Use `fetch_esummary.py`.
- **Check for publisher figure rasters before rendering anything.** Elsevier
  publishes `https://ars.els-cdn.com/content/image/1-s2.0-<PII>-gr<N>_lrg.jpg`,
  open-access, measured at **3398 × 4336 px for a single figure** against ~2500 ×
  3300 for the same figure inside a 600 dpi page render. Higher resolution, no
  re-rasterisation, free. Nature's equivalents are on `media.springernature.com`.

### Stage 00b — THE INTAKE GATE. Fill `00_value_sources.json` before writing any digitiser.

For every panel or table that supplies a value, name the tier it came from:

1. publisher Source Data file → 2. SI tables → 3. author repo / Zenodo / Dryad →
4. `pdfplumber` text table → 5. **pixels, the last resort**

Anything recorded as `5_pixels` must say **why tiers 1–4 failed**, in a sentence
a reviewer can check. "Not available" is not such a sentence.

This is the highest-value step in the whole document:

| | value source | digitisation defects |
|---|---|---|
| `NG_2025` | publisher Source Data `.xlsx` | **0** |
| `WT_2026` | the authors' own repository | **0** |
| `QL_2026`, `LL_2026`, `ZT_2026`, `XD_2026` | digitised panels | **8** (F-006, F-008, F-017, F-018, F-019, F-020, F-033, F-040) |

**Once you have numbers from tiers 1–3, use the figure to VERIFY, not to
extract.** `NG_2025`'s `stage06_heatmap.py` mapped all 310 drawn cells onto
workbook values at ρ ≥ 0.998 per panel, independently confirming block order, row
order and panel assignment. That is strictly stronger than anything available
when the figure *is* the source.

### Stages 01–07 — use the library

Import from `curation/lib`; do not write a fourth copy of a geometry helper. A
paper's own constants, panel bounds and library logic stay in the folder, next to
the figure, where a reviewer reads them.

Still true, and still the hardest-won rules:

- **Mass cannot see connectivity.** HRMS confirms the formula and therefore arm
  *counts*; it is blind to which nitrogen an arm sits on. Read the synthesis.
- **Name tokens are not chain lengths.** `KS_2024`'s `C24` is two C4 tails.
  `chains.check_name_token` reports the mismatch; it never corrects it.
- **Passing every printed statistic is necessary, not sufficient.** `QL_2026`
  reproduced every printed hit rate on bars up to 28% too high; `XD_2026` landed
  12/12 inside its stated band on tops 0–21 px wrong. Geometry must be checked
  geometrically.
- **Better than a printed statistic: a round number the authors DREW and never
  printed** — a dashed threshold, an axis limit at a round value, a colorbar's
  own tick notch. External to our pipeline in a way a recomputed marginal is not.
- **Flag on the row, not in a QC file.** A flag that lives only in a per-folder
  QC file is not a flag, it is a memo (F-001, 1,424 rows).

### Stage 08 — `08_contradictions.csv`

Every paper curated so far has contradicted itself; `NG_2025` alone yielded
seven. Record (claim, locus A, locus B, what we emitted, why). An empty file
means nobody looked.

### The anchor — `make_reference.py`, and C16

C10 only tests anything if the anchor was built by a **different code path that
is executed, not asserted**. That has been claimed twice and false twice (F-013,
a byte copy; F-037, the file written to fix F-013 re-copied the other half).
Prefer a route sharing as little as possible with the enumerator. If the paper
draws nothing, store an empty `reference_SMILES` so the absence is visible.

### `DECISIONS.md` — and the merge blocks on it

The pipeline detects that a source is silent; it does not get to decide what the
source meant. Queue these, never resolve them alone:

- which half to repair when a **name and a structure disagree** — `AC_2025`'s
  first repair was confident, internally consistent, and would have put **419
  measured values on the wrong molecules**;
- whether a `confidence = low` library belongs in the corpus;
- a value with several mutually inconsistent candidates (F-004's `N:P`);
- anything that would change a published number.

### Before merge

```bash
python3 -I provenance/build_dataset.py
python3 -I provenance/audit.py                      # C1-C14
python3 -I provenance/check_docs.py                 # C15
python3 -I provenance/check_anchor_independence.py  # C16
python3 -I provenance/check_reagent_crosscheck.py   # C17
python3 -I provenance/check_anchor_coverage.py      # C18
python3 -I provenance/test_lib.py                   # curation/lib has teeth
python3 -I provenance/test_checks.py .              # C1-C14 have teeth (~20 min)
```

---

## 4. Three things the plan got wrong

Recorded because the corrections are the useful part.

1. **"Elsevier is a hard block."** It is not. My 403 came from hitting an article
   route cold. Homepage first, then the article, and it loads. Both my attempts
   had been deep links.
2. **"The browser is for getting past institutional SSO."** The binding
   constraint is **bot-blocking**, which hits the article page long before any
   paywall. The browser unblocks Wiley/PNAS/ACS/RSC with no login at all.
3. **The reagent cross-check is error class 46, not 40.** I wrote "class 40" in
   four places. Corrected — and it is exactly the numbering hygiene the ledger
   warns about.

And one about the model eval, in `MODEL_EVAL.md`: **my clean control was not
clean and my answer key's mechanism was wrong** (now error class 50 — ground
truth transcribed from our own prose instead of from the code). Both reviewers
read the code and got it right; the key did not.

---

## 5. Open, and what to do next

### Open findings — `review/FINDINGS_2026-11-09.md`

- **F-049** `JW_2024`'s PEG-lipid: 1,873 rows call it `DMPE-PEG2000` while
  carrying the exact string 41 experiments call `DMG-PEG2000` — no N, no P, where
  a real DMPE has a phosphoethanolamine head. Both sides C122H242, so nothing
  caught it. **Needs PMID 37985700.** This is the `AC_2025` question: do not
  repair on a 41-to-1 count.
- **F-050 … F-053 are FIXED** (2026-11-09) in `lin_aromatic_2026/digitize_fig2fgh.py`,
  bar values byte-identical, `build_rows.py` confirmed not to read that file. The
  gate now runs **before** the artifacts are written — found while seeding the
  check to confirm it fires, because an `IndexError` at the CSV write was masking
  the assert.

  **F-051's repair is not the one the reviewer proposed.** It said to lower the
  light mask's `min_area` to 120 to match the dark mask. Measured: at 120 the
  light colour also matches axis text, every panel fails, **0 of 48 rows are
  emitted**. The asymmetry is load-bearing. The floor is now *auditable* instead
  — largest rejected vs smallest kept, asserted to be cleanly separated
  (measured 1,348–2,260 px² vs 7,720–18,012). **Treat reviewer findings as
  claims — including their repairs.**

### Not done

- **Item 5 of the old plan: refocus the review waves.** At the next Tier-4 round,
  point Wave 1 at the "not yet automatable" classes, the paper's own printed
  statistics, and the contradiction hunt — not at ground C1–C18 already covers.
  **CORPUS-C should not be an agent**; its job is `test_checks.py`.
- **Anchor coverage.** `ng_2025` has 6 anchors all flagged as possible copies and
  no generator — and it is the best-placed folder to fix this, because
  `stage01_vector.py` is genuinely independent of the raster enumeration.
  `ql_2026` has 1 anchor for 344 structures.
- **`amino_lipid_fg_2026`'s blocker is only partially relieved.** The five Nature
  Source Data workbooks give all 30 amino-acid blocks and **20** named members;
  only 20 of the 960 appear anywhere. Do not record it as resolved.
- **Model tiering.** `MODEL_EVAL.md`: Sonnet 5.5 matched Opus at finding a known
  defect and naming its mechanism (3/3 each), and was far behind at finding what
  nobody planted (1 extra finding vs 7). Use it for mechanical roles and for
  checking a specific suspicion; **keep Opus on Wave 1 and Wave 3.**
- The standing corpus decisions in `CLAUDE.md` §8, unchanged: re-running the
  frozen pipeline, `SB_2024`/`SY_2025`, the 48 duplicate SMILES.
