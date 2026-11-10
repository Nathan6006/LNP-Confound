# PROVENANCE — curation/songx_39664425

**PMID 39664425** · PMC11628822 · doi `10.1016/j.apsb.2024.08.032`
*Pioneering integration of combinatorial chemistry and machine learning to accelerate
the development of tailored LNPs for mRNA delivery*
He X, Wang P, Qing L, Song X — Acta Pharm Sin B 14(12):5540–5543 (2024), CC BY-NC-ND.

## Outcome

**0 rows emittable.** `article-type="letter"`, PMC subject heading **`Highlight`** — a
research-highlight commentary on someone else's paper. It contains no original
measurement, no table, no SI and no structure of its own.

The work it describes is **Li B, Raji IO, Gordon AGR, et al., "Accelerating ionizable
lipid discovery for mRNA delivery using machine learning and combinatorial chemistry",
Nat Mater 23:1002–1008 (2024), PMID 38740955** — which is **already in the corpus as
`BL_2024`** (560 in-vitro HeLa rows + 51 in-vivo mouse-muscle rows).

## Step 1 — title and endpoint, verified before any chemistry

Verified from two independent routes that agree head-to-tail, with no truncation and
no hidden qualifier:

| route | file |
|---|---|
| Europe PMC REST `resultType=core` | `source/europepmc_core_39664425.json` |
| Europe PMC JATS `fullTextXML` `<article-title>` | `source/europepmc_PMC11628822_fullTextXML.xml` |

**Endpoint discussed:** mRNA transfection / delivery — luminescent total flux in muscle
(p s⁻¹) and %tdTomato⁺ cells in liver, spleen and lung. Topically **in scope**.
**But the endpoint belongs to a different paper**; this article measures nothing.

## Evidence that there is no original data

| fact | value |
|---|---|
| `article-type` | `letter` |
| PMC `subj-group` heading | `Highlight` |
| `<abstract>` in JATS / `abstractText` in core | absent / absent |
| section titles present | Key words · Author contributions · Conflicts of interest · References · Acknowledgments |
| Methods / Results / Experimental / Data availability | **none** |
| author-contribution roles claimed | Writing – original draft · Conceptualization · Writing – review & editing (no investigation, methodology, formal analysis, validation, resources or data curation) |
| `<table-wrap>` in JATS | 0 |
| `<table>` in PMC HTML | 0 |
| `<supplementary-material>` | 0 |
| Europe PMC `hasSuppl` | `N` |
| `<fig>` / `<graphic>` | 1 / 1 |
| figure caption | ends *"Reprinted with the permission from Ref.* `<xref rid="bib1">` *. Copyright © 2024 Nature materials."* |
| `bib1` `pub-id pmid` | **38740955** |

## Acquisition

- Publisher route (recorded in `00_acquisition_raw.json`, pre-existing): the DOI
  redirected to `linkinghub.elsevier.com`, returning a "Redirecting" page with
  **0 assets**. No paywall or bot check was bypassed.
- Open-access route used instead: Europe PMC REST (`core`, `fullTextXML`) plus the
  single figure.
- Figure fetched over **two independent routes that are byte-identical**:
  - `europepmc.org/api/fulltextRepo?...&mimeType=image/jpeg` → HTTP 200, 548 920 B
  - `cdn.ncbi.nlm.nih.gov/pmc/blobs/3944/11628822/.../gr1.jpg` → HTTP 200, 548 920 B
  - md5 `8cca629044a5e7ffee0f9d91f7c320c6` on both.
- Everything reasoned over is archived in `source/`; nothing this folder asserts
  depends on the session scratchpad (`CLAUDE.md` §7, F-023).

## Tier 4b — vector geometry

Checked, and **unavailable for a structural reason**: the sole figure is a single
**flattened JPEG** (740×917 px served; EXIF records a 2960×3666 original), identical
from both routes. A flattened raster exposes no vector drawing items, so there is
nothing for `page.get_drawings()` to return. No PDF with a vector layer was obtained —
the publisher route returned 0 assets. **Reported as: zero vector items, one flattened
JPEG per figure.**

## Anchor posture — recorded by hand, as instructed

`provenance/structure_check.py`, `provenance/check_anchor_independence.py` and
`provenance/check_anchor_coverage.py` were **NOT run**. They write tracked files that
other agents read concurrently, and the brief forbids it.

- **C10 anchor: none, and none is possible.** This curation enumerates no structure, so
  there is nothing to anchor. Every structure in Figure 1 belongs to the licensed
  reprint of PMID 38740955; an anchor derived from it would anchor `BL_2024`, written
  into the wrong folder.
- `reference_structures.csv` exists with the exact required header and **one row with
  an empty `reference_SMILES`**, `expected_mismatch = no_anchor_possible`, so the
  *absence* of an anchor is visible rather than silent (`CLAUDE.md` §6).
- **No `make_reference.py`.** A generator here could only succeed; a code path that can
  only pass is the anti-pattern §6 names. See `DECISIONS.md` D-2.
- **Independence by provenance is moot** (nothing to be independent *of*), and is
  recorded as moot rather than claimed. No independence assertion is made anywhere in
  this folder — §6 records that C10's independence "was asserted twice and false twice",
  so an unearned claim is worse than none.

## Gate

`verify.py` (this folder only). Clean run **17/17, exit 0**, exit code read from the
process return code, not through a pipe. Eight seeded defects — inserted table,
article-type flipped to `research-article`, `bib1` PMID changed, reprint attribution
removed, `hasSuppl` flipped to `Y`, a `Methods` section added, figure bytes corrupted,
and a 1-row `06_rows.csv` — **each moved the agreement count to 16/17 and exited 1**.
Zero no-op seeds. Details and the control run in `07_qc.json`.

The gate also caught a defect in its own first version: a hand-transcribed prose
expectation failed because the source has U+00A0 after "Ref." and the "1" is an
`<xref>`, not a literal. It was replaced by the stronger *markup* assertion that the
caption's reference points at `bib1`.

## Scope discipline

Nothing outside `curation/songx_39664425/` was written or modified. `data/` was read
only, for the report-only corpus comparison in `02_corpus_overlap.json`. No merge into
`data/`. `curation/lib/` was neither imported nor edited — with zero structures and
zero digitised values there was no geometry, assembly or emit step to perform, and
calling those modules to produce a clean line over nothing would have been exactly the
"check that examines nothing" failure §6 warns about.
