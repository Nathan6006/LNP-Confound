# Autonomous extraction run — plan of record

Opened 2026-11-09. **This file is living**: update it as findings change the
method, and record every change in `provenance/RUN_LOG.md`.

The instruction: gather the PDFs, extract everything worth extracting, run
thorough review before anything merges, keep a detailed log, continue across
session limits until told to stop or until nothing worth collecting remains.

---

## 0. THE RULE THAT OVERRIDES EVERYTHING ELSE

**Launch, do not announce. A sentence is not an action.**

On 2026-11-10 a turn ended with "Launching batch 2 on the three largest
remaining chemistry screens" and nothing was launched. The user had to notice.
In an autonomous run that is the worst available failure mode: the session
looks busy, the log reads correctly, and no work is happening.

So, every turn, without exception:

1. **Do the thing in the same turn you say it.** If a batch is to be launched,
   the `Agent` calls go out BEFORE the sentence describing them is written.
   Never write "next I will…" as the last thing in a turn.
2. **End every turn by running the guard:**
   `python3 -I <scratchpad>/tools/check_state.py`
   It fails if any paper is marked `in_progress_agent` with no live task file,
   if any folder is marked `extracted` without `05_values.csv`, or — the case
   that caused this rule — **if zero agents are running while papers are still
   `ready_to_extract`**. An idle run is a failed run.
3. **The queue is the source of truth, not the prose.** Update
   `run_queue.json` in the same turn as the action it describes.
4. **Never stop on a question.** An unresolved question goes into that folder's
   `DECISIONS.md` and the run moves to the next paper. The only legitimate
   stops are in §7.

## 1. Scope and priority

Source of candidates: `lit_review/data.json` (104 rows, 20 verbatim PubMed
queries behind it, provenance in `lit_review/provenance/`). Plus new searches
for anything the 2026-10-25 sweep missed.

**41 rows have `status = "new"`** — not yet curated. Their composition:

| kind | meaning | count |
|---|---|---|
| `ivl` | **ionizable-lipid chemistry screen** | **26** |
| `ivv` | in-vivo screen | 9 |
| `form` | formulation screen | 5 |
| `comp` | computational / compilation | 1 |

**Priority order, as instructed:**

1. **`ivl` first, always.** The corpus already has many formulation screens;
   what it is short of is lipid chemistry. An `ivl` paper with 150 lipids
   outranks a `form` paper with 1,944 formulations.
2. Within a kind, **largest stated library first** — more rows per unit of
   verification effort, and the big ones anchor the corpus.
3. `ivv` next, `form` last. Still collect them; just not first.

**Stop condition:** a source offering **fewer than 100 usable data points** is
not worth the verification cost. Record it as `skipped_too_small` with the
count and why, so the decision is visible and not repeated.

**The `lib = -1` problem.** 15 of the 26 `ivl` rows have no library size
recorded. Those are *unknown*, not small — library size is the first thing
stage 00 establishes, and a `-1` that turns out to be 300 lipids is exactly
what this run exists to find. Do not deprioritise them on a missing number.

## 2. Phase 1 — acquire everything NOW, while the login is live

The institutional session can expire at any time, so **acquisition is front-
loaded and completed before any extraction begins.**

Per-publisher routes are measured and recorded in `BROWSER_ACQUISITION.md`.
The mix across the 41:

| publisher | n | route |
|---|---|---|
| ACS | 19 | browser, no login (may need a 2nd navigation) |
| Wiley | 9 | browser, no login |
| Elsevier | 4 | **browser + the live institutional session** |
| Springer Nature | 4 | plain request, no browser needed |
| PNAS | 3 | browser, no login |
| RSC | 1 | browser, Cloudflare self-resolves in ~8 s |
| Taylor & Francis | 1 | untested |

30 of the 41 are graded `C` ("not reachable — try CWRU library access"), but
**that grade was assigned by a plain HTTP probe on 2026-10-25, before the
browser existed.** Today's matrix says a plain probe 403s on six of eight
publishers regardless of entitlement, so grade C is expected to badly
understate what is reachable. Re-probe, do not trust the grade.

**For each paper, acquire:**
- the article PDF;
- every supplementary file (PDF, xlsx, docx, zip);
- **publisher figure rasters** where they exist — `ars.els-cdn.com/content/image/1-s2.0-<PII>-gr<N>_lrg.jpg` for Elsevier, `media.springernature.com` for Nature. Measured at ~1.4x the linear resolution of a 600 dpi page render, open-access, and they remove a whole class of rasterisation loss.

**Layout:** one folder per paper under `curation/<lastname>_<pmid>/source/`,
with `00_acquisition.json` recording URL, sha256, bytes, content type, UTC
timestamp and **the licence string from the landing page** (§8 of `CLAUDE.md`
requires CC0/CC-BY compatibility before release; far cheaper to capture now).

**Rate discipline:** one DOI at a time, sequential, no crawling, no citation-
graph following. Automated bulk retrieval is what gets an institution's whole
IP range blocked, and the subscription is the university's.

**Downloads are untrusted data**: scripts live outside the download folders,
paths are passed as arguments, readers run under `python3 -I`.

## 3. Phase 2 — the extraction loop

Per paper, following `PIPELINE_V2.md` §3. The pipeline is file-based and
resumable, which is what makes this survivable across session limits.

```
stage 00   render + pdftotext; establish the STATED library size
stage 00b  00_value_sources.json  <- THE INTAKE GATE, before any digitiser
stage 01-03  structures, via curation/lib (chains, assemble, reagents)
stage 04   values, from the tier 00_value_sources.json names
stage 05   formulation / metadata
stage 06   the paper's OWN cross-checks, measured by script
stage 07   emit, via emit.emit_rows  (refuses a literal confidence)
stage 08   08_contradictions.csv
           make_reference.py -> the C10 anchor, then C16
           DECISIONS.md for anything the source does not settle
```

**The intake gate is the single highest-value step.** The two curations that
took values from a publisher Source Data file and an author repository have
**zero** digitisation defects between them; the four that digitised produced
eight. Exhaust tiers 1–4 before writing a digitiser, and say in writing why
each failed.

**Nothing merges here.** Stage 07 emits into the folder. The merge happens
only at a batch gate (§4).

## 4. Phase 3 — the batch review gate

**Do not merge on collection.** Papers accumulate in a *staging set*; the gate
runs over the whole set; only then does `build_dataset.py` merge.

**Gate trigger — whichever comes first:**
- **3 papers** extracted and self-checked, or
- **~1,000 new rows** staged, or
- a paper whose extraction raised a structural question.

That interval is not arbitrary: `CLAUDE.md` §6 Tier 4 already sets corpus-wide
review at every 3–5 merged papers, because per-paper review structurally cannot
see cross-paper drift.

**What the gate runs, in order:**

1. **Per-paper Tier-2 adversarial review** — one briefed reviewer per paper,
   read-only, re-deriving structures and values from the PDF, not from our
   scripts. **Opus**, per `MODEL_EVAL.md`: this is the role where the measured
   gap was 7 extra findings versus 1.
2. **Cross-paper checks** — C17 (one reagent name, one structure, now across
   the new papers too), duplicate SMILES against the live corpus, experiment-ID
   collisions, normalisation-unit consistency.
3. **I adjudicate.** Reviewer findings are claims: of the Siegwart reviewer's
   12, three were real, one half-right, one flatly wrong. Reproduce the evidence
   before acting. And **validate the repair, not just the finding** — F-051's
   finding was real and its proposed repair was wrong in a way that would have
   emitted 0 of 48 rows.
4. **Which half to repair** goes to `DECISIONS.md`, never resolved alone —
   `AC_2025`'s confident, internally consistent first repair would have put 419
   measured values on the wrong molecules.
5. **Merge**, then the full gate suite:
   `build_dataset.py` → `audit.py` → C15 → C16 → C17 → C18 → `test_lib.py`,
   and `test_checks.py` at least once per batch.
6. Update `CLAUDE.md` §3 counts, `docs_claims.json`, `CHANGELOG.md`,
   `lit_review/build_data.py`'s `ADDED` registry, and `RUN_LOG.md`.

A merge is blocked by: any `OPEN` entry in a folder's `DECISIONS.md`; a missing
or copy-suspect C10 anchor; any red gate.

## 5. Working within the 5-hour usage limit

The constraint is real; the answer is **state on disk, not state in context.**

**Resumability.** `provenance/run_queue.json` is the single machine-readable
state: per paper, its stage, its status, its blockers. Any session reads that
file and continues. The per-paper stage files already make each stage
independently re-runnable — that is why this design survives a hard stop
mid-paper.

**Context discipline.**
- One paper at a time in the main thread; never hold two papers' figures.
- Read QC JSONs and summaries, not whole CSVs. Never cat a 1,000-row file.
- Push bulk reading into subagents, whose tool output stays out of my context.
- Write conclusions to disk as they are reached, not at the end of a phase.

**Model choice** (from `MODEL_EVAL.md`, measured rather than assumed):
- **Opus** — Tier-2 adversarial re-derivation, batch-gate adjudication, the
  which-half-to-repair call, any figure where the paper contradicts itself.
- **Sonnet 5.5** — stage 00 rendering and text extraction, acquisition,
  tallying, docs sweeps, mechanical diff review, checking a *specific* named
  suspicion. It matched Opus 3/3 at finding a known defect with the right
  mechanism; it was 1-versus-7 at finding what nobody planted. So: never on
  Wave 1, freely on everything mechanical.
- **Never** run a check as an agent when it is a script. `test_checks.py` is
  not a reviewer.

**The shared browser belongs to the MAIN SESSION only.** There is one Playwright
instance and the main session drives it; a subagent navigating it lands in the
middle of whatever sequence the main session is running. An agent reached for it
on 2026-11-10 and was stopped. Agents get plain clients instead, which need no
login and no shared state:
- `scratchpad/tools/pmc_fetch.py <PMCID> <outdir>` — JATS full text and captions
- PMC figure rasters at `pmc.ncbi.nlm.nih.gov/articles/instance/<id>/bin/<name>`,
  whose names are the `xlink:href` values in that XML — a sibling agent pulled
  six main-text figures at 300 dpi this way
- `static-content.springer.com/esm/...` for Springer Source Data
- `ars.els-cdn.com/content/image/1-s2.0-<PII>-...` for Elsevier assets
If an agent genuinely needs a browser-only fetch, it reports that and the main
session does it.

**Resume, do not relaunch.** When a usage limit kills agents mid-task their work
is saved; `SendMessage` to the agent id continues it with its context intact.
All six were resumed this way on 2026-11-10 rather than restarted, and one was
nearly complete.

**Subagents for concurrency.** Extraction of *different papers* is genuinely
parallel. Rules:
- each subagent owns exactly one curation folder and writes only inside it;
- it never merges, never edits shared code, never touches `data/`;
- it returns a structured summary; its findings are **claims** until I
  reproduce them;
- reviewers are **read-only** and separate from extractors — a reviewer that
  edits destroys the independent second opinion on its own change;
- **capture every report to disk the moment it lands.** In round 2, four of
  five Wave-1 reviewers died when the CLI exited and all five were blocked from
  writing their own reports, so the only copies were in context.
- cap concurrency at 3 extractors so a bad shared assumption cannot land in six
  folders before anyone looks.

**On hitting a limit:** stop cleanly — write `run_queue.json` and `RUN_LOG.md`
first, mid-paper is fine — and resume from the queue.

## 6. Logging

- **`provenance/RUN_LOG.md`** — append-only, every action with a UTC timestamp:
  what was acquired, what was extracted, what a check said, what was decided,
  what was merged, what was skipped and why. Every change, as instructed.
- **`provenance/run_queue.json`** — machine state.
- Per-folder: `00_acquisition.json`, `00_value_sources.json`, `07_qc.json`,
  `08_contradictions.csv`, `DECISIONS.md`, `PROVENANCE.md`.
- Findings continue in `provenance/review/FINDINGS_<date>.md`; error classes in
  `ERROR_CLASSES.md`, taking `max + 1` and checking it is unused.

## 7. Stop conditions

Stop and report when any of these holds:
- the queue is empty of sources with ≥100 usable points;
- the login expires and the remaining queue is Elsevier-only (**ask for a new
  login**);
- a finding implicates already-merged corpus data (stop, do not merge more);
- a check goes red for a reason not understood — never edit a number to make a
  gate green.

## 8. How this file changes

When something is learned that changes the method — a publisher route, a figure
source, a check that should exist, a reviewer brief that worked — edit the
relevant section here **and** append the change to `RUN_LOG.md` with its
evidence. The point of the run is a better pipeline as much as more rows.
