# Policy decisions — corpus admissibility and release

These are the project-level calls that close OPEN entries in many
`curation/*/DECISIONS.md` files at once. A per-folder decision that is
*an instance* of one of these should cite the P-number and close, not
re-argue it.

Opened 2026-11-10, when a merge-gate audit found **280 OPEN decisions across
46 folders**, of which 60 (in the 24 merge-candidate folders) were restatements
of six policy questions.

---

## P1 — Source licence is NOT a merge gate, and not a release gate for values

**DECIDED 2026-11-10. Supersedes the `CLAUDE.md` §8 line as written.**

§8 said: *"constituent-source licenses must be CC0/CC-BY-compatible before
corpus release."* As written this is wrong in a way that would have deleted the
corpus's own foundation. Evidence, measured not assumed:

**1. Our own baseline already fails it.** `data/Original_LNPDB.csv` is LNPDB
(Collins et al., *Nat Commun* 2026, PMID 41605942, DOI 10.1038/s41467-026-68818-1).
Europe PMC licence fields for its 42 source publications, resolved 2026-11-10:

| licence | n |
|---|---|
| none recorded (closed access) | 27 |
| CC BY-NC-ND | 7 |
| CC BY | 6 |
| CC BY-NC | 1 |
| unresolved | 1 |

**28 of 41 are `isOpenAccess: N`.** Only **6** are CC BY. Applying §8 literally
to our own inherited rows would remove ~85% of the baseline corpus — including
`AA_2008` (Akinc, *Nat Biotechnol*, closed) and `KW_2014`, which together are
57% of Result 3's pooled concordance points.

**2. LNPDB did exactly what we do, and said so in print.** Its Methods:

> "Functional data—most commonly delivery performance—were extracted from
> published heatmaps and bar plots by **digitizing the figures** and
> interpolating values based on either the color scale (heatmaps) or bar height
> (bar plots) as defined in the accompanying legends."

and

> "publications were selected from the literature based on the presence of large
> screening datasets"

It then released the result publicly and downloadable at `lnpdb.molcube.com`,
in a **CC BY 4.0** *Nature Communications* paper from MIT (Anderson, Langer) and
Lehigh, with a Data availability statement that makes no mention of source-paper
permission, and **no copyright or permissions discussion anywhere in the
article**. That is peer-reviewed, editor-passed precedent for this exact act on
this exact kind of source.

**3. The legal basis is that measurements are facts.** Under *Feist Publications
v. Rural Telephone Service*, 499 U.S. 340 (1991), facts are not copyrightable
because they owe their origin to no act of authorship. A compilation is
protected only in its **selection and arrangement**, and that protection does not
reach the underlying facts. A transfection value is a measurement of the world.
Publishers' own guidance agrees: Wiley and Elsevier both state that a new
figure or table built from factual data taken from a copyrighted work does not
require permission.

**4. Field precedent is overwhelming.** ChEMBL manually extracts bioactivity
data from closed journals (*J. Med. Chem.* and others) and releases it CC BY-SA;
BindingDB releases its own curation CC BY 4.0. Both are long-running,
institutionally hosted (EMBL-EBI, UCSD), and neither claims publisher permission.

### What P1 therefore permits and forbids

| | status |
|---|---|
| **Merging extracted numeric values into `data/`** | **ALLOWED** regardless of source licence. This is the merge gate, and it is now open. |
| **Merging structures (SMILES) we derived** | **ALLOWED.** A chemical structure disclosed in a paper is a fact about a molecule. |
| **Publicly releasing the corpus CSV** | **ALLOWED**, with full per-row citation to the source publication. Attribution is the price and we pay it anyway: every row carries `Publication_PMID`. |
| **Redistributing source PDFs, SI files, publisher Source Data files** | **FORBIDDEN.** This is verbatim copying of the protected expression. Already enforced — 110 files / 156.9 MB were removed from the repo on 2026-11-10. |
| **Publishing rendered figure images, crops, or digitisation overlays** | **FORBIDDEN.** An overlay is a derivative of the publisher's figure, which *is* protected expression. **31 such images from earlier commits are still public in git history — see OPEN ITEM below.** |
| **Reproducing a source figure in our manuscript** | Needs permission or an STM-guidelines check, per journal. Not needed for our six results, which use none. |

### Caveats that remain live

- **Contract, not copyright, is the real residual risk.** CWRU's subscription
  agreements with Elsevier/Wiley/ACS may contain terms restricting bulk text and
  data mining independently of copyright law. Our acquisition was
  article-at-a-time manual reading, not bulk API harvesting, which is the
  posture those clauses target. If a reviewer or the library asks, that
  distinction is the answer.
- **EU sui generis database right** (Directive 96/9/EC) protects substantial
  investment in a *database*, and could in principle attach to LNPDB itself
  rather than to the papers. It does not bite here: LNPDB is CC BY, which
  licenses exactly this reuse, and our extractions from individual papers are
  not extractions from any single database.
- **This is a reasoned position, not legal advice.** Before the dataset is
  released publicly, it is worth one email to the CWRU library's scholarly
  communications officer stating the plan in three sentences. That is cheap and
  it converts an inference into a sign-off.

### OPEN ITEM — git history

31 digitisation-overlay images from earlier commits are still reachable in the
public history of `github.com/Nathan6006/LNP-Confound`. Removing them requires a
history rewrite, which is the user's call and has not been made.

### Action taken

`CLAUDE.md` §8 to be rewritten from "constituent-source licenses must be
CC0/CC-BY-compatible" to the table above. **~24 folders' licence decisions close
as instances of P1**, citing this section.

---

## P2 — One-ionizable-lipid screens are ADMITTED, with a `design_axis` flag

**DECIDED 2026-11-10.**

The question: a paper that holds the ionizable lipid fixed and varies only the
formulation (helper lipid, cholesterol, PEG, ratios, dose) measures something
real, but it is not a *chemistry* screen. Does it belong in a corpus built to
study cross-laboratory *chemistry* transfer?

**Decision: admit, and label.** Reasons:

1. **Precedent already exists and is load-bearing.** `AC_2025` is a 498-row
   one-lipid screen already in `data/`. Excluding this class now means either
   removing `AC_2025` or defending an inconsistency in Methods.
2. **Exclusion is unrecoverable; labelling is reversible.** A folder refused at
   the gate costs the whole extraction to re-do. A flagged row costs one
   `WHERE` clause.
3. **Result 6 needs them.** The within-experiment learnable regime is precisely
   where formulation variance lives.
4. **Result 3 must not see them.** A screen with one ionizable lipid contributes
   no cross-lab shared *chemistry* — only a benchmark point, which is the F-064
   hazard.

**Implementation:** a new corpus column `design_axis`, one of
`chemistry` (the ionizable lipid varies) / `formulation` (it does not) /
`mixed`. Default `chemistry`. Stage 5 (Result 3) filters to `chemistry`.
This is a schema change and so must go through `provenance/schema_quality.json`,
`audit.py` C13's vocabulary check, and a seeded-failure test per §6.

**Closes:** the "does a one-lipid screen belong" decision in ~11 folders.

---

## P3 — Benchmark-only sources: DEFERRED to the analysis phase

**DEFERRED 2026-11-10 at the user's direction.** The run is in data collection;
this is a question about how Result 3 is computed, and it should be answered
when the concordance is being computed and its sensitivity can be measured
rather than argued.

The question: should a source whose *only* chemical overlap with the rest of the
corpus is the field's universal positive controls (DLin-MC3-DMA, SM-102,
ALC-0315, DOTAP) contribute points to the cross-laboratory shared-lipid
concordance in Result 3?

Why it matters: those benchmarks top nearly every screen that includes them, by
construction — they are chosen as controls *because* they work. A pair of
sources sharing only benchmarks can therefore show high rank agreement for a
reason that has nothing to do with transferable structure–activity learning.
Because Result 3 is the load-bearing result and currently sits at pooled
Spearman **0.083**, a cluster of benchmark-only pairs could move it materially
in the optimistic direction, for a spurious reason.

**What to do meanwhile:** curate and merge these sources normally (P1 and P2
permit it); do **not** decide their concordance status. Flag them so the
question is answerable later.

**Carry forward into the analysis phase as a measured sensitivity, not a
judgement call:** compute Result 3 with and without benchmark-only pairs and
report both. This already has a precedent in the project — `CLAUDE.md` §8
records the same worry for `XD_2026`, whose only shared lipids are MC3 and
SM-102 ranking 2nd and 3rd of 12.

**Affects ~13 folders.** Their benchmark-overlap decisions should be marked
RESOLVED-BY-P3 with the pointer here, not left OPEN, so they stop blocking the
merge gate.

---

## P4, P5, P6 — see the session record of 2026-11-10

P4 (`Experiment_ID` convention), P5 (schema vocabulary extension) and P6
(admissible endpoint classes) were put to the user on 2026-11-10 and the user
asked for them to be restated in plainer terms before deciding. Pending.
