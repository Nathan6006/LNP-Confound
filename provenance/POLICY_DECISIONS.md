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

## P4 — `Experiment_ID` is a SCREEN key; a new `lab` column carries the laboratory

**DECIDED 2026-11-10.** The user's words: *"it should be same screen, you can
also add a lab column. Cross laboratory should really be cross experiment."*

**(a) One tag = one screen.** Two papers from the same laboratory keep two
`Experiment_ID`s. A new column **`lab`** records the laboratory (the
corresponding author's group), so the grouping that Results 1/4 actually use can
be chosen at analysis time instead of being baked into the key.

**This renames the central claim.** The honest unit is **cross-experiment**, not
cross-laboratory: `Experiment_ID` has always been a screen key, and the
manuscript title and Results 1-5 have been describing it as a laboratory key.
With `lab` present, both are measurable and the weaker one can be reported as a
robustness check. **Flagged for the manuscript: the title
("Cross-laboratory prediction...") and every use of the phrase need to be
reconciled with what the key is.** Not a data change; a wording-and-claim change
that the user has now made deliberately.

**(b) Tag collisions resolve by date, suffixed `a`/`b`.** Two papers wanting
`YZ_2024`, and two wanting `SL_2025`: the earlier publication date takes `a`,
the later takes `b` (`YZ_2024a`, `YZ_2024b`). A bare colliding tag is never
emitted. `audit.py` must gain a uniqueness check, seeded to fail per section 6.

**(c) Multiple screening rounds inside one paper share one tag**, and the
blocking is documented. Where a paper screens in rounds (e.g.
`zhangz_40433897`: 60 formulations in Fig S4b, 18 in Fig S5), all rounds take
the same `Experiment_ID` and are z-scored together as one block. **Each such
case gets an entry in `provenance/EXPERIMENT_BLOCKING.md`** naming the paper,
the rounds, the figures, and why they are one block -- because the z-score
convention makes blocking change every value in the folder, and a pooling
choice that lives only in a curator's head is unauditable.

**Closes:** the `Experiment_ID` decision in ~17 folders and the
how-many-experiments decision in 8.

---

## P5 — Vocabulary may grow, but only by hand, and never by merging two things

**DECIDED 2026-11-10.** The user's words: *"if it is a genuinely new value, it
can take a new name, but if they are really just describing the same thing, they
should be grouped together. But be extremely careful you don't mis-group
different things together. I think before just having a script add some new
topic, you should review it and think about whether it deserves its spot."*

The rule, in order of what it forbids:

1. **No script adds a vocabulary value.** A curation may *propose* one; a human
   (or this session, reasoning explicitly in the folder's `DECISIONS.md`) decides
   whether it earns a slot, and the reasoning is recorded.
2. **Mis-grouping is the worse error and is treated as such.** Collapsing two
   genuinely different readouts into one token destroys information that no
   later analysis can recover, whereas splitting one readout into two synonyms
   is visible and repairable. **When in doubt, split.** The asymmetry is the
   point: the corpus has already been through the mirror-image failure, where
   `"A1I5R22C8 "` and `"A1I5R22C8"` were one lipid stored as two keys for the
   corpus's whole life (F-022, 915 cells).
3. **Synonyms do get merged** -- `GFP`, `GFP_fluorescence`, `gfp_normalized`
   naming one assay is three tokens for one thing and gets one.
4. `audit.py` C13 **fails the build** on any value outside the approved list, so
   a quiet addition cannot reach `data/`.

**One addition is already forced:** F-080's repair needs
**`experiment_x_model_target`** in `value_normalization_unit`. The column today
has no way to say "standardised per organ", which is part of why 378 rows
(`ZH_2023` 325, `BL_2023` 53) are mislabelled `experiment` right now.

---

## P6 — Admissible endpoints: delivery means the cargo acted

**DECIDED 2026-11-10.**

**Admitted**, because each measures whether the cargo reached the cytosol and
did its job:

| endpoint | note |
|---|---|
| luminescence (FLuc/RLuc) | the corpus's 24,908-row backbone |
| protein abundance (hEPO, FVII) | 766 rows |
| gene-editing efficiency | 141 rows, `SX_2025` -- user's explicit call |
| **knockdown** | user's explicit call. **Polarity verified, not assumed** -- see below |

**Excluded, and kept in a separate file** rather than deleted, at the user's
instruction (*"keep them in a separate file in case they are ever useful in the
future"*). Destination: **`data/excluded_endpoints.csv`**, same 79-column
schema, with an added `excluded_because` column:

| endpoint | n | why |
|---|---|---|
| `uptake` | 479 | measures where the particle went, not whether it worked |
| `diameter` | 96 | physicochemical |
| `zeta_potential` | 96 | physicochemical |
| TNS fluorescence (`whiteheadka_38888972`) | tbd | a dye assay for lipid pKa -- a property of the molecule, not a delivery outcome |
| antibacterial MIC (`pabst_37602460`) | tbd | not delivery under any reading |

The first three are the **671 `LX_2024` rows** that produced the `-585` in-vivo
delta; they were excluded by the a-priori rule of 2026-07-23 and have had no
home since. This gives them one, and makes the exclusion inspectable instead of
merely asserted.

### Knockdown polarity — MEASURED, and the corpus is already correct

A silencing readout runs backwards: less signal is better delivery. Pooling it
unflipped onto a z-scale where high = good would invert those rows. **The
corpus already handles this correctly, and here is the proof rather than the
assumption.**

`KW_2014` has 1,182 in-vitro luminescence rows (high = good) and 91 in-vivo
`protein_adundance_normalized` FVII rows. FVII silencing is measured as
*remaining* FVII protein, so if the sign had never been flipped, the two pools
would disagree. Over the **82 lipids they share**:

> **Spearman(KW_2014 vitro, KW_2014 vivo) = +0.402, p = 0.0002**

Positive, and significant. The flip was already applied upstream. So the
convention for new knockdown data is **orient so that higher = more silencing =
better delivery**, matching `LL_2012`'s `LRP6_knockdown_normalized` (112 rows)
and `KW_2014`'s FVII.

Two notes that fall out of this check:
- The column name `protein_adundance_normalized` is **doubly misleading**: it is
  misspelled (F-062, 766 rows) *and* on the `KW_2014` rows it does not hold a
  protein abundance -- it holds an inverted one. The name should say so.
- `Cargo` (mRNA 26,376 / siRNA 3,430 / pDNA 1,192) is the siRNA-bearing column.
  `Cargo_type` holds the *reporter* (FLuc, hEPO, FVII, GFP...). I searched
  `Cargo_type` for "siRNA" first and got 0, and briefly believed the corpus had
  no siRNA. It has 3,430 rows. Corrected here so the next reader does not
  repeat it.

### Ratio of means vs. mean of ratios (`dixonjb_41269090` D-4)

The user asked what this meant. A paper reporting "treated cells fluoresced 3.2x
the control" can compute that two ways: average all treated wells, average all
control wells, then divide (**ratio of means**); or divide each treated well by
its paired control and average those (**mean of ratios**). They are not equal --
the second weights low-control wells heavily and is the larger of the two
whenever the control varies. The paper's own statistic is a mean of ratios; the
corpus wants one number per formulation, which is naturally a ratio of means.
**Decision: use the ratio of means, record the paper's own figure beside it in
the curation's QC file, and note the choice in `curation_note`.** The gap is a
few percent and it is taken before a log and a z-score, so it moves ranks
negligibly; what matters is that the choice is written down rather than silent.

---

## Group decisions of 2026-11-10 (per-folder classes)

### Stereochemistry — MEASURED, and the answer is split

The user's belief was that LNPDB does not store stereochemistry. **Measured on
`Original_LNPDB.csv` (19,797 rows), it stores one kind and not the other:**

| | rows with it | share |
|---|---|---|
| **E/Z double-bond geometry** (`/`, `\`) | **7,648** | **38.6%** |
| **tetrahedral centres** (`@`) | **215** | **1.1%** |

The 38.6% is the *cis* unsaturation in oleyl/linoleyl tails, and it is
chemically load-bearing -- a *cis* tail packs differently from a *trans* one and
that is much of why unsaturated lipids work. **It must always be kept.**

The 1.1% appears in only 5 of 43 experiments (`KW_2014` 73, `JW_2024` 67,
`KL_2010` 54, `SY_2025` 11, `BL_2023` 10). **So the rule is: store E/Z always;
store a tetrahedral centre flat unless the source explicitly pins it**, which
matches the 98.9% convention and the user's intent. Closes the stereo decisions
in `dahlmanje_36701517` (cKK-E12), `peer_39912611` (2-hexyldecanoate),
`xuq_32520527` (pyrrolidine), `khanof_36652629` (six mixtures),
`maohq_39375194` (DOTAP).

### "No correct single SMILES" — recheck for intent first

The user: *"i think that is unlikely, recheck and try to find what the paper
intended."* Agreed and adopted: before any folder records a reagent as an
irreducible isomer mixture, the curation must try the supplier catalogue number,
the CAS, and the paper's own synthesis description. Only a reagent that is
*sold* as a mixture (isodecyl acrylate is the real instance) keeps that status,
and then it carries `structure_assumed = 1` with the branch named in
`curation_note`.

### Crude products — emit the intended product

The user: *"it is fine to call it the intended product, since that is what they
are trying to test anyways."* Adopted. A crude library member is emitted as the
intended structure with `positional_ambiguity = 1` (or
`structure_assumed = 1`), so the row states the convention instead of hiding it.
This matches how `NG_2025`'s 332 unpurified Mannich rows are already stored.

### Unstated metadata — one more look, then blank

The user: *"go with your recommendation of emitting blank with documentation.
However before falling back to this, go back and check for a final check to see
if you can find them."* Adopted as a two-step: each unstated field gets one
documented search pass (the paper's Methods, the SI, the author repo, the
supplier page, a sibling paper from the same lab) recorded in the folder's QC
file as *where we looked*; only then is it emitted blank. A blank with a record
of the search is different from a blank nobody tried to fill.

### Source self-contradictions — re-verify ours first, then prefer the primary

The user: *"for each one, carefully reanalyze whether the error you found is
real, and if it isn't, just fix it."* Adopted, in this order:

1. **Re-derive our own reading.** Several of these may be our digitisation, not
   the paper's error -- `zhangz_40433897`'s wrong-sign correlation is the
   likeliest candidate, and `xuq_32438474` D-2 is openly a choice between two
   defensible axis origins.
2. If ours is wrong, **fix it** and record the correction; it is not a finding.
3. If the paper really does contradict itself: **prefer the primary artifact**
   -- an author repository file or a publisher Source Data table beats a
   figure; a figure with printed numerals beats a colour scale. Where two
   artifacts of the *same* rank disagree, **emit nothing for those cells** and
   record the contradiction in `08_contradictions.csv`.
4. Every such call is written into `DECISIONS.md` with the evidence, because
   this is adjudication, not measurement.

### Group 4 scope calls — DECIDED, with the chemistry checked

The user asked for a reasoned call and guessed "no, yes, no, yes". **All four
agree with what the chemistry and the sizes say**, though for `mckinlaycj` the
reason is size, not chemistry:

| folder | call | reasoning |
|---|---|---|
| **`palczewski_41654520`** (64 rows, 30 structures) | **EXCLUDE** | Coomassie-derived: MW 1,270-1,942, 5 aromatic rings, **2 permanent sulfonates**, logP 16-26. An ionizable lipid works by being neutral at pH 7.4 and protonating in the endosome; a molecule carrying two permanent anionic charges cannot do that, so `num_protonatable_nitrogens` and `formal_net_charge` would be describing a different mechanism. Different chemical class, not a borderline case. |
| **`xuq_32438474`** (3 blocks) | **SPLIT** | `YL_2020` (36 GFP-mRNA values, 4 cell lines) is core corpus material -- admit. `YL_2020_2` (27 protein-delivery values) admit as its own `Experiment_ID` per the user's call. **`YL_2020_3` EXCLUDE**: it is Cre recombination, a downstream functional readout of the same delivery event as `YL_2020_2`, **Spearman 0.946** between them over the same 27 cells. Admitting both would put two near-copies of one ranking of 9 lipids into the concordance pool -- the double-counting that inflates Result 3. |
| **`mckinlaycj_39105751`** (5 rows) | **EXCLUDE, on size** | Chemically it is *defensible*: MW 1,019-1,064, 6-7 sp3 nitrogens, neutral -- closer to a polyamine lipidoid (cKK-E12, C12-200) than the dyes are. But **5 values cannot be z-scored meaningfully**, which is that folder's own D-2. Excluded for sample size; the chemistry objection would not have carried it. |
| **`chamolstad_41691524`** (77 rows, 1 structure) | **ADMIT** | A formulation DoE on one ionizable lipid -- exactly P2, so `design_axis = formulation` and Result 3 never sees it. |

Both exclusions go to `data/excluded_endpoints.csv` under P6's keep-don't-delete
rule, with `excluded_because` naming the reason above.
