# DECISIONS — lius_40500862 (Wu, Shi et al., JACS 2025, 147, 28665–28673)

`build_dataset.py` refuses a folder with an `OPEN` entry. Nothing here is
resolved by the pipeline: these are the places where the source is silent, or
where resolving it would change what a measurement means.

The paper is `D_si_only` — we hold `ja5c01389_si_001.pdf` and not the article.
Europe PMC reports `isOpenAccess=N`, Crossref lists no component and the ACS
standard licence, and there is no repository on GitHub or Zenodo. Every
decision below would be easier with the main text and several would be settled
by it outright.

---

### D-1
Are the spleen/liver ratio panels (Figure S5, 11 lipids; Figure S23A, 4 lipids)
corpus-emittable at all?

They are the only panels in the SI that individually measure more than four
members, so dropping them costs 11 of the 15 lipids their only value. But a
spleen/liver ratio is a **selectivity contrast between two organs**, not a
delivery magnitude, and `CLAUDE.md` §8 has already ruled that kind of quantity
non-poolable once, for `molea.csv`'s `mTP_Difference`. A model fed a ratio as
if it were a potency is fitting organ tropism, not delivery. Worse, the two
ratio panels **disagree with each other** on both lipids they share, by 2.14×
and 0.69×, with the rank inverted (`08_contradictions.csv`, row 1).

Held, not emitted, pending this decision. Both panels are digitised and stored
in `05_values.csv` with `emitted_as_corpus_candidate=False` so the decision can
be made on the numbers rather than re-done from the figures.

STATUS: OPEN

---

### D-2
Figure S24 reports an in vivo `Avg Radiance [p/sec/cm²/sr]` for four lipids,
and **never says which organ**.

The caption is "In vivo mRNA delivery comparison between carbonate-bearing LNPs
and their carboxylate-containing counterparts"; the axis says only "Avg
Radiance"; the Methods (SI p.11) describe whole-animal imaging followed by
organ harvest and imaging, so both whole-body and per-organ radiance exist in
this study. `Model_target` cannot be guessed. The adjacent panel S23B is a
single harvested-organ image, which hints at per-organ, and the paper's whole
subject is spleen tropism, which hints at spleen — but a hint is not a reading,
and assigning the wrong organ would put four measured values on the wrong
endpoint.

These four rows are digitised and held. Main-text Figure 3H, which S23B says it
is taken from, would settle it.

STATUS: OPEN

---

### D-3
Does a four-row in-vitro block deserve an `Experiment_ID`?

The only panel that is unambiguously emittable is Figure S22C: 4 bars of
IGROV1 luciferase RLU. `Experiment_value` is a per-experiment z-score of
log10 (`CLAUDE.md` §3), and a z-score over four points — two of which are
**tied at the pixel level** (see D-4) — has two effective degrees of freedom.
It would enter the corpus as a 43rd in-vitro experiment, i.e. a full extra
source key, carrying four rows.

`SL_2025` is free: the corpus holds `SL_2020` and `SL_2021` across 54 distinct
IDs and no `SL_2025`. The question is not the key, it is whether a 4-row source
is worth the source key it costs, given that source count is the denominator of
the chance baseline in `CLAUDE.md` §7 and that this project's central claim is
about source separability.

STATUS: OPEN

---

### D-4
Figure S22C cannot order `A5-CE-C7-6` against `A5-C7-6`.

Both bars' top rules are detected on the same image row (230), with full-width
coverage 0.962 and 0.963 and left/right flank delta 0.0 px on each. At 69.004
px per decade the panel resolves no better than 3.4%, so the two stored values
are identical by construction (1.487 × 10⁶ RLU each). That is a real property of
the panel, not a detector failure — `bars.bar_top` refused nothing here.

It matters because this pair is the paper's own carbonate-versus-carboxylate
comparison: in vitro the panel says they are equal, while in vivo (S24 left)
the carbonate is 2.56× higher and marked significant. Emitting two tied values
into a 4-point z-score makes half the block a tie. The alternative is to emit
one and flag the other, or to emit neither.

The tie is written onto both rows' `curation_note` by `emit.emit_rows`.

STATUS: OPEN

---

### D-5
Five lipids' tail counts are inferred, not read.

`A5-CE-C6`, `A5-CE-C6-5`, `A5-CE-C8-7`, `A5-CE-C8-8` and `A5-CE-C10-6` are
named only as x-axis labels in Figure S5. The SI gives them no HRMS, no NMR and
no drawn structure, and `reference_structures.csv` stores an empty
`reference_SMILES` for each so the absence is visible. Their headgroup, their
carbonate linkage and their alcohol identity are all pinned — the `Cm-n`
convention is fixed by three HRMS anchors *and*, independently of mass, by the
SI naming 2-hexyloctanoic acid as `A5-C7-6`'s acid. What is assumed is only
**how many** tails: 3 for the linear member and 2 for the branched ones, by
analogy with the six siblings that do carry HRMS, against a general procedure
that says all A5 members got 3 equivalents (`08_contradictions.csv`, row 2).

They carry `confidence = medium` and `structure_assumed = 1`, derived from each
record's own anchor field, never a literal. They are only reachable through D-1,
since Figure S5 is their only value.

STATUS: OPEN

---

### D-6
The source licence is not CC-compatible.

Crossref lists three licence URLs for this DOI, all `content-version: stm-asf`:
`10.15223/policy-029`, `-037` and `-045`. These are ACS standard-licence
policies, not CC0 or CC-BY. `CLAUDE.md` §8 requires every constituent source to
be CC0/CC-BY-compatible before corpus release. Four digitised numbers from a
figure are plainly not the article, but the project has not yet written down
where that line sits, and this folder should not be the place it gets decided
silently.

STATUS: OPEN

---

### D-7
The 161-member library is not curatable from the SI, and may not be curatable
at all.

The SI contains no table of any kind — its own table of contents lists only
Experimental Procedures, Supplementary Figures and References — and all 37
figures are flattened JPEGs with zero vector content, so the per-member screen
values live in main-text Figures 1–3 only. Options are: acquire the article
through a legitimate institutional route and curate properly; keep the 15
SI-namable structures with the 4 emittable values; or drop the paper and mark
the tracker row so no future session re-attempts the same SI.

STATUS: OPEN

---

### D-8
Figure S22C's formulation is not stated on the panel.

The SI gives the initial in-vitro screening formulation as
`An-CE-Cm/DOPE/cholesterol/DMG-PEG2000 = 15/20/25/2` at a lipid/mRNA weight
ratio of 20/1 (pp. 6–7), and says molar ratios were "optimized" afterwards.
Figure S22C is a comparison of "representative" lipids and does not say which
formulation it used. The screening formulation is the likely one because the
assay (IGROV1, 25 ng/well, 24 h) matches the screening protocol exactly, but
likely is not stated, so the formulation columns are left blank rather than
filled with a plausible value.

STATUS: OPEN
