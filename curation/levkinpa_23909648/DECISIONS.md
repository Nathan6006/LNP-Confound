# DECISIONS — LL_2013 (Li & Levkin, Bioconjug Chem 2013, 24(9), 1543–1551; PMID 23909648)

Queued, not resolved. `build_dataset.py` refuses a folder carrying an `OPEN`
entry, which is the point: every one of these would change what enters `data/`,
and PIPELINE_V2.md §3 is explicit that the pipeline detects that a source is
silent and does not get to decide what the source meant.

This curation produced **29 emittable rows over 3 structures** and did **not**
merge. D-1 and D-2 are the two that block a merge outright.

---

### D-1
STATUS: OPEN

AMBIGUITY: **Nine of the twelve lipidoids that carry a value have no amine
structure anywhere in this SI.** The 20 amine head groups are referenced only
as the integers 1–20 in Table S1's `Amine` column. The SI names exactly one of
them — amine 1 = 1-(2-aminoethyl)pyrrolidine, on p.3 — and draws no amine
panel at all (grep over the full text layer: `Scheme` 0 hits, `structure` 0
hits; the only three molecular drawings, on pp. 7/10/13, are the three
purified amine-1 compounds).

The twelve valued lipidoids are `1C12`, `1C14`, `1C16` (amine 1, head known)
and `3C12`, `4C16`, `4C18`, `5C14`, `10C14`, `11C9`, `12C16`, `14C6`, `16C14`
(amines 3, 4, 5, 10, 11, 12, 14, 16 — unknown). Note that the paper's **best
crude lipidoid, `4C18` at 31.7 %**, is one of the unknowns.

Each unknown amine's neutral mass IS solvable from Table S1 (`01_blocks.json`
records them), but a **mass is not a structure** and this curation does not
turn one into the other.

CANDIDATES:
 (a) Acquire the main article, read the amine panel, re-run stages 01–06.
 (b) Merge only what is structurally determined (the 29 rows, 3 compounds).
 (c) Do not merge at all — see D-8.

EVIDENCE NEEDED: the main article's amine figure. Nothing else resolves it;
no further work on this SI can, because the information is not in it. Europe
PMC, 2026-11-10: `isOpenAccess=N`, `inEPMC=N`, `inPMC=N`, `pmcid=null`,
"Subscription required". `provenance/tools/find_repo.py`, 2026-11-10:
`repo_hits: []`, all queries completed, none rate-limited.

---

### D-2
STATUS: OPEN

AMBIGUITY: **The cell line is not stated anywhere in this SI**, and neither is
the reporter gene or the assay timepoint. Grep over the full text layer: `HEK`
0 hits, `HeLa` 0 hits, `293` 0 hits as a word (the four `293` matches are
digits inside mirrored NMR peak lists), `GFP` 0, `luciferase` 0. The endpoint
is printed as "the percentage of transfected cells" with the cargo identified
only as "plasmid DNA".

`Model_target` and `Model_type` are required corpus fields. They are written
as `NOT STATED IN THE SI` in `06_rows.csv` and were **not** inferred.

CANDIDATES: (a) read the main article's methods; (b) merge with both fields
blank and a `curation_note`.

EVIDENCE NEEDED: the main article. Do not infer the cell line from the
laboratory's other papers — the sibling folder `levkinpa_29558113` is the same
group and uses HEK293T, which is exactly the kind of plausible cross-paper
assumption that becomes a fact nobody can trace.

---

### D-3
STATUS: OPEN

AMBIGUITY: Table S6 reports four **binary blends** of two ionizable lipids
(`1C16-S:1C16-D` at 2:1, 1:1, 1:2 and 1:3). The best value on the panel —
**26.1 %, the 1:2 blend** — is a blend, not a single compound. The corpus
schema has one `IL_SMILES` per row.

CANDIDATES:
 (a) Drop the 4 blend rows (what this curation did: they are in
     `06_skipped.csv` with the reason on the row).
 (b) Emit each blend twice, once per component, with a mole fraction — this
     invents a value for each component that was never measured alone.
 (c) Extend the schema for multi-IL formulations.

EVIDENCE NEEDED: a project-level convention. This is not a question about
this paper. Note that `1C16-S` and `1C16-D` ARE measured alone on the same
panel, so option (b) would duplicate a label that already carries its own
measurement — the precise situation that made `AC_2025`'s first repair assign
419 measured values to the wrong molecules.

---

### D-4
STATUS: OPEN

AMBIGUITY: **Tables S3, S4 and S6 print no DNA dose.** Only Table S5 does
(75 / 150 / 200 ng).

There is a strong inference available and it is recorded as an inference, not
applied as a fact. Table S5's L2K control is **6.9 ± 1.7** at 75 ng,
**14.0 ± 1.8** at 150 ng and **22.4 ± 4.0** at 200 ng — three clearly distinct
values. Tables S3, S4 and S6 each print L2K = **14.0 ± 1.8**, matching the
150 ng block in *both* the mean and the SD, and matching neither other dose.
Gate `Q-7` measures this and fails if it ever stops holding.

CANDIDATES: (a) set 150 ng on those rows, citing the control; (b) leave the
dose blank, as this curation did (`dna_ng = "not printed (see Q-7)"`).

EVIDENCE NEEDED: the main article's figure legends for Figures 2A, 3 and 8.
Note that the inference, while tight, is an argument from a shared control
value, and a lab that re-uses one control measurement across figures would
produce exactly this pattern without the doses being equal.

---

### D-5
STATUS: OPEN

AMBIGUITY: **No standardisation has been applied, and one value is exactly
zero.** The corpus convention is that `Experiment_value` is a per-experiment
z-score of log10. Two problems:

1. `log10(0)` is undefined. Table S5's `1C16-S:DOPE=1:2` at 150 ng is printed
   `0` with SD `0` — a non-detect. This curation flagged it
   `value_is_bound = upper` (the true efficiency is at or below the detection
   limit, not exactly zero) and it is the folder's single caveated row.
2. The 29 rows span **3 DNA doses × 3 DOPE ratios**. Pooling them into one
   per-experiment z-score mixes formulations of the same three compounds;
   splitting them is error class 22 (deleting a between-condition offset).
   This is the `WT_2026` / F-016 question in a smaller form.

`value_normalization_unit` and `Experiment_value_within_model` were left
**unassessed** rather than filled in, because writing a unit would assert a
standardisation that was never performed.

CANDIDATES: (a) one z over all 29 rows, `value_normalization_unit =
experiment`; (b) z within dose; (c) substitute a floor for the zero;
(d) drop the zero row.

EVIDENCE NEEDED: a project-level convention, plus a decision on the zero.

---

### D-6
STATUS: OPEN

AMBIGUITY: **Tables S5 and S6 disagree on the same two compounds** —
`1C16-S` is 2.0 % in Table S6 but 0.1 / 0.4 / 0 % in Table S5's 150 ng block
(**5.0×** against the best of the three, unbounded against the 1:2 cell), and
`1C16-D` is 11.3 % against 7.0 % (**1.61×**). Measured by gate `Q-8`; see
`08_contradictions.csv` C-5.

Table S6 prints **no lipidoid:DOPE ratio**, so the two panels may describe
different formulations — but nothing in the SI says so.

CANDIDATES: (a) keep both, each tagged with its own `source_table` (what this
curation did; nothing is averaged); (b) keep only Table S5, which is the one
panel that fully specifies its formulation; (c) resolve from the main
article's Figure 8 legend.

EVIDENCE NEEDED: the main article's Figure 8 legend. Until then, pooling the
two panels into one experiment silently averages a 5× disagreement.

---

### D-7
STATUS: OPEN

AMBIGUITY: **Licence.** ACS Bioconjugate Chemistry 2013 research article;
Europe PMC `isOpenAccess=N`. The SI carries "Copyright © 2013 American
Chemical Society" and is **not** CC0/CC-BY. CLAUDE.md §8 requires every
constituent source to be CC0/CC-BY-compatible before corpus release.

CANDIDATES: (a) the project's position is that extracted numeric facts are
not a derivative of the copyrighted PDF; (b) exclude on licence grounds.

EVIDENCE NEEDED: a project-level position. Identical to the sibling
`levkinpa_29558113`'s open licence question; it should be settled once, for
the corpus, not per paper.

---

### D-8
STATUS: OPEN

AMBIGUITY: **Is a 3-compound, 1-head-group experiment worth adding at all?**
If D-1 is resolved by option (b), `LL_2013` enters the corpus as 29 rows over
three structures that share one head group and differ only in tail length
(C12 vs C16) and tail count (1 vs 2).

Measured against the live corpus (`04_corpus_overlap.json`, 29,195 vitro +
1,803 vivo rows, 19,919 unique `IL_SMILES`): **0 of 3 match by full InChIKey,
0 by skeleton block, and there are 0 `IL_name` collisions.** So this
experiment would contribute **nothing** to Result 3 (shared-lipid
concordance), which is the load-bearing result — it adds one more
near-disjoint silo.

It is not nothing, though: the amine-1 head group itself appears on **638
corpus rows across 5 experiments** (`KZ_2016`, `YX_2024`, `YR_2024`,
`XH_2025`, `YX_2026_SDL`), always on far more elaborate ester/amide lipids.
So these three compounds are the *unelaborated* version of a head group the
corpus already spans — which is a genuine, if small, point of chemical
contact.

CANDIDATES: (a) merge as `LL_2013`; (b) hold until D-1 is resolved and merge
all 12 lipidoids at once; (c) do not merge.

EVIDENCE NEEDED: a project call on whether a zero-overlap 3-lipid experiment
earns its place. This is the same question `CLAUDE.md` §8 asks about
`molea.csv`, and it was answered there by keeping the data separate.
