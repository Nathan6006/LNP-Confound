# DECISIONS — li_39883839 (Xu et al., PNAS 2025, PMID 39883839)

The pipeline detects that a source is silent; it does not get to decide what the
source meant (PIPELINE_V2 §3). Every entry below is queued, not resolved.
`build_dataset.py` refuses a folder with an `OPEN` entry, which is the intended
state of this folder.

Settled questions are **not** here — they are in `08_contradictions.csv` with
the measurement that settled them.

---

### D-1 How may the three mLuc panels be blocked for z-scoring?

STATUS: OPEN

AMBIGUITY: Figures 2H, 2K and 3I each measure individual lipids as whole-animal
IVIS total flux at 0.05 mg/kg mLuc, I.M., 6 h. They overlap: H2A4B4 is in all
three, H2A4B3 in two. **Fig 2H and Fig 2K disagree by 25.8× on H2A4B3 and 2.79×
on H2A4B4, in the same direction** (`07_qc.json` V4). Fig 2K and Fig 3I agree to
**1.04×** on H2A4B4 (V2). A single `Experiment_ID` z-scored across all three
would mix two scales; three `Experiment_ID`s would give blocks of 4, 6 and 3 rows.

CANDIDATES:
1. One block over Fig 2K + Fig 3I only (9 rows, 8 distinct structures), with Fig
   2H dropped. Defensible: those two agree on their shared lipid to 4%, and Fig
   2H is the one panel that uses **crude, unpurified** Passerini product
   (Methods: "Except for the high-throughput screening and batch-testing stages,
   all other experiments utilized purified lipids").
2. Three separate `Experiment_ID`s (`YX_2025`, `YX_2025_2`, `YX_2025_3`), each
   z-scored alone. Honest, but a 3-row block carries almost no rank information.
3. One block over all 13 rows after tying Fig 2H to Fig 2K on their two shared
   lipids. **Explicitly flagged as the QL_2026 trap** (F-018/F-031): a two-anchor
   scalar offset computed from our own digitised values largely measures our own
   error, and the two anchors here disagree with each other by 9.2× on what the
   offset should be (25.79 vs 2.79), so there is no single offset.
4. Emit nothing.

EVIDENCE NEEDED: whether Fig 2H's four LNPs are crude product and Fig 2K's
purified — the Methods sentence implies it but never says so of these panels.
Failing that, a statement from the authors of the per-panel cohort. **Candidate 3
is not available on the data as it stands**, and that is a measurement, not an
opinion.

---

### D-2 Does a 9-structure, 13-row in-vivo set belong in the corpus at all?

STATUS: OPEN

AMBIGUITY: The paper's headline is a **144-member** library. That library has
**no per-member delivery value anywhere** — it was screened in pooled batches of
16 and then 4 LNPs per mouse. What is individually measurable is 9 distinct
structures, 2 of which are the benchmarks DLin-MC3-DMA and SM-102, and 7 of
which are one head group (H2 or its spacer homologs) on one tail A (A4).

CANDIDATES:
1. Merge. It adds a 2026 in-vivo I.M. muscle source and two benchmark rows.
2. Do not merge; keep the folder as a curated negative result — a documented
   case where a large published library yields no poolable per-member value.
3. Merge the 144 **structures** with no `Experiment_value`, as chemistry-only
   rows.

EVIDENCE NEEDED: a project-level call. Note the bearing on `CLAUDE.md` §8's open
question about benchmark-only contributions to Result 3: of the ≤9 structures,
MC3 and SM-102 are two, and SM-102 has only a bound (D-4).

---

### D-3 Should Figure 5's gene-editing and hEPO endpoints also be curated?

STATUS: OPEN

AMBIGUITY: Fig 5C/5D (% tdTomato+ by flow, 7 cell populations) and Fig 5F
(plasma hEPO, mIU/mL, n=8) measure SM-102, A4B4-S4 and A4B4-S3 by
**retro-orbital** injection to **liver** — a different route, organ and readout
from the I.M. mLuc panels, and only 3 lipids, 2 already covered by D-1.

CANDIDATES: (1) curate as a second `Experiment_ID` with `Model_target = liver`;
(2) leave out as a 3-lipid endpoint that cannot support a rank comparison;
(3) curate only the CD45−/CD31− panel, which the running text reproduces as
printed statistics (56% vs 48%; 34% vs 28%) and which is therefore checkable.

EVIDENCE NEEDED: whether the corpus wants a 3-row block on a third endpoint.

---

### D-4 How are the three bounded rows handled?

STATUS: OPEN

AMBIGUITY: Three of thirteen measurements are bounds, not point values, and the
corpus's `value_is_bound` column exists for exactly this.
- Fig 2H **H2A4B1** and **H2A4B2** are drawn with no resolvable height →
  `value_is_bound = upper`, bound 6.9e5 p/s (the panel's detection floor).
- Fig 2K **SM-102**'s bar top is above the axis break, and the upper segment
  carries **one** labelled tick (1e8), so no scale can be fitted to it →
  `value_is_bound = lower`, bound 4.0e7 p/s (the verified top of the calibrated
  lower segment).

CANDIDATES: (1) keep all three with their bound and let the modelling decide;
(2) drop the two upper-bounded rows and keep SM-102, since a lower bound on a
benchmark is more useful than an upper bound on a dead lipid; (3) drop all three.

EVIDENCE NEEDED: nothing from the paper. This is a corpus convention call. A
z-score of a censored value is not a measurement, and `CLAUDE.md` §3 already
warns that "a model treating these as measurements is fitting the floor of a
colour ramp".

---

### D-5 The licence is not CC0/CC-BY-compatible

STATUS: OPEN

AMBIGUITY: PNAS publishes this article under **CC BY-NC-ND 4.0**
(`00_acquisition.json`, from the JATS `<license>` element). `CLAUDE.md` §8
requires every constituent source to be CC0/CC-BY-compatible **before corpus
release**. Both the NC and the ND terms bite, and ND in particular is awkward for
a derived, re-normalised value.

CANDIDATES: (1) include and resolve at release time with the rest of the licence
audit; (2) include but mark the rows as non-redistributable; (3) seek permission;
(4) exclude.

EVIDENCE NEEDED: the project's position on whether a z-scored scalar extracted
from a figure is a "derivative" under ND. This is not a question this curation
can answer, and it applies to more sources than this one.

---

### D-6 Should the 144 library structures be recorded even though nothing is merged?

STATUS: OPEN

AMBIGUITY: `03_enumerated.csv` holds 148 named structures (147 distinct) that
reproduce 138/144 of Table S1's formulas exactly and explain the other 6 as
table defects (C-02, C-03). They have **zero** overlap with the corpus on full
InChIKey, skeleton InChIKey and `IL_name`, across both pools. That is a real,
checked asset whether or not any value is merged.

CANDIDATES: (1) keep in the folder only; (2) register in `lit_review`'s ADDED
registry as curated-but-not-merged so the sweep does not re-surface the paper;
(3) merge as chemistry-only rows (= D-2 candidate 3).

EVIDENCE NEEDED: whether `lit_review/build_data.py`'s `ADDED` registry is meant
to mean "merged" or "curated". This folder is not editing that file.
