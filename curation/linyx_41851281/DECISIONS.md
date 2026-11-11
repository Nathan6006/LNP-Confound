# DECISIONS — linyx_41851281

Lin et al., *Artificial intelligence-guided design of LNPs for in vivo targeted
mRNA delivery via analysis of the spatial conformation of ionisable lipids*,
Nature Biomedical Engineering 2026. PMID 41851281, DOI 10.1038/s41551-026-01640-8.
Proposed `Experiment_ID`: **YL_2026** (in vitro) / **YL_2026_2** (in vivo), to be
settled with the rest of the merge.

The pipeline detects that a source is silent; it does not get to decide what the
source meant. Everything below is queued, not resolved alone.

### D-1
**Which N-H carries which arm, on the four lipids where more than one answer
fits.** The Michael addition puts `n` arms on a head with `s` N-H sites. Where
`n < s` and more than one assignment gives a distinct molecule, the paper
neither draws nor states which. HRMS cannot see it — the formula is identical
either way, which is the Lin-regiochemistry lesson in `CLAUDE.md` §5 (and this
is, separately, a Lin paper).

`stage03_enumerate.py` enumerates every per-site assignment, builds each, and
counts distinct InChIKeys, so the flag is a measurement and not `n < s`. On
H4 three arms can only be 2+1 and 2+1 is 1+2, so a 3-tail H4 lipid is **not**
ambiguous and is not flagged. Exactly **four of 116** lipids are:

| ML ID | name | arms / sites | distinct structures |
|---|---|---|---|
| 7 | H4-OT1(2) | 2 / 4 | 2 |
| 10 | H4-NT1(2) | 2 / 4 | 2 |
| 96 | H4-NT12(2) | 2 / 4 | 2 |
| 106 | H3-NT12(2) | 2 / 3 | 2 |

Emitted under the convention **maximal dispersion** — one arm on each nitrogen
before any nitrogen takes a second — with `positional_ambiguity = 1`,
`confidence = medium` and a `curation_note` naming the choice. The alternative
is 2+0, i.e. one nitrogen fully substituted and one amine left free, which is a
different pKa and a different molecule.

Evidence that bears on it but does not settle it: Supplementary Fig. 12 draws
Lipid T3 and T4 (ML IDs 107 and 108, both H4) with the arms spread 2+1 and 2+2,
and the p.56 scheme separates the 3-tailed P1 from "4-tailed byproducts". Both
are cases where the distribution is forced anyway, so neither tells us what
happens at two arms.

STATUS: RESOLVED 2026-11-11 by the standing §5 rule, ratified as Group-2 policy
DECISION: **Emit as built, under the maximal-dispersion convention already
applied.** This is not a new judgement -- it is `CLAUDE.md` §5's "Flag, do not
guess -- ON THE ROW", and the folder had already done exactly that: 4 of 116
lipids carry `positional_ambiguity = 1`, `confidence = medium` and a note naming
the choice. The 2026-11-10 Group-2 policy makes it explicit: a combinatorial
library member whose regiochemistry the source does not pin is emitted as the
intended product with the ambiguity declared, which is what the corpus already
does for `NG_2025`'s 332 unpurified Mannich rows and `LL_2026`'s 840.
What makes this safe rather than lazy is that the ambiguity is MEASURED, not
inferred from `n < s`: `stage03_enumerate.py` builds every per-site assignment
and counts distinct InChIKeys, which is why 3-tail H4 lipids are correctly NOT
flagged (2+1 and 1+2 are one molecule). A reviewer can re-run that count.

### D-2
**Whether the Fig. 4c in-vitro panel is the same assay as the screen.**
The main-text legends for Fig. 3, Fig. 4 and Extended Data Fig. 2 are not in
hand — only the ESM was downloaded, and the article itself is behind the
publisher's paywall. So the assay annotation for lipids P1–P6 (cell line,
cargo, dose, timepoint) is inherited from Supplementary Fig. 10's legend rather
than read.

What licenses that: **ALC-0315's three in-vitro replicates are the same three
numbers in both panels** — 957,212.2 / 948,389.9 / 816,056.5 in
`Supplementary Fig. 10` A2:D2 and in `Fig. 4` F3:F5. Identical triplicates are
one measurement, not two, so the P rows sit in the same experiment as the
screen. The same argument covers the in-vivo side: ALC-0315's Fig. 3i and
Fig. 4c in-vivo triplicates are also identical.

Residual risk: this shows the panels share a measurement, not that the legend
says what we wrote. A reviewer with the article PDF should check the four
fields in one minute.

STATUS: RESOLVED — inherited annotation, with the identity above as the warrant

### D-3
**The in-vivo route, dose, strain and timepoint are unknown and are emitted
blank.** 19 in-vivo rows (Fig. 3i: ALC-0315 + 12 library lipids; Fig. 4c:
ALC-0315 + P1–P6) carry whole-animal luciferase total flux with no protocol
attached. The SI's Methods describe the Cy5 biodistribution study
(C57BL/6J male, IV, 0.75 mg/kg, 6 h and 24 h) and the mGFP cell-type study
(10 µg per dose, IV, days 0 and 1) — **neither is this experiment**, and
writing either one onto these rows would be a guess dressed as a protocol.

STATUS: RESOLVED 2026-11-11 -- blank, with the search recorded (Group-3 policy)
DECISION: **Emit the 19 in-vivo rows with `Route_of_administration`,
`Dose_ug_nucleicacid`, `Model_type` and the timepoint BLANK.** §3's rule is that
a blank means "not stated", and the alternative -- writing the Cy5
biodistribution protocol (C57BL/6J, IV, 0.75 mg/kg) or the mGFP protocol (10 ug,
IV, d0+d1) onto rows that are neither study -- would be a guess dressed as a
protocol, which is the F-004 failure the project already reverted once.
THE SEARCH PASS Group-3 policy requires, done 2026-11-11 before blanking:
  * Europe PMC REST, `EXT_ID:41851281` -> **no PMCID, `isOpenAccess: N`**, the
    only full-text resource listed is "Subscription required".
  * No PMC mirror exists, so `fullTextXML` is unavailable.
  * The ESM in hand carries only the two protocols named above, neither of which
    is this experiment.
So the fields are not recoverable from anything held or publicly reachable. A
reader with institutional access to Nat Biomed Eng can fill four fields in a
minute; until then the blank is the honest value and this entry says where we
looked.

### D-4
**Extended Data Fig. 2d is not emitted.** It reports ALC-0315 at 37,328.6 RLU
where Fig. 2c of the same figure reports 957,212.2 — a factor of 26 — so it is
a different assay, different cell line or different cargo. Its legend is in the
article, which is not in hand. Six triplicates (ALC-0315, Lipid T1–T6) are
therefore left out rather than labelled with the screen's annotation.

STATUS: RESOLVED 2026-11-11 -- not emitted, which is the conservative half
DECISION: **Leave Extended Data Fig. 2d out.** Its ALC-0315 reads 37,328.6 RLU
against 957,212.2 for the same compound in Fig. 2c -- a factor of **26** -- so it
is a different assay, cell line or cargo, and its legend is in the paywalled
article. Group-1 policy: where two artifacts of the same rank disagree and
nothing adjudicates them, emit nothing and record the contradiction. Labelling
six triplicates with the screen's annotation would import a 26x scale error into
a z-score. The values stay in the folder so the decision is reversible if the
article is ever acquired.

### D-5
**Lipid T1 and Lipid T2 are emitted with a measured value and no structure.**
They are the two ML-test lipids outside the H–T library (ML IDs 114 and 115),
drawn only in Supplementary Fig. 12 and built by a different chemistry
(epoxide ring-opening, not Michael addition). Their arm formulas *are* pinned,
by HRMS, and the head is H4 in both:

* Lipid T1, printed 909.8625 Da: H4 + 3 × **C18H32O**, to 0.0000 Da — a
  1,2-epoxide of an 18-carbon diene, consistent with the drawn 2-hydroxy arms
  carrying two cis alkenes.
* Lipid T2, printed 1053.9259 Da: H4 + 3 × **C19H36O3**, to 0.0000 Da — the
  same epoxide motif with a terminal ester of a branched C8 alcohol.

What is NOT established is the connectivity mass is blind to: where the alkenes
sit in T1, and where the ester sits along T2's acid chain. Those need the same
trace-plus-formula treatment the library tails got, on Supplementary Fig. 12.
Emitted with `structure_assumed = 1`, `confidence = low`, empty SMILES.

STATUS: RESOLVED 2026-11-11 -- the two rows are NOT emitted
DECISION: **Drop Lipid T1 and Lipid T2 from the emitted rows** and keep their
HRMS analysis in the folder.
The original plan was to emit them with an empty `IL_SMILES`,
`structure_assumed = 1` and `confidence = low`. Two reasons that is wrong:
  * **The pre-merge gate refuses a blank `IL_SMILES`** (added 2026-11-10), and
    rightly: the corpus is keyed on structure, every descriptor and InChIKey is
    computed from that column, and a row with no structure cannot enter Result 3
    or any fingerprint. A documented exception would be a hole in the one gate
    standing between an unresolved question and `data/`.
  * A value with no molecule attached is not a structure-activity datum. It is 2
    rows of 140.
What IS established stays recorded here and is worth keeping: both arm formulas
are pinned by HRMS to **0.0000 Da** (T1 = H4 + 3 x C18H32O at 909.8625; T2 = H4
+ 3 x C19H36O3 at 1053.9259), and the head is H4 in both. What mass cannot see
is where T1's alkenes sit and where T2's ester sits along the chain -- the
`CLAUDE.md` §5 "mass cannot see connectivity" lesson. If Supplementary Fig. 12
is ever traced the way the library tails were, these two become emittable with a
structure and should be revisited then.

### D-6
**The three commercial benchmarks are not formulation-matched to the library.**
Supplementary Table 3: the library screens at
`ionisable lipid:DOPE:Chol:DMG-PEG2000 = 45:15:38.5:1.5`, while MC3 and SM-102
use `DSPC` at `50:10:38.5:1.5` and ALC-0315 uses `ALC-0159` at
`46.3:9.4:42.7:1.6`. All eight formulations in that table sum to exactly
100 mol% (`07_qc.json` P8), so this is a deliberate design, not a transcription
slip. A cross-source comparison that treats MC3 here as the same LNP as MC3
elsewhere is comparing two formulations.

STATUS: RESOLVED 2026-11-11 by policy P3, flagged for the analysis phase
DECISION: **Emit the three benchmark rows with the formulation difference on the
row**, and carry the comparability question to the analysis phase rather than
answering it here.
The difference is real and deliberate, not a slip: the library screens at
`IL:DOPE:Chol:DMG-PEG2000 = 45:15:38.5:1.5` while MC3 and SM-102 use **DSPC** at
`50:10:38.5:1.5` and ALC-0315 uses **ALC-0159** at `46.3:9.4:42.7:1.6`, and all
eight formulations in Supplementary Table 3 sum to exactly 100 mol%
(`07_qc.json` P8).
P3 (2026-11-10) is exactly this question -- whether a source whose only shared
chemistry is the universal benchmarks should contribute to cross-source
concordance -- and the user's instruction was to defer it to the analysis phase
and measure the sensitivity both ways rather than argue it during collection.
The row carries the formulation, so the comparison is auditable either way. Note
this is the same hazard `CLAUDE.md` §8 already records for `XD_2026`.
is not mine to set

### D-7
**15 of these 116 lipids are already in the corpus, under other labs.**
Full-InChIKey matches, reported not acted on; see the overlap table in
`07_qc.json` and the run report. They are head-group/tail combinations that the
Akinc 2008 and Whitehead 2014 lipidoid libraries also contain, which is
expected for amine + acrylate Michael chemistry. Whether a re-measured lipid
enters the corpus twice, under two `Experiment_ID`s, is a corpus-level
decision — and it is exactly the shared-chemistry situation Result 3 is built
on, so it should be made deliberately.

STATUS: RESOLVED 2026-11-11 -- emit, and the overlap is now the most valuable
thing in this folder
DECISION: **Emit all 116, including the 15 that match corpus lipids by full
InChIKey.** The overlap is not a duplication problem, it is the point.
F-082 (2026-11-10) measured what this folder could not know when the question
was written: **79.1% of the in-vitro corpus's shared-lipid instances are WITHIN
one laboratory**, and the pair carrying 57% of Result 3's pooled concordance
points (`AA_2008` <-> `KW_2014`) is Anderson <-> Anderson. The entire genuinely
cross-laboratory shared chemistry of the corpus is one 55-lipid pair plus ten
pairs holding 12 lipids between them -- **67 instances in total**.
These 15 match `AA_2008`, `KW_2014` and `NC_2024` and come from a different lab,
so they are a **~22% increase in the corpus's entire genuine cross-laboratory
shared-lipid pool.** That makes this folder Result-3 material rather than
another silo, and it is the highest-value curation target the run produced.
Two conditions, both now met: the join must be on **full InChIKey**, not
`IL_name` (§3 -- 288 names carry different structures in different experiments,
and these matches are only visible on the key: `H2-OT7(2)` = `80_O10` =
`A80_C10`); and the cross-lab claim must be checkable, which the new `lab`
column makes it. The re-measurement question D-7 raised -- whether a lipid
measured twice in different labs is one datum or two -- is answered by the
corpus design: two rows, two `Experiment_ID`s, two `lab` values, and Result 3
exists precisely to compare them.
