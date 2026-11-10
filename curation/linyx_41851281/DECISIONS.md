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

STATUS: OPEN

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

STATUS: OPEN — needs the main-text Methods

### D-4
**Extended Data Fig. 2d is not emitted.** It reports ALC-0315 at 37,328.6 RLU
where Fig. 2c of the same figure reports 957,212.2 — a factor of 26 — so it is
a different assay, different cell line or different cargo. Its legend is in the
article, which is not in hand. Six triplicates (ALC-0315, Lipid T1–T6) are
therefore left out rather than labelled with the screen's annotation.

STATUS: OPEN — needs the main-text Extended Data legend

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

STATUS: OPEN

### D-6
**The three commercial benchmarks are not formulation-matched to the library.**
Supplementary Table 3: the library screens at
`ionisable lipid:DOPE:Chol:DMG-PEG2000 = 45:15:38.5:1.5`, while MC3 and SM-102
use `DSPC` at `50:10:38.5:1.5` and ALC-0315 uses `ALC-0159` at
`46.3:9.4:42.7:1.6`. All eight formulations in that table sum to exactly
100 mol% (`07_qc.json` P8), so this is a deliberate design, not a transcription
slip. A cross-source comparison that treats MC3 here as the same LNP as MC3
elsewhere is comparing two formulations.

STATUS: OPEN — flagged on the row; the corpus convention for benchmark rows
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

STATUS: OPEN
