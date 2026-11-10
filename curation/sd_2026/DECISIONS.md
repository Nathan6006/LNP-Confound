# SD_2026 — open decisions

Dong et al., *Nature Biotechnology* 2026, PMID 42806111,
DOI 10.1038/s41587-026-03298-8. Curated 2026-11-10.

Every entry below is a place where the source is silent, self-inconsistent, or
offers several mutually exclusive readings. The pipeline detects that; it does
not get to decide what the source meant. **Nothing here is resolved by this
curation, and `build_dataset.py` refuses to merge this folder while any entry
reads `STATUS: OPEN`.**

---

### D-1 The 384-lipid screen has values that cannot be read, and structures that can

STATUS: OPEN

AMBIGUITY: `03_enumerated.csv` holds all 384 structures. `05_values.csv` holds
18 values, none of them for a member of that library except LC-1. The screen's
own potencies exist only in Fig. 1c/1h (circos) and Fig. 1d/1i (Sankey), and
both panels are unreadable for two independent reasons, either of which alone
would be fatal:

1. **No numeric scale anywhere.** The heat-map colour legend is a continuous
   ramp annotated `Max` at one end and `Min` at the other — no number, no tick,
   no units. The inner bar ring has no axis and no scale rule. The Sankey's y
   axis, labelled `RLU (log2) of Ionizable lipid`, is drawn as a bare arrow
   with zero ticks. The only anchor in any of them is a dashed ALC-0315
   reference line, and one anchor cannot fit a scale — `axis.calibrate`
   requires three and would be right to refuse.
2. **No cell → lipid map.** The Fig. 1 caption states the inner rings are
   "hierarchical clusters of linkers (L1–L4), tails (Ta1–Ta4) and subtail
   (Tb1–Tb2) categories". The within-sector order is therefore a clustering's
   leaf order, and that leaf order is published nowhere — not in the caption,
   not in Methods, not in the SI. The dendrogram drawn at the top of the circle
   shows the nesting, not which leaf is which (H, Ta, Tb, L) combination.

CANDIDATES:
 (a) Merge the 384 structures with **no value**, as structure-only rows. The
     corpus schema has no such row type today.
 (b) Merge only the 18 valued rows (3 experiments) and keep `03_enumerated.csv`
     in the folder as provenance for a future value source.
 (c) Merge nothing from this paper.
 (d) Write to the corresponding author (bw.li@utoronto.ca) for the screening
     table and the circos leaf order, and hold the paper until it arrives.

EVIDENCE NEEDED: the authors' screening table, or any figure-source file that
states the leaf order. The publisher serves four assets and none is a workbook
(`00_acquisition.json → springer_source_data_probe`: n = 1..16 × seven
extensions, no `PK` magic anywhere); there is no author repository; the Data
availability statement names no accession.

---

### D-2 Which of the 384 members are SC-1…SC-5 and LC-2…LC-5?

STATUS: OPEN

AMBIGUITY: the paper draws all ten lead structures (Fig. 1g, Fig. 1l) and the
SI characterises each by formula, MW, ¹H NMR, ELSD-HPLC and MS. **Only LC-1 is
pinned to a block combination**, by Supplementary Table 5's five-step synthesis
with CAS numbers (H4 + L4 + Ta1 + Tb1, reproducing C50H97N3O4 / MW 804.34
exactly). For the other nine the printed formula matches between **2 and 12**
enumerated members (`07_qc.json → P7`), and mass cannot choose among them:
seeding LC-1's headgroup H4 → H3 leaves the molecular formula *identical* while
the InChIKey skeleton moves (`P7 seed B`, `P12`).

CANDIDATES:
 (a) Read the nine drawn structures off Fig. 1g/1l and match them by InChIKey
     skeleton. The structure strips are embedded rasters (~1378 × 280 px for
     five structures), so each structure is ~270 px wide — far below what the
     Supplementary Fig. 1 block read needed, and chain counting at that scale
     is the ±1 failure class.
 (b) Match by ¹H NMR substructure signatures: a 9H singlet near δ 1.0 is the
     pivaldehyde *tert*-butyl (L4), a 6H doublet near δ 1.0 is isobutyraldehyde
     (L3), olefinic protons near δ 5.35 are Ta4. This narrows but does not
     uniquely determine.
 (c) Leave all nine unmapped and emit their Fig. 1n values against the NAME
     only, with no structure — which is what `05_values.csv` does now, and is
     why `reference_structures.csv` carries nine rows with an empty
     `reference_SMILES`.

EVIDENCE NEEDED: a higher-resolution rendering of Fig. 1g/1l than the embedded
raster (the publisher's `media.springernature.com` figure assets), or the
authors' compound list.

---

### D-3 Supplementary Fig. 53: four more organs, by raster digitisation

STATUS: OPEN

AMBIGUITY: Supplementary Fig. 53a–d gives IV radiance in **spleen, heart, lung
and kidney** for the same 12 LNPs as Fig. 1n, n = 5, with labelled decade ticks
— so unlike the circos it *is* calibratable. But the SI page carries one
flattened raster at 1430 × 862 px for all four panels, so this is genuine
raster digitisation at roughly 150 px per decade, a weaker evidence class than
Fig. 1n's vector read. Extracting it would add 48 rows and four organs.

CANDIDATES:
 (a) Extract it and emit the four organs as four `Model_target` values inside
     one IV experiment, flagged `5b_raster`.
 (b) Extract it as a separate `Experiment_ID` so the evidence classes are never
     mixed within one block.
 (c) Leave it out, as now.

EVIDENCE NEEDED: a decision on whether a vector-read value and a raster-read
value may share an `Experiment_ID` — the corpus has no convention for this, and
`value_normalization_unit` records the normalisation unit, not the evidence
class.

---

### D-4 A second, separate lipid library in Supplementary Fig. 68

STATUS: OPEN

AMBIGUITY: Supplementary Fig. 68 reports a **different** library, made by
conventional stepwise chemistry rather than the Ugi, with its own named leads
(N-SC-1…3, N-LC-1…3 and two analogues, all with formulas and MWs in the SI) and
its own in-vivo panels. It is a second screen from the same paper, which the
corpus convention would make `SD_2026_4` or similar. Its heat maps (68c, 68d)
have the same defect as Fig. 1c/h — no numeric scale, no published member
order; its in-vivo panels (68e–g, 68i–j) do have axes.

CANDIDATES:
 (a) Curate it as a further experiment.
 (b) Leave it; the paper's own argument rests on the Ugi library.

EVIDENCE NEEDED: none to decide; this is a scope call.

---

### D-5 The licence is not Creative Commons

STATUS: OPEN

AMBIGUITY: `CLAUDE.md` §8 requires every constituent source to be
CC0/CC-BY-compatible before corpus release. The landing page for this article
carries no `creativecommons.org/licenses/` link at all; its rights statement
reads "Springer Nature or its licensor (e.g. a society or other partner) holds
exclusive rights to this article under a publishing agreement…". Checked
2026-11-10 and recorded with the method in `00_acquisition.json → licence`.

CANDIDATES:
 (a) Exclude this paper from any released corpus.
 (b) Include only facts (structures, numeric values), on the view that
     individual data points are not the copyrighted expression.
 (c) Seek permission.

EVIDENCE NEEDED: the project's position on (b), which affects far more than
this paper.

---

### D-6 Three-lipid experiments and the z-score convention

STATUS: OPEN

AMBIGUITY: `SD_2026_2` (intrathecal, brain) and `SD_2026_3` (intratracheal,
lung) each contain **three** lipids, two of which are the benchmarks ALC-0315
and LP-01 that the corpus already holds from other sources. The corpus
convention standardises `Experiment_value` within an experiment; a z-score over
n = 3 is almost entirely determined by which three were chosen, and here two of
the three are the field's universal positive controls. `CLAUDE.md` §8 already
flags exactly this hazard for the in-vivo pool ("a source whose only shared
lipids are benchmarks that top every screen can push cross-lab concordance
toward +1 for reasons that have nothing to do with transferable SAR").

CANDIDATES:
 (a) Emit all three experiments; flag the n = 3 blocks.
 (b) Emit only `SD_2026` (n = 12) and drop the two n = 3 panels.
 (c) Emit the n = 3 panels but exclude them from Stage 5's pooled concordance.

EVIDENCE NEEDED: the standing decision in `CLAUDE.md` §8 on benchmark-only
pairs, which is itself open.

---

### D-7 What Fig. 1n's radiance is a radiance *of*

STATUS: OPEN

AMBIGUITY: Fig. 1n quantifies the whole-body IVIS image; its axis says only
"Radiance intravenous injection (photons per s per cm² per sr)". The Results
text says "all formulations produced hepatic expression", and the per-organ
breakdown is a different figure (Supplementary Fig. 53), which does **not**
include liver — liver is only in Fig. 1n and in the ex-vivo DOE panels. So the
`Model_target` for these 12 rows is "whole body" on the figure's own wording
and "liver" on the text's.

CANDIDATES:
 (a) `Model_target = whole_body_luminescence`.
 (b) `Model_target = liver`, on the Results sentence.
 (c) Emit as (a) with a `curation_note` naming the text's claim — what
     `05_values.csv` does now.

EVIDENCE NEEDED: none from the paper; it is a mapping decision onto the corpus
vocabulary.
