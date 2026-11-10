# DECISIONS — andersondg_41203968

Rudra, Gupta, Reed et al., *Degradable cyclic amino alcohol ionizable lipids as
vectors for potent influenza mRNA vaccines*, Nature Nanotechnology 2025.
PMID 41203968 · DOI 10.1038/s41565-025-02044-6 · PMC13387048.

Everything here is something **the source does not settle** or something that
would change a published number. None of it is resolved in this folder.
`build_dataset.py` refuses to merge a folder while any entry below is `OPEN`.

---

### D-1 — One experiment or five? The in vivo panels are not one screen.
STATUS: OPEN
AMBIGUITY: The 100-member in vitro screen (Extended Figure 1) is a single
coherent experiment: one cell line, one dose, one timepoint, n = 3 wells. The
in vivo values are four separate panels at **four different dose/timepoint
combinations** — Figure 1f at 1 µg / 24 h, Figure 2a at 5 µg / 24 h, Figure 2b
and 2c at 1 µg / 6 h, Extended Figure 2b at 0.1 µg / 6 h — and each panel is a
different generation of the library, so they share almost no chemistry.
`Experiment_ID` is the corpus's source/environment key and a z-score is taken
within it, so this choice sets what is being standardised.
CANDIDATES:
  (a) one `DG_2025` vivo experiment pooling all 32 rows — but a 50-fold dose
      range inside one z-block deletes the dose effect, which is error class 22
      (`value_normalization_unit = experiment_x_model_type` exists because of
      exactly this);
  (b) four experiments keyed by panel (`DG_2025_vivo_gen2_tails` and so on) —
      honest, but three of the four have fewer than 11 lipids, below anything
      Stage 5 can use;
  (c) emit only the in vitro screen and drop the in vivo rows, as §4 T8 already
      does for the in vivo pool generally.
  `06_emitted_rows.csv` currently writes (b) as a PROPOSAL in a column named
  `Experiment_ID_proposed`, so nothing is committed.
EVIDENCE NEEDED: a project decision on whether a dose-series panel of six
lipids is an "environment" for Results 1–4, not more reading of this paper.

### D-2 — Four lipids whose structures this curation could not establish.
STATUS: OPEN
AMBIGUITY: AMG1441, AMG1641, AMG1741 and AMG1941 each carry a value in Figure
2b, and each has exactly one drawing in the SI. Their headgroups (Am14, Am16,
Am17, Am19) are large polyaza macrocycles drawn at 5 pt with five to seven
nitrogen labels and six hydroxyls packed inside the ring. In all four, the
reconstruction puts an arm's hydroxyl oxygen into the macrocycle — the
decomposition reports two or three distinct arms where there should be one —
and because those two defects conserve the molecular formula, **the HRMS gate
passes a structure that is wrong** (AMG1641 and AMG1741 both PASS at < 2 ppm).
That is the exact reason this folder does not treat the mass check as
sufficient.
CANDIDATES:
  (a) drop the four rows — the 32 in vivo rows become 28;
  (b) re-read those four headgroups by hand from a high-resolution render and
      store them as read-once blocks, with the loss of the ten-independent-read
      check every other block has;
  (c) take Am16/Am17/Am19 from the vendor/reagent names if the authors state
      them anywhere outside the drawing (not found in the SI text layer).
EVIDENCE NEEDED: a human reading of SI pp. 109, 118, 119 at high
magnification, or the headgroup names from the authors.

### D-3 — The two benchmarks on Figure 2c.
STATUS: OPEN
AMBIGUITY: Figure 2c measures SM-102 and cKK-E12 beside the AMG lipids, at the
same dose and timepoint. Both are already in the corpus from other sources, so
they are exactly the "identical lipid, different lab" chemistry that Result 3
is built on — and §8 already flags that a source whose only shared lipids are
the field's universal positive controls can push cross-lab concordance toward
+1 for reasons that have nothing to do with transferable SAR.
CANDIDATES: (a) emit both with corpus SMILES, adding one source to each;
(b) emit neither, as this curation currently does (`06_not_emitted.json`);
(c) emit them but mark them so Stage 5 can exclude benchmark-only pairs.
EVIDENCE NEEDED: the §8 decision on benchmark-only pairs in the pooled
concordance. This paper cannot settle it.

### D-4 — "the next best lipid AMG1010" names a lipid that is not in the panel.
STATUS: OPEN
AMBIGUITY: See `08_contradictions.csv` row 2. Screen 3 (Figure 1f) contains
Am8–Am10 × G4, G5 only. The printed "~5-fold" is reproduced exactly by
AMG104 / AMG105 = 4.71, so either the text's "AMG1010" is a slip for AMG105 or
the Source Data's "Tail 5" column header is wrong.
CANDIDATES: (a) text slip, Source Data right — the six rows stand as emitted;
(b) header wrong — three of the six rows are mislabelled AMG85/AMG95/AMG105
when they are AMG810/AMG910/AMG1010, which would put three measured values on
the wrong molecules.
EVIDENCE NEEDED: the published Figure 1f panel's own x-axis labels, or the
authors. The IVIS images in the figure may also show six distinct groups
labelled directly.

### D-5 — Two cells are exactly 0 RFU and the corpus z-scores log10.
STATUS: OPEN
AMBIGUITY: AMG41 and AMG51 are published as 0.0000 RFU with Stdev 0.0000 and
n = 3. `Experiment_value` is a per-experiment z-score of **log10**, and log10(0)
does not exist. They are emitted with `value_is_bound = upper`, which says the
truth is at or below the floor, but it does not say what number to put there.
The smallest non-zero value in the same screen is 28.91 RFU.
CANDIDATES: (a) drop the two rows; (b) substitute a floor (half the smallest
non-zero value, or the plate background) and keep the bound flag; (c) keep them
out of the z-score and re-insert them at the bottom rank.
EVIDENCE NEEDED: whichever convention the corpus already uses for censored
values — none is documented in `CLAUDE.md` §3 beyond the flag itself.

### D-6 — AMG106 fails its own HRMS, and the paper's formula is impossible.
STATUS: OPEN
AMBIGUITY: The SI prints "calculated for C114H226N4O12 [M+H]+", whose implied
neutral has a half-integer degree of unsaturation, so no such molecule exists;
the printed m/z agrees with that impossible formula to 1.1 ppm. Our structure
(C114H224N4O12) is two hydrogens lighter and re-assembles exactly from blocks
read 27 and 10 times elsewhere.
CANDIDATES: (a) keep the row at `confidence = medium` on the strength of the
block agreement; (b) drop it as the only HRMS failure inside the 10 × 10 grid.
It is currently emitted at `confidence = medium` with the failure in its note.
EVIDENCE NEEDED: nothing further from this paper — this is a judgement about
how much weight the block-consensus check carries when an external check fails.

### D-7 — Four library members the SI never draws.
STATUS: OPEN
AMBIGUITY: AMG41, AMG43, AMG73 and AMG83 carry Extended Figure 1 values but
appear nowhere in the SI's characterization section. They are emitted as
enumerated structures with `structure_assumed = 1` and `confidence = medium`.
That rests on the library rule — "10 amine headgroups × 10 tails" — being
exactly true for those four cells.
CANDIDATES: (a) keep them flagged, as now; (b) drop them, since the paper
claims all 100 were synthesised and purified but evidences only 96.
EVIDENCE NEEDED: a project rule on whether a combinatorial inference with no
drawing may enter the corpus at all.

### D-8 — The in vitro screen publishes a mean, not wells.
STATUS: OPEN
AMBIGUITY: Extended Figure 1's Source Data gives RFU, Stdev and n = 3 per
lipid; the three wells themselves are not released. Every other panel in this
curation releases the individual replicates.
CANDIDATES: (a) emit the published mean, as now; (b) require per-replicate data
for the corpus's primary in vitro pool.
EVIDENCE NEEDED: whether any corpus analysis uses within-lipid replicate
spread. Stage 5's split-half assay-noise ceiling (T2) does.

### D-9 — Licence, for the §8 release check.
STATUS: OPEN
AMBIGUITY: `00_acquisition_raw.json` records the URLs, bytes, sha256 and the
UTC timestamp but **not a licence string**, which `PIPELINE_V2.md` §3 asks for.
The PMC author-manuscript record states only "This file is available for text
mining. It may also be used consistent with the principles of fair use under
the copyright law" — that is not a CC0/CC-BY statement, and the Nature
Nanotechnology article itself is not open access.
CANDIDATES: (a) treat the extracted values as facts, not as a copyrightable
selection; (b) obtain permission; (c) exclude from any corpus redistribution
while keeping it for internal analysis.
EVIDENCE NEEDED: the licence line on the publisher's landing page, captured
with a date, and a project-level position on constituent-source licensing.
