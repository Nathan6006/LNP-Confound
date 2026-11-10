# DECISIONS -- siegwart_40644303

Questions this curation must NOT answer on its own. `build_dataset.py` refuses
to merge this folder while any entry below is `OPEN`.

WHY THIS FILE GATES THE MERGE. The pipeline is good at detecting that a source
is silent. It is not entitled to resolve it. In `AC_2025` the first repair was
confident, internally consistent, and would have assigned **419 measured values
to the wrong molecules** -- the names were wrong and the structures were right,
and nothing in the data said which. So the rule is: detect, queue, and stop.

## Format

```
### D-n  <one line>
STATUS:    OPEN | RESOLVED <date> by <who>
AMBIGUITY: what the source does not say.
CANDIDATES: each reading, and what it implies for the emitted rows.
EVIDENCE NEEDED: what would settle it.
DECISION:  (when resolved) what was chosen and why. Keep the entry.
```

## Entries

### D-1  Are Figure 2b and Figure 3b one in-vivo experiment or two?
STATUS: OPEN
AMBIGUITY: Both panels are labelled "Avg. Liver Radiance (p/sec/cm2/sr)" from
the same assay -- C57BL/6, 0.25 mg/kg, 6 h, IVIS. Figure 2b's axis tops at
1e6 and its largest bar is 7.37e5; Figure 3b's tops at 6e7 and its largest is
4.09e7. That is a ~55x offset. The two panels share no lipid and no benchmark
control, and their IVIS colour bars are 100x apart (x10^5 against x10^7), so
nothing in the paper ties them to a common scale.
CANDIDATES:
  (a) ONE experiment. Pool all 40 values and z-score together. The 3A4 family
      then sits ~1.7 log units below the rest and every 3A4-family row becomes
      a strong negative. If the offset is a cohort/imaging-session effect this
      fabricates a chemical conclusion.
  (b) TWO experiments, `SL_2025` and `SL_2025_2`, z-scored separately. This is
      what is emitted now, marked PROVISIONAL on every row. It asserts an
      offset exists and deletes any real between-family difference.
  (c) Emit neither until a third source fixes the scale.
EVIDENCE NEEDED: whether the two cohorts were imaged in one session -- a
methods statement, the raw IVIS files, or any lipid measured in both panels.
The authors' Figure 4b and 5a replot both panels against size and pK_a on
*separate axes per family*, so they do not settle it either.
IMPLICATION: every in-vivo `Experiment_value` in `09_rows_vivo.csv` changes
whichever way this goes.

### D-2  Should the siRNA half of Figure 1a be emitted at all?
STATUS: OPEN
AMBIGUITY: Figure 1a's upper stripe is 465 siRNA knockdown values at three
levels. The main text says plainly that this is a re-analysis of reference 11
(Zhou et al., PNAS 2016, PMID 26729861), which the corpus already holds as
`KZ_2016` -- same 465 lipids, four levels, 1,512 rows. The caption, however,
says the heat map is "in IGROV-1 cells", while `KZ_2016` carries `Model_type =
HeLa` on every row.
CANDIDATES:
  (a) Do not emit (current). Emitting would put the same laboratory's same
      measurements under a second `Experiment_ID` at coarser resolution --
      the "fabricate a lab" failure that kept `XD_2026`'s LipidDB out.
  (b) Emit as a separate IGROV-1 siRNA experiment, if the caption is right and
      these really are new measurements in a different cell line.
  (c) Emit nothing and additionally correct `KZ_2016`'s `Model_type`, if the
      caption is right about the cell line and the corpus label is wrong.
EVIDENCE NEEDED: PMID 26729861's own methods -- which cell line its siRNA
screen used. We did not fetch that paper.
NOTE: our read of this stripe correlates with `KZ_2016`'s stored values at
Spearman 0.899 over all 465, which says the two describe the same measurements;
it does not say which cell line they were made in.

### D-3  Two in-vivo lipids have no readable value. Drop, or record a bound?
STATUS: OPEN
AMBIGUITY: `3A4V-SC5` has NO bar drawn in Figure 2b -- zero coloured pixels
above the baseline -- while SI Figure S3 gives it a diameter, PDI, pK_a and
E.E. like every other formulation, and Figure 5a plots its marker sitting on
the x-axis rule at pK_a 6.14. `5A4-SC10` has a bar about 5 px tall and
`bars.bar_top` REFUSED it: the left flank has no covered row, which at that
height cannot be distinguished from a marker.
CANDIDATES:
  (a) Drop both rows (current). 40 of 42 in-vivo values are emitted.
  (b) Emit `3A4V-SC5` with `value_is_bound = upper` at the smallest value one
      bar-line width represents, and `5A4-SC10` likewise.
  (c) Re-digitise from Figure 5a instead -- but its marker for `3A4V-SC5` is
      merged into the axis rule, so that route does not recover it either.
EVIDENCE NEEDED: the underlying radiance numbers, or a higher-resolution
rendering of Figure 2b than the published 1500 px raster.
NOTE: lowering `bars.bar_top`'s threshold to make `5A4-SC10` pass was NOT
done. That threshold is what F-018 cost: 0.72 coverage of a trimmed window
read replicate dots as bar tops on 15 of 60 bars.

### D-4  SC3: n-propyl (this paper's Figure 1b) or neopentyl (corpus KZ_2016)?
STATUS: OPEN
AMBIGUITY: Figure 1b draws SC3 as an unbranched chain. Measured, not counted:
the raster trace finds 0 skeleton branch nodes, 3 bond segments and a bounding
box 31 x 9 px, where the two genuinely branched tails SC1 and SC2 measure 1 and
3 branch nodes and 18-19 px tall. That is n-propanethiol. The corpus's
`KZ_2016` `IL_SMILES` encodes neopentyl on all 15 SC3 lipids -- and the same
rows' `IL_tail1_SMILES` says `CCCS`, and their `Molecular.Weight` is the C5
value, so the corpus states the tail three ways and they disagree.
CANDIDATES:
  (a) n-propyl, from this paper's drawing (current). 15 structures.
  (b) neopentyl, from the upstream metadata CLAUDE.md section 9 cites.
EVIDENCE NEEDED: the reagent list of PMID 26729861 (which thiol "SC3" names),
or any HRMS/NMR for an SC3 lipid. Note HRMS would settle THIS one, because
propyl and neopentyl differ by C2H4 -- unlike the regiochemistry questions,
this is a mass difference.
IMPLICATION: 15 of our 465 structures, and a correction to `KZ_2016` if (a).

### D-5  Is Figure 1a's lowest mRNA band closed at zero, or open below?
STATUS: OPEN
AMBIGUITY: The legend prints "log 0-1", "log 1-2" and "> log2". Taken
literally the lowest band is CLOSED, so a cell in it is a measurement in
[0, 1], not a colormap floor. Taken as what a screen usually produces, "0-1"
is where everything undetectable lands, and the band is effectively "<= 1".
325 of 465 cells are in it.
CANDIDATES:
  (a) Closed band, no bound flag (current). Only the 68 cells in the open top
      band carry `value_is_bound = lower`.
  (b) Treat the lowest band as a floor and flag 325 rows `value_is_bound =
      upper`, taking this curation's caveat count from 180 to ~437.
EVIDENCE NEEDED: the raw luminescence distribution, or a statement of the
detection limit. F-040 is the warning against guessing: 65 `QL_2026` rows were
flagged as bounds and were measurements.

### D-6  3A2: piperazine (Figure 1b) or the corpus's two-nitrogen core?
STATUS: OPEN
AMBIGUITY: Figure 1b draws 3A2 with a six-ring carrying "N" on one side and
"NH" on the other plus a pendant NH2 -- 1-(2-aminoethyl)piperazine, three
nitrogens, which is also what three arms requires. The corpus's `KZ_2016`
`IL_SMILES` has two nitrogens. The molecular formulae differ: ours
C51H93N3O12S3 against the corpus's C52H94N2O12S3 for 3A2-SC6.
CANDIDATES:
  (a) the piperazine, from this paper's drawing (current), 31 structures;
  (b) the corpus structure.
EVIDENCE NEEDED: the amine list in PMID 26729861, or HRMS -- the formulae
differ, so a mass settles it.

### D-7  The paper is not CC0/CC-BY-compatible. Can it enter a released corpus?
STATUS: OPEN
AMBIGUITY: CLAUDE.md section 8 requires every constituent source to be
CC0/CC-BY-compatible before corpus release. The Europe PMC record for PMID
40644303 returns `isOpenAccess = N` and `license = null`; the Europe PMC
supplementary endpoint returns "Article with id PMC12869450 is not open access
one". The copy we hold is an NIH Public Access **author manuscript** (NIHMSID
NIHMS2133731) -- free to read, no CC licence -- and the publisher version is
subscription-only.
CANDIDATES: (a) include, on the ground that facts are not copyrightable and we
redistribute values rather than text; (b) exclude; (c) include but withhold the
figures from any release.
EVIDENCE NEEDED: the project's position on author-manuscript sources. At least
one corpus source is likely already in this category, so the answer is not
specific to this paper.

### D-8  The partial-substitution convention on 5A2 and 6A1
STATUS: OPEN
AMBIGUITY: Five cores carry fewer arms than they have N-H: 5A2 (5 arms on 7),
6A1 (6 on 7), 5A3 (5 on 6), 5A4 (5 on 6), 5A5 (5 on 6). The paper says how
many arms and never which nitrogen keeps its hydrogen. For 5A3, 5A4 and 5A5
the convention used here ("primary N-H first, then secondary from the terminus
inward") has a unique product by symmetry, and it AGREES with the corpus on all
three. For 5A2 and 6A1, on tetraethylenepentamine, it does not: our structures
and the corpus's have identical molecular formulae and different placements.
CANDIDATES: (a) our convention, 62 structures; (b) the corpus's; (c) treat both
as one structure class and accept that the stored SMILES is a representative.
EVIDENCE NEEDED: nothing in either paper distinguishes them -- a one-pot
Michael addition on a polyamine gives a mixture. This may be permanently
unresolvable, in which case the right answer is the flag, not a structure.
NOTE: all 135 affected emitted rows already carry `positional_ambiguity = 1`,
`confidence = medium` and a note naming the convention.

### D-9  Is a three-level ordinal an acceptable `Experiment_value`?
STATUS: OPEN
AMBIGUITY: The 465 in-vitro values are an ordinal with three levels, z-scored
per CLAUDE.md section 3's rule for a discretized source (z-score the BIN
ORDINAL, not a midpoint). `KZ_2016` and `SL_2021` are already in the corpus at
four levels. Three levels is coarser than anything currently there, and 325 of
the 465 fall in one level.
CANDIDATES: (a) emit as-is (current); (b) emit only the two upper bands;
(c) exclude the experiment from rank-based analyses and say so in Methods.
EVIDENCE NEEDED: a project decision about the minimum resolution a
rank-correlation input must have. This bears directly on T1/Result 3, because
of D-10.

### D-10  This experiment shares all 465 of its lipids with KZ_2016. Is that a shared-lipid pair, or a duplicate?
STATUS: OPEN
AMBIGUITY: Every one of the 465 in-vitro lipids is already in the corpus under
`KZ_2016`. If merged, `SL_2025` would create the largest shared-lipid overlap
in the corpus by a wide margin -- 465 lipids against the 231 Stage 5 currently
has in total, and 57% of its pooled points already come from one pair
(AA_2008 - KW_2014). But the two experiments differ in CARGO (siRNA in
KZ_2016, mRNA here), and the whole point of Result 3 is lipids measured for
the same endpoint in different laboratories.
CANDIDATES:
  (a) treat `SL_2025` x `KZ_2016` as a shared-lipid pair for Result 3. It
      would dominate the pooled concordance and it is a cross-CARGO
      comparison, not a cross-laboratory one -- arguably a different claim.
  (b) exclude the pair from Stage 5 and keep the rows for the source-
      classification and protocol analyses only.
  (c) exclude the experiment.
EVIDENCE NEEDED: a project decision, not a measurement. It is queued here
because merging this folder without deciding it would silently change
Result 3's headline number, and CLAUDE.md section 8 already flags the related
hazard for in-vivo benchmark-only pairs.
NOTE: this is also an opportunity, not only a risk -- the same 465 molecules,
same laboratory, different cargo, is the cleanest "one variable changes"
comparison available anywhere in the corpus, and it is a closer analogue of the
`molea.csv` within-lab cross-cell-type point than of a cross-laboratory pair.
