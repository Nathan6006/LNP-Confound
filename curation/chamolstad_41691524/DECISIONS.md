# DECISIONS — chamolstad_41691524

Kim, Xu, Lee et al., *Comprehensive Engineering of Ionizable Lipid Nanoparticles
and mRNA Elements for Next-Generation Vaccines*, ACS Nano 2026.
PMID 41691524 · DOI 10.1021/acsnano.5c19036.

Every entry below is OPEN, so `build_dataset.py` will refuse this folder. That
is the intended state: nothing here should merge on one agent's judgement.

---

### D-1
Should this paper enter the corpus at all? What is emittable is a **formulation
DoE on a single ionizable lipid** — 77 rows, 1 distinct structure — not the
96-member lipid library the paper is about. The library screen lives in
main-text Figure 1B–G, which is not on disk and is behind a paywall, and 11 of
its 12 head amines are not disclosed anywhere in the SI, so even with the
figure 88 of 96 members could not be given a structure. `ZT_2026` is the
precedent for merging a single-lipid DoE (486 rows, one lipid). The question is
whether a second such source is worth a 44th `Experiment_ID` when it adds one
chemical point and zero shared lipids.

STATUS: OPEN

### D-2
Figure S3A and Figure S4C report the same formulations, the same cell line
(C2C12), the same readout and the same timepoint, and they differ in magnitude
by a **median factor of 8.96** while rank-agreeing at **Spearman 0.935**
(n=12, p=8.1e-06). Only S4C is emitted. If this paper is merged, is that the
right choice, and are Figure S4A (DC2.4) and Figure S4C (C2C12) one
`Experiment_ID` or two? They are different cell lines, which §3's convention
would normally make one experiment with two `Model_type` values — but the
S3A/S4C factor of 9 shows the authors' own "relative to SM102" is not a stable
unit across their figures, so pooling the two panels under one z-score is not
obviously safe.

STATUS: OPEN

### D-3
**H9T6 is not a pure compound.** SI p.5 reports the MPLC-purified material as
"a mixture of several compounds", with the 4-tailed congener at 55.26% of total
HPLC–ELSD area, while the 1H NMR of the same material integrates to six tails.
We emit the designed hexa-substituted structure with `structure_assumed = 1`
and the 55.26% in `curation_note`. Does a structure that is a minority-to-
majority component of its own vial belong in a structure–activity corpus, and
if so at what `confidence`? Compare `NG_2025`, whose crude Mannich products are
carried with `positional_ambiguity = 1`.

STATUS: OPEN

### D-4
Table S1's lead-formulation summary row names **DOPE** where the abstract and
Table S1's own 1-F16 row name **DOPG** (contradiction C-01). We decline to use
the summary row. Confirm that is the right half to trust before anything in
this folder is cited for the lead formulation — the `AC_2025` precedent is that
the confident repair put 419 measured values on the wrong molecules.

STATUS: OPEN

### D-5
Table S1 prints the PEG-lipid identity on only 5 of 45 rows. We read it as a
merged cell spanning each block and assign **DMG-PEG2k** to all 1-F and 2-F
formulations (contradiction C-06). That is an inference about table layout.
The ALC-0315 row names ALC-0159 explicitly, which is consistent with the
column marking changes only — but it is not proof.

STATUS: OPEN

### D-6
**H9T2 is not emitted.** Figure S2 pairs it with H9T6 as the R-configured
counterpart, and the index gap of 4 matches the four chain lengths in the SI's
tail list, which would make T2 = (R)-glycidyl decanoate. But that is index
arithmetic on a numbering the paper never prints, and `KS_2024`'s `C24` (two C4
tails) is the standing reason not to trust a name token. Its only values are in
Figure S2B/C, whose y-axes are **broken** (two segments, different scales) and
which `axis.calibrate` cannot describe with one map.

STATUS: OPEN

### D-7
Nine bars across the two emitted panels were **refused**, not estimated: 4 in
S4A, 5 in S4C. Four of the nine have no measurable bar at all (the row above
the axis rule is 0.00 covered — the value is below one pixel, roughly <0.02);
five have a top stroke 9–16 px thick against a 2–3 px median, which is an error
bar's cap merged with the bar top on a very short bar. The refusals therefore
correlate with small values, so the emitted set is **not a random subset**. If
this paper is merged, the refused formulations must be recorded as missing, not
as absent-because-zero.

STATUS: OPEN

### D-8
The SI gives **no dose** for the in-vitro transfections in Figures S3/S4, and
no `n` per bar beyond "representative of two independent experiments". Those
fields are left blank rather than filled from the main-text methods, which are
not on disk.

STATUS: OPEN
