# DECISIONS — dixonjb_41269090 (Michalaki et al., ACS Nano 2025, PMID 41269090)

Queued, not resolved. `build_dataset.py` refuses a folder with an `OPEN` entry,
which is the point: nothing here is merged.

---

### D-1 Does the pooled SANDS barcode screen belong in the corpus at all?

The screen's own readout is a pooled DNA-barcode assay. 99 QC-passing LNPs were
co-injected intradermally at 1.5 mg/kg; LECs with high aVHH protein were sorted
by FACS; the barcodes inside those cells were sequenced; the value per LNP is a
**normalised barcode count**, and the "total LN LEC delivery" label used for the
machine-learning analysis is the **product** of the normalised PLN, ALN and BLN
readouts.

Two separate objections, and they should be settled separately.

1. **The value is relative and competitive.** A normalised barcode count depends
   on which other LNPs are in the pool: adding one strong member lowers every
   other member's share without anything about those members changing. The
   corpus's `Experiment_value` is a per-experiment z-score of a log signal, which
   is a monotone transform of a per-member quantity, not of a share.
2. **It is a colocalisation measure, not a per-LNP delivery measure.** The paper
   says so: sequencing the sorted cells identifies "LNPs colocalized with cells
   in which functional delivery occurred". The expression gate is applied to the
   CELL, and the barcode then says which particles were in that cell. An LNP that
   merely accompanies a transfecting LNP into the same cell scores.

Moot in practice: the per-LNP SANDS values are published nowhere (see
`00_value_sources.json`). But the decision should be recorded before a future
curation of a Dahlman-lab barcode paper that *does* publish them.

STATUS: OPEN

---

### D-2 Does a screen with ONE ionizable lipid belong in a cross-laboratory chemistry corpus?

Every one of the six lead LNPs carries the **same** ionizable lipid, cKK-E12.
Screening 1 is a complete 5 × 6 factorial over helper lipid and mole ratios on a
single lipomer (verified: 5 helper lipids × 6 ratio sets = 30 = the row count,
every cell present once). The only second ionizable structure among
individually-administered LNPs is the **MC3 benchmark**.

So an `EM_2025` source would contribute, to a corpus built to measure
cross-laboratory *chemistry* transfer:

- 26 rows on one structure that 3 existing experiments already hold
  (`AC_2025`, `JW_2024`, `BL_2023`), and
- 2 rows on DLin-MC3-DMA, which 11 existing experiments hold.

That makes it a **benchmark-only source** for every cross-source pair but one —
exactly the configuration `CLAUDE.md` §8 already has an open decision about,
there because benchmarks top every screen and can push concordance toward +1 for
reasons unrelated to transferable SAR. Here the benchmark does **not** top the
screen (LNP7 beats MC3 by 2.39× and 4.06×), which is the opposite hazard and
arguably more interesting, but it is still one structure against one benchmark.

The counter-argument: the corpus already contains formulation-variation rows on
a fixed ionizable lipid (`AC_2025` holds 216 cKK-E12 rows over 26 distinct
formulations), and this paper's formulations do **not** collide with any of
them.

STATUS: OPEN

---

### D-3 Should the eight non-LEC cell populations be emitted?

Tables S3 and S4 report the same aVHH+ readout for macrophages, monocytes, cDC2,
cDC1, double-negative cells, FRCs, BECs and VECs, at five sites, for MC3 and
LNP7 — 64 further cargo-expression measurements, classified inadmissible in
`05_values.csv` and counted there.

They are delivery measurements in kind. Two reasons they are held back:

- the corpus's `Model_target` vocabulary is organ-level (`liver`, `muscle`,
  `spleen`, `ear`, `brain`, `lung_epithelium`, `bone_marrow`, `whole_body`) and
  has no term for a FACS-sorted subset inside an organ, so emitting them would
  invent that convention on the strength of two ionizable lipids;
- one formulation would gain nine strongly correlated values per site inside a
  single z-block.

Against: the paper's most striking non-lymphatic result is that LNP7 reaches
~42% of cDC2 and ~28% of cDC1 in the draining node, which a reader of the corpus
would never see.

STATUS: OPEN

---

### D-4 Is a ratio of means an acceptable `Experiment_value` where the paper's own statistic is a mean of ratios?

The emitted value is `average aVHH+ cell count / average cell count` from SI
Tables S2–S4. The paper plots and quotes the mean of per-animal percentages.
Measured disagreement: Spearman 0.820641 over the 17 printed percentages, 10 of
17 within 3 percentage points, worst ratio 0.40 (Fig 2 BLN LNP1, 17.58% derived
against 44% printed) and 2.80 (Fig 4 PLN-injected MC3, 14.00% against 5%).

No per-animal data is published, so the paper's own statistic is not derivable
at any tier. The alternatives are: emit the ratio of means with the caveat on
every row (what is done), digitise the dot plots for the means (tier 5 on a
small raster, and it would still give only the 7 of 21 values whose panels we
could place), or emit nothing.

STATUS: OPEN

---

### D-5 `Model_target` needs five new terms, and `Experiment_method` is a spelling the corpus got wrong

The emitted rows introduce `lymph_node_axillary`, `lymph_node_brachial`,
`lymph_node_popliteal`, `lymph_node_popliteal_contralateral_nondraining`,
`lymph_node_axillary_nondraining`, `lymph_node_brachial_nondraining` and
`lymphatic_vessel_popliteal`. The `_nondraining` suffix is load-bearing, not
cosmetic: without it LNP7 would carry two axillary values 14.4× apart under one
key (see `08_contradictions.csv` C-5).

Separately, these rows use `Experiment_method = protein_adundance_normalized`,
which is how the corpus spells it — the misspelling is in `data_vivo.csv` today.
Matching it keeps the vocabulary closed; fixing it is a corpus-wide rename that
is not this folder's to make. Also note the readout is a **fraction of positive
cells**, not an abundance; it is the closest admissible class in the existing
vocabulary and `curation_note` says what it actually is.

STATUS: OPEN

---

### D-6 `Dose_ug_nucleicacid` is blank on all 28 rows

The paper states 1.5 mg/kg **only for the pooled SANDS screen**. The single-LNP
validation studies — which are the only source of every emitted value — state no
dose anywhere in the Methods, the figure captions or the SI. Three mutually
unsupported ways to fill it would be: carry the screen's 1.5 mg/kg across
(different experiment), take the 1 µg used in the VEGFC and Cy3 studies
(different cargo and different figure), or leave it blank.

Blank is what is emitted, following F-004's precedent.

STATUS: OPEN
