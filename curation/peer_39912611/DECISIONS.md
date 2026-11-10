# DECISIONS -- peer_39912611 (Naidu / Peer, ACS Nano 2025, PMID 39912611)

Opened 2026-11-10. The pipeline detects that a source is silent; it does not
get to decide what the source meant. `provenance/build_dataset.py` refuses a
folder that still has an `OPEN` entry, which is the intended state here: this
curation is NOT merged.

---

### D-1 The in-vitro dose series: 0.5/0.25/0.125 (Methods) or 0.25/0.15/0.025 (Figure 2)?

Methods, *LNP In Vitro Testing*: "treated with LNPs at an mLuc mRNA
concentration of 0.5 μg/mL, 0.25 μg/mL, and 0.125 μg/mL". Figure 2's x-axis
labels, under every one of the 36 triplets in all three panels: `0.25 ug/ml`,
`0.15 ul/ml`, `0.025 ug/ml`. The two series share only 0.25, and they are not
related by any single factor. The Figure 2 caption anchors its own statistics
to "the 0.25 μg/mL treatment with lipid 14", i.e. the leftmost and highest bar,
which agrees with the figure and not with the Methods. The middle label also
carries a unit typo, `ul/ml`.

What we did: carried the FIGURE labels on every row, with the unit string
`ug/mL (AXIS LABEL -- contradicts Methods, see D-2)` beside them, and emitted
all three doses rather than picking one.

One piece of evidence nobody has used yet: check `F2_2` in `07_qc.json`
identifies which dose the authors' own correlation panel used, by reproducing
the printed Pearson r against the exact pKa vector. The MIDDLE dose wins for
both cell lines that have a printed r (HeLa 0.7364 vs printed 0.7358; HepG2
0.7191 vs printed 0.7226), against 0.7534/0.6959 for the first dose and
0.5592/0.7172 for the third. That tells us which bar they used; it does not
tell us what concentration that bar was.

Not resolvable from the paper. A one-line question to the corresponding author
would settle it.

STATUS: OPEN

---

### D-2 Which `Experiment_ID`, and should the two pools share one?

Proposed `GN_2025` (first author Gonna Somu Naidu, 2025); no collision with the
54 IDs now in `data/`. Both pools are proposed under the same ID because they
are one screen of one library by one lab. The in-vitro pool then spans three
`Model_type`s and the in-vivo pool three `Model_target`s, which makes
`Experiment_value_within_model` non-trivial and is exactly the F-016 / F-028
territory: this curation leaves that column BLANK rather than compute a
within-model z whose unit it has not declared.

STATUS: OPEN

---

### D-3 The six upper-bound rows

Five lipids have no liver bar in Figure 3D (34, 35, 36, 41, 42) and one has no
lung bar in Figure 3C (39). They are not missing data: Figure 3F plots
log10(lungs/liver) for all twelve lipids including 41 and 42, so the liver
numbers exist and are simply below what panel D's axis can show. We store the
axis floor with `value_is_bound = upper` and `confidence = low`.

Two options the merge has to choose between: keep them as bounds (a model that
treats them as measurements is fitting the floor of an axis), or recover them
from the ratio panels. Recovery is possible and was NOT done here -- Figure
3F/G/H are linear-axis panels with bars running both ways from zero, and
extracting them is a second digitiser this curation did not write. Check `X2`
quantifies what the censoring costs: substituting the floor for lipid 39's lung
value is what makes our in-vitro-to-lung correlations land 0.03-0.08 below the
authors' printed ones, and putting its true value at 1e4-3e4 reproduces both.

STATUS: OPEN

---

### D-4 The Untreated control bar is not measured

Figure 3C/D/E each carry an `Untreated` bar. It is drawn HOLLOW: white fill,
black outline, and its top rule is one pixel. `curation/lib/bars.bar_top`
requires its full-width rule to hold for three contiguous rows, and the only
way to satisfy that for a hollow bar is to widen the window onto the two
outline columns it shares with its neighbours -- which puts a neighbour's bar
inside the flank test and was measured to break two other categories when
tried. All three are recorded as refusals. No corpus row depends on them; they
are a control, not a library member. If the background radiance is wanted, it
needs a detector this library does not have.

STATUS: OPEN

---

### D-5 Stereochemistry of the 2-hexyldecanoate centres

Every one of the twelve lipids has two sp3 stereocentres, at C2 of each
2-hexyldecanoate tail. The paper draws them flat in Supplementary Figure 15 and
in every per-lipid scheme, names them without a descriptor, and buys the
building block from Sigma -- i.e. racemic. Both the enumeration and the C10
anchors are emitted WITHOUT stereochemistry, so nothing asserts a configuration
the source does not state. This is a deliberate absence, not an oversight: 212
corpus rows already carry an unassigned stereocentre condition and the merge
should decide whether these 12 lipids join that baseline.

STATUS: OPEN

---

### D-6 `value_normalization_unit` and the two multi-model pools

`Experiment_value` is emitted as a per-experiment z of log10 with population sd
(`value_normalization_unit = experiment`), which is the §3 convention. For the
in-vitro pool that pools three cell lines into one z; for the in-vivo pool it
pools three organs. Whether that is right for this paper is the F-016 question
(standardising per cell type deleted a real between-assay offset in three other
experiments) in its mirror image -- here we pool, and the risk runs the other
way. The per-model split is recoverable from `05_values.csv` at any time, so
nothing is lost either way; it just has to be chosen deliberately.

STATUS: OPEN
