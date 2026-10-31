# ql_2026 — Liu et al., bispecific-linker lipidoid libraries

*"Synthesis and Screening of Bispecific Linker-based Combinatorial Libraries."*
Qiang Liu et al., **Adv Healthcare Mater 2026, 15, e02225**,
DOI `10.1002/adhm.202502225`. Not in PubMed as of 2026-10-29.

**344 datapoints → `QL_2026`. Built 2026-10-29. MERGED into `data/` 2026-10-30** (`build_dataset.py` `MERGE`; `provenance/CHANGELOG.md` §11). `data_vitro.csv` 28,541 → 28,885 rows, 41 → 42 experiments.

## The library

A one-pot three-component reaction (BL-3CR), catalyst-free, 80 °C / 4 h:

1. a primary amine attacks the linker's **NHS ester** → amide, losing NHS
2. an alkyl thiol adds across the linker's **maleimide** → thiosuccinimide
3. exception: **L2 is SPDP**, whose second site is 2-pyridyldithio, so step 2 is
   thiol-disulfide exchange → a **disulfide**, releasing pyridine-2-thione

**One arm per primary amine** (1 for 70/75/80/400, 2 for 113/114, 3 for
110/T3A), which the feed ratios in Tables S1/S2 encode exactly (1:1.05:1.1,
1:2.05:2.1, 1:3.05:3.1).

| library | lipids | composition | figure |
|---|---|---|---|
| 1 | **320** | 8 amines × 8 BLs × 5 linear thiols (C6–C14) | 2b heatmap, 2c bars |
| 2 | **24** | T3A × 8 BLs × 3 branched thiols (2C10/2C11/2C12) | 3c bars |

`2Cn` is a branched secondary thiol with **2n−1** carbons, not 2n: each arm is
n carbons counting the shared S-bearing methine, so 2C10 is C19H40S.

## Component SMILES — what the verification found

Claude chat supplied the 34 building blocks. All parse; all six reagents covered
by the corpus registry (C11) — DOPE, DSPC, DOPC, Cholesterol, DLin-MC3-DMA,
ALC-0315, SM-102 — match it exactly. `18PA` is used only in an in-vivo SORT
control, and the dioleoyl (18:1) assignment is right for that context.

**The mass check caught a real structural question.** Tables S1/S2 print MW and
observed m/z for 16 representative products. Against the first enumeration, four
rows failed badly: both L4 rows came out 28 Da (two CH2) light and both L5 rows
28 Da heavy. A swap of L4 and L5 explains that exactly, and was the obvious fix.

It is the wrong fix. Figure 2a **draws** L5 with a five-CH2 chain and L3 with
one, and the Results say "n-5CH2 alkyl chain linker (L5)". With the linkers as
drawn, all four rows match to **≤ 0.04 Da** once the tail is read as C8 rather
than the printed C6/C10. So the table mislabels those four tails, and the
chemistry is right. Fixing the half that was easier to change would have
renamed 80 lipids.

Final: **23 of 29 reference masses agree as printed, max error 1.45 Da** — all
24 rows of Table S1 and all 8 of Table S2 (the first version used only 16 and
silently dropped Table S2 entirely, which is the strongest evidence in the
paper: every one of its rows is T3A-Lx-**2C10**, so the tail is held constant
and the linker increments cannot be traded against tail length. All eight agree,
pinning L1–L8 on their own).

Of the six that do not: `113-L3-C6` is 14.07 Da light with no fit at all, a
paper error. The other five are the L4/L5 rows, and each admits **three**
mass-degenerate readings (moving a CH₂ between linker and tail changes neither
mass nor formula) — e.g. `70-L4-C6` fits `L3-C10`, `L4-C8` or `L5-C6` equally.
So which token is mistyped is **not decidable from mass**, and no single
explanation is asserted; the earlier text claiming "the table mislabels those
four tails" was over-confident. **The enumerated structures are identical under
every reading.** Table S1 also prints two rows as `400-L6-C6` and two as
`113-L6-C6`; the second of each matches the L7 product.

## Digitization

All four panels auto-detect their own geometry; every script writes a review
overlay to `review/`.

**Fig 2b, the 320-lipid heatmap** (continuous scale, so exact values, not bins).
The grid is read from the panel's own black rules — 41 vertical at 348–3298 and
9 horizontal at 1082–2022, pitch 73.78 × 117.62 — not assumed. Colour match to
the colorbar: median 1.0, max 11.3. The colorbar is linear, verified by
overlaying the predicted 2×10⁴ and 4×10⁴ positions on its printed ticks.

*Saturation is the hazard.* The colorbar stops at 6e4 while Fig 2c reaches
~4e5, so 21 cells clip to one blue. The first pass reported **zero** saturated
cells because the reference "top colour" was sampled 4 px into the bar and
picked up its black border, so nothing ever matched it. Corrected, 21 cells are
flagged; Fig 2c rescues 18 of them.

**Fig 3c, library 2** (log axis). Calibrated on the printed y-axis **tick
marks** — 2308.5 / 2459.5 / 2610.5 / 2760.5 / 2911.5 → **151.0 px per decade,
10³ at 2911.5**.

*This was wrong in the first version and Tier-2 review caught it.* The axis was
anchored on the panel **box bottom** (2959) as 10³, giving 163.75 px/decade. The
box top really does coincide with 10⁷, which is why it looked right — but the
10³ tick sits 47 px above the bottom rule, with bars running on below it. Every
Fig 3c value came out high by up to 1.76×. Three printed statistics catch it and
all three now agree exactly: Fig 3d's 2C11 rate (1/8), Fig 3e's L8 rate (0/3),
and the Results sentence "17 lipids … greater than 1 × 10⁵" (17, was 18). The
"boundary case" the first version reported at the 7×10⁵ threshold was not a
resolution limit — it was this error surfacing.

Separately, the three replicate **dots sit above the bar top** and, covering
~60–70% of the bar width, were initially mistaken for the top edge, inflating
two bars. Requiring a near-full-width rule fixes it.

**Fig 2d/2e/2f**, the printed hit-rate bars, are digitized only as a check. All
16 land on an exact k/40 or k/64, max error 0.49%.

## What validates what

`python3 -I validate.py` → 18 checks, all passing (`13_validation.json`).

- **Fig 2c is the strongest check on the heatmap.** It plots the three-tail
  lipids above 2e4, and 32 bar runs were detected: 8 for amine 110, 20 for T3A,
  4 benchmarks. Every one of the 8 lipids it names for amine 110 is exactly a
  cell the heatmap puts above 2e4 — 8 of 8 — and 19 of its 20 T3A lipids
  likewise. On the lipids both figures give an uncensored number for, the two
  independent readings (continuous colour vs log bar height) agree to a **median
  ratio of 1.029**.
- **Fig 3c reproduces every printed 3d and 3e hit rate exactly**, plus the
  Results' "17 lipids … greater than 1 × 10⁵". All three are now asserted in
  `validate.py`; the earlier version excused the 2C11 rate as a boundary case,
  and that was the one datum that would have caught the calibration error.
- **C10 anchor**: `T3A-L5-C10` written by hand from the drawn blocks with a
  different ring traversal and atom order — `exact_match`.

### Two inconsistencies inside the paper, recorded not reconciled

1. **Fig 2d/2e do not follow from Fig 2b.** They imply 36 hits above 2e4, of
   which 25 are three-tail — but the paper's own Fig 2c plots **28** three-tail
   lipids, and the heatmap gives 27. No single threshold reproduces the printed
   per-amine counts, and the disagreement is not monotone (2d has more for
   amine 75 and fewer for 80 than the heatmap), so it cannot be calibration.
   Fig 2f independently sums to 37 where 2d and 2e sum to 36.
2. **T3A-L8-C14** is plotted in Fig 2c at 5.8e4 but its Fig 2b cell is visibly
   cream and reads 7.9e3 — a 7.3× conflict between two figures of one paper.
   Flagged `value_conflict`; the heatmap value is carried, since the heatmap is
   the complete screen.

## Rows

`11_QL_2026_rows.csv`, 344 × 72. Value per lipid: Fig 2c where it plots the
lipid (27), else Fig 2b (293); Fig 3c for library 2 (24). The three-component
decomposition (`IL_head`, `IL_linker`, `IL_tail1`) is populated for all 344.

**One pooled z-score, after a scalar offset.** The four benchmark LNPs appear in
both Fig 2c and Fig 3c, giving a cross-figure ratio of **1.112× (sd 0.023)**.
Library 2 is divided by that and all 344 are z-scored **once, together**.

An earlier version standardised the two libraries separately, citing a 1.54×
offset. Both halves were wrong. The 1.54 came from the miscalibrated Fig 3c
axis. And separate standardisation does not merely remove an offset — it forces
library 2's spread onto library 1's and deletes the gap between them, which is
the paper's central claim (geometric means **4,128 vs 154,630 RLU**, ~37×).
It inverted the raw ordering of 36% of cross-library pairs while every automated
check still passed, because per-block mean 0 / sd 1 is a tautology for
concatenated standardised blocks. After the fix, library 2's mean z is **+2.12**
against library 1's −0.16, i.e. the finding survives the scaling.

Flags:

| flag | n | meaning |
|---|---|---|
| `value_censored` | 3 | saturated in Fig 2b, absent from Fig 2c — lower bound only |
| `value_floored` | 113 | at the colorbar floor, where the colormap is degenerate: runs of identical RGB span up to ~980 RLU, so the value is an upper bound within its band |
| `value_conflict` | 1 | `T3A-L8-C14`, 7.3× apart between Fig 2b and Fig 2c |

The floored flag was missing in the first version, which mattered: 113 of 320
library-1 cells sit in the degenerate cream band, the 320 cells take only 45
distinct values, and 170 are at ≤3,000 RLU. `T3A-L8-C14` now carries the
**heatmap** value, as the text always said — the first version's code silently
took the Fig 2c value instead.

Conditions, from the Methods: HeLa, 1×10⁴ cells/well in 96-well; mFluc mRNA
(1-methylpseudouridine) **20 ng/well**; 24 h; n = 3; lipid/DOPE/Chol/DMG-PEG2000
= **35/16/46.5/2.5 mol%**; citrate pH 4.0; pipette mixing (~300×). The screened
material was **crude** — blown down under N₂ and used directly; only the in-vivo
lipids were flash-purified to >95%.

`IL_to_nucleicacid_massratio` is left **blank**. The paper's "20:1" is called
"lipid/mRNA" in the Results and "LNP/mRNA" in the Methods and is never defined
as ionizable-lipid-to-mRNA. If it is total lipid, the true IL:mRNA runs 4.4–12.5
across this library, since the IL mass fraction varies 0.22–0.63 with lipid MW.
Writing 20.0 would assert something the paper does not say.

`PEG_name` is written as **DMG-PEG2000**, following the main text, Methods and
every figure. The SI Materials list instead names a phosphoethanolamine-PEG2000
from Sigma. The corpus distinguishes the two (DMG-PEG2000 21,392 rows,
DMPE-PEG2000 4,901), so the contradiction is recorded here rather than resolved.

## Pre-merge blockers, and how each was resolved

Merged 2026-10-30. Every item below was open at build time; this section records
the disposition rather than deleting the question.

- ~~`Publication_PMID` is blank (not yet indexed).~~ **Wrong premise.** PMID
  **41145968** was already indexed, with a matching DOI
  (10.1002/adhm.202502225) and first author Liu Q. Verified against NCBI
  esummary 2026-10-29 and populated on all 344 rows. *Look the PMID up before
  recording a paper as unindexed.*
- ~~Amine **114**'s degree of substitution (2 arms) ... cannot be resolved from
  this paper.~~ **Resolved; it is pinned three independent ways**, found by the
  round-1 review:
  1. the SI Materials list names it "1,5-diamino-3-azapentane (114, 97%)", whose
     hyperlink tooltip in `word/document.xml` gives the synonyms *diethylene-
     triamine, DETA, bis(2-aminoethyl)amine* — H₂N-CH₂CH₂-NH-CH₂CH₂-NH₂, exactly
     **two primary amines**;
  2. Fig 2a draws it with two labelled NH₂ and an explicit central N–H;
  3. **Table S1 supplies the control that settles it.** `400-L6-C6` is printed at
     471.7 Da; amine 400 has one primary NH₂ *and* one secondary NH, and a
     one-arm build gives 472.70 (+1.0) where two arms would give 745.0. So a
     secondary NH does not acylate under these conditions, and 114's two primary
     amines give exactly two arms.

  The 40 affected structures are correct. This folder previously understated its
  own evidence — the Table S1 control was already in the mass check and was not
  recognised as answering the question.
- **`T3A-L8-C14` is kept, carrying the Fig 2b heatmap value (7,931 RLU).**
  Decision recorded 2026-10-30, having been left open at build time. Reasoning:
  it keeps all 320 library-1 values on one uniform source, and this would
  otherwise be the only row where a contradicted source is preferred. The row is
  flagged `value_conflict` in `12_QL_2026_qc.json` — though see the open item
  below, because that flag does not reach `data/`.
- `HeLa` already exists in the corpus; no new `Model_type`. Confirmed (9,607
  pre-existing HeLa rows across 15 experiments).

## Open after merging — round-1 review findings

Recorded in `provenance/review/round_2026-10-30/FINDINGS.md`; none repaired yet.

- **117 of 344 rows (34%) carry a flag that never reaches `data/`** — 113
  floored, 3 censored, 1 conflict. The 72-column schema has no flag column, so
  the policy in `CLAUDE.md` §5/§6 is unimplementable as written. Project-wide
  (1,137 rows), not specific to this paper.
- **The floored value is above the feasible range.** The colormap is pure white
  from 0 to ≈900 RLU, so 48 cells are *upper bounds*, not measurements, and the
  assigned 1,310 RLU exceeds the whole band. Varying the assumption moves
  individual z-scores by up to ~1.0 — the dominant uncertainty in this slice,
  an order of magnitude larger than every digitization error below.
- **Fig 2c is still calibrated on y-axis label centroids, not tick marks.** The
  stored anchors give decade gaps of 406 and 408 px — unequal, which is the
  proof they are text. The ticks give exactly 203.0/203.0 native = 406.0. Every
  Fig 2c value is ≈1.7% high. **This is class 20, the same error repaired in
  Fig 3c, left live in the sibling panel of the same curation.**
- **Replicate dots are still read as the bar top on 15 of 60 bars** across
  Fig 2c and Fig 3c (up to +26%), because the `w//8` inner trim lets three
  overlapping dots exceed 85% of the *trimmed* width. This biases the
  cross-library offset from ≈1.073 to the stored 1.1116. Tell: two pairs of bars
  carry byte-identical stored values.
- Net effect of the digitization errors, measured: rms Δz 0.054, max 0.209,
  Spearman 0.99998. Real and worth fixing, but small beside the floor question.
- `113-L3-C6`'s "no fit at all, a paper error" is **wrong**: `114-L3-C6` =
  613.84 Da against the printed 613.8 is a unique fit. The degeneracy search
  moved a CH₂ only between linker and tail, never into the amine core.

## Review

Tier 2 adversarial review ran 2026-10-29 and found the Fig 3c calibration error
(E1), the conflict-row code/text mismatch (E2), the z-scoring problem (E3), the
truncated mass check (E4) and the missing floored flag (E7); all are fixed
above. It also recorded paper-internal errors worth knowing: the Results discuss
"C11 or C12 tails" for library 1, which has neither; the Methods cite "Tables
S1, S3" for feed ratios that are in S1 and S2; and the Fig 2 caption's n = 40
holds only for panel 2d (2e uses a flat denominator of 40 over groups of
160/80/80, 2f uses 64), while Fig 3's n = 8 holds for 3d but 3e's groups hold 3.

## Files

```
digitize_fig2b.py     heatmap, 320 -> 01_fig2b_320.csv, 02_fig2b_qc.json
digitize_fig2def.py   printed hit-rate bars (check only) -> 03_, 04_
digitize_fig3c.py     library 2 + benchmarks -> 05_, 06_
digitize_fig2c.py     three-tail lipids >2e4 + benchmarks -> 07_, 08_
enumerate_library.py  344 structures + mass check -> 09_, 10_
build_rows.py         corpus schema -> 11_QL_2026_rows.csv, 12_QL_2026_qc.json
validate.py           18 independent checks -> 13_validation.json
reference_structures.csv   C10 anchor
review/               overlay PNGs for every extraction
```
