# Re-measurement of three digitization defects — QL_2026 (Liu et al., Adv Healthc Mater 2026, 15, e02225)

Round 2026-10-30, independent re-measurement. Read-only against the repository: no
repository file was edited and no curation output was regenerated.

## Method common to all three findings

**Geometry re-derived from the PDF, not from the stored page PNGs.** I re-rendered
pages 4 and 5 of
`curation/ql_2026/Adv Healthcare Materials - 2025 - Liu - Synthesis and Screening of Bispecific Linker‐based Combinatorial Libraries.pdf`
with PyMuPDF at **600 dpi**, clipped to the embedded figure image's bounding box
(page 4 xref 120, bbox 58.4,76.9–533.2,658.8; page 5 xref 127, bbox 53.3,76.7–544.4,466.7).
The results are **bit-identical** to `00_pages/fig_p4_600.png` (3958×4851) and
`00_pages/fig_p5_600.png` (4093×3251) — `max |diff| = 0`. So every pixel coordinate below
is directly usable as a constant in the existing scripts, with no frame conversion.

**Resolution floor.** The figures are embedded **JPEGs at native 1976×2423 (p4) and
2046×1626 (p5)**, i.e. ≈300 dpi. The "600 dpi" page is a 2.003× upsample of that JPEG, so
every rendered row comes in identical pairs and the true resolution is **1 native pixel =
2 rendered pixels**. All uncertainties below are quoted in 600-dpi pixels; ±2 px is the
hard floor, ±4 px is realistic on a soft edge. Tick detection was therefore done on the
**native JPEG** and converted, which is where the sub-pixel precision comes from.

Overlays are in
`/private/tmp/claude-501/-Users-nathanliu-Documents-LNP-confound/ebd1fdd9-2e2e-4284-8eb1-27c940f1a93a/scratchpad/overlays/`
(six PNGs, named per finding). They are outside the repo because I am read-only here; copy
them in if they should be archived.

---

# F-017 — "Fig 2c's log axis is anchored on text centroids, not tick marks"

## VERDICT: **partly** — the diagnosis is right, the consequence is not.

The three integers in the comment on `digitize_fig2c.py:35` are **not** the tick positions,
and the "406 / 408 unequal gaps" story is real as a provenance defect. But the two constants
that the code actually uses are **correct**: `PX_PER_DECADE = 407.0` is right to 0.01 px and
`Y_1E4 = 3156.0` is low by **0.75 px**, worth **+0.42 %** on every Fig 2c value. The
previous reviewer's proposed replacements (406.0 and 3153.0) are **worse than what is
stored** — I can refute them directly.

## MY MEASUREMENT

**Page 4 (Figure 2), 600 dpi, clipped-to-figure frame. Ticks measured on the native 300-dpi
JPEG and scaled by 4851/2423 = 2.002064.**

Fig 2c's y-axis rule occupies x = 515–530 (600 dpi) / 257–264 (native). The ticks are black
stubs protruding **left** of the rule. Scanning x = 250–256 native over y = 1100–1660 with an
ink threshold finds **exactly three** ink bands and nothing else:

| tick | native y band | native centroid | 600-dpi y |
|---|---|---|---|
| 1×10⁶ | 1167–1174 | 1170.180 | **2342.78** |
| 1×10⁵ | 1370–1377 | 1373.353 | **2749.54** |
| 1×10⁴ | 1572–1581 | 1576.749 | **3156.75** |

Sub-pixel centroids are the mean over five independent x-windows (sd 0.054 / 0.005 / 0.014
native px). **There are no minor ticks** — a sensitive sweep (ink > 25/255 on the columns
immediately left of the rule, y 1100–1660 native) returns only these three bands plus the
L1-group bracket at y 1131–1135, which is above the panel. So the answer to "major and
minor" is: three major ticks, zero minor ticks, all three resolvable.

**Decade spacing implied by each adjacent pair (native px):**

| pair | centroid method | cross-correlation method |
|---|---|---|
| 10⁶→10⁵ | 203.173 | 203.215 |
| 10⁵→10⁴ | 203.396 | 203.431 |

The two gaps differ by **0.22 native px = 0.43 px at 600 dpi** — i.e. the renderer rounded
tick placement to the nearest device pixel. That is equal spacing to within the figure's own
resolution. The cross-correlation estimate is fully independent of the centroid estimate
(it shifts one tick's intensity profile onto the next on a 0.001-px grid; peak correlations
0.996 and 0.988) and agrees to 0.04 px.

**Linear fit of log10(value) against tick y (3 points, 1 dof):**

- slope = −0.00491922 decade/native px → **203.284 native px/decade = 406.99 px/decade at 600 dpi**
- residuals = +0.000183, −0.000365, +0.000183 decades; **rms 0.000258 decades**
- the same residuals in pixels: −0.037, +0.074, −0.037 native px

So the axis is linear in log10 to 2.6×10⁻⁴ decades (0.06 % in value). The ticks are
resolvable and the fit is essentially exact.

**Recommended constants (600-dpi clipped-figure frame, page 4):**

```
Y_1E4         = 3156.75     (was 3156.0;  +0.75 px)
PX_PER_DECADE = 407.0       (unchanged;   measured 406.99)
```

**Independent validation — the red threshold line.** Fig 2c carries a red dotted horizontal
line. Its intensity-weighted centroid is y = **2872.21** (rows 2867–2876). Under the three
competing calibrations:

| calibration | value of the red line | error vs a round 5×10⁴ |
|---|---|---|
| **mine** (3156.75, 406.99) | **50,020 RLU** | **+0.04 %** |
| stored (3156.0, 407.0) | 49,804 RLU | −0.39 % |
| prior reviewer (3153.0, 406.0) | 49,159 RLU | −1.68 % |

The line is obviously drawn at exactly 5×10⁴ (a round number, and the only horizontal
reference in the panel). My calibration hits it to 0.04 %. The prior reviewer's numbers
are 40× further off. **The "203.0 / 203.0 native = 406.0, 10⁴ at y = 3153.0" report is
wrong**; the true values are 203.17/203.40 native and 3156.75.

**What the three integers in the comment actually are.** I measured the y-axis *label* ink,
x = 120–240 native:

| label | ink centroid (600 dpi) | bbox midpoint (600 dpi) | tick (600 dpi) | comment says |
|---|---|---|---|---|
| 1×10⁶ | 2344.3 | 2340.4 | 2342.78 | 2342 |
| 1×10⁵ | 2750.8 | 2747.8 | 2749.54 | 2748 |
| 1×10⁴ | 3158.9 | 3155.3 | 3156.75 | 3156 |

The stored triple (2342, 2748, 3156) is not any one of these three consistently — it reads
like rounded-down label positions. So **the comment is untrustworthy and should be
rewritten** (that is an error-class-21 docstring defect, same genus as the one already
recorded for `digitize_fig3c.py`). But the *mean* of the stored anchors, 407.0, happens to
be the correct decade spacing, and the stored `Y_1E4` happens to be within 0.75 px.

**I also re-derived Fig 3c's axis the same way, to check the sibling repair.** Page 5,
native x = 102–114 left of the rule at native 115–120, five tick bands:

| tick | native centroid | 600-dpi y |
|---|---|---|
| 1×10⁷ | 1154.008 | 2307.31 |
| 1×10⁶ | 1229.870 | 2458.98 |
| 1×10⁵ | 1305.615 | 2610.43 |
| 1×10⁴ | 1380.732 | 2760.62 |
| 1×10³ | 1455.973 | 2911.05 |

Gaps 75.862 / 75.745 / 75.117 / 75.240 native; least-squares fit gives **150.913 px/decade
at 600 dpi** and **y(10³) = 2911.50**. Stored: 151.0 and 2911.5. **Fig 3c's calibration is
confirmed exactly** (0.06 % per decade, 0.00 px on the anchor); the sibling repair was
done correctly and needs no change. Fit residuals ≤0.375 native px (≤0.5 % in value).
Fig 3c's own red dotted line sits at y = 2606.17 → 1.055×10⁵, i.e. ~5 % above a nominal
1×10⁵; that is a drawn-by-eye reference line, not a calibration problem, since the five
ticks fit a straight line to 0.003 decades.

**Overlays:**
`…/scratchpad/overlays/F017_fig2c_yaxis_ticks.png` — red rules at my three tick centroids
(they land on the stubs), blue rules at the label ink centroids (visibly offset), green rule
at the red-dotted-line centroid.
`…/scratchpad/overlays/F017_fig3c_yaxis_ticks.png` — the five Fig 3c ticks.

## WHAT CHANGES

One constant, one line:

- `curation/ql_2026/digitize_fig2c.py:35` — `Y_1E4 = 3156.0` → `Y_1E4 = 3156.75`
  and replace the comment with
  `# y-axis TICK centroids (native 300-dpi JPEG, x 250..256, scaled x2.002064): 1e6 2342.78, 1e5 2749.54, 1e4 3156.75; gaps 406.76/407.21 px; linear-fit residual 2.6e-4 decades`
- `curation/ql_2026/digitize_fig2c.py:36` — `PX_PER_DECADE = 407.0` — **no change**
  (measured 406.99).
- `curation/ql_2026/digitize_fig3c.py:61–62` — **no change** (measured 2911.50 / 150.91).

Effect on the dataset, in isolation: every one of the 27 library-1 rows whose value comes
from Fig 2c rises by exactly ×1.00425 (+0.42 %). Over the 344-row slice that is
**rms Δz = 0.0006, max |Δz| = 0.0033**. The cross-figure offset moves 1.1116 → 1.1068.
No row's flag, rank or censoring status changes.

## CONFIDENCE: **certain** on the measurement (three independent estimators agree to 0.1 px
and an external round-number check confirms to 0.04 %); **certain** that the prior
reviewer's 406.0 / 3153.0 is wrong.

---

# F-018 — "the bar-top detector reads replicate dots as the bar top"

## VERDICT: **confirmed** on the mechanism and on the count; **refuted** on the corrected
cross-library offset (it is ~1.01, not 1.073) and on the claimed cause of the byte-identical
pairs being anything other than pixel coincidence on a contaminated read.

## MY MEASUREMENT

**Dpi/page:** Fig 2c = page 4, Fig 3c = page 5, both 600 dpi in the clipped-figure frame.

**Why the stored rule fails — observed directly.** Every bar is a *hollow* rectangle with
three filled replicate dots and an error whisker drawn on top in the same blue. On
`T3A-L3-C6` the stored read of y = 2713 lands where the **left dot (x 1818–1840) plus the
whisker cap (x 1848–1870)** together cover ≈85 % of the `w//8`-trimmed window
(1818–1872 of a 1809–1881 bar) while covering only ≈55 % of the full bar. The true top rule
is at y = 2755. That is exactly the failure mode described, and it is visible at 5× zoom.

**My method (dot-immune, and it is a flank test).** For each bar:

1. x-extent from the blue column profile (independent of the stored runs; 32 runs in 2c,
   28 in 3c, matching the labels).
2. The two vertical side strokes are located from a 18-px row band just above the baseline,
   where the hollow interior is unambiguously white; their centres `cL`, `cR` are recorded.
   A centred marker cluster cannot reach both.
3. Bar top = the first row, scanning down, at which a **single contiguous ink run**
   (gray ≤ 172) inside `[x0,x1]` (a) is at least 72 % of the bar width and (b) **contains
   both `cL` and `cR`** — i.e. the left flank and the right flank are joined by one
   unbroken rule — sustained for ≥2 rows (= 1 native row). The frac relaxes to 0.66/0.60
   only for the three narrowest 3c bars.

Condition (b) is the flank-agreement test. A merged dot cluster fails it because it never
reaches both outer strokes at once; it is why `T3A-L2-2C10` is read at 2519 (rule) and not
2511 (dots spanning 91 % of the width but offset 4 px right of the left stroke).

**Flank agreement, per bar.** Measuring the left and right stroke tops separately
(topmost row at which the whole stroke core is inked) and comparing:

- Fig 2c: left and right flank tops agree to **0 px on 20 of 32 bars**; they disagree on 12,
  always because a dot overlaps one stroke, and in every such case the *lower* of the two
  equals the joined-rule answer to ≤2 px. Largest disagreements: `T3A-L1-C10` 12 px,
  `110-L2-C10` 20 px, `110-L8-C10` 18 px, `T3A-L7-C14` 14 px, `DLin-MC3-DMA` 24 px.
- Fig 3c: agree to 0 px on 10 of 28, to ≤2 px on 18 of 28. Worst: `T3A-L8-2C11` 16 px,
  `T3A-L2-2C11` 14 px, `T3A-L7-2C11` 14 px — all one-sided dot overlap.
- Two bars where a naive flank rule disagrees with the joined rule by more than 2 px and the
  joined rule is right: `T3A-L5-2C12` (right stroke's outer 2 anti-aliased columns only
  become blue 22 px lower; joined rule 2569, confirmed by eye) and `ALC-0315` in Fig 3c
  (same, joined rule 2581, confirmed by eye). **Flank-only reading is not safe on this
  figure; the flanks must be required to be *joined*.**

**Cross-check that this is not just my detector.** On bars with no marker over the top edge,
my y and the stored y agree **to the pixel**: 9 bars at dY = 0 in Fig 2c and 5 in Fig 3c, and
a further 11 and 11 at dY = ±2 px (one native pixel, the resolution floor). That isolates
the cause to marker contamination rather than a different convention. All 60 bar tops were
checked against the overlays below and, for the 12 contested ones, against printed
intensity maps.

**Overlays:**
`…/scratchpad/overlays/F018_fig2c_bartops.png` and `…/scratchpad/overlays/F018_fig3c_bartops.png`
— red line/box = my bar top, orange line = the stored `y_top`, green rules = my tick
calibration. Where orange floats above red, the stored detector was on the dots.

### Per-bar results — Fig 2c (`Y_1E4 = 3156.75`, 407.0 px/decade)

| lipid | stored RLU | my RLU | ratio | dY px | marker over top |
|---|---|---|---|---|---|
| 110-L1-C10 | 27,530 | 27,336 | 0.993 | +2 | no |
| 110-L2-C8 | 68,840 | 68,354 | 0.993 | +2 | yes |
| 110-L2-C10 | 105,821 | 106,271 | 1.004 | 0 | yes |
| 110-L3-C6 | 85,350 | 85,713 | 1.004 | 0 | no |
| 110-L4-C8 | 27,530 | 27,647 | 1.004 | 0 | no |
| **110-L7-C10** | **36,119** | **28,927** | **0.801** | **+40** | no (whisker cap + dot) |
| 110-L8-C8 | 35,311 | 35,062 | 0.993 | +2 | yes |
| **110-L8-C10** | **58,755** | **54,511** | **0.928** | **+14** | yes |
| T3A-L1-C10 | 103,453 | 103,893 | 1.004 | 0 | no |
| T3A-L2-C8 | 42,801 | 42,499 | 0.993 | +2 | no |
| T3A-L2-C10 | 73,675 | 73,156 | 0.993 | +2 | no |
| T3A-L2-C12 | 61,475 | 61,736 | 1.004 | 0 | no |
| **T3A-L3-C6** | **122,589** | **97,074** | **0.792** | **+42** | yes |
| T3A-L3-C10 | 32,993 | 32,761 | 0.993 | +2 | no |
| T3A-L4-C8 | 199,413 | 198,008 | 0.993 | +2 | no |
| **T3A-L4-C10** | **109,474** | **99,295** | **0.907** | **+18** | yes |
| **T3A-L4-C12** | **170,200** | **159,704** | **0.938** | **+12** | yes |
| **T3A-L5-C8** | **24,864** | **23,597** | **0.949** | **+10** | yes |
| **T3A-L5-C10** | **44,782** | **37,952** | **0.847** | **+30** | yes |
| **T3A-L5-C14** | **24,308** | **22,809** | **0.938** | **+12** | yes |
| T3A-L6-C14 | 27,530 | 27,647 | 1.004 | 0 | no |
| T3A-L7-C10 | 47,389 | 47,590 | 1.004 | 0 | no |
| T3A-L7-C12 | 384,385 | 377,382 | 0.982 | +4 | no |
| **T3A-L7-C14** | **82,502** | **78,295** | **0.949** | **+10** | yes |
| T3A-L8-C8 | 91,346 | 89,681 | 0.982 | +4 | no |
| T3A-L8-C10 | 157,239 | 156,131 | 0.993 | +2 | no |
| T3A-L8-C12 | 132,694 | 131,759 | 0.993 | +2 | no |
| T3A-L8-C14 | 58,094 | 57,684 | 0.993 | +2 | no |
| DLin-MC3-DMA † | 48,473 | 50,933 | 1.051 | −8 | yes |
| SM102 † | 48,473 | 49,794 | 1.027 | −4 | yes |
| **ALC-0315** | **182,156** | **159,704** | **0.877** | **+24** | yes |
| Lipid 5 | 73,675 | 73,156 | 0.993 | +2 | yes |

† `DLin-MC3-DMA` and `SM102` are the only two bars in either panel whose top rule is
*crossed by the red 5×10⁴ dotted line*, which adds ink above the rule and pulls my reading
up by 4–8 px. For these two the honest statement is **y_top = 2873 ± 4** for both; both sit
on the 5×10⁴ line to within the figure's resolution. I report both variants below.

### Per-bar results — Fig 3c (`Y_1E3 = 2911.5`, 151.0 px/decade)

| lipid | stored RLU | my RLU | ratio | dY px | marker over top |
|---|---|---|---|---|---|
| T3A-L1-2C10 | 933,681 | 905,636 | 0.970 | +2 | no |
| **T3A-L2-2C10** | **492,100** | **397,501** | **0.808** | **+14** | yes |
| **T3A-L3-2C10** | **284,212** | **222,680** | **0.784** | **+16** | yes |
| **T3A-L4-2C10** | **149,795** | **124,746** | **0.833** | **+12** | yes |
| T3A-L5-2C10 | 1,346,291 | 1,305,852 | 0.970 | +2 | yes |
| T3A-L6-2C10 | 136,698 | 136,698 | 1.000 | 0 | no |
| **T3A-L7-2C10** | **933,681** | **826,454** | **0.885** | **+8** | yes |
| T3A-L8-2C10 | 293,013 | 284,212 | 0.970 | +2 | no |
| T3A-L1-2C11 | 236,686 | 222,680 | 0.941 | +4 | yes |
| T3A-L2-2C11 | 33,612 | 32,602 | 0.970 | +2 | yes |
| T3A-L3-2C11 | 97,739 | 89,193 | 0.913 | +6 | no |
| T3A-L4-2C11 | 5,560 | 5,230 | 0.941 | +4 | no |
| T3A-L5-2C11 | 731,541 | 731,541 | 1.000 | 0 | yes |
| T3A-L6-2C11 | 435,585 | 435,585 | 1.000 | 0 | yes |
| T3A-L7-2C11 | 174,470 | 174,470 | 1.000 | 0 | yes |
| T3A-L8-2C11 | 590,912 | 573,163 | 0.970 | +2 | yes |
| T3A-L1-2C12 | 839 | 814 | 0.970 | +2 | yes |
| T3A-L2-2C12 | 13,463 | 13,059 | 0.970 | +2 | yes |
| T3A-L3-2C12 | 51,513 | 51,513 | 1.000 | 0 | no |
| T3A-L4-2C12 | 12,286 | 11,559 | 0.941 | +4 | no |
| T3A-L5-2C12 | 191,186 | 185,444 | 0.970 | +2 | yes |
| **T3A-L6-2C12** | **1,346,291** | **1,155,883** | **0.859** | **+10** | yes |
| T3A-L7-2C12 | 236,686 | 236,686 | 1.000 | 0 | yes |
| T3A-L8-2C12 | 573,163 | 555,947 | 0.970 | +2 | yes |
| DLin-MC3-DMA | 54,753 | 51,513 | 0.941 | +4 | no |
| SM102 | 53,109 | 51,513 | 0.970 | +2 | no |
| **ALC-0315** | **197,107** | **154,434** | **0.784** | **+16** | yes |
| **Lipid 5** | **83,915** | **74,278** | **0.885** | **+8** | yes |

### Counts

- **|dY| ≥ 10 px (stored high by ≥ 5.8 % in 2c, ≥ 15 % in 3c): exactly 15 of 60 bars.**
  Fig 2c (10): 110-L7-C10, 110-L8-C10, T3A-L3-C6, T3A-L4-C10, T3A-L4-C12, T3A-L5-C8,
  T3A-L5-C10, T3A-L5-C14, T3A-L7-C14, ALC-0315.
  Fig 3c (5): T3A-L2-2C10, T3A-L3-2C10, T3A-L4-2C10, T3A-L6-2C12, ALC-0315.
  **The "15 of 60" claim reproduces exactly.**
- |dY| ≥ 6 px: 19 of 60. |dY| ≤ 4 px (within resolution): 41 of 60, of which 14 at exactly 0.
- Reported worst cases: `T3A-L3-C6` claimed +26 %, I measure **stored/true = 1.263 (+26 %)** ✓;
  `110-L7-C10` claimed +25 %, I measure **1.249 (+25 %)** ✓; `T3A-L3-2C10` claimed +24 %,
  I measure **1.276 (+28 %)** — the right bar, slightly larger than claimed.
- 26 of the 60 bars have a replicate marker over the top edge (≥30 % ink coverage in the
  20 px band above the rule), but only 15 of those were mis-read: a marker overlapping the
  rule is necessary, not sufficient.

### 3. Benchmark ratios and the cross-library offset

| benchmark | Fig 3c / Fig 2c, stored | Fig 3c / Fig 2c, re-measured | re-measured, red-line variant † |
|---|---|---|---|
| DLin-MC3-DMA | 1.1296 | 1.0113 | 1.0582 |
| SM102 | 1.0956 | 1.0345 | 1.0582 |
| ALC-0315 | 1.0821 | 0.9669 | 0.9669 |
| Lipid 5 | 1.1390 | 1.0153 | 1.0153 |
| **mean (the applied offset)** | **1.1116** | **1.0070** | **1.0246** |
| **sd (population, as `build_rows.py` reports it)** | **0.0234** | **0.0247** | **0.0376** |
| sd (sample) | 0.0271 | 0.0286 | 0.0435 |

**The offset is 1.007 ± 0.025, not 1.073 ± 0.034.** The claimed 1.073 is in the right
direction but ~6 % too high. With the red-line-contaminated pair pinned to the stored
`y_top = 2877`, the offset is 1.025 ± 0.038. Either way the honest reading is: **after the
marker contamination is removed there is no measurable offset between the two screens** —
the four benchmarks agree between Fig 2c and Fig 3c to within ±2.5 %, and the offset is not
distinguishable from 1.0. The stored 1.1116 was almost entirely an artifact: ALC-0315 is
dot-contaminated in *both* panels (2c +24 px, 3c +16 px, pulling in opposite directions on
the ratio), Lipid 5 in 3c (+8 px), and DLin/SM102 in 2c sit under the red reference line.

Because `build_rows.py` divides library 2 by this offset, changing 1.1116 → 1.0070 raises
every library-2 value by 10.4 %.

### 4. The two byte-identical pairs — confirmed, and explained

| pair | stored y_top | stored RLU | my y_top | my RLU | my ratio |
|---|---|---|---|---|---|
| T3A-L1-2C10 | 2463 | 933,681.3 | **2465** | 905,636 | — |
| T3A-L7-2C10 | 2463 | 933,681.3 | **2473** | 826,454 | 1.096 |
| T3A-L5-2C10 | 2439 | 1,346,291.4 | **2441** | 1,305,852 | — |
| T3A-L6-2C12 | 2439 | 1,346,291.4 | **2449** | 1,155,883 | 1.130 |

**Confirmed: both pairs are byte-identical in `05_fig3c_24.csv`, and both are spurious.**
In each pair one member is read correctly (L1-2C10 and L5-2C10, dY = +2, i.e. zero to one
native pixel) and the other is dot-contaminated (L7-2C10 read 10 px high, L6-2C12 read 10 px
high). The contaminated read happened to land on the same pixel row as the clean member's
true bar top. Mechanically: the three replicate dots sit roughly one marker-radius above
the rule, so contaminated reads cluster at a near-constant offset, which makes such
collisions more likely than chance — but the exact equality is still a pixel coincidence,
not a detector latching onto a fixed artifact.

**Caveat worth recording: exact ties are not by themselves evidence of a bug.** At 151 px
per decade on an image whose effective row pitch is 2 px, Fig 3c offers only ~75
distinguishable values per decade over a ~3.3-decade bar range, so ties among 28 bars are
unsurprising. My own corrected set still contains one exact tie (T3A-L3-2C10 and
T3A-L1-2C11 both at y = 2557 → 222,680 RLU). The two stored pairs are spurious because I
can show the contaminated member's true top is 10 px lower, not because they were tied.

### 5. Which rows change, and by how much in z

Recomputed end-to-end through `build_rows.py`'s logic (library-1 value = Fig 2c where
plotted and not in conflict, else the Fig 2b heatmap cell; library 2 = Fig 3c ÷ offset;
log10; single pooled z over all 344).

| scenario | offset applied | rms Δz | max \|Δz\| | rows with \|Δz\| > 0.05 |
|---|---|---|---|---|
| stored | 1.1116 | — | — | — |
| F-017 only | 1.1068 | 0.0006 | 0.0033 | 0 |
| F-018 only | 1.0113 | 0.0195 | 0.1515 | 13 |
| **F-017 + F-018** | **1.0070** | **0.0195** | **0.1506** | **13** |

**31 of 344 rows move by more than 0.02 z.** The complete list, F-017 + F-018 applied:

| lipid | lib | stored RLU | corrected RLU | ratio | z old | z new | Δz |
|---|---|---|---|---|---|---|---|
| T3A-L3-C6 | 1 | 122,589 | 97,079 | 0.792 | 2.039 | 1.889 | −0.151 |
| 110-L7-C10 | 1 | 36,119 | 28,928 | 0.801 | 1.247 | 1.104 | −0.143 |
| T3A-L5-C10 | 1 | 44,782 | 37,954 | 0.848 | 1.387 | 1.280 | −0.107 |
| T3A-L3-2C10 | 2 | 284,212 | 222,681 | 0.784 | 2.516 | 2.423 | −0.093 |
| T3A-L2-2C10 | 2 | 492,100 | 397,501 | 0.808 | 2.872 | 2.798 | −0.074 |
| T3A-L5-2C11 | 2 | 731,541 | 731,541 | 1.000 | 3.129 | 3.194 | +0.065 |
| T3A-L6-2C11 | 2 | 435,585 | 435,585 | 1.000 | 2.793 | 2.858 | +0.065 |
| T3A-L7-2C12 | 2 | 236,686 | 236,686 | 1.000 | 2.397 | 2.462 | +0.065 |
| T3A-L7-2C11 | 2 | 174,470 | 174,470 | 1.000 | 2.200 | 2.264 | +0.065 |
| T3A-L6-2C10 | 2 | 136,698 | 136,698 | 1.000 | 2.041 | 2.106 | +0.065 |
| T3A-L3-2C12 | 2 | 51,513 | 51,513 | 1.000 | 1.409 | 1.473 | +0.065 |
| T3A-L4-C10 | 1 | 109,474 | 99,301 | 0.907 | 1.966 | 1.903 | −0.063 |
| T3A-L4-2C10 | 2 | 149,795 | 124,746 | 0.833 | 2.101 | 2.047 | −0.054 |
| 110-L8-C10 | 1 | 58,755 | 54,514 | 0.928 | 1.563 | 1.515 | −0.048 |
| T3A-L5-2C10 | 2 | 1,346,291 | 1,305,852 | 0.970 | 3.524 | 3.570 | +0.045 |
| T3A-L1-2C10 | 2 | 933,681 | 905,636 | 0.970 | 3.287 | 3.332 | +0.045 |
| T3A-L8-2C11 | 2 | 590,913 | 573,163 | 0.970 | 2.990 | 3.036 | +0.045 |
| T3A-L8-2C12 | 2 | 573,163 | 555,947 | 0.970 | 2.971 | 3.016 | +0.045 |
| T3A-L8-2C10 | 2 | 293,013 | 284,212 | 0.970 | 2.536 | 2.581 | +0.045 |
| T3A-L5-2C12 | 2 | 191,186 | 185,444 | 0.970 | 2.259 | 2.304 | +0.045 |
| T3A-L2-2C11 | 2 | 33,612 | 32,602 | 0.970 | 1.132 | 1.177 | +0.045 |
| T3A-L2-2C12 | 2 | 13,463 | 13,059 | 0.970 | 0.539 | 0.584 | +0.045 |
| T3A-L1-2C12 | 2 | 839 | 814 | 0.970 | −1.260 | −1.216 | +0.044 |
| T3A-L5-C14 | 1 | 24,308 | 22,810 | 0.938 | 0.990 | 0.950 | −0.041 |
| T3A-L4-C12 | 1 | 170,200 | 159,715 | 0.938 | 2.252 | 2.212 | −0.040 |
| T3A-L6-2C12 | 2 | 1,346,291 | 1,155,883 | 0.859 | 3.524 | 3.490 | −0.034 |
| T3A-L5-C8 | 1 | 24,864 | 23,597 | 0.949 | 1.005 | 0.972 | −0.033 |
| T3A-L7-C14 | 1 | 82,502 | 78,299 | 0.949 | 1.783 | 1.749 | −0.033 |
| T3A-L1-2C11 | 2 | 236,686 | 222,681 | 0.941 | 2.397 | 2.423 | +0.025 |
| T3A-L4-2C12 | 2 | 12,286 | 11,559 | 0.941 | 0.479 | 0.504 | +0.025 |
| T3A-L4-2C11 | 2 | 5,560 | 5,230 | 0.941 | −0.035 | −0.010 | +0.025 |

Note the sign pattern in library 2: bars whose raw value is unchanged (ratio 1.000) still
gain +0.065 z, purely because the offset divisor drops from 1.1116 to 1.0070. That is the
dominant library-2 effect, not the per-bar correction.

Per-library mean z: library 1 −0.1591 → −0.1610; library 2 +2.1213 → +2.1466; the gap
widens slightly, 2.280 → 2.308. Geometric-mean RLU: lib 1 4,128 → 4,113; lib 2 154,630 →
145,480; the headline 37.5× six-tail-over-three-tail gap becomes 35.4×. **The paper's
between-library conclusion is untouched.**

## WHAT CHANGES

Three edits, all in the two digitizers:

- `curation/ql_2026/digitize_fig2c.py:75–77` and `curation/ql_2026/digitize_fig3c.py:110–117`
  — replace the `w // 8` inner-trim plus `cov >= 0.85` rule with the joined-flank test:
  locate the two side-stroke centres from an 18-px row band just above the baseline, then
  take the first row carrying one contiguous run of ink (gray ≤ 172) that spans ≥72 % of the
  bar width **and contains both stroke centres**, sustained ≥2 rows. Do not substitute a
  stricter coverage threshold alone — a three-dot cluster reaches 91 % contiguous coverage
  on `T3A-L2-2C10` and 97 % on `110-L1-C10`, so coverage cannot separate them; the
  both-flanks-joined condition can.
- `curation/ql_2026/digitize_fig2c.py:35` — the F-017 constant above.
- No change to `build_rows.py`; it recomputes the offset from the CSVs, so it will pick up
  1.0070 automatically. Consider recording in `build_rows.py`'s docstring (line 13–17) that
  the offset is **not** significantly different from 1, since the prose currently asserts
  "It is a SCALAR ~1.1x".

Per-row deltas: the 31 rows above. Flags, censoring and conflict status are unchanged for
every row (the T3A-L8-C14 Fig 2b/2c conflict still triggers: corrected Fig 2c value 57,684
vs heatmap 7,931, ratio 7.3).

## CONFIDENCE: **certain** that the stored detector mis-reads 15 of 60 bars and that the
offset is ≈1.01 rather than 1.11. **Probable** on the individual corrected pixel rows at
the ±2 px level — the figure is a 300-dpi JPEG and the top rules carry a 2–4 px
anti-aliased ramp; the two benchmarks under the red line in Fig 2c are ±4 px.
**Certain** that the two byte-identical pairs are spurious; **speculative** on the exact
claim that two equal values in one panel is "the signature of a detector latching onto
something that is not the rule" — ties are also expected from the pixel grid alone, and my
corrected set contains one.

---

# F-019 — "the Fig 2b colormap floor sits ABOVE the feasible range"

## VERDICT: **confirmed**, and the underlying cause is a second instance of error class 20 —
the colorbar's zero is anchored on the **outer edge** of its black border instead of the
border's centre.

## MY MEASUREMENT

**Page 4, 600 dpi, clipped-figure frame.**

### 1. The colorbar's own calibration, re-derived from its ticks

The colorbar box runs x 3349–3444 (black borders) with the painted ramp in x 3359–3432. It
carries two tick notches, at 4×10⁴ and 2×10⁴, protruding into the bar from both sides.
Darkness-weighted centroids, averaged over the left (x 3361–3375) and right (x 3417–3431)
notches:

- 4×10⁴ at **y = 1393.98**
- 2×10⁴ at **y = 1680.60**
- → **143.31 px per 1×10⁴**, → **0 RLU at y = 1967.23**, 6×10⁴ at y = 1107.36

Independent check from the border geometry: the bottom black border spans y 1962–1974
(centre **1968.0**) and the top border y 1100–1109 (centre **1104.5**), giving 143.9 px per
1×10⁴ and 0 RLU at y = 1968.0. The two estimates agree to 0.8 px. I adopt
**0 RLU at y = 1967.6, 143.6 px per 1×10⁴**.

`digitize_fig2b.py:40` uses `BAR_TOP, BAR_BOT = 1103, 1973`. **1973 is the outer edge of the
bottom border, 5.4 px below the true zero.** The slope is nearly right (145.0 vs 143.6 px
per 1×10⁴, 1 % high) but the zero point is displaced, so **every Fig 2b value is biased high
by a roughly constant ≈ +370 RLU**: +1 % at the 2×10⁴ hit threshold, −0.02 % at 4×10⁴,
and **+38 % at the floor**.

Overlay `…/scratchpad/overlays/F019_colorbar_calibration.png` draws my 0/1e4/…/6e4 rules on
the colorbar: they land on the two notches, on the centre of the bottom border beside the
printed "0", and on the centre of the top border beside "6×10⁴". The script's `BAR_BOT`
line is drawn in blue, visibly below the zero.

### 2. The LUT, and where the ramp dies

LUT = per-row median RGB of the ramp interior (x 3365–3428), tick rows excluded. The bottom
of the ramp, with values on the tick-anchored scale:

| y range | RGB | value range (RLU) | RGB distance to the entry above |
|---|---|---|---|
| 1896–1905 | (252, 249, 244) | 4,986 … 4,359 | — |
| 1906–1921 | (251, 250, 245) | 4,290 … 3,245 | 1.73 |
| 1922–1927 | (252, 251, 247) | 3,176 … 2,827 | 2.45 |
| 1928–1937 | (253, 252, 248) | 2,758 … 2,131 | 1.73 |
| 1938–1947 | (253, 253, 251) | 2,061 … 1,435 | 3.16 |
| 1948–1953 | (254, 254, 252) | 1,365 … 1,017 | 1.73 |
| **1954–1955** | **(255, 255, 255)** | **947 … 877** | **3.32** |
| 1956–1961 | (253,253,253)/(254,254,254)/(245,245,245) | — | border anti-aliasing (neutral grey, no cream tint — not ramp) |

**The ramp reaches pure white (255,255,255) at y = 1954, i.e. at ≈ 947 RLU, and there is no
gradient at all below it.** The claim's "pure white at about 900 RLU" is right.

**Resolution limit, quantified.** Adjacent LUT entries near the floor differ by an RGB
Euclidean distance of **1.73 to 3.46** (1–3 levels out of 255) while each entry spans
**349 to 1,045 RLU**. For comparison, the measurement noise is far smaller: within-cell
pixel sd over the sampled interior of the 113 floor cells is **mean 0.06, max 0.34**
levels, and the cell-to-cell sd of the sampled median among nominally identical cells is
**0.00–0.50** levels. So the limit is the **colormap's own quantisation**, not JPEG noise
or sampling — the figure simply does not encode anything below ~950 RLU.

### 3. Counts and names

From `01_fig2b_320.csv` (320 cells; `floored = 1` on 113):

- **48 cells sit at the colormap white point** (stored value 1,310.3 RLU). Their sampled
  colours are (254,254,254) on 43 and (255,255,255) on 5 — all neutral grey/white with no
  cream tint, i.e. whiter than the last ramp entry (254,254,252), hence at or below the
  white point:

  `110-L5-C6, 113-L2-C6, 113-L2-C8, 113-L3-C8, 113-L3-C10, 113-L5-C8, 113-L6-C10, 113-L8-C6,
  113-L8-C8, 114-L2-C6, 114-L3-C8, 114-L3-C10, 114-L6-C6, 114-L8-C6, 400-L2-C6, 400-L2-C8,
  400-L3-C8, 400-L3-C10, 400-L3-C12, 400-L5-C8, 400-L7-C6, 400-L7-C8, 400-L8-C6, 70-L2-C6,
  70-L2-C8, 70-L2-C14, 70-L3-C8, 70-L5-C6, 70-L5-C8, 70-L5-C10, 70-L8-C6, 70-L8-C8,
  75-L2-C6, 75-L2-C8, 75-L3-C8, 75-L5-C6, 75-L5-C8, 75-L5-C10, 75-L6-C10, 80-L2-C6,
  80-L2-C10, 80-L3-C8, 80-L3-C10, 80-L3-C12, 80-L5-C6, 80-L5-C8, 80-L7-C10, 80-L8-C6`

- **65 cells are exactly one LUT step above**, at (254,254,252), stored 1,724.1 RLU:

  `110-L3-C14, 110-L4-C14, 110-L6-C6, 113-L3-C12, 113-L4-C10, 113-L5-C6, 113-L5-C10,
  113-L6-C6, 113-L7-C6, 114-L2-C8, 114-L2-C14, 114-L3-C12, 114-L4-C10, 114-L5-C6, 114-L6-C8,
  400-L1-C12, 400-L2-C10, 400-L4-C12, 400-L5-C6, 400-L5-C10, 400-L6-C6, 400-L6-C8,
  400-L7-C12, 400-L8-C8, 70-L1-C12, 70-L2-C10, 70-L2-C12, 70-L3-C10, 70-L3-C14, 70-L4-C6,
  70-L5-C12, 70-L6-C6, 70-L6-C10, 70-L6-C12, 70-L6-C14, 70-L7-C8, 70-L7-C12, 70-L8-C14,
  75-L1-C12, 75-L2-C10, 75-L2-C12, 75-L2-C14, 75-L3-C10, 75-L3-C12, 75-L3-C14, 75-L4-C6,
  75-L4-C14, 75-L6-C6, 75-L6-C8, 75-L6-C14, 75-L7-C6, 75-L7-C12, 75-L8-C6, 75-L8-C8,
  80-L1-C12, 80-L4-C14, 80-L5-C10, 80-L6-C6, 80-L6-C8, 80-L6-C10, 80-L6-C12, 80-L8-C10,
  T3A-L2-C6, T3A-L3-C12, T3A-L3-C14`

  48 + 65 = 113, matching the stored `floored` count exactly.

  Honest qualification: 43 of the 48 render (254,254,254), which is a neutral grey between
  the white entry (255,255,255, distance 1.73) and the last cream entry (254,254,252,
  distance 2.00). Those two are separated by **3.32 RGB units**. The separation is real but
  only just above the per-cell reproducibility, so the white/one-step-above split should be
  treated as ±a few cells, not exact.

  Overlay `…/scratchpad/overlays/F019_fig2b_floored_cells.png` boxes the 48 in red, the 65
  in orange and the 21 saturated cells in blue. The red and orange cells are visually
  indistinguishable from each other and from plain white — which is the finding.

### 4. The narrowest defensible interval, and where 1,310.3 sits

The boundary between two LUT entries is the midpoint of their value ranges.

- **White-point cells (n = 48):** anything from 0 up to the midpoint between the white band's
  top (947.1) and the next entry's bottom (1,016.7), i.e. **[0, 982] RLU**. There is no
  information below 947 at all, so the honest answer to "what IS recoverable" is:
  **only the upper bound. The interval is 0 to ~980 RLU, and nothing narrower is defensible.**
- **One-step-above cells (n = 65):** **[982, 1,400] RLU** (midpoint of 1,016.7–1,364.9 and
  its neighbours), i.e. a genuine measurement to ±18 %.

**Arithmetic for the stored values:**

- Stored 1,310.3 vs the white-cell interval [0, 982]:
  1,310.3 − 982 = **+328 RLU above the top of the entire feasible band**; 1,310.3 / 982 =
  **1.33×** the upper bound, and 1,310.3 / 491 (the interval midpoint) = **2.67×**.
  In log10, log10(1310.3) − log10(982) = **0.1254 decades**, which at the slice's pooled
  sd of 0.6699 decades is **0.187 z** — simply replacing 1,310.3 by the defensible upper
  bound already moves those 48 rows by nearly a fifth of a standard deviation, in the same
  direction for all of them.
  **So yes: the stored value is not merely imprecise, it is outside and above the band,
  and the bias is systematic.**
- Stored 1,724.1 vs the one-step-above interval [982, 1,400]: 1,724.1 / 1,400 = **1.23×**,
  also above the band. Same cause (the −5.4 px zero error), smaller in relative terms.
- The next group up is affected too: stored 2,413.8 for the (253,253,251) cells (n = 57) vs
  a true band of [1,400, 2,096] — 1.15× above.

### 5. Sensitivity of the slice's z-scores

All 113 floored cells set to the scenario value, everything else held at stored, pipeline
otherwise identical to `build_rows.py` (offset 1.1116, pooled z over 344):

| scenario | floor value | rms Δz | max \|Δz\| | mean z lib 1 | mean z lib 2 | z gap | pooled sd (decades) |
|---|---|---|---|---|---|---|---|
| (a) stored | 1,310.3 / 1,724.1 | — | — | −0.1591 | +2.1213 | 2.280 | 0.6699 |
| (b) white point | 923 | **0.1277** | **0.2098** | −0.1530 | +2.0398 | 2.193 | 0.7322 |
| (c) half the white point | 462 | **0.2508** | **0.4790** | −0.1436 | +1.9150 | 2.059 | 0.8315 |
| (d) 100 RLU | 100 | **0.4384** | **0.9877** | −0.1252 | +1.6699 | 1.795 | 1.0843 |
| (b*) white point, 48 white cells only | 923 | 0.0711 | 0.1621 | −0.1566 | +2.0879 | 2.244 | 0.6908 |
| (d*) 100 RLU, 48 white cells only | 100 | 0.4156 | 0.8623 | −0.1333 | +1.7773 | 1.911 | 0.8873 |
| (e) LUT-correct (923 / 1,167) | per group | 0.0887 | 0.1541 | −0.1546 | +2.0617 | 2.216 | 0.7151 |

**The claimed sensitivity reproduces**: max |Δz| = **0.21 at the white point** and **0.99 at
100 RLU**, exactly the reported figures; my rms at 100 RLU is 0.438 against the claimed
0.415. This is an order of magnitude larger than F-017 (rms 0.0006) and F-018
(rms 0.0195, max 0.151) combined, and **F-019 is therefore the dominant digitization
uncertainty in this slice**, as claimed.

**Does the paper's between-library conclusion survive?** Yes, in every scenario — and it
strengthens:

| scenario | geo-mean RLU lib 1 | geo-mean RLU lib 2 | ratio | P(random lib-2 > random lib-1) |
|---|---|---|---|---|
| (a) stored | 4,128 | 154,630 | 37.5× | 0.9223 |
| (b) 923 | 3,450 | 154,630 | 44.8× | 0.9223 |
| (c) 462 | 2,701 | 154,630 | 57.3× | 0.9370 |
| (d) 100 | 1,574 | 154,630 | 98.3× | 0.9370 |

The *z gap* shrinks from 2.28 to 1.80 only because pushing 113 library-1 cells downward
inflates the pooled sd from 0.67 to 1.08 decades. On the raw scale the six-tail advantage
grows from 37× to 98×. No ordering is inverted. So the conclusion survives — but the
z-scores that the corpus actually stores are not a stable representation of the effect
size, because 35 % of library 1 sits on an unmeasurable floor whose assumed depth sets the
scale.

### 6. Tie structure of `Experiment_value` — confirmed

`11_QL_2026_rows.csv`: 344 rows, **83 distinct `Experiment_value`s**. Full structure:

| tie size | number of distinct values | source |
|---|---|---|
| 65 | 1 (z = −0.72506) | Fig 2b LUT entry (254,254,252), stored 1,724.1 RLU |
| 57 | 1 (z = −0.50691) | Fig 2b LUT entry (253,253,251), stored 2,413.8 RLU |
| 48 | 1 (z = −0.90298) | Fig 2b colormap **white point**, stored 1,310.3 RLU |
| 22 | 1 (z = −0.34399) | stored 3,103.4 RLU |
| 16 | 2 (z = −0.26285, −0.08594) | stored 3,517.2 and 4,620.7 RLU |
| 10 | 1 | |
| 8 | 1 | |
| 5 | 1 | |
| 4 | 3 | |
| 3 | 4 | |
| 2 | 6 | |
| 1 | 61 | |

**Confirmed exactly: 83 distinct values across 344 rows, three largest ties 65, 57, 48.**
Worth stating plainly: **192 of the 344 rows (56 %) carry one of only six values**, and all
six come from the bottom six LUT entries of the Fig 2b colorbar. The slice's apparent
resolution is an artifact of a colormap that quantises its bottom 4,000 RLU into six steps.

## WHAT CHANGES

- `curation/ql_2026/digitize_fig2b.py:40` —
  `BAR_TOP, BAR_BOT = 1103, 1973` → anchor on the colorbar **ticks**, not the border edges:
  `BAR_TOP, BAR_BOT = 1107.4, 1967.6   # 6e4 and 0 from the 4e4/2e4 tick notches at y=1393.98/1680.60, 143.31 px per 1e4`
  (or equivalently keep two anchors and set `PX_PER_1E4 = 143.6`). This alone removes a
  systematic +370 RLU on all 320 cells and brings 1,310.3 → 947 and 1,724.1 → 1,165.
- `curation/ql_2026/digitize_fig2b.py:86` — the LUT window `np.arange(BAR_TOP+18, BAR_BOT-17)`
  must stop where the ramp stops (y = 1955), not 6 px into the border's grey skirt; rows
  1956–1961 are border anti-aliasing and currently sit in the LUT as (253,253,253) /
  (254,254,254) / (245,245,245), which is why a pure-white cell nearest-matches at y = 1954
  rather than being flagged as off-scale.
- **The real change is semantic, not numeric.** The 48 white-point cells are **upper bounds**,
  not measurements. They should carry `value_is_bound = upper` with the bound at **~980 RLU**
  (not a point value of 1,310.3 and not 947), exactly as the saturated cells carry a lower
  bound. `build_rows.py:36–39` already describes `value_floored` as "an upper bound within
  its band" — but `row_flags` then stores 1,310.3 as the value, which is **above** the band,
  so the stored number contradicts its own flag. Either store the bound (982) or, better,
  store the band midpoint (491) with an explicit interval; the two differ by 0.3 z for those
  rows, which is the honest width of the ignorance.
- Per-row deltas: with the colorbar simply re-anchored and the floor cells set to the LUT
  value (scenario e), **113 rows move**, rms Δz 0.089, max |Δz| 0.154. If the 48 are instead
  carried as upper bounds at 982 RLU, those 48 move by −0.187 z each.

## CONFIDENCE: **certain** that the ramp reaches pure white at ≈947 RLU, that 113 cells sit
in the colormap's dead or near-dead zone, that 1,310.3 is above the feasible band for 48 of
them, and that the tie structure is 83/65/57/48. **Certain** on the colorbar zero-point
error (two independent anchors agree to 0.8 px and both disagree with 1973 by >5 px).
**Probable** on the exact 48/65 split, which depends on a 3.3-unit RGB separation.
**Speculative** on any attempt to put a point value on a white cell — the data are not there.

---

# IF I HAD TO PICK ONE

**F-019.** It is the only one of the three that changes what the dataset *means*.

The arithmetic decides it. F-017 moves the slice by rms 0.0006 z — three orders of
magnitude below the others; it is a documentation defect (the comment on
`digitize_fig2c.py:35` describes anchors that are not the ones the figure has) attached to a
constant that happens to be right, and the prior reviewer's proposed "fix" would have made
it worse. F-018 is a real detector bug that mis-reads a quarter of the bars and inflates the
cross-library offset from a true ~1.01 to 1.1116, but it moves rms 0.0195 z and tops out at
0.151 z on one row, and it touches only 51 of 344 rows because Fig 2c supplies values for
just 27 library-1 lipids. Both are worth fixing; neither changes a conclusion.

F-019 is different in kind. One hundred and thirteen of 320 library-1 cells — 35 % of the
larger library, 33 % of the whole slice — carry values the figure does not contain. Forty-
eight of them are not measurements at all: the colormap is flat white from 0 to ~947 RLU, so
those cells are upper bounds, and the stored 1,310.3 sits 1.33× above the top of the entire
band they could occupy and 2.67× above its midpoint. That is a one-directional bias on a
sixth of the slice. Depending only on what one assumes those cells really are, the slice's
z-scores move by up to 0.99 — larger than the combined effect of every other digitization
error in this curation by a factor of six, and larger than the z separation between many
adjacent lipids. Any model trained on `Experiment_value` is being told that 192 of 344
LNPs take one of six values, all six set by the bottom six quantisation steps of a
colorbar, and the lowest of them placed above its own feasible range.

Two caveats on scope. First, the *paper's* headline — six-tail lipids beat three-tail ones —
survives all three findings and in fact strengthens under every F-019 scenario (37× → 98×
on raw geometric means, with no rank inversions). F-019 is a defect in the corpus slice, not
a challenge to Liu et al. Second, and this is why I have not leaned on any printed statistic
anywhere above: the Fig 3d/3e hit rates and the "17 lipids > 1×10⁵" sentence are computed
from the *top* of the dynamic range and are blind to everything that happens at the floor,
so they cannot discriminate F-019 at all, and they also survived F-018's 25 % errors. Every
verdict here rests on re-derived geometry with an overlay, plus two external round-number
checks that the paper did not print as statistics: Fig 2c's red reference line landing on
5.00×10⁴ (0.04 %) and the colorbar's zero landing on the centre of its own border.

---

# Addendum, 2026-10-30 — F-018 patch, bar-by-bar verification, and knock-on numbers after F-017 + F-019

Written at the coordinator's request after F-017 and F-019 were applied. Still read-only on
the repository except this file: the patched scripts were built and run in
`…/scratchpad/patchtest/`, which contains copies of `digitize_fig2c.py`, `digitize_fig3c.py`,
`00_pages/` and `01_fig2b_320.csv` taken from the repo at 18:34, so the diffs below apply
cleanly to the current files.

## A1. The geometry question — every bar is hollow; nothing is solid

The reverted attempt concluded "Fig 2c's first bar is 12 px wide and the stroke band shows
all 12 columns fully coloured, i.e. solid, not hollow". That reading was one level of
nesting off. **The 12 columns are the left side stroke, not the bar.** Measured directly, in
the band y = [baseline−26, baseline−8):

| | Fig 2c | Fig 3c |
|---|---|---|
| bar width | 66–72 px | 46–50 px |
| coloured column groups per bar | **exactly 2** on all 32 | **exactly 2** on all 28 |
| stroke widths | 8–14 px | 4–10 px |
| white interior fraction of the band | 0.67–0.76 | 0.57–0.70 |
| solid bars | **0** | **0** |

For `110-L1-C10` specifically: bar x = 557–627 (70 px), left stroke cols 557–568 (12 px),
right stroke cols 619–626 (8 px), 48 white columns between them, interior fraction 0.69.

I also ran the reverted `_bar_top_UNUSED` as it stands in the repo. **It does not fail on
bar 1** — it returns 2977, which is within 2 px of my hand measurement of 2979. Across all
60 bars it agrees with my independent measurement to ≤2 px on **56**, and differs on four:

| bar | `_bar_top_UNUSED` | mine | cause |
|---|---|---|---|
| DLin-MC3-DMA (2c) | `None` | 2869 | red 5×10⁴ dotted line crosses the top rule; a **blue-only mask** is punched through there |
| SM102 (2c) | 2877 | 2873 | same |
| T3A-L2-2C10 (3c) | 2511 | 2519 | the `i - cur[-1] > 2` **gap bridging** welds the merged dot cluster into one run that reaches both stroke centres |
| ALC-0315 (3c) | 2573 | 2581 | same |

So the reverted code was close, and the two things it needed were a **luminance mask instead
of a blue mask** and **strictly contiguous runs instead of 2-px gap bridging** — not a
different model of the geometry. Reverting was still the right call, since neither of us had
validated it bar by bar; what follows is that validation.

## A2. The patch

Both files get the same two helpers (`_runs`, `bar_top`) and the same call-site change, so a
12-px stroke on a wide hollow bar and a hypothetical solid bar go through one code path: the
`len(groups) >= 2` branch takes the outermost two stroke centres, the `len(groups) == 1`
branch takes the quarter points of a solid fill, and the `else` falls back to just inside the
bar edges. On this figure only the first branch ever fires, and the function says so.

Full unified diffs are saved as
`…/scratchpad/overlays/digitize_fig2c.py.patch` (234 lines) and
`…/scratchpad/overlays/digitize_fig3c.py.patch` (163 lines). The substance:

### Shared helpers (inserted in both files, replacing `_bar_top_UNUSED` in `digitize_fig2c.py`)

```python
INK_THR = 172            # "anything drawn": see bar_top() below


def _runs(idx, gap=0):
    """Split a sorted index list into runs, allowing `gap` missing entries.

    gap=0 means STRICTLY contiguous. That matters: bridging even a 2 px gap is
    enough to let a merged replicate-dot cluster pass the test below (it does
    on T3A-L2-2C10 and on ALC-0315 in Fig 3c).
    """
    out, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] > gap + 1:
            out.append(cur)
            cur = []
        cur.append(i)
    out.append(cur)
    return out


def bar_top(ink, x0, x1, top_search, base_y, frac=0.72, nrow=2,
            band=(24, 8), min_stroke=3):
    """Topmost row of a bar's top rule, immune to the replicate dot markers.

    THE GEOMETRY, measured rather than assumed. All 60 bars in Figs 2c and 3c
    are HOLLOW: in a band just above the baseline each one shows exactly two
    fully-coloured column groups (the side strokes, 4-14 px wide) with 57-76%
    white between them. Fig 2c bars are 66-72 px wide, Fig 3c bars 46-50 px.
    None is solid. An earlier attempt at this fix read the 12-column LEFT
    STROKE of the first bar (110-L1-C10, cols 557-568) as if it were the bar
    and concluded the bars were 12 px wide and solid; the bar is cols 557-627.
    The branch below still handles a solid bar, so the two cases share one code
    path, but on this figure that branch never fires.

    THE RULE. Take the first row carrying ONE STRICTLY CONTIGUOUS run of ink
    inside [x0, x1) that spans at least `frac` of the bar width AND contains
    both side strokes' centres, sustained for `nrow` rows (2 = one native pixel
    of the 300 dpi source). A dot cluster cannot satisfy it: it does not reach
    both strokes without a gap. The bar's own top rule does, because it is the
    thing that joins them.

    `ink` must be a LUMINANCE mask ("anything drawn"), not a blue mask. Fig 2c's
    red dotted 5e4 reference line crosses the top rule of DLin-MC3-DMA and
    SM102; on a blue-only mask the rule is punched through there and
    DLin-MC3-DMA returns no top at all. The consequence is that `top_search`
    must now clear any BLACK furniture too -- see the assert in each caller.

    Returns (y_top, y_bot, frac_used), or (None, None, None).
    """
    w = x1 - x0
    # 1. side strokes, from a band just above the baseline where the hollow
    #    interior is unambiguously white.
    colfrac = ink[base_y - band[0]:base_y - band[1], x0:x1].mean(0)
    on = [i for i, v in enumerate(colfrac) if v > 0.95]
    groups = [g for g in _runs(on) if len(g) >= min_stroke] if on else []
    if len(groups) >= 2:               # hollow bar: outermost two groups
        cL = (groups[0][0] + groups[0][-1]) / 2.0
        cR = (groups[-1][0] + groups[-1][-1]) / 2.0
    elif len(groups) == 1:             # solid bar: quarter points of the fill
        a, b = groups[0][0], groups[0][-1]
        cL, cR = a + (b - a) * 0.25, a + (b - a) * 0.75
    else:                              # nothing resolvable: just inside the edges
        cL, cR = 1.0, w - 2.0

    # 2. scan down for the joined top rule. The fallback fracs never fire on
    #    this figure; `frac_used` is written to the CSV so that if one ever
    #    does, the bar is visible as needing review rather than silently rescued.
    for f in (frac, 0.66, 0.60):
        need = f * w
        ok = np.zeros(base_y - 6 - top_search, dtype=bool)
        for k in range(len(ok)):
            idx = np.flatnonzero(ink[top_search + k, x0:x1])
            if idx.size == 0:
                continue
            for sg in _runs(list(idx)):
                if (sg[-1] - sg[0] + 1) >= need and sg[0] <= cL and sg[-1] >= cR:
                    ok[k] = True
                    break
        for k in range(len(ok) - nrow + 1):
            if ok[k:k + nrow].all():
                j = k
                while j < len(ok) and ok[j]:
                    j += 1
                return top_search + k, top_search + j - 1, f
    return None, None, None
```

### `digitize_fig2c.py` — mask construction (after `coloured = …`)

```python
    # Bar RUNS are found on the blue mask (the ink mask would merge the black
    # axis into them); bar TOPS are read on the ink mask (the blue mask breaks
    # where the red 5e4 reference line crosses a top rule).
    ink = img.mean(2) <= INK_THR
```

### `digitize_fig2c.py` — the measuring loop

```python
    # TOP_SEARCH must clear every black element above the bars, because the
    # top rule is now read on a luminance mask. Guard it rather than trust it.
    for (x0, x1) in runs:
        idx = np.flatnonzero(ink[TOP_SEARCH, x0:x1])
        assert idx.size == 0 or max(r[-1] - r[0] + 1 for r in _runs(list(idx))) \
            < 0.72 * (x1 - x0), "TOP_SEARCH starts inside a full-width rule"

    rows, boxes = [], []
    for (x0, x1), lab in zip(runs, LABELS):
        # Corrected 2026-10-30, review finding F-018. The old rule took the
        # first row where >=85% of a w//8-trimmed window was coloured; three
        # overlapping replicate dots fill >=85% of THAT while covering only
        # 55-70% of the full bar, so it fired on the dots on 15 of the 60 bars
        # across Figs 2c and 3c, by up to +28% in value, and inflated the
        # cross-library offset from a true ~1.01 to 1.1116.
        ytop, ybot, fused = bar_top(ink, x0, x1, TOP_SEARCH, int(BASE_Y))
        assert ytop is not None, f"no bar top for {lab}"
        assert ytop < int(BASE_Y) - 26, f"bar top for {lab} overlaps the stroke band"
        value = 10 ** (4.0 + (Y_1E4 - ytop) / PX_PER_DECADE)
        rows.append({"lipid_name": lab, "is_benchmark": int(lab in BENCH),
                     "RLU": round(float(value), 1), "x0": x0, "x1": x1,
                     "y_top": ytop, "rule_thickness_px": ybot - ytop + 1,
                     "frac_used": fused, "source_figure": "Fig 2c",
                     "extraction_method": "bar_chart_log_axis"})
        boxes.append((x0, ytop, x1))
```

Also in `digitize_fig2c.py`, the QC `anchors` string is still the stale pre-F-017 text and
should go with it:

```python
-        "calibration": {"y_1e4": Y_1E4, "px_per_decade": PX_PER_DECADE,
-                        "anchors": "label centroids 1e6=2342, 1e5=2748, 1e4=3156"},
+        "calibration": {"y_1e4": Y_1E4, "px_per_decade": PX_PER_DECADE,
+                        "anchors": "y-axis TICK centroids 2342.78/2749.54/3156.75 (F-017)"},
```

### `digitize_fig3c.py` — a new constant, and why it is not optional

```python
INK_THR = 172            # "anything drawn": see bar_top() below
# The panel BOX TOP rule occupies y = 2301-2314. The old search start,
# int(Y_1E3 - 4*PX_PER_DECADE) + 4 = 2311, is INSIDE it. That was harmless
# while the top was read on a blue mask (the box is black), but bar_top() reads
# a luminance mask, and at 2311 every bar returns 2311. Hence an explicit
# constant, with a guard in main().
TOP_SEARCH = 2320
```

with `top = int(Y_1E3 - 4.0 * PX_PER_DECADE) + 4` replaced by `top = TOP_SEARCH`, the same
`ink = img.mean(2) <= INK_THR`, the same guard and loop as above, and
`y0 = TOP_SEARCH - 20` in `overlay()`.

**This is the one trap in the patch.** I verified it: with `top_search = 2311` the patched
detector returns 2311 for every Fig 3c bar. Switching from a hue mask to a luminance mask
makes black furniture visible to the detector, so the search window has to be tightened and
guarded, not inherited. The guard asserts rather than assumes.

## A3. Robustness of the one free parameter

`INK_THR` is the only tuned number. Sweeping it over the whole plausible range, measuring
every one of the 60 bars at each value:

| `INK_THR` | bars returning `None` | frac fallback used | max shift vs 172 | bars moved |
|---|---|---|---|---|
| 150 | 3 | 0 | 4 px | 9 |
| 160 | 1 | 0 | 4 px | 6 |
| 166 | 0 | 0 | 4 px | 4 |
| **172** | **0** | **0** | — | — |
| 178 | 0 | 0 | 2 px | 5 |
| 185 | 0 | 0 | 8 px | 11 |

172 sits in the middle of a stable plateau. It breaks down below ~160 (bars stop resolving)
and above ~180 (at 185, ALC-0315 in Fig 3c jumps 8 px back onto the dot answer). Within
[166, 178] no bar moves by more than 4 px and only one moves by 4. The `frac` fallback chain
(0.66, 0.60) **never fires at any threshold**; `frac_used` is written to the CSV so that if it
ever does, the bar shows up as needing review instead of being silently rescued.

## A4. Per-bar verification — all 60 bars

Produced by running the patched scripts in place; `repo` columns are the current
`07_fig2c_bars.csv` / `05_fig3c_24.csv` (F-017 applied, old bar rule). "My independent y" is
the hand measurement from the main report, made with different code. Bars whose top moves by
≥10 px are bolded.

| # | panel | lipid | repo `y_top` | patched `y_top` | Δ px | my independent y | agrees ≤2 px | repo RLU | patched RLU | ratio | rule px |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2c | 110-L1-C10 | 2977 | 2979 | +2 | 2979 | yes | 27,647 | 27,336 | 0.989 | 12 |
| 2 | 2c | 110-L2-C8 | 2815 | 2817 | +2 | 2817 | yes | 69,132 | 68,354 | 0.989 | 10 |
| 3 | 2c | 110-L2-C10 | 2739 | 2739 | +0 | 2739 | yes | 106,270 | 106,270 | 1.000 | 12 |
| 4 | 2c | 110-L3-C6 | 2777 | 2777 | +0 | 2777 | yes | 85,713 | 85,713 | 1.000 | 8 |
| 5 | 2c | 110-L4-C8 | 2977 | 2977 | +0 | 2977 | yes | 27,647 | 27,647 | 1.000 | 12 |
| 6 | 2c | **110-L7-C10** | 2929 | 2969 | +40 | 2969 | yes | 36,273 | 28,927 | 0.797 | 10 |
| 7 | 2c | 110-L8-C8 | 2933 | 2935 | +2 | 2935 | yes | 35,461 | 35,062 | 0.989 | 8 |
| 8 | 2c | **110-L8-C10** | 2843 | 2857 | +14 | 2857 | yes | 59,004 | 54,511 | 0.924 | 8 |
| 9 | 2c | T3A-L1-C10 | 2743 | 2743 | +0 | 2743 | yes | 103,893 | 103,893 | 1.000 | 8 |
| 10 | 2c | T3A-L2-C8 | 2899 | 2905 | +6 | 2901 | **NO (+4)** | 42,983 | 41,548 | 0.967 | 10 |
| 11 | 2c | T3A-L2-C10 | 2803 | 2805 | +2 | 2805 | yes | 73,989 | 73,156 | 0.989 | 6 |
| 12 | 2c | T3A-L2-C12 | 2835 | 2835 | +0 | 2835 | yes | 61,736 | 61,736 | 1.000 | 8 |
| 13 | 2c | **T3A-L3-C6** | 2713 | 2755 | +42 | 2755 | yes | 123,110 | 97,074 | 0.789 | 6 |
| 14 | 2c | T3A-L3-C10 | 2945 | 2947 | +2 | 2947 | yes | 33,134 | 32,761 | 0.989 | 6 |
| 15 | 2c | T3A-L4-C8 | 2627 | 2629 | +2 | 2629 | yes | 200,261 | 198,008 | 0.989 | 10 |
| 16 | 2c | **T3A-L4-C10** | 2733 | 2751 | +18 | 2751 | yes | 109,940 | 99,295 | 0.903 | 10 |
| 17 | 2c | **T3A-L4-C12** | 2655 | 2667 | +12 | 2667 | yes | 170,923 | 159,704 | 0.934 | 10 |
| 18 | 2c | **T3A-L5-C8** | 2995 | 3005 | +10 | 3005 | yes | 24,970 | 23,597 | 0.945 | 10 |
| 19 | 2c | **T3A-L5-C10** | 2891 | 2921 | +30 | 2921 | yes | 44,973 | 37,952 | 0.844 | 10 |
| 20 | 2c | **T3A-L5-C14** | 2999 | 3011 | +12 | 3011 | yes | 24,411 | 22,809 | 0.934 | 12 |
| 21 | 2c | T3A-L6-C14 | 2977 | 2977 | +0 | 2977 | yes | 27,647 | 27,647 | 1.000 | 12 |
| 22 | 2c | T3A-L7-C10 | 2881 | 2881 | +0 | 2881 | yes | 47,590 | 47,590 | 1.000 | 10 |
| 23 | 2c | T3A-L7-C12 | 2511 | 2515 | +4 | 2515 | yes | 386,020 | 377,382 | 0.978 | 14 |
| 24 | 2c | **T3A-L7-C14** | 2783 | 2793 | +10 | 2793 | yes | 82,852 | 78,295 | 0.945 | 12 |
| 25 | 2c | T3A-L8-C8 | 2765 | 2769 | +4 | 2769 | yes | 91,734 | 89,682 | 0.978 | 12 |
| 26 | 2c | T3A-L8-C10 | 2669 | 2671 | +2 | 2671 | yes | 157,908 | 156,131 | 0.989 | 6 |
| 27 | 2c | T3A-L8-C12 | 2699 | 2701 | +2 | 2701 | yes | 133,258 | 131,759 | 0.989 | 8 |
| 28 | 2c | T3A-L8-C14 | 2845 | 2847 | +2 | 2847 | yes | 58,340 | 57,684 | 0.989 | 10 |
| 29 | 2c | DLin-MC3-DMA | 2877 | 2869 | -8 | 2869 | yes | 48,680 | 50,933 | 1.046 | 8 |
| 30 | 2c | SM102 | 2877 | 2873 | -4 | 2873 | yes | 48,680 | 49,794 | 1.023 | 8 |
| 31 | 2c | **ALC-0315** | 2643 | 2667 | +24 | 2667 | yes | 182,930 | 159,704 | 0.873 | 10 |
| 32 | 2c | Lipid 5 | 2803 | 2805 | +2 | 2805 | yes | 73,989 | 73,156 | 0.989 | 10 |
| 33 | 3c | T3A-L1-2C10 | 2463 | 2465 | +2 | 2465 | yes | 933,681 | 905,636 | 0.970 | 6 |
| 34 | 3c | **T3A-L2-2C10** | 2505 | 2523 | +18 | 2519 | **NO (+4)** | 492,100 | 373,979 | 0.760 | 6 |
| 35 | 3c | **T3A-L3-2C10** | 2541 | 2557 | +16 | 2557 | yes | 284,212 | 222,680 | 0.784 | 8 |
| 36 | 3c | **T3A-L4-2C10** | 2583 | 2595 | +12 | 2595 | yes | 149,795 | 124,746 | 0.833 | 6 |
| 37 | 3c | T3A-L5-2C10 | 2439 | 2441 | +2 | 2441 | yes | 1,346,291 | 1,305,852 | 0.970 | 6 |
| 38 | 3c | T3A-L6-2C10 | 2589 | 2589 | +0 | 2589 | yes | 136,698 | 136,698 | 1.000 | 4 |
| 39 | 3c | T3A-L7-2C10 | 2463 | 2471 | +8 | 2471 | yes | 933,681 | 826,454 | 0.885 | 4 |
| 40 | 3c | T3A-L8-2C10 | 2539 | 2541 | +2 | 2541 | yes | 293,013 | 284,212 | 0.970 | 8 |
| 41 | 3c | T3A-L1-2C11 | 2553 | 2557 | +4 | 2557 | yes | 236,686 | 222,680 | 0.941 | 10 |
| 42 | 3c | T3A-L2-2C11 | 2681 | 2683 | +2 | 2683 | yes | 33,612 | 32,602 | 0.970 | 2 |
| 43 | 3c | T3A-L3-2C11 | 2611 | 2617 | +6 | 2617 | yes | 97,739 | 89,193 | 0.913 | 4 |
| 44 | 3c | T3A-L4-2C11 | 2799 | 2803 | +4 | 2803 | yes | 5,560 | 5,230 | 0.941 | 8 |
| 45 | 3c | T3A-L5-2C11 | 2479 | 2479 | +0 | 2479 | yes | 731,541 | 731,541 | 1.000 | 8 |
| 46 | 3c | T3A-L6-2C11 | 2513 | 2515 | +2 | 2513 | yes | 435,585 | 422,501 | 0.970 | 6 |
| 47 | 3c | T3A-L7-2C11 | 2573 | 2573 | +0 | 2573 | yes | 174,470 | 174,470 | 1.000 | 8 |
| 48 | 3c | T3A-L8-2C11 | 2493 | 2495 | +2 | 2495 | yes | 590,912 | 573,163 | 0.970 | 8 |
| 49 | 3c | T3A-L1-2C12 | 2923 | 2925 | +2 | 2925 | yes | 839 | 814 | 0.970 | 6 |
| 50 | 3c | T3A-L2-2C12 | 2741 | 2743 | +2 | 2743 | yes | 13,463 | 13,058 | 0.970 | 4 |
| 51 | 3c | T3A-L3-2C12 | 2653 | 2653 | +0 | 2653 | yes | 51,513 | 51,513 | 1.000 | 6 |
| 52 | 3c | T3A-L4-2C12 | 2747 | 2751 | +4 | 2751 | yes | 12,286 | 11,559 | 0.941 | 8 |
| 53 | 3c | T3A-L5-2C12 | 2567 | 2569 | +2 | 2569 | yes | 191,186 | 185,444 | 0.970 | 8 |
| 54 | 3c | **T3A-L6-2C12** | 2439 | 2449 | +10 | 2449 | yes | 1,346,291 | 1,155,883 | 0.859 | 8 |
| 55 | 3c | T3A-L7-2C12 | 2553 | 2553 | +0 | 2553 | yes | 236,686 | 236,686 | 1.000 | 8 |
| 56 | 3c | T3A-L8-2C12 | 2495 | 2497 | +2 | 2497 | yes | 573,163 | 555,947 | 0.970 | 4 |
| 57 | 3c | DLin-MC3-DMA | 2649 | 2653 | +4 | 2653 | yes | 54,753 | 51,513 | 0.941 | 8 |
| 58 | 3c | SM102 | 2651 | 2653 | +2 | 2653 | yes | 53,109 | 51,513 | 0.970 | 8 |
| 59 | 3c | **ALC-0315** | 2565 | 2581 | +16 | 2581 | yes | 197,107 | 154,434 | 0.784 | 4 |
| 60 | 3c | Lipid 5 | 2621 | 2629 | +8 | 2629 | yes | 83,915 | 74,278 | 0.885 | 8 |

**Every bar is accounted for: 60 rows, 58 agreeing with my independent measurement to ≤2 px
(one native pixel), 2 differing by 4 px.** Both exceptions are the threshold-sensitive bars
identified in A3:

- `T3A-L2-C8` (2c): patched 2905, hand 2901. 4 px = **1.0 %** in value. Immaterial.
- `T3A-L2-2C10` (3c): patched 2523, hand 2519. 4 px = **6.2 %** in value, because Fig 3c has
  only 151 px per decade. The row at y = 2519–2520 sits exactly on the ink threshold (its
  darkest pixels land in the 153–178 grey band that `INK_THR = 172` bisects). The honest
  value is **y = 2521 ± 2**, i.e. 373,979–449,074 RLU. Either way the repo's 2505 is wrong by
  14–18 px; this residual is small against the defect being fixed, but it should be recorded
  rather than presented as exact.

Rule thicknesses come out at 6–14 px in Fig 2c and 2–10 px in Fig 3c, consistent with a
stroked rectangle outline at the two panels' different scales — a sanity check that the
detector is finding rules and not dots, since a dot cluster would not produce a band of that
thickness bounded below by the hollow interior.

Patched-run overlays (the scripts' own output, not hand-drawn):
`…/scratchpad/overlays/F018_PATCHED_fig2c_overlay.png` and
`…/scratchpad/overlays/F018_PATCHED_fig3c_overlay.png`.

## A5. Knock-on numbers, recomputed on top of F-017 + F-019

Re-derived through `build_rows.py`'s value-assembly, offset and pooled-z logic, reading the
repo's **current** `01_fig2b_320.csv` (F-019 applied: floor cells now 948.6 / 1,367.1 /
2,064.6 RLU, 113 floored, 21 saturated) and `09_enumerated_344.csv`.

### Cross-library offset

| benchmark | repo now (F-017+F-019, old bar rule) | with F-018 patch |
|---|---|---|
| DLin-MC3-DMA | 1.1248 | 1.0114 |
| SM102 | 1.0910 | 1.0345 |
| ALC-0315 | 1.0775 | 0.9670 |
| Lipid 5 | 1.1342 | 1.0153 |
| **mean (applied)** | **1.1069** | **1.0071** |
| **sd (population)** | **0.0233** | **0.0247** |
| sd (sample) | 0.0270 | 0.0286 |

This reproduces the coordinator's current `build_rows.py` figure of 1.1069 ± 0.0233 exactly,
and lands within 0.0001 of the 1.0070 I reported before F-019 — as expected, since F-019
touches only the Fig 2b heatmap and the offset is measured on the two bar panels. The
conclusion is unchanged and is worth stating in the prose at `build_rows.py:13-17`, which
currently asserts "It is a SCALAR ~1.1x": **after the marker contamination is removed the two
screens agree to within ±2.5 % and the offset is not distinguishable from 1.0.**

### Slice statistics

| | repo now | with F-018 patch |
|---|---|---|
| pooled log10 mean | 3.6637 | 3.6631 |
| pooled log10 sd | 0.7094 | 0.7084 |
| mean z, library 1 | −0.1566 | −0.1583 |
| mean z, library 2 | +2.0885 | +2.1106 |
| **z gap** | **2.2451** | **2.2689** |
| geo-mean RLU, library 1 | 3,569 | 3,556 |
| geo-mean RLU, library 2 | 154,630 | 144,926 |
| **geo-mean ratio** | **43.3×** | **40.8×** |

Effect of applying F-018 on top: **rms Δz = 0.0188, max |Δz| = 0.1419, 31 rows move by more
than 0.02 z, 12 by more than 0.05 z.** Essentially identical to the 0.0195 / 0.1506 I
projected before F-019 landed. The six-tail advantage stays at ~41× and the between-library
conclusion is untouched.

Worth noting how much F-019 already moved the baseline on its own: the geometric mean of
library 1 fell from 4,128 to 3,569 RLU and the geo-mean ratio rose from 37.5× to 43.3×,
entirely from re-anchoring the colorbar. That remains the largest of the three effects.

### Rows that change when F-018 is applied (|Δz| > 0.02, 31 of 344)

| lipid | lib | repo RLU | patched RLU | ratio | z now | z after | Δz |
|---|---|---|---|---|---|---|---|
| T3A-L3-C6 | 1 | 123,110 | 97,074 | 0.789 | 2.011 | 1.869 | −0.142 |
| 110-L7-C10 | 1 | 36,273 | 28,927 | 0.797 | 1.263 | 1.127 | −0.136 |
| T3A-L2-2C10 | 2 | 492,100 | 373,979 | 0.760 | 2.797 | 2.692 | −0.105 |
| T3A-L5-C10 | 1 | 44,973 | 37,952 | 0.844 | 1.395 | 1.293 | −0.101 |
| T3A-L3-2C10 | 2 | 284,212 | 222,680 | 0.784 | 2.461 | 2.374 | −0.087 |
| T3A-L5-2C11 | 2 | 731,541 | 731,541 | 1.000 | 3.040 | 3.103 | +0.063 |
| T3A-L7-2C12 | 2 | 236,686 | 236,686 | 1.000 | 2.349 | 2.411 | +0.062 |
| T3A-L7-2C11 | 2 | 174,470 | 174,470 | 1.000 | 2.162 | 2.224 | +0.062 |
| T3A-L6-2C10 | 2 | 136,698 | 136,698 | 1.000 | 2.013 | 2.075 | +0.062 |
| T3A-L3-2C12 | 2 | 51,513 | 51,513 | 1.000 | 1.416 | 1.476 | +0.061 |
| T3A-L4-C10 | 1 | 109,940 | 99,295 | 0.903 | 1.942 | 1.883 | −0.059 |
| T3A-L4-2C10 | 2 | 149,795 | 124,746 | 0.833 | 2.069 | 2.019 | −0.050 |
| 110-L8-C10 | 1 | 59,004 | 54,511 | 0.924 | 1.561 | 1.515 | −0.045 |
| T3A-L5-2C10 | 2 | 1,346,291 | 1,305,852 | 0.970 | 3.413 | 3.458 | +0.045 |
| T3A-L1-2C10 | 2 | 933,681 | 905,636 | 0.970 | 3.189 | 3.234 | +0.045 |
| T3A-L8-2C11 | 2 | 590,912 | 573,163 | 0.970 | 2.909 | 2.954 | +0.044 |
| T3A-L8-2C12 | 2 | 573,163 | 555,947 | 0.970 | 2.891 | 2.935 | +0.044 |
| T3A-L6-2C11 | 2 | 435,585 | 422,501 | 0.970 | 2.723 | 2.767 | +0.044 |
| T3A-L8-2C10 | 2 | 293,013 | 284,212 | 0.970 | 2.480 | 2.523 | +0.044 |
| T3A-L5-2C12 | 2 | 191,186 | 185,444 | 0.970 | 2.218 | 2.262 | +0.043 |
| T3A-L2-2C11 | 2 | 33,612 | 32,602 | 0.970 | 1.154 | 1.196 | +0.042 |
| T3A-L2-2C12 | 2 | 13,463 | 13,058 | 0.970 | 0.594 | 0.635 | +0.041 |
| T3A-L5-C14 | 1 | 24,411 | 22,809 | 0.934 | 1.020 | 0.981 | −0.039 |
| T3A-L1-2C12 | 2 | 839 | 814 | 0.970 | −1.105 | −1.067 | +0.039 |
| T3A-L4-C12 | 1 | 170,923 | 159,704 | 0.934 | 2.212 | 2.174 | −0.038 |
| T3A-L5-C8 | 1 | 24,970 | 23,597 | 0.945 | 1.034 | 1.002 | −0.032 |
| T3A-L7-C14 | 1 | 82,852 | 78,295 | 0.945 | 1.769 | 1.737 | −0.031 |
| T3A-L6-2C12 | 2 | 1,346,291 | 1,155,883 | 0.859 | 3.413 | 3.384 | −0.030 |
| T3A-L1-2C11 | 2 | 236,686 | 222,680 | 0.941 | 2.349 | 2.374 | +0.025 |
| T3A-L4-2C12 | 2 | 12,286 | 11,559 | 0.941 | 0.538 | 0.560 | +0.022 |
| T3A-L4-2C11 | 2 | 5,560 | 5,230 | 0.941 | 0.052 | 0.074 | +0.022 |

The library-2 rows with ratio 1.000 still gain ~+0.06 z: their raw readings do not change,
only the offset divisor does. That is the dominant library-2 effect.

## A6. Three corrections to my own earlier report

1. **Wrong.** I wrote that `digitize_fig2b.py:86`'s LUT window runs "6 px into the border's
   grey skirt" and that rows 1956–1961 "currently sit in the LUT". They do not:
   `np.arange(BAR_TOP + 18, BAR_BOT - 17)` = `arange(1121, 1956)` ends at **1955**, which is
   exactly where the painted ramp ends. That recommendation was unnecessary and should be
   ignored. The ramp/LUT boundary analysis it sat next to is unaffected — the white point is
   still y = 1954–1955.
2. **Scope.** My F-018 patch note said "no change to `build_rows.py`". Still true for the
   code, but the docstring claim "It is a SCALAR ~1.1x" becomes false once the patch lands.
3. **Not evidence.** After the patch, Fig 3c's hit check reproduces the printed 3d and 3e
   rates exactly (2C10 3/8, 2C11 1/8, 2C12 1/8; L1 1/3, L5 2/3, L6 1/3, L7 1/3), and the
   borderline case the old `06_fig3c_qc.json` had to explain in prose — a bar sitting within
   one pixel of the 7×10⁵ threshold — is no longer borderline. **I am not offering this as
   evidence that the patch is right**, for the reason recorded in the main report: these same
   rates also reproduced on the uncorrected values, so they cannot discriminate. I mention it
   only because that stale `resolution_note` in `06_fig3c_qc.json` (which still quotes
   701,135 RLU, a value from a superseded calibration) can now be deleted rather than
   rewritten.

## A7. What I did not do

I did not touch `build_rows.py`, `11_QL_2026_rows.csv` or `12_QL_2026_qc.json`; the numbers
in A5 are re-derived, not produced by running the pipeline. After the patch is applied the
pipeline should be re-run and `12_QL_2026_qc.json`'s `cross_figure_offset` block checked
against **1.0071 / sd 0.0233→0.0247**; if it reports anything else, the patch did not land as
tested.
