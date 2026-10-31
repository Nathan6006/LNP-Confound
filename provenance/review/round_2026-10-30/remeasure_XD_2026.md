# Re-measurement of XD_2026 Fig 5C — finding F-020 (bar-top detector fires inside the bar)

Reviewer: independent re-derivation, 2026-10-30. Read-only on the repository; nothing in
`curation/` or `data/` was edited or regenerated. All geometry below was derived from a fresh
600-dpi render of page 13 of `curation/xd_2026_lipidai/1-s2.0-S2211383525007221-main.pdf`
made with PyMuPDF in the scratchpad, **not** from the curation's stored render.

Sanity note on that render: my 600-dpi rasterisation of page 13 is **bit-identical** to the
curation's `00_pages/png/p13_600.png` (`np.array_equal` → True, max abs diff 0). So pixel
coordinates quoted here and pixel coordinates in `04_fig5c.json` live on the same grid and are
directly comparable. The figure is a flattened raster: a single embedded 1541x1864 PNG
(xref 491) placed at 108.51,240.02–478.30,687.29 pt, i.e. **300 dpi native**. The 600-dpi page
render is an exact 2x of it (scale 1.99974 in x, 1.99957 in y). Every measurement below was
repeated on the native 300-dpi image as a cross-check.

Review overlay (ticks, stored tops, re-measured tops, flank windows drawn on the figure):
`/Users/nathanliu/LNP/confound/provenance/review/round_2026-10-30/overlay_XD_2026_fig5c.png`

---

## VERDICT

**CONFIRMED.** F-020 is real, the stated mechanism is correct, and the stated magnitude is right
to within ~0.01 decades. Three details in the finding as written need correcting, all of which
make the finding *stronger*, not weaker:

1. The 60%-rule reproduction is **12 of 12 at exactly 0 px**, not "10 of 12 to ±1 px". There are
   no exceptions to explain. The mechanism is pinned, not merely correlated.
2. The flank-vs-stored gap is **0.3–11.3 native px (0.5–21.4 px at 600 dpi)**, not "0.5–10.9 px";
   the quoted range is right if read in native 300-dpi figure pixels, which is also the unit in
   which "90.00 px per decade" is correct. Stating the unit would avoid a 2x ambiguity.
3. Mean bias is **+0.051 decades** (+0.053 with sub-pixel edges), not +0.060. Max is **+0.117**
   (E4), close to the quoted +0.122. The z-score rms of **0.034** and all-12-rows-change are
   exactly as reported.

The companion claim is also **CONFIRMED**: the axis calibration is sound. This is **not** error
class 20. The two errors do not compound, because there is only one error.

---

## AXIS CALIBRATION

Method: the y-axis spine is the longest vertical dark run in the panel, at page x = 1063
(600 dpi). Ticks are dark runs immediately left of the spine. Major ticks are ~11 px long,
minor ticks ~3 px, so the two populations separate cleanly by run length and I located majors
in the window x ∈ [spine−13, spine−7] (which only long ticks reach) and all ticks in
[spine−5, spine−2]. Centroids are ink-intensity-weighted (sub-pixel), not threshold midpoints.

### Major ticks (600-dpi page coordinates)

| decade | tick y (600 dpi) | spacing to next | px/decade |
|---|---|---|---|
| 10^9 | 4931.339 | 179.824 | 179.824 |
| 10^8 | 5111.163 | 179.971 | 179.971 |
| 10^7 | 5291.133 | 180.081 | 180.081 |
| 10^6 | 5471.215 | — | — |

Linear fit (log10 value vs y): `log10 = -0.0055567918 * y + 36.40205401`

- **px per decade = 179.9600 at 600 dpi = 89.9800 in native 300-dpi figure pixels.**
- Residuals of the 4 majors, in decades: `+0.000367, -0.000387, -0.000327, +0.000347`
  (max 3.9e-4 decades = 0.07 px). Not literally zero, but 3 orders of magnitude below the
  defect being measured.
- The curation's stored slope is -0.005557093 (179.9502 px/decade). Difference from mine:
  3e-7 decades/px, i.e. <2e-4 decades over the whole 740-px plot. Immaterial.

Repeating the whole fit on the **native 300-dpi** image independently:
major ticks at native y = 1465.438, 1555.208, 1645.184, 1735.403; spacings 89.771, 89.976,
90.219; **89.9873 px/decade**; residuals ≤1.3e-3 decades.

### Minor-tick validation

The fit predicts the log minor ticks (2,3,…,9 within each decade). I detected **35 tick
centroids** in the panel; 31 of them are minors falling inside the plotted range, and the fit
predicts all 31:

| decade block | minors checked | residual range (px, 600 dpi) |
|---|---|---|
| 10^9 → 10^8 | 9e8 … 2e8 (8) | −0.49 … +0.48 |
| 10^8 → 10^7 | 9e7 … 2e7 (8) | −0.48 … +0.84 |
| 10^7 → 10^6 | 9e6 … 2e6 (8) | −0.44 … +0.52 |
| below 10^6 | 9e5 … 3e5 (7) | −0.48 … +1.45 |

Minor-tick residual: **rms 0.53 px, max 1.45 px = 0.008 decades**. (2e5 is predicted at
y=5596.9, below the plotted range; there is no tick there and it is excluded.)

**Confirms 90.00 px/decade (89.98 measured) with effectively zero residual, validated against
31 unlabelled minor ticks rather than five.** The calibration is sound; the error is entirely in
bar-top detection.

---

## BAR TOPS

Method, immune to the replicate discs: for every column of a bar I find the first scan row at
which the lavender-fill test holds for 6 consecutive rows. Occlusion can only *remove* fill, so
an occluded column reports a top that is too low; it can never report one that is too high.
I take the modal value over the **5 outermost columns on each side** (skipping 1 px of
anti-aliased edge) and require the two sides to agree. Columns in the middle are never used.

Flank agreement, all 12 bars: **left and right integer tops differ by 0 px on 11 bars and 1 px
on SM-102**. Sub-pixel (50%-crossing of the blue-minus-red fill ramp), left vs right: max
disagreement **0.56 px** (E5), median 0.11 px. The "flanks agree to ≤1 px on all 12" claim holds.

`my pixel top` below is the integer flank top, directly comparable to the stored integer top.
Values use my own calibration. (Using the curation's stored calibration instead changes every
z-score by <1e-5.)

| lipid | stored top_y | my top_y (L / R) | Δ px (600 dpi) | Δ px (native) | stored value | my value | Δ decades |
|---|---|---|---|---|---|---|---|
| ALC-0315 | 4966 | **4953** (4953 / 4953) | −13 | −6.5 | 6.4064e+08 | 7.5729e+08 | **+0.0727** |
| SM-102   | 4971 | **4964** (4965 / 4964) | −7  | −3.5 | 6.0093e+08 | 6.5787e+08 | **+0.0393** |
| MC3      | 5089 | **5082** (5082 / 5082) | −7  | −3.5 | 1.3277e+08 | 1.4536e+08 | **+0.0393** |
| E1       | 5008 | **4990** (4990 / 4990) | −18 | −9.0 | 3.7429e+08 | 4.7170e+08 | **+0.1004** |
| E2       | 5039 | **5037** (5037 / 5037) | −2  | −1.0 | 2.5174e+08 | 2.5852e+08 | **+0.0115** |
| E3       | 5189 | **5180** (5180 / 5180) | −9  | −4.5 | 3.6930e+07 | 4.1483e+07 | **+0.0505** |
| E4       | 5280 | **5259** (5259 / 5259) | −21 | −10.5| 1.1526e+07 | 1.5097e+07 | **+0.1172** |
| E5       | 5229 | **5214** (5214 / 5214) | −15 | −7.5 | 2.2136e+07 | 2.6850e+07 | **+0.0838** |
| E6       | 5509 | **5491** (5491 / 5491) | −18 | −9.0 | 6.1533e+05 | 7.7573e+05 | **+0.1006** |
| E7       | 5527 | **5527** (5527 / 5527) | 0   | 0    | 4.8874e+05 | 4.8940e+05 | +0.0006 |
| E8       | 5462 | **5462** (5462 / 5462) | 0   | 0    | 1.1228e+06 | 1.1242e+06 | +0.0006 |
| E9       | 5494 | **5494** (5494 / 5494) | 0   | 0    | 7.4553e+05 | 7.4652e+05 | +0.0006 |

Δ decades: **mean +0.0514, max +0.1172 (E4), rms 0.0658, min +0.0006**. Every Δ is positive —
the bias is strictly one-sided, as claimed.

Two independent checks on these tops:

- **Native 300 dpi.** Repeating the flank measurement on the embedded 1541x1864 PNG and mapping
  back to 600-dpi page coordinates reproduces every top to **≤1.1 px**: 4951.6, 4963.6, 5081.5,
  4989.6, 5035.6, 5179.5, 5257.5, 5213.5, 5491.5, 5525.4, 5461.5, 5493.5. The 2x upsample
  introduces no bias.
- **Replicate markers.** The bar height should be the mean of the three plotted replicates. The
  ink-intensity centroid of the three discs per bar (a proxy for their mean) agrees with my
  flank tops to **rms 0.021 decades**, and with the stored tops only to **rms 0.057 decades** —
  and the sign of the stored disagreement is positive on 9 of 12 bars (the biased ones), i.e.
  the stored bars sit systematically *below* their own replicate cloud.

Sub-pixel variant (50%-crossing instead of first-fill-row): tops are a uniform 0.40–0.62 px
*above* the integer tops. This is a constant edge convention, so it cancels exactly in the
z-score (max |Δz| between the two variants is 0.0018) and shifts Δ decades only from +0.0514
mean / +0.1172 max to **+0.0534 mean / +0.1193 max**.

---

## MECHANISM

**Yes, completely.** Re-applying the curation's own rule — "first scan row where more than 60%
of the inner bar width (`x0+6 : x1-6`) is lavender fill" — to my independent render reproduces
`04_fig5c.json`'s `top_y` on **12 of 12 bars, to 0 px**. Not 10 of 12; not ±1 px. There are no
exceptions, so the question "which two do not, and why" has no cases: the rule *is* the error.

The reason the 9 biased bars and the 3 clean bars both fall out of the same rule is the
occluded fraction at the true top row:

| lipid | inner width | fill at true top | fraction | >60%? | fraction at stored top |
|---|---|---|---|---|---|
| ALC-0315 | 35 | 15 | 0.43 | no | 0.69 |
| SM-102 | 33 | 7 | **0.21** | no | 0.73 |
| MC3 | 34 | 17 | 0.50 | no | 0.62 |
| E1 | 35 | 20 | **0.57** | no | 0.63 |
| E2 | 34 | 15 | 0.44 | no | 0.62 |
| E3 | 35 | 13 | 0.37 | no | 0.63 |
| E4 | 34 | 15 | 0.44 | no | 0.62 |
| E5 | 35 | 9 | **0.26** | no | 0.66 |
| E6 | 35 | 17 | 0.49 | no | 0.63 |
| E7 | 35 | 34 | 0.97 | **YES** | 0.97 |
| E8 | 34 | 23 | 0.68 | **YES** | 0.68 |
| E9 | 35 | 25 | 0.71 | **YES** | 0.71 |

The three replicate discs are ~18 px across at 600 dpi on a ~46 px bar, they are drawn in the
middle of the width, and they straddle the top edge — so at the true top row they knock the fill
fraction down to 0.21–0.57 on nine bars. The detector walks down the bar until enough of the
discs are behind it, which is 2–21 px (0.011–0.117 decades) too far. On E7/E8/E9 the discs
happen to sit low enough, or far enough off the top edge, that the fraction at the true top is
already 0.68–0.97, so the rule fires at the right row by luck.

E1 is the warning: 0.57 versus a 0.60 threshold, a 3-percentage-point miss, yields an 18 px
(0.100 decade) error. The threshold has no margin.

Per-column top profiles make the picture unmistakable — e.g. E4, left to right across the bar:
`5259 ×14, then 5266 5266 5267 5268 5272 5279 5282 ×12, then 5259 ×15`. Flat flanks, a sagging
middle, nothing else.

One trap worth recording for whoever writes the fix: a naive `min` over *all* columns is **not**
safe. On E1, columns 23–25 report a top of 4982 — 8 px *above* the real top — because the
anti-aliased halo of the disc that floats above the bar passes the lavender colour test. Flank
columns, or a mode, avoid this; a global minimum does not.

---

## PROPAGATION

`stage07_emit.py` computes `Experiment_value = z(log10(value))` over the 12 lipids with
population sd (ddof=0). I reproduced the stored z-scores from `04_fig5c.json` exactly — all 12
agree with `07_XD_2026_rows.csv` to 0.0e+00 — then recomputed with my values.

Stored: mean log10 = 7.338150, pop sd = 1.174654.
Mine:   mean log10 = 7.389581, pop sd = 1.186097.

| lipid | stored Experiment_value | my Experiment_value | Δz |
|---|---|---|---|
| ALC-0315 | +1.250122703619 | +1.255954 | **+0.0058** |
| SM-102 | +1.226468537042 | +1.204419 | **−0.0220** |
| DLin-MC3-DMA | +0.668230205813 | +0.651597 | **−0.0166** |
| E1 | +1.051427704368 | +1.082611 | **+0.0312** |
| E2 | +0.904771871588 | +0.862419 | **−0.0424** |
| E3 | +0.195146874264 | +0.192473 | **−0.0027** |
| E4 | −0.235358957446 | −0.177637 | **+0.0577** |
| E5 | +0.005913541644 | +0.033185 | **+0.0273** |
| E6 | −1.318719786694 | −1.264543 | **+0.0542** |
| E7 | −1.403874786373 | −1.433200 | **−0.0293** |
| E8 | −1.096370620866 | −1.128680 | **−0.0323** |
| E9 | −1.247757286961 | −1.278598 | **−0.0308** |

- **rms |Δz| = 0.0335** (the finding's "rms 0.034" — exact).
- **max |Δz| = 0.0577, on E4.**
- **All 12 rows change**, including E7, E8 and E9, whose *pixel tops do not change at all*.
  This is the point the finding makes and it is correct: the bias is one-sided and non-constant,
  so it moves the experiment mean and sd and therefore reaches every row. Note the sign flip —
  the three unbiased bars move *down* in z (−0.029 to −0.032) purely because the nine biased
  bars moved up. A constant pixel offset would cancel; this one does not.

Sensitivity: using sub-pixel tops instead of integer gives rms 0.0331 / max 0.0578; using the
curation's stored calibration instead of mine gives rms 0.0331 / max 0.0578. The result does
not depend on either choice.

---

## PAPER'S OWN NUMBERS

### Fig 5D reality bands

`stage07_emit.py` asserts `band(v) == CLASS[n]` with A > 5e7, B > 5e6, C otherwise.
**All 12 still fall in their stated band under my values — 12/12, the assertion still passes.**

| lipid | Fig 5D class | stored value → band | my value → band | margin to nearest boundary |
|---|---|---|---|---|
| ALC-0315 | A | 6.406e+08 → A | 7.573e+08 → A | 1.18 dec |
| SM-102 | A | 6.009e+08 → A | 6.579e+08 → A | 1.12 dec |
| DLin-MC3-DMA | A | 1.328e+08 → A | 1.454e+08 → A | 0.46 dec |
| E1 | A | 3.743e+08 → A | 4.717e+08 → A | 0.97 dec |
| E2 | A | 2.517e+08 → A | 2.585e+08 → A | 0.71 dec |
| **E3** | **B** | 3.693e+07 → B | **4.148e+07 → B** | **0.081 dec** |
| E4 | B | 1.153e+07 → B | 1.510e+07 → B | 0.52 dec |
| E5 | B | 2.214e+07 → B | 2.685e+07 → B | 0.27 dec |
| E6 | C | 6.153e+05 → C | 7.757e+05 → C | 0.80 dec |
| E7 | C | 4.887e+05 → C | 4.894e+05 → C | 1.01 dec |
| E8 | C | 1.123e+06 → C | 1.124e+06 → C | 0.65 dec |
| E9 | C | 7.455e+05 → C | 7.465e+05 → C | 0.82 dec |

No band is crossed. **But flag E3**: correcting the bias moves it from 0.132 to **0.081 decades**
below the B/A boundary — 14.6 px at 600 dpi, the thinnest margin in the panel. The two-figure
cross-check that `stage07_emit.py` relies on as its proof of correct digitisation ("two figures
agreeing is the check that the digitisation is right") is, for E3, worth about 15 px. It is a
weak check, and it passed equally well on the *wrong* values — which is precisely why it did not
catch F-020. That sentence in the `stage07_emit.py` docstring overstates what the check buys.

### Printed Fig 5D / 5E statistics

Both reproduced by reading the Fig 5D grid algorithmically (lavender-cell detection over the
12-row x 3-column Reality and Prediction blocks, page x 2117–2591 and 2664–3138, y 4928–5662):

- Reality row = `AAAAABBBCCCC` — **exactly the `CLASS` string hard-coded in `stage07_emit.py`.**
- Prediction row = `AAAAACCBCCBC` (E3→C, E4→C, E5→B, E6→C, E7→C, E8→B, E9→C).
- **Overall accuracy 9/12 = 75.0%, incorrect 3/12 = 25.0%** — matches the Fig 5E pie (75% / 25%)
  and the text on page 5066: *"we achieved an overall prediction accuracy of 75 percent."*
- **Class B accuracy 1/3** — matches the text: *"the predictive accuracy for class B ionizable
  lipids is subpar, with only one-third of the predictions being correct"*, and the two class-B
  misses are both B→C, matching *"the consistent misclassification from class B lipids to class
  C."* Class A 5/5 (text: *"a high degree of concordance ... for class A"*), class C 3/4.

The paper's own numbers are internally consistent and the curation transcribed the Reality
column correctly.

---

## THE EXACT CHANGE

**File:** `/Users/nathanliu/LNP/confound/curation/xd_2026_lipidai/stage04_fig5c.py`

**Lines to replace: 65–67** (inside `for nm,(x0,x1) in zip(names,segs):`), currently

```python
    col=fill[:, x0+6:x1-6]
    rows=np.nonzero(col.sum(1) > 0.6*(x1-x0-12))[0]
    top=rows.min()
```

**New rule** — measure the top from the flank columns only, and require the two flanks to agree:

```python
    # F-020: the three replicate discs are drawn in the MIDDLE of the bar and
    # straddle its top edge, so a "60% of the width is fill" test fires 2-21 px
    # INSIDE the bar. Occlusion can only remove fill, never add it, so the
    # unoccluded flanks carry the true edge. Measure there and cross-check them.
    def _coltop(x, K=6):                      # first row of K consecutive fill px
        c = fill[:, x]; run = 0
        for y, v in enumerate(c):
            if v:
                run += 1
                if run == K: return y - K + 1
            else: run = 0
        return None
    def _mode(v):
        v = [t for t in v if t is not None]
        c = Counter(v); m = max(c.values())
        return min(k for k, n in c.items() if n == m)
    L = _mode([_coltop(x) for x in range(x0+1, x0+6)])   # 5 leftmost columns
    R = _mode([_coltop(x) for x in range(x1-6, x1-1)])   # 5 rightmost columns
    assert abs(L - R) <= 1, ('flank disagreement', nm, L, R)
    top = min(L, R)
```

plus `from collections import Counter` next to the existing imports on line 7
(`import numpy as np, json`). The 1-px inset at each end skips the anti-aliased bar edge; the
`assert` is the guard that makes the method self-checking, and it passes on all 12 bars
(max flank disagreement 1 px, on SM-102 only).

Do **not** replace the rule with `np.nonzero(fill[:, x0:x1].sum(1) > 0)[0].min()` or any other
global minimum: on E1 the halo of the disc floating above the bar passes the lavender test and
would put the top 8 px too high.

**What it gives on all 12 bars** (verified by running the replacement standalone against the
same render, repo file untouched):

| lipid | top_y | value (p/sec/cm^2/sr) | log10 | new Experiment_value |
|---|---|---|---|---|
| ALC-0315 | 4953 | 7.5729e+08 | 8.879 | +1.255954 |
| SM-102 | 4964 | 6.5787e+08 | 8.818 | +1.204419 |
| MC3 | 5082 | 1.4536e+08 | 8.162 | +0.651597 |
| E1 | 4990 | 4.7170e+08 | 8.674 | +1.082611 |
| E2 | 5037 | 2.5852e+08 | 8.412 | +0.862419 |
| E3 | 5180 | 4.1483e+07 | 7.618 | +0.192473 |
| E4 | 5259 | 1.5097e+07 | 7.179 | −0.177637 |
| E5 | 5214 | 2.6850e+07 | 7.429 | +0.033185 |
| E6 | 5491 | 7.7573e+05 | 5.890 | −1.264543 |
| E7 | 5527 | 4.8940e+05 | 5.690 | −1.433200 |
| E8 | 5462 | 1.1242e+06 | 6.051 | −1.128680 |
| E9 | 5494 | 7.4652e+05 | 5.873 | −1.278598 |

Band assertion in `stage07_emit.py`: **12/12 still agree with Fig 5D**, so the downstream
assert continues to pass unchanged. `07_XD_2026_rows.csv`, `07_XD_2026_values.csv` and
`07_qc.json` would all need regenerating (not done here — read-only review).

Two optional follow-ups, both out of scope for the defect itself:

1. The docstring on lines 1–6 says "bar tops are the topmost row of the lavender bar fill",
   which is what the code *should* do and is not what it does. Fix the prose with the code.
2. The claim in `stage07_emit.py`'s docstring that Fig5C-vs-Fig5D agreement "is the check that
   the digitisation is right" should be softened: the check passed on the defective values, and
   its tightest margin (E3) is 0.08 decades.

---

## CONFIDENCE

**Certain** on: the axis calibration (89.98 px/decade native, residuals ≤4e-4 decades on the
majors and ≤0.008 decades on 31 minors, confirmed independently at 300 and 600 dpi); the
mechanism (the stored 60% rule reproduces all 12 stored tops to 0 px, and the fill fraction at
the true top is <0.60 on exactly the nine biased bars and >0.60 on exactly the three unbiased
ones); the 12 re-measured tops (left/right flanks agree to ≤1 px integer and ≤0.56 px sub-pixel
on every bar, and the native-resolution re-measurement reproduces them to ≤1.1 px); the
propagated z-deltas (rms 0.0335, max 0.0577 on E4); and both printed Fig 5D/5E statistics
(75%/25% and class-B 1/3).

**Probable** on: the exact third decimal of the bar values, which depends on the edge convention
(first-fill-row vs 50%-crossing) to about 0.5 px = 0.003 decades. This does not propagate: the
offset is constant and cancels in the z-score to <0.002.

**Speculative** on: nothing load-bearing. The replicate-disc centroid cross-check (rms 0.021 dec
against my tops vs 0.057 against the stored tops) is a supporting argument, not a measurement —
the three discs overlap and I used an ink-weighted centroid rather than resolving them
individually. It points the right way but I would not quote it as a number.
