# Round-2 adversarial review — `QL_2026` (Liu et al., Adv Healthc Mater 2026, PMID 41145968)

**Provenance.** Returned as text by a read-only reviewer (the harness blocks
subagents from writing here) and saved verbatim. **Its findings are claims.** My
adjudication — what reproduced, what I applied, what I held — is in `FINDINGS.md`
under F-031…F-037, not here. Note it independently re-derived two items already
logged (its E-4/E-5 ≈ the round-1 QL E-2/E-1) and **disagrees with the stored
magnitude on one of them**, which is recorded in my adjudication.

---

## SUMMARY — verdict on the 1.0071 offset

**Stop quoting 1.0071. It is not a measurement, and my independent re-measurement gives the reciprocal.** Re-reading all four benchmark bars in both panels with my own full-width-rule detector, the four 3c/2c ratios are 1.0053 / 0.9974 / 0.9605 / 1.0090 — geometric mean **0.9929**, log-sd 0.0226, s.e.m. 1.1%. The curation's 1.0071 and my 0.9929 straddle 1.0 and differ by more than the correction either one applies. Three further facts close it: (i) the applied 0.71% correction is **0.2 of one native pixel** on Fig 3c, where one native pixel is 1.5%; (ii) the bar-top **reading convention** biases the ratio by **+1.8% to +3.1%** on its own, because Fig 2c's pen is thicker than Fig 3c's (top rule median 10 px vs 6 px, side strokes 8 px vs 7 px) and both panels are read at the rule's *top edge* rather than its centre — so the convention artefact is 2.5–4.4× larger than the "offset" measured through it; (iii) removing the division entirely changes the emitted corpus by **rms Δz 0.0009, max |Δz| 0.0047, Spearman 1.000000, zero rank changes**. A 4-anchor offset *is* defensible as a **scale check** — it has the power to detect a real inter-run offset of 1.5× and it returned 1.00, which is what licenses pooling the two panels — but it cannot measure a 0.7% number. The honest statement is "the four shared benchmarks place the two screens on a common scale to within ±2.5%, consistent with 1.0", and the division should either be dropped or kept with that sentence attached. **One concrete hazard: `validate.py` line 154 asserts `1.0 < st.mean(rr) < 1.3`. It passes by 0.0071 on a quantity with sd 0.0247; my equally careful re-measurement would fail it.** That `1.0 <` floor is inherited from the superseded 1.1116 era and was never re-derived after F-018 removed the effect — error class 35 one level up, at the threshold instead of the window. Beyond the offset: the data chain is in excellent shape — I reproduced all 344 structures from an independent enumeration, both bar panels to ≤1 native pixel, the heatmap's saturated and floored sets exactly, and the corpus z-scores to 9e-16 — and the three round-1 repairs (F-017, F-018, F-019) are all genuinely present in the emitted data.

---

## 1. CONFIRMED CORRECT

Method note: my 600-dpi page renders and the curation's `00_pages/*.png` are **the same raster**, offset by exactly (dx 486, dy 640) for p4 and (486, 639) for p5 — established from identical tick centroids, so every number below is directly comparable to the stored constants.

**C-1. Fig 3c log-axis calibration (F-017 / class-20 repair holds).** Five y-axis tick marks, darkness-weighted centroids in my render: 2946.46 / 3098.04 / 3249.59 / 3399.73 / 3550.20 → **150.94 px/decade**, 10³ at 3550.20. Stored: 151.0 and 10³ at 2911.5 (= 3550.5 in my frame). Agreement 0.3 px and 0.04%. Anchored on ticks, not the panel box.

**C-2. Fig 3c bar tops (F-018 repair holds).** All 28 bars re-measured with a detector that requires one contiguous ink run spanning the **full** bar width (x0+3 → x1−4) for 2 rows. **22 of 28 identical to the pixel; 6 differ by exactly one native (300 dpi) pixel.** Max value difference 3.5%, median |diff| 0.2%, no systematic sign. **No bar is read on the replicate dot cluster** (the old failure was up to +28%). Visual overlay confirms every red box sits on the rule, not the dots.

**C-3. Fig 3c library 2 — all 24 re-derived independently** (mine vs stored RLU):

| lipid | mine | stored | lipid | mine | stored |
|---|---|---|---|---|---|
| T3A-L1-2C10 | 904,148 | 905,636 | T3A-L1-2C11 | 222,180 | 222,681 |
| T3A-L2-2C10 | 384,786 | 373,979 | T3A-L2-2C11 | 32,502 | 32,602 |
| T3A-L3-2C10 | 222,180 | 222,681 | T3A-L3-2C11 | 88,958 | 89,193 |
| T3A-L4-2C10 | 124,435 | 124,746 | T3A-L4-2C11 | 5,054 | 5,231 |
| T3A-L5-2C10 | 1,303,912 | 1,305,852 | T3A-L5-2C11 | 708,327 | 731,541 |
| T3A-L6-2C10 | 136,362 | 136,698 | T3A-L6-2C11 | 434,732 | 422,501 |
| T3A-L7-2C10 | 825,063 | 826,454 | T3A-L7-2C11 | 174,060 | 174,470 |
| T3A-L8-2C10 | 283,603 | 284,212 | T3A-L8-2C11 | 554,917 | 573,163 |
| T3A-L1-2C12 | 810 | 814 | T3A-L5-2C12 | 185,013 | 185,444 |
| T3A-L2-2C12 | 13,013 | 13,059 | T3A-L6-2C12 | 1,154,105 | 1,155,883 |
| T3A-L3-2C12 | 51,365 | 51,513 | T3A-L7-2C12 | 236,160 | 236,686 |
| T3A-L4-2C12 | 11,518 | 11,559 | T3A-L8-2C12 | 554,917 | 555,947 |

**`T3A-L2-2C10` specifically** (round 1's 6.2%-off case, "honest value y = 2521 ± 2"): my detector returns **y = 2521** in the curation's frame, exactly round 1's figure. The stored value is y = 2523, one native pixel low, 2.8% in value — i.e. that item is now closed at the resolution limit, not 6.2% off.

**C-4. Fig 2c calibration, validated by a round number the authors drew but never printed.** Tick centroids measured in the curation's own raster: 2343.87 / 2750.16 / 3157.02 (stored comment: 2342.78 / 2749.54 / 3156.75; ≤1.1 px). The red dotted reference line's centroid is 2872.74 (stored 2872.21). It evaluates to **50,019 RLU on my constants and 50,021 on theirs — 0.04% off a round 5×10⁴ either way.** This is a genuinely external check and it is decisive; the 0.1% disagreement between 406.58 and 407.0 px/decade is immaterial.

**C-5. Fig 2c bars.** 32/32 runs detected at the same windows; **28 identical to the pixel, 4 differ by one native pixel**; max value difference 2.0%, no systematic sign.

**C-6. Fig 2b geometry.** Grid auto-detected independently: **41 vertical lines** (833.5 → 3784.5, pitch 73.775) and **9 horizontal** (1721.5 → 2662.5, pitch 117.625) — identical to the stored constants after the frame offset. Row labels read off the figure are 70 / 75 / 80 / 400 / 113 / 114 / 110 / T3A in that order, label centroids on row centres (no class-16 off-by-one). Column brackets L1…L8 with C6/C8/C10/C12/C14 inside each confirmed from the printed axis.

**C-7. Fig 2b cell reading.** All 320 cells re-read with my own LUT and my own tick-derived value scale. **Saturated set identical: 21 of 21, same lipids. Floored set identical: 113 of 113, same lipids.** 45 distinct values with the identical tie structure (65 / 57 / 48 / 22 / 21 / 16 / 16 / 10). Values differ only by the calibration offset discussed in E-3.

**C-8. Fig 2b ↔ Fig 2c agreement.** Over the 10 shared lipids that are neither saturated nor floored, median heatmap/bar ratio **1.004** (the curation reports 0.997 over n=9, which is the reciprocal with the conflict row excluded — consistent). The `T3A-L8-C14` conflict reproduces exactly: heatmap 7,587 vs Fig 2c 57,874 = **7.6×**.

**C-9. All 344 structures.** I wrote an independent enumerator (my own per-linker SMILES templates, my own core templates, no shared code with `enumerate_library.py`) and it reproduces **344 of 344 canonical SMILES exactly**. Thiosuccinimide regiochemistry is correct and in fact has no freedom: maleimide's two ring carbons are equivalent, so thia-Michael addition can only give 3-(alkylthio)succinimide, which is what is stored. L2 correctly gives a **disulfide** (`R-S-S-CH₂CH₂C(=O)NH-`) via pyridyl-dithio exchange, not a thioether.

**C-10. Linkers L1–L8 confirmed against a drawing the curation never cites.** **Figure 4a** (p. 6) draws all eight BLs in full *and* draws the complete six-tail product. Read at 600 dpi: L1 = SMCC (NHS-CO-cyclohexane-1,4-diyl-CH₂-maleimide), L2 = SPDP, L3 = 1 CH₂, L4 = 3 CH₂, L5 = 5 CH₂, L6/L7/L8 = 1/2/3 ethylene-glycol units. **All eight match the stored block SMILES.** Fig 4a's product also independently confirms one amide per primary amine on T3A with the central tertiary N untouched.

**C-11. Tail carbon counts, geometrically.** 2C10 = C19 is pinned by mass (all eight Table S2 rows are T3A-Lx-2C10 and my independent build matches: e.g. T3A-L1-2C10 1747.83 vs printed 1747.8; two C10 arms would be +84 Da). Using that to fix the bond pitch from Fig 3a (50.7 px/bond), the three drawn upper arms measure **9.03 / 9.96 / 10.99 bonds** → 2C10/2C11/2C12 = C19/C21/C23, exactly as stored. Counted by geometry, not by eye.

**C-12. Independent mass check vs Tables S1/S2.** My own build: **23 of 31 reference masses agree (max 1.45 Da)**; the eight that do not are the six the curation already records plus the two from the duplicated label block. I independently resolve two of them: the "113-L3-C6 613.8" row is **114**-L3-C6 (613.85, a unique fit — confirming PROVENANCE's own correction of round 1), and the second "L6" block in Table S1 is **L7** (400-L7-C6 = 516.75 vs printed 516.7; T3A-L7-C6 = 1260.73 vs 1260.7). The five L4/L5 anomalies all fit a single reading — the **tail** token is mistyped and should be C8 in every one — which is consistent with Fig 4a's drawn linker lengths and leaves every enumerated structure unchanged.

**C-13. Arm counts pinned.** Table S1's `400-L6-C6` at 471.7 vs my one-arm 472.70 (two arms would be 745.0) confirms a secondary NH does not acylate; 114's two primary amines therefore give exactly two arms. Reproduced independently.

**C-14. The four benchmark LNPs are absent from the corpus** — 0 of 344 `IL_name`s match, and none of the three benchmark SMILES the corpus holds elsewhere appears in the slice. Excluding them is right: they are calibration standards, not library members, and emitting them under `QL_2026` would put the field's universal positive controls into a screen whose chemistry is otherwise entirely novel.

**C-15. Fig 2d/2e/2f do not follow from Fig 2b — and it is not our calibration.** With *my* values: at the 2×10⁴ threshold the per-amine counts are [0,1,3,2,3,5,8,19] against the printed [0,2,2,1,4,2,8,17]. The disagreement is **non-monotone** (114 is over by 3, 75 and 113 are under by 1), so no threshold reconciles them — I swept 1.5×10⁴ to 4×10⁴. Tail counts [3,7,12,8,11] vs printed 2f [3,6,13,8,7]. The curation's "paper-internal inconsistency, recorded not reconciled" is independently confirmed.

**C-16. The C10 anchor is real.** `T3A-L5-C10`: reference string 125 bytes / sha `1a0781…`, enumeration 123 bytes / sha `06fe7e…` — **not a byte copy** — identical InChIKey `DKFSQTNBZODLDT-UHFFFAOYSA-N`, and its MolWt 1290.98 matches Table S1's printed 1290.9 to 0.08 Da. (See E-9 for what it does not cover.)

**C-17. Metadata, against the Methods on p. 11.** HeLa 10⁴ cells/well in 96-well; lipid/DOPE/Chol/DMG-PEG2000 = 35/16/46.5/2.5; LNP/mRNA 20:1 w/w; 20 ng FLuc mRNA per well; 24 h; citrate pH 4.0; ~300× pipette mixing; crude product for in vitro, flash-purified >95% only for in vivo. All match the emitted rows. Leaving `IL_to_nucleicacid_massratio` **blank** is correct — the paper says "LNP/mRNA", not ionizable-lipid/mRNA.

**C-18. The value chain reproduces end to end.** Rebuilding z from `01_`/`05_`/`07_` + the offset gives `data_vitro.csv`'s `Experiment_value` to **8.9e-16**; all 344 SMILES identical between `11_QL_2026_rows.csv` and the corpus; `Experiment_value_within_model` == `Experiment_value` (single `Model_type`), `value_normalization_unit = experiment` on all 344, molar percents sum to 100.

---

## 2. ERRORS FOUND

### E-1 — 1.0071 is a convention artefact, not a measured offset  *(this is the manuscript number)*

| | |
|---|---|
| file | `curation/ql_2026/build_rows.py` (offset), `12_QL_2026_qc.json`, PROVENANCE.md |
| severity | **value** (negligible in the data) / **high** as a reported quantity |
| confidence | certain |

Four independent lines of evidence:

1. **Re-measurement flips the sign.** My four ratios: 1.0053 / 0.9974 / 0.9605 / 1.0090, geometric mean **0.9929** (log-sd 0.0226). Stored: 1.0114 / 1.0345 / 0.9670 / 1.0153, mean 1.0071 (sd 0.0247). The single largest per-benchmark difference is **SM102**, where the stored Fig 3c value (51,513.4) is byte-identical to MC3's; I read it one native pixel lower (49,822). The ALC-0315 anchor disagrees by ~4% in *both* readings, so the anchor set is not homogeneous either.
2. **Below resolution.** 0.71% is 0.2 of one native pixel on Fig 3c (1 native px = 1.5%).
3. **Dominated by the reading convention.** Both panels are read at the **top edge** of the bar's top rule, but the data value sits at the stroke **centre**. Fig 2c's pen is thicker than Fig 3c's: top-rule thickness medians **10 px vs 6 px**, side-stroke medians **8 px vs 7 px** (600 dpi). Half-thickness over px-per-decade gives +1.23% … +2.29% for Fig 2c and +1.99% … +5.48% for Fig 3c, so centre-reading would multiply the 3c/2c ratio by **1/1.018 … 1/1.031**, taking 1.0071 to **0.989 … 0.977**. The artefact is 2.5–4.4× the offset it produced.
4. **No effect on the data.** Setting the offset to 1.0: **rms Δz 0.0009, max |Δz| 0.0047, Spearman 1.000000.** Switching the whole slice to stroke-centre reading: rms Δz 0.006, max 0.037, library-2 mean z 2.111 → 2.095.

**Recommendation.** Keep the benchmark comparison, drop the number. Report it as "the four shared benchmarks agree between the two panels to within ±2.5%, consistent with a common scale (1.00); no rescaling applied", and remove the division or state that it is below the digitization resolution. **Verify by:** re-measuring any one benchmark bar top in both panels at 600 dpi and noting that a one-native-pixel choice moves the ratio by 1.5%.

### E-2 — `validate.py`'s offset gate has no margin and is pinned to the superseded value

| | |
|---|---|
| file | `curation/ql_2026/validate.py` line 154 |
| severity | **the check, not the data** |
| confidence | certain |

```python
check(st.pstdev(rr) < 0.1 and 1.0 < st.mean(rr) < 1.3, ...)
```
Measured mean 1.0071 passes the `1.0 <` floor by 0.0071, against a per-anchor sd of 0.0247. My independent measurement (0.9931) **fails this assertion**. The window `[1.0, 1.3]` was plainly chosen when the offset was 1.1116 and was never re-derived after F-018 reduced it to ~1.0 — the same shape as error class 35 ("every threshold calibrated under the old mask has to be re-derived, not inherited"), here at the threshold rather than the search window. A test that an equally valid re-measurement would fail is not testing what it claims. Replace with a two-sided consistency test, e.g. `abs(log10(mean(rr))) < 3*sd/sqrt(4)`.

### E-3 — `Y_VMIN` is not the value its own comment derives; the floor bound is 2.8–4.1% high

| | |
|---|---|
| file | `curation/ql_2026/digitize_fig2b.py`, `Y_VMAX, Y_VMIN = 1107.4, 1967.6` |
| severity | **value**, small, but on exactly the 48 rows F-019 exists for |
| confidence | certain on the arithmetic; probable on which anchor is right |

The comment states the anchors come from the colorbar's 4e4 and 2e4 tick notches at **1393.98 and 1680.60, 143.31 px per 1e4**. Those give y(0) = 1680.60 + 2×143.31 = **1967.22** — and the comment itself says "the tick-derived 1967.23". The code uses **1967.6**, 0.38 px lower; that figure matches the *border-centre* cross-check (1967.5) the comment presents as an independent route, not the ticks it attributes the constants to. This is error class 21, in the same file the class was first logged against.

My own notch centroids, measured darkness-weighted **in the curation's own raster**, are 1394.11 and 1680.58 → px/1e4 = 143.235, **y(0) = 1967.05**.

Consequence, all at the pure-white LUT row y = 1954:

| anchor | y(0) | floor value |
|---|---|---|
| code as written | 1967.6 | **948.6** (stored) |
| code's own stated ticks | 1967.22 | 922.5 |
| my tick centroids | 1967.05 | 911.1 |

So the upper bound carried on the 48 pure-white rows is 2.8–4.1% too high. Negligible at the top of the scale (<0.1%). **Verify by:** recomputing `Y_VMIN` from the two notch centroids in the file's own comment and re-running the LUT.

### E-4 — 65 of the 113 `value_is_bound = upper` rows are not on the floor, and their note is false for them

| | |
|---|---|
| file | `curation/ql_2026/build_rows.py` `row_flags()`; 65 rows of `data_vitro.csv` |
| severity | **metadata**, with direct analysis consequences |
| confidence | certain |

(Independently re-derived; already logged as **QL E-2** in round 1, so treat this as corroboration with the numbers re-measured.)

The 113 split into exactly two LUT bands:

| n | cell RGB | stored RLU | the band the truth occupies |
|---|---|---|---|
| 48 | (254,254,254) / (255,255,255) | 948.6 | **[0, 948.6]** — open downward. Genuine upper bound. |
| 65 | (254,254,252) | 1367.1 | **[1018, 1367]** — bounded below. A measurement with ±15% quantization. |

Their shared `curation_note` says the cell "sits on the colormap floor … every true value from 0 to there renders the same colour". For the 65, that is false: their colour is measurably not white (distance 3.32 from (255,255,255), versus the next band up at 4.90), and a true value of 0 would render white, which they are not. The flag rule is `‖col − white‖ < 4.0`, which straddles the two bands. `upper` should be on **48** rows; the 65 need either no flag or a different one. As it stands, a consumer filtering on `value_is_bound` discards 65 real measurements, or treats a 1,367 RLU reading as possibly zero.

### E-5 — every Fig 2b value is the **top** of its quantization band (median +12.8%)

| | |
|---|---|
| file | `digitize_fig2b.py`, `i = int(d.argmin())` |
| severity | **value**, 186 rows |
| confidence | certain on the mechanism and the bias; the repair is a decision |

(Also already logged as **QL E-1** in round 1; re-derived, and I disagree with the magnitude quoted there.)

Within a band of identical LUT rows every distance ties, so `argmin` returns the lowest index = smallest y = **largest value**. Band widths are wide at the bottom and 1 px at the top, so the bias is one-sided and strongly non-uniform — it does not cancel in a z-score:

| band top | band bottom | n cells | bias of top-of-band |
|---|---|---|---|
| 1,367.1 | 1,018.4 | 65 | +14.6% |
| 2,064.6 | 1,436.9 | 57 | +17.9% |
| 2,762.1 | 2,134.4 | 22 | +12.8% |
| 4,296.7 | 3,250.4 | 16 | +13.9% |
| ≥ 1.2×10⁴ | | 1 each | ≤2.3% |

186 unsaturated, unfloored cells move; median bias **+12.8%**.

**On the magnitude.** Using band **midpoints for the genuinely measured cells only** (leaving the 113 floored rows as bounds, which is what they are): **rms Δz 0.038, max |Δz| 0.066, Spearman 1.000000** (the map is monotone, so no rank changes). Pushing midpoints onto the floored rows too — replacing a bound with a point estimate of 474 RLU — gives rms Δz 0.116, max 0.275, and is the wrong thing to do. The 0.149 / 0.321 figure in the open list looks like the second variant; I would not use it to justify a decision.

### E-6 — `PROVENANCE.md` still states the superseded offset

| | |
|---|---|
| file | `curation/ql_2026/PROVENANCE.md` line 142 |
| severity | **metadata / documentation**, high leverage |
| confidence | certain |

> "The four benchmark LNPs appear in both Fig 2c and Fig 3c, giving a cross-figure ratio of **1.112× (sd 0.023)**. Library 2 is divided by that…"

The applied value is **1.0071 (sd 0.0247)** and has been since 2026-10-30. The file's own round-1 section (line 243) mentions the correction to ≈1.073, which is itself superseded. This is class 21 on precisely the number the manuscript is said to carry — a reader of this folder gets 1.112.

### E-7 — `PROVENANCE.md` line 5 contradicts the corpus and the rest of its own file

> "DOI `10.1002/adhm.202502225`. **Not in PubMed as of 2026-10-29.**"

PMID 41145968 is populated on all 344 corpus rows, and line ~205 of the same file records that the "not indexed" premise was wrong. The header is the first thing read. (Overlaps the already-logged QL E-5 about the stale "Open after merging" section; the specific lines are 5, 142, 217, 229–250.)

### E-8 — the PEG identity contradiction reaches the corpus unflagged

| | |
|---|---|
| severity | **metadata**, 344 rows |
| confidence | certain on the contradiction, unresolvable on which is right |

The Methods and every figure say **DMG-PEG 2000**; the SI Materials list (my own extraction of `word/document.xml`) says *"…cholesterol phosphoethanolamine-N-[methoxy(polyethyleneglycol)-2000 (PEG-2000) were purchased from Sigma-Aldrich"*. The corpus distinguishes the two (DMG-PEG2000 21,392 rows; DMPE-PEG2000 4,901). All 344 rows assert `PEG_name = DMG-PEG2000` with **blank `confidence` and blank `curation_note`**. The contradiction is recorded only in `PROVENANCE.md` — the exact shape §6 names: *a flag that lives only in a per-folder file is a memo*. The schema now has a column that can carry it.

### E-9 — C10 coverage: one typed anchor, no generator, and the paper's best drawing unused

| | |
|---|---|
| severity | **process**, not a current structural error |
| confidence | certain |

The anchor is not a copy and its mass checks out (C-16). But: it is a **typed literal** in `reference_structures.csv` with no `make_reference.py`, so the "different code path" claim is asserted rather than executed — the standard CLAUDE.md §6 now sets. More materially, it was built from "Fig 2a blocks + Methods", i.e. from our own reading of the components, so it does not test that reading against any drawn whole molecule. **Figure 4a draws the complete six-tail product and all eight BLs in full** and is cited nowhere in the anchor record. An anchor derived from Fig 4a would, for free, cover **L1 and L2** — currently the two linkers with no anchor at all, and the two most structurally distinctive (cyclohexane; disulfide rather than thiosuccinimide).

### E-10 — the floored rows' note carries two different numbers for one quantity

The note reads "pure white at about **948** RLU … the stored value is the UPPER BOUND of a **[0, ~982]** RLU band". 982 is the pre-F-019 figure; the stored bound is 948.6 (and per E-3 the tick-derived value is 911–922). Cosmetic, but it is the sentence a `data_vitro.csv` reader sees.

### E-11 — latent error class 30 in `build_rows.py`

`comp.setdefault(name, smiles)` over the live corpus — first row wins, resolved by file order. This is the exact pattern repaired in `curation/build_new_curation.py` (which now raises on an ambiguous name); the repair was not propagated to its sibling. **Currently harmless, verified**: DOPE, Cholesterol and DMG-PEG2000 each map to exactly one SMILES corpus-wide. It becomes live the moment any of those names acquires a second structure.

### E-12 — the 7.3× conflict row is not counted as a caveat (already handled)

`T3A-L8-C14` carries only a `curation_note`; `rows_with_a_caveat` keys on `confidence`/`positional_ambiguity`/`structure_assumed`/`value_is_bound`, so the most-qualified value in the slice is counted as clean. I re-derived this independently before finding that `audit.py` already records it and resolves it with `rows_with_only_a_note`. Noting only to confirm the measurement: **117 rows carry a note, 116 carry `value_is_bound`**, and the 117th row's `curation_note` **is** populated — the parenthetical in `CLAUDE.md` §6 and in round-2 FINDINGS ("carries no flag column at all") is factually wrong; what it lacks is `value_is_bound`. A `confidence = low` on that one row would both flag and count it.

---

## 3. UNVERIFIABLE

- **The true values of the 48 pure-white cells.** Anywhere in [0, ~911–949]. Only the authors' source data settles it.
- **The true values of the 3 censored rows** (`400-L8-C12`, `114-L4-C14`, `114-L7-C14`), stored at 59,051 and flagged `lower`. The saturated cells that Fig 2c *does* rescue reach 377,382, so the censoring is plausibly up to ~6× — the note says "lower bound" without saying how far.
- **Which token is mistyped in the five L4/L5 Table S1 rows.** Mass-degenerate per row; a single coherent reading exists across rows (tail = C8) but mass cannot prove it. Immaterial: all 320 combinations are enumerated, so no stored structure depends on it.
- **DMG-PEG2000 vs a phosphoethanolamine-PEG2000** (E-8). Needs the authors.
- **The thiosuccinimide C3 configuration.** A catalyst-free thia-Michael on a crude, unpurified product is racemic; 301 of 344 molecules carry 1–3 unassigned stereocentres there (the 43 L2 disulfide lipids have none). Stored unspecified, which is correct — flagging it would be noise, not information.
- **The 40 L1 lipids' cyclohexane-1,4-diyl cis/trans.** Genuinely undetermined from the paper and not covered by C9 (which only sees C=C). Arguably a `structure_assumed` candidate; stated rather than pressed.

---

## 4. SEVERITY

| id | changes a structure | changes a value | metadata only | cosmetic |
|---|---|---|---|---|
| E-1 offset is not a measurement | — | yes (rms Δz 0.0009 — negligible in data; **high as a reported number**) | — | — |
| E-2 `validate.py` gate has no margin | — | — | — | the check itself |
| E-3 `Y_VMIN` ≠ stated derivation | — | **yes**, 320 rows (≤4.1%, concentrated on the 48 floor rows) | — | — |
| E-4 `upper` flag wrong on 65 rows | — | — | **yes**, 65 rows | — |
| E-5 top-of-band convention | — | **yes**, 186 rows (rms Δz 0.038) | — | — |
| E-6 PROVENANCE asserts 1.112× | — | — | **yes** | — |
| E-7 "Not in PubMed" header | — | — | **yes** | — |
| E-8 PEG contradiction unflagged | — | — | **yes**, 344 rows | — |
| E-9 C10 coverage | — | — | process | — |
| E-10 948 vs ~982 in one note | — | — | — | yes |
| E-11 latent class 30 | — | — | latent | — |
| E-12 conflict row not a caveat | — | — | yes (already handled) | — |

**Nothing I found changes a structure.** The 344 enumerated molecules survived a from-scratch re-derivation, a 31-row mass check, and a figure (Fig 4a) the curation never used.

---

## 5. CONFIDENCE

| id | confidence | basis |
|---|---|---|
| E-1 | **certain** | four independent arguments; the decisive one (pen-thickness asymmetry, 10 px vs 6 px) is measured, not argued |
| E-2 | **certain** | the literal in the source, against my measured 0.9931 |
| E-3 | **certain** on arithmetic; **probable** on which anchor is right | the comment's own two numbers disagree with the constant; my notch centroids measured in the curation's own raster |
| E-4 | **certain** | the two bands are separable in RGB and their value intervals do not overlap |
| E-5 | **certain** on mechanism and bias; the **repair is a decision** | band structure computed directly from the LUT |
| E-6, E-7, E-10 | **certain** | file text vs emitted data |
| E-8 | **certain** that the paper contradicts itself; **speculative** on which lipid | SI `document.xml` vs Methods |
| E-9 | **certain** | folder contents and the anchor's own `note` |
| E-11 | **certain** that the pattern is present; **certain** that it is currently harmless | checked all three names across both corpus files |
| E-12 | **certain** | direct count of the emitted rows |

---

### Priority if you only act on three

1. **E-1** — rewrite how the offset is described wherever it appears (`build_rows.py` docstring, `12_QL_2026_qc.json`, `PROVENANCE.md` line 142, and the manuscript). The data need not change.
2. **E-2** — replace the `1.0 < mean < 1.3` gate; it currently certifies a number my re-measurement contradicts.
3. **E-4** — the flag fix is unambiguous and cheap: 48 rows `upper`, 65 rows re-described. It can be applied without E-5's decision.

E-3 should ride along with any E-5 decision, since both move the same values and both are concentrated at the bottom of the heatmap's range.
