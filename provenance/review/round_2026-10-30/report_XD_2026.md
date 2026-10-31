# Wave-1 adversarial review — XD_2026 (Duan et al., *Acta Pharm Sin B* 2026, PMID 42592089)

Reviewer: Tier-3 adversarial agent, 2026-10-30. Read-only.
Target: `curation/xd_2026_lipidai/` → 12 rows of `XD_2026` in `data/data_vivo.csv`.
Scope: (a) whether the first Tier-2 pass's nine fixes are actually right, (b) the ground it
did not cover. Previous conclusions were treated as claims throughout.

**Method, stated once.** Everything below was re-derived from the PDF, not from the
curation's scripts. The figure on p.5065 is a single flattened raster, embedded at
**1541 x 1864 px in a 369.8 x 447.3 pt box = exactly 300 dpi native**; I worked on that
raster directly (extracted with `fitz.Document.extract_image`), because a 600-dpi page
render is a 2x interpolation of it and adds no information. Structures were counted by eye
at 8-16x magnification and cross-checked with sub-pixel centreline profiling (a different
estimator from `stage02_trace.py`'s skeletonise -> RDP -> vertex count). Fig 5C was
re-digitised with my own tick calibration and my own bar-top estimator. Fig 5D was read by
pixel-sampling its cell grid. Benchmark SMILES were written from their chemical names by
fragment assembly and compared by InChIKey.

**Headline:** the curation's central and most aggressive claim — that the paper draws
ALC-0315 one carbon short — is **correct**, and I reproduce it independently. All twelve
structures and all metadata are correct. One real value-changing defect was found, in the
Fig 5C bar-top detector. Four documentation defects were found, one of which (the
"299 vs 47" claim) is recorded in two places that contradict each other and is
unreproducible from the data.

---

## 1. CONFIRMED CORRECT

### 1.1 The ALC-0315 claim — verified independently, and it holds

This was the claim most worth attacking, because on its strength the corpus overrides a
published structure. I counted the carbons myself.

Method: column-by-column ink profile of the native raster, taking the stroke centreline
y(x) and locating its turning points. Bond pitch in this panel is 11.5-12.0 px in x.

| | glyph/vertex x (native px) | carbons |
|---|---|---|
| **ALC-0315 head** | O ~= 113 -> **125 (valley), 137 (peak), 148.5 (valley)** -> N ~= 159 | **3** |
| **E1 head** (390 px below, same panel, same renderer) | O ~= 83 -> **95 (peak), 106.5 (valley), 118 (peak), 130 (valley)** -> N ~= 142 | **4** |

Both chains are measured by the identical procedure; both have the same bond pitch; both
terminate on an `N` glyph with a single bond of normal length (12.6 px for ALC, no
intervening turning point between the last carbon and the nitrogen). The real compound is
((4-hydroxybutyl)azanediyl)bis(hexane-6,1-diyl) bis(2-hexyldecanoate) — four. **The paper
draws three.** Confirmed.

The claim that the error is *localised*, not a convention, also holds, and I tested it
harder than the curation did:

- **ALC-0315's own spacers are drawn correctly.** Profiling the upper arm gives
  N(159) - C(159.5, a near-vertical bond) - C(170.5) - C(182) - C(193) - C(204.5) -
  C(216) - O(226.5): **six** carbons, i.e. hexane-6,1-diyl, exactly right. So the renderer
  is not systematically dropping a carbon; the defect is confined to the head.
- **SM-102 is drawn correctly.** N-(CH2)5-C(=O) hexanoyl arm (vertices at native x ~ 123,
  136, 148.75, 161, 173, carbonyl 185); undecyl ester (11 vertices from the ester O);
  octanoate arm (7 CH2 + carbonyl); heptadecan-9-yl alcohol.
- **MC3 is drawn correctly**: acyclic ester (no dioxolane ring), 4-(dimethylamino)butanoate
  (3 CH2 + carbonyl), both arms cis,cis-dienes drawn peak-to-peak.

The repair chosen is also the right half. The name `ALC-0315` has an external authority
(a defined compound), the paper uses it as a benchmark, and two other benchmarks in the
same panel reproduce the corpus exactly. Emitting the corpus structure keeps the value with
the compound the authors say they used. The drawn (3-carbon) form is preserved in
`reference_structures.csv` and surfaces in `provenance/structure_check.json` as
`expected_mismatch_paper_drawing_is_wrong` — the absence is visible, not buried in prose.

### 1.2 The geometric traces of E1-E9 — spot-checked by eye, all agree

I counted **nine** of the nine by eye at 8-16x on the native raster (the brief asked for
four). The specific failure mode the brief names — a vertex-counting algorithm that misses a
vertex at a particular bond angle — is a real risk in this panel, because every tail
approaches its nitrogen through a **near-vertical bond** whose turning point is easy to drop.
I checked that geometry explicitly in E4, E5, E7 and in the ALC-0315 spacer. The trace
counts it correctly every time; it is my *low-magnification* eye that dropped it first (I
initially read SM-102's hexanoyl arm as 5 carbons at 8x before resolving the N-adjacent
vertex at 16x). That is a point in the tracer's favour, not against it.

| lipid | what I counted by eye | agrees with emitted SMILES |
|---|---|---|
| E1 | head HO-(CH2)4-N; two hexane-1,6-diyl linkers (6 regular steps each); two ureas; tail = CH2-CH(hexyl)(octyl), i.e. 8-carbon main + 6-carbon branch | yes |
| E2 | tail CH2-CH(OH)-C12 (12 vertices from the carbinol outward); N(tail)2-CC-NH-CO-CC-piperazine-CC-CO-NH-CC-N(tail)2 | yes |
| E3 | tail 14 C, OH on C2 (12 vertices beyond the carbinol); propyl-piperazine-propyl | yes |
| E4 | all four tails **14 C** (upper-left: 14 vertices x 437->601; lower-left: 13 in-band + 1 at the N-bend); propyl linkers (3 vertices) | yes |
| E5 | tails **14 C** (13 in-band 417->564 plus the N-adjacent vertex at (563, 475)); backbone N(tail)2-CC-NH-CO-CC-N(Me)-CC-N(Me)-CC-CO-NH-CC-N(tail)2 | yes |
| E6 | tail 14 C, OH on C2; backbone = hydroxy analogue of E7 | yes |
| E7 | tails **14 C** (13 zig-zag vertices + the N-bend); bis-amide piperazine backbone | yes |
| E8 | tails 14 C, OH on C2; ethylenediamine-bis-amide backbone | yes |
| E9 | all four tails **12 C** (upper-left 956->1108, upper-right 1329.5->1480.5, both exactly 12 vertices) | yes |

Formula arithmetic independently confirms each: e.g. E4 = 4 x 14 + 2 x 3 + 4 (piperazine) =
C66 ok; E9 = 4 x 12 + 4 x 2 + 2 x 3 + 4 = C66, N8, O4 ok; E1 = 2 x 16 + 2 (ureas) + 12 + 4 =
C50, N5, O3 ok.

### 1.3 The three benchmarks, written from their chemical names

Written by me by fragment assembly from the name alone, canonicalised, compared by InChIKey:

| name | my string -> InChIKey | corpus InChIKey | match |
|---|---|---|---|
| ALC-0315 — ((4-hydroxybutyl)azanediyl)bis(hexane-6,1-diyl) bis(2-hexyldecanoate) | QGWBEETXHOVFQS-UHFFFAOYSA-N | QGWBEETXHOVFQS-UHFFFAOYSA-N | yes; C48H95NO5, 765.7210 |
| SM-102 — heptadecan-9-yl 8-[(2-hydroxyethyl)(6-oxo-6-undecoxyhexyl)amino]octanoate | BGNVBNJYBVCBJH-UHFFFAOYSA-N | BGNVBNJYBVCBJH-UHFFFAOYSA-N | yes; C44H87NO5, 709.6584 |
| DLin-MC3-DMA — (6Z,9Z,28Z,31Z)-heptatriaconta-6,9,28,31-tetraen-19-yl 4-(dimethylamino)butanoate | NRLNQCOGCKAESA-**KWXKLSQISA**-N | NRLNQCOGCKAESA-KWXKLSQISA-N | yes; C43H79NO2, 641.6111 |

The MC3/KC2 trap is cleared explicitly: the corpus string is the **acyclic ester**
(`Rings = 0`), not the KC2 dioxolane, and the stereo layer is the all-cis tetraene.
Both are 641.6111 Da, so only the InChIKey and the ring count discriminate — both do.

### 1.4 Fig 5C axis calibration — sound, and the class-20 failure did *not* occur

I located the major ticks myself as the long horizontal marks left of the spine (x 73-80,
8 px; minor ticks are only 3-4 px):

```
10^9 -> y = 1465.5     10^8 -> 1555.5     10^7 -> 1645.5     10^6 -> 1735.5
```

Exactly **90.00 px/decade**, zero residual, so `log10(v) = (2275.5 - y)/90`. This was then
checked against five *unlabelled minor* ticks, which the curation did not use:

| minor tick | predicted y | observed y |
|---|---|---|
| 8 x 10^8 | 1474.2 | 1474 |
| 7 x 10^8 | 1479.4 | 1479 |
| 6 x 10^8 | 1485.5 | 1485.5 |
| 5 x 10^8 | 1492.6 | 1492 |
| 4 x 10^8 | 1501.3 | 1501 |

The `ql_2026` failure mode (class 20 — a log axis anchored on the panel box) is **not**
present: the calibration is anchored on tick marks and is correct to better than 0.01
decades across the whole axis. The bar->name mapping is also confirmed — I read the rotated
x-axis labels (ALC-0315, SM-102, MC3, E1...E9) and the twelve axis tick centres match my
twelve bar centres to < 1 px.

### 1.5 Fig 5D and both printed numbers

Read by sampling the 12 x 6 cell grid (rows y 1464-1830, 30.58 px pitch; Reality columns
centred at x 644/724/804, Prediction at 920/1000/1080):

| | ALC | SM | MC3 | E1 | E2 | E3 | E4 | E5 | E6 | E7 | E8 | E9 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| Reality | A | A | A | A | A | B | B | B | C | C | C | C |
| Prediction | A | A | A | A | A | **C** | **C** | B | C | C | **B** | C |

- overall correct **9/12 = 75.0 %** — matches the printed 75 % (text p.5066 and the Fig 5E pie)
- class B correct **1/3** — matches "only one-third of the predictions being correct"
- additionally: class A 5/5, class C 3/4, and "the absence of misclassifications as class A"
  (p.5067) — no non-A row is predicted A.

Class thresholds independently confirmed from two places: section 2.1.1 ("class A ... exceeding
5 x 10^7, class B ... between 5 x 10^7 and 5 x 10^6, ... class C ... below 5 x 10^6") and the xlsx
Sheet1 header (`Three-part method A.(>5*10^7) B.(5*10^6~5*10^7) C.(<5*10^6)`).

### 1.6 All twelve rows, all seventy-two fields

No sample. Every field of every row, end to end.

- **Structures**: 12/12 parse, 12/12 canonical, 0 duplicate SMILES within the experiment,
  0 disconnected fragments, 0 radical electrons.
- **Descriptors**: all 14 RDKit columns + all 5 ML columns recomputed from each row's own
  SMILES — **0 mismatches** across 12 x 19 = 228 values. `has_ester/has_carbonate/
  has_disulfide` re-derived by SMARTS — 12/12 correct (esters only on the three benchmarks).
- **Formulation**: `50 / 10 / 38.5 / 1.5` summing to 100, DSPC / Cholesterol / DMG-PEG2000,
  `microfluidics`, `Aqueous_buffer = citrate`, `Dialysis_buffer = citrate` — all verbatim
  from Methods 2.3.3 ("0.01 mol/L citric acid buffer (pH 3)... 50 %, 10 %, 38.5 %, 1.5 %...
  microfluidic device, dialyzed using 0.01 mol/L citric acid buffer at pH ~ 6 for 2 h").
- **Model**: `Mouse_BALBc` from Methods 2.3.2 ("BALB/c mice (6-8 weeks)... Male animals");
  `Cargo = mRNA`, `Cargo_type = FLuc` from 2.3.1 ("mRNA encoding luciferase (mLUC)");
  6 h timepoint from 2.3.4.
- **Identifiers**: PMID 42592089 and DOI 10.1016/j.apsb.2025.10.040 present and
  consistent on all 12 rows; `Publication_link` resolves to the same article.
- **Pipeline fidelity**: `07_XD_2026_rows.csv` -> `data/data_vivo.csv` differs **only** in
  `Index`, `LNP_ID` and the 22 descriptor columns, all of which the emit file leaves blank
  for `build_dataset.py` to fill. Every curated field transfers byte-exactly.
- **z-score**: `Experiment_value` equals `z(log10(04_fig5c value))` with population sd to a
  max absolute difference of **0.0** — the stored value inverts exactly back to its source.
  (Mean 0 / sd 1 is noted as a tautology and is not offered as evidence.)
- **Correctly left blank**: `Dose_ug_nucleicacid` and `IL_to_nucleicacid_chargeratio`. I
  searched the whole paper: the only dose given is d-luciferin 150 mg/kg (the substrate),
  and N:P appears nowhere. Not imputed — right call.
- **Shared lipids**: SM-102 with SB_2024/SS_2018/SY_2025, DLin-MC3-DMA with
  BL_2023/SB_2024/ZC_2023/ZH_2023 — identical SMILES in every case, so they are genuine
  cross-source anchors, not name collisions.

### 1.7 C10 and C11 for this folder — both present and both pass

The packet text I was issued says in section 5 "No anchor records found for this folder" and in
section 6 lists `{}` for C11. **Both statements are false of the repository**, and the live packet
on disk (`round_2026-10-30/packet_XD_2026.md`, regenerated 01:46:31 after `brief.py` was
fixed at 01:46:16) already contains the correct JSON. The copy pasted into my prompt is a
stale snapshot; this is not a repo defect, and I flag it only so the coordinator knows the
briefing text was out of date. For the record:

- `provenance/structure_check.json` holds **14** anchors, 12 of them for `xd_2026_lipidai`:
  SM-102 and DLin-MC3-DMA `exact_match`, ALC-0315 `expected_mismatch_paper_drawing_is_wrong`
  (reference `WESYLBGVEJKNRX` = the 3-carbon drawn form vs enumerated `QGWBEETXHOVFQS`),
  E1-E9 `NO_DRAWN_STRUCTURE` with the absence spelled out.
- The anchor strings are written in a **different atom order** from the emitted ones
  (ALC-0315 and SM-102 head-first, MC3 amine-first), so the comparison is not a copy test.
  Residual caveat: they were still transcribed by the same session from the same figure, so
  for SM-102/MC3 they test transcription, not independent derivation. My name-derived
  InChIKeys (1.3) close that gap.
- `provenance/reagent_registry.json` (C11) covers ALC-0315, SM-102, DLin-MC3-DMA, DSPC and
  Cholesterol, all with `status: ok`, and `XD_2026` is listed among the experiments for each.

### 1.8 The LipidDB exclusion — endorsed, with every supporting number re-derived

I agree the 485-row LipidDB must not carry one `Experiment_ID`, and I verified each of the
curation's stated reasons from the primary files rather than from its summary:

| claim | my independent result |
|---|---|
| `descriptors_1107.csv` carries `inject_way` with four levels | 961 / 85 / 50 / 11, and that feature is absent from the released xlsx |
| 7 duplicate structures with conflicting class labels | reproduced, spreads up to 5 orders of magnitude |
| one paper supplies 124 of 485 rows | Lu/Percec *Pharmaceutics* 2023 = 124 rows |
| top 10 sources supply 318 rows | 124+40+31+30+25+18+14+13+12+11 = 318 |
| 83 lipids shared with our corpus | 83 **unique lipids** (= 89 rows: 81 vitro, 69 vivo) |
| Table 1 vs xlsx: 3 luminescence disagreements, 2 also class | row 4 printed 3.16e7/B vs xlsx 9.0e3/C; row 5 printed 5.8e7/A vs xlsx 6.62e6/B; row 15 printed 2.0e8 vs xlsx 2.0e10 (both A) |
| Table 1 row 4 cites a different source than the xlsx | Table 1 cites ref 27; xlsx row 14 attributes it to Zhang D/Percec |

Two additions of my own, both supporting the exclusion:

- **SM-102 is in the training set** (xlsx Sheet1 row 354, exact canonical match). The paper
  validates its model on a lipid it trained on. ALC-0315 and MC3 are not present.
- **Table 1 row 8's structure appears four times in the xlsx**, from three different source
  papers, with values 1.5e9, 1.5e9, 1.4e10 and 1.0e8 — a 140x spread under
  one label. All four happen to be class A, so this case does not appear in the "7
  conflicting class" list, but it is the same pathology.

### 1.9 "No external anchor for E1-E9" — I tried to break this and could not

Methods 2.3.1 says "**As shown in Table 1**, Lipid E 1-9 and mRNA encoding luciferase
(mLUC) were provided by WestGene Biopharmaceutical." That cross-reference looked like a
live lead, so I followed it: **Table 1 is the LipidDB sample table** (26 printed rows plus
"27 ..."), and contains no E lipid. I checked all six tables in the paper, found no
supplementary information is offered anywhere in the text, and searched the cloned repo
(`repo/`, 53 files) — the only data files are the LipidDB xlsx, `descriptors_1107.csv`,
`Model_Refinement_1.xlsx`, the per-model contrast workbooks and a pickled model. **None of
E1-E9 appears in the LipidDB under any canonical or stereo-blind form.** The curation's
"no name, formula, HRMS, synthesis or SI" claim survives. The authors' own "Table 1"
pointer is a citation error in the paper.

---

## 2. ERRORS FOUND

### E-1 — Fig 5C bar tops are detected *inside* the bars, biased low on every bar

**File / line.** `curation/xd_2026_lipidai/stage04_fig5c.py`:

```python
col=fill[:, x0+6:x1-6]
rows=np.nonzero(col.sum(1) > 0.6*(x1-x0-12))[0]
top=rows.min()
```

**What is wrong.** The rule declares the bar top at the first row where **more than 60 % of
the bar's width** is lavender fill. The three replicate scatter circles are drawn as white
discs with thin dark rims and in most bars they **straddle the bar's top edge**, occluding
the middle of the width. The detector therefore cannot fire until it has descended past the
circles. The bias is one-sided — always low, never high — and its size is whatever fraction
of the width the circles happen to cover, which is why it is *not* a constant offset and
does not cancel in the z-score.

**Evidence.**

1. *Direct pixel evidence.* For ALC-0315, the bar's antialiased top edge is a full-width
   non-white row at **y = 1476**, with solid fill from 1477 (L = bar fill, o = intermediate,
   . = white, # = ink):
   ```
   1475 .......ooooooooo...........   <- circles only; both flanks pure white
   1476 ooooooooo##o#oooooooooo..     <- full-width edge of the rectangle
   1477 LLLLoooooo#o#ooLLLLLLLLo.     <- fill, interrupted by the circles
   1483 LLLLLoooooLoooLLLLLLLLLL.     <- first row that passes the 60 % test
   ```
   The curation's reported value, 8.8066, back-solves to y = 1483.1 — **seven pixels inside
   the bar.**
2. *Independent estimator.* I measured each bar top from its **flank columns only**
   (x0+1...x0+3 and x0+20...x0+22), which the centred circle clusters never occlude. Left and
   right flanks agree to <= 1 px on **every one** of the twelve bars — a bar is a rectangle,
   so this is its top edge, unambiguously.
3. *Mechanism confirmed by reproduction.* Applying the curation's own 60 %-of-width rule to
   my native raster reproduces its implied tops to within 1 px on **10 of 12** bars. The
   discrepancy is the rule, not the raster or the calibration.
4. *Pattern.* The four bars where the curation agrees with me (delta <= 1.4 px: E2, E7, E8, E9)
   are exactly the four whose scatter circles sit clear above the bar top. The five worst
   (delta 6.9-10.9 px: ALC-0315, E1, E4, E5, E6) are exactly those whose circles straddle it.

**Correct values.** Calibration `log10(v) = (2275.5 - y)/90` on the native raster; `y` =
the bar's top antialiased row.

| lipid | curation value | curation log10 | **correct value** | **correct log10** | corpus z | **correct z** | dz |
|---|---|---|---|---|---|---|---|
| ALC-0315 | 6.406e+08 | 8.8066 | **7.644e+08** | **8.8833** | +1.2501 | +1.2523 | +0.0022 |
| SM-102 | 6.009e+08 | 8.7788 | **6.901e+08** | **8.8389** | +1.2265 | +1.2148 | -0.0116 |
| DLin-MC3-DMA | 1.328e+08 | 8.1231 | **1.487e+08** | **8.1722** | +0.6682 | +0.6527 | -0.0155 |
| E1 | 3.743e+08 | 8.5732 | **4.823e+08** | **8.6833** | +1.0514 | +1.0837 | +0.0322 |
| E2 | 2.517e+08 | 8.4009 | **2.610e+08** | **8.4167** | +0.9048 | +0.8588 | -0.0460 |
| E3 | 3.693e+07 | 7.5674 | **4.137e+07** | **7.6167** | +0.1951 | +0.1843 | -0.0109 |
| E4 | 1.153e+07 | 7.0617 | **1.525e+07** | **7.1833** | -0.2354 | -0.1811 | +0.0542 |
| E5 | 2.214e+07 | 7.3451 | **2.747e+07** | **7.4389** | +0.0059 | +0.0344 | +0.0284 |
| E6 | 6.153e+05 | 5.7891 | **8.046e+05** | **5.9056** | -1.3187 | -1.2585 | +0.0602 |
| E7 | 4.887e+05 | 5.6891 | **4.948e+05** | **5.6944** | -1.4039 | -1.4366 | -0.0327 |
| E8 | 1.123e+06 | 6.0503 | **1.151e+06** | **6.0611** | -1.0964 | -1.1274 | -0.0310 |
| E9 | 7.455e+05 | 5.8725 | **7.644e+05** | **5.8833** | -1.2478 | -1.2773 | -0.0295 |

Mean log10 shift **+0.060 decades** (range +0.005 to +0.122). Effect on the stored
`Experiment_value`: **max |dz| = 0.060, rms 0.034**, and **E6 and E9 swap rank** (they
differ by 2 px, so that swap is itself inside my uncertainty and should not be read as a
finding about those two lipids).

**How to verify this is the correct value.** Three routes, all cheap: (i) print the pixel
map around any bar top and confirm the full-width antialiased edge row sits above the
detected one — the ALC-0315 map above is reproducible in four lines; (ii) measure the top
from the flank columns only and check left and right agree; (iii) replace the `> 0.6 * w`
rule with `np.nonzero(col.any(1))[0].min()` restricted to the flanks, and confirm the four
unoccluded bars (E2, E7, E8, E9) are unchanged while the occluded five move up.

**Severity: changes a VALUE** (all 12 rows). **Confidence: certain** on the mechanism,
direction and sign; **probable** on the exact magnitudes, which carry +/- 1 px = +/- 0.011
decades.

### E-2 — `confidence = medium` never reached the corpus; three documents say it did

**Files.** `CLAUDE.md` section 9 (2026-10-29 entry); `provenance/CHANGELOG.md:399`;
`provenance/ERROR_CLASSES.md` rows "Commercial reagent that is an isomer mixture" and
"A structure with no external anchor at all".

**What is wrong.** `data/data_vivo.csv` and `data/data_vitro.csv` are both **72 columns with
an identical header, and there is no `confidence` column, no `positional_ambiguity` column,
no `*_assumed` column and no free-text note column.** I grepped both headers. The
`confidence = medium` flag for E1-E9 exists **only** in
`curation/xd_2026_lipidai/07_XD_2026_values.csv` and in the C10 note. CLAUDE.md section 9 states
"**E1-E9 carry `confidence = medium`**" inside a sentence whose subject is the merge into
`data_vivo.csv`, which a reader will take as a property of the corpus. It is not.

The same gap is corpus-wide, not specific to this paper: `ERROR_CLASSES.md` prescribes
"flag `confidence = medium` and say so" and "flag `positional_ambiguity` / `*_assumed`,
confidence `low`" as the standing mitigations for three separate error classes, and the
schema cannot hold any of them. `CHANGELOG.md:236` makes the same claim for 24 `C1-C4` rows
(`quaternization_site_assumed=1` and confidence `low`) and `:149` for 140 `Ac8-1/Ac10-1`
rows.

**Consequence.** A consumer reading only `data/` cannot tell that nine of these twelve
structures have **no external anchor of any kind** and were read off a figure with a
demonstrated one-carbon error on its own controls. Exhaustively, for this slice: *nothing*
qualitative reached the corpus — not `confidence`, not the ALC-0315 override, not the
`Model_target` / `Route` / `Dialysis_buffer` readings. Only the blank `Dose_ug_nucleicacid`
and blank `IL_to_nucleicacid_chargeratio` are visible, and those are absences, not flags.

**Correct state.** Either add the column to the schema and the builder and populate it, or
amend section 9 / CHANGELOG / ERROR_CLASSES to say plainly that the flag lives in the curation
folder and the corpus cannot carry it. The second is cheap and immediate; the first is the
real fix.

**Severity: METADATA only.** **Confidence: certain.**

### E-3 — the "299 vs 47" claim is unreproducible, and its two records contradict each other

**Files.** `CLAUDE.md` section 9: "the LipidDB publication count was 299, not 47".
`provenance/ERROR_CLASSES.md`, "De-duplicating reference strings on an assumed identifier
format": "cutting only at `doi:` gave 299 'publications' against a true 47". These read
in opposite directions: the first says 299 is the answer, the second says 299 was the bug.
My brief, derived from section 9, told me the true number *is* 299. It is not.

**What I did.** Pulled all 485 `References` strings from the xlsx and de-duplicated under
nine different rules:

| rule | distinct |
|---|---|
| raw strings, whitespace-normalised | **399** |
| strip a trailing `doi:` token | 286 |
| **cut at the first `doi:`** (the naive rule the registry blames) | **116** |
| cut at the last `doi:` | 116 |
| delete any `10.xxxx/...` substring | 226 |
| cut at the first identifier of **any** kind (`doi:` / `http` / `PMID` / `PMCID` / `YYYY;`) | **56** |
| ...then letters-only | **50** |
| everything before the first 4-digit year, letters-only | 49 |
| delete all digits | 50 |

**299 is not produced by any of them** (closest: 286). The naive `doi:`-only cut gives 116,
not 299.

**What is actually right.** Cutting at the first identifier of any kind gives **56**
distinct source strings, **50** after letters-only normalisation, against the paper's stated
**47**. The paper's own number is independently corroborated by its citation span — section 2.1.1
cites "24-40,41-70", and 70 - 24 + 1 = **47** exactly. So 47 is right, the residual 50 - 47
is three strings that survive normalisation (formatting variants), and the honest statement
is "47 stated, ~50 recoverable". **The curation folder already has this right**:
`06_lipiddb_qc.json` records `distinct_reference_strings_after_normalisation: 56`,
`..._letters_only: 50`, `distinct_source_publications_stated_by_paper: 47`, and never
mentions 299. The defect is entirely in the two documents.

**Correct state.** Delete "299" from `CLAUDE.md` section 9 and `ERROR_CLASSES.md`, or replace it
with the reproducible pair (116 under the naive `doi:` cut -> 50 under the correct cut -> 47
stated). Keep the error *class*, which is real and valuable; it is only the number that is
wrong.

**Severity: METADATA / documentation.** **Confidence: certain** that 299 is unreproducible
from the released file; **probable** as to which of the two contradictory readings was
intended.

### E-4 — a printed summary statistic was not reproduced, and it resolves an "Unexplained"

**File.** `curation/xd_2026_lipidai/README.md`, "Known defects *in the published dataset*":
"The paper states 475 lipid samples; the xlsx has 485 rows. **Unexplained** ... 485 + 720 =
1205 != the stated 1107 augmented points either."

**What is missing.** The paper prints a *third* count that the curation never reproduced.
p.5067 (section 3.4): "After data amplification, the dataset comprised **472 class A samples, 447
class B samples, and 287 class C samples**." Counting the three-part class column across
all three xlsx sheets:

| sheet | rows | A | B | C |
|---|---|---|---|---|
| Sheet1 | 485 | 201 | 171 | 113 |
| sheet1+tail | 360 | 135 | 138 | 87 |
| sheet1-tail | 360 | 135 | 138 | 87 |
| **total** | **1205** | **471** | **447** | **287** |
| paper | 1206 | 472 | 447 | 287 |

**B and C match exactly; A is one short.** This is close to decisive: the dataset the paper
describes as "after data amplification" **is** the 1205-row xlsx, and section 2.1.2's "1107 valid
data points" is the outlier — confirmed by `descriptors_1107.csv`, which has 1107 rows and a
*different* class split (A 437 / B 402 / C 268). So the README's framing is half right: the
1107 mismatch is real, but 1205 is corroborated by the paper itself and should not sit in an
"Unexplained" bucket beside it. This is exactly the rule the review program exists to
enforce — reproduce **every** printed summary statistic.

**Severity: METADATA only** (LipidDB is not merged; nothing in `data/` moves).
**Confidence: certain.**

### E-5 — `04_fig5c.json` cannot be re-derived from its own stored numbers

**File.** `curation/xd_2026_lipidai/04_fig5c.json`.

**What is wrong.** `calibration.tick_rows` = `[81.0, 261.5, 441.0, 621.0]` are **crop-relative**
(the script's `Y0 = 4850` has been subtracted), while every `bars[].top_y` = `[4966 ... 5527]`
is in **full-page** coordinates (`top + Y0`). Applying the stored `slope`/`intercept` to the
stored `top_y` yields `log10 ~= -18.1`, not 8.81. The two coordinate frames in one artifact
cannot be composed. The values themselves were computed consistently inside the script, so
no number is wrong — but the artifact is not independently checkable, which is the whole
point of storing it. (This is a cousin of class 21: the stored provenance does not invert
back to its source.)

**Correct state.** Store both in the same frame, and store `Y0`/`X0` explicitly.

**Severity: COSMETIC** (no value changes). **Confidence: certain.**

### E-6 — "600 dpi" overstates the evidence; the figure is 300 dpi at source

**Files.** `README.md`, `stage04_fig5c.py` docstring and the `p13_600.png` path;
`CLAUDE.md` section 5 ("Render with `pymupdf` at 600 dpi").

**What is wrong.** Page 13's figure is a single embedded raster of **1541 x 1864 px** placed
in a **369.79 x 447.26 pt** box — that is **300.0 dpi in both axes**. A 600-dpi page render
is a 2x interpolation and contains no information the native raster does not. The curation
(and the first Tier-2 reviewer, which worked "at 300 dpi native") reached the right answers
anyway, and the tick calibration survives the resampling (the 600-dpi tick spacings
180.5/179.5/180.0 are the interpolation of a true 180.0). But the claim invites the reader
to believe the measurement has twice the precision it has, and the general advice in section 5 is
worth qualifying: *check the embedded raster's native dpi first; rendering above it buys
nothing.*

**Severity: COSMETIC.** **Confidence: certain.**

### E-7 — the "Fig 5C <-> Fig 5D 12/12" cross-check has no discriminating power

**File.** `README.md`: "**All 12 digitised values fall in their own Fig 5D band** — two
figures agreeing is what makes the digitisation trustworthy, not the fit alone."

**What is wrong.** The claim is true but carries almost no information, and it was relied on
as the validation that the digitisation is right. The A/B/C bands are a **full decade wide**;
the digitisation error found in E-1 is **0.06 decades**. My values — which differ from the
curation's by up to 0.122 decades — also land **12/12**. A test that both the right and the
wrong answer pass is not evidence. The nearest value to a boundary is E3, which sits 0.077
decades below the A/B threshold under my reading and 0.13 under the curation's; even a 3x
larger error would not have moved it.

This matters beyond this paper: the same sentence pattern ("two figures agree, therefore the
digitisation is trustworthy") will recur, and it should be stated with the band width next
to the error, so the reader can see whether the check could have failed.

**Severity: METADATA / documentation.** **Confidence: certain.**

### Minor notes (not counted as findings)

- `06_lipiddb_qc.json`'s `overlap_with_corpus_vitro / _vivo / _union` (81 / 69 / 89) are
  **row** counts, while `README.md`'s "83 lipids shared with our corpus" is the **unique
  lipid** count. Both are correct — I reproduced 81/69/89 and 75/67/83 — but the key names
  do not say "rows".
- `07_qc.json`'s `assumptions.Route_of_administration` describes `intravenous` as a reading
  because Methods 2.3.4 says only "post-vaccination". In fact section 3.4 of the body states
  "6 h after tail vein injection" and the Fig 5 caption says "intravenously administered"
  and "intravenous injection" — three printed statements that agree. The flag is honest but
  over-cautious; this field is as well attested as any other in the row.

---

## 3. UNVERIFIABLE

These are genuinely unknowable from the paper, and a flag is the correct deliverable. They
are **not** padding for the findings list.

1. **The identity of E1-E9.** The paper gives no name, formula, HRMS, NMR, synthesis or
   supplementary information; the only "Table 1" pointer is a citation error; the cloned
   repository contains none of the nine. Fig 5A is the sole source, and that panel has a
   verified +/-1-carbon error on one of its three controls. *Settled by:* the WestGene
   Biopharmaceutical product specification, or the authors' characterisation data.
   The current posture (emit from the trace, flag `confidence = medium`, store an empty
   `reference_SMILES` row per lipid) is right — subject to E-2, which says the flag does not
   actually reach the corpus.
2. **Whether the ALC-0315 head error is a drawing error or a real 3-carbon analogue.**
   Methods 2.3.1 says the benchmarks were "synthesized in our laboratory" and gives no
   characterisation at all. The evidence for a drawing error is strong (E1 drawn with four in
   the same style; both other benchmarks drawn exactly; a defined compound name used as a
   literature benchmark), and I endorse the override — but it is an inference, not a fact.
   *Settled by:* the authors' HRMS/NMR for their ALC-0315 batch.
3. **`Model_target = whole_body`.** Methods 2.3.4 describes *both* whole-animal imaging at
   6 h *and* ex-vivo heart/liver/spleen/lung/kidney on the same instrument, and neither the
   Fig 5C caption ("Luciferase expression histogram") nor the body text says which ROI the
   histogram plots. The Fig 5C y-axis unit (p/sec/cm^2/sr) is identical for both.
   `whole_body` is the defensible reading for one number per lipid, and 22 corpus rows
   already use that category. *Settled by:* the source data or the authors.
4. **mRNA dose and N:P ratio.** Absent from the paper. Correctly blank; the only dose
   printed anywhere is d-luciferin 150 mg/kg.
5. **n per group in Fig 5C.** Three scatter points are drawn per bar but the paper never
   states n, and never says whether they are individual animals or technical replicates.
   `Experiment_batching = individual` is a reasonable reading; it is not stated.
6. **Whether the bar heights are arithmetic or geometric means of the three points.** I
   attempted to settle this by detecting the replicate discs and comparing, but at 300 dpi
   the discs are ~7 px across with faint antialiased rims and overlap each other in most
   bars; centroid precision is +/-2-3 px, which is the same order as the effect. Not
   determinable from this raster. It does not affect E-1, which is a property of the bar's
   own top edge.
7. **The one-row discrepancy in E-4** (class A 471 vs 472). Could be the authors' arithmetic,
   a row dropped between versions, or the five blank rows in `sheet1+tail`. Not determinable.

---

## 4. SEVERITY

| # | finding | severity |
|---|---|---|
| E-1 | Fig 5C bar tops detected inside the bars (60 %-width rule vs. overlying scatter discs) | **changes a VALUE** — all 12 rows; max abs(dz) 0.060, rms 0.034 |
| E-2 | `confidence = medium` asserted of the corpus; no such column exists in either 72-column file | **METADATA only** |
| E-3 | "299 vs 47" unreproducible; `CLAUDE.md` section 9 and `ERROR_CLASSES.md` contradict each other | **METADATA only** |
| E-4 | the paper's printed class counts (472/447/287) never reproduced; they resolve an "Unexplained" | **METADATA only** |
| E-5 | `04_fig5c.json` mixes crop-relative and page coordinates; cannot be re-derived from itself | **cosmetic** |
| E-6 | "600 dpi" asserted of a figure that is 300 dpi native | **cosmetic** |
| E-7 | the 12/12 Fig 5C<->5D agreement is offered as validation but cannot discriminate | **METADATA only** |

**No structure changes.** All 12 emitted structures are correct as drawn-and-corrected, and
I verified each one by an independent route.

## 5. CONFIDENCE

| # | confidence | basis |
|---|---|---|
| E-1 | **certain** on mechanism, direction and sign; **probable** on magnitudes (+/- 1 px = +/- 0.011 decades) | pixel maps show the full-width edge row above the detected row on every biased bar; left and right flank estimates agree to <= 1 px on all 12 bars; applying the curation's own rule to my raster reproduces its tops to +/- 1 px on 10/12 |
| E-2 | **certain** | both corpus headers read and compared; no such column exists |
| E-3 | **certain** that 299 is unreproducible (9 de-duplication rules tried, none yields it); **probable** as to which reading was intended | the curation's own QC file records 56/50/47 and never 299 |
| E-4 | **certain** | class counts recomputed from all three xlsx sheets; B and C match the printed numbers exactly |
| E-5 | **certain** | applying the stored calibration to the stored `top_y` gives log10 ~= -18 |
| E-6 | **certain** | embedded image 1541 x 1864 px in a 369.79 x 447.26 pt box |
| E-7 | **certain** | my own values, differing by up to 0.122 decades, also land 12/12 |

---

### Assessment of the previous Tier-2 pass

Its structural work holds up completely. The ALC-0315 finding is real and I reproduce it by
an independent method; the E1-E9 traces are correct including the N-adjacent near-vertical
vertex, which is exactly where a vertex counter would be expected to fail; the benchmark
SMILES are the right compounds; the metadata is faithful to the paper; and the decision not
to merge LipidDB is sound, with every supporting statistic reproducible. The claim that none
of its nine fixes touched an emitted structure or value is true as far as it goes.

Where it stopped short is the part no structural check covers: it validated the Fig 5C
digitisation against a cross-check (the A/B/C bands) too coarse to fail, and so did not
examine the bar-top estimator itself; it did not reproduce the paper's printed class counts;
and the one number it is most celebrated for — 299 — does not come out of the released file
under any rule I could construct, while the curation folder it produced records the correct
figures throughout. The lesson is the one already in `ERROR_CLASSES.md`, applied one level
up: a cross-check is only worth its width, and a number recorded in a changelog needs the
same re-derivation as a number recorded in the corpus.
