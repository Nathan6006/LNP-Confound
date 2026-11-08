# Round-2 Wave-1 review — `LL_2026` / `LL_2026_2` (Lin et al., JACS 2026, 148, 6960−6969, PMID 41689530)

**Provenance.** Returned as text by a read-only reviewer and saved verbatim.
**Its findings are claims.** My adjudication is in `FINDINGS.md`. This is the
report that settles F-008, F-009 and the `10A2-3B-EP10` anchor.

---

Read-only. Nothing in the repo was written. Python run as `python3 -I`; pages rendered at 600 dpi with pymupdf and looked at with the Read tool; overlays built and inspected.

**Starting point:** my own 600-dpi render of PDF page index 2 is **byte-identical** (max abs channel diff 0) to `curation/lin_aromatic_2026/00_pages/p2_600.png`, so every pixel measurement below is on the same raster the curation used, independently produced.

---

## 0. THE THREE NAMED ITEMS — VERDICTS FIRST

### F-008 — **REAL, exactly as claimed. Currently harmless. Fix it.**

`digitize_fig2a.py:42` sets `BODY = (640, 679, 3161, 1073)`.

Measured from the raster, three independent ways:

| quantity | script | measured |
|---|---|---|
| left edge | 640 | **641.1** (full colour from x=644, transition 640–643) |
| right edge | **3161** | **3103.6** (full colour to 3102, transition 3103–3106, white from 3107) |
| column pitch | 120.048 | **117.26** |
| top / bottom | 679 / 1073 | 677 / 1075 (pink header band 617–674, white gutter 675–676) |

Right-edge overrun = **3161 − 3103.6 = 57.4 px**. The round-1 reviewer's "57 px" and "pitch 117.33" are both correct.

I recovered the 22 column boundaries by x-gradient summed over the 7 row-interiors: 641.1, 755.2, 874.4, 991.8, 1107.7, 1226.8, 1342.9, 1460.1, 1578.4, 1696.2, 1812.9, 1930.7, 2047.5, 2165.0, 2282.9, 2399.6, 2518.7, 2636.6, 2752.7, 2871.1, 2988.1, 3103.6 — uniform to ±1.5 px.

**(b) Exactly which columns it mis-samples.** Window centre drift = `0.35 + 2.79·c` px; sampling half-width 30 px; true cell half-width 58.6 px.

| col | label | fraction of window in its own cell | in the NEXT cell | off the heatmap (white) |
|---|---|---|---|---|
| 0–10 | 1A1…10A1 | 1.00 | 0 | 0 |
| 11 | 10A2 | 0.95 | 0.05 | 0 |
| 13 | 11A2 | 0.86 | 0.14 | 0 |
| 15 | 12A1(−2) | 0.77 | 0.23 | 0 |
| 17 | — | 0.68 | 0.32 | 0 |
| 19 | 12A1 | 0.59 | **0.41** | 0 |
| 20 | **12A2** | 0.54 | 0 | **0.46** |

So columns **11–19** sample up to 41 % of the *neighbouring* lipid's cell, and column **20 (12A2)** samples 46 % pure white off the right edge. Measured directly: max white fraction in any sampled window = **0.367** under the script's geometry, **0.000** under mine.

**(c) Do any bin assignments change? No — 0 of 147.** I re-digitized Fig 2a under both geometries with the same legend swatches: both give the identical 147 bins, and both match `04b_fig2a_147.csv` exactly (histogram 19 / 41 / 87). The per-pixel median survives because contamination never exceeds 50 %; **column 20 clears it by 4 percentage points.**

The diagnostic that *should* have caught it: `max colour distance` is **2.29** under the script geometry and **0.00** under the measured one. CLAUDE.md §5's ">~5 means geometry, not colours" gate is calibrated for catastrophic offsets and is blind to a drift that eats 46 % of the sampling window.

**One thing round 1 did not notice: the defect is visible in the curation's own `review/fig2a_overlay.png`.** Cropping its right-hand end shows the green body box running past 12A2 into white and the red sampling boxes in the last columns sitting visibly off-centre, straddling the cell boundary. The protocol's overlay was generated and saved; the right-hand end was not looked at.

**Proposed correction:** `BODY = (641, 677, 3104, 1075)`. (Round 1's suggested `(640, 675, 3104, 1073)` is also fine; the difference is sub-pixel.) *Verify by re-running the digitization and requiring `max_distance == 0.00` and `max white fraction == 0.000`, not just "bins unchanged".*

**Severity:** cosmetic today, latent off-by-N. One more column, or a 2.5 % larger pitch error, flips column 20 to white → the lightest bin.
**Confidence:** certain.

**Checked and clean by contrast: Fig 2e.** Measured body x ∈ [3385.5, 4624.5], y ∈ [734, 4141.5] vs the script's `(3385, 733, 4627, 4142)` — within 2.5 px everywhere; drift 0.3 px/column against a 31 px half-window. Re-digitizing under my geometry and under the script's gives **identical 840 bins**, both matching `04_fig2e_840.csv` 840/840. Max colour distance **2.24**, minimum runner-up margin **49.73** (palette separation 50.2) — no cell is anywhere near a bin boundary.

---

### F-009 — **All three sub-questions resolved. The evidence supports a FLAG, not a structure change, and nothing stronger.**

**(a) Does the SI specify (S) for 5B, and is any other block specified?**
Yes, and no. SI p.2, Materials, verbatim: *"2,3-epoxypropyl-4-methoxyphenyl ether (4B), **(S)-glycidyl trityl ether (5B)**, benzyl acrylate (6B)…"*. A case-insensitive regex over the whole SI and the whole main text returns **exactly one** occurrence of `(S)` and **zero** of `(R)-`. All 29 other building blocks — including every EP epoxide and 1B/2B/3B/4B — are named without a stereodescriptor. 5B is the only one.

**(b) Is the regiochemistry such that C2 survives?**
Yes, and the paper demonstrates it rather than my assuming it. The general method (SI p.3) is uncatalysed aminolysis of a monosubstituted epoxide in IPA at 75 °C, which goes at the less-hindered terminal CH₂. The paper's own ¹H NMR of 10A2-3B-EP10 (Fig S21) pins it for the same glycidyl-ether class: *δ 4.11–3.80 (m, 6H)* = 2 × ArO**CH₂** (4H) + 2 × **CH**OH (2H), with all 24 N-CH₂ protons at 2.76–2.23. Attack at C2 instead would give ArOCH₂–**CH**(N)–**CH₂**OH and an entirely different pattern (a 1H N-CH near 2.8 and a 2H CH₂OH near 3.5). So the terminal-CH₂ regiochemistry is observed, and no bond to C2 is broken — its configuration is retained.

**(c) Any stereochemical characterisation of the products?**
**None.** Zero hits across main text and SI for optical rotation, `[α]`, polarimetry, circular dichroism, chiral HPLC, Chiralpak, "enantio-", "ee", "diastereo-". The only product characterisation is ¹H NMR plus HRMS on four compounds, none of which is a 5B lipid. No enantiomeric excess is stated for the reagent either.

**Scope in the corpus.** 141 rows (120 `LL_2026` = 5 amines × 2 equiv × 12 tails; 21 `LL_2026_2` = the 5B row of Fig 2a × 21 columns) — matches the packet. Those rows carry **210** 5B carbinols in total (one per arm; nA2 lipids carry two). Across the whole 987-row slice RDKit finds **3,955 unassigned tetrahedral stereocentres and 0 specified**; all but those 210 are genuinely racemic (EP tails, 1B–4B, and the 7B methacrylate α-carbon).

**My verdict: carry it as a flag. Do not rewrite the SMILES.** Reasons, in order of weight: (i) no product-level stereochemical data exists, (ii) no ee is stated for the reagent, (iii) the library was isolated as crude product with no chromatography, (iv) every lipid already carries 4–6 racemic carbinols, so asserting one centre converts an unspecified structure into a *specific diastereomer the paper never claims* — a strictly worse misstatement than leaving it blank. Round 1's reviewer reached the same conclusion and was right.

**Exact flag text.** Set `structure_assumed = 1` on the 141 rows and append to `curation_note`:

> `Aromatic module 5B is supplied as (S)-glycidyl trityl ether — the only building block this paper gives a stereodescriptor — and amine opening at the terminal CH2 leaves that carbinol untouched, so each 5B arm is in principle configurationally defined. Stored without stereochemistry: the paper reports no optical rotation, chiral HPLC or other stereochemical characterisation of any product, states no enantiomeric excess for the reagent, and the remaining carbinols are racemic, so the lipid is a diastereomer mixture regardless.`

If a `confidence` change is wanted too, `high → medium` on those 141 rows is defensible and consistent with how Ac8-1/Ac10-1 are treated; I would not go to `low`, because the connectivity is not in doubt.

**One trap for whoever ever does encode it.** I checked the CIP outcome with `rdCIPLabeler` on correctly-connected trityl SMILES written with an identical neighbour ordering at the stereocentre: the descriptor does **not** flip on opening — (S)-reagent gives an (S)-product. In the corpus's own SMILES ordering the (S) arm is

```
C[C@H](O)COC(c1ccccc1)(c1ccccc1)c1ccccc1       # (S)  — verified by rdCIPLabeler
C[C@@H](O)COC(c1ccccc1)(c1ccccc1)c1ccccc1      # (R)  — wrong enantiomer
```

(My first hand analysis of this got it backwards by following the epoxide ring closure to the real atom instead of to a CIP duplicate. Quote the computed token, not a derivation.)

**Severity:** structure (stereochemistry only; connectivity unaffected). **Confidence:** certain on (a) and (c), certain on (b) given the NMR.

---

### The C10 anchor `10A2-3B-EP10` — **derived from scratch; it matches. The anchor is correct.**

Derived entirely from Fig S17 and the SI synthesis text, by plain string assembly — no `RWMol`, no `AddBond`, no `RemoveAtom`, and `enumerate_library.py` never read or imported.

**My SMILES (as written, before canonicalisation):**
```
N(CC(O)COc1ccccc1)(CC(O)CCCCCCCC)CCCN(CC(O)CCCCCCCC)CCCCN(CC(O)CCCCCCCC)CCCN(CC(O)CCCCCCCC)CC(O)COc1ccccc1
```
**Canonical:**
```
CCCCCCCCC(O)CN(CCCCN(CCCN(CC(O)CCCCCCCC)CC(O)COc1ccccc1)CC(O)CCCCCCCC)CCCN(CC(O)CCCCCCCC)CC(O)COc1ccccc1
```
- Formula **C68H126N4O8**, ExactMolWt 1126.9576, **[M+H]⁺ = 1127.9648**
- **InChIKey `JUAPTPKSYXADRC-UHFFFAOYSA-N`** — skeleton `JUAPTPKSYXADRC`, stereo layer `UHFFFAOYSA` (no centres asserted)
- **Identical to the stored `reference_structures.csv` anchor and to the corpus row** (canonical-SMILES equality, not just InChIKey skeleton).

**Evidence, feature by feature:**

| feature | source | evidence |
|---|---|---|
| backbone = spermine, N-(CH₂)₃-N-(CH₂)₄-N-(CH₂)₃-N | Fig S17 (SI p.28), product drawing | vertex spacing measured at 33.0 px; N-to-N distances 130 / 160 / 126 px → 4 / 5 / 4 bonds → **3 / 4 / 3** carbons. Also SI Materials: "spermine (10A)" |
| 2 aromatic arms, **one on each terminal (primary) N** | Fig S17 drawing; synthesis route; Fig S3 caption | route is self-proving: **P1** MS [M+H]⁺ calc 395.19 = C14H24F6N4O2 = spermine + **2** CF₃CO (primary-selective) · **P3** 403.33 = spermine + **2** Boc (secondaries) with free primaries · **P5** 503.36 = spermine + **2** × 3B (C9H10O2) and nothing else. Each mass I recomputed from my structures: 395.188 / 403.328 / 503.359 |
| aromatic arm = PhO-CH₂-CH(OH)-CH₂-N | Fig S17, zoomed at 600 dpi | unsubstituted phenyl, ether O, 3-carbon chain, OH on the middle carbon; blocks.csv "glycidyl phenyl ether (3B)" matches SI Materials verbatim |
| 4 EP10 tails, N-CH₂-CH(OH)-C₈H₁₇ | Fig S17 | geometric vertex trace of one tail: 9 turning points above the N-CH₂ at 33 px spacing → C2 (carbinol) … C10 (CH₃) → **10 carbons**, i.e. 1,2-epoxydecane opened at the terminal CH₂ |
| all four N tertiary, no N-H left | Fig S21 ¹H NMR | *δ 2.76–2.23 (m, **24H**)* = 4 N × 3 CH₂ = 12 N-CH₂. Full integral accounting: 4+6+6+4+24+16+48+12 = **120 of 126 H**, the 6 missing being the exchangeable OH |
| no stereocentres asserted | Fig S17; SI Materials | no wedges anywhere in the figure; 3B and EP10 both named without a descriptor |
| overall | SI p.32 | **HRMS calc 1127.9648, found 1127.9687; mine 1127.9648 (+0.04 ppm)** |

**So the anchor is right — it simply cannot be *proved* right by C10 as stored.** Write the generator from the recipe above; the one-line route is: `ARM = "CC(O)CO" + "c1ccccc1"`, `TAIL = "CC(O)" + "C"*8`, then `"N(" + ARM + ")(" + TAIL + ")CCCN(" + TAIL + ")CCCCN(" + TAIL + ")CCCN(" + TAIL + ")" + ARM`. Pure string concatenation, no RDKit graph surgery — maximally disjoint from the enumerator. **Caveat for the generator:** my own string route tripped the ring-closure collision hazard twice during this review (the `6A` benzene and the `7A` piperazine both reusing ring number `1` with the arm's phenyl). Use a high ring index in any template that has its own ring, and assert the parsed formula.

---

## 1. CONFIRMED CORRECT (re-derived independently, with numbers)

1. **Raster provenance.** My 600-dpi render of PDF page index 2 == `00_pages/p2_600.png`, max channel diff **0**.
2. **Fig 2e geometry and bins.** 840/840 bins reproduced under two geometries; histogram **45 / 134 / 321 / 340**; max colour distance **2.24**; min runner-up margin **49.73**. The QC file's four stored legend RGBs `(233,238,244) (194,212,226) (123,174,215) (69,142,204)` reproduce exactly.
3. **Fig 2e labels, read at 600 dpi.** Column headers `8A1 8A2 9A1 9A2 10A1 10A2 11A1 11A2 12A1 12A2` confirmed by eye; 84 row labels `1B-EP6 … 7B-Ac10-1` confirmed, matching `ROW_LABELS`' aromatic-major / tail-minor order. Legend layout matches the order of `LEGEND[]`.
4. **Fig 2a bin indexing 2–4 is confirmed by colour identity, not just by the legend text.** The three Fig 2a swatches are Fig 2e's bins 2, 3, 4 at RGB distance **0.0 / 2.0 / 2.0**, against runner-up distances 50–65.
5. **The 48 Fig 2f/2g/2h bars — fully re-derived, tick-anchored, and they hold.** Scale from a least-squares line through tick centres only: panels f and g share an axis and give **5.12 px/%** (5 ticks, max residual 0.8 px); panel h gives **5.346 px/%** (6 ticks) and its fitted 0-tick lands at y = 4057.5, which is exactly the independently detected x-axis — a free consistency check. Bar tops measured from the **left and right flanks** (error class 32), required to agree within 2 px, every bar passing. Result: **48/48 within 0.65 pp** of my own Fig 2e marginals (panel f max 0.57, g max 0.64, h max 0.50). **The F-007 bar is present and exact: EP6 `>1e6` = 1.40 % vs 1.4 %.** The F-006 fix is genuinely implemented (`rep_l` appears nowhere in the scale; tick spacing is least-squares fitted and asserted).
6. **Printed statistics reproduced — all of them, none excused.** "1B−4B achieved an overall hit rate (RLU > 100,000) of over 85 %" → 87.5 / 90.8 / 87.5 / 89.2, pooled **88.75 %** ✓. "single aromatic ring … higher than triple aromatic rings" → trityl 5B (61.7) and anthracenyl 7B (58.3) are the two lowest of seven ✓. "chains of moderate lengths (C8−C14) yielded optimal delivery" → EP8 98.6, EP10 90.0, EP12 81.4, EP14 84.3 vs EP6 30.0, EP16 78.6, EP18 70.0 ✓. "a moderate increase in hydrophobic tail numbers slightly enhanced delivery" → 8A 69.6, 9A 69.6, 10A 76.8, 11A 83.9, 12A 93.5 against 5/6/6/6/7 N-H sites ✓.
7. **30/30 reagent names in `blocks.csv` match the SI Materials list verbatim**, including the three easy to get wrong: 7B is a **meth**acrylate, 2A is `NCCCNCCCN(C)C`, and Ac8-1/Ac10-1 are correctly left structureless.
8. **Twenty structures rebuilt from scratch by my own string-template route, 20/20 match the corpus.**
9. **All four HRMS anchors**, each computed from *my* structure: 10A2-3B-EP10 1127.9648 (+0.04 ppm), 9A2-3B-EP10 1113.9492 (−0.01), 9A2-4B-EP10 1173.9703 (+0.02), 10A2-3B-EP12 1240.0900 (+0.04). *HRMS is blind to connectivity — the regiochemistry below is what actually carries these.*
10. **Regiochemistry is printed, per-amine, and the enumerator matches it.** Fig S1 caption, verbatim: *"For 1B-5B bearing epoxy groups, the amines were first reacted with aromatic analogs (1 or 2 eq), followed by coupling of the remaining sites with alkyl chains. For 6B and 7B with acrylate groups, alkyl chains were first introduced, followed by coupling with the aromatic analogs."* Fig S3 caption, verbatim: *"**The primary amines of 8A-12A first reacted with the representative aromatic analog 3B (1 or 2 eq).** The remaining reactive sites were then coupled with hydrophobic alkyl chains."* Figures S2 and S3 then *draw* nA1 and nA2 products for **every one of 1A–12A**. I checked each drawing at 600 dpi; the enumerator's "primary before secondary, one per N" rule reproduces all of them.
11. **Fig 2a's column set is predicted by the paper's own rule.** 1A, 2A and 3A appear only as `A1`; those are exactly the three precursors for which two arms on two distinct primary nitrogens is impossible. Corroborated by the main text.
12. **z-scores exact.** Re-derived from the bin ordinals with population sd: max |z_recomputed − z_stored| = **0.00e+00** for both experiments. `LL_2026` 4 distinct values (equal spacing 1.1475); `LL_2026_2` 3 (equal spacing 1.4043).
13. **Metadata, all traced to source.** Fig 2e formulation 26.3/26.3/43.9/3.5 at IL:mRNA 10:1 ✓; Fig 2a formulation 24.2/32.3/40.3/3.2 at 20:1 ✓; dose 25 ng/well → 0.025 µg ✓; IGROV1, Fluc mRNA, citrate pH 4.4, 1:3 rapid mixing ✓.

---

## 2. ERRORS FOUND

### E1 — F-008, confirmed and reproduced. *(see §0)*
**File:** `digitize_fig2a.py:42`. **Correct value:** `BODY = (641, 677, 3104, 1075)`. **Severity:** cosmetic now, latent off-by-N. **Confidence:** certain.

### E2 — F-009, confirmed on all three sub-questions; deliverable is a flag. *(see §0)*
**File:** 141 rows in `data/data_vitro.csv`; flag emitted via `build_rows.py:quality()` and `enumerate_147.py`. **Severity:** structure (stereochemistry only). **Confidence:** certain.

### E3 — `positional_ambiguity = 1` on all 987 rows, and the `curation_note` on them is wrong on two counts
**New this round.**

The note reads, on every one of the 987 rows:
> *"One-pot library, no chromatographic separation of regioisomers: which **ring nitrogen** bears which arm is **a curation convention**, not a measurement."*

Both emphasised clauses are false.

(a) **"a curation convention"** — it is the paper's own assignment. Fig S3's caption states it in words and Figures S1/S2/S3 draw the product for every amine. I enumerated every way of distributing the arms over the N-H sites and canonicalised each product: **under the paper's stated rule the product is unique for every (amine, equiv) pair in the library, all 24 of them.** Stronger still, for `1A1`, `4A1`, `5A1`, `6A1`, `7A1` and `11A1` the product is unique even if you *ignore* the rule and allow arms anywhere — **126 corpus rows** are flagged ambiguous where no ambiguity exists on any reading.

(b) **"ring nitrogen"** — only `3A` (N-aminoethylpiperazine) has a reactive ring nitrogen, and it contributes 7 rows. The phrase is wrong on **980 of 987** rows, and this note is the user-facing text a downstream reader sees in `data_vitro.csv`.

What *is* defensible is a flag saying the real material is a crude one-pot product — the SI says the final products were obtained by evaporating solvent with no purification — so the actual sample is a mixture of regioisomers and incompletely-substituted species around the drawn structure. That is a different claim from "we chose".

**Correct value:** either (i) keep `positional_ambiguity = 1` and reword to *"One-pot library, crude product with no chromatographic purification: the stored structure is the regiochemistry the paper draws (Figs S1–S3 — aromatic arms on the primary amines, one per nitrogen), and the real material is a mixture around it"*; or (ii) set `positional_ambiguity = 0` on the 126 rows where the product is unique under any assignment and reword the rest. I prefer (i), with (ii) folded in.
**Severity:** metadata only. **Confidence:** certain on the two factual errors in the note; the flag *value* is a judgement call.

### E4 — the two sibling experiments label the same molecules with two different vocabularies
`build_rows.py:191-193` writes `amine_name` / `aromatic_name` / `tail_name`; `enumerate_147.py:98-100` writes `amine` / `aromatic` / the literal `"EP10"`. Result in the live corpus:

| | `IL_head_name` | `IL_tail1_name` | `IL_tail2_name` |
|---|---|---|---|
| `LL_2026` | `Spermine`, `Norspermine`, … | `Glycidyl phenyl ether`, … | `1,2-epoxydecane`, … |
| `LL_2026_2` | `10A`, `9A`, … | `3B`, … | `EP10` |

Same paper, same chemistry, same reagents — two naming systems, so any grouping or join on `IL_head_name` splits each component in two. This is the column set F-010's repair existed to populate. Secondary: the aromatic module goes into `IL_tail1_name` while `blocks.csv` classifies 1B–7B as `block_class = linker`, and `IL_linker_name` is blank on all 987 rows.
**Severity:** metadata only. **Confidence:** certain.

### E5 — `enumerate_147.py:111` hard-codes `confidence: "high"` (the class-29 shape)
`row.update(BR.quality({**r, "confidence": "high"}))`.

**It is currently correct, and the comment above it is honest about why.** `04b_fig2a_147.csv` has no confidence column; the 840-row rule assigns `medium` only for tails `Ac8-1`/`Ac10-1`; `EL.build(..., "EP10", ...)` hard-wires EP10 for all 147, so `high` is exactly what the rule yields. Verified: `LL_2026_2` is 147/147 `high`, 0 `structure_assumed`.

But the construct is class 29 in form: a literal at a call site standing in for a per-item assessment, correct only by a property of the current input. **This is the line that will make F-009's 21 `LL_2026_2` 5B rows carry `high` after you add the 5B flag**, because the literal overwrites whatever a per-item record would say.
**Correct value:** derive it, e.g. `"confidence": "medium" if r["tail"] in EL.SUPPLIED else "high"`, or have `enumerate_library.build()` return the confidence alongside the SMILES.
**Severity:** cosmetic today; **changes structures' metadata** once the F-009 flag lands. **Confidence:** certain.

### E6 — error class 30 (`comp.setdefault` by name) is still live in both of this folder's emitters
`build_rows.py:149-153` and `enumerate_147.py:61-65` both do `comp.setdefault(r[f"{pre}_name"], r[f"{pre}_SMILES"])` over the live corpus — first row wins, over a file in alphabetical experiment order. Round 1 hardened `build_new_curation.py` and left these two.

**Currently harmless:** `DOPE`, `Cholesterol` and `DMG-PEG2000` each resolve to exactly **1** distinct SMILES. Becomes order-dependent the moment any experiment introduces a second structure under one of those names.
**Severity:** cosmetic today, latent structure error. **Confidence:** certain.

### E7 — `DMG-PEG2000` is emitted with a discrete SMILES that `blocks.csv` says cannot exist
`blocks.csv` for DMG-PEG2000: *"Polydisperse PEG conjugate; no discrete SMILES can be assigned, so left empty."* The emitted rows nevertheless carry a specific fixed-length PEG chain, pulled out of the corpus by the E6 lookup. CLAUDE.md §3 takes the folder's side (*"PEG-lipids are deliberately absent — the block length is a distribution"*). Class 21.
This is a corpus-wide convention rather than something `LL_2026` invented, so I flag it as an inconsistency to resolve in one place.
**Severity:** metadata. **Confidence:** certain that the two disagree; no view on which should win.

### E8 — the recorded "one-sided −0.52 pp bias" in the Fig 2f/g/h check is a property of the detector, not of the figure
`11_LL_2026_qc.json` records *"max |delta| 1.0 pp, mean −0.52 pp … the residual −0.52 pp one-sided bias is recorded, not tuned away"*, and `08_fig2fgh_qc.json` shows **48/48 deltas negative**.

My independent measurement gives **max |delta| 0.65 pp, mean −0.10 pp, 27/48 negative** — a two-sided residual. Where the two scales agree almost exactly (panel h: theirs 5.347, mine 5.346 px/%), their values are still ~0.5 pp lower, so the gap is in the bar-top detector, not the calibration. I profiled one top edge pixel-by-pixel: fill `(242,218,218)`, background `(254,254,254)`, ramp over ~4 px; the 50 % crossing is at y ≈ 3531.7, while a threshold detector lands at 3532–3534. **The whole disputed quantity is 1–3 px of edge definition, i.e. 0.2–0.6 pp.** That is the resolution floor of this check; neither number is "wrong", but recording −0.52 pp as a one-sided property of the figure overstates what the measurement can resolve.
**Severity:** metadata. **Confidence:** probable.

### E9 — the 48-bar check is narrower than the phrase "48 printed hit rates" suggests

(a) **They are not printed.** The paper prints exactly one hit-rate statistic ("over 85 %", which I reproduced). The 48 are bar heights in Figures 2f/2g/2h — graphically encoded, recovered by our own pixel code. They are still genuine external ground truth, but calling them "printed" invites the reader to treat them as digits read off the page.

(b) **What they constrain.** The 48 numbers are 24 group marginals × 2 thresholds. They pin the *count* of cells ≥ bin 3 and = bin 4 within each group. They say nothing about *which* cell inside a group, and **nothing whatsoever about the bin-1 vs bin-2 distinction**, since both sit below 1e5. That is **179 of 840 cells (21 %)** whose assignment this validation does not touch at all. Those rest entirely on the colour match, which is strong here (palette separation 50.2, observed max distance 2.24, min runner-up margin 49.73) — but the claim should say so.
**Severity:** metadata / claim scope. **Confidence:** certain.

### E10 — the Fig 2a defect was visible in the curation's own review overlay
`review/fig2a_overlay.png`, right-hand end: the green body box runs past 12A2 into white and the red sampling boxes in the last columns straddle their cell boundaries. The overlay was produced; its most informative region was not looked at. **An overlay only helps if the check includes *where to look*** — here, both extreme ends of every axis.
**Severity:** process. **Confidence:** certain.

### E11 — `LL_2026` and `LL_2026_2` share 70 lipids with byte-identical SMILES, across two `Experiment_ID`s. This lands in Stage 5.
Not a curation error — the two really are different experiments (different formulation and IL:mRNA ratio) and the separation is correct. But it has a live consequence that nothing in the folder flags:

- 987 rows, 917 unique SMILES → **70 duplicates**, all the same lipid measured in both screens.
- Concordance over those 70: **Spearman 0.496** (p = 1.3e-5), Pearson 0.672, 52/70 same bin.

When the frozen pipeline is re-run on the current corpus, this becomes a **70-point cross-"source" pair in Stage 5's pooled shared-lipid concordance** — likely the second-largest contributor after AA_2008↔KW_2014 — and it is a *same-lab, same-paper, same-cell-line* pair that differs only in formulation. Pooling it with genuine cross-laboratory pairs would pull Result 3's headline upward for a reason that has nothing to do with cross-laboratory transfer. Same class as the open §8 decision about benchmark-only pairs in the in-vivo pool.

Read the other way, 0.496 is a *useful* new rung for the paper's own argument: same lab, same cells, only the formulation changes → 0.50, against molea's same-lab cross-cell-type 0.30 and the cross-lab 0.083.
**Severity:** changes a reported value, if Stage 5 is re-run without a decision. **Confidence:** certain on the numbers; the pooling question is a judgement.

### E12 — `digitize_fig2a.py`'s stated reason for indexing the bins 2–4 is vacuous
Docstring: *"Bins are therefore indexed 2-4 … rather than 1-3 which would silently shift every value by a decade."* The stored `Experiment_value` is a per-experiment z-score, invariant to an additive shift in the ordinal — 1-3 and 2-4 give bit-identical z-scores (confirmed; spacing 1.4043 either way). The *choice* is nonetheless correct, for a reason the docstring does not give: it makes the bin index directly comparable to Fig 2e's, which `07b_fig2a_qc.json`'s 70-lipid overlap table depends on. And item 4 above is independent evidence for it.
**Severity:** cosmetic. **Confidence:** certain.

---

## 3. UNVERIFIABLE — and what would settle each

1. **The branch isomer of iso-octyl / iso-decyl acrylate** (Ac8-1, Ac10-1; 140 rows). The SI gives only trade names; both are commercial oxo-alcohol esters sold as mixtures. Correctly flagged `confidence = medium` + `structure_assumed = 1`. *Settled only by:* the supplier's CoA or a GC-MS of the lot.
2. **Which N-H actually survives to take the aromatic arm for 6B and 7B** (282 rows). Fig S1 reverses the order for these two, so in reality the leftover site is whichever is least reactive, plausibly an internal secondary nitrogen rather than the terminal one Fig S1 draws. We follow the drawing. *Settled only by:* NMR or LC-MS/MS of an actual 6B/7B product.
3. **Regiochemical and substitutional purity of any library member.** The SI is explicit that final products are crude. Every stored structure is the idealised fully-capped single regioisomer. *Settled only by:* characterisation the paper does not contain.
4. **Enantiomeric excess of the 5B reagent and of the 5B products.** "(S)-glycidyl trityl ether" states a configuration, not a purity. *Settled only by:* the reagent CoA plus chiral HPLC. This is why F-009's deliverable is a flag.
5. **Whether any 5B product retains its trityl group.** The general method has no acidic step, so retention is expected and the corpus stores it retained. No 5B compound is characterised.

---

## 4. SEVERITY

| # | finding | changes a structure | changes a value | metadata only | cosmetic |
|---|---|---|---|---|---|
| E1 | F-008, 57 px window overrun | — | no (0/147 bins move) | — | **yes**, latent off-by-N |
| E2 | F-009, 5B stereocentre unflagged | **stereochemistry** (flag, not rewrite) | — | — | — |
| E3 | `positional_ambiguity` + wrong note on 987 rows | — | — | **yes** | — |
| E4 | two naming vocabularies across the sibling IDs | — | — | **yes** | — |
| E5 | hard-coded `"high"` at `enumerate_147.py:111` | metadata, **once F-009 lands** | — | — | yes, today |
| E6 | class-30 `setdefault` lookup | latent | — | — | yes, today |
| E7 | DMG-PEG2000 SMILES vs `blocks.csv` | — | — | **yes** | — |
| E8 | "−0.52 pp one-sided bias" overstated | — | — | **yes** | — |
| E9 | "48 printed hit rates" — scope of the claim | — | — | **yes** | — |
| E10 | overlay generated, right-hand end not inspected | — | — | — | process |
| E11 | 70 lipids shared across the two IDs | — | **yes, if Stage 5 is re-run undecided** | — | — |
| E12 | vacuous docstring rationale | — | — | — | **yes** |

**No stored structure and no stored `Experiment_value` in either experiment is wrong.** 840/840 Fig 2e bins, 147/147 Fig 2a bins, 20/20 spot-checked structures, 4/4 HRMS anchors, 3/3 intermediate MS values, and both z-score vectors all reproduce from scratch.

---

## 5. CONFIDENCE

| # | confidence | basis |
|---|---|---|
| E1 | **certain** | edge positions and 22 boundaries measured three ways; both geometries re-digitized; visually confirmed in my overlay and in the curation's own |
| E2 | **certain** | single `(S)` in the whole SI; zero hits for any chiral assay; regiochemistry confirmed by the paper's own NMR integrals; CIP outcome computed, not reasoned |
| E3 | **certain** on both factual errors in the note; the flag *value* is a judgement | 24 exhaustive product enumerations; two verbatim captions; 980/987 have no ring nitrogen |
| E4 | **certain** | read straight out of `data_vitro.csv` |
| E5 | **certain** on the construct and on its current correctness | code read + 147/147 `high` in the corpus |
| E6 | **certain** | code read + distinct-SMILES count over the live corpus |
| E7 | **certain** that the folder's text and its output disagree | — |
| E8 | **probable** | two honest detectors differ by ~0.5 pp; the true edge is definitional at 1–3 px |
| E9 | **certain** | the paper prints one hit-rate statistic; 179/840 cells lie below both thresholds |
| E10 | **certain** | looked at the stored overlay |
| E11 | **certain** on the numbers; the pooling decision is yours | 70 identical SMILES, Spearman computed from the corpus |
| E12 | **certain** | z-score is affine-invariant; verified numerically |

---

### Note on my own errors, since the program asks for it
Two intermediate results were wrong before I caught them, both classes this project already tracks. (1) My first independent structure build reported a mismatch for `6A2-3B-EP10` and a spurious second product for `7A1`; both were **ring-closure collisions in my own SMILES templates** — the arm's `c1ccccc1` colliding with a benzene and a piperazine numbered `1` — the identical hazard that fused piperazine to cholesterol in the Siegwart curation. The corpus was right both times. (2) My hand CIP analysis of the 5B carbinol concluded the descriptor flips on ring opening; it does not. I had followed the epoxide ring closure to the real stereocentre instead of to a CIP duplicate atom. `rdCIPLabeler` settled it. Both are instances of the standing rule: **run the primary object, do not reason about it.**
