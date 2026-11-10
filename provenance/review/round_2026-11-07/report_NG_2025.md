# Tier-2 adversarial review — `curation/ng_2025/` (Gong et al., *Nat Biomed Eng* 2025, 9:2181–2195)

**Provenance.** Returned as text by a read-only reviewer and saved verbatim.
**Its findings are claims.** My adjudication — what reproduced, what I applied,
what I held — is in `FINDINGS.md` under F-048…, not here.

---

## VERDICT: **MERGE WITH FIXES**

The 333 emitted rows are, as far as I can break them, right. I re-derived both experiments' `Experiment_value` vectors from the publisher's Source Data workbook and reproduced the in-vitro block to max |Δ| = 5e-7 and **both in-vivo blocks exactly (|Δ| = 0.0)**; I re-derived the phenol tail lengths, the amine cores, the thiol tail lengths, the C-a16 connectivity and the B-a1 anchor by my own geometry and my own NMR assignment, and found no emitted structure and no emitted value wrong. The B-a1 anchor claim is **genuinely executed** (I confirmed both that `vector_smiles()` reads `01_fig1d_graph.json` at run time *and*, independently, that the Fig. 1d vector bond segments encode exactly that structure).

**The single most important finding is negative, and it is about the flags, not the data: the README's stated justification for flagging all 104 dialkyl-ketone lipids `confidence = low` does not reproduce.** The README asserts "15 of 52 [1:2:1 lipids] are drawn as the open-chain mono-adduct — I 6 (b4, b5, b7, b9, b10, b11), J 1 (b5), K 4 (b1, b2, b7, b12), L 4 (b1, b2, b4, b10)". I rendered and looked at nine of the fifteen named lipids in Supplementary Figs 8 and 9 — **I-a-b4, I-a-b5, I-a-b7, I-a-b9, I-a-b10, I-a-b11, K-a-b2, K-a-b12, L-a-b4** — and **every single one is drawn as the cyclic 4-oxopiperidine-3,5-dicarboxylate**, with a closed six-membered ring, exactly as the scheme and exactly as enumerated. The count that error class 25 was written around is itself wrong. This does not change a structure (the emitted ring form is *more* supported than the curation believed), but it means a `confidence = low` on 104 rows would be carrying a reason that is false — the precise failure mode §6 warns about when it says provenance text must not assert what the code/data does not.

---

## 1. CONFIRMED CORRECT

**C1. In-vitro `Experiment_value` (310 rows) — re-derived from the workbook.**
Read `41551_2025_1422_MOESM3_ESM.xlsx` sheet `Fig1` myself, mapped the six blocks to (tail, amine) by my own rule, and computed `z = (x − mean)/sd(ddof=0)` on the **already-log10** workbook values. Max |mine − emitted| = **4.999e-7** (= 6-dp rounding in the CSV); ρ = 1.000000. Control: z of the *linear* luminescence differs by up to 9.04, and ddof=1 differs by up to 4.75e-3 — so the population-sd-of-log10 convention is confirmed, not merely consistent.

**C2. In-vivo `Experiment_value` (23 rows) — exact.**
`z(log10(arithmetic mean of the three replicates))`, ddof=0, per block. Max |mine − emitted| = **0.000000** for both `liver` (n=12) and `muscle` (n=11). (The alternative `mean(log10(·))` differs by up to 0.015, so the aggregation step is pinned too.)

**C3. The workbook→lipid mapping, independently (this is the big one).**
The Fig. 1e–j heat maps are **vector rects with exact fill colours**, not raster. I extracted all 310 filled cells from page 3's drawing list, binned them into six grids by position, reduced each panel's RGB cloud to its first principal component (variance fraction 1.0000 in every panel — the colormap is exactly 1-D), and tested *every* column permutation against the workbook:

| panel | identity-permutation \|ρ\| | best competing permutation | rows reversed |
|---|---|---|---|
| Fig.1e | **0.9996** | 0.9124 | 0.4390 |
| Fig.1f | **0.9995** | 0.8666 | −0.1556 |
| Fig.1g | **0.9981** | 0.8845 | −0.0035 |
| Fig.1h | **0.9997** | 0.7213 | 0.0601 |
| Fig.1i | **0.9999** | 0.6394 | −0.2622 |
| Fig.1j | **0.9999** | 0.7143 | 0.8296 |

The identity mapping is uniquely best in all six, and the row order (a1…a16 / b1…b13 top-to-bottom) is likewise pinned. Panel letters are fixed by position; tail letters by the PDF text layer, whose x-coordinates I read: `A B C D` at x 321/340/359/378 under panel e, `E F G H` under f, **`I1 J1 K1 L1`** under g, **`I2 J2 K2 L2`** under h, `M N O` under i, `P Q R` under j. This reproduces the curation's C7 numbers by a completely different route.

**C4. Phenol tails A/B/C/D = 0/6/8/9 carbons — geometric, three independent ways.**
(a) Fig. 1b is vector; connected-component analysis of its bond-length segments gives **7 / 13 / 15 / 16** segments = (6 ring + 1 C–O) + 0/6/8/9 chain bonds. (b) Supplementary Fig 17c's x-axis ticks read **0, 6, 8, 9** ("Number of carbon in the alkyl chain") — a part of the paper not used to derive them. (c) Methods print "**p-octylphenol**" for the C-a16 synthesis.

**C5. Tail C is *n*-octyl, not *tert*-octyl.** 4-*tert*-octylphenol is at least as common a catalogue item as 4-*n*-octylphenol, and formula/mass cannot tell them apart. The printed C-a16 NMR settles it: **δ 0.90 (t, J = 6.8 Hz, 12H)** is a linear terminal CH₃; *tert*-octyl would give a 9H singlet near 0.70 and no triplet. This test did not exist in the curation's QC and it closes a real hole.

**C6. C-a16 connectivity pinned by a full NMR assignment** (my own, beyond the curation's "integrals sum to 106"): 6.99–6.82 (8H) + 6.72–6.67 (4H) = 12 ArH = four 1,2,4-trisubstituted rings; **4.05–4.14 (2H)** = the two diketopiperazine α-CH; **3.86–3.43 (8H)** = four Ar–CH₂–N (⇒ exactly **4 arms**, not 2 and not 6); 2.47–2.66 (12H) = 4×ArCH₂ (8H) + 2×ε-NCH₂ (4H), which requires the lysine ε-N to be **tertiary**; 2.06–1.61 (8H) + 1.60–1.39 (12H) = 4×octyl β-CH₂ (8H) + 2×lysine β,γ,δ (12H); 1.29 (40H) = 4×5 CH₂; 0.90 (12H) = 4 CH₃. Every multiplet is accounted for with zero residue. Independently, Methods run the synthesis at **4.5 phenol : 4.2 HCHO : 1.0 amine** = 4 arms. `C72H112N4O6`, MW 1129.71 vs the printed 1129.71.

**C7. Amine cores a1–a16 — 16/16 by eye off a 900-dpi render of Fig. 1c**, matching `02_blocks_hand.csv` exactly, including a16 = cyclo(Lys-Lys). a2 = `NCCN(CC)CC` is pinned by the printed sentence "N,N-diethylethylenediamine (referred to as amine a2 … and b10 …)".

**C8. Amines b1–b13 — 13/13** off Supplementary Fig 7c; b10 = `NCCN(CC)CC` confirmed by the same printed sentence.

**C9. Degree of substitution.** Supplementary Fig 1 draws all five product types; Supp Figs 4a, 7a, 10a, 15a print the stoichiometries ("2 eq. thiol : 2 eq. HCHO : 1 eq. amine"; "1 eq. ketone : 2 eq. HCHO : 1 eq. amine"). I verified on the drawn products that **secondary N–H is left intact** (A-a12 and A-a15 both show the central N–H unsubstituted with 2 arms on each terminal NH₂) and that the **6-arm case is real** (A-a13, tris(2-aminoethyl)amine, three N's each bearing two 2-hydroxybenzyl arms). Phenol attachment is **ortho to OH** in every drawing I opened, including tail A where para is chemically available.

**C10. Thiol tails E/F/G/H = C6/C8/C10/C12 — my own trace** of the native 876×825 raster: midline-crossing counts 6 and 8 for E and F, and chain-bottom y at 179/225/269/314 px (Δ = 46/44/45 px = exactly 2 carbons at ~23 px/bond) fixing G and H at 10 and 12. These are the four standard commercial n-alkanethiols.

**C11. Both ketone ring systems** — Supplementary Fig 1 draws the 1:2:1 product as 1-R-4-oxopiperidine-3,5-dicarboxylate and the 1:4:2 product as the 1,5-bis(ester) **bispidinone** with two N–R; the per-lipid panels I opened agree (see §2, E1). Ester chains I/J/K/L = ethyl/butyl/hexyl/octyl, read off Supp Fig 7b and confirmed on the drawn products.

**C12. The B-a1 C10 anchor is real.** `stage08_reference.py` calls `vector_smiles('01_fig1d_graph.json', 20)` at import time; the JSON carries genuine page coordinates (`pdf`, `page: 2`, `clip`, 69 bond segments, 11 glyphs) and its largest component has **35 atoms**, which is exactly B-a1's heavy-atom count. Separately I rebuilt Fig. 1d from the raw PDF line items myself: the reagent phenol gives 13 bond-length segments (= hexyl), each product arm gives 15 (6 ring + ring–O + ring–CH₂ + CH₂–N + 6 chain = hexyl), and the linker component gives 4 segments between the two nitrogens (= propylene) plus two N–CH₃. The anchor is derived, not asserted, and it is correct.

**C13. Supplementary Fig 17c transposition — reproduced.** The bars stand at **0 / 3 / 6 / 2** over ticks 0/6/8/9, matching the Results text and contradicting the Source Data (0/2/6/3) and the paper's own 11-lipid in-vivo panel (2 B-lipids, 3 D-lipids). Confirmed in two independent places, as the README says.

**C14. The hit rule.** Threshold log10 = 5.938 over the 64 phenols returns exactly the 11 in-vivo lipids; the near-miss that sets the window's floor is B-a5 at 5.93774. Over all 310 it also catches G-a12 at 6.07403 (thiol), which was never taken into mice. Both claims reproduce from the workbook.

**C15. Metadata, field by field** (P = printed fact, R = reading):

| field | value | P/R | source |
|---|---|---|---|
| `Model_type` (vitro) | MC38 | **P** | "A mouse colon cancer cell line MC38 was used for in vitro screening" |
| `Dose_ug_nucleicacid` (vitro) | 0.01 | **P** | "a dose of 10 ng per well mRNA" |
| `Mixing_method` (vitro) | handmixed | **P** | "were mixed using pipette mixing" |
| `Mixing_method` (vivo) | microfluidics | **P** | "animal experiments … Micro&Nano microfluidic chip" |
| `Dialysis_buffer` (vivo) | PBS | **P** | "dialysing against 1× PBS … for 2 h" |
| `Dialysis_buffer` (vitro) | blank | **R** | the dialysis sentence attaches to the microfluidic clause; defensible |
| `Aqueous_buffer` | citrate | **P** | "10 mM citrate buffer (pH = 3)" |
| molar ratios 50/10/38.5/1.5 | DOPE | **P** | "ionizable lipid, DOPE, cholesterol and DMG-PEG … 50/10/38.5/1.5" |
| MC3 row `HL_name` = **DSPC** | — | **P** | "For the preparation of DLin-MC3-DMA LNPs, DLin-MC3-DMA, **DSPC**, cholesterol and DMG-PEG … 50/10/38.5/1.5" ✔ the curation got this right |
| `Model_type` (vivo) | Mouse_B6 | **P** | "6–8-week-old female C57BL/6J mice … (around 20 g)" |
| `Dose_ug_nucleicacid` (vivo) | 2.50 | derived | 0.125 mg kg⁻¹ (**P**) × 20 g (**P**) |
| `Model_target` = liver, `Route` = IV | — | **P** | "intravenously injected … luciferase signal was mainly detected in the liver"; Fig 1m caption "luminescence levels in livers" |
| `Model_target` = muscle, `Route` = IM | — | **P** | "the luminescence level **at the injection site** was assessed 6 h post-injection (Fig. 1n–p…)" — `muscle` is correct, not `whole_body` |
| `Cargo`/`Cargo_type` | mRNA / FLuc | **P** | firefly luciferase mRNA throughout |
| `Experiment_batching` = individual | — | **R** | correct: Supp Fig 17a is the *pooled* comparison, Fig 1e–j are not |
| `IL_to_nucleicacid_massratio` blank | — | **R** | "lipid and mRNA at a mass ratio of 20:1" does not say total vs ionizable; leaving it blank is right |

**C16. Component identities.** `DOPE`, `DSPC`, `Cholesterol`, `DLin-MC3-DMA` all canonicalize **identically** to `provenance/reagent_registry.py` and to the corpus's existing strings; `DMG-PEG2000` matches the corpus's single existing variant and is deliberately not in the registry. C11 will pass. All vocabulary values already exist in the corpus; `MC38` is a genuinely new `Model_type` (vitro goes 16 → 17 cell types).

**C17. Table hygiene, re-run independently.** 0 unparseable, 0 non-canonical, 0 radicals, 0 disconnected fragments, **0 padded cells** (class 37), 0 duplicate `(Experiment_ID, IL_name, Model_target, Model_type)` keys, 0 names carrying two SMILES, 0 SMILES carrying two names. Corpus overlap is **exactly 1** — `DLin-MC3-DMA`, which matches the corpus string bit-for-bit and uses the majority spelling. All descriptor columns are blank, which is correct: `build_dataset.py` recomputes them.

**C18. The held-back row is correctly held back.** I confirmed the contradiction from the primary objects, not the README: the PDF text layer carries **two** rotated `DLin-MC3` axis labels on page 3, at x = 202.5 (panel m) and **x = 477.8 (panel p)**, while the workbook's Fig.1m header reads `Dlin-MC3-DMA` and its Fig.1p header reads `SM-102`. The two headers were typed differently, so the difference is deliberate, not a copy. The paper used both benchmarks and associates SM-102 with the intramuscular route elsewhere — which tips the balance toward the workbook header and *away* from the axis, but does not settle it. Holding the row back is the right call and I would not reverse it.

---

## 2. ERRORS FOUND

**E1 — the open-chain ketone count does not reproduce.** *(certain; documentation / changes a flag's stated reason, not a structure)*
File: `README.md` §"Three things the paper contradicts itself about" item 2, and `stage10_flags.py` docstring.
Claim: "In the 1:2:1 library **15 of 52** are drawn as the open-chain mono-adduct — I 6 (b4, b5, b7, b9, b10, b11), J 1 (b5), K 4 (b1, b2, b7, b12), L 4 (b1, b2, b4, b10)".
Evidence: I extracted the native rasters of Supplementary Figs 8 and 9 (`995×1454` and `1376×1429`) and rendered nine of the fifteen named structures at 3× with inversion. **All nine are drawn with the closed six-membered 4-oxopiperidine ring**: I-a-b4, I-a-b5, I-a-b7, I-a-b9, I-a-b10, I-a-b11, K-a-b2, K-a-b12, L-a-b4. That is 6/6 of the tail-I list and 3/3 of the others I sampled.
Correct value: unknown, but **much lower than 15, and 0 for tail I**. How to verify: classify all 104 drawn products mechanically — the cyclic form encloses one extra hole in the ink, so a per-structure connected-component hole count, compared across the four tails for the *same* amine, decides it without eyes.
Why it matters: this is the stated reason 104 rows would carry `confidence = low`. A flag whose reason is wrong is worse than no flag, by the project's own F-011 logic.

**E2 — `stage10_flags.py`'s docstring and the README disagree with each other** *(certain; cosmetic)*. The docstring says "of tail I's 13 lipids at 1:2:1, **five** are drawn as the open-chain mono-adduct, while tail J's are drawn cyclic for 12 of 13"; the README says six for I and 12 of 13 cyclic for J. Both are secondary artifacts of a count nothing implements — the exact pattern of error class 39. Neither survives E1.

**E3 — README says all rows carry DOPE; one carries DSPC.** *(certain; cosmetic, data is correct)*
`README.md`: "Every row here carries 50 / 10 / 38.5 / 1.5 (IL / **DOPE** / cholesterol / DMG-PEG)". The emitted table correctly has `HL_name = DSPC` on the single `NG_2025_Fvivo_MC3` row, matching the printed Methods sentence. Fix the prose, not the data.

**E4 — `Formulation_ID` prefix does not match `Experiment_ID` on the 23 in-vivo rows.** *(certain; cosmetic but mechanically checkable)*
Rows with `Experiment_ID = NG_2025_2` carry `Formulation_ID` = `NG_2025_Fvivo` / `NG_2025_Fvivo_MC3`. Across the live corpus only **292 of 30,665** rows violate "`Formulation_ID` starts with `Experiment_ID` + `_`", and those are a single upstream oddity (`AP_2025` carrying `JGW_2024_*`). Rename to `NG_2025_2_Fvivo` and `NG_2025_2_Fvivo_MC3` before merge.

**E5 — the in-vivo blocking has no vocabulary term, and will poison `Experiment_value_within_model`.** *(certain; metadata/convention)*
`NG_2025_2` is z-scored per (`Experiment_ID` × **`Model_target`**), which is *not* one of `schema_quality.json`'s three allowed `value_normalization_unit` values. The precedent the README cites is real — I checked `ZH_2023` and `BL_2023` in `data_vivo.csv` and **each organ sub-block is exactly mean 0 / sd 1**, so those two are per-organ-standardised too — but both are labelled `experiment`, which the corpus should probably not repeat. Separately: `NG_2025_2` has one `Model_type` (`Mouse_B6`) spanning two already-standardised blocks, so `build_dataset.py` will compute `Experiment_value_within_model` **over the union of two blocks whose offset has already been deleted** — error class 41, exactly.
Two facts to inform the decision: the between-block offset here is **recoverable** (both panels are in the same unit, radiance p s⁻¹ cm⁻² sr⁻¹, with the ×10⁶ / ×10⁵ multipliers printed), so the §3 "only legitimate where unrecoverable" exemption does not apply on its face; and **liver-vs-muscle Spearman over the 11 shared lipids is 0.982** (p = 8.4e-8), so pooling would cost almost nothing in rank terms. That 0.982 is also a free manuscript datum: a same-lab, same-lipid, cross-route, cross-organ concordance point, sitting above molea's 0.30 and far above the cross-lab 0.083.

**E6 — `dialkylate()` bonds to atom index 0 of the fragment.** *(probable; process, no current error)*
`stage04_enumerate.py` does `m.AddBond(n, off, …)` where `off` is the insertion offset, i.e. **atom 0 of the arm SMILES**. CLAUDE.md §6 bans exactly this. I checked all 14 `ARM` strings and atom 0 *is* the intended attachment carbon in every one, so nothing is wrong today — but the anchor builder in the same folder already does it correctly with a `*` dummy, so the fix is one line and the inconsistency is a live trap for the next editor.

**E7 — no C10 anchor for the thiol library.** *(certain; missing coverage)*
The six anchors cover phenol (×2), pyrrole, acetophenone and both ketone rings. The **64 thiol rows have none**, although Supplementary Figs 5 and 6 draw all 64 products in full.

**E8 — the paper's own stoichiometry is written in two orders, and only one gives a sensible product.** *(certain; paper-level, curation chose correctly)*
Main Fig. 1g/1h panel titles read "**Amine/formaldehyde/dialkyl ketones** = 1:2:1" and "= 1:4:2"; the Fig. 1 caption and Supplementary Figs 7/10 read "**dialkyl ketone : formaldehyde : amine** ratio of 1:2:1 / 1:4:2". For 1:2:1 the orders coincide; for **1:4:2 they do not** — the panel title would mean 1 amine : 4 HCHO : 2 ketones, which closes nothing, while the caption/SI order means 1 ketone : 4 HCHO : 2 amines, which is the bispidinone the SI actually draws. The curation followed the SI and is right.

**E9 — Supplementary Fig 2 draws A-a2 with the a1 backbone.** *(certain; paper-level, curation chose correctly)*
I measured the drawn chain in Supp Fig 2a: **A-a1 and A-a2 both show three carbons between the two nitrogens**, but a2 is `N,N`-diethyl**ethylene**diamine by the printed main-text sentence and by Fig. 1c. The product panel is a template artefact. The curation emitted the Fig. 1c reading, which is right.

---

## 3. UNVERIFIABLE

**U1 — amines b3 and b4: aminal or hydrazine (28 lipids).** **Aminal** is drawn in the four amine-list panels **and also in the four ketone product panels** — eight panels, not four. **Hydrazine** is drawn unambiguously in the pyrrole and acetophenone product panels. Chemistry pulls the other way: an aminal bearing a free NH₂ is not an isolable, purchasable reagent, whereas 1-amino-4-methylpiperazine is an ordinary catalogue item — though 1,1-dimethylhydrazine is a restricted/high-hazard chemical. **No mass, formula or NMR is given for any b3/b4 lipid**, and the two options differ by CH₂ = 14 Da. Keep the curation's choice, keep both strings, flag the rows.

**U2 — the true open-chain count among the 104 ketone lipids.** E1 refutes the stated list; I did not classify all 104.

**U3 — all 310 in-vitro lipids are crude.** Printed, verbatim: "the solvents were evaporated to obtain **310 different crude lipid compounds**. These crude compounds were subsequently utilized to prepare LNPs for in vitro screening." The emitted SMILES is the nominal major product of a one-pot reaction that was never purified or characterised (only C-a16 has any characterisation at all). **The curation flags it nowhere.** It is the largest single gap in the flag design.

**U4 — tail D, 4-nonylphenol.** Drawn linear n-nonyl; commercial "4-nonylphenol" is predominantly a **branched technical isomer mixture**, and unlike tail C there is no NMR to break the tie. 16 lipids.

**U5 — the "~15-fold" claim.** A 15× threshold implies MC3 ≈ 3.5e5 and only four hits, while reproducing the eleven needs MC3 in (1.73e5, 2.18e5]. Supplementary Fig 17f plots MC3 at roughly 2.3e5 — the two numbers most likely come from different runs.

**U6 — O-b10.** Source Data make it the highest value in the pyrrole panel; Supplementary Fig 17b re-measures it ~24× lower while five companions agree within 1.9×. Not a column-assignment artefact (Fig 1i pinned at |ρ| = 0.9999).

**U7 — IM timepoint.** Main text "6 h"; Supplementary Fig 18 caption "After 24 h". Same experiment by the main text's own cross-reference. No schema column.

**U8 — mRNA chemistry.** Results: "5-methoxyuridine-modified"; Methods: "1-methyl pseudouridine-modified". No schema column.

**U9 — the Fig. 1p benchmark.** See C18. Genuinely unresolvable; held back correctly.

---

## 4. SEVERITY

| # | severity |
|---|---|
| E1 | **changes a flag's stated reason** (and, if acted on, the `confidence` of 104 rows). No emitted structure or value changes. |
| E2 | cosmetic |
| E3 | cosmetic (prose only; the data is right) |
| E4 | cosmetic, but it breaks a near-universal corpus invariant — fix before merge |
| E5 | **metadata / corpus-convention**; could change how 23 rows are consumed |
| E6 | cosmetic today, latent structural hazard |
| E7 | missing verification coverage on 64 rows |
| E8, E9 | cosmetic (record-keeping); both confirm the curation's choices |
| U1 | would **change 28 structures** if reversed |
| U3 | **changes the confidence of all 310 in-vitro rows** |
| U4 | would change 16 structures if reversed |

## 5. CONFIDENCE

| # | confidence |
|---|---|
| E1 | **certain** — nine rendered structures, each looked at |
| E2, E3, E4, E8, E9 | certain |
| E5 | certain on the facts; the *decision* is the caller's |
| E6 | certain that the construct is there; certain that it is currently harmless |
| E7 | certain |
| U1 | the drawings are certain; **which reagent is real is speculative either way** |
| U3 | certain (printed verbatim) |
| U4 | probable; unknowable for this paper |
| U5–U9 | certain as *contradictions*; unresolvable as facts |

---

## 6. RECOMMENDATION FOR THE 7 QUALITY COLUMNS

### `value_normalization_unit` — all 333 rows, required
- **310 in-vitro rows** → `experiment`.
- **23 in-vivo rows** → the E5 decision. *(recommended)* **pool the two organs into one z** and set `experiment`: the offset is recoverable, the rank cost is nil (ρ = 0.982), and it makes `Experiment_value_within_model` meaningful. *(status quo)* keep per-organ blocks and set `experiment_x_model_type` with a note.

### `confidence`
- **`medium` on all 310 in-vitro rows** (U3) — nominal products of unpurified crude mixtures. Same for the **11 crude in-vivo rows**; **MC3 row `high`**.
- **`low` on the 28 b3/b4 lipids** (U1).
- **Do *not* put `low` on the 104 ketone lipids on the "drawn open-chain" grounds** (E1).
- Net: **310 vitro = 28 low + 282 medium**, **23 vivo = 22 medium + 1 high**.

### `positional_ambiguity`
- **1 on all 310 in-vitro rows and on the 11 crude in-vivo rows (321 rows)**, `0` on MC3. A Mannich on a phenol gives ortho *and* para products, both ortho positions of a 4-alkylphenol are open, and nothing was purified.

### `structure_assumed`
- **1 on the 28 b3/b4 lipids**; **1 on the 16 tail-D lipids** (U4). Tail C needs **no** flag because the NMR triplet rules out *tert*-octyl.

### `value_is_bound`
- **Blank on all 333 rows.** Every in-vitro value comes from the Source Data spreadsheet; the minimum is 0.845098 (= log10 7) with no repeated clamp value. Nothing is censored.

### `curation_note` — per group
1. **All 310 in-vitro:** *"Structure is the nominal product of an unpurified crude Mannich reaction; the paper screened 310 crude compounds and characterised only C-a16. Value is the authors' own Source Data for Fig. 1e–j, not a digitisation."*
2. **The 28 b3/b4** (append): the aminal/hydrazine contradiction.
3. **The 16 tail-D** (append): branched technical isomer mixture.
4. **O-b10 only:** the Supplementary Fig 17b disagreement.
5. **The 11 crude in-vivo:** same crude lipids; dose derived from 0.125 mg/kg × ~20 g.
6. **The MC3 in-vivo row:** commercial benchmark, DSPC not DOPE; the matching Fig. 1p benchmark is held back.
7. **The 104 ketone lipids:** no ring-system note, or one saying every per-lipid drawing checked agrees with the scheme.

### Other things to do before merge
- Rename the two in-vivo `Formulation_ID` values (E4).
- Correct the README's DOPE/DSPC sentence (E3) and retract the "15 of 52" count (E1/E2).
- Add one thiol C10 anchor (E7).
- Record E8 and E9 in the curation's notes.
- For CLAUDE.md §8: this merge takes **DLin-MC3-DMA from 5 to 6 in-vivo sources**, and here MC3 ranks **7th of 12** in the liver block (z = +0.246) rather than at the top — a useful counter-example to the "benchmark-only pairs push concordance toward +1" worry.
