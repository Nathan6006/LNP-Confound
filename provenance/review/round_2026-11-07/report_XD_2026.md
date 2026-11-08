# Wave-1 review — `XD_2026` (Duan et al., *Acta Pharm Sin B* 2026, PMID 42592089, LipidAI)

**Provenance.** Returned as text by a read-only reviewer and saved verbatim.
**Its findings are claims.** My adjudication is in `FINDINGS.md`.

---

**Summary.** **None of the nine E-lipid structures disagrees with ours.** I re-traced all twelve structures from the embedded native Fig 5A raster with my own geometry, assembled SMILES by my own string construction, and compared by **full InChIKey**: 12/12 exact, stereo included. **The ALC-0315 override is correct and should stand** — the drawn head is unambiguously 3 carbons where the compound has 4, the bond lengths rule out a drawing convention, and the alternative hypotheses all fail against internal controls in the same panel. My independently assembled "as drawn" ALC-0315 reproduces the stored C10 anchor's InChIKey (`WESYLBGVEJKNRX-UHFFFAOYSA-N`) exactly. I also re-digitised all 12 Fig 5C bars from my own tick calibration (rms 0.18 px over 35 ticks) and got the stored values back to **rms Δz = 0.0027, max 0.0050**; re-read Fig 5D programmatically and reproduced both printed statistics; and confirmed no LipidDB row reached `data/`. Four findings, **none touching a structure or a value**: a retracted "299 publications" figure still asserted in `stage07_emit.py`, an unguarded corpus-lookup-by-name of exactly error class 30 in the same file, a `structure_assumed` flag used outside its own schema definition in a way that contradicts F-030 from this same round, and an unexecuted provenance claim in `stage03_structures.py`.

---

## Method (what I actually did, so you can judge the independence)

I worked from the **embedded image**, not a page render: page 13 carries a single 1541 × 1864 RGB PNG (xref 491), which I extracted losslessly. The curation digitised Fig 5C from a 600-dpi page render, which is a 2× upsample of this same PNG — so my measurements are on the original pixels and carry no interpolation. All inspection was at 8–18× NEAREST magnification with the Read tool; all counting was by a per-column ink-extremum detector I wrote (`profile → local extrema → merge runs`), cross-checked against the magnified images. Annotated overlays were built and looked at for the bar tops and for ALC-0315's full carbon skeleton.

Scales came only from axis furniture. Nothing in my pipeline reads a curation output.

---

## 1. CONFIRMED CORRECT

### 1.1 All nine E structures — independent re-trace, full InChIKey match

I measured every chain in every molecule. Bond pitch is 12.3–12.5 px for E2–E8, 13.8 px for E9, 11.0–11.5 px for the column-1 molecules; every chain I report has uniform spacing within ±1 px, and I resolved the near-vertical bonds adjacent to each tertiary N by magnified inspection (the extremum detector cannot see them — that is the one systematic trap in this panel, and it costs exactly one carbon per arm if missed).

| lipid | what I measured from Fig 5A | agrees with corpus |
|---|---|---|
| **E1** | 4-hydroxybutyl head; two **(CH₂)₆** spacers; **urea** (NH–CO–NH) each arm; tails = CH₂–CH(**hexyl**)(**octyl**) = 2-hexyldecyl | ✓ |
| **E2** | 4 × 2-hydroxytetradecyl (**12 C beyond the carbinol**, all four); amine–CH₂CH₂–NH–CO–CH₂CH₂–N(piperazine) on both sides, i.e. symmetric about the ring | ✓ |
| **E3** | 4 × 2-hydroxytetradecyl; **propyl**–piperazine–**propyl** | ✓ |
| **E4** | 4 × **tetradecyl**; propyl–piperazine–propyl; no oxygen | ✓ |
| **E5** | 4 × tetradecyl; CH₂CH₂–NH–CO–CH₂CH₂–N(CH₃)–CH₂CH₂–N(CH₃)–CH₂CH₂–CO–NH–CH₂CH₂ | ✓ |
| **E6** | 4 × 2-hydroxytetradecyl; CH₂CH₂–N(CH₃)–CH₂CH₂–CO–**N(piperazine)**…**N**–CO–CH₂CH₂–N(CH₃)–CH₂CH₂ (bis-acylpiperazine) | ✓ |
| **E7** | as E6 with tetradecyl tails | ✓ |
| **E8** | 4 × 2-hydroxytetradecyl; CH₂CH₂–N(CH₃)–CH₂CH₂–CO–NH–CH₂CH₂–NH–CO–CH₂CH₂–N(CH₃)–CH₂CH₂ (no ring) | ✓ |
| **E9** | 4 × N–CH₂–CO–NH–**dodecyl**; propyl–piperazine–propyl | ✓ |

I then built these as SMILES in my own notation and compared InChIKeys against `data/data_vivo.csv`:

```
E1 NMAVUZZYJAKLPH  E2 CMQHGIZKEZAQLN  E3 FTGMRNQBRRDPKP  E4 KSYAQVXWXINEEZ
E5 YKSYJAIFMNBXSQ  E6 KQZORBMPXVHRNP  E7 PATHEZXQHLBCLT  E8 HSSNQSBFSXNWAC
E9 XQVZJKYDENCKEO
```
**9/9 full-InChIKey identical** to the emitted rows, formulas identical to `03_structures.csv`.

**This does not justify upgrading `confidence` above `medium`, and it should not be read that way.** Two concordant readings of the same drawing drive *transcription* risk to ~0; they say nothing about *drawing* risk, which is the stated basis for `medium` and which this panel has demonstrably realised once. Keep `medium`.

### 1.2 Benchmarks — and the ALC-0315 verdict

| | drawn in Fig 5A (my trace) | published compound | |
|---|---|---|---|
| SM-102 | 2-hydroxyethyl head; N–(CH₂)₅–CO–O–**undecyl**; N–(CH₂)₇–CO–O–CH(**octyl**)(**octyl**) | identical | ✓ exact |
| DLin-MC3-DMA | (CH₃)₂N–(CH₂)₃–CO–O–CH(arm)(arm); **both arms C18**, *cis* double bonds between arm-C9=C10 and arm-C12=C13 counting from the carbinol-attached carbon | identical | ✓ exact |
| ALC-0315 | **HO–(CH₂)₃–N**; both arms correct **(CH₂)₆**; both acyls correct **2-hexyldecanoyl** | head is **(CH₂)₄** | ✗ one carbon |

Measured head geometry, native pixels: O(113.0) → 125.5 → 137.0 → 147.5 → N(158.0). Bond lengths **11.3 / 11.3 / 10.9 / 10.7 px**. Three vertices.

**The "drawing convention" alternative fails five ways:**
1. The heteroatom→first-vertex bond is a **full-length bond** (11.3 px, same as every interior bond). A convention that let a label absorb a CH₂ would show a short first bond — and would have to apply to every label in the panel.
2. **SM-102's head in the same panel reads 2** (vertices at 98.5, 111.2 → N at 122), which is correct. A carbon-swallowing convention would make it read 1.
3. **E1's head in the same panel, same style, 390 px below, reads 4** (95.2, 106.5, 118.5, 130.0 → N at 141), which is correct for the structure it must be. A convention would make it 3 and ALC 2.
4. The rest of ALC-0315 is drawn exactly right — both hexamethylene arms (I counted 6 on each, including the vertex at the top/bottom of the steep N bond that a naive extremum scan misses) and both branched acyl groups. A convention would not be localised to one substituent of one molecule.
5. There is **no authorial structure file that could contradict the name**: ALC-0315 does not appear in the authors' own LipidDB (I canonicalised all 485 rows and searched for both the 3-C and 4-C forms — zero hits; SM-102 is the only one of the twelve present, as serial 353), and the released repo contains no SMILES for any of the twelve validation lipids.

So the name denotes the published compound, the measurement belongs to the physical reagent, and emitting the registry structure keeps the value on the right molecule — the `AC_2025` test, answered correctly. `reagent_registry.json` holds exactly the emitted string.

One small point in the paper's favour and against the registry framing: **Methods 2.3.1 says "ALC-0315, SM-102, and MC3 are synthesized in our laboratory"**, not purchased, while the Results call them "commercially accessible ionizable lipids". The packet lists all three as C11 *commercial* reagents. The conclusion is unchanged, but the C11 note for this slice should not say "commercial reagent" on the authority of the paper.

### 1.3 Values — Fig 5C re-digitised from scratch

**Calibration from axis furniture only.** The y-axis rule sits at x = 79; tick marks protrude to its left. I detected **35 tick marks** at x = 77 and fitted nothing: anchoring at y = 1465.5 → 10⁹ with 90.00 px/decade (the four major ticks at 1465.5 / 1555.5 / 1645.5 / 1735.5, spacing **exactly 90.0, 90.0, 90.0**), all 35 ticks land on log₁₀(k) positions at **rms 0.178 px = 0.0020 decades**, major-tick residual **0.00**.

**Bar tops from flank columns.** Bars are a lavender fill (211,212,232); replicate discs are white with black rims and occlude the middle of the bar width. I took the 5 leftmost and 5 rightmost columns of each bar and required a run of ≥8 filled rows. **Left-flank and right-flank medians agreed to 0 px on 12/12 bars.**

| lipid | my top (native y) | curation top (converted) | Δ | my log₁₀ | stored log₁₀ | Δ dec |
|---|---|---|---|---|---|---|
| ALC-0315 | 1476 | 1476.0 | 0 | 8.8833 | 8.8789 | +0.0044 |
| SM-102 | 1482 | 1481.5 | 0.5 | 8.8167 | 8.8177 | −0.0010 |
| DLin-MC3-DMA | 1541 | 1540.5 | 0.5 | 8.1611 | 8.1620 | −0.0009 |
| E1 | 1495 | 1494.5 | 0.5 | 8.6722 | 8.6732 | −0.0010 |
| E2 | 1518 | 1518.0 | 0 | 8.4167 | 8.4121 | +0.0046 |
| E3 | 1590 | 1589.5 | 0.5 | 7.6167 | 7.6174 | −0.0007 |
| E4 | 1629 | 1629.0 | 0 | 7.1833 | 7.1784 | +0.0049 |
| E5 | 1607 | 1606.5 | 0.5 | 7.4278 | 7.4285 | −0.0007 |
| E6 | 1745 | 1745.0 | 0 | 5.8944 | 5.8891 | +0.0053 |
| E7 | 1763 | 1763.0 | 0 | 5.6944 | 5.6891 | +0.0053 |
| E8 | 1731 | 1730.5 | 0.5 | 6.0500 | 6.0503 | −0.0003 |
| E9 | 1746 | 1746.5 | 0.5 | 5.8833 | 5.8725 | +0.0108 |

Recomputing z(log₁₀), population sd, over my twelve values against the stored `Experiment_value`: **rms Δz = 0.0027, max |Δz| = 0.0050 (E8)**. The post-F-020 values are confirmed to half a native pixel everywhere. I also reproduced the F-020 mechanism incidentally: `inner_fill_fraction_at_true_top` in `04_fig5c.json` is 0.21–0.57 on the nine shifted bars and 0.68–0.97 on the three that did not move, so the 0.60 threshold is exactly the discriminator.

### 1.4 Fig 5D — read by cell colour, and all three printed statistics

I located the two 3 × 12 blocks geometrically (Reality x 604–845, Prediction x 878–1118, rows y 1463–1831, cell 80.5 × 30.75 px) and sampled each cell centre. Highlight = (211,212,232), background = (226,244,248); the separation is absolute.

- **Reality**: A A A A A | B B B | C C C C. My digitised values give exactly this under the printed thresholds (A > 5×10⁷, B 5×10⁶–5×10⁷, C < 5×10⁶): **12/12**. The thinnest margin is E3 at 7.6167 against the A boundary 7.69897, i.e. **0.082 decades ≈ 7.4 native px**.
- **Prediction**: A A A A A | **C C B** | C C **B** C.
- Agreement **9/12 = 75%** — matches "we achieved an overall prediction accuracy of 75 percent" and the 75/25 pie.
- Class B: E5 correct, E3 and E4 both predicted C → **1/3** — matches "only one-third of the predictions being correct".
- And the text's third, unquantified claim, **"the consistent misclassification from class B lipids to class C"**, is reproduced exactly: both B misses go to C. That sentence is an independent check on my cell read that neither the curation nor round 1 used.

### 1.5 Metadata, field by field, against the paper

| field | emitted | paper | verdict |
|---|---|---|---|
| molratios 50 / 10 / 38.5 / 1.5 | ✓ | Methods 2.3.3 verbatim | **printed fact** |
| `Mixing_method` microfluidics | ✓ | 2.3.3 | printed fact |
| `Aqueous_buffer` citrate | ✓ | 0.01 M citric acid pH 3 | printed fact |
| `Dialysis_buffer` citrate | ✓ | 0.01 M citric acid pH ≈ 6, 2 h | printed fact; the category is unique to this experiment (12 of 1,780 vivo rows) |
| `Model_type` Mouse_BALBc | ✓ | BALB/c 6–8 wk, **male**, random grouping | printed fact |
| `Route` intravenous | ✓ | Fig 5 caption "intravenously administered" / "after intravenous injection", Results "6 h after tail vein injection" | **printed fact in two places** — `07_qc.json` calls this a reading; it is more strongly supported than that |
| timepoint 6 h | n/a | 2.3.4 and caption | printed fact |
| `Model_target` whole_body | ✓ | **genuinely unresolved** — see §3 |
| `Dose_ug_nucleicacid` blank | ✓ | **absent.** The only dose is d-luciferin 150 mg/kg. I grepped for dose/dosage/µg/mg·kg⁻¹/encapsulat* — nothing. | correct to leave blank |
| N:P | blank | **absent.** The only ratio besides the lipid mol% is the aqueous:ethanol 3:1 mixing ratio. | correct to leave blank |
| PMID 42592089, DOI 10.1016/j.apsb.2025.10.040 | ✓ | | fine |

Derived columns are self-consistent (`num_protonatable_nitrogens` = 4 on the tetra-amines and 1 on the benchmarks and E1, correctly excluding amide N; `has_ester` = 1 only on the three benchmarks; ALC-0315's `ExactMolWt` 765.721 = the 4-carbon head, not the drawn 751.705).

### 1.6 LipidDB did not leak

`build_dataset.py`'s `MERGE` takes exactly one file from this folder. `XD_2026` is **12 rows in `data_vivo.csv`, 0 in `data_vitro.csv`**. I canonicalised all 485 LipidDB rows plus the 720 MTA rows (986 unique structures) and intersected with the 12 emitted canonical SMILES: **one hit, SM-102**, which is the same molecule legitimately appearing in both. Per the brief I did not attempt to settle the LipidDB publication count.

---

## 2. ERRORS FOUND

### E-1 — `stage07_emit.py` still asserts the retracted "299 publications"

**File:** `curation/xd_2026_lipidai/stage07_emit.py`, line 12 (module docstring):
> `LipidDB is a compilation across 299 publications, not an experiment`

Round 1 retracted that number explicitly: `ERROR_CLASSES.md` records *"299 against a true 47" was a reviewer claim propagated unchecked*, and the folder's own `06_lipiddb_qc.json` has `distinct_source_publications_stated_by_paper: 47`, `..._after_normalisation: 56`, `..._letters_only: 50`. The repairs reached `CHANGELOG.md`, `CLAUDE.md` and `ERROR_CLASSES.md` but not this file — and this is the file that emits the rows, so it is the one a future reader opens first. (`stage06_lipiddb.py` line 41 also contains "299", but correctly, as a description of the *wrong* naive `doi:` cut it argues against; leave that one alone.)

**Correct value:** no single count; say "the 47 studies §2.1.1 names" as the folder's own `README.md` already does.
**Severity: cosmetic.** **Confidence: certain.**

### E-2 — `comp()` and `il()` resolve components by name from the live corpus with `.iloc[0]` — error class 30, unguarded

**File:** `curation/xd_2026_lipidai/stage07_emit.py`, lines 53–66.

```python
def comp(col_name, col_smiles, name):
    s = corpus.loc[corpus[col_name]==name, col_smiles].dropna()
    assert len(s), name
    return s.iloc[0]
...
def il(name):
    s = corpus.loc[corpus.IL_name==name,'IL_SMILES'].dropna()
    if not len(s):
        s = pd.read_csv(.../'data_vitro.csv')...
    return Chem.MolToSmiles(Chem.MolFromSmiles(s.iloc[0]))
```

This is the exact construct round 1 wrote up as error class 30 / F-014 and fixed in `curation/build_new_curation.py`, whose `load_corpus_components` now returns an `ambiguity` set and whose `smi()` **refuses** rather than picking. The remedy was scoped to the file where the bug was found. Here, six names (`DSPC`, `Cholesterol`, `DMG-PEG2000`, `ALC-0315`, `SM-102`, `DLin-MC3-DMA`) are resolved by first-row-wins over a corpus written in alphabetical experiment order, with no ambiguity check.

Two aggravating details and one mitigating one:
- **`DLin-MC3-DMA` is the specific name that was ambiguous in round 1** (2,190 ester rows + 203 `AC_2025` dioxolane rows, with `AC_2025` sorting first). It is unambiguous *today* only because the `AC_2025` name repair landed — precisely the "fixed by luck of row order, not by the guard" situation `build_new_curation.py`'s docstring warns about.
- The lookup is **circular**: `stage07_emit.py` reads `data/data_vivo.csv`, which `build_dataset.py` produces from `stage07_emit.py`'s own output. A regeneration therefore depends on the corpus state at the time it is run, and nothing records which state that was.
- Mitigating: I checked every one of the six names today and each resolves to exactly **1** distinct SMILES in both pools, matching `reagent_registry.json`. **No emitted value or structure is wrong.**

**Correct repair:** have `comp()`/`il()` raise on >1 distinct canonical SMILES for a name, and take the three benchmark structures from `provenance/reagent_registry.json` directly rather than from the corpus.
**Severity: metadata / latent.** **Confidence: certain on the mechanism, certain that it is currently benign.**

### E-3 — `structure_assumed = 1` on E1–E9 is outside the flag's own definition, and contradicts F-030 from this same round

`schema_quality.json` defines `structure_assumed` as *"1 where **a specific structural feature** was assumed rather than read from the paper … Pairs with curation_note, which **must say WHICH feature**."*

For E1–E9 **no specific feature was assumed** — the whole structure was read from Fig 5A, correctly, as I have now confirmed atom by atom. The paired note says *"Structure recovered only by geometric trace of the Fig 5A drawing … there is no external anchor against which to check it"*, which names no feature because there is none to name. The flag is carrying "this structure has no external anchor", which is a different proposition.

What makes this worth raising rather than nitpicking is that **F-030, applied in this same round, declined to set `structure_assumed` on ALC-0315 on exactly this reasoning**: *"Nothing was assumed; it was overridden … Stretching `structure_assumed` to cover it would make the flag mean two things."* The same test applied to E1–E9 says the flag is being stretched there too. One of the two decisions should move.

**I do not recommend clearing the flag.** It is currently the only machine-readable marker that these nine rows have no anchor, and losing that is worse than the definitional drift. The clean repair is to **widen the definition in `schema_quality.json`** (e.g. "…or where the entire structure rests on a single reading with no external anchor") and say so in the `meaning` string, or add a distinct term and migrate these nine. Either way the decision should be recorded once, because C13 enforces the vocabulary and not the semantics, so nothing can catch this drifting again.

**Severity: metadata, 9 vivo rows.** **Confidence: certain on the definitional mismatch; the repair is a judgement call.**

### E-4 — `stage03_structures.py`'s provenance claim is not executed (error class 24, on the enumerator rather than the anchor)

Its docstring: *"E1-E9 + benchmarks read from Fig 5A by geometric vertex tracing (stage02). **Every chain length is counted from the stage02 polylines, not by eye.**"*

The twelve SMILES are **typed literals** in the `DRAWN` dict. The script imports `pandas` and `rdkit` only; **no code in this folder reads `02_*.json`**, and nothing asserts that a typed chain length equals the traced polyline's vertex count. The trace was genuinely run — the `02_*.json` files exist and the vertex counts are right — but a human transcribed them, and the transcription is the step the docstring claims is machine-made. This is class 24 one slot over from where it was first found.

The cost here turned out to be zero: I verified all twelve independently. But the claim is the thing that would stop a future reader from re-checking.

**Correct repair:** either weaken the docstring to "chain lengths transcribed from the stage02 polylines", or add an `assert` that loads each `02_*.json` and compares its vertex count per path against the typed SMILES.
**Severity: cosmetic.** **Confidence: certain.**

### E-5 — two errors in the paper, neither ours, worth one line in the folder README

- Methods 2.3.1: *"**As shown in Table 1**, Lipid E 1-9 and mRNA encoding luciferase were provided by WestGene…"* — **Table 1 contains none of E1–E9.** It is the 26-row LipidDB sample table (serial / structure / SMILES / luminescence / class / ref), which I read on page 4. The cross-reference is wrong, and a future curator chasing it for an E-lipid identifier will waste the trip. Worth recording precisely *because* it looks like a lead to an anchor and is not.
- Methods 2.3.1 says the three benchmarks were "synthesized in our laboratory"; the Results (§3.4) call them "commercially accessible ionizable lipids". Internally inconsistent; immaterial to us.

**Severity: cosmetic.** **Confidence: certain.**

### Retracted finding — the packet's `flag_columns` block

When I was handed `packet_XD_2026.md` its `flag_columns` block listed **only `confidence`**, while the corpus slice has four populated quality columns. I was about to report this as a fresh instance of F-029. **On re-reading the file at the end of the session it has been regenerated and now shows all four.** So I was briefed from the pre-regeneration copy. Recording it only so you know the packet I worked from differed from the one in the repo — and noting that none of my conclusions came from the packet.

---

## 3. UNVERIFIABLE

- **`Model_target` = `whole_body`.** I could not break the tie and I do not think it can be broken from the paper. Methods 2.3.4 describes *both* whole-animal imaging at 6 h *and* ex-vivo heart/liver/spleen/lung/kidney on the same instrument; Fig 5B shows both rows; the Fig 5C caption says only "Luciferase expression histogram"; and the unit is identical for either ROI. I checked the one thing that could have discriminated — whether the Fig 5C rank order matches the whole-body panel better than the ex-vivo liver panel — and both orderings agree with Fig 5C, including the three lipids (E6, E7, E9) with no visible signal in either. `whole_body` remains the defensible reading for one number per lipid. **Would be settled by:** the authors' source data or their Living Image ROI definitions.
- **The identity of E1–E9.** I went looking for an anchor the earlier passes may not have: the released repo contains 6 Python scripts, `LipidDB-Ionizable lipid Library.xlsx` (Sheet1 485 / sheet1+tail 366 / sheet1-tail 361), `descriptors_1107.csv` (1,107 rows, descriptor columns only, plus an `inject_way` feature), `Model_Refinement_1.xlsx`, eight base-model contrast workbooks and a `.pkl`. **There is no SMILES, name, formula or identifier for any of the twelve validation lipids anywhere in it**, and I grepped the scripts for SMILES/ALC/SM-102/MC3/validation. The "no external anchor of any kind" statement is now confirmed against the released code, not only against the PDF. **Would be settled by:** WestGene Biopharmaceutical, or a later paper from this group naming E1–E9.
- **Stereochemistry.** E2, E3, E6 and E8 each carry **four** undefined carbinol stereocentres; E1 carries two undefined 2-hexyldecyl branch centres. Fig 5A draws all of them flat, with no wedges anywhere in the panel. The emitted SMILES are correctly stereo-unspecified. Genuinely unknowable; worth knowing if anyone uses stereo-aware fingerprints on this slice.
- **What the bar represents.** Three replicate discs per bar; the paper states no *n* and no error-bar convention, so whether the bar is an arithmetic or geometric mean of the three is unrecoverable. On a log axis the difference is real but small relative to the between-lipid spread.
- **LipidDB publication count** — not pursued, per the brief.

---

## 4. SEVERITY

| # | finding | severity |
|---|---|---|
| E-1 | "299 publications" in `stage07_emit.py` | **cosmetic** (documentation) |
| E-2 | unguarded corpus-lookup-by-name, class 30 | **metadata / latent** — changes nothing today; could change a structure on any future regeneration |
| E-3 | `structure_assumed` outside its definition on 9 rows | **metadata** (9 vivo rows) |
| E-4 | unexecuted provenance claim in `stage03_structures.py` | **cosmetic** |
| E-5 | the paper's own Table 1 cross-reference and the commercial/synthesised contradiction | **cosmetic** |

**No finding changes a structure. No finding changes a value.** Every emitted structure and every emitted `Experiment_value` in `XD_2026` reproduced independently.

---

## 5. CONFIDENCE

| # | confidence | basis |
|---|---|---|
| **ALC-0315 override upheld** | **certain** | head geometry measured to ±0.3 px; three internal controls in the same panel; as-drawn InChIKey reproduces the stored anchor; no authorial structure file exists to contradict the name |
| **9/9 E structures correct** | **certain** | full-InChIKey match from an independently constructed SMILES set, after per-chain geometric measurement at native resolution |
| **12/12 values correct** | **certain** | tick calibration rms 0.0020 decades from 35 ticks, flank tops agreeing L/R to 0 px on 12/12, rms Δz 0.0027 against the stored values |
| E-1 | certain | the retraction is in `ERROR_CLASSES.md` and the number is in the file |
| E-2 | certain on mechanism; certain it is benign today | read the code; measured the name→structure multiplicity in both pools |
| E-3 | certain on the mismatch; the repair is judgement | the schema text and F-030's own reasoning |
| E-4 | certain | the script's imports and the absence of any read of `02_*.json` |
| E-5 | certain | read Table 1 on page 4 and both sentences in full |

---

## Suggested disposition

1. Fix E-1 and E-4 in place — one-line docstring edits, no regeneration needed.
2. Decide E-3 once, for the whole corpus, and record it in `schema_quality.json` rather than in a findings file. It is the kind of thing that will be re-litigated at the next curation otherwise.
3. Treat E-2 as the real one. It is not wrong today, but it is a known class sitting unguarded in an emitter, and the thing that makes it safe (the `AC_2025` name repair) is unrelated to this folder and could be undone by any future merge that reuses one of those six names. Porting `build_new_curation.py`'s ambiguity guard is a small change with a permanent payoff, and the three benchmarks should come from `reagent_registry.json`, not from the corpus the builder writes.
4. Record in the folder README that the repo contains no structures for the twelve, so the next reviewer does not re-search it — and that Methods' "as shown in Table 1" is a dead end.
