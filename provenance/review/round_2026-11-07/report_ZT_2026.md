# Round-2 Wave-1 review — `ZT_2026` (Tian et al., Adv Mater 2026, PMID 42634137)

**Provenance.** Returned as text by a read-only reviewer and saved verbatim.
**Its findings are claims.** My adjudication is in `FINDINGS.md`. **E1 is
confirmed** — I reproduced the byte-identity myself before writing this file.

---

**Round 2, read-only. Nothing in the repository was modified.**

**Summary.** I re-derived the `ZT_2026` slice from the PDFs by my own route and found **seven errors plus one packet bug — none of which changes a structure or an endpoint value**. The corpus content is in good shape: I independently re-digitized all 486 Fig 2c cells at 600 dpi with my own panel detection, my own grid and my own legend swatches and got **0/486 disagreements**; I re-read all 51 Fig 1d printed numerals by eye and got **51/51 agreement**; and the C6O2B2 structure, the six amine heads, both alkyls and all three cholesterol tails reproduce against Fig 1b/2a and the SI NMR. **The most important finding is E1: the C10 anchor is still not independent of the enumeration.** `make_reference.py`, written in round 1 specifically to replace an anchor that was a byte copy of `03_enumerated.csv` (F-013), hard-codes an `AMINE_HALF` template that I show is *character-for-character* the enumerator's own canonical SMILES with the sterol substring replaced by `%` and the ring-closure digit renumbered 1→5 — while its docstring says the amine half was "typed by hand from Fig 2a". Only the sterol half (the half F-013 was about) is genuinely independent. The anchor therefore cannot detect the connectivity and regiochemistry errors C10 exists for (classes 14/15); it is a regression lock, not a verification. This is error class 24 recurring inside the file written to retire it. I separately confirmed, by reading Fig 2a myself, that the C6O2B2 connectivity is in fact **correct** — so this is a verification-integrity finding, not a structural one.

---

## 1. CONFIRMED CORRECT

### 1.1 The four repairs I was asked to re-verify in one pass each — all four are in the emitted data

| finding | check I ran | result |
|---|---|---|
| **F-013** (537 rows carried epi-cholesterol) | RDKit ester-cleavage reaction over **all 51 unique `ZT_2026` `IL_SMILES`**; compared the liberated C27H46O sterol's **full** InChIKey against cholesterol typed by me from the name *cholest-5-en-3β-ol* | **51/51 → `HVYWMOMLDIMFJA-DPAQBDIFSA-N`** (3β). The 3α epimer key `…-VEIPTCAHSA-N` appears **0** times. |
| **F-011** (24 assumed structures relabelled `high`) | `data/data_vivo.csv` | 24 rows `confidence = low`, `structure_assumed = 1`, note *"ASSUMED: no ring N and the paper draws no structure"*. **I verified the note is true of all 24**: C1–C4 are `CNCCN`, `CN(C)CCN`, `NCCN(CCN)CCN`, `CNCCN(CCNC)CCNC` — no ring in any of them (confirmed against Fig 1b at 600 dpi) — and Fig 2a draws only C6O2B2. |
| **F-014** (`MC3-CB*` built on the KC2 dioxolane) | substructure search on the three `MC3-CB*` rows | **no 1,3-dioxolane**; the core is the linoleyl **ester** `…CH(C18:2)2-O-C(=O)CCC[N+]…`, C73H128NO4⁺. Correct. |
| **F-025** (applied today) | `data/data_vivo.csv` | `Lipid 14-CB1/2/3` → **`low`**, `TT3-CB*`/`OF-2-CB*`/`MC3-CB*`/`SM-102-CB*` → **`medium`**, with **two** differentiated notes (the `low` one names the 2-way amine choice and cites F-015). Applied as the validator recommended, not as originally proposed. |
| **F-030** (OF-2 tails C17→C18) | rebuilt the corrected core SMILES and compared to `LL_2025`'s independently curated `OF-02` | **C84H152N4O6 / `UMSCHHWDFWCONK-RGLFHLPNSA-N`, full InChIKey identical including all stereo.** The live `OF-2-CB*` rows carry C18 tails. |
| F-024's 6 note-less rows | `C5O1B*`/`C6O1B*` | now carry a note: *"Tail reagent is isodecyl acrylate (branched terminus), a commercial branched isomer mixture…"*. Closed. |

### 1.2 Fig 2c — fully independent re-digitization, 486/486

Method, sharing nothing with `digitize_fig2c.py`: rendered page 5 at **600 dpi** (the curation worked at 300); located the six panels by thresholding on colour-swatch distance and taking smoothed row/column bands; derived each panel's 9×9 grid by equal division of its own measured extent; sampled a 19×19 px median at each cell centre; matched to **four swatches I sampled myself** from the legend bar in the same raster (`[165,212,226] [238,240,248] [238,210,208] [225,167,169]`, reproducing the stored values to ≤2 counts).

- Bin histogram **17 / 138 / 238 / 93** — identical to `04b_fig2c_486.csv`.
- **0 of 486 cells disagree** on `bin_index`, and `panel_row`/`panel_col` agree on all 486.
- Max colour distance **5.83**, minimum runner-up margin **40.76** — no ambiguous cell.
- Marginals of my own `>1e6` set: C6O2B2 **41/21/31**, DMG-PEG **35/29/29**, DODAP **8/31/54**, helper **DOPC 3 / DMPC 17 / DOPG 16 / DOPE 9 / EggSM 13 / DSPC 35**. All match Figs 2d–g exactly.
- **Axis orientation verified directly from the figure**, not inherited: I read the DOPC panel's furniture at 600 dpi — left axis `C6O2B2 35/50/65` (outer rows), right axis `DODAP 10/30/50` (inner rows), top `10/30/50` under the panel title `DOPC (A)` (outer columns), bottom `DMG-PEG 1 2 4` (inner columns). That is exactly `r//3, r%3, c//3, c%3`.
- **Legend ordering verified**: blue `<10⁴`, pale `10⁴–10⁵`, pink `10⁵–10⁶`, red `>10⁶`.

### 1.3 The 486-point DoE grid reproduces exactly

Recovering raw levels from the stored mol-percents: all **486** rows map to a grid point, **0** unmatched, **0** percentage-degenerate grid points, and each of the six helpers carries the **full 81-point 3×3×3×3 design with no repeats**. Stored ratios sum to 100 ± 0.01 once `fifthcomponent_molratio` (DODAP) is included.

### 1.4 Fig 1d — all 51 printed values re-read by eye at 600 dpi

I located the heatmap band programmatically (y 4698–5030, x 866–3880, 17 columns), cropped each row in overlapping chunks, and read the numerals from the rendered images. **All 51 agree with `04c_fig1d_51.csv` exactly**, including the three I would most expect to be misread (`C2O1B1 = 1.7×10³`, `SM-102-CB2 = 7.0×10³`, `C6O2B2 = 3.0×10⁶`). Column header order (MC3-C, SM-102-C, Lipid 14-C, OF-2-C, TT3-C, C1O1…C6O2) confirmed from the figure.

### 1.5 `Experiment_value` — re-derived, both pools

- **vitro:** population z (ddof=0) of the **bin ordinal** 1–4 with counts 17/138/238/93 gives −2.394584 / −1.091372 / +0.211839 / +1.515050 — matching the four stored values to 1e-15. The ordinal was used, **not** `rlu_midpoint`; the four values are exactly equally spaced, as the convention requires for a discretized source.
- **vivo:** population z of `log10(value_printed)` reproduces all 51 to **2.2e-16**.
- `Experiment_value_within_model` equals `Experiment_value` on all 537 rows (single `Model_type` each).

### 1.6 Four printed statistics the curation never used, all reproduced from the extracted 51 values

This is the strongest independent evidence for the Fig 1d read, because it comes from four other panels.

| source | paper says | from our 51 values |
|---|---|---|
| **Fig 1e** (bars + replicate dots) | — | **Fig 1e's bars are the arithmetic mean of the three B-tail values, and its dots are the individual Fig 1d numbers.** Calibrating only from the y-axis **tick marks** (y = 3800.5 / 3906 / 4010.5 / 4115 px → **104.8 px per 10⁶**, three spacings within 1 px), the four tall bars measure 1.79 / 0.37 / 1.73 / 2.09 ×10⁶ against Fig 1d means of 1.767 / 0.360 / 1.800 / 2.100 ×10⁶ — **mean error 2.4 %, no fitted parameter**. The dot patterns are diagnostic: C6O2's three dots sit at ~3.0, ~2.2, ~1.1 ×10⁶ (Fig 1d B2/B3/B1) and C5O1's at 2.7, 1.3, 1.3 ×10⁶ (B2/B1/B3). |
| **Fig 1f** (`**** p < 0.0001`, piperazine vs not) | piperazine ≫ no-piperazine | 12 vs 39 lipids; medians 1.44×10⁶ vs 6.3×10⁴; **t on log10 p = 4.1e-12**, Mann-Whitney p = 1.9e-7. Dot counts in the two violins are consistent with 12 / 39. |
| **Fig S8a/b** (`ns`) | O1 vs O2 ns; B1/B2/B3 ns | **p = 0.982** and **p = 0.585** (36 synthesized-core lipids). |
| **Fig S9a/b/c** (hit rate >10⁶ p/s) | piperazine higher; C6 > C5 on both expression and hit rate | piperazine **9/12 = 75 %** vs non-piperazine **0/39 = 0 %**; C6 **6/6** vs C5 **3/6**; C6 mean 10^6.27 vs C5 10^5.87, p = 0.052. |

### 1.7 Structures re-derived from the figures and the SI NMR

- **C6O2B2 against Fig 2a (600 dpi):** I traced the drawing independently. Far-left N bears 2 dodecyl-acrylate arms → –CH₂CH₂– to the plain piperazine N; the **N⁺** is the other ring N, carrying the cholesteryl chain and –CH₂CH₂–N(1 arm)–CH₂CH₂–N(2 arms). **That is exactly the stored connectivity**, including which of the two inequivalent C6 piperazine nitrogens is quaternised. Five arms, five nitrogens, C3–O drawn as a bold wedge (3β).
- **Linker length from Fig 2a:** N⁺–C–C–C–C–C(=O)–O = **5-carbon**, matching B2.
- **B1/B2/B3 from Fig 1b:** drawn in full as Br–(CH₂)ₙ₋₁–C(=O)–O–cholesteryl with a bold 3β wedge. Chain vertices counted on the rendered image: **B1 = 3 C, B2 = 5 C, B3 = 7 C** — matching the text ("3-, 5-, or 7-carbon linkers") and the corpus (measured on all 51 SMILES: B1 → 3, B2 → 5, B3 → 7, no exceptions).
- **Amine heads from Fig 1b:** C1 MeNH-CH₂CH₂-NH₂ (3 N–H), C2 Me₂N-CH₂CH₂-NH₂ (2), C3 tris(2-aminoethyl)amine (6), C4 tris(2-(methylamino)ethyl)amine (3), C5 1,4-bis(3-**aminopropyl**)piperazine (4), C6 H₂N-CH₂CH₂-pip-CH₂CH₂-NH-CH₂CH₂-NH₂ (5). All six match `HEADS` including the propyl-vs-ethyl distinction between C5 and C6.
- **Alkyls from Fig 1b:** O1 is 8-methylnonyl acrylate (branch carbon at 8.1 bond-pitches from the ester O, isopropyl terminus) = isodecyl, C10; O2 is linear C12. Matches.
- **NMR, re-read from the SI text layer:** C6O2 δ4.04 (t, **10H**) + δ0.87 (t, **15H**) ⇒ 5 arms; **C5O2 δ4.04 (t, 8H) + δ0.95–0.75 (12H) ⇒ 4 arms** (an arm-count anchor for C5 that `07_qc.json` does not mention); C6O1 10H ⇒ 5 arms. B1's δ3.58 (t, BrCH₂) / δ2.89 (t, CH₂CO) with nothing between is diagnostic of 3-bromopropanoate; the aliphatic envelope grows by ~2 CH₂ from B1→B2→B3.
- **The 73.9 ppm ¹³C claim is TRUE but could not be checked from the text.** `enumerate.py` cites "an isolated peak at 73.9 ppm" as 3β evidence. **The authors' own printed ¹³C list has no such peak** — it jumps 76.74 → 64.94 — so a text-layer check would have called the claim false. I measured the spectrum trace instead: calibrating on the 23 axis ticks (220→0 ppm, 9.355 px/ppm), there **is** a clean isolated peak at **73.76–73.97 ppm, height 224 px against a 151–162 px baseline**, and **nothing above baseline anywhere in 70–72 ppm** where a 3α ester would sit. The docstring is right and the paper's peak list is incomplete.

### 1.8 Other things I checked and found correct

- Helper-lipid SMILES, all six, by InChIKey: DMPC `CITHEXJVPOWHKC-UUWRZZSWSA-N`, DOPC `SNKAWJBJQDLSFF-NVKMUCNASA-N`, DOPE `MWRBNPKJOOWZPW-NYVOMTAGSA-N`, DOPG `DSNRWDQKZIEDDB-SQYFZQSCSA-M`, DSPC `NRJAVPSFFCBXDT-HUESYALOSA-N`, EggSM `RWKUXQNLWDTSLO-GWQJGLRPSA-N` — each is the right molecule for its name (but see E7 on EggSM).
- All 36 enumerated SMILES in `03_enumerated.csv` canonicalize identically to their corpus rows, and every `confidence` transfers unchanged; the vitro and vivo C6O2B2 strings are byte-identical.
- Top-14: all 14 Table S1 / Fig 2j formulations are present in the 486 and **all 14 land in my own `>10⁶` bin**. Fig 2j rows 1–12 re-read from the figure match `04e_top14.csv` exactly, including the lead DSPC-3 = 35/50/30/2.
- `make_reference.py` contains **no** `RWMol`, `AddBond`, `InsertMol` or `RemoveAtom`, takes the sterol from `reagent_registry.json`, and re-roots it with `MolToSmiles(rootedAtAtom=)` with an InChIKey assertion on the re-rooted string. That part of the independence claim is genuine and well built. (The amine half is not — E1.)
- In-vivo metadata against the Methods: C57BL/6 ✓, intravenous ✓, 0.1 mg/kg Luc mRNA ✓, 6 h ✓, `Model_target = brain` ✓ (Fig 1d is explicitly whole-brain flux), `CHL_*` deliberately empty ✓ (the paper's central claim is no free cholesterol).
- The 93-vs-91 discrepancy **is** disclosed, in `07_qc.json` and `PROVENANCE.md`, with the correct argument (the paper's own four Fig 2d–g marginals each sum to 93, and 238+93 = 240+91 = 331).
- A single `Formulation_ID` across the 51 vivo rows is the corpus-wide convention (FG_2026, LL_2026, MO_2026, YX_2026 all do the same), not an anomaly.

---

## 2. ERRORS FOUND

### E1 — The C10 anchor's amine half is the enumeration's own output, not a hand reading of Fig 2a
**File:** `curation/siegwart_42634137/make_reference.py`, `AMINE_HALF` → propagates into `reference_structures.csv`'s `note` and into the C10 status in the packet.

**What is wrong.** The docstring and the emitted note both state the amine half was *"typed by hand from Fig 2a (p.5)"* and that the anchor is *"Independently derived"*. It is not. I took `03_enumerated.csv`'s C6O2B2 string, replaced the sterol substring with `%`, and renumbered the ring-closure digit 1→5. The result is **character-for-character identical** to `AMINE_HALF`:

```
enumeration, sterol → '%', ring digit 1 → 5:
CCCCCCCCCCCCOC(=O)CCN(CCC(=O)OCCCCCCCCCCCC)CCN(CCC(=O)OCCCCCCCCCCCC)CC[N+]5(CCCCC(=O)O%)CCN(CCN(CCC(=O)OCCCCCCCCCCCC)CCC(=O)OCCCCCCCCCCCC)CC5
AMINE_HALF in make_reference.py:
CCCCCCCCCCCCOC(=O)CCN(CCC(=O)OCCCCCCCCCCCC)CCN(CCC(=O)OCCCCCCCCCCCC)CC[N+]5(CCCCC(=O)O%)CCN(CCN(CCC(=O)OCCCCCCCCCCCC)CCC(=O)OCCCCCCCCCCCC)CC5
IDENTICAL: True
```

A hand-typed SMILES cannot coincide with RDKit's canonical traversal — which happens to start at a dodecyl chain terminus and places every branch in canonical order — in 157 characters. The `%` placeholder sits exactly where the sterol's canonical substring begins, and the docstring's own explanation for choosing ring digit 5 ("so it cannot collide with the 1–4 the registry's sterol string uses") is the edit you make after pasting.

**Why it matters.** C10 exists to test connectivity and regiochemistry against external ground truth (error classes 14 and 15). For those classes this anchor can only ever agree, because it *is* the enumeration. Only the sterol half — class 31, the one F-013 was about — is genuinely checked. The anchor is a useful regression lock against *future* changes to `enumerate.py`, but it is not retrospective verification, and the note says it is.

**Correct value / how to verify.** The structure itself is right — I confirmed C6O2B2's connectivity, arm placement and N⁺ position directly from Fig 2a at 600 dpi (§1.7). The repair is to the provenance, not the molecule: either re-type the amine half in a different, non-canonical traversal (e.g. head-first from the piperazine) and let RDKit canonicalize it before comparison, or state plainly in the note that the amine half is transcribed from the enumeration and that the anchor therefore tests the sterol only. To verify my claim, run the four-line transformation above.

### E2 — Three of `digitize_fig2c.py`'s four marginal checks are broken, record themselves as failing, and are reported as passing
**Files:** `curation/siegwart_42634137/digitize_fig2c.py` lines 303–312 and `07b_fig2c_qc.json` → contradicted by `07_qc.json`.

`marg(key)` groups the outstanding cells by `r["C6O2B2_molpct"]` / `"DMGPEG_molpct"` / `"DODAP_molpct"` and then does `got.get(k, 0)` against keys `35/50/65`, `1/2/4`, `10/30/50`. The 2026-10-28 `column_naming_fix` renamed the ratio columns to `*_molratio` and put genuine percentages in `*_molpct`, so those three lookups now miss **every** key and fall through to the literal `0`. The emitted QC file records exactly that:

```json
"C4_fig2d_by_C6O2B2": {"digitized": {"35":0,"50":0,"65":0}, "reported": {"35":41,"50":21,"65":31},
                       "delta": {"35":-41,...}, "rank_order_matches": false}
```
`C7_fig2g_by_helper` passes only because `helper_lipid` was not renamed.

Meanwhile `07_qc.json` states *"ALL FOUR EXACT … 15 of 15 bins match the digitized counts with zero deviation."*

**The substance still holds** — I recomputed the marginals from `04b` using `*_molratio` and got 41/21/31, 35/29/29, 8/31/54 and the six helper counts exactly, and `digitize_bars.py`/`07d_bars_qc.json` already do it that way. So nothing in the data is wrong. What is wrong is that the panel-orientation argument quoted in `PROVENANCE.md` check 3 now rests on a check that no longer runs, and a re-run today would emit a QC file declaring three failures beside a summary saying they pass.

This is **error class 28 verbatim** (`else 0` / `.get(k, 0)` turning a detector miss into a number) plus **error class 35's general form** (a rename invalidating everything keyed on the old name). Fix: key on `*_molratio`, and make the lookup `got[k]` with a `KeyError`, not `got.get(k, 0)`.

### E3 — The in-vitro dose/n contradiction reaches no row (error class 34, in this curation, after round 1)
**Files:** `curation/build_new_curation.py` ~line 420 (`"Dose_ug_nucleicacid": 0.02`) and `curation/new_curation_qc.json` → `paper_internal_contradictions_recorded.in_vitro_dose_and_n`.

Fig 2c's caption says **20 ng mRNA per well, n = 3**; Methods §4.8 says **15 ng per well, n = 5**, and the curation established (correctly) that §4.8 is the only in-vitro subsection and describes the same wells, so the paper contradicts itself. The curation chose the caption and wrote `UNRESOLVED` into the QC JSON.

All **486** vitro rows carry `Dose_ug_nucleicacid = 0.02`, `confidence = high`, and **blank `curation_note`, blank `structure_assumed`, blank `value_is_bound`**. A consumer of `data_vitro.csv` sees a settled 20 ng dose at the project's highest confidence. The disclosure exists only in a per-folder QC file — *"a flag that lives only in a per-folder QC file is not a flag, it is a memo."* This is the same shape as the open `WT_2026` E5 item this round.

`PROVENANCE.md` compounds it twice: its conditions block states *"20 ng mLuc per well, 24 h, n = 3"* as fact, and its section headed **"Two places the paper disagrees with itself"** lists only the 240/91 count and the Fig 1d max — not this, the third and the only one that reached a stored value.

**Repair:** a `curation_note` on the 486 rows naming both readings and which was taken. No value need change.

### E4 — All six C5 lipids record their quaternization site as "drawn in Fig 2a", and the paper draws no C5 structure
**Files:** `03_enumerated.csv` (`quaternization_site`, `quaternization_site_assumed`) for `C5O1B1-3` and `C5O2B1-3`; produced by `enumerate.py::_quaternization_site`.

The function returns the literal string `"drawn in Fig 2a (piperazine N)"` whenever the chosen nitrogen is in a ring, and `quaternization_site_assumed` is set to `0`. **Fig 2a draws exactly one molecule, C6O2B2.** No C5 conjugate appears in Fig 1b, Fig 2a or anywhere in the 23-page SI (Fig 1b draws the C5 *head*, unquaternized). So six rows assert figure provenance that does not exist, and declare the site not-assumed when it was chosen by heuristic.

**The structures are nevertheless right.** I ran `enumerate.py`'s own `michael()` on the C5 head and checked symmetry with `Chem.CanonicalRankAtoms(breakTies=False)`: the two piperazine nitrogens of the tetra-alkylated C5 core fall in **one** symmetry class (`{80}`), so the choice is forced and cannot be wrong. For C6 the two ring nitrogens are **distinct** (`{98, 99}`) and Fig 2a does settle it — I verified the chosen one carries the three-arm branch, as drawn.

**Repair:** make the note reflect the actual basis — "forced by symmetry: the two piperazine N are equivalent" for C5, "drawn in Fig 2a" for C6 — which is exactly the symmetry-class computation `_site_confidence()` already implements for the commercial cores. The two functions currently answer the same question two different ways.

### E5 — "B1/B2/B3 ester connectivity is inferred and flagged medium" is wrong in both halves
**Files:** `07_qc.json` → `structures.confidence_flags`; `enumerate.py`, the comment above `TAILS`.

1. **It is not inferred.** Fig 1b draws B1, B2 and B3 in full, each as Br–(CH₂)ₙ₋₁–C(=O)–**O–C3 of cholesterol with a bold 3β wedge**, at roughly 4× the scale of Fig 2a's sterol. I rendered and read all three.
2. **No row is flagged medium for it.** `TAILS` is a 2-tuple `(n_link, description)` with no confidence field, and the confidence expression (`enumerate.py:309`) consults only the head key and the alkyl key. The six `C5O2B*`/`C6O2B*` rows are `high` with blank `structure_assumed`; the six `*O1B*` rows are `medium` **because of O1**, not because of the tail.

**Consequence beyond the docs:** `reference_structures.csv` holds one anchor and `make_reference.py`'s note states the sterol stereocentres "come from the registry, not from a fresh read of the figure, so this anchor cannot catch an error in the registry entry itself." That limit is self-imposed — **Fig 1b draws the complete sterol three times**, so three figure-derived anchors are available and would close exactly the gap the note describes.

### E6 — `Dose_ug_nucleicacid = 2.0` on the 51 vivo rows encodes an undocumented body-weight assumption
**File:** `curation/build_new_curation.py` ~line 535.

The paper gives **0.1 mg/kg** and Methods §4.11 gives the mice as **18–20 g**, i.e. **1.8–2.0 µg**. The stored 2.0 µg is the top of that range. The mg/kg → µg conversion and the 20 g assumption appear in no file. Up to an 11 % overstatement on a stored absolute dose. Severity is low but it is a *value* in a column, derived by an undeclared convention. Repair: either store 1.8–2.0 with a note, or keep 2.0 and say in `curation_note` that it assumes a 20 g mouse.

### E7 — Three named catalogue reagents in this slice are outside the C11 registry, and one is a natural mixture stored as a single species

| reagent | rows | column | note |
|---|---|---|---|
| **DODAP** | **537** (every `ZT_2026` row) | `fifthcomponent_name` | C11 *does* scan `fifthcomponent`, so adding the entry costs one line. Stored as C41H77NO4 / `NYDLOCKCVISJKK-JKUUOIIVSA-N`, which is the right molecule; it is simply unchecked. |
| **DMPC** | 81 | `HL_name` | unchecked |
| **EggSM** | 81 | `HL_name` | unchecked **and** a convention |

**EggSM is the sharper half.** Egg sphingomyelin is a natural mixture (~84 % N-palmitoyl plus C18:0, C22:0, C24:1…). The corpus stores the single species SM 16:0 with `confidence = high` and no note, on 81 rows. CLAUDE.md §3 already states the principle for the identical case — PEG-lipids are *deliberately absent* from the registry "because the block length is a distribution, so no single SMILES is right" — and `18PG` (error class 19) is the precedent for the harm. EggSM appears in exactly one experiment, so this is cheap to fix now.

### E8 — Packet bug: `brief.py` drops the fifth component, so every `ZT_2026` row in the packet appears to violate C5
`provenance/review/brief.py:145–151` builds `molratios` from `("IL", "HL", "CHL", "PEG")` only, and `sample_rows`'s `keep` list likewise omits `fifthcomponent_molratio`. Every sampled `ZT_2026` row therefore shows ratios summing to ~79 % with the missing 20.55 % never named. For an experiment whose fifth component is 10–50 mol % of the particle, that is materially misleading and would send a reviewer after a non-existent C5 violation. The same function's C11 block *does* read `fifthcomponent_name`, so the omission is an oversight. One-line fix in both places.

Separately — **not a bug** — the packet's vivo `confidence` tally is stale against the live corpus; the packet was generated before F-025 was applied. The difference is exactly the six `Lipid 14-CB*`/`TT3-CB*` rows.

---

## 3. UNVERIFIABLE

- **91 vs 93 "outstanding" LNPs.** The paper contradicts itself. I attempted to settle it from Fig 2h, the 486-point scatter: calibrating on the shaded band boundaries and counting ink above the 10⁶ line gave **79–83 points** by two area methods, with ±15 % uncertainty from dot overlap. **That cannot discriminate 91 from 93.** The curation's decision to use the discrete bins is correct. Needs the authors' source data.
- **Which tertiary amine is quaternised in C1–C4** (24 rows). Already flagged `low` + `structure_assumed = 1` with the reason on the row. Correctly handled; unknowable without the authors' NMR of a C1–C4 conjugate.
- **The conjugation site on the 15 commercial-core lipids.** The SI has NMR for C6O2, C6O1, C5O2, B1, B2, B3 and C6O2B2 **only** — no spectrum, MS or structure for any `-C` conjugate. F-025's resolution is the right posture.
- **The isodecyl isomer distribution in O1.** Flagged today; nothing further knowable.
- **`n` per group in vitro (3 vs 5).** Not representable in the schema. Needs the authors.
- **Whether `Dose_ug_nucleicacid` should be 1.8 or 2.0** (E6). The Methods give a weight range.

---

## 4. SEVERITY

| # | what it changes |
|---|---|
| **E1** | **nothing in the data** — but it removes the retrospective verification C10 is the project's only source of. Verification integrity; highest-leverage item here. |
| **E2** | **nothing in the data.** Validation integrity: three checks silently return 0 and a summary reports them as passing. |
| **E3** | **metadata only**, 486 rows — a contested dose stored as settled, caveat invisible to any consumer. |
| **E4** | **metadata only**, 6 rows. Structures right (symmetry-forced); recorded evidence false. |
| **E5** | **metadata only** — plus three unexploited C10 anchors. |
| **E6** | **changes a value** in `Dose_ug_nucleicacid` on 51 rows by up to 11 %. No endpoint moves. |
| **E7** | **metadata**, 537 + 81 + 81 rows. C11 coverage gap, and one convention (EggSM) unflagged. |
| **E8** | **cosmetic**, but a reviewer-facing artifact that can manufacture a false finding. |

**No finding changes a structure. No finding changes an `Experiment_value`. Row counts stay 486 vitro / 51 vivo.**

---

## 5. CONFIDENCE

| # | confidence | basis |
|---|---|---|
| **E1** | **certain** | exact string equality after a two-token mechanical transformation, reproducible in four lines. |
| **E2** | **certain** | the failure is written into `07b_fig2c_qc.json` in the repository. |
| **E3** | **certain** | both paper statements read from the PDF; blank `curation_note` read from `data_vitro.csv`. |
| **E4** | **certain** on the false provenance string and on the symmetry equivalence; **certain** the structure is right. |
| **E5** | **certain** — rendered Fig 1b and read all three drawn tails; `TAILS` has no confidence field. |
| **E6** | **certain** on the arithmetic and the absence of any record; the repair is a judgement call. |
| **E7** | **certain** on the registry gap; **probable** that EggSM should carry a mixture flag. |
| **E8** | **certain** — read from `brief.py` and visible in the packet's own JSON. |

---

### One line for the ratchet

E2 and E4 are the same shape as classes 28 and 24 and need no new class. **E1 does:** *a replacement anchor written to cure a copied anchor, which re-copies the half the original bug did not touch.* The detectable signature is cheap — if the anchor string, after any placeholder substitution and ring-renumbering, is a substring match to the enumeration's canonical SMILES, the anchor is a transcription. That test would have fired on this file the day it was written, and it generalises to the six `ng_2025` anchors and the one `LL_2026` anchor still flagged `anchor_may_be_a_copy`.
