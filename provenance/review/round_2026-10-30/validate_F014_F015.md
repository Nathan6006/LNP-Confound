# Validation review — F-014 and F-015

Reviewer: validation pass, read-only. Paper: Tian et al., *Adv Mater* 2026,
PMID 42634137, DOI 10.1002/adma.74762 (`curation/siegwart_42634137/`).
All measurements below were re-derived from the PDFs and from `data/` with
independently written rdkit scripts; nothing was taken by re-running the
project's own build or audit code.

> **State-of-the-repo warning.** Files moved underneath this review.
> `curation/build_new_curation.py` (17:23:50), `curation/new_curation.csv`
> (17:23:56), `curation/new_curation_qc.json` (17:23:56) and both corpus files
> (17:24:28) were rewritten by another session *during* this review, after I had
> already measured the pre-repair state. Both states are reported. F-014 was
> applied mid-review; F-015 was not.

---

## 1. VERDICT F-014 — rewrite the STRUCTURE, keep the name

**The name `MC3-C` is right and the stored structure was wrong.** This is the
opposite half from the `AC_2025` repair, and the paper settles it directly.

### What settles it

**Main text p. 3, Figure 1b, row labelled "Commercially available cores".**
The panel draws five parent cores. The one captioned **"Dlin-MC3-DMA (MC3)"**
is drawn as:

> (CH₃)₂N–CH₂CH₂CH₂–C(=O)–O–CH(R)(R′), where R and R′ are two straight
> dilinoleyl chains, each carrying two *cis* alkene marks in the
> 9,12 / 6,9 pattern.

That is the **acyclic 4-(dimethylamino)butanoate ester of dilinoleylmethanol**
— DLin-MC3-DMA. **There is no ring of any kind in the drawing, and in
particular no 1,3-dioxolane.** Read as an image at 700 dpi from
`page 3, clip (300,278)–(612,345)`. The figure caption does not enumerate the
commercial cores (it lists only "amine head groups (C1–C6), alkyl groups (O1
and O2), and cholesterol tails (B1, B2, B3)"), but the panel itself has a
separate left-hand band labelled "Commercially available cores" containing
SM-102, OF-2, TT3, Dlin-MC3-DMA (MC3) and Lipid 14.

**Figure 1d** (same page) is the 17 × 3 heatmap. Its column labels are
`MC3-C, SM-102-C, Lipid 14-C, OF-2-C, TT3-C, C1O1 … C6O2`; its rows are
`B1/B2/B3`. So the measurement is labelled by the paper as the cholesterol
conjugate **of the molecule Fig 1b names and draws as MC3**. The measurement
and the name travel together; the SMILES was generated afterwards by a script.

**Figure 1e / p. 4** confirms `MC3-C` is a conjugate and not the benchmark:
Fig 1e has separate bars for `MC3` and `SM-102` under a bracket labelled
"Benchmark", and `MC3-C, SM-102-C, Lipid 14-C, OF-2-C, TT3-C` under
"Cholesterol-modified". The text (p. 4) reads: "conventional MC3 LNPs and
SM-102 LNPs **formulated with unmodified ionizable lipids** exhibited minimal
brain delivery."

Because the label is the paper's own and the structure is a derived artifact of
a lookup, rewriting the structure keeps the measurement with the molecule the
paper measured. Rewriting the name would invent a KC2 conjugate the paper never
made — `KC2` / DLin-KC2-DMA appears nowhere in the paper or the SI.

### Measurements — confirmed, with one correction

Re-derived with rdkit over `data/data_vivo.csv` + `data/data_vitro.csv`
(ring count = `CalcNumRings`, dioxolane = SMARTS `[#6]1-[#8]-[#6]-[#8]-[#6]-1`):

| | rows | rings | dioxolane | esters | exact mass |
|---|---|---|---|---|---|
| `MC3-CB1/2/3` **as stored at review start** | 3 | **5** | **1** | 1 | 1082.9838 / 1111.0151 / 1139.0464 |
| `MC3-CB1/2/3` **as stored now (17:24)** | 3 | **4** | **0** | 2 | *identical masses* |
| corpus `DLin-MC3-DMA` | **3,086** | 0 | 0 | 1 | 641.6111 |
| corpus `KC2` | 203 | 1 | 1 | 0 | 641.6111 |

The table in the review packet is confirmed, except the `DLin-MC3-DMA` row
count: I measure **3,086**, not 3,097, scanning `IL/HL/PEG/CHL/fifthcomponent`
name columns in both corpus files. The 3,097 figure comes from
`provenance/reagent_registry.json`, which additionally folds in the alias
`Dlin-MC3-DMA` (`reagent_registry.py:93`, `ALIASES`). Both numbers are right
for their own definition; neither affects the verdict.

Note the trap worked exactly as warned: **the wrong and the right structure have
the same molecular formula and the same exact mass to four decimals** at every
linker length (C73H128NO4⁺ / 1082.9838 etc.). Mass and formula discriminate
nothing here. Ring count and the dioxolane SMARTS do.

### Root cause — verified, and sharper than reported

`curation/build_new_curation.py`:

```
159      comp.setdefault(n, s)            # load_corpus_components()
322      PARENT = {"MC3-C": "DLin-MC3-DMA", "SM-102-C": "SM-102"}
339-341  if core in PARENT and comp.get(PARENT[core]):
             parent_smiles = comp[PARENT[core]]
             conf, how = "medium", f"corpus {PARENT[core]}"
```

The lookup is by name into the live corpus, with `setdefault` — **first row in
file order wins**. That detail is load-bearing and was missing from the reported
story.

At build time the name `DLin-MC3-DMA` was *ambiguous*, not simply wrong. In the
pre-patch snapshot `provenance/superseded/LNPDB_vitro.csv` the name carries
**two different structures**: 2,190 rows of the correct ester, and **203 rows
(the `AC_2025` block) of the dioxolane**. Which one `setdefault` returns depends
only on row order:

* `provenance/superseded/LNPDB_vitro.csv` is in **chronological** experiment
  order (`AA_2008, KL_2010, … AC_2025` at row 17050). Its first
  `DLin-MC3-DMA` row is #4972, `SP_2020`, the **ester**. Running the builder's
  lookup against this snapshot returns the ester — the reported story does
  **not** reproduce here.
* The live `data/data_vitro.csv` is in **alphabetical** experiment order, which
  puts `AC_2025` second (row 497). Pre-patch, that block's `DLin-MC3-DMA` rows
  were the **dioxolane**.

Reconstructing that state exactly — today's alphabetically-sorted corpus with
`AC_2025`'s `KC2` ↔ `DLin-MC3-DMA` name patch un-applied — the builder's lookup
returns the **dioxolane**, and the `max`-index quaternization of it reproduces
the stored `MC3-CB1/CB2/CB3` **bit-for-bit (canonical SMILES identical, all
three)**. Root cause proven.

So the mechanism is: *name lookup, into a corpus that was sorted alphabetically,
where one name transiently covered two structures and the mislabelled block
happened to sort to the front.* The `AC_2025` name patch (CHANGELOG §10,
2026-10-29) is confirmed as the upstream event; the build (`new_curation.csv`,
2026-10-28 23:51) predates it.

### Current state — the repair was applied mid-review, and applied correctly

As of 17:23–17:24 today, `MC3-CB1/2/3` in both `curation/new_curation.csv` and
`data/data_vivo.csv` carry **4 rings, no dioxolane, 2 esters**, and match my
independent rebuild from the real DLin-MC3-DMA parent bit-for-bit. Example:

```
MC3-CB1  CCCCC/C=C\C/C=C\CCCCCCCCC(CCCCCCCC/C=C\C/C=C\CCCCC)OC(=O)CCC[N+](C)(C)
         CCC(=O)O[C@@H]1CC[C@@]2(C)C(=CC[C@H]3[C@@H]4CC[C@H]([C@H](C)CCCC(C)C)
         [C@@]4(C)CC[C@@H]32)C1
```

`Experiment_value` did **not** move: the three z-scores (−0.360215491,
−0.699161915, −0.457890866) are unchanged and still reproduce from Fig 1d's
printed cells 7.0×10⁴ / 4.1×10⁴ / 6.0×10⁴ via `log10` + population z over the
51 values. Structure moved, value stayed — the right shape for this repair.

**Caveat on durability.** The repair was effected only by regenerating from the
now-patched corpus. `PARENT` at line 322 and `setdefault` at line 159 are
unchanged, so the mechanism that produced the error is still in place: the
correctness of these six rows remains a silent function of the corpus's name
table and its row order at the moment of the next rebuild. Nothing in the
builder or in C10/C11 would catch a recurrence (see §3).

**Confidence: certain.**

---

## 2. VERDICT F-015 — CONFIRMED as a code defect; the SITE is UNVERIFIABLE

Both halves of the reviewer's claim are true. A third point, which changes the
recommended repair, is added.

### 2.1 Second copy exists — `curation/build_new_curation.py:199-204`

```python
199    cands = [a for a in mol.GetAtoms()
200             if a.GetSymbol() == "N" and a.GetTotalNumHs() == 0
201             and sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1) == 3]
202    if not cands:
203        return ""
204    t = max(cands, key=lambda a: a.GetIdx()).GetIdx()
```

Line 204 is `max(atom index)`. Still present at 17:23:50, after the F-014
repair.

### 2.2 The docstring claims parity — `curation/build_new_curation.py:164-170`

```python
164  def quaternize(core_smiles, n_link):
165      """Cholesteryl-O2C-(CH2)n- onto the most substituted tertiary N.
166
167      Same construction as the synthesized library (see siegwart_42634137/
168      enumerate.py); used here only for MC3-C and SM-102-C, whose parent SMILES
169      come from the corpus. The paper does not draw these two conjugates, so they
170      carry a medium-confidence flag in the QC.
171      """
```

Two false statements, not one:

* **"Same construction as the synthesized library"** — it is not.
  `enumerate.py:123-166` (`_quaternization_site`) explicitly names
  `max(atom index)` as "an artifact of SMILES writing order" that it replaced.
* **"used here only for MC3-C and SM-102-C"** — also false. Line 327 routes
  `TT3-C`, `Lipid 14-C` and `OF-2-C` through the same `quaternize` via
  `COMMERCIAL_CORES` at line 324. This matters: the docstring asserts that the
  very molecule the finding is about is not handled by this function.

### 2.3 Where the two rules actually diverge — my own run

Both rules re-implemented from scratch and run over all five commercial parents:

| parent | tertiary N | `max(index)` picks | repaired rule picks | same atom | **same product** |
|---|---|---|---|---|---|
| DLin-MC3-DMA | [43] | 43 | 43 | yes | yes |
| KC2 (dioxolane) | [42] | 42 | 42 | yes | yes |
| SM-102 | [19] | 19 | 19 | yes | yes |
| TT3 | [6, 40, 73] | 73 | 6 | **no** | **yes** (C₃-symmetric arms) |
| OF-2 | [8, 53] | 53 | 8 | **no** | **yes** (C₂-symmetric arms) |
| **Lipid 14** | **[1, 11]** | **11 — central trialkyl N** | **1 — terminal N(CH₃)₂** | **no** | **NO** |

`Lipid 14` is the only molecule where the divergence survives into a different
product. Atom 11 is the internal nitrogen bearing the two
`–(CH₂)₈–O–C(=O)–CH(C₆H₁₃)(C₈H₁₇)` arms; atom 1 is the terminal dimethylamino.
Checked against `curation/new_curation.csv`: all three `Lipid 14-CB*` rows match
the `max(index)` build exactly and the repaired-rule build not at all — i.e.
**the stored structures are alkylated at the hindered central nitrogen.**
`SM-102-CB*`, `TT3-CB*` and `OF-2-CB*` are unaffected. **3 rows, exactly as
reported.**

As with F-014, every mass-like descriptor is blind to this: the two candidate
`Lipid 14-CB1` structures agree on **all nine** descriptors I computed —
formula C86H159N2O8⁺, exact mass 1348.2091, heavy atoms 96, rings 4, TPSA
108.44, LogP 23.474, rotatable bonds 60, Fsp3 0.9302, HBA 9. They are
constitutional isomers that differ only in which nitrogen carries the charge.

### 2.4 Which nitrogen does the paper draw? — it does not

I went through the figures. **The paper and SI draw no cholesterol conjugate of
any commercial core.**

* **Fig 1b (p. 3)** draws Lipid 14 as the **parent**:
  (CH₃)₂N–(CH₂)₃–C(=O)O–CH₂CH₂–N[(CH₂)₈–O–C(=O)–CH(C₆H₁₃)(C₈H₁₇)]₂. Both
  nitrogens are neutral and unalkylated. Read at 700 dpi. (This does confirm
  the hard-coded `COMMERCIAL_CORES["Lipid 14"]` parent SMILES is right.)
* **Fig 2a (p. 5)** draws the only quaternized structure in the paper,
  **C6O2B2**, with the `N⁺` unambiguously on a **piperazine ring nitrogen**,
  cholesteryl-ester linker hanging off it. This vindicates `enumerate.py`'s
  ring-preference rule for C5/C6 — and says nothing about Lipid 14, which has
  no ring nitrogen.
* **SI Figures S1–S7** are ¹H/¹³C NMR of **C6O2, C6O1, C5O2, B1, B2, B3 and
  C6O2B2 only**. **S13** is the structure of **C6O2**. There is no NMR, no MS
  and no drawing for `MC3-C`, `SM-102-C`, `Lipid 14-C`, `OF-2-C` or `TT3-C`
  anywhere in the 23-page SI.
* **Methods, p. 11** gives only the reaction, not the regiochemistry:
  "cholesterol tails (1.0 equiv) and commercially available cores or synthesized
  cores (1.0 equiv) were dissolved in anhydrous acetonitrile … heated at 95 °C
  for 24 h."

So: **the quaternization site of `Lipid 14-CB*` is unverifiable from the
paper.**

### 2.5 The important extra finding — porting the repaired rule does NOT fix it

The obvious repair ("make the second copy match `enumerate.py`") is **not a
fix**. For Lipid 14 *both* rules fall through to the atom-index tiebreak,
because neither nitrogen is in a ring and `beta_ester_arms` returns 0 for both.
The repaired rule only lands on the terminal N(CH₃)₂ because of how the SMILES
happens to be typed. I wrote the identical molecule tail-first and re-ran both:

| Lipid 14 SMILES written… | `max(index)` picks | repaired rule picks |
|---|---|---|
| head-first (as in the builder) | central trialkyl N | terminal N(CH₃)₂ |
| **tail-first (same molecule)** | **terminal N(CH₃)₂** | **central trialkyl N** |

The rules simply swap answers. The repaired rule is a SMILES-ordering artifact
here too — the exact defect its own docstring says it was written to eliminate.

The reason is a gap in the detector: `beta_ester_arms`
(`enumerate.py:143-158`) walks `N–C–C–C(=O)` and so cannot see Lipid 14's
central nitrogen arm, which is `N–CH₂CH₂–O–C(=O)R` (β-**acyloxy**, the oxygen
intervenes). Counting β-acyloxy as electron-withdrawing alongside β-ester makes
the rule order-independent — I tested
`[NX3;H0]-[CH2]-[CH2]-[$([CX3](=O)[OX2]),$([OX2][CX3]=O)]`: it picks the
terminal N(CH₃)₂ for Lipid 14 from **either** SMILES spelling, and is a clean
regression (still the piperazine N for the C5O2 and C6O2 cores, matching
Fig 2a; unchanged for MC3, SM-102, TT3, OF-2).

### 2.6 Verdict and recommended repair

**CONFIRMED** on every checkable element: the second copy exists, it is
`max(atom index)`, its docstring falsely claims parity (and falsely denies
covering Lipid 14), and it alkylates the hindered central nitrogen of
`Lipid 14-CB1/CB2/CB3`.

**UNVERIFIABLE** on the chemistry: the paper draws no `-C` conjugate, so no
figure can say which nitrogen is the site.

Recommended deliverable, in order:

1. **Flag the 3 rows**, do not silently re-point them. Whatever is written,
   `Lipid 14-CB*` carries an inferred conjugation site, and the honest record is
   a flag. Drop their QC confidence from `high` to `medium` (see §3).
2. **Fix the docstring at 164-170** regardless — both of its claims are false
   today, and the second one actively misdirects.
3. **De-duplicate the rule.** Two copies of this logic is the real defect; the
   builder should import `_quaternization_site` from `enumerate.py` rather than
   hold a second copy.
4. If a site must be chosen, choose the **terminal dimethylamino**, on the
   order-independent β-acyloxy criterion above, and record it as **inferred**.
   The chemistry supports it (dimethylalkylamine, unhindered, no β-EWG, versus
   a trialkylamine flanked by two C8 chains with a β-acyloxy arm — a
   Menshutkin alkylation at 95 °C would strongly favour the former), but that is
   chemical inference, not evidence from this paper.

**Confidence: certain** on the code defect and on which rows diverge;
**certain** that the site is undrawn; **probable** that the terminal NMe₂ is
the real site, on chemical grounds alone.

---

## 3. OTHER STALE NAME-LOOKUP ROWS

**Complete enumeration.** Two mechanisms in `build_new_curation.py` read the
live corpus by name: `PARENT` (line 322, ionizable-lipid parents) and
`smi(name)` (line 269, `comp.get(name) or SUPPLIED.get(name, "")`, helper
components). Everything else is hard-coded or enumerated from the figure.
Provenance determined by reconstructing the corpus **excluding `ZT_2026`'s own
rows** — the naive check is circular, because four of the helper lipids now
appear in the corpus *only* as the rows this file contributed.

| lipid / component | rows | parent source | verdict |
|---|---|---|---|
| `MC3-CB1/CB2/CB3` | 3 | **corpus name lookup** `DLin-MC3-DMA` | **WAS STALE — F-014. Repaired 17:24; now correct (ester, 4 rings). Mechanism unchanged.** |
| `SM-102-CB1/CB2/CB3` | 3 | **corpus name lookup** `SM-102` | **CLEAN.** Corpus `SM-102` = C44H87NO5, 709.6584, one tertiary N; matches heptadecan-9-yl 8-((2-hydroxyethyl)(6-oxo-6-(undecyloxy)hexyl)amino)octanoate and the Fig 1b drawing. Name unambiguous in the corpus (1 distinct structure, 2,295 rows). Both quaternization rules agree. |
| `DSPC` (HL) | 81 | **corpus name lookup** | **CLEAN.** C44H88NO8P, 789.6248, 18:0/18:0 PC. |
| `DOPE` (HL) | 132 | **corpus name lookup** | **CLEAN.** C41H78NO8P, 743.5465, 18:1/18:1 PE. |
| `DMG-PEG2000` (PEG) | 537 | **corpus name lookup** | **CLEAN as a representative.** C122H242O50, PEG-45 dimyristoylglycerol. PEG block length is a distribution; per CHANGELOG §10 PEG-lipids are deliberately outside the C11 registry. No single SMILES is strictly truthful, but this is the corpus convention and it is applied consistently. |
| `TT3-CB*`, `OF-2-CB*` | 6 | hard-coded `COMMERCIAL_CORES`, read off Fig 1b | not name-lookup. Structures unaffected by the rule divergence (symmetric arms). |
| `Lipid 14-CB*` | 3 | hard-coded `COMMERCIAL_CORES` | not name-lookup. **F-015 applies** — parent correct, site wrong/unverifiable. |
| 36 `CxOyBz` | 36 | `03_enumerated.csv` from Fig 1b blocks | not name-lookup. |
| `DOPC`, `DMPC`, `DOPG`, `EggSM` (HL) | 324 | `SUPPLIED` literals | not name-lookup — **absent from the pre-merge corpus**, so the `comp.get` arm never fired. Structures check out (C44H84NO8P / C36H72NO8P / C42H78O10P⁻ / C39H79N2O6P). |
| `DODAP` (5th) | 537 | `SUPPLIED` literal | not name-lookup. C41H77NO4, 647.5853. |

**`MC3-CB*` was the only stale one.** I also checked every name the builder
looks up for a second ambiguity of the kind that caused F-014: in today's
`data_vitro.csv`, `DLin-MC3-DMA`, `SM-102`, `DSPC`, `DOPC`, `DMPC`, `DOPG`,
`DOPE`, `EggSM`, `DMG-PEG2000` and `Cholesterol` each resolve to **exactly one
distinct structure**. The lookup table is currently unambiguous. That is a
property of this moment, not a guarantee.

### Three further defects found while enumerating

1. **Confidence is laundered at the merge boundary.**
   `03_enumerated.csv` marks the 24 `C1*/C2*/C3*/C4*` rows **`low`**, with
   `quaternization_site = "ASSUMED: no ring N and the paper draws no
   structure"`. `build_new_curation.py:336` overwrites this with a flat
   `"high"` for anything found in the enumerated file. The per-row confidence
   that `enumerate.py` carefully computed is discarded; `new_curation_qc.json`
   reports `low: []`.

2. **`Lipid 14-CB*` and `TT3-CB*` are labelled `high` against the QC's own
   stated rule.** `new_curation_qc.json → confidence_basis` says *medium* means
   "the cholesterol CONJUGATION SITE is inferred (the paper does not draw any -C
   conjugate)". The paper draws **no** `-C` conjugate for any of the five
   commercial cores — I verified this against the full main text and the 23-page
   SI. So all 15 should be at most `medium`; six are `high`.

3. **Two module-docstring claims are false.**
   `build_new_curation.py:50-52` says DOPC/DMPC/DOPG/EggSM and DODAP "are not in
   the corpus" — true at build time for all five, and `new_curation_qc.json →
   supplied_component_SMILES` lists all five, which is correct. But the QC's
   `how` field (line 449) says Fig 1b "DRAWS for all five (MC3, SM-102, TT3,
   Lipid 14, OF-2)" — true of the **parents**, and the sentence is used to
   justify the conjugates' structures, which Fig 1b does not draw.

### Why neither finding was caught

* **C10** (`structure_check.json`) compares an enumerated structure against a
  drawn one. Its only `siegwart_42634137` anchor is **`C6O2B2` / Fig 2a**,
  which passes. Fig 1b draws five further structures (MC3, SM-102, TT3,
  Lipid 14, OF-2) and **none of them is an anchor** — which is exactly the gap
  F-014 fell through. Adding those five as C10 anchors would have caught it on
  the day.
* **C11** (`reagent_registry.json`) pins 13 catalogue reagents by name. It
  reports `DLin-MC3-DMA ok, 3,097 rows`. `MC3-CB1` is not a registered name, so
  the derived conjugate is invisible to it. C11 checks the parent; nothing
  checked that the derived child was built from that parent.

---

## 4. EVIDENCE FROM THE PAPER

All page numbers are PDF pages of
`Advanced Materials - 2026 - Tian - Esterified Cholesterol Conjugates…pdf`
and `adma74762-sup-0001-suppmat.pdf`. Figures were rendered as images at
160–700 dpi with pymupdf; the structure panels carry no text layer.

**P. 3, Fig 1b, band "Commercially available cores"** — five parent structures.
`Dlin-MC3-DMA (MC3)`: acyclic chain `(CH₃)₂N–CH₂CH₂CH₂–C(=O)–O–CH<` branching
into two dilinoleyl chains, each with two *cis* alkene marks. **No ring.**
`Lipid 14`: `(CH₃)₂N–(CH₂)₃–C(=O)O–CH₂CH₂–N[(CH₂)₈–O–C(=O)–CH(hexyl)(octyl)]₂`,
both nitrogens neutral. `SM-102`, `OF-2`, `TT3` also drawn. *This panel is the
single piece of evidence that settles F-014.*

**P. 3, Fig 1b, band "Synthesized cores"** — amine heads C1–C6 and alkyl
acrylates O1/O2; band "Cholesterol tails" — B1/B2/B3 as Br–(CH₂)ₙ–C(=O)–O–
cholesteryl, with the ester on the steroid C3 (β-wedge drawn), n = 2/4/6.

**P. 3, Fig 1d** — heatmap, columns
`MC3-C | SM-102-C | Lipid 14-C | OF-2-C | TT3-C | C1O1 … C6O2`, rows `B1/B2/B3`.
`MC3-C` cells read 7.0×10⁴ (B1), 4.1×10⁴ (B2), 6.0×10⁴ (B3); `Lipid 14-C` reads
4.3×10⁴, 3.6×10⁴, 8.6×10⁴. These match `04c_fig1d_51.csv` exactly and reproduce
the stored z-scores.

**P. 3, Fig 1e** — bar chart; bracket "Benchmark" spans `MC3` and `SM-102`,
bracket "Cholesterol-modified" spans `MC3-C, SM-102-C, Lipid 14-C, OF-2-C,
TT3-C, C1O1 … C6O2`. Establishes that `MC3-C` ≠ `MC3` and that `MC3-C` is a
conjugate.

**P. 3, Fig 1 caption (quoted)** — "(b) Chemical structures of the lipid library
components, including amine head groups (C1–C6), alkyl groups (O1 and O2), and
cholesterol tails (B1, B2, and B3). … (e) Total brain luciferase activity
mediated by BBB-targeting LNPs formulated with cholesterol-conjugated lipids
(synthesized core structures are shown), and benchmark MC3 and SM-102 LNPs
formulated with unmodified ionizable lipids." *The caption omits the commercial
cores; the panel shows them.*

**P. 2 (quoted)** — "…or by directly utilizing commercially available cores such
as SM-102, DLin-DMA-MC3, and many more (Figure 1b; Figures S1–S7). These
tertiary amino lipid cores were subsequently functionalized with
cholesterol-derived alkyl bromides to generate a cholesterol-based quaternary
amino lipid structure." (The text's "DLin-DMA-MC3" is a typo for
DLin-MC3-DMA; the figure label is "Dlin-MC3-DMA (MC3)".)

**P. 4 (quoted)** — "conventional MC3 LNPs and SM-102 LNPs formulated with
unmodified ionizable lipids exhibited minimal brain delivery following
intravenous administration."

**P. 5, Fig 2a** — structure of **C6O2B2**, the only quaternized structure drawn
anywhere. The `N⁺` is on a **piperazine ring nitrogen**, with the
cholesteryl-ester C5 linker on it and the cholesterol C3-O ester clearly drawn.
Confirms `enumerate.py`'s ring-preference rule for C5/C6 cores.

**P. 11, Methods (quoted)** — "cholesterol tails (1.0 equiv) and commercially
available cores or synthesized cores (1.0 equiv) were dissolved in anhydrous
acetonitrile in a sealed vial … heated at 95 °C for 24 h." No regiochemistry
given.

**SI pp. 2–8, Figs S1–S7** — ¹H NMR of C6O2, C6O1, C5O2, B1, B2, B3, and
¹H/¹³C of C6O2B2. S7's C6O2B2 δ 4.54–4.51 (1H) is the cholesteryl-ester C3-H
(free cholesterol C3-H is ~3.5), confirming the C3 ester. **SI p. 14, Fig S13** —
structure of C6O2 (unconjugated control). **No commercial-core conjugate is
characterized or drawn anywhere in the SI** (pp. 9–23 are biology figures,
Table S1 formulations, Table S2 primers).

---

## 5. CONFIDENCE

| claim | confidence |
|---|---|
| Fig 1b (p. 3) draws MC3 as the acyclic ester, not the dioxolane | **certain** — read as an image at 700 dpi; the drawing has no ring at all |
| F-014: the **structure** was the wrong half; the name `MC3-C` is right | **certain** — the paper's own label plus its own drawing of the molecule under that label |
| Stored `MC3-CB*` were the KC2 dioxolane conjugate (pre-17:24) | **certain** — 5 rings, dioxolane SMARTS = 1, and a bit-for-bit rebuild from the KC2 parent |
| Root cause = `setdefault` name lookup × alphabetical corpus order × transiently ambiguous `DLin-MC3-DMA` | **certain** — reconstructed the pre-patch state and reproduced the stored SMILES exactly, all three |
| The applied repair (17:24) is correct and moved no values | **certain** — matches my independent rebuild; z-scores unchanged and re-derivable from Fig 1d |
| The repair is **not durable**: `PARENT`/`setdefault` unchanged, no check covers it | **certain** — lines 159/322 verified in the current file; C10 anchors and C11 registry inspected |
| F-015: second copy exists, is `max(index)`, docstring falsely claims parity | **certain** — quoted at `build_new_curation.py:199-204` and `:164-170` |
| F-015: exactly `Lipid 14-CB1/2/3` diverge, at the central vs terminal N | **certain** — both rules re-implemented and run over all five parents; stored rows match `max(index)` |
| The paper draws no conjugate of any commercial core | **certain** — full main text and all 23 SI pages checked |
| Porting `enumerate.py`'s rule verbatim does **not** fix Lipid 14 | **certain** — demonstrated by rewriting the same molecule tail-first; the two rules swap answers |
| The true site is the terminal dimethylamino | **probable** — chemical inference only (sterics + β-acyloxy EWG); no evidence in this paper |
| `SM-102-CB*`, `DSPC`, `DOPE`, `DMG-PEG2000` name lookups are clean | **certain** for the first three (formula, exact mass, connectivity, single corpus structure); **probable** for `DMG-PEG2000`, where PEG length is a distribution and no single SMILES is strictly truthful |
| Confidence laundering (24 `low` → `high`; 6 conjugates `high` against the QC's own rule) | **certain** — compared `03_enumerated.csv` against `new_curation_qc.json` |
