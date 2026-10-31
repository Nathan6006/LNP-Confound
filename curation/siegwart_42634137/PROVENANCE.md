# Where every number in this curation comes from

Tian et al., *Esterified Cholesterol Conjugates Enable LNP-Mediated mRNA Delivery
to the Blood-Brain Barrier*, Adv Mater 2026, PMID 42634137, doi 10.1002/adma.74762.

Both PDFs are **flattened raster** — one image per page, no vector art and no
extractable text inside any figure (the embedded images are 217 dpi native in the
main text, 220–600 dpi in the SI). Nothing here was read from a PDF text layer or
a vector path; every figure number was recovered from pixels, and every such step
has a script and a review overlay.

## The two things this paper counts, which are not the same thing

| | what varies | what is fixed | n |
|---|---|---|---|
| **51 lipids** | the lipid structure | formulation held constant | 51 |
| **486 LNPs** | the formulation | **one** lipid, C6O2B2 | 486 |

The abstract's "screening of 486 LNPs" is a formulation DoE, not a chemical
library. The paper states it directly (p. 4): *"the relative molar ratios of lipids
within LNPs were adjusted across the following ranges: C6O2B2 (35, 50, 65), DODAP
(10, 30, 50), helper phospholipid (10, 30, 50), and DMG-PEG (1, 2, 4), yielding a
library of **486 C6O2B2 LNP formulations**"* — 3×3×3×3 = 81 ratio combinations ×
6 helper lipids (DOPC, DMPC, DOPG, DOPE, EggSM, DSPC) = 486. Every one of the 486
carries `IL_SMILES` for C6O2B2 and nothing else.

Table S1 proves these levels are molar **ratios**, not percentages: it prints both,
and ratio/Σratio×100 reproduces all 56 printed percentages to 0.00 (`top14.py`).
Reading them as percentages would misstate every formulation.

## Source of each field

| data | n | where | how | script | overlay |
|---|---|---|---|---|---|
| lipid structures (6 amine heads, 2 acrylates, 3 chol tails) | — | Fig 1b | structure diagram, read as building blocks | `enumerate.py` | — |
| C6O2B2 full structure | 1 | Fig 2a | structure diagram | — | — |
| arm stoichiometry | — | SI NMR | δ4.04 (t,10H) + δ0.87 (t,15H) ⇒ 5 acrylate arms | `enumerate.py` | — |
| lipid SMILES | 36 | derived | RDKit from the blocks above | `enumerate.py` | — |
| **in vivo brain flux** | **51** | **Fig 1d** | **printed numerals, glyph clustering** | `digitize_fig1d.py` | `review/fig1d_overlay.png` |
| **in vitro screen** | **486** | **Fig 2c** | **4-bin discrete heatmap colour** | `digitize_fig2c.py` | `review/fig2c_overlay.png` |
| outstanding-LNP marginals | 15 | Figs 2d–g | bar height, tick-calibrated | `digitize_bars.py` | `review/fig2dg_overlay.png` |
| top-14 formulations | 14 | Fig 2j **and** Table S1 | printed numerals, double-sourced | `top14.py` | — |
| group means | 15 | Fig S11 | log-axis bar top | `digitize_figS11.py` | `review/figS11_overlay.png` |

Experimental conditions, from the Methods and figure captions, attached to every row:

- **in vivo (the 51)** — C57BL/6 mice, I.V., luciferase mRNA 0.1 mg/kg, 6 h,
  n = 3, whole-brain total flux in photons/s. Formulation held at
  47.6/23.8/23.8/4.8 (lipid/DODAP/DOPE/DMG-PEG), ethanol dilution.
- **in vitro (the 486)** — IGROV-1 cells, 20 ng mLuc per well, 24 h, n = 3
  biologically independent samples, RLU; LNPs made by pipette mixing.
  Total lipid/mRNA 40:1 wt/wt (Table S1).

## The 486: why a colour read is exact here

The Fig 2c legend is **discrete** — `<1e4 / 1e4–1e5 / 1e5–1e6 / >1e6` — so each
cell is one of exactly four printed colours and the read is nearest-neighbour
matching against swatches sampled from the legend in the same raster. No
interpolation, no inferred scale. Max distance from any cell to its matched swatch
is **4.72** on a 0–255 scale, with **zero** cells within 10 of a second swatch.

`bin_index` (1–4) is the defensible quantity. `rlu_midpoint` is a geometric
half-decade placeholder for numeric convenience and must not be used as a measured
value.

## Seven independent checks, all passing

1. **Four marginals** — Figs 2d–g counts vs the Fig 2c digitization: **15 of 15
   bins exact**. Measured by script, not by eye (`digitize_bars.py`), so no human
   sits on both sides of the comparison. Only Fig 2g has real degrees of freedom
   (6 bars, 2 fit parameters); its max residual is **0.08 LNPs**.
2. **Absolute tick calibration** — px-per-tick ÷ px-per-LNP gives 10.09, 10.15,
   20.01, 9.97 LNPs per tick across the four panels: all within 1.5% of a round
   interval, independent of the Fig 2c counts.
3. **Orientation** — the row-axis marginals are one non-monotone (41/21/31) and one
   monotone (8/31/54); a transposition or outer/inner swap breaks at least one.
   Both match, so the axis assignment is verified by the data, not assumed.
4. **Fig 2h distribution** — ink per decade band in the 486-point scatter: 17.1 /
   50.0 / 28.6 / 4.4 %, against 19.1 / 49.0 / 28.4 / 3.5 % from the bins. A seventh
   figure panel, same answer.
5. **Top-14 containment** — all 14 formulations in Table S1 are present among the
   486 **and all 14 sit in the >1e6 bin**, as they must, having been selected as
   top performers. By chance that is ~1e-10.
6. **Fig 1d colour bar** — cell fill mapped through the continuous white-to-blue
   bar, never touching the glyphs: Spearman **0.986** against the parsed numerals
   overall, **0.995** above 2e5 (the bar saturates below ~1e5).
7. **Fig S11 rank order** — matches wherever the authors report significance
   (DODAP 50 > 30 > 10, ****; DSPC top, ****). The orders that differ are exactly
   the comparisons they mark "ns", and are near-ties.

## Two places the paper disagrees with itself

Both are recorded rather than silently resolved.

- **240 / 91 vs 238 / 93.** The text says 240 acceptable and 91 outstanding.
  Digitization gives 238 and 93 — and the paper's own four bar charts each
  independently sum to **93**. The totals reconcile exactly: 238+93 = 240+91 = 331
  cells at or above 1e5. The prose is off by two at the 1e6 boundary; the figures
  and this digitization agree.
- **Fig 1d max.** C6O2B2-B2 at 3.0e6 is the largest value, and the paper
  independently names C6O2B2 as the lead compound. Consistent.

## Method notes worth keeping

- **OCR was rejected.** Tesseract misread the raised exponent on Fig 1d
  (`1.2×10⁵` → `109`) and dropped a decimal point. Glyph connected-components in a
  rigid 7-glyph layout, binarized and aspect-preserved, then complete-linkage
  clustered, gave 51/51 cells segmenting correctly with all furniture assertions
  passing. Aspect preservation matters: stretching glyphs to a square collapsed
  0/9, 6/8 and 5/6, and 5/6 in an exponent is a tenfold error.
- **Resolution matters more than algorithm.** At 300 dpi the mantissa clusters were
  impure no matter the metric; at 600 dpi, k=14 gives 14 pure clusters.
- **Log-axis bars have no meaningful height.** Fig S11's quantity is the bar top
  against the decade ticks. Overlaid scatter dots split the bar fill into
  components, and a bounding-box top reported the lower fragment — reversing
  DODAP 50 and 30. Median of per-column topmost fill fixes it.

## Deliberately not extracted

- **Per-point values from Fig 2h.** 486 dots at ~4.5 px radius overlap heavily
  (194 connected components for 486 points, largest ~10 dots merged), and the
  x-ordering carries no formulation labels, so even perfect deconvolution would
  give a distribution rather than per-formulation values. Used as a distribution
  check only (#4 above).
- **Fig S12 cytotoxicity.** Viability screen, not transfection.
- **Figs 3–5.** Mechanism, BBB integrity, protein corona — not library data.

## Files

```
digitize_fig1d.py    04c_fig1d_51.csv      07c_fig1d_qc.json
digitize_fig2c.py    04b_fig2c_486.csv     07b_fig2c_qc.json
digitize_bars.py     04d_fig2dg_bars.csv   07d_bars_qc.json
digitize_figS11.py   04f_figS11_means.csv  07f_figS11_qc.json
top14.py             04e_top14.csv         07e_top14_qc.json
enumerate.py         03_enumerated.csv     07_qc.json
review/              fig1d_overlay.png  fig2c_overlay.png
                     fig2dg_overlay.png figS11_overlay.png
                     fig1d_glyphs_mantissa.png  fig1d_glyphs_exponent.png
```

Every overlay draws the exact pixels each number was computed from. The two glyph
montages are the labelling surface: 14 mantissa + 4 exponent cluster
representatives were assigned by eye, and those 18 decisions are the only
human-supplied values in the Fig 1d read (plus one bottom-tick exponent per Fig
S11 panel).

---

## Corrections, 2026-10-30 (round-1 corpus review)

Three structural errors in this folder, all found by review rather than by a
check, and all now carrying a check of their own. Full detail in
`provenance/CHANGELOG.md` §12 and `provenance/review/round_2026-10-30/`.

### 1. Every structure in this folder stored epicholesterol (F-013)

All 36 enumerated lipids, all 537 corpus rows, and the 15 commercial-core
conjugates built in `../build_new_curation.py` carried the **3α** epimer.

`enumerate.py::_cholesteryl_ester` bonded the acid's oxygen to the sterol's C3
and then called `RWMol.RemoveAtom` on cholesterol's own hydroxyl oxygen, which
removes that oxygen from the middle of C3's neighbour list **without correcting
C3's chiral tag** — so the tag came to mean the opposite of what it meant.
Reproduced from scratch: genuine cholesterol through that code path gives C3 =
`R`; the same ester written as a plain SMILES string gives `S`.

Nothing in the project could see it. The epimers share molecular formula
`C117H218N5O12+`, exact mass 1885.6597, the InChIKey **skeleton** block, and all
17 descriptor columns bit-for-bit.

Now assembled the way the chemistry does it: the **acid** loses its OH and the
**alcohol's** oxygen is retained as the ester oxygen, so no atom adjacent to a
stereocentre is ever deleted. Verified against an independent plain-string
oracle that uses no `RWMol` at all — exact match for all three linkers, and the
**full CIP label vector** matches, not just C3.

That 3β is correct is settled from the paper, not assumed: Fig 1b (p.3) draws
B1/B2/B3 with the C3–O bond as a bold wedge **cofacial with the C18/C19 angular
methyls**; Fig 2a shows the same relation in a different orientation; ¹H C3-H is
4.71–4.58, the broad axial *tt* of a 3β-acyloxy sterol (3α would be near 5.0 and
far narrower); and the Fig S7 ¹³C spectrum — whose typed peak list omits C3
entirely — has an isolated peak at 73.9 ppm against ~70.9 for a 3α ester. No
epi, allo or 3α sterol is mentioned anywhere in the paper or SI.

**Note the history, because it is the point.** Fragment assembly was *mandated*
for this molecule after string concatenation fused the piperazine to cholesterol
through a SMILES ring-closure collision. The prescribed remedy then introduced a
silent stereochemical failure in the same molecule. Neither route is the safe
one; both need an external, stereo-aware check. That check is now `audit.py`
**C14**, which cleaves every ester inside `IL_SMILES` corpus-wide and pins the
liberated sterol's full InChIKey against `provenance/reagent_registry.json`.

### 2. The C10 anchor was a copy of the thing it checked (F-013/E2)

`reference_structures.csv`'s `reference_SMILES` was **byte-for-byte identical**
to `03_enumerated.csv`'s `C6O2B2` — same sha256, 223 bytes — while its own note
claimed it had been "assembled from fragments, independently of enumerate.py …
through the C3 oxygen (3-beta)". The string it annotated was 3α. So
`structure_check.py` reported `exact_match` **while comparing a value to
itself**, on this project's flagship curation, for the one molecule C10 exists
to protect. `ERROR_CLASSES.md` already listed this exact hazard.

Replaced by `make_reference.py`, which **builds** the anchor rather than
asserting its provenance, and shares as little as possible with the enumerator:
no `RWMol`, no `AddBond`, no `RemoveAtom` — plain string substitution — with the
amine half typed by hand from Fig 2a (p.5), the 5-carbon linker from Fig 1b
(p.3), and the sterol taken **verbatim** from `reagent_registry.json`, so the
note's claim is mechanically true instead of a promise. It self-checks before
emitting (one N⁺, five N, six esters, no radical, and its own sterol must cleave
to cholesterol) and it refuses to write an anchor that is itself wrong.

Confirmed able to disagree: it matches the corrected enumeration and does **not**
match the old 3α string. C10 now fails a stereo-only difference when the
reference's stereochemistry came from a pinned authority, and reports
`anchor_may_be_a_copy` for a byte-identical reference in a folder with no
generator script.

**Stated limit:** the sterol's eight stereocentres come from the registry, not
from a fresh read of the figure, so this anchor cannot catch an error in the
registry entry itself. That is why the registry pins a chemical name
("cholest-5-en-3beta-ol") and not only a SMILES string.

### 3. The quaternization site, and the second copy of this code (F-014/F-015)

`../build_new_curation.py` carried its own copies of both the ester assembly and
the site-selection rule, with a docstring claiming parity with this folder's
`enumerate.py`. It had no parity, and that cost two further errors:

- `MC3-CB1/2/3` were built on the **KC2 dioxolane** under an MC3 name (5 rings,
  1 dioxolane, against MC3's 0 and 0). The name is right and the structure was
  wrong — the *opposite* of the `AC_2025` repair — settled by Fig 1b, which
  draws "Dlin-MC3-DMA (MC3)" as the acyclic ester with no core ring. The cause
  was a corpus lookup resolved by **file order**: 292 component names map to
  more than one structure, `DLin-MC3-DMA` was one of them, and the corpus is
  written in alphabetical experiment order. That lookup now raises instead of
  choosing.
- `Lipid 14-CB1/2/3` were alkylated at the hindered central nitrogen. The fix
  was **not** to port this folder's repaired rule verbatim: its β-arm counter
  walked only `N-C-C-C(=O)` and could not see Lipid 14's `N-C-C-O-C(=O)`
  β-**acyloxy** arm, so both nitrogens scored 0 and the rule fell through to the
  atom-index tiebreak — the SMILES-ordering artifact it had been rewritten to
  eliminate. β-acyloxy is now counted too, and the choice is invariant across
  four spellings of each parent.

Both copies are gone: `build_new_curation.py` imports this folder's
implementation. `test_enumerate.py` tests both functions at the function level —
the string oracle, the full CIP vector, and the ordering-invariance of the site
rule — rather than 537 rows later.

### 4. 24 disclosed `low`-confidence rows were labelled `high` (F-011)

`03_enumerated.csv` flags 24 of the 36 lipids `confidence = low` with the
verbatim reason *"ASSUMED: no ring N and the paper draws no structure"*.
`../build_new_curation.py` read only `lipid_name → IL_SMILES` and then wrote the
literal `"high"` for every lipid it found, discarding both that column and
`quaternization_site_assumed`. So the rows whose regiochemistry the paper draws
nothing for were labelled this project's highest confidence, and the only
surviving provenance asserted the opposite of the truth.

They now read `low`, with that sentence as the row's `curation_note`, and the 15
commercial-core conjugates read `medium` with `structure_assumed = 1` — the
paper draws no `-C` conjugate for any of the five, so every one of those sites
is inferred. The site remains **unverifiable** from this paper; the deliverable
is the flag.

**No `Experiment_value` changed in any of this.** The 486 in-vitro and 51
in-vivo z-scores are unchanged.
