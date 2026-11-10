# PROVENANCE — `mitchellmj_41603124`

Hamilton AG, Thatte AS, Xu J, … Weissman D, **Mitchell MJ**.
*"High-Throughput In Vivo Screening Using Barcoded mRNA Identifies Lipid
Nanoparticles With Extrahepatic Tropism for In Situ Immunoengineering."*
**Adv Mater** 38(13):e14370, Mar 2026. PMID **41603124** · PMCID **PMC12957869**
· DOI 10.1002/adma.202514370 · **CC BY**.

Proposed `Experiment_ID` **`AH_2026`** (first author Hamilton). Collision check
executed: 54 `Experiment_ID`s in the live corpus, **0** beginning `AH`.

**Rows emitted: 0.** `06_rows.csv` is header-only, 79 columns, so the zero is
visible. Four entries in `DECISIONS.md` are `OPEN`, which is itself a merge
block.

---

## Why zero — two independent blockers, either one sufficient

**1. The screen's endpoint is not delivery.** 3,840 attributable values of
*normalized intracellular b-mRNA accumulation* — barcoded mRNA recovered from
sorted cells and sequenced. The authors sampled deliberately **before** peak
translation and say so, then ran Figure 3 separately *"to confirm
transfection"*. `data/`'s `Experiment_method` vocabulary is delivery-only and
this is none of its four values. Pooling is **not** the objection: `LR_2023`
contributes 407 of 1,803 in-vivo rows at `Experiment_batching = barcoded` — but
with `Cargo_type = peptide_barcode` and `Experiment_method =
protein_adundance_normalized`, a **translated** readout. (D-1)

**2. The validation panels are admissible endpoints and attributable, but the
structures are not there.** Of the 7 formulations in Figure 3, **1** ionizable
lipid is determined from this paper (`18D4`) and **1** comes from the registry
(`SM-102`); the other 5 are not determinable. `data/` has 0 blank `IL_SMILES`
in 30,998 rows. The two that drop out carry the panel's two lowest liver
values, so a partial emit would truncate the range every z-score is computed
over. (D-2)

---

## The three numbers

| | |
|---|---|
| **distinct ionizable lipids** | **133** named (120 `nDm` library + 13 catalogue), of which **1** is fully determined from this paper, **2** more are drawn but unreadable, **6** composable in principle, **111** not determinable, **13** catalogue (6 pinned by the registry) |
| **distinct formulations** | **133**, all in Supplementary Table 3, three recipes, **133/133 mol% sum to exactly 100** |
| **individually measured** | **7** (LNP 9/37/94/96/97/112 + SM-102) × 2 organs in Fig 3b; the same 7 × 9 populations in Fig 3f/3j |
| **what the design could support** | 120 × 16 populations × 2 strains = **3,840** per-LNP screen values |
| **what the publication attributes** | **3,840** — Supplementary Figure 14 labels every axis. This paper does *not* have the sibling's attribution failure; it has an endpoint failure |

Supplementary Table 3's `Index` is **not** monotone in the amine number
(…8, 10, 11, 13, 12, 14 … 20, 9), so the index was **read**, never computed.

---

## Anchor posture — recorded BY HAND

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run**: they write tracked files other
agents read. What those checks would find:

- `reference_structures.csv` has the exact 7-column header and **9 rows, 1
  carrying a SMILES**.
- **`18D4` is a genuine, independent C10 anchor.** `03_enumerate.py` composes
  it by plain SMILES string interpolation, **arm-first**. `make_reference.py::build`
  builds it with `assemble.attach` on atom-mapped atoms, **head-first from the
  piperazine**. Both give `MFMRADGRKUHRDG-UHFFFAOYSA-N`. The two routes share
  no code.
- **Independence is enforced at run time by AST, not by a substring scan** —
  a text scan would flag this file's own prose. `_assert_independent()` walks
  `build()`'s AST and refuses on any `JoinedStr`, `+`/`%` on a string, or
  `.join`/`.format`, and requires at least one `assemble.attach` call.
  **The guard is made to fail on purpose**: `--selftest` runs it against a
  deliberately string-built function and requires it to raise. My first
  `build()` *failed* that guard, correctly — it used `%` to set an atom map.
- **The gate is proven to test connectivity**, by seeding a regioisomer of
  `18D4` (the secondary N keeps its H; one primary N takes both arms):
  formula `C46H85N5O12`, average mass 900.209 and monoisotopic mass
  **899.619473** all identical, **Δ 0.000000 Da**; only the InChIKey separates
  them (`MFMRADGRKUHRDG` vs `UHVAORHCRFDQOO`).
- **`17D9.2` and `20D8i` carry empty `reference_SMILES` with the reason on the
  row**, so the absence is visible. They are *drawn* and still not
  determinable — the distinction is the point.
- **The six excipients/benchmarks carry empty `reference_SMILES` too, and this
  is not the sibling's excipient-anchor case.** That sibling's paper *drew*
  Cholesterol and DOPE, so it could anchor them. This paper only **names**
  them in Supp Table 3 and draws no excipient anywhere, so there is nothing of
  its own to compare against. Recorded instead as a check: all six resolve
  identically in `provenance/reagent_registry.py` and in the live corpus,
  **6/6 by full InChIKey**.
- `C14-PEG2000`, `DMG-PEG2000` and `ALC-0159` are deliberately unpinned — a PEG
  block length is a distribution (CLAUDE.md §3).

**C18 would report 8 of 9 anchor rows empty for this folder.** That is the
honest state and it should be baselined with the reason, not closed.

---

## What was measured, in numbers

Full detail in `07_qc.json` (23 checks). The ones that carried weight:

- **133/133** formulations sum to exactly 100 mol% (trap 8 clean, 0 offenders).
- **6/6** agreement of the lead-LNP identity across three independent loci.
- **18D4**: `a+b = 8.017` carbons from the printed ions against the drawn 4+4;
  Δ 0.612 Da on `[M+H]⁺` with the nearest alternative 28.053 Da away —
  a **45.8×** discrimination margin.
- **17D9.2**: `a+b = 18.993` across three independently labelled ions. So the
  token `9.2` is **not** a chain length — 19 total, not 2 × 9.2 = 18.4 — and the
  drawn maleate is unsymmetrical, 19 being odd. Trap 1, live.
- **Fig 3b calibration**: x residual **0.3988 px** over 8 decade ticks (LOO max
  0.571), y residual **0.2895 px** over 6 (LOO max 0.48). Decade ticks are
  separated from the 5× minors by **measured tick length** — 4 px vs 3 px,
  alternating **12/12** — not by taking every other one.
- **The printed `18.1×` gate**: measured **18.123×**, px error **0.01**.
  Seeded (SM-102 shifted 4 px) it moves to **10.993×** and the command exits 1
  against 0 clean — each run its own process, verdict from `returncode`.
- **Fig 2 column count**: **124** (sweep 123–125) against 122 printed in the
  abstract and 133 in the Conclusions. Neither printed number fits; 133 is
  excluded by 8 bins.
- **Precision honesty**: 18.42 px per decade ⇒ one pixel is a factor **1.133**,
  and `axis.TICK_RESIDUAL_PX = 2.5 px` — set for 600 dpi renders — would permit
  a factor-**1.37** error here and pass. The gate is not wrong; it is being
  asked to protect a panel it was not scaled for. Native dpi is **not
  measurable** (`axis.native_tol` refuses a guess) because no PDF rendition of
  the article is obtainable; the px-per-decade is reported instead.

## Corpus overlap

- **0** of 120 `nDm` names match any corpus `IL_name`.
- **4** catalogue names are present: `C12-200`, `DLin-MC3-DMA`, `SM-102`,
  `ALC-0315`.
- **No structure sweep was run, and that is deliberate**: 111 of 120 members
  have no determinable structure, so there is no InChIKey to compare and a
  substructure "superset" against structures we do not hold would be
  meaningless.
- **Token-pattern test, before any name join.** Amine tokens are the contiguous
  run **1–20 with no gaps** — the signature of local indices, not of the sparse
  shared catalogue, in which **0** bare numerics 1–20 appear at all. The damage
  a naive substring join would have done was *measured*: token `"1"` captures
  **20,837 rows / 14,437 distinct structures / 50 experiments**. **No name join
  on an amine token is performed anywhere in this curation.**

## Already-merged arms, and the two staged same-lab siblings

No arm of this screen is already merged. Three Mitchell-lab folders, three
unrelated numbering systems — **same lab ≠ same numbering**:

| folder | grammar | overlap with this paper |
|---|---|---|
| `mitchellmj_38982196` | `<amine>T-O<tail>`, amines 1–25 | **0** shared names, no shared chemistry |
| `mitchellmj_41555650` | `C<tail>-<core>`, cores 114/200/306/490/494 | **3** shared names: `C12-200`, `C12-494`, `C14-494` |
| this paper | `<amine>D<maleate>`, amines 1–20 | — |

`C12-200` matches the corpus (`AC_2025`) by full InChIKey; `C12-494` and
`C14-494` are not in the corpus. **But this paper attributes no value to any of
the three**: the 13 catalogue reagents appear only in Figure 2's *unlabelled*
heat map, Supplementary Figure 14 is the 20×6 `nDm` grid and excludes them, and
Figure 3b's only benchmark is SM-102. So the shared chemistry yields no usable
cross-source pair.

**Do benchmarks top the screen?** In Figure 3b, SM-102 has the highest liver
luminescence of the seven (2.83e11 p/s) — **18.1× LNP 112**, the paper's own
printed figure. On spleen the order reverses: LNP 112 at 6.90e8 ≈ SM-102's
6.89e8 while the paper's lead claim is tropism, not magnitude. Relevant to
CLAUDE.md §8's open question on benchmark-only pairs, and noted because the
SM-102 comparison **confounds lipid with formulation** — SM-102 uses the
clinical recipe, all six `nDm` LNPs share one different recipe.

## Hazards hit in this run

- **The session scratchpad is SHARED with the sibling curation agents.** A glob
  copy out of it put **17 figure files from other papers** into `source/figs/`.
  Caught by listing the directory immediately after the copy; the 6 survivors
  were re-hashed against my own download record, **6/6** identical. A file-count
  or size check would have passed.
- **Prose lagged its own number, twice** — once in a templated verdict sentence
  that went on asserting "the figure sides with the abstract" after the measured
  column count had moved, and once when I let my own probe values (123, 124) be
  reported as "printed numbers". Both are in `08_contradictions.csv` as
  retracted suspicions with the measurement that retired each, rather than
  edited away.
- **Error class 20, live, in my own first border detector**: "the outermost
  full-height dark column" anchored on a saturated purple heat-map **cell**,
  because a high-value cell is darker than the panel's grey border. Fixed by
  selecting on saturation. The colour ramp's floor is **white**, so a panel's
  span can never be taken from where colour starts.
