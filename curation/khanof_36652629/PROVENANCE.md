# PROVENANCE — khanof_36652629

Tilstra, G.; Couture-Senécal, J.; Lau, Y. M. A.; Manning, A. M.; Wong, D. S. M.;
Janaeska, W. W.; Wuraola, T. A.; Pang, J.; Khan, O. F.
**Iterative Design of Ionizable Lipids for Intramuscular mRNA Delivery.**
*J. Am. Chem. Soc.* 2023. PMID **36652629**, DOI **10.1021/jacs.2c10670**.
Corresponding: Omar F. Khan, University of Toronto.

Curated 2026-11-10. **Nothing in this folder has been merged into `data/`.**

---

## 1. Title and endpoint, verified first

The brief warned that a sibling burned its budget on a paper whose truncated
title hid a different endpoint. Both were checked before any digitiser was
written.

- **Title**, from the Europe PMC core record for PMID 36652629 queried
  2026-11-10: *"Iterative Design of Ionizable Lipids for Intramuscular mRNA
  Delivery."* The SI cover page prints the same title. Not truncated.
- **Endpoint**: total luminescent flux (p/s) from firefly luciferase in the
  **hindlimb muscle** of live **C57BL/6** mice after **intramuscular**
  injection, imaged on an IVIS Spectrum. Read verbatim off the caption of every
  value panel. This is an in-vivo delivery endpoint and it is the same endpoint
  in all six panels.
- Explicitly **not** the endpoint, and not emitted: apparent pKa (TNS),
  predicted LogD (MarvinSketch), size/PDI/EE, molecular-dynamics N–P contacts,
  tdTomato gene-editing fluorescence (Fig S6, different cargo and units), and
  mRNA copies per particle (Tables S5/S6).

## 2. Sources held, and what could not be obtained

| | |
|---|---|
| `source/ja2c10670_si_001.pdf` | the Supporting Information, 39 pages, 3.0 MB, full text layer on 38 pages |
| main article | **NOT HELD.** Europe PMC reports `isOpenAccess=N`, `inEPMC=N`, `hasSuppl=N`, `pmcid=null`. No paywall or bot-check bypass was attempted. |
| author repository | **NONE.** `provenance/tools/find_repo.py --doi 10.1021/jacs.2c10670 --pmid 36652629 --title "iterative design ionizable lipids intramuscular" --author Khan` returned `repo_hits: []`, 0 candidates; 17 unrelated repositories were filtered out. The SI names no repository, Zenodo deposit or Dryad DOI. |

Each failure is a checkable sentence, not "not available".

## 3. Value tier: 5 (pixels), and why 1–4 and 4b all fail

Recorded in `00_value_sources.json` before any digitiser existed, and
re-runnable with `python3 -I stage00_scan_for_values.py`, which writes
`00_scan_for_values.json`.

- **tier 1** — ACS publishes no Source Data product for JACS, and there is no
  PMC asset set to enumerate.
- **tier 2** — the SI has exactly six tables and none carries a flux value:
  S1 and S3/S4 are sequences, S2 is physicochemical characterisation of all 28
  formulations, S5/S6 are mRNA copies per particle for two formulations.
- **tier 3** — 0 repository hits, see above.
- **tier 4** — measured: across all 39 pages and **43,853 extracted characters**
  there are **0** lines carrying luminescence vocabulary together with a number
  in a plausible flux range. Every large numeric token in the text layer is an
  NMR shift, a coupling constant, a reagent mass, a yield or an HRMS *m/z*.
- **tier 4b `pdf_vector_geometry`** — checked explicitly, because a sibling read
  54 values straight out of filled-rectangle coordinates. **It does not apply.**
  Figure pages 2–8 carry **0 drawings, 0 rectangle items and 0 filled paths**,
  and 1–3 embedded raster images each (0.6–2.2 megapixels per page). The panels
  were flattened before deposit. By contrast pages 9–17, the synthesis schemes,
  **do** carry vector line and curve geometry (12–83 paths per page), which is
  why the structures are not a pixel problem even though the values are.

## 4. The three numbers

| | count | what it is |
|---|---|---|
| **distinct structures** | **11** | 8 ionizable lipids synthesised here + 3 benchmarks (SM-102, DLin-MC3-DMA, ALC-0315) that the paper formulated but did not make |
| **formulations** | **28** | Table S2's formulation IDs. One lipid, many batches: SM-102 is 7 of them, OC2-K3-E10 is 8, C3-K2-E14 is 2. The suffix-free names hide this; the ID does not. |
| **individually measured** | **18** | bar means in the six delivery panels. 10 distinct structures (ALC-0315 appears in no flux panel), across 6 rounds. |

Not a combinatorial library: 3 cores × 2 linkers × 5 tails is 30 possible
lipids and 8 exist. "Iterative design" is literal here.

## 5. Corpus overlap — report only, nothing acted on

`stage04_corpus_overlap.py` → `04_corpus_overlap.json`. Three keys, because
they fail differently: full InChIKey, InChIKey skeleton block, and `IL_name`.

- **0 of the 8 novel lipids** hit anything, on any key, in either pool.
- **3 of 3 benchmarks** hit, and on **all three keys at once**:
  - SM-102 — vitro `AC_2025`, `WT_2026`; vivo `SB_2024`, `SS_2018`, `SY_2025`, `XD_2026`
  - DLin-MC3-DMA — vitro `AC_2025`, `SP_2020`, `WT_2026`, `YZ_2022`, `YZ_2024`; vivo `BL_2023`, `NG_2025_2`, `SB_2024`, `XD_2026`, `ZC_2023`, `ZH_2023`
  - ALC-0315 — vitro `AC_2025`, `LZ_2025`; vivo `XD_2026`
- That triple agreement is a **free external check on our structures**: the
  benchmark SMILES came from `provenance/reagent_registry` (C11) and reproduce
  the corpus exactly, independently of this paper.
- The only shared chemistry with the corpus is therefore benchmarks. See
  `DECISIONS.md` D-2 and `CLAUDE.md` §8.

## 6. Checks, with numbers

All in `07_qc.json`; run `python3 -I stage07_qc.py`.

| check | result |
|---|---|
| **A** stroke-centreline convention | worst disagreement between the baseline rule's centre and the lowest tick's centroid, over 6 panels: **0.013 px** (0.000/0.003/0.002/0.000/0.000/0.013). Neither quantity comes from a bar. |
| **B1** our formula vs SI formula | **8/8** exact |
| **B2** SI formula vs SI's own "calc'd" *m/z* | **8/8** within 0.016 Da — but on **two different mass scales**: average on C3-K2-E14 and C3-K2-E10, monoisotopic on the other six |
| **B3** ours vs SI "found" *m/z* | **8/8**; on C3-K2-E14 the calc and the found are on **different scales from each other** (C-2) |
| **C** Fig S11 pKa vs Table S2 | **10/11** exact to the printed 2 dp, across two separately typeset sources. The miss is ALC-0315, which is absent from Table S2 entirely (C-6). 28/28 Table S2 rows parsed. |
| **D** cross-round ratios | 5 lipids in more than one panel: **1.053x, 1.662x, 1.751x, 1.861x, 5.011x** |
| **E** Fig S9A, unprinted | drawn 2.06 and 5.93; ours 3.162x and 7.251x. **Does not agree** — recorded as C-5 and D-4, not as a pass |
| **G** the paper's own significance marks | **12/12** sorted by a single fold threshold at 2.074x. Largest "ns" 1.945x, smallest starred 2.212x, margin **0.056 decades**. A t-test is not a function of the fold alone, so this is a real test of the bar tops. |
| **F** seeded defects | **3/3 fired**, each in its own process with `returncode` read from the process object, never through a pipe |

### What the seeds were, and what moved

1. **`dots_on_both_top_corners`** — two saturated dots of the bar's own colour
   pasted onto *both* top corners of Fig S1B's SM-102 bar. This is the one shape
   `bars.bar_top`'s flank test cannot see, because two symmetric dots agree with
   each other perfectly. Caught by the rule-thickness test: the false rule was
   thinner than the side wall drawn with the same pen. Exit 1.
2. **`tick_error_only_LOO_can_see`** — this folder applies `axis.LOO_TOL_PX` in
   the pixels of the *original* raster, which is 4.07 render px on S4B rather
   than 1.5 (see §7). A widened number has to be shown to still have teeth, so
   the seed searches for a displacement the max-residual test **passes** and
   leave-one-out **rejects**. One exists: tick 0 moved **−3.25 px**. Exit 1.
3. **`halo_widened_window`** — the bar window widened by 8 px each side, i.e.
   `trim_to_walls` undone. The top rule then covers less than 90% of the
   inflated width and a flank fails. Exit 1.

Each is credited because the **agreement count moved**: the control run passes
and the seeded run refuses.

## 7. Two tolerances were widened, and both are measured, not chosen

Every panel is a flattened raster at a **measured** native resolution —
219.88, 220.26, 220.27, 221.25, 221.25 and 188.49 dpi, computed from each
embedded image's pixel dimensions divided by its placement rectangle. Rendering
at 600 dpi does not create information.

- `axis.LOO_TOL_PX` (1.5 px) and `bars.FLANK_TOL_PX` (2.0 px) are tolerances in
  **real** pixels. At 600 dpi over a 220 dpi image, 1.5 px is **0.55 of a real
  pixel** — below the raster's own quantisation, so no correct detection can
  meet it.
- Both are therefore applied as `tol × 600 / native_dpi`, and **both numbers are
  recorded per panel** in `05_digitize.json` (`loo_worst_px`,
  `loo_worst_native_px`, `loo_tol_native_px`, `flank_tol_px_used`,
  `flank_tol_native_px`).
- Worst leave-one-out residual in native pixels across all six panels:
  **0.689**, against the library's 1.5. Worst tick-fit residual: 1.25 px.
- The widening was then seeded against (seed 2 above). This is the one judgement
  call in the folder and it is the reason C-11 is recorded as a *retracted*
  suspicion rather than quietly fixed.

Nothing else was loosened. `bars.bar_top`'s 0.90 coverage rule and `nrow=3` are
untouched; six refusals were fixed by correcting the **mask and the window**,
never the threshold. The full sequence of refusals and what each one turned out
to be is in the `stage05_digitize.py` docstring.

## 8. Anchor posture — recorded by hand

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run** (out of scope for this task), so the
posture is stated here instead.

- `make_reference.py` builds **8 anchors**, one per synthesised lipid, and all
  **8 agree with the enumeration by full InChIKey**.
- **The code path is different, and executed, not asserted.** The enumerator
  (`stage01_enumerate.py`) uses RDKit `RWMol` + `CombineMols` + `AddBond` with
  atom-mapped attachment atoms, builds core-first, and imports `curation/lib
  assemble`. `make_reference.py` uses **plain Python string concatenation**, no
  atom maps, writes **tail-first**, and imports no project module; RDKit is used
  only to parse the finished string and print an InChIKey.
- **String assembly is normally banned** because it once fused piperazine to
  cholesterol through a ring-closure digit collision. That hazard needs rings.
  Every structure here is acyclic, and the script **asserts** that no anchor
  contains a ring-closure digit before writing anything. Where the hazard cannot
  occur, the route sharing least with the enumerator is the better anchor.
- **The anchors earned their keep immediately.** The first three runs disagreed,
  at 0/8, 1/8 and 5/8, for three instances of one class — a SMILES fragment has
  a direction:
  1. both arms written outward, so the core nitrogen bonded to a **tail**;
  2. the tails reversed inside their branches, so each nitrogen sat **past the
     carbinol** at the far end of the chain;
  3. `OC2`/`OC3` written oxygen-first, giving **N–O–CH2CH3** in place of
     N–CH2CH2–OH, an amine-N-oxide ether.

  **Every one of the three had the correct molecular formula**, parsed,
  canonicalised and carried no radical. C3-K2-A14 passed through the second of
  them because a plain alkyl tail reads the same in either direction. Nothing
  but the InChIKey comparison saw any of it. Mass cannot see connectivity, and
  neither can a formula check.
- **2 of the 8 anchors are byte-identical to the enumeration** (C3-K2-A14,
  OC2-K2-A10) — the two with direction-symmetric tails, where a correct
  re-derivation canonicalises to the same string. `structure_check.py`'s
  `anchor_may_be_a_copy` would flag these; the folder **has** a generator
  script, so the posture is the correct one. Byte-identity is not the signal
  (F-037), which is exactly why the other 6 differ byte-wise and agree
  chemically.
- **3 declared gaps**, with empty `reference_SMILES` so the absence is visible:
  SM-102 (drawn in Fig S1A but not traced — structure from C11 and corpus-matched
  3 ways), DLin-MC3-DMA and ALC-0315 (never drawn in the SI).
- **Stereochemistry is deliberately absent** from all 8 anchors and all 8
  enumerations. The carbinols come from racemic epoxides, so the compounds are
  stereoisomer mixtures. C10's stereo comparison therefore carries **no
  information** for this folder, which each anchor's note says rather than
  letting it look like agreement.

## 9. Files

| file | what |
|---|---|
| `00_value_sources.json` | the intake gate, written FIRST |
| `stage00_scan_for_values.py` → `00_scan_for_values.json` | the measurements the gate cites |
| `stage01_enumerate.py` → `01_blocks.csv`, `03_enumerated.csv` | 10 blocks; 8 lipids + 3 benchmarks |
| `stage04_corpus_overlap.py` → `04_corpus_overlap.json` | three-key overlap, report only |
| `stage05_digitize.py` → `05_digitize.json` | the six panels, tier 5 |
| `stage06_values.py` → `05_values.csv` | 18 rows, 38 columns, flags via `emit.emit_rows` |
| `stage07_qc.py` → `07_qc.json` | checks A–G and the seeds |
| `08_contradictions.csv` | 8 open, 3 retracted with the measurement that retired each |
| `DECISIONS.md` | D-1 … D-7, all `STATUS: OPEN` |
| `make_reference.py` → `reference_structures.csv` | 8 anchors + 3 declared gaps |

Rerun order: `stage00_scan_for_values` → `stage01_enumerate` →
`stage04_corpus_overlap` → `stage05_digitize` → `stage07_qc` →
`stage06_values` → `make_reference`. All `python3 -I`.

## 10. What is NOT verified

- **The dose**, on every row. Not stated anywhere in the SI (D-3).
- **The timepoint** on Figs S4B, S4C and S7A. Stated only on S1B/S2B/S3B (C-8).
- **Fig S8**, the dose–response series — its lipid identity is contradicted
  between its caption and its own legend (C-1), so no row is emitted from it.
- **Fig S9A's fold increases**, which disagree with the panels they derive from
  by 1.53x and 1.22x (C-5). Unresolved, and recorded rather than reconciled.
- **Which formulation Fig S4D's OC2-K3-E10 column describes** (C-4, D-5).
- **The regiochemistry of epoxide opening** is taken from the paper's own
  drawings (Figs S1A/S2A/S3A/S4A show the hydroxyl on the second carbon) and
  from the standard anti-Markovnikov outcome; HRMS cannot see it, and no NMR
  assignment in the SI distinguishes the two regiochemistries directly.
- **Absolute stereochemistry**: unknowable, and genuinely absent — the products
  are mixtures (D-7).
