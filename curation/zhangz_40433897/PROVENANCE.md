# PROVENANCE — zhangz_40433897

Cui K, Zhang F, Shi T, Zhao S, Zhou Y, Liu X, Hu Y, Hu Z, Kong L, **Zhang Z**.
*Iterative Screening of Vitamin E-Based Functional Lipid Nanoparticles for mRNA
Delivery.* **ACS Nano** 2025;19(22). PMID **40433897**,
DOI **10.1021/acsnano.5c01378**.

**Endpoint: nucleic-acid delivery** — mRNA firefly-luciferase luminescence, mRNA
and plasmid GFP⁺ fraction by flow cytometry, and siRNA (siLuc) knockdown of a
luciferase reporter. In scope.

**Original research, not a commentary.** Europe PMC `core` record, fetched
2026-11-10: `abstractText` length **1199** (non-zero), and the SI carries 17
figures of measurements, a 60-row orthogonal formulation design and 8 HRMS
characterisations. The discriminator used was abstract length plus the presence
of Methods-grade material, not `pubType` (which reads only "Journal Article").

---

## 1. What was acquired, and what was not

| | |
|---|---|
| held | `source/nn5c01378_si_001.pdf`, 1,176,524 bytes, sha256 `c5a3248e3ec19d6e9532c546cdf325a71b338a8e27b5db5163f6126809b509d2`, 11 pages |
| | Pre-supplied in `source/` by the caller; its own watermark reads "Downloaded from pubs.acs.org/… by PHILLIPS EXETER ACADEMY user on 10 November 2026". |
| **not held** | **the main article.** `isOpenAccess=N`, `inEPMC=N`, `inPMC=N`, `pmcid=null`; the only full-text URL is the DOI, "Subscription required". NCBI `elink dbfrom=pubmed&db=pmc&id=40433897` returns **only** `pubmed_pmc_refs` (7 citing PMCIDs) and **no `pubmed_pmc` self-link**, so no PMC deposit exists and there is no PMCID to fetch. ACS needs a browser and this run was forbidden the shared browser. No paywall or bot check was bypassed. |
| not held | no author repository. `provenance/tools/find_repo.py --pmid 40433897 --title "vitamin E lipid nanoparticle mRNA" --author Zhang` returned nothing above noise; no Zenodo/Dryad DOI appears in the SI. |

**Consequence, stated plainly:** main-text Figures 1–6, which carry the full
iterative screen over the 35-amine library, were never seen. Everything in this
folder is SI-only, and the screen's main results are **not** in it.

---

## 2. Tier ladder — `00_value_sources.json`

| tier | status |
|---|---|
| 1 publisher Source Data | **absent.** One supplementary asset exists for this DOI. Recorded as UNVERIFIED-ABSENT: the landing page could not be checked without the browser. |
| 2 SI numeric table | **present, but carries no endpoint.** `find_tables()` recovers Table S1 (flow antibodies) and Table S2 (60 prescriptions, mol%). Table S2 is formulation metadata, not a measurement. |
| 3 author repo | absent |
| 4 PDF text/vector | absent for values. The only numerals in the text layer are Table S1/S2 cells, the 8 HRMS pairs, and figure-caption *n*'s. |
| **4b vector geometry** | **impossible, with a checkable reason.** `pymupdf.get_drawings()` returns exactly **2 paths on each of pages 4–11**, and those two are the ACS download-watermark rectangles, not figure content — pages 1–3, which carry no figure at all, have 2, 7 and 3. **Zero** vector bond paths and **zero** vector axis paths exist anywhere in the document; every figure is one flattened DCTDecode JPEG. |
| **5 pixels** | the only route to any value here. |

`triage.py`'s `table_pages` reported `[1]` — one page — while
`pdfplumber.find_tables()` succeeds on all 11. Its text heuristic was wrong
again, in the direction the brief warned about. Its token count was 322 under
`numbers_regex = \d+\.\d+|\d{3,}`.

---

## 3. Native DPI, and both tolerance figures

Measured per embedded raster as `pixels / rect_points × 72`, on the JPEG itself
and not on a page render.

**Range: 173.8 – 220.5 dpi.** The seven digitised panels:

| panel | px | dpi x / y |
|---|---|---|
| Fig S4 (incl. S4b) | 1146 × 583 | 219.9 / 219.9 |
| Fig S5 (a and b) | 1267 × 305 | 219.9 / 219.5 |
| Fig S7 | 436 × 438 | 219.5 / 219.6 |
| Fig S8 (a and b) | 803 × 412 | 216.2 / 216.2 |
| Fig S13 | 648 × 377 | 189.9 / 190.0 |

`curation/lib` tolerances are written for a 600 dpi render, so at ~220 dpi the
scale factor is **×2.731** (×3.159 at Fig S13's 190 dpi):

| | 600-dpi figure | native-px equivalent |
|---|---|---|
| `axis.TICK_RESIDUAL_PX` | 2.5 px | **6.827 px** (6.827 / 7.895) |
| LOO residual tolerance | 1.5 px | **4.096 px** (4.092 / 4.099 / 4.163 / 4.739) |

**Measured worst leave-one-out residual, per panel, in native px:**
Fig S5a 0.429 · Fig S5b 0.277 · Fig S4b 0.500 · Fig S7 0.000 ·
Fig S8a 0.500 · Fig S8b 0.250 · Fig S13b 0.500.

**The rescale is not load-bearing here**, stated two ways so neither framing
hides the other:

* taken as **native px**, the worst residual anywhere is 0.500 px, which clears
  the unrescaled 1.5 px number on **7 of 7** panels — this is what
  `07_qc.json`'s `loo_passes_unrescaled_tol` reports;
* taken as **600-dpi equivalents** (×600/native), they are 1.170 / 0.757 /
  1.364 / 0.000 / 1.388 / 0.694 / **1.580**, so six of seven clear 1.5 px and
  Fig S13b (the 190 dpi panel) misses it by 0.08 px.

Both numbers are recorded per panel in `07_qc.json` check G.

**The widened tolerance was seeded anyway** (`07_qc.json` seed
`G_axis_loo_native_tolerance`, exit 1): on Fig S5a's real tick ladder
(63.0/128.0/194.0/259.0 for 120/80/40/0), displacing one tick gives worst-LOO
1.333 px at +1.0, 2.667 at +2.0, **4.333 at +3.25 — which still passes the
4.096 px rescaled tolerance by 0.24 px** — and 5.000 at +4.0, which fails. So
the rescale costs real sensitivity between 2 and 3.5 px and is stated as such,
not hidden.

---

## 4. The C10 anchor posture — recorded by hand, as instructed

> **There is no anchor for any VNxCyOz lipid, and that is the honest posture,
> not an omission.**

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run** — they write tracked files other
agents read. This section is the hand record.

**Why there is no anchor.** C10 compares an enumerated structure against the
structure the paper **draws**. This paper draws **no VNxCyOz structure
anywhere**: Figure S1b draws the product with a generic `R`, and Figure S2 is
NMR and HRMS traces with no structure diagram. `reference_structures.csv`
therefore carries **8 rows with an empty `reference_SMILES`** and
`expected_mismatch = no_anchor_paper_draws_no_structure`, so the absence is
visible rather than silently missing.

**The one real anchor is VE-AC**, the tocopheryl-acrylate intermediate, which
Figure S1b **does** draw in full and for which Figure S2(a) prints an exact
mass. It is the right fragment to pin because it appears **twice in every
library member**: a wrong VE-AC is a wrong corpus, and nothing else in this
paper could detect it.

| | |
|---|---|
| anchor literal | `CC(C)CCCC(C)CCCC(C)CCC[C@]1(C)CCc2c(C)c(OC(=O)C=C)c(C)c(C)c2O1` |
| canonical | `C=CC(=O)Oc1c(C)c(C)c2c(c1C)CC[C@@](C)(CCCC(C)CCCC(C)CCCC(C)C)O2` |
| formula | C32H52O3 |
| InChIKey | `KODIAYXWITYQQW-YFCIGXFPSA-N` |
| [M+Na]⁺ calc | 507.3809 — Figure S2(a) prints calc **507.3809**, found 507.3804. **Δ = −0.03 mDa** |

**Independence is by provenance, not by byte-difference.** The canonical form of
this anchor is byte-identical to `stage03b_structures.py`'s VE-AC. That is
expected and is **not** the signal — F-037 is precisely the lesson that a
correct re-derivation canonicalises to the same string, so testing byte-identity
is a trap. Independence here rests on four things that were executed, not
asserted:

1. **`make_reference.py` imports no assembler.** Checked on the **AST**, not by
   grepping text: the module's imports are `ast, csv, hashlib, os, rdkit, sys`
   and it makes no call to `RWMol`, `AddBond`, `RemoveAtom`, `CombineMols`,
   `SetChiralTag`, `esterify` or `attach`. The enumerator
   (`stage03b_structures.py`) builds VE-AC through `curation/lib/assemble.esterify`.
2. **The anchor is a source string literal**, typed from the drawing.
3. **It is written in the reversed traversal.** The enumerator writes tocopherol
   from the aromatic 5-methyl inward and closes the pyran last
   (`Cc1c(C)c2c(c(C)c1O)CC[C@@](C)(…)O2`). The anchor starts at the far end of
   the phytyl tail, walks inward to the C2 quaternary centre, and closes the
   **aromatic** ring last. The acrylate is vinyl-first here, acid-first there.
4. **Direction is checked explicitly, because mass is blind to it.**
   `check_direction()` asserts the acryloyl is an **aryl** ester written
   `C(=O)O-aryl` (1 match), the chromanol oxygen is **in the ring** (1 match),
   the vinyl is terminal (1 match), and **no free phenol survives** (0 matches).

**The gate was made to fail on purpose.** Seed `anchor_direction_gate` in
`07_qc.json` rewrites the ester backwards — `Ar-O-C(=O)-vinyl` → `Ar-C(=O)-O-vinyl`
— which is **isomeric**: formula C32H52O3 both ways, exact mass 484.3916 both
ways, **Δ = 0.000000**. A patched copy of `make_reference.py` run as its own
process **exits 1** with *"the acrylate is not written as an aryl ester (0
matches). A reversed O-C(=O) has the same formula and the same exact mass."*
Verdict read from `returncode`, not from stdout.

**Stereo caveat on the anchor.** Only tocopherol's C2 is assigned (2*R*). The
4′ and 8′ side-chain centres are left **unassigned on both sides**, because the
paper names no isomer and reports no optical rotation; assigning RRR-α-tocopherol
would be invention. So C10 must compare these on a stereo-tolerant basis for
those two centres while C2 is asserted. Recorded in the anchor row's note.

---

## 5. Chemistry — what is known and what is not

Derived in `stage01_amines.py` (reads Fig S1a's 35 CAS labels and Fig S2's 8
HRMS pairs out of the text layer with subscript x-positions resolved against
the baseline glyphs) and built in `stage03b_structures.py`.

* The scaffold is **R-NH₂ + 2 × VE-acrylate → R-N(CH₂CH₂CO-O-tocopheryl)₂**, drawn
  in Fig S1b. Every VNxCyOz has exactly **two** arms, which the O₆ (or O₈) count
  forces.
* Core amine formulas follow exactly from the printed ion formula minus
  2 × C32H52O3.
* **Three structures are PROPOSED** (VN2C2, VN3C2, VN3C1), reproducing the
  printed `[M+H]⁺` calc to **−0.08, −0.03 and −0.03 mDa**.
* **Four are not** (VN2C4, VN3C2O, VN2C1 and the disputed panel (g)).
* **No structure is pinned.** See `DECISIONS.md` **D-1**. The decisive fact is
  that the cores of **VN2C1 — the paper's own lead compound — and VN3C2O have no
  amine of the right formula among the 34 distinct drawn amines at all**, so
  "the only drawn amine with this formula" is an argument whose premise the
  paper itself falsifies.

**An unprinted test the paper's own numbers pass 3/3:** the calc masses put
VN2C1→VN3C1, VN2C2→VN3C2 and VN2C4→panel (g) **exactly one CH₂ apart, to
≤ 0.05 mDa**. The paper never states that VN2→VN3 is a one-carbon step.

---

## 6. Verification state

`07_qc.json` holds 10 checks and **10 seeded defects, 10/10 caught**, each seed
in its own process with the verdict read from `returncode`. Nothing is credited
on a clean run alone.

`08_contradictions.csv` holds 7 entries — 4 OPEN against the paper, 3 RETIRED,
two of which were **my own** defects (benchmark SMILES written from memory; the
inset table's rule read as bar LNP43's top) kept with the measurement that
retired them.

**Not run, by instruction:** `provenance/structure_check.py`,
`check_anchor_independence.py`, `check_anchor_coverage.py`,
`provenance/build_dataset.py`, `audit.py`. Nothing outside this folder was
written. Nothing was merged into `data/`.
