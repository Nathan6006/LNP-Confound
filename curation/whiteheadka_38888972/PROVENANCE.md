# whiteheadka_38888972 — provenance

| | |
|---|---|
| Paper | Strelkova Petersen DM, Weiss RM, Hajj KA, Yerneni SS, Chaudhary N, Newby AN, Arral ML, Whitehead KA. *Branched-Tail Lipid Nanoparticles for Intravenous mRNA Delivery to Lung Immune, Endothelial, and Alveolar Cells in Mice.* Adv Healthc Mater 2024 |
| DOI | 10.1002/adhm.202400225 |
| PMID | **38888972** (indexed; no `esummary_cache` claim of unindexed status is made or needed) |
| PMCID | PMC11368637 — NIH author manuscript `nihms-2006250` |
| Proposed `Experiment_ID` | **`KW_2024`** — not settled, see `DECISIONS.md` **D-7** |
| Library | **580** = 29 amine heads × 20 acrylate tails, verified three ways |
| Rows emittable | **49** (in vivo). 531 members have no delivery measurement |
| Value source | **tier 5, pixels**, for every value — `00_value_sources.json` records why tiers 1–4 do not exist |
| Merge status | **BLOCKED**: 9 OPEN entries in `DECISIONS.md`; `build_dataset.py` refuses the folder |
| Licence | **not CC** — see `DECISIONS.md` **D-6** |

## The one number this paper is about

The abstract says *"a library of 580 ionizable lipidoids"*, of which *"over 40
enabled protein expression"*. Both are true and they are about different sets.

- **580** is exact and was counted three independent ways: Table S1's 29 rows ×
  Table S2's 20 rows; the Figure 2A heatmap's grid, measured as 19 internal
  column boundaries and 28 internal row boundaries fitted to a uniform grid
  (residual 1.10 px and 1.07 px); and 580 distinct canonical SMILES enumerated
  in `03_enumerated.csv`, with 580 distinct InChIKeys.
- **49** is how many ever went into an animal. The other **531** have a TNS
  fluorescence value and nothing else.
- **"over 40"** is of the 49, not of the 580: 49 of 49 clear the untreated bar,
  48 clear twice it, 44 clear 1e6 p/s.

So the extractable *chemistry* here is 580 structures and the extractable
*delivery data* is 49 rows. The intake gate says so before the first digitiser
runs, rather than at the end.

## Where the values come from

`00_value_sources.json` is the record. In short: Advanced Healthcare Materials
published **no Source Data file**; the sole SI item is a `.docx` whose every
live table is Table S1 (heads), Table S2 (tails) or Table S3 (antibody panels),
none of which carries a measurement; the JATS full text names **no repository
and no accession**; and the figures are flattened 300 dpi rasters with no text
layer, so there is no `pdfplumber` route either. Each of those four statements
is checkable against the files in `source/`.

**"Figure S2. Raw GFP transfection efficacy data" is not the screen.** Its own
caption says it is the pre-normalisation flow-cytometry data for `514O6,10`
formulated three ways — one lipidoid. The per-lipidoid values live only in
main-article Figure 2, which was not in `source/` at all until it was fetched
from `pmc.ncbi.nlm.nih.gov/articles/instance/11368637/bin/` (see
`00_acquisition.json`).

## The structural claim that mass cannot see

Each lipidoid is the Michael adduct of one amine head and one acrylate tail, so
**the number of tails is the head's N-H count**. Molecular formula and exact
mass are *consequences* of that count, not evidence for it — PIPELINE_V2's
"mass cannot see connectivity", in its purest form.

The paper pins the rule in one place only: Figure 1A draws `306 + 4 Oi10 →
306Oi10` with all four tails, on a head with exactly 4 N-H. Measured over all
29 heads, 15 have 4 N-H, 9 have fewer and **5 have more** (174 and 200 have 5,
110 and 122 have 6, 111 has 7) — and Methods says only 4 equivalents of
acrylate were used. Those 100 products cannot be what the paper names them.
They are emitted fully substituted, following the corpus's own precedent for
this chemistry (`AA_2008`, `KL_2010`, `KW_2014` all store full substitution),
with `confidence = low` and `structure_assumed = 1`, and the question is
queued as **D-2**.

## The paper's own cross-checks

Full numbers in `07_qc.json`. **Passing these is necessary and not sufficient**
— `QL_2026` reproduced every printed hit rate on bars up to 28 % too high — so
the geometric checks are listed first.

| | check | result |
|---|---|---|
| P3 | the **drawn** 10⁸ "nominally high potency" rule, a round number the authors drew and never printed as a coordinate | **0.99748 × nominal, −0.0011 decades** — a tenth of a pixel |
| P4 | Figure 2C's point cloud vs the (Fig 2A TNS, Fig 2B flux) pairs it must be, matched one-to-one | **30 markers**: mean Δ −80 AU on x (1.2 px), **+0.0094 decades** on y (0.8 px) |
| P5 | Figure 3A's stacked totals vs Figure 2B's bars — the same number drawn twice by the authors | **44 compared**, median −0.004, **rms 0.027 decades**; 6 positions refused, not filled in |
| P1 | library size 580 | **exact, three ways** |
| P9 | "an additional 13 LNPs over 10⁸ p/s" | **exactly 13**, on mean+SD |
| P8 | "402O8,12 and 500O6,10 … over 10⁹ p/s" | the three named lipidoids and no fourth, on mean+SD |
| P7 | "over 40 enabled protein expression" | 49 / 48 / 44 under three readings |
| P10 | "O14–O18 lowest TNS, O10/Oi10/O8,12 highest" | O18, O16, O14 are ranks 1, 2, 4 of 20; Oi10 and O10 are 20 and 19; **O8,12 is 15** |
| P6 | "514O6,10 produced 35 % lung, 35 % spleen, 30 % liver" | **FAILS: 40.6 / 34.5 / 24.9** by angular span, 0 unresolved rays, sum 100.000 % |
| P2 | Figures 2B and 3A show "the same 49 LNPs" | **FAILS: 51 vs 50 detected x ticks**; `306O12` is missing from 3A |
| P11 | "~75 % of the most efficacious LNPs … branched, unsaturated tails" | not checkable as written — branched is 90–92 % at every cut, *branched **and** unsaturated* is one lipidoid |
| P12 | "29 amine heads (80 through 602)" | count right, range wrong — Table S1 starts at **60** |

P8 and P9 together are the useful finding: two independent printed counts (3 at
10⁹, 13 at 10⁸) both land **exactly** on the upper-error-bar reading and
neither lands on the bar-top reading. The paper's potency thresholds are
mean+SD thresholds while the corpus stores means. The whisker detector that
establishes this is deliberately separate from `bars.bar_top` and its output is
never a value — error class 32 is a bar-top detector reading a whisker.

## Contradictions

`08_contradictions.csv`, **15 entries**. The three that change structures:

- **C-1/C-2/C-3** — heads **113**, **200** and **303** are drawn in Figure 1B
  as different molecules from the ones Table S1 names, each one methylene
  short; in 200 and 303 the drawn structure is an **aminal** (N–CH₂–N), which
  is not what Enamine or Sigma sells. For 113 the corpus settles it
  independently: `KL_2010`'s head `113` is `CN(CCN)CCN`, the name's structure.
  The names are emitted and the question is queued as **D-1**, because this is
  the `AC_2025` shape — a name and a structure disagreeing — and that is not
  ours to decide.
- **C-11** — the 1:4 stoichiometry versus five heads with 5–7 N-H (**D-2**).
- **C-13/C-14** — `Oi10` is a commercial branched-isomer mixture drawn as one
  isomer; `X8` is the only X tail whose name omits cis/trans while the drawing
  shows trans. Both flagged `structure_assumed = 1` on the row.

**C-15 is recorded although it is not a contradiction**: head `516`'s but-2-ene
is drawn with a **wavy bond**, an explicit declaration that the geometry is
undefined, and its name gives no E/Z either. The stored structure has no
geometry, which is exactly what both loci say, so the row carries no caveat
flag. Without this entry an unassigned C=C would read as an omission.

## Reproducibility

`sh run_all.sh` runs every stage in dependency order, `python3 -I` throughout
(the files in `source/` are untrusted input). A clean run from deleted outputs
takes **7.6 s** and reproduces every artifact **byte-identically**: all nine of
`01_blocks.csv`, `03_enumerated.csv`, `04_fig2b_bars.csv`, `05_values.csv`,
`06_tns.csv`, `07_fig2c_points.csv`, `08_pies.csv`, `07_qc.json` and
`reference_structures.csv` match sha256.

## The C10 anchor, and its teeth

`make_reference.py` builds **306Oi10** from Figure 1A, which draws the entire
reaction with every bond — the only place in the paper where the substitution
*count* is drawn rather than implied. The route shares nothing with the
enumerator: four short fragment strings substituted into a tail-first template
by `str.replace`, **no RWMol, no AddBond, no RemoveAtom, no CombineMols, no
reaction**, with RDKit used only to parse and canonicalise at the end. It
agrees with `03_enumerated.csv` by full InChIKey
(`BPPZRZOKGADTQE-UHFFFAOYSA-N`).

C16 was **made to fail on purpose** on this folder: pasting the enumeration's
canonical SMILES into `make_reference.py` as a literal made
`provenance/check_anchor_independence.py` exit 1, naming
`whiteheadka_38888972/make_reference.py:SEEDED_PASTE`; reverting returned it to
exit 0. `09_seeded_checks.py` does the same for this folder's own assertions —
**8 of 8** seeded defects are caught, including a positive control that must
*not* raise.

## What could not be verified

- The three disputed head structures (**D-1**) — no HRMS, no NMR, no elemental
  analysis anywhere in the paper or SI. Nothing internal decides it.
- The substitution state of the 100 lipidoids built on heads with >4 N-H
  (**D-2**) — no mass spectrum for any of them.
- The branch of `Oi10` (**C-13**) — a commercial isomer mixture.
- The in-vitro HeLa values (**D-4**) — they exist, as an unlabelled scatter.
- `Dose_ug_nucleicacid` (**D-5**) — 0.75 mg/kg with no body weight.
- The licence (**D-6**) — the Wiley landing page was not fetched.
