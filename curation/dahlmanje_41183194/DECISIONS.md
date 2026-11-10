# DECISIONS — `KG_2025` (Gentry 2025, PNAS, PMID 41183194)

Paper: *Glycolipid nanoparticles target the spleen and detarget the liver
without charge.* Gentry K, Lian L, …, Dahlman JE. Proc Natl Acad Sci U S A
122(45) e2409569122. DOI 10.1073/pnas.2409569122.

`Experiment_ID` = **`KG_2025`**, from the FIRST author (Kara Gentry). The folder
is named for the PI, which is correct as a folder name and is not the key.

Nothing here has been merged into `data/`. Entries marked `STATUS: OPEN` block a
merge by design.

---

### D-1. Does a source whose only chemical variable is the helper lipid belong in this corpus?

The screen contains exactly **two** ionizable lipids, SM-102 and cKK-E12, and
both are already in the corpus — SM-102 on **2,312 rows across 6 experiments**,
cKK-E12 on **241 rows across 3**. Every one of the ten rows emitted here carries
an ionizable lipid the corpus already has. What varies is the helper lipid, the
PEG-lipid and the molar ratio.

This is not disqualifying on precedent: `ZT_2026` (486 rows, **1** ionizable
lipid, 6 helper lipids), `YZ_2022` and `YZ_2024` (1,080 rows each, 1 IL, 6 HLs),
`LZ_2025` (96 rows, 1 IL, 8 PEGs), `LL_2025`, `SP_2020` and `SB_2024` are all
already in `data/` on the same shape, and `ZT_2026` is this project's own
flagship curation. So the corpus has repeatedly accepted single-IL,
excipient-varying experiments.

The question is whether it should, for THIS project specifically, whose Results
1 and 3 are about ionizable-lipid chemistry. A source contributing no new
ionizable chemistry cannot inform Result 1, and its contribution to Result 3 is
described in D-7.

STATUS: OPEN

---

### D-2. The pooled 109-LNP screen publishes ratios and unlabelled points, and no attributable per-LNP value.

The design supports a per-LNP normalized delivery value for each of 109
formulations in each of up to 25 sorted cell types. The publication attributes
**none** of them:

- **Fig. 1G/1H** aggregate normalized barcode count into TWO groups per library
  (C8 vs C16 glycolipid). A violin is a distribution; no LNP identity, molar
  ratio or PEG identity survives it.
- **Fig. S11** plots 52 points for the SM-102 pool, and they are disqualified
  twice over: the y quantity is `Relative Normalized BC Count (spleen/liver)`,
  a **ratio**, and every point is **unlabelled**.
- **Fig. S4**, captioned as the ranking of GlycoLNPs by spleen-to-liver ratio,
  is a cartoon containing no LNP name and no number at all (see
  `08_contradictions.csv`).

So the barcode screen contributes **zero rows**, and the reason is attribution,
not pooling. Pooling alone is not disqualifying — `LR_2023` contributes 407
in-vivo rows as `barcoded`. (`LR_2023` is Rhym/Anderson at MIT, PMID 37127709,
**not** a Dahlman-lab experiment; that was checked against Europe PMC rather
than assumed.)

The ratio question is the one that needs a ruling, because it recurs: `CLAUDE.md`
§8 already holds that a selectivity contrast is not a poolable readout, and this
paper's entire thesis is a spleen-to-liver contrast. Does a screen that
publishes only ratios contribute nothing, or should the corpus gain a way to
carry a contrast?

STATUS: OPEN

---

### D-3. Six of the ten emitted rows rest on a pixel read of Fig. 2B.

`00_value_sources.json` records the whole paper at tier 5 for anything not
printed: there is no publisher Source Data file, no author repository
(`find_repo.py`: 5 queries, 0 hits), and the only deposit is SRA PRJNA1150220,
which holds scRNA-seq reads and no delivery value. `fullTextPDF` returns HTTP
404, so even the main figures exist only as flattened 708×720 and 796×671 JPEGs.

Fig. 2B was digitised anyway because it is the only locus attributing an
absolute per-organ value to GlycoLNP45, 59 and 62, and because it carries its
own external anchor: its GlycoLNP61 bars measure the quantity Fig. 3C PRINTS as
54 and 78. Measured deltas **+0.67** and **−0.16** percentage points, against a
panel resolution of 0.5666 and 0.5682 percent per pixel. Those two numbers took
no part in the calibration or the detection.

Accepted, with the geometry gated: scale and origin fitted from detected ticks
only (LOO residual 0.5 px liver, 0.0 px spleen, against the **strict** 1.5 px
tolerance — `axis.native_tol` was not used and no widening is relied on); the
zero tick, excluded from the fit because it is drawn fainter than the detector
threshold, is predicted at 208.25 px against a drawn baseline at 208 and a faint
tick measured at 208.333; and the two panels, calibrated from disjoint tick sets
in disjoint pixel columns, give px-per-point figures agreeing to a ratio of
1.00284.

STATUS: SETTLED — accepted, and the rows carry `confidence = medium`.

---

### D-4. Fig. 3A/3C and Fig. 3D disagree, in opposite directions, on the same condition.

Same formulations, same 0.8 mg/kg, same organ and cell type, same 16 h:

| | Fig. 3C printed | Fig. 3D digitised | ratio |
|---|---|---|---|
| GlycoLNP61, splenic macrophages | 78 | 88.03 | ×1.129 |
| Control LNP61, splenic macrophages | 66 | 46.37 | ×0.703 |

The paper's own effect size, Glyco ÷ Control, is **1.1818** from Fig. 3C and
**1.8984** from Fig. 3D — a ratio of ratios of **1.6064**.

The digitisation is not the problem: Fig. 3D's ticks sit at 384.032 / 480.0 /
576.077 px with a LOO residual of **0.109 px**, every bar passed `bars.bar_top`
with a flank delta of 0.0, and the panel's own dose response is strictly
decreasing (88.03 > 77.10 > 26.59) — a consistency test nothing in this curation
imposed. The cohorts differ (Fig. 3A shows ~12 points per bar, Fig. 3D shows 3).

The printed Fig. 3C values are emitted because printed numerals outrank a bar
read, and the two were **not** averaged. But which cohort is the better estimate
of this formulation's delivery is the paper's question, not ours, and the answer
moves the effect size by 1.6×.

STATUS: OPEN

---

### D-5. `Control LNP61`'s liver value has two printed candidates, 81 and 82.

Fig. 3C's pentagon prints `Li (54,81)`; the Results text says 82%. The figure
legend states that the parenthesised values are the percent delivery, which
makes Fig. 3C the primary locus, and that is what is emitted. On this
experiment's scale the two candidates differ by 0.051 z. The row carries
`confidence = medium` and names the conflict in `curation_note`; they were not
averaged.

STATUS: OPEN

---

### D-6. The glycolipid's sphingoid base length is not pinned by the raster.

`C16 Galactosyl(α) Ceramide` is **absent from the corpus** under every name
tried, so its structure had to be established rather than looked up.

What IS established: the N-acyl chain measures **16 carbons** off Fig. 1A, with
all three independent methods agreeing (turns 15, mid-line crossings 15, apex
vertices 15 bonds; +1 for the carbonyl carbon), and `chains.check_name_token`
confirms the token "C16" means the acyl chain — the `KS_2024` `C24` trap tested
and passed. The galactose is α and the base carries one trans double bond, both
drawn.

What is NOT: the sphingoid base's length. Mid-line crossing gives 14 bonds in
both panels; apex counting gives 14 for the C8 panel and 13 for the C16 panel.
At ~9 px per bond on a flattened JPEG a ±1 ambiguity is real, and
`chains.count_carbons`' contract is to refuse rather than average or prefer.

The base is therefore taken from the reagent's chemical name (d18:1
sphingosine), and **all eight rows carrying this helper lipid are flagged
`structure_assumed = 1`** with the reason in `curation_note`. The flag is on the
row, not in this file (F-001).

STATUS: OPEN

---

### D-7. What `KG_2025` would actually add to the in-vivo shared-lipid set — measured, not assumed.

`CLAUDE.md` §8 carries an open decision about benchmark-only pairs. A sibling
curation raised the same question from the other side (**F-064**: merging
`dahlmanje_36701517` would give the in-vivo pool a **third** shared chemistry,
cKK-E12, and unlike MC3 and SM-102 it is *not* a universal benchmark). So I
measured what `KG_2025` contributes rather than arguing it.

**SM-102 — a fifth source, with NO organ-matched pair.** SM-102 sits in four
in-vivo experiments today, and their targets are:

| source | `Model_target` |
|---|---|
| `SB_2024` | ear (13 rows) |
| `SY_2025` | ear (2) |
| `SS_2018` | whole_body (1) |
| `XD_2026` | whole_body (1) |

`KG_2025`'s six SM-102 rows are **liver and spleen**. The intersection with
every existing SM-102 source is **empty**. So merging raises SM-102's source
count from 4 to 5 while contributing **zero** organ-matched cross-source points
— it inflates the shared-lipid *count* without adding a comparable pair.

**cKK-E12 — in the live corpus, also no organ-matched pair.** cKK-E12 is in the
vivo pool once, `BL_2023`, 6 rows under the name `MD-1`, targets **muscle** and
**lung_epithelium**. `KG_2025`'s four cKK-E12 rows are liver and spleen.
Intersection again **empty**.

**But against the staged sibling, the pair IS organ-matched — and that is the
problem.** `curation/dahlmanje_36701517/06_rows.csv` stages 37 rows, all
cKK-E12, with targets liver 13, lung_epithelium 9, kidney 6, heart 5 and
**spleen 4**. If both folders merge, cKK-E12 becomes a shared lipid between
`KG_2025` and `dahlmanje_36701517` on **both liver and spleen**, matched on
organ, on a non-benchmark research lipid — exactly the rare and valuable
cross-source point F-064 argues for.

Except that **both papers are from the Dahlman laboratory.** `Experiment_ID` is
the project's environment key, so two papers from one lab are two environments
and the pair would be scored as cross-source agreement. The manuscript's claim
is about cross-**laboratory** prediction. See D-11.

STATUS: OPEN

---

### D-8. The lung, kidney and heart values are printed but not emitted.

Fig. 3C prints ten more numerals: L (13,0), K (1,1), H (0,0). Four of the six
are **0%**, which cannot enter the corpus's per-experiment z-score of log10, and
emitting one half of a pair would put a one-sided organ into the
standardisation. The two kidney values are printed as `1` at the very bottom of
the range, where 1 versus 0.5 versus 1.4 moves log10 by 0.45 decades. No cell
type is stated for any of these three organs, where spleen and liver both name
one. All six are catalogued in `05_values.csv` with
`corpus_admissible = no` and the reason on the row.

STATUS: SETTLED — excluded, with every excluded value recorded.

---

### D-9. The article licence fails the §8 release requirement.

Europe PMC's `license` field for PMC12625924 is **`cc by-nc-nd`**, corroborated
by the JATS permissions block. `CLAUDE.md` §8 requires constituent-source
licences to be CC0/CC-BY-compatible before corpus release. **NonCommercial and
NoDerivatives both fail**, so this is two independent failures, not one. There
is no repository and therefore no second licence to consider.

STATUS: OPEN

---

### D-10. `GlycoLNP60` has a complete formulation and no delivery value.

It is characterised in Fig. S5 and Fig. 2A (132 nm, PDI 0.22, EE 18%) and the
text then excludes it from dosing on its encapsulation efficiency. It is carried
in `03_enumerated.csv` with `dosed = no` and appears in `05_values.csv` with
`corpus_admissible = no`, so the design/attribution gap is visible rather than
silently dropped.

STATUS: SETTLED — enumerated, not emitted.

---

### D-11. Two papers from one laboratory are two `Experiment_ID`s but one laboratory.

This curation and `curation/dahlmanje_36701517` are both first-authored out of
the Dahlman lab at Georgia Tech/Emory, and both are keyed on their own first
author, so they would enter the corpus as two independent sources. On the
measurements in D-7 they would also supply the in-vivo pool's only
**organ-matched, non-benchmark** shared-lipid pair.

That pair would be a within-laboratory, cross-paper comparison presented to
every downstream analysis as a cross-source one. For Result 3 — the load-bearing
result, which holds chemistry constant and varies the laboratory — this is a
distinct hazard from the benchmark one §8 already records, and it points the
opposite way: a benchmark-only pair can push concordance toward +1 for reasons
unrelated to transferable SAR, and a same-lab pair can do so for reasons
unrelated to *cross-laboratory* transfer specifically.

Nothing in the corpus currently records which laboratory an `Experiment_ID` came
from. Whether it should, and whether same-lab pairs should be excluded from
Result 3 or reported separately as a within-lab rung (the role `molea.csv`
already plays in §8), is not a call this curation can make.

STATUS: OPEN

---

### D-12. `KG_2025` is the sixth instance of F-067's one-lipid formulation screen, and it is already curated the way F-067 recommends.

F-067 records that five papers in the current run are one-lipid formulation
screens carrying ~74% of the run's staged rows, and recommends **option 3** as
the default: admit only the individually-measured validation arms, where the
chemistry actually varies, rather than the large pooled screen whose variance is
formulation.

`KG_2025` is a *two*-lipid formulation screen and lands on option 3 by
construction rather than by choice: the 109-formulation pooled screen
contributes **zero** rows because it publishes nothing attributable (D-2), and
the ten emitted rows are exactly the individually-dosed validation arm. If the
user adopts option 1 (exclude one-lipid screens) this folder is excluded too; if
option 3, it is already in the recommended shape.

STATUS: OPEN — governed by the F-067 policy decision, not separately decidable here.
