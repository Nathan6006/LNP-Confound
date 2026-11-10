# DECISIONS — `zhangn_42569578` → proposed `MY_2026` (PMID 42569578)

Queued, not resolved. `build_dataset.py` refuses a folder with an `OPEN` entry,
which is the correct state for this paper: **two of the three open items each
independently block the merge**, and one of them decides which molecule 12
measured values belong to.

---

### D-1 Which binding is right when the workbook and the figures disagree about 12 formulations?

`mmc2.xlsx` and SI Table S2 agree **perfectly** on which ionizable lipid has
which composition — a 25/25 composition bijection. They disagree only on the
**numbering**: the workbook orders the blocks ALC-0315(1-7), CKK-E12(8-13),
MC-3(14-19), SM-102(20-25); Table S2 and every figure x-axis order them
ALC-0315(1-7), SM-102(8-13), MC-3(14-19), CKK-E12(20-25). It is a clean
two-block swap with within-block order preserved, so 13 formulations sit at the
same index and **12 are permuted**.

That alone would be harmless bookkeeping. It is not harmless because Figs 4c
and 5c plot the workbook **in sheet order** while labelling their x-axis groups
**in Table S2 order**. Verified: every one of the 25 bars in Fig 5c and Fig 4c
matches the workbook row of the same number. So the two sources bind the lipid
name to different things, and for 12 formulations they disagree about which
molecule a measured value belongs to.

**Evidence that the name travels with the COMPOSITION** (workbook + Table S2):

- The 25/25 composition bijection itself: the lipid name is attached to the
  same composition in both documents, independently.
- The ILC discriminator, derived from the paper's own Fig 3a and its own stated
  selection rule, and using none of the disputed labels. The paper reports that
  ionizable-lipid concentration correlates **positively** with size for SM-102
  (ρ = 0.64, the strongest of the four) and **negatively** for CKK-E12, and the
  screen prioritised **small** size. So the selected CKK-E12 set should sit at
  **high** concentration and the selected SM-102 set at **low**. Measured: the
  CKK-E12-named block has mean ILC **8.54 mM with 5 of 6 at the 10 mM
  maximum**; the SM-102-named block has mean ILC **3.13 mM with 0 of 6 at
  10 mM**. This matches the workbook naming.

**Evidence that the name travels with the VALUE** (figures + abstract + text):

- The abstract gives two values to 4 dp — "SM-102- and CKK-E12-based
  formulations … reduced to 0.0744 and 0.0611, respectively". Both appear
  exactly in column V. Under the figure-label reading the pairing is exactly
  right; under the workbook reading it is reversed.
- "SM-102-based LNPs exhibited the highest uptake efficiency … MC3-based the
  lowest" (§3.3). Measured block means for high-Cy3 uptake: the block the
  figures call SM-102 is **1.441e7**, the highest; the one they call CKK-E12 is
  6.259e6; ALC-0315 4.480e6; MC-3 4.223e6, the lowest. Matches the figure
  labels.
- "CKK-E12 … demonstrated the strongest silencing effect" (§3.4) is true under
  the figure-label reading and false under the workbook reading.

Both sides are real and neither is clerical. This is the `AC_2025` shape — a
name and a structure disagreeing — where the first repair was confident,
internally consistent, and would have put **419 measured values on the wrong
molecules**. `PIPELINE_V2` §3 names this case explicitly as one to queue.

Resolving it needs the authors, or the `merged_quantified_results_with_SizePDI.xlsx`
that the repo references by local path and does not commit.

Note the two readings are **not symmetric in cost**: under either one the 13
ALC-0315 and MC-3 rows are unaffected, but the 12 SM-102/CKK-E12 rows swap
between two lipids whose corpus footprints differ by an order of magnitude
(SM-102 2,312 rows across 6 experiments; cKK-E12 241 across 3).

STATUS: OPEN

---

### D-2 Does a source whose shared chemistry is *only* universal benchmarks belong in the concordance pool?

This paper contributes **zero new ionizable lipids**. All four — ALC-0315,
SM-102, DLin-MC3-DMA and cKK-E12 — are commercial catalogue reagents, bought
from Avanti and MedChemExpress, and **all four are already in the corpus**
(313 / 2,312 / 3,098 / 241 rows by full InChIKey).

`CLAUDE.md` §8 already carries this question for `XD_2026`, where *two* of the
shared lipids were benchmarks: "a source whose only shared lipids are benchmarks
that top every screen can push cross-lab concordance toward +1 for reasons that
have nothing to do with transferable SAR." Here the fraction is **4 of 4**. If
merged, `MY_2026` would enter Stage 5's shared-lipid pool contributing nothing
but universal positive controls, against in-vitro experiments that carry the
same four (`AC_2025`, `WT_2026`, `YZ_2022`, `YZ_2024`, `LZ_2025`, `SP_2020`,
`JW_2024`).

This is the §8 decision in its purest available form, and it should be settled
on the §8 ledger rather than inside this folder.

STATUS: OPEN

---

### D-3 Is a 4-lipid formulation screen a chemistry source at all?

F-067 records that five one-lipid formulation screens accounted for ~74% of
staged rows, and that a policy is needed. This is the same class with 4 lipids
instead of 1: the 25 rows vary helper lipid, three molar ratios and
concentration while the **chemistry is held fixed by construction** — exactly
4 distinct structures over 25 rows.

Precedent exists both ways: `ZT_2026` is a 486-formulation DoE on **one** lipid
and was merged. So the question is not whether such a source can be merged but
whether it should count toward a corpus whose unit of analysis is the ionizable
lipid.

STATUS: OPEN

---

### D-4 Which transform carries a knockdown readout onto the corpus scale?

`Experiment_value` is a per-experiment z-score of log10 with population sd, and
for delivery **higher must mean better**. This paper's endpoint is *relative*
S100P expression, where **lower** is better. The corpus's one in-vitro knockdown
precedent, `LL_2012` (`LRP6_knockdown_normalized`, 112 rows), stores the
knockdown so that higher = better, confirmed: mean 0.000000, sd(ddof=0)
1.000000, floor-heavy with its minimum repeated.

The natural monotone choice is `z(-log10(relative expression))`, which is
exactly "z of log10 of a delivery-proportional signal". It is computed and
recorded in `07_qc.json` (`would_be_z`) so the choice is visible and reversible,
but it is a convention, not a measurement, and no row was emitted on it.

A new `Experiment_method` term would also be needed —
`S100P_knockdown_normalized`, by analogy with `LRP6_knockdown_normalized`.

STATUS: OPEN

---

### D-5 The per-well siRNA dose for the gene-silencing arm is not stated

§2.4 states the uptake doses precisely (9 pmol/well high, 3 pmol/well low).
§2.5 states only that S100P siRNA was diluted to 2 nmol/mL and that 250 µL of
a fivefold-diluted LNP preparation went onto each well of a 12-well plate. The
delivered pmol depends on the 3:1 aqueous:organic ratio and on an encapsulation
efficiency the paper never reports.

A number could be arithmetically constructed. It would be a guess, so
`Dose_ug_nucleicacid` stays blank and `05_values.csv` records
`dose = "not stated for this arm"` on all 25 admissible values.

STATUS: OPEN
