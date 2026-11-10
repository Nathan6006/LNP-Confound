# KK_2015 — decisions queued for the user

Kauffman KJ et al., *Nano Letters* 15(11):7300–7306 (2015), PMID 26469188,
DOI 10.1021/acs.nanolett.5b02497. Folder named for the PI (`andersondg_26469188`);
the corpus key derives from the **first author**, so it is `KK_2015`.

Nothing below is resolved by this curation. `PIPELINE_V2.md` §3: *the pipeline
detects that a source is silent; it does not get to decide what the source meant.*

---

### D-1 — Admit a one-ionizable-lipid formulation screen at all? (F-067)

This paper contributes **exactly one ionizable lipid: C12-200.** Every one of
its 38 library formulations plus the `Original` benchmark uses it. The factors
that vary are the phospholipid identity (DSPC / DOPE / DOPC, and a DSPE arm that
never formulated), the four molar percentages, and the C12-200:mRNA weight ratio
— all of which are existing corpus columns (`HL_name`, `IL_molratio`,
`HL_molratio`, `CHL_molratio`, `PEG_molratio`, `IL_to_nucleicacid_massratio`),
so the screen is fully representable. 36 rows are staged in `06_rows.csv`.

For a corpus keyed on ionizable-lipid chemistry the contribution is **1 lipid,
however many rows exist.** Against that:

* the 36 arms are **individually measured** (n=3 mice each, printed values with SD),
  not a pooled or modelled summary;
* **chemistry does vary across arms** — three distinct phospholipids — just not
  the ionizable lipid;
* there is direct precedent: **`AC_2025` is itself a C12-200 formulation screen**
  and carries 498 in-vitro rows under this same structure;
* it would be the corpus's **first in-vivo C12-200**, bridging a 518-row
  in-vitro structure into the in-vivo pool.

And the honest cost: 36 in-vivo rows is **2.0 % of the current 1,803-row in-vivo
pool** for one lipid, which is exactly the dilution F-067 was opened about.

**Recommendation:** admit, *or* admit only the 6-point Library C weight-ratio
ladder plus `Original` (7 rows), *or* decline. Not a call this curation should make.

STATUS: OPEN

---

### D-2 — A-03's cholesterol: emit the printed 31.5, or the over-determined 21.5?

`06_rows.csv` currently carries **21.5**, with the reason on the row in
`curation_note` and in `08_contradictions.csv` C-01.

The printed 31.5 makes A-03's four components sum to **110.0 mol%**, which is
physically impossible. 21.5 is forced three ways: 13 of 14 Library A rows have
printed cholesterol *exactly* equal to `100 − (IL + phospholipid + PEG)`;
A-03's other three values are triple-attested (Table S1, the p10 design table's
numeric columns, and the 174-assertion coded-level bijection); and **A-12
carries the identical (60, 16, 2.5) levels and prints 21.5.**

This is not a guess among candidates — it is the unique solution — but it
**changes a published number**, which `PIPELINE_V2.md` §3 says to queue.
A-03's measured EPO value (176 ± 21 ng/mL) is unaffected either way.

STATUS: OPEN

---

### D-3 — The PEG-lipid is never identified

The Table S1/S3/S4 column header is the bare word **`PEG`**. The SI never names
the PEG-lipid, and the main text — which would — is closed-access and not held.
Emitted as `PEG_name = "Unknown"` with a blank `PEG_SMILES`, which is an
existing corpus convention (1,798 rows across `KZ_2016` and `RG_2023`).

Guessing was **specifically refused**. A 2015 Anderson-lab C12-200 formulation
would conventionally use C14-PEG2000, but the briefing records that
`C14-PEG2000` is drawn in the literature as a phosphoethanolamine-carbamate-PEG
and stored in *this* corpus as the glycerol ether — so a plausible guess puts a
wrong structure on all 36 rows. If the user can supply the main text, this
resolves immediately.

STATUS: OPEN

---

### D-4 — Licence: this source fails §8

Off the Europe PMC core record, not the landing page: `license = null`,
`isOpenAccess = N`, `fullTextUrlList availabilityCode = S` ("Subscription
required"), and `authMan = epmcAuthMan = nihAuthMan = N`, so there is no PMC
author manuscript either. The article is **All Rights Reserved, © 2015 American
Chemical Society.** Funders recorded are Shire and Skoltech/MIT with no NIH
grant, consistent with there being no public-access deposit. No repository
exists, so there is no second licence to assess.

`CLAUDE.md` §8 requires every constituent source to be CC0/CC-BY-compatible
before corpus release. **This one is not.** If admitted, only derived numerical
values could be redistributed, never the SI PDF, and §8 needs an explicit ruling
on whether a closed-access source may contribute values at all.

STATUS: OPEN

---

### D-5 — `confidence = medium` or `low` for a structure this paper never attests?

All 36 rows carry `confidence = medium` and `structure_assumed = 1`. The paper
provides **no structural attestation of C12-200 whatsoever**: no drawing, no
molecular formula, no HRMS, no molecular weight, no synthesis. The structure is
carried over from the corpus, where `AC_2025` (498 rows) and `JW_2024` (20 rows,
under the name `other`) agree on InChIKey `RVHYPUORVDKRTM-UHFFFAOYSA-N`.

`low` is defined as "an element of the structure is an assumption the paper does
not pin", and here *the whole structure* is unpinned by this paper — which
argues for `low`. Against: C12-200 is a named catalogue compound from this same
lab, and the amine-200 core is independently corroborated by three corpus
experiments and by an arm-count derivation (`stage04_structures.py`). `medium`
was chosen, with the reason stated verbatim on every row.

STATUS: OPEN

---

### D-6 — `Cargo_type`: `EPO` or the corpus's `hEPO`?

Emitted as **`EPO`**. The SI says only "Serum EPO" and "erythropoietin"; it
never states the species, and the main text is not held. The corpus token for
the two existing erythropoietin experiments (`RG_2023`, `LL_2025`) is **`hEPO`**,
which asserts human. Emitting `EPO` is honest but will not join with those two
experiments; emitting `hEPO` joins but asserts something the source does not say.

STATUS: OPEN

---

### D-7 — `Model_type`: the mouse strain is not stated

Emitted as **`Mouse_unspecified`**, a new token following the corpus's
`Mouse_<strain>` pattern. The Table S1 footnote says only "into mice"; no strain
appears anywhere in the SI. `Model_type` is never blank anywhere in the current
corpus (0 of 30,998 rows), so a blank was not an option. If the main text names
the strain this becomes `Mouse_B6` or similar.

STATUS: OPEN

---

### D-8 — Tables S3 and S4: 8 admissible values, not emitted

Both readouts are **admissible delivery endpoints** (luminescence; target-protein
knockdown) and are recorded in `05_values.csv` with `corpus_admissible = yes`.
Neither is emitted:

* **Table S3**, 2 rows, luciferase mRNA, total flux `(p/s)×10⁹`. A 2-point z is
  exactly ±1 by construction and carries no information; and pooling `p/s` with
  `ng/mL` into the EPO block's single z would standardise two incommensurable
  quantities across a different cargo.
* **Table S4**, 6 rows, FVII siRNA, serum FVII as % of PBS — a different cargo,
  a readout that is floored at 0 and inverted in direction, and only 2
  formulations measured at 3 doses each.

Together they would be 8 rows on **2 formulations of the same single ionizable
lipid**. Padding the row count was refused.

STATUS: OPEN

---

### D-9 — F-061: this paper cannot adjudicate it, and `data/` is untouched

The corpus holds two incompatible structures for this lab's catalogue amine 200:

| where | name | formula | N | arms (derived) |
|---|---|---|---|---|
| `AC_2025` (498 rows) | `C12-200` | C70H145N5O5 | 5 | 5 |
| `KW_2014` (10 rows) | `A200_C10…C14` | C75–C95 H…N5O10 | 5 | **5.0** |
| `NC_2024` (13 rows) | `A200_O*`, `A200_X*` | …N5O10 | 5 | 5 |
| `KL_2010` (4 rows) | `200C8/C10/C12/C14` | C20–C26 H…N3O | 3 | **1.0** |

The arm counts are *derived*, not read: within a series holding the amine fixed
and walking the tail length, `dC/d(tail)` **equals the number of arms**.
`KW_2014` gives a constant 5.0 across four steps; `KL_2010` gives a constant 1.0
across three. Both series are internally self-consistent, which is why this is a
real conflict rather than a typo — two curations each coherently encode a
different molecule for the same catalogue amine.

**This paper supports neither form by its own evidence.** It cites C12-200 by
name only; no structure, formula, HRMS, MW or synthesis appears, and the weight
ratio and molar percentages are independent recipe parameters from which no
molecular mass is derivable. Nothing in `data/` was changed. The `AC_2025`
precedent — a confident repair that put 419 measurements on the wrong molecules
— is exactly why.

STATUS: OPEN

---

### D-10 — `Model_target = liver`, but the measurement is serum

Emitted as **`liver`**, on the authority of the abstract's own words:
"optimize lipid nanoparticle formulations for mRNA delivery **to the liver** in
vivo". The measurement itself is **circulating serum EPO**, not an organ-resolved
readout, so the organ attribution is the authors' claim about where translation
occurs, not something this assay measures. Recorded so a consumer joining on
`Model_target` knows which it is.

STATUS: OPEN
