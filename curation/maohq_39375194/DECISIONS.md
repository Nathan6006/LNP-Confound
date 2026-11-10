# DECISIONS — maohq_39375194 (Cheng et al., ACS Nano 2024, PMID 39375194)

Queued, not resolved. `provenance/build_dataset.py` refuses a folder with an
`OPEN` entry, which is the intended state: none of these can be settled from
the paper, and the pipeline detects that a source is silent — it does not get
to decide what the source meant.

---

### D-1 — Licence: neither the article nor the repository is CC0/CC-BY compatible

The article is closed (Europe PMC `isOpenAccess = N`, `license = null`, DOI
marked "Subscription required"). The PMC deposit `PMC11512640` is an NIH author
manuscript whose `<permissions>` reads *"This file is available for text mining.
It may also be used consistent with the principles of fair use under the
copyright law."* — not a Creative Commons licence.

The GitHub repository `MaoResearchGroup/ml_lnp_design_rules` carries
**GPL-3.0**. That is a copyleft *software* licence; applied to `Raw_Data/` it is
not CC0/CC-BY and would propagate to any redistribution.

`CLAUDE.md` §8 requires every constituent source to be CC0/CC-BY-compatible
before corpus release. This source satisfies neither side. It is a **release**
blocker, not a curation blocker: the numbers can be used for analysis now, but
the "Data & code availability" section cannot claim a clean licence chain for
them.

Note the asymmetry worth deciding on generally: the already-merged `YZ_2022`
(Nat Commun, PMID 35879315) **is** open access and covers the same library's
HepG2 arm. So the licence-clean and licence-unclear versions of one experiment
are in different papers.

STATUS: OPEN

---

### D-2 — `18PG` is SOPG in one paper from this lab and DOPG in the other, for the same library

Measured, both ways, in `08_contradictions.csv` C-2:

| | says | evidence |
|---|---|---|
| ACS Nano 2024 (this paper) | **DOPG**, 18:1/18:1 | SI Table 2 draws a cis C=C in *both* tails; SI Table 3 prints `Double_bonds = 2`; `1080_LNP_OHE.csv` carries 2 on all 180 `18PG` rows |
| Nat Commun 2022 (PMID 35879315) | **SOPG**, 18:0/18:1 | names it verbatim: *"1-stearoyl-2-oleoyl-sn-glycero-3-phospho-(1'-rac-glycerol) (18PG)"* — describing **this same 1,080-formulation library** |

`provenance/reagent_registry.py` pins `18PG` to SOPG on the Nat Commun
authority, over **1,474** corpus rows in `WT_2026`, `YZ_2022` and `YZ_2024`.

What we emitted: the structure **this paper draws** (dioleoyl PG, 1'-rac head,
`DSNRWDQKZIEDDB-GCMPNPAFSA-M`) on all **972** `18PG` rows, with
`confidence = medium` and the conflict in `curation_note` on every one.

Why it is not repaired here. This is the `AC_2025` shape — a name and a
structure disagreeing — where the first confident repair was internally
consistent and assigned **419 measured values to the wrong molecules**. And the
blast radius is not symmetric:

- if Nat Commun is right, only **our 972 rows** are wrong;
- if this paper's drawing is right, **360 already-merged rows** (`YZ_2022` +
  `YZ_2024`, same library) carry the wrong helper lipid today, and `WT_2026`'s
  1,114 `18PG` rows are a separate question again because that is a different
  purchase.

Cheapest way to settle it: the Avanti catalogue number from either lab's order
records, or an author query. Avanti sells both (`840503` SOPG, `840475` DOPG),
both as the 1'-rac head, and both are plausible for a six-lipid charge panel.

Do not merge until this is answered. Whichever way it goes, `C17` will have
something to say.

STATUS: OPEN

---

### D-3 — The dose and timepoint for this screen are stated twice, differently

Same measurements, two Methods sections (`08_contradictions.csv` C-1):

| | per-well dose | incubation | final conc. | plasmid |
|---|---|---|---|---|
| ACS Nano 2024 | 4 µL into **100 µL** = **0.1 µg**/well | **48 h** | 1 µg/mL | Aldevron |
| Nat Commun 2022 | 8 µL into **200 µL** = **0.2 µg**/well | **24 h** | 1 µg/mL | Translate |

The corpus's `YZ_2022` rows already carry `Dose_ug_nucleicacid = 0.2`, correctly
for that paper. We emit **0.1 µg** and 48 h on the 4,680 training rows, because
the five cell lines we emit are described *only* by this paper — but the same
paragraph demonstrably misdescribes the HepG2 arm, whose values are provably
identical to `YZ_2022`'s, so it is not a clean authority for anything.

Three readings, none decidable from the papers:

1. The five new cell lines really were run at 0.1 µg / 48 h and HepG2 at
   0.2 µg / 24 h; the ACS Nano Methods retro-describes everything with the
   newer protocol.
2. Everything was run at the Nat Commun protocol and the ACS Nano Methods is
   wrong.
3. Everything was run at the ACS Nano protocol and the Nat Commun Methods is
   wrong — least likely, since `YZ_2022` is the primary report of that arm.

Under reading 1 the emitted dose is right. Under 2 it is wrong on 4,680 rows.
An author query settles it; a guess does not.

STATUS: OPEN

---

### D-4 — DOTAP: the paper draws a stereocentre the corpus does not store

SI Table 2 draws a bold wedge to H at DOTAP's propane C2. The reagent registry's
`DOTAP` has no stereocentre, and **2,564** corpus rows across six experiments
carry that unspecified structure.

We emitted the registry's unspecified form with `structure_assumed = 1` and the
reason on the row, because emitting a stereo-defined DOTAP would give one
reagent name two structures across experiments — exactly what `C17` gates.

The decision is whether the registry should gain the stereocentre (changing
2,564 existing rows) or the drawing should be recorded as a depiction
convention. This affects an excipient's stereochemistry only: no ionizable lipid
and no measured value moves either way.

STATUS: OPEN

---

### D-5 — `RLU_HepG2` is already in the corpus; confirm the exclusion

Not a judgement call about chemistry, but it changes 1,080 rows, so it is
recorded rather than assumed.

All 1,080 formulations join `YZ_2022` one-to-one on
(helper lipid, four molar percents, N/P). Fitting the affine map from our ln
values to the corpus z and expressing the residual back in **ln units**:

- rms residual **2.8531e-4**
- expected rms for ±0.0005 three-decimal rounding: **2.8868e-4**
- ratio **0.988**; max |residual| 5.33e-4

So the discrepancy *is* the repository file's rounding and nothing else, and the
corpus holds these values at **higher** precision than the repository. We
therefore exclude all 1,080 `RLU_HepG2` values from `06_rows.csv`
(`corpus_admissible = no` in `05_values.csv`, with that reason).

What should be confirmed is the *consequence*: `LC_2024` as emitted is five cell
lines, and the sixth lives under a different `Experiment_ID` from a different
paper. Any analysis that treats `Experiment_ID` as "a lab's screen" will see
this one screen as two environments. For a project whose central claim is about
source separability (`CLAUDE.md` §1), that is a modelling-relevant fact and not
a bookkeeping detail.

The alternative — re-emitting HepG2 under `LC_2024` and dropping `YZ_2022` — is
*not* recommended: it would lose precision and would discard the open-access,
licence-clean version of the data (D-1).

STATUS: OPEN

---

### D-6 — Is a 1,080-formulation screen of ONE ionizable lipid wanted in the corpus at all?

This paper contributes **one** ionizable-lipid structure, DLin-MC3-DMA, which
the corpus already holds in 3,098 rows across eleven experiments. What is new is
**formulation space**: 1,080 + 72 compositions over six helper lipids, five new
cell lines.

There is precedent — `ZT_2026` is a 486-formulation DoE on one lipid — so this
is not a new category. But it is worth deciding deliberately, because the
corpus's source key is `Experiment_ID` and the analyses in `scripts/` are about
*chemical* transfer across sources. An experiment with one lipid and 1,152
formulations contributes:

- **nothing** to Result 3 (shared-lipid concordance needs lipids shared across
  labs; MC3 is already shared, and adding 4,680 more MC3 rows changes which
  pairs dominate the pooled estimate);
- **a great deal** to Result 6 (the within-experiment learnable regime) and to
  the within-lab cross-context question in `CLAUDE.md` §8, for which see below.

**The measurement that makes this interesting.** The same 1,080 formulations,
one lab, one library, one cargo, only the target cell changing, give
within-library Spearman against the HepG2 arm of:

| cell | n | Spearman |
|---|---|---|
| Neuro2a | 720 | 0.494 |
| B16F10 | 1,080 | 0.488 |
| HEK293T | 1,080 | 0.476 |
| PC3 | 1,080 | 0.447 |
| ARPE-19 | 720 | 0.330 |

That is a five-fold replication of the `molea.csv` result `CLAUDE.md` §8 relies
on (hepatocyte vs chondrocyte, ρ = 0.30, n = 1,152) — on a different lab, a
different cargo (pDNA, not mRNA) and a different axis of variation (formulation,
not head-group chemistry). It lands in the same place: **even within one lab and
one library, changing only the target cell costs most of the rank agreement.**
`molea`'s 0.30 is the low end of this range, not an outlier.

Whether that belongs in `data/` or stays a standalone companion point, as
`molea.csv` did, is the decision.

STATUS: OPEN
