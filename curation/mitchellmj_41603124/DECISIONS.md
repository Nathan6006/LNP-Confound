# DECISIONS — `mitchellmj_41603124`

Hamilton et al., *Adv Mater* 38(13):e14370 (2026). PMID **41603124**, PMCID
PMC12957869, DOI 10.1002/adma.202514370. **Michael J. Mitchell lab (UPenn)**,
verified from the Europe PMC `core` record's `authorString` — Mitchell MJ is
last author, Hamilton AG first, so the proposed `Experiment_ID` would be
**`AH_2026`**, not an `MM_*`.

Licence **CC BY** (`license: "cc by"`, `isOpenAccess: "Y"`). No repository
exists, so no repository licence.

`build_dataset.py` refuses to merge a folder with an `OPEN` entry. That is the
intended state of this folder.

---

### D-1
STATUS: OPEN

**QUESTION:** Is the b-mRNA screen — 3,840 attributable values, the largest
single block of in-vivo LNP measurements this project has seen in one paper —
a delivery endpoint?

**WHAT THE PAPER MEASURES.** Normalized intracellular b-mRNA accumulation:
total RNA is extracted from FACS-sorted cells, reverse-transcribed, and the
barcode region sequenced. The quantity is **how much of the delivered mRNA is
present inside the cell**, not how much protein it made.

**THE AUTHORS DRAW THE SAME DISTINCTION, EXPLICITLY.** Results, verbatim:
*"As intracellular exogenous mRNA is vulnerable to enzymatic degradation, we
opted to analyze b-mRNA accumulation at a relatively short time following
treatment to avoid potential screening bias introduced by extensive
degradation. However, peak mRNA translation occurs relatively late after
treatment."* They sampled deliberately **before** translation, and then ran
Figure 3 — single-plex NanoLuc and EGFP — *"to confirm transfection of lead
mRNA LNP candidates identified by our high-throughput b-mRNA screen."* Figure 3
exists because accumulation is not transfection.

**WHY POOLING IS NOT THE OBJECTION.** `data/` has `Experiment_batching` with
vocabulary `{barcoded, individual}`, and `LR_2023` contributes **407 of 1,803**
in-vivo rows at `barcoded`. But `LR_2023` carries `Cargo_type =
peptide_barcode` and `Experiment_method = protein_adundance_normalized` — a
**translated** readout, detected as peptide. The corpus precedent is for a
pooled screen whose readout still requires expression. This one does not.

**THE VOCABULARY.** Measured from `data/data_vivo.csv` on 2026-11-10,
`Experiment_method` takes exactly four values: `luminescence_normalized` (836),
`protein_adundance_normalized` (766), `editing_efficiency_normalized` (141),
`luminescence_relative_to_Spikevax` (60). Sequenced cargo abundance is none of
them.

**CANDIDATES:**
1. Exclude. The readout is cargo biodistribution into sorted cells, which §
   the project's own endpoint list puts with cellular uptake and particle
   biodistribution. This is what `05_values.csv` records.
2. Admit as a fifth `Experiment_method`, e.g. `cargo_abundance_normalized`,
   and let the modelling decide whether it pools with transfection. The paper
   itself provides the evidence needed to judge that: Figure 2c–e gives the
   full correlation structure, and Figures 3b/3f/3j give a 7-LNP paired
   comparison of accumulation against both luminescence and protein.
3. Admit only the populations where the paper's own validation shows
   accumulation tracking transfection, and exclude the rest.

**EVIDENCE NEEDED:** a project decision on whether the corpus endpoint is
"delivered and expressed" or "delivered". Note the direction of the risk:
candidate 2 would add 3,840 rows from one source to a 1,803-row in-vivo pool —
**more than doubling it, and tripling it as a single `Experiment_ID`** — which
bears directly on Result 4's leave-one-source-out and on §8's open question
about benchmark-dominated pairs.

---

### D-2
STATUS: OPEN

**QUESTION:** The validation panels (Figure 3b, 3f, 3j) **are** delivery
readouts and **are** attributable. Should a partial experiment be emitted from
them when 5 of the 7 ionizable lipids have no determinable structure?

**WHAT IS DETERMINABLE.** The paper draws 3 of its 120 library members in full
(Supplementary Figure 2: `17D9.2`, `18D4`, `20D8i`), with a mass spectrum and a
¹H NMR for each. Of the seven formulations in Figure 3:

| LNP | lipid | structure | why |
|---|---|---|---|
| 97 | `18D4` | **determined** | drawn; MS confirms a+b = 8.017 carbons, Δ 0.612 Da against a nearest alternative 28.053 Da away — a **45.8×** margin |
| 125 | `SM-102` | **determined** | `provenance/reagent_registry.py`, which agrees with the live corpus by full InChIKey |
| 96 | `17D9.2` | **no** | drawn, and the three printed ions pin a+b = 18.993 carbons on TREN with three arms — but mass cannot see connectivity, and 11+8 and 10+9 are isomeric, as is every branch position, on a 137 dpi raster |
| 112 | `20D8i` | **no** | drawn; M = 1878.58 ± 0.41 Da, but amine 20 is a branched polyamidoamine whose core and arm count were not resolved at 137 dpi, so the arithmetic has two unknowns |
| 94 | `17D8i` | **no** | composable in principle (amine 17 from `17D9.2`, maleate `8i` from `20D8i`) but inherits `20D8i`'s unresolved maleate |
| 9 | `2D8` | **no** | amine 2 and maleate 8 are drawn nowhere in this paper |
| 37 | `7D4` | **no** | amine 7 is drawn nowhere in this paper |

Amines 1–16 and 19, and maleates `6.2`, `8` and `18`, are defined only in
**reference 34** — Xu, Han, Xu, Ricciardi et al., *"Plug-and-play assembly of
biodegradable ionizable lipids for potent mRNA delivery and gene editing in
vivo"*, bioRxiv Mar 2025 — a **different paper**, not in `source/`.

**WHY A PARTIAL EMIT IS WORSE THAN NONE HERE.** `data/` has **0 blank
`IL_SMILES` in 30,998 rows**, so the 5 undeterminable lipids cannot be emitted.
That leaves 2 of 7 — and the two that drop out at the bottom, `2D8` and `7D4`,
carry the **two lowest liver values in the panel** (7.22e7 and 6.22e7 p/s
against SM-102's 2.83e11). `Experiment_value` is a per-experiment z-score of
log10, so truncating the bottom of a 7-point range changes every z that
survives. A 2-row experiment built from the top of a truncated distribution is
not a cheaper version of the real one; it is a different quantity.

**AND THE 2 THAT SURVIVE ARE CONFOUNDED.** The six `nDm` LNPs share one recipe
exactly (IL 35 / DOPE 16 / Chol 46.5 / C14-PEG2000 2.5). `SM-102` uses the
clinical recipe (50 / DSPC 10 / Chol 38.5 / DMG-PEG2000 1.5). So the one
surviving pair differs in **both** ionizable lipid and formulation.

**CANDIDATES:**
1. Emit nothing. `06_rows.csv` is header-only. This is the current state.
2. Emit the 2 determinable lipids × 2 organs as a 4-row in-vivo experiment,
   with `confidence = low` and a note that the range is truncated.
3. Obtain reference 34, determine amines 1–20 and all six maleates from it,
   and revisit — which would unlock not 7 lipids but potentially all 120.

**EVIDENCE NEEDED:** reference 34. Note this is the `amino_lipid_fg_2026`
situation in reverse: there the blocker was which members exist; here the
members are fully enumerated (Supp Table 3, 133 rows, Index contiguous, 133/133
mol% summing to exactly 100) and the blocker is purely what the building blocks
*are*.

---

### D-3
STATUS: OPEN

**QUESTION:** If D-1 ever resolves toward admitting the screen, is a value that
is its LNP's **share of its own pool** poolable across sources?

**THE FORMULA, STATED IN PROSE AND NEVER INSTANTIATED.** Methods 4.6, verbatim:
*"Normalized accumulation was calculated as the ratio of output read fraction
to input read fraction. Specifically, read fractions were first computed by
dividing the read count for each barcode by the total number of barcode counts
within each sample. Each sample read fraction (output) was then divided by the
corresponding mean read fraction from the uninjected pool samples (input) to
yield normalized accumulation."* That is

    nd_i(s) = [ c_i(s) / Σ_k c_k(s) ] / mean_p[ c_i(p) / Σ_k c_k(p) ]

**COMPOSITIONAL: YES, in the numerator.** `c_i(s)/Σ_k c_k(s)` is LNP *i*'s
share of its own pool, so every value depends on its pool-mates.

**BUT ONE STEP LESS SO THAN THE SIBLING'S, AND THAT MATTERS.** The sibling's
formula carried an outer renormalisation, `/ Σ_k[same] × 100%`, forcing the
values onto a unit simplex. This paper has **no such step**, and its
denominator is a **fixed reference** — the mean input fraction from uninjected
pool samples, the same number for every output sample. The pool-mate dependence
is removed from the input side. That is the materially better case.

**THE CONSTANT SUM SURVIVES ANYWAY, AND THE PAPER'S OWN DESIGN IS WHY.**
`Σ_i nd_i(s) = Σ_i f_out,i(s)/f_in,i`. The barcode pool was narrowed from 200
to 134 precisely to make the `f_in,i` uniform — Figure 1g, verbatim: *"b-mRNAs
demonstrating 20% enrichment or depletion (outside the shaded area) were
excluded from later use."* With every retained `f_in,i` within ±20% of 1/N,
`Σ_i nd_i` is pinned to **[N/1.2, N/0.8] = [110.8, 166.3]** for N = 133. The
values still lie on a simplex; it is scaled by N instead of by 100.

**RECONSTRUCTION RESIDUAL: NOT COMPUTABLE, and that is the finding.** This
paper publishes no read counts, no input fractions, no per-LNP normalized
accumulation values and no worked example. Data availability reads *"available
from the corresponding author upon reasonable request"*, and a search of the
article XML and all 22,296 characters of SI text returns **zero** accessions
(GEO, SRA, PRJNA, E-MTAB, Zenodo, Dryad, GitHub, figshare). The formula is
fully stated and never instantiated, so there is nothing numeric to reconstruct
against. The sibling reached residual 0.0 because its paper tabulated a worked
example; this one does not.

**WITHIN THIS PAPER THE PROBLEM DOES NOT BITE.** One pool of 133 LNPs serves
every cell population and both strains, so all 3,840 values sit on the *same*
simplex. It bites across sources — which is exactly the comparison this project
measures.

**CANDIDATES:** (1) never pool a compositional readout across sources; (2) pool
only after a per-source rank transform, which is what Stage 5 already does and
which is invariant to the scaling; (3) treat pool size as a covariate.

**EVIDENCE NEEDED:** a project decision, and it is not specific to this paper —
it applies to `LR_2023`'s 407 merged rows too, which were merged without it
being asked.

---

### D-4
STATUS: OPEN

**QUESTION:** Figure 2 plots neither 122 nor 133 LNPs. Which number describes
the screen?

Measured from the figure's own geometry: borders detected as neutral-dark
columns at x = 48 and x = 530, 481 px of cells, dominant FFT component **124
cycles identical in all five strips**, and **123–125** over a ±2 px sweep of
both borders. 133 and 134 are 8–9 bins outside; 122 is one column below.
Supplementary Table 3 lists 133 formulations; the abstract claims 122 were
screened simultaneously.

What is certain is that **Figure 2 does not plot 133 LNPs**, so the
Conclusions' "library of 133 distinct mRNA LNPs" is the number *formulated*.
Roughly nine to eleven formulations are in Supp Table 3 and not in Figure 2,
and the paper never says which or why — Supplementary Figure 10 reports
physicochemical characterisation for the library with `n ≥ 4` but prints no
exclusion criterion.

**AFFECTS:** the denominator N in D-3's constant-sum bound, and any future row
count. **EVIDENCE NEEDED:** the authors' QC criterion, or the underlying counts.
