# DECISIONS -- whiteheadka_30637934

Hajj, Ball, Deluty, Singh, Strelkova, Knapp, Whitehead.
*Branched-Tail Lipid Nanoparticles Potently Deliver mRNA In Vivo due to
Enhanced Ionization at Endosomal pH.* Small 15(6), 2019.
PMID 30637934, DOI 10.1002/smll.201805097. Experiment_ID **KW_2019**.

Questions this curation must NOT answer on its own. `build_dataset.py` refuses
to merge this folder while any entry below is `OPEN`.

WHY THIS FILE GATES THE MERGE. The pipeline is good at detecting that a source
is silent. It is not entitled to resolve it. In `AC_2025` the first repair was
confident, internally consistent, and would have assigned **419 measured values
to the wrong molecules** -- the names were wrong and the structures were right,
and nothing in the data said which. So the rule is: detect, queue, and stop.

Evidence for every entry is in `08_contradictions.csv`, `07_qc.json` and
`07_seeded_checks.json`.

## Entries

### D-1  Which excipient takes 16 mol% and which takes 46.5 mol%
STATUS: OPEN

AMBIGUITY: SI Methods p.2 reads "mixing lipidoid, cholesterol (Sigma Aldrich),
1,2-dioleoyl-sn-glycero-3-phosphoethanolamine (DOPE, Avanti Polar Lipids),
C14-PEG2000 (Avanti Polar Lipids), at a molar ratio of 35 : 16 : 46.5 : 2.5".
Read positionally that is cholesterol 16 / DOPE 46.5. The four values sum to
exactly 100 either way, so the sum cannot discriminate and nothing else in this
SI can. It is the mapping that is in doubt, not the numbers.

CANDIDATES:
  (a) **cholesterol 16, DOPE 46.5** -- the sentence's own component order.
      Implies 11 rows whose excipient ratios match NO other experiment in the
      corpus.
  (b) **DOPE 16, cholesterol 46.5** -- the corpus convention and the field's.
      Implies the paper's number order is transposed relative to its component
      list.

EVIDENCE: 20 corpus experiments use this exact `35 / x / y / 2.5` ratio and
every one of them assigns helper lipid 16 and cholesterol 46.5. **Zero rows**
in the 29,195 vitro + 1,803 vivo corpus have `CHL_molratio` 16 or
`HL_molratio` 46.5. `AC_2025` varies the pair across four formulations
(28/34.5, 8/54.5, 16/46.5, 50/12.5), always summing to 62.5, with cholesterol
the larger in three of four. `BL_2023_3CR` uses the identical printed reagent
name `C14-PEG2000` with DOPE 16 / cholesterol 46.5. Separately, the same lab's
2024 paper -- staged, unmerged, at `curation/whiteheadka_38888972` -- records
its own Methods as "lipidoid : cholesterol : DOPE : PEG2000 = 35 : 46.5 : 16 :
2.5", i.e. candidate (b) with the component order this paper uses.

EMITTED PENDING A RULING: candidate **(b)**, with the discrepancy written into
`curation_note` on all 11 rows. Reason: an excipient field that differs from
every other experiment in the corpus for a prose-ordering reason would
manufacture precisely the kind of source-separability signal this project
exists to measure. That is a reason to prefer (b) as the default, not a reason
to close the decision.

EVIDENCE NEEDED: the main article's Methods (behind a subscription), or any
other Whitehead-lab paper stating this formulation with an unambiguous pairing.

### D-2  Whether the in-vivo panels are digitised and emitted at all
STATUS: OPEN

AMBIGUITY: three SI panels carry in-vivo values and none is a clean screen.
`Figure S3D` is a time course for 3 lipidoids on a log axis with only 3 labelled
decade ticks. `Figure S4` is 2 lipidoids x 3 doses where each bar is a SIX-ORGAN
STACKED TOTAL (measured, see `08_contradictions.csv` X-4), so the drawn top is a
sum and the per-organ values are differences of near-coincident cumulative
rules. `Figure S5` is 2 lipidoids x 3 timepoints of serum **EPO protein by
ELISA**, a different cargo and a different readout class, on two panels with
different y limits (0-2000 vs 0-200).

CANDIDATES:
  (a) emit nothing in vivo -- the current state.
  (b) digitise `Figure S4`'s six stacked totals at 0.75 mg/kg and emit 2 rows
      as a `KW_2019_2` in-vivo experiment.
  (c) digitise all three.

EVIDENCE: `CLAUDE.md` section 8 already records that the in-vivo pool's only
cross-source chemistry is MC3 and SM-102, and flags benchmark-only pairs as a
live problem for Result 3. A 2-lipidoid in-vivo environment adds an
`Experiment_ID` and no new shared chemistry: both 306Oi10 and 306O10 are already
in the corpus in vitro (NC_2024) and 306O10's skeleton is in KW_2014 in vivo
under the name `A306_C10`. Nothing is lost by waiting and something is risked by
adding a 2-member environment.

EVIDENCE NEEDED: a ruling on whether 2-member in-vivo environments are wanted
at all, which is the same question as the open `SB_2024`/`SY_2025` item.

### D-3  Whether Figure S6B becomes a separate Experiment_ID
STATUS: OPEN

AMBIGUITY: Figure S6B's 100 ng series shares four lipidoids and the identical
nominal dose with Figure S6A and disagrees with it by **0.50x to 1.28x**, with
the rank of those four **inverted** (Spearman rho = **-0.20**). Both are HeLa /
unmodified firefly-luciferase mRNA / RLU / 24 h / n=4. So S6B is a second plate,
and RLU is gain-dependent -- the SI itself warns of exactly this for its TNS
readout on p.3 ("the absolute value of fluorescence values will depend on plate
reader settings, including the gain").

CANDIDATES:
  (a) emit S6A only -- the current state. 11 rows, 1 environment.
  (b) emit S6B as `KW_2019_2`, 12 rows over 4 lipidoids x 3 doses, with its own
      per-experiment z. Honest about the two plates, but creates a 4-lipid
      environment whose dose varies within it.
  (c) pool both under `KW_2019` with one z. **Rejected outright**: that deletes
      a between-assay offset the data says is real, which is error class 22 and
      the F-016 / `WT_2026` defect.

EVIDENCE: measured ratios and ranks in `02_fig_s6_geometry.json`
(`cross_panel_S6A_vs_S6B_100ng`). Note which way the disagreement cuts: S6B
makes 306Oi10 the best at 100 ng, agreeing with the paper's headline, while S6A
makes 306O10 the best by 1.48x. The abstract calls both "the two most potent
materials" without ordering them in vitro and locates the tenfold 306Oi10
advantage **in vivo**, so the paper does not adjudicate between its own panels.

EVIDENCE NEEDED: a ruling on (a) vs (b). Nothing further is recoverable from
the source.

### D-4  The isodecyl branch in 306Oi10
STATUS: OPEN

AMBIGUITY: the SI names only "alkyl-acrylate tails of varying lengths
(Sartomer, TCI Chemicals)". Isodecyl acrylate is sold as a BRANCHED-ISOMER
MIXTURE, so "Oi10" does not denote one molecule. `CLAUDE.md` section 6 already
lists "the exact branch of iso-decyl acrylate" among the things that are
unverifiable from the paper alone, and `ZT_2026` carries 6 rows at
`confidence = medium` for this same reagent.

CANDIDATES: (a) store 8-methylnonyl, the corpus convention; (b) store a
different branch; (c) refuse to store 306Oi10 at all.

EMITTED PENDING A RULING: (a), flagged `confidence = medium`,
`structure_assumed = 1`, with the reason on the row. The corpus's `NC_2024`
`A306_Oi10` is **full-InChIKey identical** to what this curation built, so the
convention is at least consistent across two independent curations of the same
reagent name -- which is agreement on a convention, not evidence about the
molecule.

EVIDENCE NEEDED: the Sartomer/TCI catalogue number, or an SI NMR. Neither is in
the source.

### D-5  Whether all 11 rows should carry structure_assumed = 1 for crude product
STATUS: OPEN

AMBIGUITY: SI Methods p.2 states "In vitro experiments were conducted with
**crude lipidoids**. For in vivo experiments, the lipidoids were purified ... to
isolate the fully substituted lipidoid product." Figure S6A is in vitro, so
every emitted value was measured on a crude 1:4 reaction product -- a mixture of
substitution states -- while the stored structure is the nominal
tetra-substituted species.

CANDIDATES:
  (a) `structure_assumed = 1` on all 11 with a note -- the current state.
  (b) leave the flag blank, matching `KW_2014` and `NC_2024`, whose A306 rows
      carry no flag for the same chemistry. But a blank means "not assessed",
      not "clean", so (b) would be withholding an assessment we actually made.
  (c) `positional_ambiguity = 1` as well. **Not emitted**: the 306 core's four
      N-H sites are two symmetric pairs, so the tetra-product has no
      regiochemical choice. The uncertainty is substitution DEGREE, which
      `structure_assumed` carries. Stated on the row so the distinction is not
      lost.

EVIDENCE NEEDED: a ruling on whether "crude one-pot product" is a corpus-wide
flag. If it is, it applies to far more than these 11 rows -- most lipidoid
libraries in the corpus are crude -- and should be set by the builder, not by
one curation.
