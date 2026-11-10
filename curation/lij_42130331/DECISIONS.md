# DECISIONS -- lij_42130331 (Guo et al., ACS Nano 2026, PMID 42130331)

Questions this curation must NOT answer on its own. `build_dataset.py` refuses
to merge this folder while any entry below is `OPEN`.

WHY THIS FILE GATES THE MERGE. The pipeline is good at detecting that a source
is silent. It is not entitled to resolve it. In `AC_2025` the first repair was
confident, internally consistent, and would have assigned **419 measured values
to the wrong molecules** -- the names were wrong and the structures were right,
and nothing in the data said which. So the rule is: detect, queue, and stop.

## Format

```
### D-n  <one line>
STATUS:    OPEN | RESOLVED <date> by <who>
AMBIGUITY: what the source does not say.
CANDIDATES: each reading, and what it implies for the emitted rows.
EVIDENCE NEEDED: what would settle it.
DECISION:  (when resolved) what was chosen and why. Keep the entry.
```

## Entries

### D-1  What assay produced Fig. 1C and Fig. 1E -- i.e. what the 300 values ARE
STATUS:    OPEN
AMBIGUITY: The SI never describes the screen. Table S1 and Table S2 are
           captioned only as "Source data for the heatmap in Fig. 1C / Fig.
           1E"; the SI has no methods section, no figure legend for Fig. 1, and
           no units anywhere. The main text is paywalled: the article page
           resolves to `/article-abstract/`, is marked "Available to Purchase",
           and its "Figures & tables" control yields no caption text to an
           unentitled session (tested 2026-11-10; `document.body.innerText`
           stayed at 7,470 characters and contains no "Figure 1"). Europe PMC
           reports the paper as not open access, not in PMC, with no
           supplementary deposit. So **cell line, readout, units, dose,
           timepoint, and even in vitro versus in vivo are unknown.**
CANDIDATES:
  (a) An in-vitro luciferase reporter screen in a cell line, values in RLU.
      Supported, but NOT established, by three things in the SI: Fig. S3 is a
      cell-viability panel covering *both* libraries, which only makes sense if
      both were screened in cells; the magnitudes (1.7e2 - 9.8e4) are typical
      of a plate-reader luminescence well; and Fig. S4 reports an in-vitro
      luciferase comparison of the lead against DLin-MC3-DMA, SM-102 and
      ALC-0315 in **HeLa** cells.
  (b) The same, but in a different cell line from Fig. S4's HeLa. Fig. S4 is a
      benchmark comparison of the LEAD, which is downstream of the screen, so
      its cell line is not evidence about the screen's.
  (c) An in-vivo screen. The paper is spleen-tropic and in vivo throughout, but
      a 300-member in-vivo screen is implausible and Fig. S3's whole-library
      viability panel points the other way.
  Under (a) or (b) the rows belong in `data_vitro.csv`; under (c) in
  `data_vivo.csv`. `Model`, `Model_type`, `Model_target`, `Cargo` and
  `Experiment_method` are emitted BLANK on all 300 rows until this is settled;
  a blank is "not stated", not "clean".
EVIDENCE NEEDED: the main-text PDF, or the Fig. 1 caption, or the Methods. One
           page settles all five fields. An institutional session on
           pubs.acs.org would reach it; this curation had none.

### D-2  Are Library 1 and Library 2 one experiment or two?
STATUS:    OPEN
AMBIGUITY: Table S1 (50 values, range 17 - 1,980, median 380) and Table S2 (250
           values, range 150 - 98,000, median 1,100) are two panels of one
           figure, but they share **no compound at all** -- Library 2 is
           Library 1 esterified, so the two sets are disjoint by construction.
           Nothing anchors one table's scale to the other's.
CANDIDATES:
  (a) ONE experiment, `QG_2026`, with a single pooled z-score of log10 over all
      300 values. This is the corpus convention (one z per `Experiment_ID`) and
      it PRESERVES the ~50x between-library difference, which is the paper's
      central structure-activity claim -- esterifying the beta-hydroxy is what
      produces the lead. **This is what is currently emitted.**
  (b) TWO experiments, `QG_2026` and `QG_2026_2`, each standardised separately.
      This DELETES the between-library difference. `QL_2026` is the precedent:
      standardising its two libraries separately erased a 37x difference and
      inverted 36% of cross-library pairs while every check still passed. There
      the fix was a scalar offset tied through four shared benchmarks; here
      there are no shared compounds, so no such tie exists.
EVIDENCE NEEDED: whether Fig. 1C and Fig. 1E share a colour scale or an axis
           label (visible in the main-text figure), or a Methods sentence
           saying the two libraries were screened together. Note that this
           question is downstream of D-1: if the two panels turn out to be
           different assays, (b) follows immediately.
           The raw values are kept in `05_values.csv`, so either choice is a
           one-line change to `POOL_LIBRARIES` in `stage07_emit.py`.

### D-3  What formulation was the screen run at?
STATUS:    OPEN
AMBIGUITY: The SI names the excipients only for the LEAD compound's
           optimisation. Fig. S11 optimises "the ratios of the ionizable lipid,
           DOPE, and cholesterol while maintaining a constant ratio of
           DMG-PEG2000", and Fig. S12 runs a second round -- both downstream of
           the Fig. 1C/1E screen, so their ratios are not the screen's. Molar
           ratios, N:P, mRNA dose, mixing method and buffers are stated nowhere
           in the SI for any experiment.
CANDIDATES:
  (a) Emit the excipient IDENTITIES (DOPE, cholesterol, DMG-PEG2000) with blank
      ratios, on the grounds that Fig. S11 names them.
  (b) Leave all formulation columns blank, because Fig. S11 describes the
      optimisation of the lead and says nothing about what the 300-member
      screen used. **This is what is currently emitted.**
  F-004 is the precedent for refusing to reconstruct: `WT_2026`'s N:P had three
  mutually inconsistent candidates and was blanked on the user's decision.
EVIDENCE NEEDED: the main text's Methods.

### D-4  Does a library whose connectivity has no independent check belong in the corpus?
STATUS:    OPEN
AMBIGUITY: Every one of the 300 structures rests on two drawings, Fig. S1A
           (the scheme) and Fig. S1B / Fig. S2 (the blocks). The SI contains no
           HRMS, no NMR table, no SMILES, no InChI and no deposited structure
           file -- a grep of the full extracted text for
           `smiles|inchi|HRMS|ESI-MS|calcd` returns **0**. The only
           characterisation of any member is Fig. S9, a mass spectrum and a 1H
           NMR of the lead published as IMAGES with no printed m/z and no shift
           list. The two bonds most at risk are exactly the two that a
           molecular formula cannot see:
             * where the thiolactone ring opens (drawn: at the carbonyl, giving
               an amide to the amine and a free thiol);
             * which epoxide carbon takes the sulfur (drawn: the terminal CH2,
               leaving the hydroxyl on the carbon that carries R).
           Both alternatives are C- and H-identical. This is the Lin
           regiochemistry situation with the SI synthesis route removed.
CANDIDATES:
  (a) Merge all 300 with `confidence = medium` and the regiochemistry stated in
      `curation_note`. **This is what is currently emitted.**
  (b) Merge with `confidence = low`, on the grounds that "no external check
      exists at all" is weaker than the `medium` rows elsewhere in the corpus,
      which mostly have HRMS or a named reagent behind them.
  (c) Hold the whole library out until the main text's synthetic procedures can
      be read.
EVIDENCE NEEDED: the main text's synthesis section, or the numbers behind
           Fig. S9. The m/z of the lead alone would confirm the arm COUNT (and
           therefore the secondary-amine rule that gives 96 and 400 two arms);
           it would still say nothing about either regiochemistry.

### D-5  The source is not licensed for redistribution
STATUS:    OPEN
AMBIGUITY: `CLAUDE.md` section 8 requires every constituent source to be
           CC0/CC-BY-compatible before corpus release. This article is
           "(c) 2026 American Chemical Society", marked "Available to Purchase"
           on its landing page, and Europe PMC reports `isOpenAccess = N`. No
           Creative Commons licence is offered on any part of it, including the
           SI these 300 values come from.
CANDIDATES:
  (a) Merge for internal analysis and exclude from any released corpus.
  (b) Merge and release only the derived values (a z-score of a number the
      authors printed), treating the extracted numbers as facts rather than as
      an expression. That is a judgement for the user, not for this curation.
  (c) Do not merge.
EVIDENCE NEEDED: the user's call, consistent with whatever is decided for the
           other subscription-only sources already in the corpus.
