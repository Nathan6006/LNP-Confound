# DECISIONS -- liu_42585110 (KJ_2026)

Questions this curation must NOT answer on its own. `build_dataset.py` refuses
to merge this folder while any entry below is `OPEN`, and two of them are.

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

### D-1  The article is All Rights Reserved, which fails CLAUDE.md section 8
STATUS:    OPEN
AMBIGUITY: None about the licence itself -- it is not ambiguous. What is
           undecided is whether this project accepts a non-CC source.
           CLAUDE.md section 8 requires every constituent source to be
           CC0/CC-BY-compatible before corpus release, and this one is not.
           Three independent records agree and none of them is the publisher's
           landing page: Europe PMC's core record for EXT_ID:42585110 returns
           `isOpenAccess: "N"`, `license: null` and availability
           "Subscription required"; Crossref lists exactly two licences for
           the DOI, Wiley's own `#vor` terms and `tdm_license_1.1`, neither
           of them Creative Commons; and page 1 of the PDF reads verbatim
           "All rights reserved, including rights for text and data mining and
           training of artificial intelligence technologies or similar
           technologies. (c) 2026 Wiley-VCH GmbH".
           **The CC line in this PDF is a trap.** Every page carries Wiley's
           download stamp ending "...OA articles are governed by the
           applicable Creative Commons License". That is conditional
           boilerplate printed on every Wiley PDF regardless of access status
           and it does not assert that THIS article is OA. The briefing for
           this curation read it as a CC licence; the three records above are
           what settled it.
CANDIDATES:
  (a) Merge and flag, as the corpus already does for 14 CC BY-NC-ND rows.
  (b) Merge but exclude this source from any released corpus.
  (c) Do not merge.
EVIDENCE NEEDED: A ruling from the user on section 8's scope -- specifically
           whether it governs the released corpus only, or also the working
           one. Nothing further can be learned about the licence itself.

### D-2  What IS the quantity Figure 1a's -1.0 .. +1.0 colour scale encodes?
STATUS:    OPEN
AMBIGUITY: The caption says "Luciferase expression in DC2.4 cells after 6 h
           incubation", and that is the whole of what the paper says. The
           colorbar has NO axis title, NO tick marks and no units; the article
           and the SI never state what quantity was scaled or how. This
           matters because the corpus's `Experiment_value` is defined as a
           per-experiment z-score of **log10** signal.
WHAT IS MEASURED, and is not in doubt:
  * The scale is a MIN-MAX RESCALING onto [-1, +1] of something monotone in
    expression. The brightest cell's RGB equals the colorbar's top interior
    colour and the darkest equals its bottom interior colour, so the colour
    limits were set to the data's own extremes; nothing is clipped, so no cell
    is a bound.
  * The read itself is sound on three independent tests: the paper's three
    named top performers are our top three in the paper's own order; its
    stated ATX-RA > ATX-c-RA series effect reproduces at p = 0.0045 with 15/20
    matched-pair wins; and the diverging map's white point lands at +0.0094,
    i.e. 0.47% of the half-range from where the +-1.0 labels imply. Max LUT
    colour distance 3.46 against the 5.0 gate.
  * The three unlabelled rows ARE replicates: inter-row Spearman 0.983-0.996,
    within-column sd 0.037 against between-column sd 0.410.
CANDIDATES:
  (a) The plotted quantity is LOG-scaled. Then a per-experiment z of the panel
      value is EXACTLY the corpus convention, because a z of an affine
      function of log10 signal is the z of the signal -- the identity
      CLAUDE.md section 3 already relies on for
      `Experiment_value_within_model`. The 41 emitted in-vitro rows are
      correct as they stand.
  (b) The plotted quantity is LINEAR in luminescence. Then the emitted
      in-vitro `Experiment_value` is a z of a linear quantity, which is off
      the corpus convention. Ranks are unaffected either way; values are not.
  (c) Drop the in-vitro arm and merge only the 82 in-vivo rows, which are in
      stated physical units and carry no such question.
EVIDENCE TRIED, AND IT DID NOT SETTLE IT: the shape of the distribution. If
           panel a were linear in luminescence it should be strongly
           right-skewed, as this paper's own in-vivo radiances are (skew +0.71
           to +1.73, median at 0.19-0.29 of range). In log10 those same data
           are left-skewed (skew -0.59 to -1.16, median at 0.61-0.66). Panel a
           is nearly symmetric: skew +0.10, median at 0.42 of range. That is
           not the linear shape, which weakly favours (a), but it matches
           neither reference cleanly and the in-vitro distribution need not
           resemble the in-vivo one. **Not determinable from what the paper
           publishes.**
EVIDENCE NEEDED: a statement from the authors, or the underlying DC2.4
           luminescence values.

### D-3  `Dose_ug_nucleicacid` is blank on all 123 rows
STATUS:    RESOLVED 2026-11-10 by this curation, following F-004's posture
AMBIGUITY: The column is in micrograms. In vitro the SI gives a CONCENTRATION,
           0.5 ug/mL for 6 h in a 96-well plate, and never states the well
           volume. In vivo the caption gives 0.1 mg/kg IV and no animal mass
           is stated anywhere.
CANDIDATES: in vitro, 0.05 ug at 100 uL or 0.1 ug at 200 uL -- both plausible,
           neither stated. In vivo, 2.0 ug at a nominal 20 g mouse, which is
           what several corpus experiments appear to carry, but the conversion
           is ours and not the paper's.
DECISION:  BLANK, with the printed figure in `curation_note` on every row.
           F-004's `N:P` had three mutually inconsistent candidates and was
           blanked on the user's instruction; a dose that requires an
           unstated volume or an unstated body mass is the same shape. The
           information is not lost -- it is on the row, in the units the paper
           prints -- and a future reader can convert knowingly.

### D-4  `N:P` and the lipid-to-mRNA mass ratio are blank
STATUS:    RESOLVED 2026-11-10 by this curation
AMBIGUITY: None. They are not stated.
DECISION:  BLANK. Measured rather than assumed: a regex sweep for `N/P`,
           `N:P`, "nitrogen...phosph", "weight ratio", "mass ratio", `w/w`,
           `wt/wt` and "encapsulation" over the SI's full text and the main
           PDF's text layer returns ZERO hits in the SI and one hit in the
           main text, which is an Introduction sentence about what ionizable
           lipids do. There is nothing to recover.

### D-5  The benchmark is emitted under the corpus spelling, not the paper's
STATUS:    RESOLVED 2026-11-10 by this curation
AMBIGUITY: The paper writes `Dlin-MC3` throughout. The corpus holds the same
           molecule under `DLin-MC3-DMA` (3,081 vitro rows) and `Dlin-MC3-DMA`
           (17 vivo rows). The paper's own spelling matches **0** corpus rows
           by name.
DECISION:  Emit `DLin-MC3-DMA` with the structure from
           `provenance/reagent_registry.py`, which states it from the
           chemical name rather than inheriting whatever the corpus contains.
           The SI's own expansion, "Dilinoleylmethyl-4-dimethylaminobutyrate",
           is that compound. Recorded here because the alternative -- keeping
           the paper's spelling -- would add a THIRD spelling of a molecule
           that already has two, and C17 pins one reagent name to one
           structure. The paper's spelling is preserved in
           `02_corpus_overlap.json`.

### D-6  The PEG-lipid's SMILES comes from the corpus, not from a registry
STATUS:    RESOLVED 2026-11-10 by this curation
AMBIGUITY: `provenance/reagent_registry.py` deliberately holds no PEG-lipid,
           because a PEG block length is a distribution and no single SMILES
           is right. The SI names "methoxy-PEG2000-carbamate-dimyristoyl
           glycerol (PEG2k-C-DMG)".
DECISION:  Take it from the corpus through `corpus_components.smi()`, the
           guarded loader that RAISES on an ambiguous name, under the name
           `DMG-C-PEG2000`. The match is not a coincidence of spelling: that
           corpus entry is C119H235NO49 with two myristoyl esters and one
           carbamate, which IS methoxy-PEG-carbamate-dimyristoyl glycerol.
           The corpus's far more common `DMG-PEG2000` is the ETHER-linked
           C122H242O50 and would have been the wrong molecule -- the
           `C14-PEG2000` trap named in EXTRACTION_BRIEF section 5, where a
           name unambiguous in the corpus is simply wrong for the paper in
           hand. The PEG length remains a convention inherited from the
           corpus and is not a claim about this paper's material.
