# DECISIONS -- whiteheadka_38888972

Strelkova Petersen, Weiss, Hajj, Yerneni, Chaudhary, Newby, Arral, Whitehead.
*Branched-Tail Lipid Nanoparticles for Intravenous mRNA Delivery to Lung
Immune, Endothelial, and Alveolar Cells in Mice.* Adv Healthc Mater 2024.
PMID 38888972, DOI 10.1002/adhm.202400225.

Questions this curation must NOT answer on its own. `build_dataset.py` refuses
to merge this folder while any entry below is `OPEN`.

WHY THIS FILE GATES THE MERGE. The pipeline is good at detecting that a source
is silent. It is not entitled to resolve it. In `AC_2025` the first repair was
confident, internally consistent, and would have assigned **419 measured values
to the wrong molecules** -- the names were wrong and the structures were right,
and nothing in the data said which. So the rule is: detect, queue, and stop.

Evidence for every entry below is in `08_contradictions.csv` and `07_qc.json`.

## Format

(each real entry begins with a level-3 Markdown heading; the `###` is omitted
in this example on purpose -- `build_dataset.py`'s `open_decisions()` parses
line by line and does NOT skip fenced code blocks, so a literal `### D-n` here
would register as a tenth OPEN decision. Reported, not worked around.)

```
D-n  <one line>
STATUS:    OPEN | RESOLVED <date> by <who>
AMBIGUITY: what the source does not say.
CANDIDATES: each reading, and what it implies for the emitted rows.
EVIDENCE NEEDED: what would settle it.
DECISION:  (when resolved) what was chosen and why. Keep the entry.
```

## Entries

### D-1  Three amine heads: the Table S1 NAME and the Figure 1B DRAWING disagree
STATUS:    OPEN
AMBIGUITY: Heads **113**, **200** and **303** are drawn in Figure 1B as
           different molecules from the ones Table S1 names (contradictions
           C-1, C-2, C-3). In all three the drawing is one methylene short:
           113 loses its N-methyl, 200 and 303 each lose two carbons from an
           arm, which in those two makes the drawn structure an **aminal**
           (N-CH2-N) -- not a stable isolable amine and not what Enamine or
           Sigma sells under those names.
CANDIDATES:
  (a) **The NAME is right** (what this folder emitted). Supported by: these are
      catalogue reagents with vendor CAS numbers; two of the three drawings are
      chemically implausible; and for head 113 a SECOND, INDEPENDENT curation
      already in the corpus -- `KL_2010`'s `IL_head_SMILES` for head `113` --
      is `CN(CCN)CCN`, the name's structure. 60 of 580 structures, 5 of the 49
      emitted rows (200O10, 200O12, 200Oi10, 303X2, 303Oi10).
  (b) **The DRAWING is right.** Then head 113 has 5 N-H rather than 4, so every
      113 lipidoid gains a fifth tail, and heads 200 and 303 lose two carbons
      each. This changes formula and mass, not just connectivity.
  (c) Emit nothing for those three heads and drop the 5 affected rows.
EVIDENCE NEEDED: the vendor catalogue numbers (TCI for 113, Enamine for 200,
           Sigma for 303), or any characterisation of a 113/200/303 lipidoid in
           this group's earlier papers. `LNP_ML/data_files/Whitehead_siRNA/`
           uses the same head numbering and may settle 113 directly.
NOTE:      This is the `AC_2025` shape exactly -- a name and a structure
           disagreeing, with no internal evidence for which half is wrong -- so
           it is queued rather than decided, even though (a) looks strong.

### D-2  What to store for the five heads that 4 equivalents cannot fully substitute
STATUS:    OPEN
AMBIGUITY: Methods says "Stoichiometric ratio of 1:4 was used" and Figure 1A
           draws `306 + 4 Oi10`, on a head with exactly 4 N-H. Five heads have
           more: **174** and **200** have 5 N-H, **110** and **122** have 6,
           **111** has 7 (counted by RDKit, not by eye). Four equivalents
           cannot fully substitute them, so those 100 products are necessarily
           under-substituted mixtures, and the paper names each as one
           compound. (Contradiction C-11.)
CANDIDATES:
  (a) **Full substitution** (what this folder emitted), with
      `confidence = low` and `structure_assumed = 1` on all 100. Matches the
      corpus's own precedent for this chemistry -- `AA_2008`, `KL_2010` and
      `KW_2014` all store the fully substituted lipidoid -- so the column stays
      comparable across experiments.
  (b) **Exactly four tails**, placed... where? The stoichiometry fixes the
      count but not the positions, so this needs a positional convention and
      `positional_ambiguity = 1` on all 100.
  (c) **Exclude those 100** from the corpus. Affects 8 of the 49 emitted rows:
      110X7, 111O10, 122X4, 122Oi10, 174O11, 200O10, 200O12, 200Oi10.
EVIDENCE NEEDED: any mass spectrum, elemental analysis or stated MW for a
           lipidoid built on 110/111/122/174/200 from this group. None is in
           the paper or its SI.

### D-3  Is TNS fluorescence an admissible corpus endpoint?
STATUS:    OPEN
AMBIGUITY: 531 of the 580 members have a TNS fluorescence value (Figure 2A,
           extracted to `06_tns.csv`) and no delivery measurement at all. TNS
           reports LNP **surface ionization**, which is a physicochemical
           property of the formulation in the same class as diameter, zeta
           potential and uptake -- and the corpus's in-vivo definition
           explicitly excludes those (CLAUDE.md changelog 2026-07-23). It is
           also measured on a polyU-mRNA formulation, i.e. on particles that
           carry no reporter and were never dosed.
CANDIDATES:
  (a) **Not an endpoint** (what this folder assumed): emit 49 rows; keep
      `06_tns.csv` in the folder as a non-corpus artifact. The paper then
      contributes 49 rows and 580 structures.
  (b) **An endpoint under a new `Experiment_method`** such as
      `TNS_fluorescence_pH5`: 580 rows from one paper, more than any existing
      vitro experiment except `WT_2026`. This would be the first non-delivery
      readout in the corpus and would change what "pooling delivery
      measurements" means in the manuscript's own framing.
EVIDENCE NEEDED: a decision about the corpus's endpoint definition, not about
           this paper. Bears directly on CLAUDE.md section 8's open
           "endpoint-comparability grouping" item.

### D-4  The in-vitro HeLa values exist but carry no identity
STATUS:    OPEN
AMBIGUITY: Figure 2C (bottom) plots in-vitro HeLa luminescence against in-vivo
           flux for the same 49 LNPs, and Methods describes the assay fully
           (15,000 cells/well, BrightGlo, N=4). The points are **unlabelled**,
           so no in-vitro value can be attributed to a lipidoid. The paper has
           in-vitro delivery data for 49 lipidoids and no way to name them.
CANDIDATES:
  (a) **Drop it** (what this folder did). The paper contributes in-vivo rows
      only.
  (b) **Attribute through the shared y axis.** Both Figure 2C panels and
      Figure 2B use in-vivo flux, so in principle each 2C point maps back to a
      2B bar. In practice the 49 bars crowd within ~0.01-0.02 decades of each
      other near 1e7 and the markers in 2C overlap -- only 30 of 49 are even
      resolvable as separate blobs (`07_fig2c_qc.json`). A non-injective map
      would put measured values on the wrong molecules, which is the AC_2025
      failure.
  (c) **Ask the authors.**
EVIDENCE NEEDED: a labelled version of Figure 2C, or the underlying table.

### D-5  `Dose_ug_nucleicacid` is not derivable
STATUS:    OPEN
AMBIGUITY: The dose is stated only as **0.75 mg/kg total mRNA**. The paper
           gives the strain, sex and age of the mice (C57BL/6 female,
           6-8 weeks) but no body weight, so micrograms cannot be computed.
CANDIDATES:
  (a) **Blank** (what this folder emitted), following the F-004 precedent where
      an unrecoverable `N:P` was blanked rather than reconstructed.
  (b) Assume a body weight (a 6-8 week female C57BL/6 is ~17-20 g, giving
      12.8-15 ug) and record the assumption on the row.
  (c) Add a `Dose_mg_per_kg` column to the schema, since mg/kg is what in-vivo
      papers actually report and the corpus keeps losing it.
EVIDENCE NEEDED: a stated body weight, or a schema decision.

### D-6  Licence: this paper is not CC
STATUS:    OPEN
AMBIGUITY: The PMC record's `<permissions>` block reads *"This file is
           available for text mining. It may also be used consistent with the
           principles of fair use under the copyright law."* That is the NIH
           Author Manuscript Collection statement, **not** a Creative Commons
           licence. CLAUDE.md section 8 requires every constituent source to be
           CC0/CC-BY-compatible before corpus release.
CANDIDATES:
  (a) Treat the extracted values as facts, which are not copyrightable, and
      release.
  (b) Check the Wiley landing page for a CC-BY licence on the version of
      record (the article may be open access there even though the PMC deposit
      is an author manuscript).
  (c) Exclude the paper from any released corpus while keeping it internally.
EVIDENCE NEEDED: the licence string off `onlinelibrary.wiley.com`. It was not
           fetched here -- see `00_acquisition.json`.

### D-7  The `Experiment_ID`
STATUS:    OPEN
AMBIGUITY: `KW_2024` is proposed and is free. But the corpus already carries
           `KW_2014` (Whitehead 2014, PMID 24969323), which is the same
           laboratory and the same lipidoid chemistry, and `curation/` contains
           a `whiteheadka_30637934` folder for a third Whitehead paper. The
           `INITIALS_YEAR` key does not distinguish laboratory from paper, and
           §2 of CLAUDE.md treats `Experiment_ID` as the ENVIRONMENT key for
           every cross-source analysis.
CANDIDATES:
  (a) `KW_2024`, treating each paper as its own environment (the convention).
  (b) Something that marks the shared laboratory, so Result 1 and Result 4 can
      distinguish "a different lab" from "the same lab, a later paper".
EVIDENCE NEEDED: a decision about what the source key is for. This is the same
           question CLAUDE.md section 8 raises under endpoint-comparability,
           and it is sharper here than anywhere else in the corpus: `KW_2014`
           and this paper share the head-code namespace (both use numbered
           amine heads from the same series) but not the tail namespace.

### D-8  Per-organ rows: liver / lungs / spleen are drawn but not recoverable
STATUS:    OPEN
AMBIGUITY: The corpus has in-vivo `Model_target` values of `liver`, `spleen`
           and `lung_epithelium`, and this paper measured all three for all 49
           LNPs. They appear only as **stacked bars on a log axis**
           (Figure 3A), which do not decompose: a segment's drawn height is
           log(cumulative) minus log(cumulative below it), so a one-pixel error
           in the bottom segment propagates multiplicatively into the top one.
           Six lipidoids additionally have a pie (Figure 3B/3C) giving exact
           fractions -- and for one of those six the printed percentages and
           the drawn pie disagree by 5 pp (C-6).
CANDIDATES:
  (a) **No per-organ rows** (what this folder did). 49 whole-body rows.
  (b) Per-organ rows for the six lipidoids with pies, as fraction x total.
      Six lipidoids x three organs = 18 rows, on a panel whose printed
      statistic does not match its own drawing.
  (c) Attempt the stacked decomposition, with the error propagation stated.
EVIDENCE NEEDED: whether the manuscript needs organ-resolved rows at all.
           Note Figure 3A also has one FEWER lipidoid than Figure 2B (C-5).

### D-9  `Model_target` for the emitted rows
STATUS:    OPEN
AMBIGUITY: Figure 2B plots the SUM over six organs **excised and imaged ex
           vivo** (liver, spleen, lungs, pancreas, heart, kidneys). The nearest
           corpus vocabulary term is `whole_body`, which elsewhere (`XD_2026`,
           `SX_2025`) means a whole-animal IVIS image of a live mouse. These
           are not the same measurement: an ex-vivo organ sum excludes
           everything not in those six organs and is not attenuated by skin.
CANDIDATES:
  (a) `whole_body` (what this folder emitted), with the distinction recorded in
      `05_values_qc.json`.
  (b) A new term such as `six_organ_sum`, which no existing row uses.
EVIDENCE NEEDED: a vocabulary decision. If (a) stands, the note should move
           onto the row, because a flag that lives only in a per-folder QC file
           is not a flag (F-001).
