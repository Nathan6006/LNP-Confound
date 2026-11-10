# DECISIONS — dahlmanje_36701517

Radmand et al., *The Transcriptional Response to Lung-Targeting Lipid
Nanoparticles in Vivo*, Nano Letters 23(3), 2023. PMID 36701517, PMC9912332,
doi 10.1021/acs.nanolett.2c04479. Licence **CC BY** (Europe PMC core record,
`license: "cc by"`, `isOpenAccess: "Y"`).

`provenance/build_dataset.py` refuses to merge a folder that has an `OPEN`
entry here. Nothing below is resolved alone.

---

### D-1 Does a paper with ONE ionizable lipid belong in a corpus keyed on ionizable-lipid structure?

All 144 formulations use **cKK-E12**. What varies is the helper lipid (18), the
PEG-lipid (2) and the molar ratio (4). The 37 rows in `06_rows.csv` therefore
carry **one distinct `IL_SMILES`** and three distinct formulations.

That is the `SB_2024` situation (2 unique lipids, flagged, kept) made sharper.
A per-experiment z-score over 37 rows that share one lipid is a contrast across
organ and cell type, not across chemistry, so the experiment contributes no
structure–activity signal to Results 1–5 and would add an `Experiment_ID` whose
source-classification features are constant.

Against that: cKK-E12 **is already in the corpus three times** — `AC_2025`
(vitro, 216 rows, `IL_name` = `CKK-E12`), `JW_2024` (vitro, 19 rows, `IL_name`
= `other`) and **`BL_2023` (VIVO, 6 rows, `IL_name` = `MD-1`)** — all three at
the identical full InChIKey `ZISVTYVLWSZJAL-UHFFFAOYSA-N`. `CLAUDE.md` §8 notes
that the in-vivo pool's only cross-source chemistry is MC3 and SM-102. Merging
this paper would add **cKK-E12 as a third shared in-vivo lipid**, and
`BL_2023`'s cKK-E12 rows are DOTAP/DOPE formulations in lung epithelium — the
same chemistry and the same target as Cat-LNP. That is a shared-lipid pair with
chemistry held constant, which is what Result 3 is made of.

**STATUS: OPEN**

---

### D-2 The pooled barcode screen: admissible in principle, not transcribed here

Supplementary Figure 6 is not a plot. It is **three data tables**, embedded as
rasters on SI pages 15, 16 and 17, giving for every LNP its barcode, LNP number,
per-organ average normalized delivery, the average of those, the ionizable
lipid, cholesterol, PEG-lipid, helper lipid, charge and all four molar ratios.
This is the barcode-to-formulation map, and it is the paper's best value source.

The values are **pool shares**. Fig S3a states it as a worked example — 3, 1 and
1 of five sequenced barcodes are printed as 60, 20 and 20 — which reconstructs
as 3/5, 1/5, 1/5 at residual **0.000** and sums to **100.000**. On real data the
cationic table's lung, liver and spleen columns sum to 97.364, 89.527 and
93.450 (shortfalls 2.636, 10.473, 6.550), all in the same direction and
accounted for by the one pooled LNP the text says exists and the table omits
(see `08_contradictions.csv` X-1). So each value depends on the other members
of its pool, which no corpus column records.

Pooling is not itself disqualifying — `LR_2023` already contributes 407
barcoded in-vivo rows. Two things are unresolved: whether these pool shares are
poolable with absolute readouts, and that **only the 22-row cationic table has
been transcribed here**. The neutral (58 rows) and anionic (50 rows) tables are
not. Transcribing them is a real job and it must carry its own validation; the
one used here is that the table's `average normalized delivery` column is the
arithmetic mean of its three organ columns to **1.7e-9 on 22 of 22 rows**, which
a single misread digit in any of 88 numerals would break.

**STATUS: OPEN**

---

### D-3 The aVHH arm (Fig S17) is a second cargo-expression readout and is not extracted

Fig S17 measures **% aVHH+ cells** and **aVHH MFI** by flow cytometry for the
same three LNPs — a protein-abundance endpoint, admissible in kind. It is a
different cargo (aVHH mRNA, not Cre), a different mouse (C57BL/6, not Ai14), a
different timepoint (1 day, not 4) and a different dose path. Emitting it
beside the Fig S13 rows would put two assays under one `Experiment_ID`, which is
the shape of error class 22.

**STATUS: OPEN**

---

### D-4 The schema has no cell-type column, so these rows duplicate on their key

The measurement is resolved to 21 cell types across 5 organs. The corpus
resolves to `Model_target` (organ) only, so 37 rows collapse onto 15 distinct
`(Experiment_ID, Formulation_ID, Model_target)` triples and the cell type lives
in `curation_note`. That means duplicate `(Experiment_ID, IL_name,
Model_target, Model_type)` keys — the count `audit.py` C7 tracks, currently
baselined at 40.

`Model_target` also needs two values the vivo corpus does not yet have:
**`kidney`** and **`heart`** (existing set: bone_marrow, brain, ear, liver,
lung_epithelium, muscle, spleen, whole_body).

**STATUS: OPEN**

---

### D-5 `Cargo_type = "Cre"` is a new vocabulary value

The existing set is FLuc, FVII, GFP, RLuc, base_editor, hEPO, peptide_barcode.
Cre recombinase is none of them. `SX_2025` pairs `Mouse_Ai14` with
`base_editor`, but that is a different chemistry and reusing the label would
assert something false. `Cre` is declared here rather than forced.

**STATUS: OPEN**

---

### D-6 cKK-E12's stereochemistry is not pinned by anything in this paper

cKK-E12 is a lysine-derived diketopiperazine with two defined α-carbons. The
paper draws no structure for it anywhere; it names a vendor and a catalogue
number (Organix Inc., O-8744). The stored string is the corpus's, which is
**flat** — no stereocentres. Every emitted row therefore carries
`confidence = medium` and `structure_assumed = 1` with a note saying so. The
alternative — writing the stereoisomer from the chemical name — would be
inventing a structure the paper does not state, and `provenance/reagent_registry.py`
has no cKK-E12 entry to take it from.

**STATUS: OPEN**

---

### D-7 Six An-LNP liver values are unreadable and are not reconstructed

The Fig S13c Liver bars are drawn as parallelograms with a ~2 px horizontal
shear over their height, and `bars.bar_top` refuses 5 of 6 on flank
disagreement. A de-shear was attempted from the panel's own furniture and
rejected (see `08_contradictions.csv` X-8). Those six values exist in Fig S15d
in background-subtracted form and are recorded in `05_values.csv` as
inadmissible, because mixing two normalisations inside one `Experiment_ID` is
error class 22. An-LNP therefore contributes 6 rows where Cat-LNP contributes 16
and Neu-LNP 15, and the imbalance is a refusal, not a finding about An-LNP.

**STATUS: OPEN**
