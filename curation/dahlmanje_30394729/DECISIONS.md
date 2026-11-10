# DECISIONS — dahlmanje_30394729 (Sago et al., JACS 2018, PMID 30394729)

Proposed `Experiment_ID`: **CS_2018** (first author Sago CD; the folder is named
after the senior author, Dahlman JE).

`build_dataset.py` refuses a folder with an `OPEN` entry. D-1 is open, so this
folder cannot merge until a human settles it. That is the intended state.

---

### D-1 — 7C1 has no single SMILES. Does a polymeric ionizable lipid belong in a corpus keyed on structure?

STATUS: OPEN

The paper's only ionizable lipid is **7C1**: branched polyethylenimine, Mn ≈ 600
(1 molar eq), ring-opened with **14 molar eq** of a C15 1,2-epoxide (R = C13),
90 °C / 72 h / EtOH (SI Fig S1B). The product is a polydisperse, statistically
substituted polymer. Fig 3A and Fig 5A draw **one representative arm** with the
PEI backbone nitrogen replaced by a diethylamino group; that cartoon is
`CCCCCCCCCCCCCC(O)CN(CC)CC` (C19H41NO), and it is not the lipid.

Facts that bear on it:

- `data/data_vitro.csv` and `data/data_vivo.csv` have **0 rows with a blank
  `IL_SMILES`** — 29,195 and 1,803 rows respectively. The schema has no way to
  say "polymer".
- The project already takes this position for PEG-lipids: `reagent_registry.py`
  excludes them deliberately, because "the block length is a distribution, so no
  single SMILES is right". 7C1 is the same situation on the ionizable lipid.
- Writing a representative arm, or an arbitrary n-mer, into `IL_SMILES` would be
  guessing a structure, which the protocol forbids outright.

The three live options: (a) emit nothing — what this folder does; (b) extend the
schema with an explicit polymer representation and a flag, which is a corpus-wide
change and not a curation's call; (c) emit the five individually measured
formulations with the representative arm as `IL_SMILES` and
`structure_assumed = 1`, which puts a 299 Da small molecule where a ~2 kDa
polymer belongs and would be read by every descriptor column as fact.

I recommend (a) now and (b) as a separate corpus decision. Note that resolving
this unlocks only **five** rows (D-2), so the cost/benefit is poor unless the
schema question is being answered for its own sake.

---

### D-2 — the pooled screen: 146 measured formulations, zero attributable values

STATUS: OPEN

Library 1 pooled 115 of 120 formulated LNPs; library 2 pooled 31 of 40. Each
pooled LNP has a `Normalized DNA Delivery (%)` value. **Not one of them is
attributable**, for two independent reasons:

1. Fig 3E / 5D plot them as unlabelled scatter points (115 and 31 identical grey
   circles), and Fig 3F/3G/5E/5F re-plot the same values against diameter, also
   unlabelled. There is no key and no stated ordering.
2. **Which** LNPs were pooled is never published. 5 of 120 and 9 of 40 failed the
   20–200 nm stability gate and are not named, so even a perfect digitisation
   could not be aligned to the composition tables.

Separately, the readout itself is a question: the quantity is DNA barcode
abundance recovered from FACS-sorted ICAM-2-low cells. The functional element is
the **gate**, not the measurement. The corpus's one barcoded experiment,
`LR_2023` (407 of 1,803 in-vivo rows), carries `Cargo_type = peptide_barcode`
with `Experiment_method = protein_adundance_normalized` — a *translated* readout.
Pooling is not the objection; the readout is. If a DNA-barcode-recovery class is
ever admitted, it needs its own `Experiment_method` value and the decision is
corpus-wide, not per paper.

Nothing here is actionable without new information from the authors (the
barcode-to-LNP map and the stable subsets).

---

### D-3 — library 2 LNP 36 and LNP 40 sum to 110 mol%. Which number is wrong?

STATUS: OPEN

Both read `60 / 45 / 5 / 0` (7C1 / cholesterol / PEG / DSPC) in SI Fig S2F — the
same ratio in both PEG blocks, at a stride of 4, i.e. one error in a block pasted
twice (`08_contradictions.csv` C-3). Two repairs are supported by different
evidence and they disagree:

- **Series pattern.** Rows 33/34/35 hold cholesterol + PEG = 40 exactly
  (38.5 + 1.5, 37.5 + 2.5, 36.5 + 3.5). Row 36 continuing that series is
  `60 / 35 / 5 / 0`: the cholesterol figure is the typo, PEG stays at 5.
- **Fig 5H.** The 'Formulated' pie prints 40 / 40 / 20 where the table gives
  45 / 20 / 35 (or 45 / 40 / 15 inclusive). Reconciliation needs exactly **two**
  LNPs to move from the low PEG band into >20 — and these are exactly two LNPs.
  That implies the PEG figure is the typo and should be above 20.

The first is a local pattern, the second is an independent rendering of the whole
table that is exact on five other panels. I am not choosing: this is the
`AC_2025` shape, where a confident internally-consistent repair put 419 measured
values on the wrong molecules. It changes nothing in `data/` today (no rows are
emitted), so it costs nothing to leave open.

---

### D-4 — the PEG-lipid identities

STATUS: OPEN

The paper names `C14PEG2000`, `C16PEG2000`, `C18PEG2000`, `C14PEG3000`,
`C18PEG3000` and names no catalogue product. Fig 3B draws a diacyl
glycero-phospho-ethanolamine–PEG with `(OCH2CH2)45OCH3`; 45 × 44.05 = 1982.4 Da,
0.9 % from the printed 2000 (never multiplied out in the paper). By that skeleton
C14/C16/C18 would be DMPE-/DPPE-/DSPE-PEG2000, and `DMPE-PEG2000` and
`DSPE-PEG2000` already exist in the corpus `PEG_name` vocabulary. **That is an
inference from a drawn skeleton, so it is recorded and not applied.** Only
matters if D-1 is resolved in favour of emitting.

---

### D-5 — licence

STATUS: OPEN

Europe PMC: `isOpenAccess = N`, `license = null`. The PMC copy is an NIH author
manuscript whose permissions read "available for text mining ... fair use". That
is **not** CC0/CC-BY-compatible, so under `CLAUDE.md` §8 this source would block
corpus release even if D-1 were resolved. Flagged for the release decision, not
for this curation.
