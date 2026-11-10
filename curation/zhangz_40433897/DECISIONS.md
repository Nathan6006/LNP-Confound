# DECISIONS — zhangz_40433897

Cui, Zhang, Shi, Zhao, Zhou, Liu, Hu, Hu, Kong, Zhang.
*Iterative Screening of Vitamin E-Based Functional Lipid Nanoparticles for mRNA
Delivery.* ACS Nano 2025;19(22). PMID 40433897. DOI 10.1021/acsnano.5c01378.

Every entry below is for the user, not for me. `build_dataset.py` refuses a
folder with an `OPEN` entry, which is the intended state of this folder: the
75 rows in `06_rows.csv` must not reach `data/` until D-1 is settled.

---

### D-1 — No VNxCyOz structure in this paper is PINNED. Three are PROPOSED. Do they go into the corpus?

The paper draws **no** VNxCyOz structure anywhere. Figure S1b draws the reaction
with a generic `R-NH2`; Figure S1a draws 35 commercial amines, each labelled
with a CAS registry number; Figure S2's caption prints one exact mass per
compound. So every structure here is an inference from three measured things:

1. the scaffold `R-NH2 + 2 x VE-acrylate -> R-N(CH2CH2CO-O-tocopheryl)2`, whose
   intermediate VE-AC reproduces the printed `[M+Na]+` calc of 507.3809 to
   **-0.03 mDa** from our own build;
2. the exact formula of each product, printed in Figure S2;
3. the 34 distinct amines drawn in Figure S1a.

(1) + (2) fix the core amine's **formula** exactly. (3) then says which drawn
reagents carry that formula. The result, per compound:

| compound | core formula | drawn amines with that formula | connectivity determinate? | proposed |
|---|---|---|---|---|
| VN2C2 | C6H16N2 | **1** — CAS 100-36-7 | yes (one N with 2 H) | yes |
| VN3C2 | C7H18N2 | 2 — 104-78-9, 10061-68-4 | 104-78-9 only | yes, by the CH2 ladder |
| VN3C1 | C5H14N2 | **1** — CAS 109-55-7 | yes | yes |
| VN2C4 | C6H14N2 | 2 — 7154-73-6, 694-83-7 | neither uniquely | **no** |
| panel (g) | C7H16N2 | 1 — 51387-90-7 | yes, but its NAME is disputed (X-1) | **no** |
| VN3C2O | C7H18N2O2 | **0** | — | **no** |
| VN2C1 | C4H12N2 | **0** | — | **no** |

**What makes the three proposals strong.** The measured calc masses put
VN2C1 -> VN3C1, VN2C2 -> VN3C2 and VN2C4 -> panel (g) exactly one CH2 apart,
3/3 within **0.05 mDa**, which the paper never states. That is consistent with
`VNx` = the number of carbons between the two nitrogens and `Cy` = the head
amine type, and under that reading C1 = dimethylamino and C2 = diethylamino,
with **both** rungs of each ladder (109-55-7/104-78-9 and 100-36-7/104-78-9)
actually drawn in Figure S1a. Our builds of the three reproduce the printed
`[M+H]+` calc to **-0.08, -0.03 and -0.03 mDa**.

**What makes them only proposals.** Two of the seven cores — including the
paper's own LEAD compound VN2C1 — have **no** amine of the right formula in the
drawn set at all (X-4). So "the only drawn amine with this formula" is an
argument whose premise the paper itself falsifies: the drawn set is not a
complete inventory of the cores used. A C6H16N2 core could equally be
N,N'-diethylethylenediamine or 1,6-hexanediamine, neither drawn, and the exact
mass cannot tell them apart. Mass cannot see connectivity.

**Options.**

- **(a)** Merge the 75 rows with the proposed structures, `structure_assumed=1`,
  `confidence=low`, and the curation note already on each row. The chemistry is
  a three-way inference, not a reading.
- **(b)** Merge nothing from this paper. Zero structures are read from a
  drawing, which is the standard C10 exists to enforce.
- **(c)** Merge only the Figure S4b block (58 rows, one lipid, VN2C2 — the one
  compound with a unique formula match, a determinate connectivity, and a
  drawn +CH2 partner), and drop the Figure S5 blocks where four of six lipids
  have no structure.

I have not chosen. Emitting (a)'s rows while leaving this OPEN is deliberate:
the rows exist so the decision can be made against them, and the merge gate
stops them.

STATUS: OPEN

---

### D-2 — Figure S4b's printed DMG-PEG correlation has the wrong sign. Report it, or treat our digitisation as suspect?

The inset table of Figure S4b prints four coefficients. Against Table S2's mol%
(PDF text layer, tier 2) and our 58 readable bars (tier 5):

| factor | printed | our Spearman | signed delta | delta on magnitude |
|---|---|---|---|---|
| VN2C2 | -0.636 | -0.6275 | +0.0085 | 0.0085 |
| DDAB | -0.501 | -0.5170 | -0.0160 | 0.0160 |
| Chol | +0.727 | +0.7215 | -0.0055 | 0.0055 |
| DMG-PEG | **-0.519** | **+0.5324** | +1.0514 | **0.0134** |

All four magnitudes reproduce to **<= 0.016**. Three signs reproduce. The
fourth does not, and no construction from Table S2 gives DMG-PEG a negative
coefficient (mol% +0.532, rank +0.532, share-of-block +0.299, Chol/DMG-PEG
ratio -0.296 — wrong magnitude). The design is collinear in the direction that
forces the positive sign.

This also settles something the paper never states: the metric is **Spearman**,
not Pearson (Pearson gives -0.623 / -0.391 / +0.760 / +0.581, only two of which
land within 0.05).

The question for the user is whether a published sign error gets written into
Methods as a reason to treat this source's reported statistics as unreliable,
or whether it is noted and otherwise ignored. It does **not** affect any value
we emit.

STATUS: OPEN

---

### D-3 — Ten of 118 bars are refused. Accept the gaps, or find another route?

`curation/lib/bars.bar_top` refuses 9 bars outright and one more is a
zero-height control; all refusals were respected and none was fixed by lowering
a threshold. Causes, each diagnosed by looking at the pixels:

| panel | bar | cause |
|---|---|---|
| FigS5a | VN2C3 @ PEG-0.5 | a replicate ring straddles the LEFT top corner; flanks disagree by 3 px |
| FigS5a | VN2C2 @ PEG-1 | same, 4 px |
| FigS4b | LNP30, LNP54 | the right flank's top rule is 2 rows where the full-width rule is 3 |
| FigS7 | Lip3000 | replicate ring on the left corner, 3 px |
| FigS8a | **VN2C1** | a replicate dot cluster spans the whole right flank; flanks disagree by 6 px |
| FigS8a | Lip3000 | 3 px |
| FigS8b | VE+SM102 | 5 px |
| FigS8b | PBS, FigS13b PBS | zero-height control bars; the detector's baseline guard fires, correctly |

The one that costs something is **FigS8a's VN2C1**, because Figure S8a is the
VN2C1-vs-Lip3000 comparison and Figure S8b is the VN2C1-vs-SM-102 comparison;
losing S8a's VN2C1 makes the S8a/S8b repeat ratio uncomputable (it is the only
same-axis repeat in the whole SI). Reading it would mean taking the right
flank's row over the left's, which is precisely "returning the better-looking
flank" — the move the module refuses and the move that produced F-018.

None of the ten is corpus-admissible on other grounds anyway (controls, or
lipids with no structure), so nothing currently depends on this.

STATUS: OPEN

---

### D-4 — Figure S4b's 60 formulations and Figure S5's 18 are different rounds. One block or six?

`01_blocks.csv` emits **six** blocks, one per round, rather than pooling:

- `R1_formulation_DoE` — 60 prescriptions of ONE lipid (VN2C2), 4T1-Luc, siLuc.
- `R2_structure_x_PEG_4T1` — 6 structures x 3 PEG, 4T1-Luc, siLuc.
- `R3_structure_x_PEG_MC38` — the same 6 x 3, MC38-Luc.
- `R4_lipid_to_pDNA_ratio`, `R5_VE_TPGS_substitution`, `R6_storage`.

R1 and R2 share a cell line, a cargo and a readout, and share the lipid VN2C2 —
but **no formulation is shared between them**: R1's compositions are in Table
S2 in mol%, while R2's three levels are labelled only `PEG-0.5 / PEG-1 / PEG-2`
and no mol% is given for any component. So no repeat ratio is computable
between the two rounds and they must not be pooled into one standardisation
block.

If the user wants a single `Experiment_ID`, that is a decision about what an
"experiment" means here, not something the source settles.

STATUS: OPEN

---

### D-5 — `Dose`, `timepoint` and all formulation mol% for Figures S5, S7, S8 and S13 are not recoverable from the SI.

The SI states no dose of nucleic acid, no assay timepoint, and no component
mol% for any figure except S4b (Table S2). The main article is behind the ACS
paywall and this run was forbidden the browser, so those fields are blank, not
guessed. `05_values.csv` carries `NOT STATED IN THE SI` rather than an empty
cell, so the absence is visible.

If the main text is obtained later, these fields can be filled — but a row
merged now would carry blanks that look like "not applicable".

STATUS: OPEN
