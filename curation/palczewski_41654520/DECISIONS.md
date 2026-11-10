# DECISIONS — palczewski_41654520

Zhang et al., *A combinatorial synthetic strategy for developing genome-editing
protein-delivery agents targeting mouse retina*, Nat Commun 2026.
PMID 41654520 · PMC12992809 · doi 10.1038/s41467-026-69077-w · CC BY 4.0

Queued for the user. The pipeline detected that a source is silent or
self-contradictory; it does not get to decide what the source meant. Every
entry below is `OPEN`, so `build_dataset.py`'s merge gate will refuse this
folder until they are settled — which is the intended state: **nothing here has
been merged into `data/`.**

---

### D-1
STATUS: OPEN

AMBIGUITY: The CPP5 pentapeptide's **sequence**, which fixes the structures of
`CBB18` and `CBB19`. The SI synthesis heading (p29) writes **KLVPM**; the main
text Results writes **KLPVM**. Both are permutations of the same five residues,
so both give C27H50N6O6S and identical exact masses — the HRMS check in
`stage07_qc.py` P1 passes on *both* candidates to within 1e-4 (`CBB18_KLVPM`
and `CBB18_KLPVM` are both `exact`). This is the textbook case of mass being
blind to connectivity.

CANDIDATES:
- `KLVPM` — the SI, which is where every other structure in this curation came
  from and which is the document closest to the bench.
- `KLPVM` — the main text, which is the version a reader of the paper proper
  sees, and which is the only place the word "sequence" appears.

EVIDENCE NEEDED: anything that orders the Val and Pro — the peptide supplier's
catalogue entry, an MS2 fragmentation series, or a prior paper defining "CPP5".
The authors' own earlier work on CBB-protein delivery (ref 31 in this paper)
may name it. Failing that, the honest outcome is to emit `CBB18`/`CBB19` with
no SMILES and keep only their values, which is what this folder currently does.

Current state: **no SMILES emitted** for `CBB18` or `CBB19`. Both candidates
are in `03_enumerated.csv` as `CBB18_KLVPM`, `CBB18_KLPVM`, `CBB19_KLVPM`,
`CBB19_KLPVM`. Their four value rows carry `confidence = low` and
`structure_assumed = 1`.

---

### D-2
STATUS: OPEN

AMBIGUITY: Where `CBBox11` is hydroxylated. It is the DMSO-storage oxidation
product of `CBB11` (bis-linoleoyl). The SI says only that the material is **two
chromatographically separated isomers** whose MS2 "suggested hydroxylation
events on one of the polyunsaturated fatty acid moieties" (Supplementary
Figure 2). Which arm, and which of the several allylic carbons, is not
determined — and "two isomers" means the isolated material is a mixture, not a
compound.

CANDIDATES: any mono-hydroxylated linoleoyl arm. The allylic positions flanking
the 9,12-diene (C8, C11, C14) are the chemically likely ones, but the paper
assigns none of them and we will not pick.

EVIDENCE NEEDED: the MS2 fragment assignment behind Supplementary Figure 2c/2f,
or an NMR of the separated isomers. Neither is in the paper.

Current state: **no SMILES emitted.** `CBBox11` appears in `03_enumerated.csv`
with an empty `canonical_SMILES`, `confidence = low`, `structure_assumed = 1`
and `positional_ambiguity = 1`. Its Fig 2d value (mean 2.057% RFP conversion,
n=2) is kept, because the value is real whatever the isomer ratio is.

---

### D-3
STATUS: OPEN

AMBIGUITY: `CBB7`'s HRMS entry contradicts itself and we do not know which half
is wrong. The SI gives the formula `C95H150N9O7S2 3+` and the m/z `531.3706;
found 531.3704`. That formula computes to **531.0361**, so the printed m/z is
0.3345 high — exactly 1.0035 Da spread over a charge of 3.

CANDIDATES:
- The **formula is right and the printed m/z is a typo.** This is what we
  emitted. It is supported independently of the formula string by homolog
  spacing: `CBB4` and `CBB7` differ by exactly C12, so their 3+ m/z must differ
  by 48.0000, and 579.0361 − 48.0000 = 531.0361. `CBB2`, `CBB4` and `CBB14`,
  the other three quaternised members, are all internally consistent.
- The **m/z is right and the structure has one more hydrogen than we think.**
  Uncomfortable, because the *found* value tracks the printed calculated one to
  0.4 ppm, which is not what a typo in a calculated value usually looks like.
  We could not construct a chemically sensible species at 531.3706.

EVIDENCE NEEDED: the raw CBB7 spectrum. The Source Data workbook's "Synthesis
of CBB" sheet carries LC-MS traces per compound and may contain it; we did not
mine 92,168 rows of m/z-versus-intensity for this, and that is the obvious next
step if the decision matters.

Current state: `CBB7` emitted on the formula, flagged in
`08_contradictions.csv` row 4. **No tolerance was widened to make it pass**;
`stage07_qc.py` P1's gate names CBB7 explicitly as the one exception.

---

### D-4
STATUS: OPEN

AMBIGUITY: Whether this paper belongs in the corpus at all, and if so, which of
its four assay blocks become `Experiment_ID`s.

The corpus's in-vitro pool is **mRNA / nucleic-acid delivery**. **There is no
mRNA anywhere in this paper.** The cargoes are:

| block | cargo | cell line | readout as printed | n |
|---|---|---|---|---|
| Fig2c | Cre recombinase **protein** | CS = HEK293-loxP-GFP-RFP (human) | RFP conversion (%) | 4 |
| Fig2d | Cre recombinase **protein** | same | RFP conversion (%) | 2 |
| Fig4c | **ABE8e RNP** (base-editor protein + sgRNA) | rd12 reporter, NIH/3T3 (mouse) | eGFP/mCherry (%) | 4 |
| Fig4d | **ABE8e RNP** | same | eGFP/mCherry (%) | 3 |

Fig2c/4c are flow cytometry; Fig2d/4d are fluorescence microscopy, and the
*same* compound disagrees across the pair by 0.55× and 1.86× (`08_contradictions.csv`
row 9). So they are four blocks, not two, and certainly not one.

CANDIDATES:
- **Exclude entirely.** Protein and RNP delivery is a different endpoint class
  from mRNA transfection; pooling it would mean the corpus's "delivery" label
  spans cargoes that nothing says are comparable.
- **Include all four as separate `Experiment_ID`s** (e.g. `KP_2026`,
  `KP_2026_2`, `KP_2026_3`, `KP_2026_4`), with `Model_type` carrying the cargo
  and the instrument.
- **Include only the two flow-cytometry blocks** (Fig2c, Fig4c), which are the
  22-compound full-library panels at n=4, and leave the microscopy panels out.

EVIDENCE NEEDED: a ruling on endpoint comparability — the same open question
`CLAUDE.md` §8 already lists as "Endpoint-comparability grouping". Note what
the paper would buy if admitted: Fig2c and Fig4c screen the **identical 22
compounds** with two different cargoes in the same lab, which is a clean
within-lab, same-chemistry, cargo-swap concordance point — the same kind of
rung as `molea.csv`'s hepatocyte-vs-chondrocyte 0.30, and directly relevant to
Result 3.

Current state: `05_values.csv` holds all 64 library-member rows with the panel,
cargo, instrument, cell line, ratio, dose and timepoint on every row. Nothing
is merged.

---

### D-5
STATUS: OPEN

AMBIGUITY: Which salt form to store for `CBB22`. It is commercial Coomassie
Brilliant Blue G-250, used as a library member in its own right (it has values
in Fig 2c and Fig 4c). We emit the **neutral free acid**, C47H49N3O7S2,
832.06 g/mol. The SI weighs it as "854 mg, 1.0 mmol", which is the **mono-sodium
salt** — 832.06 + (Na − H) = 854.04, agreeing to 0.005%.

CANDIDATES:
- **Free acid** (emitted). Consistent with every other free-sulfonic-acid
  member of this library — `CBB5`, `CBB16` and all ten `CBBZ` compounds are
  reported by HRMS as `[M+H]+` of the neutral free acid, and those eleven
  structures all reproduce exactly on that basis.
- **Mono-sodium salt.** What was actually in the bottle, and what the only
  independent mass evidence for `CBB22` points at.

EVIDENCE NEEDED: the corpus convention for a reagent sold as a salt. This
changes `formal_net_charge` and `ExactMolWt` for one row in each of two blocks,
and nothing else. `stage07_qc.py` P2 carries a per-row `weighed_form` field so
the weighing check is correct either way — it does not depend on this decision.

---

### D-6
STATUS: OPEN

AMBIGUITY: How `Experiment_value` should be normalised if this folder is
merged. This curation emits **raw printed units** — "RFP conversion (%)" and
"eGFP/mCherry (%)" — and leaves `value_normalization_unit` blank, because no
z-scoring was done here.

CANDIDATES:
- One z of log10 **per assay block** (the §3 convention, `value_normalization_unit
  = experiment`), giving four blocks at mean 0 / sd 1.
- One z across Fig2c+Fig2d, and one across Fig4c+Fig4d, i.e. per cargo — this
  would pool flow with microscopy, and `08_contradictions.csv` row 9 is the
  argument against it.

EVIDENCE NEEDED: the D-4 ruling first; the normalisation question only arises
if the blocks are admitted. Noted here so that a blank
`value_normalization_unit` on 64 rows is read as "this folder did not z-score",
not as "nobody looked".

Every value is strictly positive (minimum 0.0003 on `CBBZ4` in Fig 4d), so
log10 is defined throughout and no flooring is needed.
