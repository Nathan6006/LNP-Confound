# DECISIONS — Sago et al. 2018, PMID 30275336

Queued, never resolved alone (PIPELINE_V2 §3). Every `### D-n` carries a plain
`STATUS:` line.

---

### D-1. Is a POOLED COMPETITIVE barcode readout comparable in kind to the corpus's per-LNP measurements?

**What the readout actually is, in the paper's own words.** Methods, verbatim:
*"Counts for each particle, per tissue, were normalized to the barcoded LNP
mixture we injected into the mouse. This 'input' DNA provided the DNA counts and
was used to normalize DNA counts from the cells and tissues."* Results:
*"we ranked LNPs by calculating the 'normalized delivery' of each barcode.
Normalized delivery is analogous to counts per million in RNA-sequencing data."*
Units as printed: **Normalized DNA Delivery (%)**.

**Measured, not inferred.** `stage04_readout.py` reconstructs the formula from
the paper's own worked example (SI Fig S1b, p.2) and reproduces every printed
cell: **6/6 column totals, 18/18 table-2 cells, 3/3 input shares, 18/18 table-3
cells**. The formula is

    normalized_delivery_i = [ (count_i / Σ_j count_j) / (input_i / Σ_j input_j) ]
                            / Σ_k[ same ]  × 100 %

The paper prints `Total (%) = 100` under **both** normalised tables, and the
unrounded reconstruction gives column sums of **exactly 100.000000, residual
0.0**. Three seeded defects each move the agreement count
(`wrong_raw_count` → 15/18; `drop_input_normalisation` → 1/18;
`no_renormalisation` → 0/18 and residual 213.37).

**So pooling is not incidental — the value is compositional.** Each LNP's number
is its *share* of a pool that is constrained to sum to 100%. Change the pool
membership and every member's value changes. The three in vivo screens used
three different pools (71, 78 and 106 members), so their values live on three
different simplices and are not on a common scale even within this paper.

**But the corpus has already admitted a pooled screen, and that must be faced
honestly.** `data/` carries an `Experiment_batching` column whose vocabulary is
`{barcoded, individual}`, and `LR_2023` (PMID 37127709) contributes **407 vivo
rows** at `Experiment_batching = barcoded`, `Cargo_type = peptide_barcode`.
A pooled, multiplexed in vivo screen is therefore already in the corpus and
declared as such. **Pooling alone is not grounds for refusal here**, and a
refusal on that basis would contradict the project's own precedent.

**Recommendation.** Treat pooling as a *declared property*, not a disqualifier,
and decide separately whether the compositional (sum-to-100%) normalisation is
a further step beyond `LR_2023`'s — it may be, since a share-of-pool is
rank-preserving within a pool but not comparable across pools. This question is
**moot for this paper** (D-2 blocks it outright) but will recur for the next
Dahlman-lab screen, and it bears on whether `LR_2023`'s 407 rows belong in
Result 3's shared-lipid concordance.

STATUS: OPEN

---

### D-2. ZERO structures are determinable. This, not pooling, is what blocks the paper.

Eight lipid-amine compounds are named: **7C1, 73, 78, 80, 83, 88, 90, 92**.
None has a determinable structure.

| compound | amine | tail | why not determinable |
|---|---|---|---|
| 7C1 | PEI600 | C13 epoxide | polydisperse branched polymer |
| 73, 78, 80 | PEI1200 | C12 epoxide / acrylate / acrylamide | polydisperse branched polymer |
| 83, 88, 90 | triethylenetetramine | C12 epoxide / acrylate / acrylamide | amine is determinable; **degree of substitution never stated** |
| 92 | macrocyclic tetraamine | C10 epoxide | ring size unconfirmed **and** degree of substitution never stated |

SI Fig S2b and S1j draw each product as a **single repeat-unit cartoon**: one
tertiary nitrogen bearing one arm, with the backbone truncated to two stubs.
That is not a molecule. The paper itself attributes 7C3's batch-to-batch
variability to *"the polymeric backbone of 7C3 and the resulting heterogenous
chemical structure"*, and concludes that *"future FIND screens should focus on
compounds with defined chemical structures"* — the authors are telling us these
are not defined structures.

`data/data_vitro.csv` + `data_vivo.csv` contain **0 blank `IL_SMILES` in 30,998
rows**. A row with no structure cannot be emitted, and inventing one is
forbidden. `reference_structures.csv` therefore carries **8 rows with an empty
`reference_SMILES`**, so the absence is visible rather than silent.

STATUS: OPEN

---

### D-3. Even the individually-measured subset carries zero ionizable structural variance.

The admissible-looking subset is the lead-LNP validation: **7C2 and 7C3**,
injected one LNP per mouse, with **25 printed numerals** (14 in the SI Fig S4a
table, 11 in body prose). These are genuinely individual measurements with a
declared n, dose, route and timepoint, and the endpoints map onto corpus
vocabulary (`editing_efficiency_normalized`, `protein_adundance_normalized`).

They are still inadmissible, for the D-2 reason — and for a second reason worth
recording: **7C2 and 7C3 share the same ionizable compound, 7C1.** They differ
only in helper lipid (18:1 Lyso PC vs DOPE) and mol ratio. So even if 7C1 had a
SMILES, this paper would contribute an experiment of 2 rows with **identical
`IL_SMILES`** — zero structural variance, and therefore nothing at all for a
structure–activity or cross-lab concordance analysis.

STATUS: OPEN

---

### D-4. The screens' values are unattributable, independently of D-2.

Normalized delivery for the 71 + 78 + 106 pooled in vivo LNPs (and the 54
in vitro) is published **only** as unlabelled grey points (Fig 2E, SI Fig S6a–c,
whose x-axis is the *bioreplicate*, not the LNP) and as an unlabelled
Euclidean-clustered heat map (Fig 2F). Fig 2F does carry a numeric colour scale
(`Log2 Normalized Delivery (%)`, 0–3) — so the refusal is **not** "no numeric
scale"; it is that the columns bear no LNP labels and their order is the
clustering's, which is unrecoverable.

**The design supports 309 per-LNP values. The publication attributes 0 of them.**
No digitiser was written and none should be: there is nothing to key a measured
pixel onto.

STATUS: OPEN

---

### D-5. The licence is not release-compatible.

The article is **CC BY-NC-ND 4.0** (PMC JATS `<permissions>`, verbatim).
`CLAUDE.md` §8 requires constituent-source licences to be CC0/CC-BY-compatible
before corpus release; **NC and ND both fail that test**. No author repository
exists, so there is no second licence to fall back on. Recorded even though zero
rows are proposed, because the next Dahlman-lab PNAS paper will hit the same
clause.

STATUS: OPEN

---

### D-6. `C14-PEG2000` means a different molecule here than it does in the corpus.

MEASURED (`stage03_chaincheck.py`): the PEG-lipid drawn in SI Fig S2b has a
**14-carbon acyl chain** (13 contiguous bond quads, 5.92 pt stride, uniformity
gate passed) on a **phosphoethanolamine-carbamate-PEG2000** head — i.e. 14:0
PEG2000 PE / DMPE-PEG2000. The corpus resolves the name `C14-PEG2000`
unambiguously to the **glycerol-ether DMG-PEG2000**, a different molecule.

Nothing is emitted from this paper, so nothing is wrong today. It is queued
because it is error class 30 in a column nobody guards: a future curation that
keys a PEG name against the corpus would attach the wrong structure, and the
guarded lookup in `corpus_components.py` would *not* catch it — the name is
unambiguous in the corpus, just wrong for this paper.

STATUS: OPEN

---

### D-7. Seven formulations do not sum to 100 mol%.

418 of 424 sum to exactly 100. The exceptions are the in vitro library's LNP #1
(101.5) and Screen 3's LNPs 65/81/97/113 (101, the **same ratio replicated
across four helper lipids at a stride of 16** — one error in a block pasted four
times) and 143 (99); LNP 141's helper cell is blank while its other three sum to
100. Verified to be the paper's, not the extraction's: the table grid and an
independent text extraction of the same page return identical digits.

Printed values stand, unaltered, with `molpct_sum` carried in
`01_formulations.csv`. No correction is applied — the parallel blocks make a
typo likely but that is an inference.

STATUS: OPEN
