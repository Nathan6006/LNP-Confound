# DECISIONS — khanof_36652629

Tilstra, Couture-Senécal et al., *Iterative Design of Ionizable Lipids for
Intramuscular mRNA Delivery*, J. Am. Chem. Soc. 2023. PMID 36652629,
DOI 10.1021/jacs.2c10670.

`build_dataset.py` refuses a folder with an `OPEN` entry. Every item below is
one the pipeline detected and does not get to settle on its own.

---

### D-1 Six rounds, or one experiment?

The six value panels are the six rounds of an iterative design. They share an
endpoint (total flux, p/s, hindlimb, IM, C57BL/6, firefly luciferase) and a lab,
which is what makes pooling tempting.

They should not be pooled, and the paper's own Table S2 is the evidence rather
than our numbers: the same ionizable lipid carries a **different formulation ID**
in different panels — C3-K2-E14 is row 3 in Fig S1B and row 4 in Fig S2B,
C3-K2-E10 is row 6 then row 7, C3-K3-E10 is row 8 then row 10. Those rows differ
in size, PDI and pKa, so they are different batches of particle.

Measured, for the five lipids that appear in more than one panel:

| lipid | panels | max/min |
|---|---|---|
| C3-K2-E14 | S1B 1.711e6 · S2B 3.415e5 | **5.011x** |
| C3-K3-E10 | S3B 2.912e6 · S4B 1.565e6 | 1.861x |
| OC2-K3-E10 | S4B 1.135e7 · S4C 1.988e7 · S7A 1.184e7 | 1.751x |
| SM-102 | S1B 8.487e6 · S4B 1.135e7 · S4C 1.411e7 · S7A 1.090e7 | 1.662x |
| C3-K2-E10 | S2B 4.265e5 · S3B 4.490e5 | 1.053x |

5.011x between two panels of the same molecule, while another re-formulated
lipid moves only 1.053x, is not a spread that one scale absorbs.

Proposal: **six blocks, one per round**, which is how `05_values.csv` is
written. The alternative — one block of 18 rows — would put a 5x batch effect
inside an experiment and call it structure–activity.

STATUS: OPEN

---

### D-2 The only shared chemistry with the corpus is benchmarks

All 8 novel lipids are absent from the corpus on all three keys. The only
overlap is SM-102, DLin-MC3-DMA and ALC-0315, and all three match on full
InChIKey, skeleton and `IL_name` together.

This lands directly on the open question in `CLAUDE.md` §8: *"a source whose
only shared lipids are benchmarks that top every screen can push cross-lab
concordance toward +1 for reasons that have nothing to do with transferable
SAR."* Here it is sharper than in `XD_2026`, because this paper's benchmarks do
**not** top the screen — in Fig S7A the four bars sit within 1.945x of each
other and every comparison is marked "ns", and in Fig S4B SM-102 and the
paper's own lead read identically.

Whether these rows may enter a shared-lipid pair in Stage 5 is not ours to
decide.

STATUS: OPEN

---

### D-3 The mRNA dose is nowhere in the material we hold

Fig S8 varies the dose from about 0.1 to 1.25 µg, so the other six panels used
some fixed dose — and no caption, no table and no methods paragraph in the SI
states it. The main article is not open access (Europe PMC: `isOpenAccess=N`,
`inEPMC=N`, `pmcid=null`) and was not obtained; no paywall bypass was attempted.

`dose` is blank on all 18 rows. The options are to emit with a blank dose, or to
acquire the main article first. Reconstructing it from the lipid:RNA ratio in
Table S2 would need the injected lipid mass, which is also unstated — that is
the F-004 shape and it is not attempted.

STATUS: OPEN

---

### D-4 Fig S9A disagrees with the panels it is derived from

The authors drew two fold-increase bars, Tails 2.06 and Core 5.93. Formed from
our reads of the same panels those ratios are 3.162x and 7.251x — 1.534x and
1.223x higher.

A within-panel ratio cancels any constant row offset and any calibration-origin
error, so this is not an origin problem in our reader; and the two gaps are not
equal, so it is not a scale error either. Our S1B and S4B tops are independently
supported by check G, where all 12 of the paper's own significance marks sort
correctly on a single threshold. The remaining explanations are that the authors
computed these from something other than the plotted arithmetic means — a
geometric mean or a per-animal ratio would both run lower — or that Fig S9A
disagrees with Figs S1B and S4B.

Fig S9A supplies no row either way. What is open is whether this is enough
doubt to caveat the S1B and S4B values.

STATUS: OPEN

---

### D-5 Fig S4D's OC2-K3-E10 characterisation matches no OC2-K3-E10 formulation

Fig S4D reads (95 nm, 0.15, 97%, 6.8). Table S2 row 11, the only unmodified-mRNA
OC2-K3-E10 at N/P 5, is (102 nm, 0.12, 97%, 6.76). The size and PDI instead
match row 8 exactly, which is **C3-K3-E10** — a different lipid, and the bar
immediately to its left in the same panel.

We store row 11 and mark `formulation_match` as partial. The N/P, lipid:RNA,
size, PDI, EE and pKa carried on that row therefore rest on a two-of-four match.
Which half is wrong is a question for the source, not for us — this is the
`AC_2025` shape, where the first confident repair moved 419 measured values onto
the wrong molecules.

STATUS: OPEN

---

### D-6 Fig S4C uses a different cargo

Fig S4C is nucleoside-**modified** firefly luciferase mRNA; the other five
panels are unmodified. Table S2 agrees, giving modified-mRNA formulations 13 and
14 for exactly these two bars. It is kept as its own round, `R5_core_modmRNA`,
with `cargo` recorded on the row.

Whether a modified-mRNA block belongs in the same corpus pool as unmodified-mRNA
blocks is a corpus-level question.

STATUS: OPEN

---

### D-7 Six lipids are stereoisomer mixtures stored as flat SMILES

The E-series tails come from ring-opening a **racemic** 1,2-epoxyalkane, four
times per molecule. The product is a mixture and no single configuration is
correct, so the stored SMILES has 4 unassigned stereocentres (0 for the two
A-series lipids, which have none). That is the corpus norm, and it is stated in
`curation_note` on every affected row rather than left to be inferred.

It also means C10's stereo comparison carries no information for this folder,
which `reference_structures.csv` says in each anchor's note. Flagging it as a
`confidence` caveat was considered and rejected: the connectivity is pinned by
the synthesis, the drawn structure and an 8/8 formula match, and lowering
`confidence` would mislabel a connectivity-certain structure.

STATUS: OPEN
