# DECISIONS — `dahlmanje_29489381`

Paunovska, Sago, Monaco et al., *A direct comparison of in vitro and in vivo
nucleic acid delivery mediated by hundreds of nanoparticles reveals a weak
correlation*, Nano Lett 18(3):2148–2157 (2018). PMID 29489381, PMC6054134.

**This curation emits 0 corpus rows.** The decisions below are queued, not
resolved. The merge gate blocks on any `OPEN` entry, which is correct here:
nothing should merge from this folder until a human has ruled on D-1 and D-2.

---

### D-1
STATUS: OPEN

AMBIGUITY: Whether this paper can enter a corpus keyed on ionizable-lipid
structure at all. Not one of its 15 lipomers has a determinable SMILES.
Seven of fifteen (`7C1`, `104`, `68`, `69`, `78`, `79`, `106`) are built on
PEI 600 or PEI 1200 — branched poly(ethylenimine) of nominal MW 600 or 1200,
which is a polydisperse mixture, not a molecule. The SI's own schemes concede
this by drawing `PEI 600` as diethylamine with the label written underneath.
For the remaining eight the amine is a defined reagent (spermine, spermidine,
triethylenetetramine, 1,4,8,12-tetraazacyclopentadecane) but the **degree of
substitution is never stated** — spermine offers 6 N–H and the paper never
says how many carry a tail, so the product ranges over roughly 442–1,653 Da.

CANDIDATES:
(a) exclude the paper entirely;
(b) include it as a *formulation* source with a blank `IL_SMILES` and the
    lipomer id as `IL_name`, accepting that it cannot contribute to Result 3
    (shared-lipid concordance) or to any fingerprint-based analysis;
(c) include only the eight non-PEI lipomers with an assumed fully-substituted
    structure flagged `confidence = low` + `structure_assumed = 1`.

EVIDENCE NEEDED: (c) is a guess, not a flag — the molecular formula itself is
open, which is a different thing from an uncertain regiochemistry, so my
recommendation is against it. Settling (a) vs (b) is a project-level question
about whether the corpus admits rows with no structure. Note that `7C1` and
`104` are *not* currently in the corpus under any name (checked both pools),
so nothing here collides with existing rows.

---

### D-2
STATUS: OPEN

AMBIGUITY: Whether the 87 rank rows of SI Fig S4e/f/g are usable as a readout
even in principle. They are the only per-LNP-identified numbers in the paper
we can reach, but they are (i) **ranks**, not delivery values; (ii) **censored**
— the in vivo columns list only top-ranked particles, so a blank means
"outside the listed set", not "not measured"; and (iii) **selected on the
outcome** — the Fig S4 caption says the listed rows are the particles coloured
in Fig S4d, i.e. chosen because they did well in vivo.

CANDIDATES:
(a) do not use them at all;
(b) use the *in vitro* columns only (iMAEC / RAWs), which are uncensored ranks
    within the full stable pool (observed 1–95 of 112, 1–53 of 105, 1–64 of 64)
    and so are a genuine partial ranking — but still only over an
    outcome-selected subset of particles;
(c) use them as a rank-only auxiliary dataset outside `Experiment_value`.

EVIDENCE NEEDED: whether the project wants a rank-valued, outcome-selected
source in the corpus at all. Moot unless D-1 resolves in favour of inclusion.

---

### D-3
STATUS: OPEN

AMBIGUITY: One LNP's lipomer mole percent is printed two ways. SI Fig S4f
(p20, endothelial sub-table) prints **36** for lipomer 58 / cholesterol 60 /
PEG 5 / helper 0 / diameter 31.2 nm; SI Fig S3d LNP #54 and Fig S4f's own
macrophage sub-table both print **35** for the same particle.

CANDIDATES:
(a) 35 — two of three printings, and the only value consistent with the
    sum-to-100 constraint that all 564 design rows satisfy (36+60+5+0 = 101);
(b) 35 in the design table and 36 verbatim in the rank table, each kept as its
    own locus prints it, which is what is on disk now;
(c) 36.

EVIDENCE NEEDED: none beyond what is here — the invariant settles it on the
merits. The decision is one of **policy**: whether a curation may correct a
published number when an unprinted design constraint decides it, or must carry
the disagreement forward. Recorded because the project's own rule is that the
pipeline does not get to decide what a silent or self-contradicting source
meant. Low stakes: it changes one formulation field on a row we do not emit.

---

### D-4
STATUS: OPEN

AMBIGUITY: The Fig S4e/f/g panels state neither the in vitro dose nor the
timepoint for the ranks they list. Methods fixes 20 ng total DNA per well and
72 h for every experiment after the dose-response test, and the corresponding
library panels (Fig S3j–l, S3m–o) are labelled "72 hours" — but library 2 also
has 4 h and 48 h datasets (Fig S3bb–ll, mm–uu), so the inference is not free.

CANDIDATES:
(a) record 20 ng / 72 h, citing Methods;
(b) leave dose and timepoint blank on these rows;
(c) record them with `structure_assumed`-style provenance text saying the
    panel itself states neither.

EVIDENCE NEEDED: nothing further is obtainable from the paper. `05_values.csv`
currently does (c) — it carries the Methods values with the sentence "NOT
printed on the Fig S4 panel itself" attached. Moot unless D-2 resolves in
favour of use.

---

### D-5
STATUS: OPEN

AMBIGUITY: Library 4 (the 144-formulation spleen screen, SI Fig S5a; 85 stable)
is a **separate experiment** from the 281, run on a different question with a
different readout set (8 splenic cell types, Figs 5–6). It is fully extracted
here (144/144 design rows) but has no per-LNP readout we can reach at all —
Fig 5B is a clustering heatmap in the main text, which we do not hold.

CANDIDATES:
(a) treat library 4 as a second `Experiment_ID` if this paper is ever merged;
(b) drop it, since it contributes no reachable value;
(c) keep the design table as reference material only.

EVIDENCE NEEDED: moot unless D-1 resolves in favour of inclusion. Flagged so
that a later merge does not silently pool 420 + 144 formulations into one
experiment — the paper never pooled them and the tracker's 281 excludes
library 4.
