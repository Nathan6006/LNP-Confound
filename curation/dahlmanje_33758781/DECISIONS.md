# DECISIONS — `dahlmanje_33758781` (Gan et al. 2020, prospective `ZG_2020`)

Queued, not resolved alone. `provenance/build_dataset.py` refuses a folder with
an `OPEN` entry, which is correct here: **this folder emits zero rows and should
not be added to `MERGE`.**

---

### D-1 Zero rows, for two independent reasons

The paper measures functional mRNA delivery for 109 LNPs across 16 cell types
in vivo, and none of it is emittable.

**Reason 1 — one ionizable lipid.** Every one of the 109 LNPs carries
`cKK-E15`. The library axes are the helper phospholipid (10), the sterol (3)
and the molar ratio (4). A corpus keyed on ionizable-lipid chemistry gains
exactly **one** structure from this paper no matter how many rows are written,
and any rows written would differ only in `HL_*`, `CHL_*` and `*_molratio`.

**Reason 2 — nothing is attributable.** The design supports 109 × 16 = **1,744**
per-LNP values and the authors computed them; the paper attributes **0**. Fig 3c
plots ~109 unlabelled points; Fig 3f/3g collapse to 10 and 3 group means over
~11 LNPs each; Fig 3h is a derived enrichment statistic. The 11 of 120
formulations dropped at QC are never identified, so the 109 members cannot even
be enumerated.

Either reason alone is sufficient. The closest call is Fig 3i (LNP1, injected
alone, individually measured in 3 liver cell types) — blocked by X-6: the paper
never states LNP1's sterol or molar ratio.

STATUS: OPEN

---

### D-2 `cKK-E15`'s tail length is an inference, not a reading

The paper names `cKK-E15` twice, draws it once (Fig 2a), cites **no** source for
it, and abbreviates each of its four tails to a `CnH2n+1` label whose subscripts
measure 2–3 px in the only figure we can obtain (709 px web JPEG; the
publisher's full-resolution figure sits behind a browser this task does not
use). The **topology** is readable and matches the corpus `cKK-E12` exactly — a
2,5-diketopiperazine, two `(CH2)4` linkers to two tertiary nitrogens, each
bearing two `CH2-CH(OH)-R` arms. The integer *n* is not.

The inference `E12 → 12-carbon arms` is **measured from the corpus**, not
remembered: `AC_2025`'s `CKK-E12` has four arms of exactly 12 carbons
(`stage02_corpus_overlap.py` walks them). Applying that convention gives
`cKK-E15` 15-carbon arms.

**Why this matters more than it looks.** `cKK-E12` is already in the in-vivo
pool (`BL_2023`, as `MD-1`, 6 rows). The in-vivo pool has only **two** lipids
shared across experiments. So if the convention inference were wrong and
`cKK-E15` were in fact `cKK-E12`, this paper would supply a **third** shared
in-vivo structure. Under the inference it supplies none: 0 rows match the
inferred full InChIKey `FGDZRFPFIDLJAJ-UHFFFAOYSA-N`, and 0 match even its
skeleton block.

`03_enumerated.csv` therefore stores `cKK-E15` with an **empty** SMILES,
`confidence = low`, `structure_assumed = 1`, and `reference_structures.csv`
stores an empty `reference_SMILES` so the gap is visible. It is not written as
a structure anywhere.

STATUS: OPEN

---

### D-3 The regiochemistry of the ten phospholipids

Eight product names put the adamantylacetyl on the primary (C3) hydroxyl;
**A-10 and A-11 print it on C2**, with the fatty acyl on C3 (X-1). We emit all
ten with adamantylacetyl on C3, from the synthesis route and the uniform NMR
pattern. The two regioisomers are identical in formula and in exact mass to
**0.000000 Da**, so no mass check can adjudicate; seed `S5` confirms the HRMS
gate runs green on the wrong one.

This affects no emitted value (there are none) but it does affect `A-11`, the
paper's lead, and it would affect `HL_SMILES` if this folder were ever revived.
Flagging rather than silently preferring the majority spelling.

STATUS: OPEN

---

### D-4 Would these rows be admissible *in kind*, if they were attributable?

Worth settling for the next barcoded Dahlman-lab paper, since four more sit in
`curation/`.

The readout **does require expression**: barcodes are sequenced out of
tdTomato+ FACS-sorted cells, and a cell is tdTomato+ only once Cre mRNA has been
translated into functional Cre protein that excises a Lox-Stop-Lox cassette.
That clears the F-066 bar, which a screen sequencing cargo abundance *before*
translation fails. `Cre→tdTomato` recombination is explicitly admissible.

So the answer is **yes in kind, no in fact**. Pooling is not the obstacle —
`LR_2023` contributes 407 of 1,803 in-vivo rows from a barcoded pool. The
obstacles are D-1's two reasons.

One caveat for whoever settles this: Fig 3e is `QUANT`/ddPCR **biodistribution**
(cargo abundance, translation-independent) and is *not* admissible, even though
it sits in the same figure as the admissible panels and is measured on the same
sorted cells.

STATUS: OPEN

---

### D-5 `Experiment_ID` would be `ZG_2020`, not `JD_2020`

Europe PMC's `core` record gives the author string as
`Gan Z, Lokugamage MP, Hatit MZC, Loughrey D, Paunovska K, Sato M, Cristian A,
Dahlman JE`. **Gan Z is first author**; Dahlman is last and corresponding, and
is what the folder is named after. The project's convention keys on the first
author, so the id would be `ZG_2020`. No `ZG_2020` exists in either corpus file
(43 vitro + 20 vivo ids enumerated 2026-11-10). Recorded because a sibling in
this run had the same folder-vs-first-author mismatch.

STATUS: OPEN
