# DECISIONS — `dahlmanje_30216729` (Sago CD et al., *Nano Letters* 2018, PMID 30216729)

Queued, not resolved alone. `provenance/PIPELINE_V2.md` §3: the pipeline detects
that a source is silent; it does not get to decide what the source meant.

---

### D-1 — Is cargo-abundance biodistribution admissible as a delivery endpoint?

The whole paper measures one thing: how much of an administered DNA barcode is
recoverable from sorted cells, by ddPCR (absolute counts) and by Illumina
sequencing (compositional "normalized delivery"). The cargo is a 91-nt
chemically modified single-stranded DNA barcode. It encodes nothing, is never
transcribed and is never translated.

My reading is that this is inadmissible, and that the case is not close:

- the brief's not-delivery list names **particle biodistribution** and
  **cellular uptake** explicitly;
- the paper's own title, abstract, Discussion and axis labels call the quantity
  biodistribution, accumulation and uptake — "Normalized Nanoparticle
  Biodistribution (%)", "Nanoparticle Accumulation (Au)", "% Normalized
  Nanoparticle **Uptake**";
- the `LR_2023` barcoded precedent (407 of 1,803 in-vivo rows) is narrower than
  "pooling is allowed": it requires a **translated** readout, and a screen
  sequencing cargo abundance before translation matches none of `data/`'s
  in-vivo `Experiment_method` values.

The one place this could be argued the other way is that two panels are labelled
"Absolute Delivery (AU)" (Fig 5f, Fig 6b). That is a label on the product of a
compositional share and a ddPCR count; the measured molecule is still the
barcode. I did not treat the word on the axis as changing what was measured.

STATUS: OPEN — recommendation: inadmissible, emit nothing. This is a
policy call about a readout class and it will recur: this laboratory's barcoding
papers are numerous, and a ruling here sets the precedent for all of them.

---

### D-2 — Does a one-ionizable-lipid formulation screen belong in the corpus at all?

This is the seventh one-lipid screen this run (F-067 records five being 74% of
staged rows). Here the ionizable component is `7C1` on **248 of 248**
formulations; what varies is the excipient molar ratios, the PEG alkyl tail, the
PEG molecular weight, the presence of DSPC, and the mouse genotype.

A corpus keyed on ionizable-lipid chemistry gains **one** lipid from this paper
however many rows exist — and in this case it cannot even gain that one, because
7C1 has no writable structure (D-3).

STATUS: OPEN — no action needed for this folder, since D-1 already closes it
at zero. Flagged because F-067's policy question is now answered by a seventh
instance and the pattern is not incidental: the Dahlman barcoding design holds
chemistry fixed by construction, so these papers are *structurally* one-lipid.

---

### D-3 — `7C1` has no atom-resolved structure. Confirm it must never be written.

SI p14 draws 7C1 as a **generic scheme**: an epoxide labelled `R=C13` (hence a
C15 epoxide) plus `PEI600`, 90 °C, 72 h, EtOH. PEI600 is branched polydisperse
polyethylenimine of nominal MW 600, and the degree of epoxide substitution is
never stated. The paper gives no molecular formula, no HRMS and no NMR for it.

So 7C1 is the ionizable-component analogue of a PEG-lipid, which `CLAUDE.md` §3
keeps out of `provenance/reagent_registry.py` precisely because the block length
is a distribution and no single SMILES is right.

Measured, so the claim is not a guess: exact name match for `7C1` in the corpus
returns **0** rows; `PEI600` and `PEI` return **0**. A naive substring join on
`7C1` returns 36 rows and 35 structures, and every one is an Anderson-lab
lipidoid (`117C14`, `A7C10`, `B7C10`) where the digits are tail lengths.

STATUS: OPEN — recommendation: if any future curation of a 7C1 paper is
admitted, `IL_SMILES` must be **blank** with `confidence = low` and a
`curation_note` naming the polydispersity, never a representative monomer. A
made-up monomer here would be indistinguishable from a measurement downstream.

---

### D-4 — No bond trace of the drawn cholesterol and DSPC was attempted.

SI p14 draws cholesterol and DSPC in full, and they are vector: in the structure
band, **391 filled paths carrying 1,598 segments** against only **4 stroked
paths**, so bonds here are filled polygons with `width=None` and any trace needs
`include_filled=True`.

A full connectivity trace would have produced the folder's only genuinely
independent C10 anchor. I did not do it, for a stated reason rather than a
pretended impossibility: this paper emits zero rows, so the anchor would protect
nothing, and cholesterol's eight stereocentres are not recoverable from a 2D
drawing without wedge analysis — which is the exact shape of the measured
failure where a hand-built cholesterol came out a diastereomer with 4 of 8
centres inverted, identical in formula, in exact mass to six decimals and in the
skeleton InChIKey block.

Consequence, stated plainly: **this folder has zero independent anchors.** The
two structures in `reference_structures.csv` are registry pins and are marked
`anchor_is_a_registry_pin_NOT_independent`; they are not counted as a C10
result.

STATUS: OPEN — recommendation: leave untraced unless this paper is ever
admitted, in which case trace DSPC first (one stereocentre, not eight).

---

### D-5 — Library 2: the paper says 115 pooled, the diameters allow at most 110.

The paper's two inclusion criteria are a 20–200 nm hydrodynamic diameter and a
single-inflection autocorrelation curve. Only the diameter is printed per LNP,
so the diameter gate is a **superset** of the pooled set and cannot be smaller
than it.

Library 1 behaves: **119** in band vs **111** stated pooled, a superset margin of
+8. Library 2 does not: **110** in band vs **115** stated pooled, a margin of
**−5**. Nine library-2 diameters sit outside the band (LNP 3 = 290 nm,
11 = 12.5, 23 = 15, 29 = 264, 37 = 270, 50 = 16, 68 = 17, 92 = 297, 103 = 545).

Either five printed diameters are wrong, or five out-of-band LNPs were pooled
anyway, or the stated 115 is wrong. The SI never marks which LNPs were pooled,
so it is not resolvable from the paper.

STATUS: OPEN — no value depends on it here. Recorded because it is the kind
of self-inconsistency that would silently set a pool size if this design were
ever admitted.

---

### D-6 — The one admissible readout class in the paper has no data.

Methods, "Endothelial RNA interference", describes injecting 7C1 carrying
1 mg/kg siRNA and measuring **ICAM2 protein MFI** by flow cytometry at 72 h.
Protein-silencing MFI **is** corpus-admissible.

But the only siRNA ever injected is the non-targeting control **siLuciferase**,
and its ICAM2 MFI is *defined* as 100%. There is no targeting siRNA and
therefore no knockdown. Measured: `ICAM` occurs exactly **twice** in the
165,357-character article XML, both inside that one Methods paragraph, and
**zero** times in the 25,876-character SI text.

Separately, Methods also says cells were isolated "96 (for in vivo gene editing)
hours after injection", and the paper contains no gene-editing experiment at
all — no Cas9, no sgRNA, no TIDE, no indel measurement.

STATUS: OPEN — recommendation: treat both as Methods boilerplate carried from
a sibling paper. Flagged rather than dismissed because these two sentences are
the only places in the paper where an admissible endpoint could have hidden, and
a reader grepping Methods for "gene editing" or "protein expression" would
conclude wrongly.
