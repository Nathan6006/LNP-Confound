# ng_2025 — Gong et al., *Nat Biomed Eng* 2025, PMID 40681859

Mannich-reaction combinatorial libraries; lead ionizable lipid **C-a16**.
DOI 10.1038/s41551-025-01422-8. Six libraries, **310 crude lipids screened in
vitro**, 11 taken into mice.

**Status: extracted and checked, NOT merged.** `provenance/build_dataset.py`
does not reference this folder.

## What the paper gives us, and why this one was cheap

The publisher's Source Data workbook (`41551_2025_1422_MOESM3_ESM.xlsx`, sheet
`Fig1`) holds **the numbers behind the Fig. 1 heat maps**, so not one in-vitro
value is digitised from pixels. The structures are the work; the values are a
spreadsheet read.

| output | contents |
|---|---|
`03_enumerated.csv` | the 310 lipids, name / library / tail / amine / SMILES |
`04_values.json` | 310 in-vitro log10 luminescence + both in-vivo panels with replicates |
`07_NG_2025_rows.csv` | **333 rows** in the 72-column corpus schema |
`07_flags.csv` | per-row confidence flags (116 of 310 lipids) |
`07_held_back.csv` | rows deliberately not emitted, with the reason |
`09_qc.json` | every check below |
`reference_structures.csv` | six C10 anchors |

Rows: `NG_2025` 310 in-vitro (MC38); `NG_2025_2` 23 in-vivo in two endpoint
blocks — liver after IV (Fig. 1m, 12 rows incl. DLin-MC3-DMA) and muscle after
IM (Fig. 1p, 11 rows). Blocking by `(Experiment_ID, endpoint)` and z-scoring
each block separately follows `BL_2023` and `ZH_2023`, and avoids inventing two
"sources" out of one lab for the source-classification analysis.

## Structures: read by machine wherever the page allowed it

**Fig. 1b and 1c are vector line art, not raster.** `stage01_vector.py` reads
the page's own bond coordinates, so the four phenolic tails and all sixteen
amine cores come out of the PDF as graphs rather than as something counted by
eye — this project's documented ±1 failure mode (ERROR_CLASSES row 19). Atom
identity falls out of each label's path signature, which is exactly separable
on this page (`N` = 10 line items and no curves; `H2N`/`NH2` = 26 + 18; the
stacked `N/H` = 22 + 0; `O` = curves only; `HO` = 12 + 20).
`stage02_blocks.py` compares those twenty machine-read structures against the
twenty written by hand: **20/20 identical as canonical SMILES.**

The Supplementary figures are raster at ~106 dpi effective, so tails E–R and
amines b1–b13 were read from the images. Chain lengths there still are not
eye-counted: `stage03_zigzag.py` takes the ink centroid per pixel column and
counts turning points. On the thiol tails it agrees with an independent
skeleton trace (13/11/9/7 vertices → C12/C10/C8/C6).

**Degree of substitution is not one rule** (`stage04_enumerate.py`):

* phenol / thiol / pyrrole / acetophenone — two tail arms per **primary** NH2;
  secondary N-H, tertiary N and the a16 amide N-H do not react.
* dialkyl ketone **1:2:1** — a double Mannich that closes a ring:
  1-R-4-oxopiperidine-3,5-dicarboxylate.
* dialkyl ketone **1:4:2** — the same again: the bicyclic bispidinone, two
  amine residues per lipid.

The ketone products are therefore *not* `R-N(CH2-tail)2`, and the two printed
ratios are exactly the one-ring and two-ring closures.

## Checks (`09_qc.json`)

* **C1 mass.** The paper characterises one compound. Enumerated C-a16 is
  **C72H112N4O6, MW 1129.71** against the **1129.71** printed in the mass
  spectrum's own inset (Supplementary Fig 19a), with [M+H]+ 1130.7 vs the
  observed 1130.4 and [M+2H]2+ 565.4 vs 565.8. Its 1H NMR integrals sum to 106
  H against 112 in the structure: the missing 6 are the 4 phenol OH and 2 amide
  NH, which is itself a structural statement.
* **C2 the paper's own hit rule, within the phenol library.** "Hits" are lipids
  above 5× DLin-MC3-DMA. One threshold (log10 = 5.938) reproduces **all 11**
  lipids the paper then took into mice, exactly — no extras, none missing —
  *among the 64 phenol lipids*, which is the scope the paper itself states the
  rule in. Applied to all 310 the same threshold also catches **G-a12** (thiol,
  6.074), which clears the stated criterion yet was never taken into mice, so
  no single threshold over the whole library reproduces exactly those 11. The
  phenol window is also narrow: (5.937736, 6.036851].
* **C3** 310 enumerated = 310 stated. **C4** 310 unique structures, no
  duplicates. **C5** all parse, all contain nitrogen, no radicals.
* **C6** molar ratios sum to 100 for every row.
* **C7 heat map vs workbook.** Every one of the 310 drawn cells is matched to
  its workbook number through the figure's own axis labels; Spearman between
  drawn colour and stored value is **≥ 0.998 in all six panels** — the block
  order, the row order and the panel assignment are all confirmed, not assumed.
* **C8** zero of the 310 novel lipids already exist in the corpus; the one
  benchmark row matches the corpus DLin-MC3-DMA string exactly.
* **Chain lengths confirmed by the paper's own axis.** Supplementary Fig 17c
  plots hits against "number of carbon in the alkyl chain" with ticks at
  **0 / 6 / 8 / 9** — exactly the phenol tail lengths read geometrically for
  A / B / C / D, from a part of the paper that was not used to derive them.
* **C9** each z-score block is mean 0 / sd 1 (population sd).
* **C10** six drawn-structure anchors, one per chemistry plus both ketone ring
  systems, all **exact matches**. The strongest is **B-a1**, whose product the
  paper draws in Fig. 1d — vector, so its reference string was machine-read
  from bond coordinates rather than typed by someone who already knew the
  answer.

## Three things the paper contradicts itself about

1. **Amines b3 and b4** (28 lipids, flagged `low`). The four amine-core panels
   draw one CH2 between the two nitrogens — aminals, which are not isolable
   reagents. The pyrrole and acetophenone product panels draw the same two
   amines with a direct **N–N** bond: 1,1-dimethylhydrazine and
   1-amino-4-methylpiperazine, both catalogue chemicals. Emitted as the
   amine-core panels draw them, since those four panels agree and are the
   paper's statement of what was bought; the alternative SMILES is in
   `07_flags.csv`, one line to reverse.
2. **The ketone products** (104 lipids, flagged `low`). Enumerated from the
   schemes. In the 1:2:1 library **15 of 52** are drawn as the open-chain
   mono-adduct instead — I 6 (b4, b5, b7, b9, b10, b11), J 1 (b5), K 4 (b1, b2,
   b7, b12), L 4 (b1, b2, b4, b10) — so the same amine is drawn both ways under
   different tails. The 1:4:2 panels are mixed too, roughly 18 of 52, making
   about **a third of all 104** ketone lipids drawn against their own scheme.
   The decisive detail is that **every open-chain drawing gives the Mannich
   nitrogen an extra methyl it cannot have** (K-a-b1's amine is methylamine,
   yet it is drawn as N(CH3)2 plus one glutarate arm): they are template
   artefacts, not an alternative product, which is why all 104 are enumerated
   as rings.
3. **The Fig. 1p benchmark (held back).** The x-axis labels the twelfth bar
   DLin-MC3-DMA; the Source Data column for those same three numbers is headed
   SM-102. Different molecules, so the row is not emitted — a wrong benchmark
   identity would land measured values on the wrong structure in exactly the
   shared-lipid comparison this corpus exists for.

Four smaller ones, recorded but not blocking:

* The phenol hit counts are printed as **0/3/6/2** for tails A–D in the Results
  text **and again, graphically, as the bars of Supplementary Fig 17c**, where
  the paper's own Source Data and its own in-vivo panel give **0/2/6/3** — B
  and D transposed, in two independent places, not a sentence-level typo.
* "~15-fold higher luminescence" for C-a16 over DLin-MC3-DMA cannot hold at the
  same time as the hit rule: 15× implies MC3 ≈ 3.5e5, at which only **four**
  lipids clear 5× rather than eleven. Reproducing the eleven needs MC3 in
  (1.7e5, 2.2e5], i.e. a fold-change of 24–30×. Most likely the two numbers
  come from different runs (Supplementary Fig 17f is a separate n = 3
  experiment), but the paper does not say so.
* The main text puts the IM readout at 6 h; Supplementary Fig 18 says 24 h.
* The Results call the screening mRNA 5-methoxyuridine-modified; Methods say
  all mRNAs were 1-methylpseudouridine-modified.

And one value disagreement inside the paper, flagged in `07_flags.csv`:
**O-b10** is the highest value in the whole pyrrole panel (4.07e4) in the
Fig. 1 Source Data, but Supplementary Fig 17b, which re-measures six library
members individually, plots it ~24× lower; the other five agree with Fig. 1 to
within 1.9×.

## Not extracted / left blank, on purpose

* `IL_to_nucleicacid_massratio` — Methods give "lipid and mRNA at a mass ratio
  of 20:1" without saying whether that is total lipid or ionizable lipid, and
  the column means the latter.
* `Dialysis_buffer` for the in-vitro screen — the dialysis sentence follows the
  microfluidic sentence, and the screen was pipette-mixed; it is not stated
  whether screening LNPs were dialysed. Filled for the in-vivo rows (PBS).
* Zeta potential and apparent pKa — never reported.
  **Size, PDI and encapsulation efficiency ARE printed**, as a table in
  Supplementary Fig 17d covering all 11 in-vivo hits plus DLin-MC3-DMA (size
  75.9–89.6 nm, PDI 0.12–0.189, EE 82.4–96.5 %). They are not emitted only
  because the 72-column schema has no columns for them; the numbers are in this
  README and in `09_qc.json` so they are not lost.
* `Dose_ug_nucleicacid` for the in-vivo rows is **derived**: 0.125 mg kg-1 ×
  the ~20 g body weight Methods state = 2.5 µg.

**Two molar ratios exist in this paper.** Every row here carries
50 / 10 / 38.5 / 1.5 (IL / DOPE / cholesterol / DMG-PEG), which is what Methods
state for the screening and in-vivo-screening LNPs. Fig. 2a separately gives
the *optimised* C-a16 LNP as **35 : 16 : 46.5 : 2.5**, size 80.6 ± 9.6 nm, PDI
0.095 — that formulation belongs to the downstream C-a16 experiments (Figs 2-5)
and none of those are emitted here. The Methods sentence is general, so the
in-vivo screen of the 11 hits is taken to be at 50/10/38.5/1.5; the paper does
not restate the ratio for that panel, which is the one soft spot in the
formulation metadata.

## Reproduce

```bash
cd curation/ng_2025
python3 -I stage01_vector.py s41551-025-01422-8.pdf 2 85 103 297 300 3.6 > 01_fig1bc_graph.json
python3 -I stage02_blocks.py 01_fig1bc_graph.json T_A,T_B,T_C,T_D,H_a1,H_a2,H_a3,H_a4,H_a5,H_a8,H_a9,H_a6,H_a7,H_a11,H_a10,H_a12,H_a13,H_a14,H_a15,H_a16
python3 -I stage04_enumerate.py && python3 -I stage05_values.py
python3 -I stage06_heatmap.py && python3 -I stage07_emit.py
python3 -I stage08_reference.py && python3 -I stage09_qc.py && python3 -I stage10_flags.py
```

Reviewed adversarially on 2026-10-30 (Tier 2, read-only). It re-derived every
structural claim independently and found **no emitted structure and no emitted
value wrong**. Its nine findings were all verified here: seven were real and are
fixed above (the B-a1 anchor was a typed literal while the note claimed it was
machine-read — it is now derived at run time from `01_fig1d_graph.json`; the
ketone open-chain count; the unscoped hit-rule claim; the Supplementary Fig 17d
table; the 15-fold and Supp Fig 17c contradictions; the benchmark
`Formulation_ID`; the O-b10 disagreement). One correction to the reviewer: the
vector route does work on Fig. 1d — it needs the glyph merge radius at 3.0 pt
rather than 3.6, which is why that radius is a script argument.

`NG_2025_blocks.csv`, `NG_2025_conditions.csv` and `NG_2025_checks.csv` are a
prior extraction done in Claude chat, kept for comparison. Its 20 Fig. 1b/1c
structures are confirmed 20/20; its b3/b4 reading is the one flagged above; and
its stated degree-of-substitution rule ("2 arms per primary NH2") is wrong for
the 104 dialkyl-ketone lipids, which are ring systems.
