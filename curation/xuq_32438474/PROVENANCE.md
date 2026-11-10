# PROVENANCE -- curation/xuq_32438474

**Paper.** Li Y, Jarvis R, Zhu K, Glass Z, Ogurlu R, Gao P, Li P, Chen J, Yu Y,
Yang Y, Xu Q. *Protein and mRNA Delivery Enabled by Cholesteryl-Based
Biodegradable Lipidoid Nanoparticles.* Angew Chem Int Ed 2020.
PMID 32438474 - PMCID PMC7679290 - DOI 10.1002/anie.202004994.
Proposed `Experiment_ID` **YL_2020** (see DECISIONS D-6).

**Nothing here has been merged into `data/`.** The corpus is unchanged at
29,195 x 79 vitro / 1,803 x 79 vivo.

## Run order

```bash
python3 -I stage01_blocks.py       # 01_blocks.csv, 01_tails.json
python3 -I stage03_enumerate.py    # 03_enumerated.csv   (18/18 formulas agree)
python3 -I make_reference.py       # reference_structures.csv  (C10 anchors)
python3 -I stage05a_si_vector.py   # 05a_*  Figures S24C/E/F, PDF vector
python3 -I stage05b_fig3b.py       # 05b_*  Figure 3B, raster
python3 -I stage05_values.py       # 05_values.csv
python3 -I 07_seed_checks.py       # 07_seed_checks.json  (structure seeds)
python3 -I stage07_qc.py           # 07_qc.json           (everything + seeds)
```

## Where the values are

| panel | what | n | tier |
|---|---|---|---|
| Figure 3B | GFP **mRNA**, % GFP+ cells, 12 samples x 4 cell lines | 48 (36 library) | 5_pixels |
| Figure S24C | (-30)GFP-Cre **protein**, % GFP+, 12 x 3 doses | 36 (27 library) | 4b vector |
| Figure S24F | recombination, % DsRed+, 12 x 3 doses | 36 (27 library) | 4b vector |
| Figure S24E | mean fluorescence, 12 bars, **no stated dose** | 12 | not emitted |
| in vivo | nothing numeric | **0** | - |

`00_value_sources.json` records why tiers 1-4 fail for Figure 3B, with a
sentence per tier that a reviewer can check.

## The library

**18 distinct ionizable structures**, not 18 formulations: 9 OCholB (heads 75,
76, 77, 78, 80, 81, 87, 90, 304), 3 O16B and 3 OChol (heads 75-77), 3
single-tailed OCholB-1 (heads 75-77). `75-OCholB-F` and friends are the SAME
structures formulated 1/1 w/w with DOPE and are counted as formulations. The
paper states no library size anywhere; it was established from the SI's two
[M+H]+ tables (Figures S13 and S16) and cross-checked against the Results'
own "six out of the nine examined OCholB LNPs".

## Claims this folder makes that `check_docs.py` could gate

- `03_enumerated.csv` has 18 rows and 18 distinct full InChIKeys.
- All 18 agree with the SI's [M+H]+ formula tables.
- `05_values.csv` has 120 rows, 90 of them on library members.
- `emit.caveat_count` over the library rows is 20 (the 77- and 87- series).
- Corpus overlap is 0 by full InChIKey, 0 by skeleton InChIKey, 0 by `IL_name`.
