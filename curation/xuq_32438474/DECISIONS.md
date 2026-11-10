# DECISIONS -- curation/xuq_32438474 (proposed `YL_2020`)

Li, Jarvis, Zhu, Glass, Ogurlu, Gao, Li, Chen, Yu, Yang, Xu.
*Protein and mRNA Delivery Enabled by Cholesteryl-Based Biodegradable Lipidoid
Nanoparticles.* Angew. Chem. Int. Ed. 2020. PMID 32438474, PMCID PMC7679290,
DOI 10.1002/anie.202004994.

`provenance/build_dataset.py` refuses to merge a folder with an `OPEN` entry.
Nothing here is for this curation to answer alone: the pipeline detects that a
source is silent, it does not get to decide what the source meant
(`PIPELINE_V2.md` s3).

**Nothing in this folder has been merged into `data/`.** The corpus is
unchanged.

---

### D-1
Protein delivery is a different cargo from mRNA, and this paper's largest
value block is protein. Figure 3B is 36 library values of **GFP mRNA**
transfection across four cell lines. Figures S24C and S24F are 27 + 27 library
values of **(-30)GFP-Cre recombinase protein** delivery and of the
Cre-mediated recombination it causes, in one cell line at three protein
concentrations. They share no cargo, no dose unit (ug/mL of mRNA against nM of
protein) and no readout definition, so they are written as three separate
blocks (`YL_2020`, `YL_2020_2`, `YL_2020_3`) in `05_values.csv`.

The question is which of the three, if any, belong in a corpus whose other 43
in-vitro experiments are overwhelmingly nucleic-acid delivery. Pooling the
protein blocks with mRNA screens would put a protein-delivery rank next to an
mRNA-delivery rank under one `Model_target`; keeping only Figure 3B discards
two thirds of the paper's quantitative content. Also note S24F is a
*downstream functional* readout of the same delivery event as S24C
(Spearman 0.946 between them over the 27 shared cells), so admitting both
would put two near-duplicate rankings of the same 9 lipids into the pool.

STATUS: OPEN

### D-2
Figure 3B's zero point. The y scale is fitted on the four detected tick marks
only, as `curation/lib/axis.calibrate` requires. That fit puts value 0 at row
**1220.78**, while the first row of the drawn axis rule is **1223** -- a
**2.22 px, 0.48 percentage-point** gap. Both of the paper's printed anchors
come out low by almost exactly that much: Lpf2k/HeLa is 97.56 against a stated
98 (-0.44 pp) and naked mRNA/HeLa is 1.46 against a stated 2 (-0.54 pp).

One systematic explains both, and shifting the origin to the axis rule would
make both land. It has not been done. A calibration moved so that the printed
statistics land is a calibration fitted onto its own validation target --
error class 27, which hid a uniform ~1.1% bias across all 48 of `LL_2026`'s
printed hit rates (F-006) -- and `axis.calibrate` deliberately offers no way
to supply an origin separately from the scale (F-033).

So the choice is between two defensible zeros, and it is a research decision:
accept a bounded -0.5 pp systematic on every Figure 3B value, or re-anchor on
the axis rule and say so in Methods. It matters most at the floor, where the
corpus takes log10 before z-scoring: 1.03 against 1.51 is 0.17 of a decade.

STATUS: OPEN

### D-3
Figure S24E gives 12 mean-fluorescence-intensity values, and neither the panel
nor its caption states the protein concentration. Panels C and F on the same
page both carry a 25/50/100 nM legend; E carries none, and there is no fourth
legend anywhere on the page. The values are extracted into
`05a_si_vector_values.csv` (0-800 a.u. axis, held-out minor ticks reproduce to
0.29 a.u.) and are deliberately **not** carried into `05_values.csv`, because a
delivery value with no dose is not comparable to anything.

Emit them with a blank dose, or leave them out?

STATUS: OPEN

### D-4
Head 87 is N,N-bis(2-hydroxyethyl)-1,3-propanediamine: it offers the acrylate
a primary amine **and two primary alcohols**. Aza-Michael and oxa-Michael
addition give the identical molecular formula, so the SI's own [M+H]+
(C77H133N4O10S4) cannot distinguish 2 tails on the amine from 1 on the amine
and 1 on an alcohol. The paper draws no structure for 87-OCholB.

Both tails are placed on the amine -- the stated chemistry ("Michael addition
reaction", amine heads plus acrylate tails, 70 C, 48 h, no catalyst) and the
standard lipidoid route. The row carries `structure_assumed = 1` and
`confidence = medium` with that reason in `curation_note`. 87-OCholB is one of
the two members the paper itself calls consistently ineffective, so its value
is near the floor either way, but the structure is a genuine assumption.

Confirm the assignment, or exclude 87-OCholB?

STATUS: OPEN

### D-5
Nine of the 18 distinct structures carry **no delivery value in this paper**:
75/76/77-O16B and 75/76/77-OChol appear only in the Figure 1C cytotoxicity
panel, and 75/76/77-OCholB-1 only in Figure S23's Cas9 panels (which report a
knockout that did not happen rather than a delivery efficiency). They are
enumerated in `03_enumerated.csv` and formula-confirmed 9/9, but they would
enter the corpus as structures with no `Experiment_value`.

Include them as value-less library members, or restrict the merge to the nine
OCholB lipidoids?

STATUS: OPEN

### D-6
`Experiment_ID`. The corpus convention is `INITIALS_YEAR`, and the five most
recent curations all use the **first** author (`ZT_2026` Tian, `LL_2026` Lin,
`WT_2026` Toh, `XD_2026` Duan, `QL_2026` Liu) while this folder is named after
the **corresponding** author (`xuq_32438474`, Qiaobing Xu). First author here
is Yamin Li, giving `YL_2020`; corresponding gives `QX_2020`. Neither collides
with the 43 vitro / 20 vivo ids in the live corpus. `YL_2020` is used
throughout this folder as a placeholder.

STATUS: OPEN

### D-7
**No in-vivo value exists in this paper.** Figures 5, 6, 7 and Supplementary
Figures S27-S45 report every in-vivo result as the presence or absence of
tdTomato signal in a micrograph or an IVIS image. There is no axis, no
intensity scale and no table; Figure 7B is a GFAP/NeuN colocalisation
percentage of tdTomato+ cells, which is a cell-type attribution and not a
delivery measurement. Three facts are recoverable and are recorded in
`07_qc.json` rather than invented as numbers: 80-OCholB gave lung signal,
76-OCholB gave spleen signal (protein and Cre mRNA both), and no construct gave
liver signal.

If the corpus wants a presence/absence in-vivo row for those, it needs a
convention it does not currently have.

STATUS: OPEN
