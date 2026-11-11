# liu_42585110 -- provenance

**Status: EXTRACTED 2026-11-10. NOT MERGED.** `DECISIONS.md` carries two OPEN
entries and `build_dataset.py` refuses a folder with one, which is the intended
behaviour: D-1 (an All-Rights-Reserved source against CLAUDE.md section 8) and
D-2 (what Figure 1a's unitless scale is) are questions this curation is not
entitled to answer.

| | |
|---|---|
| Paper | Jiang K, Shao B, Li R, Liang T, You C, Wang L, Cen J, Tan J, Zhang Y, Chen M, Liu S. *Retinoic-Acid-Derived Ionizable Lipids Enable Spleen-Selective Antigen Expression for mRNA Cancer Vaccination.* Angew. Chem. Int. Ed. 2026; 65:e3823572 |
| DOI | 10.1002/anie.3823572 |
| PMID | **42585110** (indexed; Europe PMC core record retrieved 2026-11-10) |
| PMCID | none. `inEPMC: "N"`, `inPMC: "N"`, `hasSuppl: "N"` |
| Experiment_ID | **KJ_2026** -- from the FIRST author, Jiang Kai. The folder is named for the PI (Liu) and the folder name is **not** the key. Checked against all 54 distinct corpus `Experiment_ID`s **by PMID, not by the year token**: no `KJ_*` exists and no corpus row carries PMID 42585110, so there is no collision of the F-075 kind |
| Licence | **All Rights Reserved, Wiley-VCH GmbH. NOT Creative Commons.** See `00_acquisition.json` and DECISIONS.md D-1 |
| Rows emitted | **123** = 41 in vitro (DC2.4) + 82 in vivo (41 formulations x liver and spleen) |
| Distinct ionizable lipids | **40** novel RAILs + the Dlin-MC3 benchmark = **41 distinct structures**, all 41 by full InChIKey |
| Value source | all tier **5_pixels**; see `00_value_sources.json` |

## The endpoint

Firefly-luciferase reporter expression, twice:

* **in vitro** -- luminescence in DC2.4 cells after 6 h, 0.5 ug/mL mRNA
  (Figure 1a);
* **in vivo** -- whole-organ Avg Radiance in liver and spleen 6 h after IV
  administration of 0.1 mg/kg FLuc mRNA to BALB/c mice, ex vivo IVIS,
  normalised to organ area (Figure 1c).

Both require translation of the delivered mRNA, which is the corpus's delivery
endpoint. Nothing pooled, nothing barcoded, so the compositional-readout
question does not arise.

## Why every value is pixels

`00_value_sources.json` holds the checkable sentences. In short, and each of
these is a measurement rather than an assertion:

* **No Source Data file exists.** Wiley publishes one supplementary asset for
  this article, the SI .docx.
* **The SI contains no tables at all.** `<w:tbl>` element count = 0 AND the
  string "Table S" occurs **0** times in 959,060 characters of
  `word/document.xml`. The briefing warned that a `<w:tbl>` count of 0 can mean
  the tables are EMF images -- a sibling curation in this run recovered three
  that way -- so all **83** files in `word/media` were enumerated and typed:
  48 PNG + 28 SVG (Schemes S1-S12 and Figures S1-S29: synthesis schemes and
  NMR/ESI-MS spectra), 4 TIFF (Figures S30-S33), 2 EMF (Figures S34-S35) and
  the Figure S37 pair. **The two EMF metafiles were opened**; they are bar
  charts of cell viability and endosomal escape. Neither is a table and
  neither is a delivery endpoint. There were no hidden tables.
* **No repository.** Data Availability, verbatim: "The data that support the
  findings of this study are available from the corresponding author upon
  reasonable request." `find_repo.py` returns `repo_hits: []` for the first
  author AND for the PI. An accession sweep over **3,126,538 characters**
  (40,301 main-PDF text + 959,060 SI document.xml + 2,127,177 SI SVG text
  layers) returned 5 apparent hits -- `GSE1`, `SRP3`, `SRP7`, `SRX1`, `IPX6` --
  **every one of which was checked by printing its surrounding context and is
  base64-encoded embedded-font data inside an SVG.** True count of real
  accessions: **0**.
* **No tier-4b vector geometry.** `page.get_drawings()` returns 187 path
  groups on page 1 (the masthead) and 1-3 on every other page; the five figure
  pages carry 3, 2, 2, 2, 2. Those are the Wiley download watermark and the
  page rule.
* **Native dpi measured, not guessed:** Figure 1 is embedded at 2051 x 2352 px
  on a 492.24 x 564.27 pt rectangle = **299.97 x 300.13 dpi**. The raster is
  EXTRACTED, never re-rendered, so there is no second rasterisation.

## The 40 lipids, and the factorial

`stage02_blocks.py` proves the library size instead of restating it:
**4 ATX tail/linker motifs x 2 linker architectures x 5 retinoids = 40**,
and refuses if that product is not 40. The paper words the same factorial the
other way round -- "two molecular scaffolds (ATX-c-RA and ATX-RA), four linker
motifs, and five retinoic-acid-derived hydrophobic domains" -- and what it
calls a linker motif is our tail motif. See `08_contradictions.csv` row 4;
AT4 in particular is not a linker-length variant of AT1 at all, it is a
reversed ester on a branched acid.

Structures are built in `stage03_enumerate.py` through `curation/lib/assemble`
only: `esterify` once and `attach` twice per lipid, with the retinoid put on
the linker BEFORE the head arm so that `esterify` always sees exactly one free
hydroxyl. **No RWMol, no RemoveAtom, no string concatenation anywhere.**

## The paper's own cross-checks, measured by script

Full numbers in `07_qc.json`. Passing every printed statistic is necessary and
not sufficient -- `QL_2026` reproduced every printed hit rate on bars up to 28%
too high -- so the load-bearing checks here are the ones that are not
recomputed marginals.

| check | result |
|---|---|
| **V1** Figure 1d is a redundant encoding of Figure 1c; read independently and compared | Spearman **0.993**, Pearson 0.996, slope **1.015**, median relative difference 4.7% over 26 comparable formulations |
| **V2** the zero of each bar axis, PREDICTED from a fit that excluded it | Figure 1c: drawn baseline 1507.361 px, predicted **1507.734** -> error **0.373 px = 0.0034 units**. Figure 1d: error **-1.004 px = 0.018 units** |
| **V3** the white point of the Figure 1a diverging colormap | **+0.0094**, i.e. **0.47%** of the half-range; max LUT colour distance **3.46** against the 5.0 gate |
| **V4** printed in-vivo statistics | AT4-13CRA spleen/liver **4.193** against the printed "approximately 4.1"; **rank 1 of 41**; splenic share **80.74%** against the paper's 78.9/19.2 renormalised to **80.43%** |
| **V5** printed in-vitro prose | the three lipids the paper names are our top **3 of 41 in its own order**; ATX-RA beats ATX-c-RA by **0.337** at Mann-Whitney **p = 0.0045**, **15/20** matched-pair wins |
| **V6** every structure against the SI's own ESI-MS | **39/39** within 1.5 mDa, worst **0.50 mDa** -- which IS the 3-decimal rounding floor. The 40th, AT3-Bex, has no printed ESI-MS and is the only one |
| **V7** the head-arm carbon count | 4-carbon reading misses by at most **0.50 mDa**; 3-carbon reading misses by at least **14,015.32 mDa** on its best case, a factor of **28,031** |
| **V8** the formulation | derived from the SI's own rule; DSPC:chol:PEG ratio preserved to **0.0**; the control's four percentages are Onpattro's published 50 : 10 : 38.5 : 1.5 to **0.0 mol%** |
| **V9** Figure 2c's PRINTED NUMERALS | mean **+1.77 pp**, worst **5.31 pp**, ordering of all five preserved |

### Round numbers the authors drew or used and never printed

These are external to this pipeline in a way a recomputed marginal is not:

1. **The zero of both bar axes.** Each axis is calibrated on its ticks for
   3/2/1 (panel c) and 6/4/2 (panel d) ONLY; the 0 tick is fused with the drawn
   baseline and excluded, and `axis.calibrate` is passed
   `forbid_target=[baseline]` so that fitting onto the validation target is a
   **TypeError**, not a silent class-27 absorption. The fit then predicts where
   0 falls: **0.373 px** and **-1.004 px** from the rule the authors drew.
2. **The white point of the Figure 1a colormap.** A diverging red-white-blue
   map has its white point at the value its middle label names. The figure
   draws no tick for it. Measured **+0.0094**.
3. **The Onpattro ratio.** The SI lists four control percentages and never says
   what they are. They are Onpattro's published composition, 50 : 10 : 38.5 :
   1.5, to **|delta| 0.0 mol%**.

### Measured sensitivity limit, stated honestly

Our Figure 1c splenic shares sit **+1.77 pp** on average and **5.31 pp** at
worst (AT4-Bex) from Figure 2c's printed numerals. The tall bars are read to
~1 px on a 110.18 px-per-unit axis, i.e. 0.009e8, so that residual is **not**
read precision on the tall bars -- it lives in the short bars, where a 1 px
error is a large fraction of the value and both organs enter a ratio.

Figure 1d exceeds the Figure 1c ratio in 18 of 26 comparisons, which could
mean 1d is the mean of per-animal ratios while 1c gives per-organ means --
those are not algebraically equal. **Tested and NOT established**: sign test
p = 0.122, median ratio 1.015. Reported as unsettled rather than asserted.

## Corpus overlap

`02_corpus_overlap.json`. **Every query carries a positive control that asserts
a specific named molecule known to be present**, because a non-zero count is not
evidence that a control works (F-070). Controls: `DLin-MC3-DMA` by InChIKey
**3,098** rows, cholesterol **30,397**, DSPC **7,038**, plus a round-trip on our
own lead.

* **ZERO** of the 40 RAILs exists anywhere in the corpus's 19,919 structures.
  A measured zero, on a join that demonstrably finds MC3.
* A SUPERSET substructure sweep, labelled as such: **0** corpus structures
  carry the retinoid trimethylcyclohexene-polyene motif and **0** carry the
  bexarotene tetrahydronaphthalene -- each SMARTS asserted to match one of our
  own molecules first. My first version of the retinoid SMARTS returned a clean
  corpus 0 **while not matching AT1-ATRA**, a molecule built three files
  earlier; the control is what caught it. That is F-070 reproduced from scratch.
* **Token-pattern test.** The paper's own spelling `Dlin-MC3` matches **0**
  corpus rows by name while the molecule is present **3,098** times under two
  other spellings. Every join here is on full InChIKey.
* **No arm of this paper is already merged**: 0 corpus rows carry PMID
  42585110 and no `KJ_*` experiment exists.

### The in-vivo contribution is the load-bearing part

The in-vivo pool has exactly **two** structures shared across experiments,
`DLin-MC3-DMA` and `SM-102`. This paper measures MC3 in **liver and spleen**:

* **Liver** -- MC3 currently has three in-vivo liver sources (`NG_2025_2`,
  `ZC_2023`, `ZH_2023`). KJ_2026 would be the **fourth**, adding 3 new
  organ-matched cross-source pairs on the one benchmark the in-vivo pool
  actually shares.
* **Spleen** -- the **first MC3 spleen measurement anywhere** in the in-vivo
  corpus, and only the **second spleen source of any kind**: all 162 existing
  spleen rows come from `ZH_2023` alone. `ZH_2023` carries MC3 only in liver,
  so the two spleen sources would still share no lipid; the gain is organ
  coverage, not a new shared-lipid pair.

### The benchmark does not top the screen

`Dlin-MC3` ranks **32 of 41** in vitro and has the **lowest** spleen-to-liver
ratio of all 41 in vivo. This is not CLAUDE.md section 8's hazard of a
universal positive control topping a screen for reasons unrelated to
transferable SAR.

## Contradictions

`08_contradictions.csv`, six rows. The one that changes structures: the SI text
and Materials both say **4-amino-1-butanol** while Schemes S4, S5 and S9 all
**draw three carbons**. Mass settles it at a factor of 28,031 and four carbons
is what is built -- the same class as `XD_2026`, where the paper drew ALC-0315
with a 3-carbon hydroxyalkyl head where the compound has 4.

## C10 anchors

`make_reference.py` **builds** two anchors; `reference_structures.csv` stores
what it built. Per the run's hard rules this folder did **not** run
`provenance/structure_check.py`, `check_anchor_independence.py` or
`check_anchor_coverage.py` -- they write tracked files other agents read -- so
the posture is recorded here by hand and the comparison is made inside
`make_reference.py` itself.

**Two anchors, of different chain lengths and different linker architectures**,
because one length cannot catch an off-by-one that scales with n:

| anchor | drawn in | result |
|---|---|---|
| `AT1-c-ATRA` | SI Scheme S5, row 2 | **exact match**, full InChIKey `ALAWTIHWVGEAKA-VDLJIEFOSA-N` |
| `AT4-Bex` | SI Scheme S9, bottom row | **exact match**, full InChIKey `CZMRBDVYRKKGAS-UHFFFAOYSA-N` |

The generator shares as little as possible with the enumerator: it imports
nothing but `rdkit`, `csv` and `ast`; it builds by writing one SMILES string
per molecule rather than by fragment assembly, so there is no `RWMol`
anywhere; and it writes each molecule in the **opposite traversal**, starting
at the retinoid's ring and finishing at the lipid tail, where the enumerator
starts at the retinoid acid and finishes at the head arm.

**Independence is executed, not asserted.** Two guards run on every build:
an **AST guard** that walks this module's own `Import`/`ImportFrom`/`Call`/
`Attribute` nodes only -- so the ban list, which is itself a set of strings in
that file, cannot match itself -- and a **template-provenance guard** that
refuses if any string literal of >=25 characters is a substring of an
enumerated SMILES. Byte-identity of the emitted anchor is **not** the signal
and testing for it is a trap (F-037): a correct re-derivation is byte-identical
to a correct enumeration once canonicalised, and both of these are.

**The anchor earned itself on its first run.** `AT1-c-ATRA` came back
`exact=False` with a different InChIKey **skeleton** block -- C54H93NO7 against
the enumeration's C55H95NO7. The error was in the anchor, not the enumeration:
the hand-written T1 arm had five carbons from the nitrogen to the carbonyl
where 6-bromohexanoate gives six, and the SI's printed [M+H]+ of 882.719 says
the enumeration is the correct half. An anchor that was a copy could not have
disagreed.

## Seeds -- every check made to fail on purpose

`09_seeds.json`, written by `seed_checks.py`. **8 of 8 credited**, plus a
clean control that must come back clean before any verdict is read.

The harness obeys two rules learned expensively. **F-069:** it copies the
folder into a temp directory, mutates the copy, runs there, and then verifies
every tracked artifact in the real folder is byte-identical to its pre-suite
sha256 -- an earlier harness elsewhere left a seeded regioisomer anchor in a
real folder, written by a run that exited non-zero, so the non-zero exit
*concealed* the write. **Section 10:** verdicts are read from the child's
`returncode`, one process per seed, and a seed is credited only when its
**metric moves**; a seed whose metric line is absent is refused, because an
exception on the way to a check looks identical to the check firing.

| seed | route | result |
|---|---|---|
| S1 AST guard | refusal | seeded `import assemble`; fires with the right message |
| S2 template guard | refusal | seeded a 60-char literal lifted from an enumerated SMILES; fires |
| S3 C10 anchor regioisomer | metric | exact matches **2 -> 1**. The defect has **identical formula C56H91NO5 and delta exact mass 0.000000 Da** and is visible only to the InChIKey -- which is what proves the gate tests connectivity |
| S4 retinoid E/Z | refusal | 13CRA written all-trans; the geometry assertion fires |
| S5 Z-bond position | refusal | 13CRA and 9CRA swapped; both have 3E/1Z so a count alone would pass |
| S6 extra CH2 in the AT1 tail | metric | HRMS agreement **39 -> 29** |
| S7 heat-map grid shift | refusal | half-cell shift; `heatmap.match`'s distance diagnostic fires |
| S8 bar-window widened 6 px | metric | bars read **77 -> 0** |

Two seeds were **not credited on the first run and the seeds were retargeted,
not the guards weakened**. S1's `import assemble` sat at module level and
raised `ModuleNotFoundError` before the guard ran, so the expected message
never appeared; it now sits inside a function that is never called, which the
AST still sees. S6 tripped a *different* check -- the ESI-MS convention
determination, which chose between the protonated and neutral-H readings on
the **worst** residual, so one wrong structure made both conventions fit badly
and the choice flipped. That was a real defect in the check: "which convention
does the SI use" is a property of the column, while "is any structure wrong"
is the HRMS gate. The convention is now chosen on the **median** (0.32 mDa for
the neutral-H reading against 0.75 for the protonated, margin 0.43 mDa) and
the strict worst-case bound remains the gate. S6 then moved the metric it was
aimed at.

## Reproducibility

Every stage reads only from `source/` and from files this folder writes; no
script reads an input from a session scratchpad (CLAUDE.md section 7). Run
order:

```
python3 -I stage01_si_text.py        # SI .docx -> 00_pages/si_text.txt + media census
python3 -I stage02_blocks.py         # -> 01_blocks.csv
python3 -I stage03_enumerate.py      # -> 03_enumerated.csv   (39/39 HRMS)
python3 -I stage03b_corpus_overlap.py# -> 02_corpus_overlap.json
python3 -I stage04_values.py         # -> 04_*.csv, 04_geometry.json
python3 -I stage06_verify.py         # -> 07_qc.json
python3 -I stage05_formulation.py    # -> 05_values.csv   (reads 07_qc.json)
python3 -I stage07_emit.py           # -> 06_rows.csv
python3 -I make_reference.py         # -> reference_structures.csv
python3 -I seed_checks.py            # -> 09_seeds.json
```
