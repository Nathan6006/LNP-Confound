# PROVENANCE — maohq_38583366

Zhu Y, Cai SS, Ma J, Cheng L, Wei C, Aggarwal A, Toh WH, Shin C, Shen R, Kong J,
Mao SA, Lao Y-H, Leong KW, Mao H-Q. *Optimization of Lipid Nanoparticles for
Gene Editing of the Liver via Intraduodenal Delivery.* **Biomaterials** 308
(2024) 122559. PMID **38583366**, PMCID **PMC11099935**,
DOI **10.1016/j.biomaterials.2024.122559**, PII **S0142961224000930**.

Identity, first author and endpoint resolved from Europe PMC's REST `core`
record in one call. Original-data paper by the stated discriminator:
`abstractText` is **1,473 characters** (not 0) and the JATS carries
*2. Materials and Methods* and *3. Results and Discussion* with named
subsections. `pubType` was not used.

**First author is Zhu Yining**, so the key derives to `YZ_2024` — which is
already a *different* paper in the corpus (PMID 38082180). See `DECISIONS.md`
D-3.

---

## THE ENDPOINT

In vivo firefly-luciferase expression from **plasmid DNA** LNPs, read as total
bioluminescence flux (photons/s) by **ex vivo IVIS** after intrahepatic,
intravenous or **intraduodenal** injection; plus an in-vitro HepG2 luciferase
screen reported as ln(RLU) at 24 h. The gene-editing arm reads Cas9 excision →
tdTomato (confocal and FACS) and TIDE indel %.

## HOW MANY DISTINCT IONIZABLE LIPIDS: **ONE**

`DLin-MC3-DMA`, already in the corpus on **3,098** rows under two spellings.

| | |
|---|---|
| distinct ionizable-lipid **structures** | **1** |
| distinct helper-lipid structures | 6 (DOTAP, DDAB, DOPE, DSPC, 14PA, 18PG) |
| **formulations** in the library | **1,080** = 3 × 5 × 4 × 3 × 6, *proved*, equal to the deposited row count |
| **individually measured** in vitro | 192 printed in SI Tables S1–S6 (of 1,080) |
| **individually measured** in vivo | **20** (Fig 4B), plus 30 organ values (Fig 4D), 12 dose points (Fig 4E), 12 route/organ (Fig 4G) |
| design-supported but **not attributed** in vivo | 68 cluster/sub-cluster measurements (Fig 3A/3B/3C/3D) over pools of 8 or 4 |
| **emitted** | **0** |

The 1,080 are formulations of one lipid. The varying one-hot group is the
*helper*.

## Value sources

`00_value_sources.json` is the gate and was written before any digitiser.
Summary: **tier 2** for the 192 SI values (parsed as a real `w:tbl`, header
asserted cell by cell); **tier 5** for everything else, with the reason tiers
1–4 fail stated per panel. mmc2/mmc3/mmc4 were probed on the Elsevier CDN in
four extensions each — **twelve HTTP 404s**, so `mmc1.docx` is the only
supplement. The Europe PMC PDF route returned **HTTP 403 as a 5,522-byte HTML
page**, rejected on magic bytes (`3c21444f43`), so there is **no page rect and
no derivable native dpi**; `axis.native_tol` refuses a guessed dpi by design and
was therefore not used. All tolerances are stated in each panel's own pixels and
converted to a factor per pixel.

Figures came from `ars.els-cdn.com/.../gr<N>_lrg.jpg`, 3059×3180 to 3591×3475 px.

**Accession sweep, explicit null:** 0 hits for `GSE\d+`, SRP, PRJNA, SRX, SRR,
SRA, accession, deposit, github, zenodo, dryad, figshare, osf.io, "data
availab", "Source Data", "code availab" across **78,475** characters
(60,875 JATS + 17,600 SI). `find_repo.py` run twice, `--author Zhu` and
`--author Mao`: no repository for this PMID.

**Untrusted-data posture.** Every download went to a uniquely named,
*asserted-empty* directory in the session scratchpad, and the docx was unzipped
into a second one. No glob copy was used (F-065). Every parser lives in a
different directory from the extracted archive, takes its input path as an
argument, and runs under `python3 -I`.

## The headline finding: the in-vitro arm is already merged

See `DECISIONS.md` D-1 and `02_corpus_overlap.json`. Composition-only join,
**1,080/1,080**, Spearman **1.000000**, ln-unit residual rms **2.8531e-4**
against **2.8868e-4** predicted from three-decimal rounding alone. The paper's
192 SI values are in that column as a multiset **192/192**, each helper's top 32.

## C10 anchor posture — recorded by hand

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` write tracked files other agents read, so **none of
them was run**. The local equivalent is `scripts/stage05_c10.py` →
`04_c10_anchor_check.json`.

Fig 1 draws three components in full. Two are anchored:

| anchor | posture | result |
|---|---|---|
| DLin-MC3-DMA | **DERIVED** from the IUPAC name, head first | exact match, **full InChIKey** `NRLNQCOGCKAESA-KWXKLSQISA-N` |
| DMG-PEG2000 (n_EO = 45) | **DERIVED** from the chemical name, head first from the methyl cap | exact match, **full InChIKey** `NLDCRRFPGOSAHV-UHFFFAOYSA-N` |
| Cholesterol | **NO ANCHOR**, empty `reference_SMILES` | gap recorded, see below |
| DOTAP, DDAB, DOPE, DSPC, 14PA, 18PG | **NO ANCHOR**, empty `reference_SMILES` | this paper draws none of them |

**2 of 2 anchors are exact on the FULL InChIKey.** Both are byte-identical to
the corpus strings once canonicalised — **that is expected and is not the
signal** (F-037). Independence is established at run time, two ways, in
`make_reference.py`:

* `ast_independence()` walks *this module's own parsed AST*, inspecting
  `Import`/`ImportFrom`/`Call`/`Attribute` nodes only, so the ban list cannot
  match itself. Measured: `forbidden_imports_found: []`,
  `forbidden_calls_found: []`.
* `template_provenance()` requires that no string literal of ≥25 characters in
  the module be a substring of any target.
* **Both were made to fail on purpose**, in child processes, against a clean
  control: a forbidden import placed in a never-called function (so the AST sees
  it and nothing executes) and a 40-character slice of the target pasted in as a
  literal. **2 of 2 credited.**

The traversal is **opposite** to the corpus's: the corpus writes these lipids
tail first; this file composes head first, from the dimethylamino nitrogen
outward for MC3 and from the methyl cap inward for the PEG lipid, and every
chain is a one-atom token multiplied by a number, so a chain length exists in
one place and is never a hand-typed run of `C`s.

**The gate is proved to test connectivity, not mass.** Each anchor was also
compared against a seeded defect with **identical formula and identical exact
mass**: an MC3 regioisomer with the ester moved from C19 to C18 of the same C37
tetraene (C43H79NO2, Δ 0.000000 Da) and a reversed myristoyl ester on the PEG
lipid (C122H242O50, Δ 0.000000 Da). Both change the full InChIKey. **2 of 2
credited.**

### What Fig 1 settles about cholesterol, and why there is still no anchor

Fig 1 draws cholesterol with **the C3 oxygen on a bold wedge (3β)** — the exact
feature F-013 inverted on 537 `ZT_2026` rows — plus the Δ5 olefin and an
iso-octyl side chain. That reading is recorded here as prose and **not** turned
into a SMILES. Hand-writing an eight-stereocentre sterol is a worse anchor than
none: one agent's hand-written cholesterol came out a diastereomer with 4 of 8
centres inverted at identical formula, identical exact mass to six decimals and
an identical skeleton InChIKey block. The registry string is **not** counted as
an anchor either — comparing the corpus to the corpus is comparing a value to
itself. **Against the seven empty anchors the guards are VACUOUS, not passing.**

### The PEG bracket

Fig 1 subscripts the polyether bracket **44**, and draws **one further
ethyleneoxy unit outside the bracket** before the methyl cap, so the drawn total
is **45**. The corpus `DMG-PEG2000` string carries exactly 45. Taking the
subscript at face value would have been a 44 Da error that only the InChIKey
catches. Verified by zooming the bracket to 5× native.

## Digitisation

Two panels, both through `curation/lib` (`axis`, `bars`); no geometry was
written here.

**The scale is fitted on major ticks only.** The minor ticks are held out and
used as the test. They are **not** at log10(2..9): this figure's minor ticks are
equally spaced in **linear** value inside each decade (2.8, 4.6, 6.4, 8.2 — an
Origin "4 minor ticks" log axis), so they are *unequally* spaced in the fitted
coordinate and a wrong scale cannot reproduce them by accident.

| | Fig 4B | Fig 4D |
|---|---|---|
| px per decade | 298.992 | 111.196 |
| 1 px, as a factor | 1.00773 | 1.02092 |
| decade gaps differ by | **0.077 px** | — (4 majors, LSQ) |
| held-out minor ticks | **10**, max **0.48 px**, rms 0.257 | **16**, max **0.528 px**, rms 0.298 |
| bars read / refusals | **20 / 0** | **30 / 0** |
| max flank delta | 2.0 px | 1.0 px |

**Three levels, all read from the panel, none assumed**: the white page, the
bar's dotted fill and the dark border stroke. The ink mask is the midpoint of
the first two, and every bar asserts that its three levels are separated by ≥30
`min(RGB)` units before any threshold is taken — the floor is *auditable*, not
tuned (the shape of F-051's repair).

**The value is the stroke CENTRE, by two-sided half maximum** against the page
above and the bar's own fill below. `bars.bar_top` is used as the detector and
the refusal gate; reporting its row — the stroke's upper edge — against a scale
fitted on tick *centroids* would be F-019. The correction is **+0.54 px on
average** (a factor of 1.00418), min −0.153, max +1.745.

### Unprinted anchors

1. **Fig 4B's drawn x-axis rule lands on a round 5,000**, to **0.26 px** — an
   axis limit the authors chose and never printed. The rule is **detected**, not
   hard-coded, so the check can actually see its own input.
2. **Fig 4D's six printed pie percentages** (96.7 / 91.4 / 95.5 / 95.2 / 96.0 /
   87.5 % liver) reproduce from liver ÷ sum of the five organ bars at
   **max |Δ| 0.318 pp, rms 0.140 pp, 6 of 6**. The authors drew those pies and
   never computed them from these bars, so this is external to the read, and one
   misread bar of thirty breaks it. It is also the strongest available evidence
   that Fig 4B — same method, three times finer axis — is read correctly.
3. The text's *"over 87% of the observed bioluminescent flux localized in the
   liver"*: smallest computed share **87.51%**.

### Checks with teeth

`scripts/seed_checks.py`, **6 of 6 seeds credited** against a clean control, one
child process per seed, verdicts read from `returncode`:

| seed | route | metric |
|---|---|---|
| `tick_shift_4b_subtle` (+3 px, residual 2.0 px — *passes* `calibrate`) | gate | decade gap 0.0774 → 6.0774 px |
| `tick_shift_4b_gross` (+6 px) | `axis.calibrate` refusal | — |
| `baseline_shift_4b` (whole tick ladder +8 px) | gate | round-5,000 offset 0.26 → 7.74 px |
| `occlude_bar_4b` | `bars.bar_top` refusal | — |
| `liver_bar_4d` (−30 px on one liver bar) | gate | pie share 0.318 → 76.292 pp |
| `minor_tick_4d` (one held-out minor +5 px) | gate | minor residual 0.528 → 4.825 px |

Two seeds were **retargeted rather than crediting them weakly**: moving the
drawn baseline rule hit the level-separation guard before reaching the gate it
was aimed at, and a major-tick shift in Fig 4D moved the minor residual only
0.528 → 1.057 px, under the gate — widening the gate to credit it would have
been widening the gate, not the seed. Both are recorded in
`01_seed_results.json`.

**F-069:** the harness writes only into a temp directory and asserts that
`01_fig4b.json` and `02_fig4d.json` are byte-identical before and after the
suite. Both stages were then run twice more and produced byte-identical output.

## Licence — **FAILS §8**

NIH author manuscript. The JATS permissions block reads verbatim: *"This file is
available for text mining. It may also be used consistent with the principles of
fair use under the copyright law."* Europe PMC: `isOpenAccess = N`,
`license = null`. Crossref lists only Elsevier TDM and accepted-manuscript user
licences. **Not CC0/CC-BY.** No repository exists, so there is no second licence
to check.
