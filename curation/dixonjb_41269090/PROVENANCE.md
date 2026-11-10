# PROVENANCE — dixonjb_41269090

Michalaki E, Jeong K, Chin R, Qi Z, Liebman L, González-Vargas Y, Schrader
Echeverri E, Paunovska K, Muramatsu H, Pardi N, Tamburini BJ, Jakus Z,
Dahlman JE, Dixon JB. *Lipid Nanoparticle Library Screen Reveals Lymphatic
Endothelial Cell-Targeting Lipid Nanoparticle for Delivering Vascular
Endothelial Growth Factor C mRNA after Lymphatic Injury.* **ACS Nano**, 2025.
PMID **41269090** · PMCID **PMC12874650** · DOI **10.1021/acsnano.5c12080**.

The title in our tracker was truncated at "Reveals Lymphatic…". Resolved from
the Europe PMC REST `core` record (`EXT_ID:41269090 AND SRC:MED`, hitCount 1) on
2026-11-10, together with the research-article discriminator: `abstractText`
length **1,873** and both a Methods and a Results section, so this is ordinary
research with original data, not a commentary.

**Nothing in this folder has been merged into `data/`.** No file outside this
folder was written.

## 1. The endpoint

**Delivery readout:** the **percentage of lymphatic endothelial cells
(CD45−/CD31+/PDPN+) expressing the protein translated from the delivered aVHH
mRNA**, measured by flow cytometry 12–16 h after intradermal injection of a
single named LNP formulation. That is cargo expression, not particle
accumulation.

**Cargo expression does exist in this paper**, and this matters, because
"lymphatic targeting" here is reported *both* ways:

| | what it measures | admissible |
|---|---|---|
| Figures 2B–D, 4A, 4B, S6 | % of cells positive for the aVHH **protein** | **yes** |
| Figure 5 | IVIS radiant efficiency of **Cy3-conjugated miRNA** | no — label biodistribution, no expression step |
| SANDS screen | normalised **DNA-barcode counts** in sorted cells | no — pooled, relative, and unpublished |
| Figure 7 | VEGFC **ELISA** | cargo protein, but normalised to a different-cargo control group and plotted only |

The screen proper is a **pooled barcode assay** (SANDS): 99 QC-passing LNPs
co-injected, aVHH-high LECs sorted, barcodes sequenced. Its per-LNP values are
published nowhere — not in a table, not in text, and not recoverable from
Figure S3a, which is a histogram of the delivery *distribution* whose bin
heights are counts of LNPs. The admissible subset is the single-LNP validation,
exactly as the brief anticipates. The admissibility argument is `DECISIONS.md`
**D-1**.

## 2. Where the values are, and the tier

`00_value_sources.json` is the intake gate and was written before any extraction
code in this folder.

Every emitted value comes from **SI Tables S2, S3 and S4** as
`aVHH+ cell count / cell count`. Those three tables are **pasted into the SI as
bitmaps**: their pages carry 286, 337 and 344 characters of text — caption and
ACS download watermark only — and `pdfplumber.find_tables()` finds nothing on
them. Measured native resolutions: **198.5, 249.7 and 249.8 dpi**
(`px / rect_pt × 72`, never guessed). So the tier is **5, printed numerals**:
the mildest form of pixels, with no axis, no scale, no origin and no colour
ramp, so error classes 20, 23, 27, 28, 32 and 43 cannot arise — nothing is
measured, only read. `stage00_scan_for_values.py` re-derives that claim and
exits non-zero if it is false.

Each table was read **twice**, the second time from a different crop at a
different magnification: **42/42** cells of Table S2 and **40/40** cells of the
LECs columns of S3/S4 agree.

The composition tables are tier 4: SI page 2 (Screening 1) and SI page 21
(Table S1) are native text tables and are extracted, not read.

## 3. The three numbers

| | |
|---|---|
| **Distinct ionizable-lipid structures, individually administered** | **2** — cKK-E12 and DLin-MC3-DMA; **1** excluding the benchmark |
| **Distinct ionizable-lipid structures, whole library** | **NOT DETERMINABLE.** Five are named (cKK-E12 + the stereopure lipomers C12S, C12R, C13S, C13R, "provided from the Dahlman lab"). Figure 1A *draws* all of them, but the only obtainable copy is 708 × 391 px for twelve panels — about 110 × 57 px each — with no SI structure, formula, HRMS or citation. |
| **Distinct formulations** | 150 claimed formulated, 99 claimed administered, **30 verifiable** (Screening 1's native table; Screenings 2 and 3 are bitmaps at 122.7/108.6/99.2 dpi and are not transcribed because no emitted value depends on them), **7 individually administered** |
| **Individually measured delivery values emitted** | **28** = 18 (six leads × three nodes) + 10 (two LNPs × five sites) |

**A library of one lipid.** Screening 1 is a complete factorial, proved rather
than asserted: 5 helper lipids × 6 mole-ratio sets = 30 = the row count, every
cell present exactly once, one lipomer throughout. That is why all six leads
carry cKK-E12 and why `DECISIONS.md` **D-2** exists.

**A name token is not a chain length, again.** Figure 1A's "R Chains" key reads
C12S = C₁₂H₂₅ and C13S = C₁₃H₂₇. The arm is CH₂–CH(OH)–R, so an arm of "C12"
carries **14** carbons. The number is the R-chain length, not the arm length.

## 4. Corpus overlap

Full InChIKey, skeleton key and `IL_name` against the live
29,195 × 79 / 1,803 × 79 corpus (`02_corpus_overlap.json`).

- **cKK-E12** is already present, exactly, in **three** experiments — and under
  **three different names**: `CKK-E12` in `AC_2025` (216 rows), `MD-1` in
  `BL_2023` (6 vivo rows), and `other` in `JW_2024` (19 rows). C17 pins one name
  to one structure; the reverse is not checked, and keying on `IL_name` would
  have missed this lipid in two of the three. Our rows therefore carry the
  SMILES the corpus already uses, not a fresh one.
- **DLin-MC3-DMA**: exact hits in 11 experiments.
- **Skeleton-only hits: 0** for both — no stereoisomer-only match anywhere.
- **Near-miss sweep, labelled a SUPERSET**: a deliberately loose
  `[NX3;!$(N=*)](CC(O)[CH2])` matches **1,584** distinct corpus structures
  across 19 experiments. That is a whole chemotype, not a near-miss of cKK-E12;
  it is reported so the exact answer is not a silent zero.

**Already-merged arm test.** `AC_2025` (COMET) varies the *same* four formulation
parameters on the *same* ionizable lipid, so it is the natural place for this
screen to be a duplicate. Matched formulation-by-formulation on
(helper lipid, IL/CHL/PEG/HL mole percents): **26 corpus formulations vs our 6,
zero collisions.** The endpoints differ anyway — `AC_2025` is in-vitro DC2.4
luminescence — but the formulation test is the stronger one and it is clean.

**Token-pattern test.** The library labels are `LNP1…LNP30`, **contiguous with
no gaps**. A contiguous 1..N numbering is the `AC_2025` failure mode: the shared
lipidoid amine catalogue is sparse and non-contiguous, so a corpus name lookup
on "LNP7" would return a plausible wrong molecule. **No corpus name lookup is
performed on these tokens**; they live in `Formulation_ID` only, never in
`IL_name`. The structural attestation instead is that both ionizable lipids are
catalogue reagents with a supplier and catalogue number printed in the Methods
(BroadPharm BP-29590; MedKoo 555308), resolved through
`provenance/reagent_registry.py` and `curation/corpus_components.py`'s guarded
`smi()`, which raises on an ambiguous name.

**Does the benchmark top the screen? No.** LNP7 beats MC3 by **2.39×** at the
draining popliteal node and **4.06×** at the popliteal lymphatic vessel on our
derived values (6.2× and 3.7× on the paper's printed ones). Relevant to
`CLAUDE.md` §8's open decision on benchmark-only pairs — this source would be
benchmark-only for every cross-source pair but one.

## 5. One z-block, and why

All 28 rows are one experiment, `EM_2025`, with
`value_normalization_unit = experiment`.

The first draft used two blocks, one per cohort. That is error class 22:
standardising sub-blocks separately **deletes** the offset between them, which
is what F-016 cost 5,040 rows. Here the offset is measurable rather than
hypothetical — LNP7 is in both cohorts, and at the one site that drains the
injection in both (the popliteal node) they agree to within **11.4%** (37.34%
vs 33.52%, ratio 1.114). So the cohorts sit on the same scale and one z is
right.

What *does* differ between cohorts is which nodes drain the injection: study A
injected both paws, study B only one hindlimb. That is a fact about the site,
not a scale offset, so it is carried in `Model_target` as a `_nondraining`
suffix rather than deleted by a separate standardisation. Without it LNP7 would
carry two axillary values **14.4×** apart under one key. After the fix exactly
**one** duplicate `(IL_name, Formulation_ID, Model_target)` key remains, and it
is the true replicate pair above.

## 6. Anchor posture — recorded by hand

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run**: they write tracked files other
agents read. This section is the hand record.

`reference_structures.csv`, header exactly
`lipid_name,reference_SMILES,source_figure,read_on,read_by,expected_mismatch,note`:

| lipid | reference_SMILES | posture |
|---|---|---|
| `DLin-MC3-DMA` | present | **derived**, exact InChIKey match |
| `cKK-E12` | **empty** | the absence is visible on purpose |

**The MC3 anchor is derived, and honest about what it tests.** It is not read
off the drawing. `make_reference.py` walks the IUPAC name that
`provenance/reagent_registry.py` records as that reagent's authority —
*(6Z,9Z,28Z,31Z)-heptatriaconta-6,9,28,31-tetraen-19-yl
4-(dimethylamino)butanoate* — and composes the 37-carbon chain **one carbon at a
time as a string**, reading the alkene locants and the ester locant. No molecule
object is edited anywhere. It agrees with the enumerated structure:
`NRLNQCOGCKAESA-KWXKLSQISA-N`. What that tests is the registry string against
the chemical name by a second route; it does **not** test the paper's drawing.

**Independence is verified by provenance, not by bytes** (F-037). The emitted
anchor *is* byte-identical to the registry string, and that is not the signal —
a correct re-derivation is byte-identical to a correct enumeration once
canonicalised, which is exactly why C10's own `anchor_may_be_a_copy` missed
F-037. Instead `_ast_guard` walks **this module's own AST at run time** and
refuses on any import of `assemble` / `build_new_curation` or any call to
`RWMol`, `AddBond`, `RemoveAtom`, `CombineMols`, `ReplaceSubstructs`,
`EditableMol`, `esterify` or `attach`. It walks the AST rather than scanning
text, because a substring scan flags the module's own docstring.

**The guard was made to fail on purpose.** `--prove-guard` writes a temporary
copy with `import assemble` added and requires rejection: it reports
`REJECTED (['import assemble'])` and exits 0 only on that.

**The gate was shown to test connectivity.** `--seed-regioisomer` moves the
ester from C19 to C18. Formula `C43H79NO2` **both**; exact mass **641.611081**
both, Δ **0.000000**; InChIKey `NRLNQCOGCKAESA-KWXKLSQISA-N` vs
`DPHDIRYUGZIOOZ-PNOCIVJDSA-N`. Only the key sees it. Exit code 1.

**Why cKK-E12 has no anchor.** Figure 1A *does* draw it — my first pass recorded
"NO structure disclosed" from the Methods sentence alone, and that was wrong,
which is itself the lesson about reading a paper from one sentence. But the
panel is **130 × 70 px** inside a 708 × 391 px figure, the only copy obtainable:
Europe PMC's `fulltextRepo` returned an HTML error page and the PMC OA service
returned 404 / an HTML shell on both hosts. No bot check was bypassed and no
paywall was crossed. At that size the 2,5-diketopiperazine with its two ring
carbonyls and two NH, the two tertiary nitrogens, the four CH₂–CH(OH)–R arms and
their `C10H21` labels are legible and **consistent** with the emitted structure
(C₆₀H₁₂₀N₄O₆, four 2-hydroxydodecyl arms) — but the methylene count of each
lysine side chain, the one thing a C10 anchor exists to test, is not resolvable.
Counting it by eye is what this project forbids. **No anchor is better than a
guessed one.**

**Stereo assertions, stated rather than counted.** `DLin-MC3-DMA` has **no**
stereocentre, so "stereocentres preserved" is **vacuous** for it; its four C=C
bonds do carry explicit Z geometry and that is asserted separately, because CIP
centre counting is blind to double-bond geometry. `cKK-E12` has **2**
stereocentres and the stored corpus string leaves **both unassigned** — not
vacuous, a real gap, and the reason every cKK-E12 row carries
`confidence = medium` with the reason in `curation_note`.

## 7. Flags, on the row

Written by `curation/lib/emit.emit_rows`, the one audited writer, with per-item
callables. 28 rows: `confidence` **medium ×23** (cKK-E12, undeclared
stereochemistry) and **high ×5** (MC3, from the registry);
`positional_ambiguity` 0; `structure_assumed` 0; `value_is_bound` blank;
`value_normalization_unit` `experiment`; `curation_note` populated on **28/28**.

Both refusals were exercised, **each in its own process, verdict read from the
return code** (F-012's shape), and both returned **3**:

- a literal `confidence="high"` → refused (error class 29, F-011/F-025/F-046);
- a caveat with an empty note → refused (F-024).

The control run returns 0.

## 8. Checks that move a number

`07_qc.json`. Five further seeds each perturb one input and the **agreement
count must move**, not merely raise; each runs in its own process and the
verdict is the return code. All five moved exactly their own metric and failed
exactly their own check:

| seed | metric | moved |
|---|---|---|
| `s2_digit` | Table S2 double-read agreement | 42 → **41** |
| `s3_lec_digit` | S3/S4 LECs double-read agreement | 40 → **39** |
| `screening1_helper` | Screening-1 distinct factorial cells | 30 → **29** |
| `mc3_ratio` | MC3 vs Onpattro max |Δ| mol% | 0.0 → **0.5** |
| `z_value` | z round-trip observed/predicted ratio | 0.9966 → **1.9414** |

**A round number the authors used but never stated.** Table S1 gives the MC3
arm as 50 : 38.5 : 10 : 1.5 (ionizable : cholesterol : DSPC : PEG). That is the
published patisiran/Onpattro molar ratio exactly, max |Δ| **0.0** mol%. The
paper calls its benchmark "an FDA-approved hepatocyte-targeting LNP formulation
in Onpattro" but never says the ratio *is* Onpattro's. External to anything we
derived. (`mole_percents_sum_to_100` is recorded beside it and labelled
**WEAK**, because it is close to a tautology — the authors composed the table to
sum to 100 and nothing we derived enters it.)

**Two independent composition tables.** The six leads appear on SI page 21 and
SI page 2; **48/48** field comparisons exact, and `stage03_enumerate.py` raises
rather than picking if they disagree.

**Round-trip tolerance derived, not rounded.** `Experiment_value` is stored to
6 dp, so the z carries a 5 × 10⁻⁷ half-ulp. Observed max reconstruction error
**2.0164 × 10⁻⁵** pp against **2.0233 × 10⁻⁵** pp predicted from that rounding
alone — ratio **0.9966**. A round 10⁻⁹ bound was tried first and could only ever
fail.

## 9. Licence

| | |
|---|---|
| Article | **CC BY** (Europe PMC `license: "cc by"`, `isOpenAccess: Y`, `inPMC: Y`). **Compatible** with `CLAUDE.md` §8. |
| SI PDF | same article licence. It carries an ACS "Downloaded from pubs.acs.org … by PHILLIPS EXETER ACADEMY user" institutional watermark on most pages; that is an access stamp, not a separate licence. |
| PMC OA assets (`source/pmc_oa/`) | same article licence, fetched through the PMC Open Access asset route after NCBI `elink` confirmed a real `pubmed_pmc` **self**-link (`12874650`), not merely `pubmed_pmc_refs`. |
| Author repository | **none exists.** `find_repo.py` returned 0 candidates on every author term tried (Dixon, Michalaki, Dahlman, Paunovska); its single code hit is an unrelated laser-physics repo matched on the literal PMID. The article names a "custom-built Python script for SANDS analysis" with no URL and carries no Data Availability statement. So there is no second licence to clear — and no repository means **tier 3 was unavailable**, which is why the values are tier 5. |

Caveat worth keeping: the licence is read off the Europe PMC record, not off the
ACS landing page, which was not visited.

## 10. What could NOT be verified

1. **The library's distinct-structure count.** Four of the five ionizable lipids
   have no readable structure. Figure 1A draws them at ~110 × 57 px.
2. **The paper's own plotted statistic.** No per-animal data is published, so
   the mean of per-animal percentages is not derivable at any tier. The emitted
   ratio of means disagrees with the 17 printed percentages at Spearman
   **0.820641**, 10/17 within 3 points, worst ratio **0.40**
   (`DECISIONS.md` D-4, `08_contradictions.csv` C-3 and C-4).
3. **Dose.** Not stated for any single-LNP study; blank on all 28 rows (D-6).
4. **The cKK-E12 stereocentres.** The commercial reagent is a single
   diastereomer; the corpus string is flat. Kept flat for corpus consistency
   (C17), flagged `confidence = medium` with the reason on the row.
5. **Screenings 2 and 3** of Figure S1 (120 of the 150 formulations) are
   bitmaps and were not transcribed; no emitted value depends on them.
6. **Figure 2A's per-lead diameters and PDIs**, and every Figure 5–9 panel, were
   not read; none is a delivery endpoint.
7. **The 150 and 99 counts** are the paper's claims. 30 are verifiable.
