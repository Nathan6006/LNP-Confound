# Wave-1 adversarial review — `curation/maohq_39375194` (LC_2024 / LC_2024_2, 5,112 rows)

Reviewer: read-only Wave-1 agent, 2026-11-10. Target is 71% of everything this run staged.
Method: values, structures and metadata re-derived from the primary sources (a fresh clone of the
authors' repository, the SI PDF, and the Europe PMC/NCBI records), never by re-running the folder's
scripts. Every corpus query carried a positive control.

Controls used, all reproducing the brief's figures exactly: `DSPC` by name **7,038** rows
(vitro 6,401 + vivo 637); `DLin-MC3-DMA` by full InChIKey `NRLNQCOGCKAESA-KWXKLSQISA-N`
**3,098** rows across two spellings (`DLin-MC3-DMA`, `Dlin-MC3-DMA`); cholesterol by
`CHL_name` **30,397** rows.

**Bottom line.** I could not break the numbers. All 5,112 emitted values, all four molar percents,
the helper assignment and the N:P re-derive from the authors' deposit to float noise, and I found
**four external anchors in the publication that the curation did not use — all four pass**. No
emitted value and no emitted structure is wrong. Eight findings follow, all in documentation,
flagging posture and verification integrity; the largest is a flag that was computed and then died
at the emit step (the F-001 shape), on 748 rows.

---

## CONFIRMED CORRECT

**C1 — Every value re-derives from the deposit.** I built my own long table from
`Raw_Data/1080_LNP_OHE.csv` and `Compiled_validation.xlsx`, computed my own mean/sd, and joined on
`Formula label` × cell.

| | n | my μ(ln) | my σ(ln) | max abs Δz vs emitted | rms |
|---|---|---|---|---|---|
| `LC_2024` | 4,680 | 4.341323 | 2.750764 | **1.465e-14** | 4.84e-15 |
| `LC_2024_2` | 432 | 5.546311 | 2.451561 | **6.66e-16** | 3.77e-16 |

Helper name and N:P agree on all 5,112. Molar percents: exact (max diff 0.0) on the validation set,
max 4.55e-5 on the training set — that is 4-dp rounding, not an error. 0 null values, 0 duplicate
`(Experiment_ID, Formulation_ID, Model_type)` keys.

**C2 — Repository provenance verified upstream.** Fresh `git clone` into a uniquely named,
verified-empty scratch directory (`exists_before_count = 0`), nothing executed:

- commit **`82cee58057a6168b61739c5863c6ac8a4ff19a7e`** is the current `HEAD`, dated
  **2024-10-31 14:22:46 -0400** — matching `PROVENANCE.md` to the second;
- `Raw_Data/1080_LNP_OHE.csv` md5 **`98baa22d8d8e9a92648e1c9b9904198f`** ✓;
- `Raw_Data/Compiled_validation.xlsx` md5 **`5e1ea0cd1e8f7c114ff355d58c4206b4`** ✓;
- `source/nn4c07615_si_001.pdf` sha256 **`27ef09af…c5f39e8a`** ✓.

Note for the record: this folder's source assets are **not** in
`provenance/deleted_sources_2026-11-10.json` (that manifest covers the 16 zero-yield folders; the
nine `maohq` entries in it belong to `maohq_38583366`, a different paper). `source/github/` is still
on disk and byte-identical to upstream.

**C3 — The 72 `HepG2` rows are NOT a `YZ_2022` duplication.** They are the de-novo validation arm.
My first join returned 0 overlap everywhere, which was *my* bug: `YZ_2022` stores
`IL_to_nucleicacid_massratio` ∈ {8.48, 16.96, 25.44} and no charge ratio, while this folder stores
`IL_to_nucleicacid_chargeratio` ∈ {4, 8, 12}; the map is mass ratio = 2.12 × N:P. After reconciling
that, joining on helper + four percents to 1 dp + N:P:

| arm | n | distinct keys | overlap with `YZ_2022`'s 1,080 |
|---|---|---|---|
| `LC_2024` B16F10 / HEK293T / PC3 | 1,080 ea | 1,080 | **1,080** (same library, expected) |
| `LC_2024` Neuro2a / ARPE-19 | 720 ea | 720 | **720** |
| `LC_2024_2` all six cells | 72 ea | 72 | **0** |

The validation library is a complete **6 helpers × 3 × 2 × 2 factorial** at levels disjoint from the
training grid (`HL/(IL+HL)` {1.5, 20, 40}, `(IL+HL)` {30, 70}, `PEG/(Chol+PEG)` {0.5, 7}, N:P **10**
against the training N:P {4, 8, 12}). Nothing merged is duplicated; the exclusion was applied to
exactly the 1,080 training HepG2 values and nothing else.

**C4 — D-5 reproduced independently.** 1,080 one-to-one composition join (`validate='one_to_one'`
passed). Spearman(repo `RLU_HepG2`, corpus `YZ_2022` `Experiment_value`) = **1.000000**. Affine fit
residual expressed in ln units: rms **2.85314e-4**, max 5.335e-4, against **2.88675e-4** expected
for ±0.0005 three-decimal rounding — ratio 0.988. The curation's numbers are exact.

**C5 — One-hot inversion verified, not taken on trust.** `Ionizable_Lipid` has exactly **one**
distinct value (`Dlin-MC3-DMA`) on all 1,080 training and all 72 validation rows. The one-hot group
`{18PG, 14PA, DOPE, DSPC, DOTAP, DDAB}` sums to exactly 1 on **1,080/1,080** and **72/72** rows;
each column sums to 180 (training) and 12 (validation); `idxmax` inversion equals the file's own
`Helper_lipid` categorical with **0 mismatches** in both files. The varying group is the helper, as
claimed.

**C6 — The 792 vs 1,152 gap is in the primary file.** `RLU_N2a` and `RLU_ARPE19` are blank on all
180 `14PA` and all 180 `DDAB` formulations — whole-arm, 1,080 − 360 = 720. Not a curation artefact.

**C7 — The z block is right (not error class 22).** Each `Experiment_ID` is separately mean 0 /
pop sd 1 (`LC_2024` n=4,680, μ=−1.82e-17; `LC_2024_2` n=432, μ=+6.58e-17). The packet's single
pooled "mean −1.57e-17, sd 1.000" is uninformative — pooling two standardised blocks *always* gives
0/1. What matters is that cell types are **pooled**, so the between-assay offsets survive: per-cell
mean z runs HEK293T **+1.165** → PC3 **−0.668**. `Experiment_value_within_model` is blank, which is
correct — `provenance/build_dataset.py:868` computes it at merge, and the corpus column is populated
on 29,195/29,195 rows.

**C8 — The two-ID split is defensible, and I tested the alternative.** Separate z does delete a
level difference between the libraries, but that difference is **not** a scalar: validation minus
training mean ln is B16F10 **+4.07**, PC3 +2.57, HepG2 +2.28, Neuro2a +1.52, ARPE-19 +0.17,
HEK293T **−0.95**. Pooling would not have been more faithful. The maxima are near-identical across
libraries (e.g. Neuro2a 10.44 vs 10.72), so the validation library raises the floor, not the
ceiling — consistent with a new design region rather than a selected subset.

**C9 — Descriptors.** All 14 published descriptor columns plus the 5 ML-change-log columns recompute
exactly from the emitted `IL_SMILES` and match the **3,081** in-vitro corpus rows carrying the
identical MC3 SMILES on every column. (`03_structures_qc.json`'s "3,081 rows / 5 experiments" is the
in-vitro scoping and is correct; 3,098 / 11 experiments is vitro+vivo. I checked this before
reporting it as a discrepancy — it is not one.)

**C10 — Anchors reproduce.** All seven stored `reference_SMILES` give exactly the InChIKeys recorded
in `04_c10_anchor_check.json` when I re-canonicalised them myself, and all seven are byte-identical
to the emitted structures (expected, per F-037). See E-7 for the posture problem.

**C11 — Excipient consistency (C17).** Seven of eight emitted excipients are byte-identical to the
single structure the corpus already stores: DSPC (7,038 rows), DOPE (17,737), DOTAP (2,564),
DDAB (360), 14PA (360), Cholesterol (30,397), DMG-PEG2000 (22,973) — each with exactly one distinct
corpus SMILES. **18PG is the only divergence**, which is the flagged conflict.

**C12 — Vocabulary and schema.** 79/79 columns in corpus order. No new value in any controlled
column except four new `Model_type` cell lines (`ARPE-19`, `B16F10`, `Neuro2a`, `PC3`).
`Aqueous_buffer = acetate`, `Mixing_method = handmixed`, `Experiment_method =
luminescence_normalized`, `Cargo = pDNA` all already exist. Quality columns conform to
`schema_quality.json` with 0 violations; **0** caveat rows have an empty `curation_note`.
`LC_2024` and `LC_2024_2` are both free (54 corpus IDs checked).

**C13 — Licence (F-059) confirmed from primary records, not landing pages.** Europe PMC core record
for PMID 39375194: `isOpenAccess = N`, `license = null`, PMCID `PMC11512640`, abstract 1,770 chars.
The repository `LICENSE` is the verbatim GNU GPL v3 ("GNU GENERAL PUBLIC LICENSE / Version 3,
29 June 2007", 553 non-empty lines). Neither is CC0/CC-BY. The asymmetry `PROVENANCE.md` notes is
real: PMID 35879315 (`YZ_2022`) returns `isOpenAccess = Y`, `license = cc by`.

**C14 — D-3's dose and timepoint are right for what is emitted.** I pulled the full JATS via NCBI
efetch (Europe PMC's `fullTextXML` 500s on this closed-access record; `efetch.fcgi?db=pmc&id=11512640`
works). Methods, verbatim: *"4 μL of a particle suspension at 25 μg pDNA/mL was pipetted into 100 μL
of culture media"* = **0.1 μg/well**, and *"After 48 h of incubation"*. Also *"with pipette mixing
**without further purification**"* → `handmixed` ✓ and blank `Dialysis_buffer` ✓, and *"25 mM
magnesium acetate buffer (pH 4.0)"* → `acetate` ✓. The contradiction with PMID 35879315 is real but
the emitted values are this paper's own numbers. For a z-scored value the dose matters only as
metadata; it would matter if anyone pooled these rows with `YZ_2022`'s.

**C15 — Four external anchors in the publication, all passing, that the curation did not use.**
`CLAUDE.md` §5 asks for the paper's own cross-checks. `07_qc.json` contains no entry for any of
these; I ran them.

1. **A printed point value.** The paper: *"The two exceptions were N2a and ARPE-19 cells, which were
   the most highly correlated cell lines with a Spearman's rank correlation coefficient of **0.71**."*
   Measured on the deposit: Spearman(`RLU_N2a`, `RLU_ARPE19`) = **0.7144** over the 720 common
   formulations, and it is the **maximum of all 15 pairs**. (Rank-invariant to their floor/min-max,
   so the check is robust: floored at 1.5 it is 0.7136 and still the max.)
2. **A printed band claim.** *"low (r < 0.35) or moderate rank correlation (0.36 < r < 0.67)"* for
   the other pairs → **14 of 15** pairs land inside those bands, none in the (0.35, 0.36) gap, and
   the 15th is the 0.714 exception the text names.
3. **A printed qualitative ordering.** *"highly conserved while using DDAB, moderately conserved
   using DOTAP and DOPE, and weakly conserved in DSPC, 18PG, and 14PA."* Mean cross-cell Spearman
   per helper: DDAB **0.699** > DOTAP **0.523** ≈ DOPE **0.511** > DSPC **0.334** ≈ 18PG **0.322** >
   14PA **0.268**. A **6/6** rank match. This one independently validates the *helper assignment* —
   a mislabelled helper would scramble the grouping.
4. **A figure, for the validation arm.** Supplementary Figure 3 ("In vitro validation of refined
   models (n = 72)") plots ground-truth transfection per cell on a 0–12 ln axis. My eyeball x-ranges
   off a 200-dpi render vs the xlsx min/max: B16F10 (3.4, 11.0) vs 3.36–11.60; HepG2 (1.0, 8.7) vs
   1.06–8.63; PC3 (1.4, 9.5) vs 1.29–9.57; HEK293 (3.3, 11.8) vs 3.32–11.67; N2a (2.9, 10.7) vs
   2.90–10.72; ARPE19 (0.6, 8.6) vs 0.60–8.58. All six match, and the ranges are distinct enough
   (HepG2 max 8.6 vs HEK 11.7) that a column swap would show.

Together these validate the values, the **cell-line → column mapping** and the **helper assignment**
against the publication, independently of the repository. `DECISIONS.md` D-6's cross-cell Spearmans
also reproduce exactly (Neuro2a 0.494, B16F10 0.488, HEK293T 0.476, PC3 0.447, ARPE-19 0.330).

**C16 — The seeds are genuine.** All 13 file seeds name a specific metric, record a baseline and a
seeded value that differ, and a non-zero returncode, with a passing control; all 5 emit seeds raise
with a named refusal. This satisfies "a new check is not done until it has been made to fail on
purpose". See E-8 for the one class not covered.

**C17 — The 18PG drawing, re-read by me.** I extracted SI page 2's image `xref 68` at native
resolution (786×145, paired to the `18PG` label by y-overlap: label y 639.1–651.4, image rect y
627.6–662.4) and counted carbons on both tails: C1 carbonyl, seven vertices to the cis C9=C10, seven
vertices to C17, terminal C18 — **18:1 cis-Δ9 in each tail**, against a DSPC control (`xref 67`)
drawn fully saturated. The head glycerol carries no wedge (rac); the backbone C2 carries a bold
wedge to H. The curation's reading of the drawing is correct, and the emitted
`DSNRWDQKZIEDDB-GCMPNPAFSA-M` is C42H78O10P⁻, 773.5338, two C=C, backbone CIP **R**, head glycerol
unspecified — exactly as documented.

---

## ERRORS FOUND

### E-1 — The packet is wrong: it names 3 of the 11 experiments that already hold MC3

**File:** `provenance/review/brief_staged.py:177`, and every packet in `round_2026-11-10/`.

`where = sorted({f"{e}:{nm}" for e, nm in ck[k]})[:3]` truncates to three and prints the list with
**no ellipsis and no count**, so it reads as complete. The packet says MC3 is already in the corpus
as `AC_2025:DLin-MC3-DMA, BL_2023:DLin-MC3-DMA, NG_2025_2:DLin-MC3-DMA`. It is in **11 experiments,
3,098 rows**: AC_2025, BL_2023, NG_2025_2, SB_2024, SP_2020, WT_2026, XD_2026, **YZ_2022**,
**YZ_2024**, ZC_2023, ZH_2023. Alphabetical truncation drops exactly the two experiments that are
*the same library as this curation* — the central fact of the folder. The folder's own
`07_qc.json::near_miss_sweep` lists all 11 correctly, so the repository is right and the packet is
wrong.

**Severity:** cosmetic for the data; material for review quality. A reviewer briefed only from the
packet would not learn the `YZ_2022` relationship from it. **Affects all 29 packets** whose
structure appears in more than three experiments. **Confidence: high.**

### E-2 — D-2's evidence table points the reader the wrong way, in both directions

**File:** `curation/maohq_39375194/DECISIONS.md` D-2, and `08_contradictions.csv` C-2.

The brief asked which half is wrong — the name or the structure. I could not settle it (see
UNVERIFIABLE), but I can correct the balance of evidence, which D-2 gets wrong on both sides.

*The SOPG side is stronger than D-2 says.* D-2 presents it as one sentence in Nat Commun 2022. The
verbatim name *"1-stearoyl-2-oleoyl-sn-glycero-3-phospho-(1'-rac-glycerol) (18PG)"* appears in at
least **three** Mao-lab papers, two published *after* this one:

| PMID | PMCID | journal, year | names 18PG as SOPG |
|---|---|---|---|
| 35879315 | PMC9310361 | Nat Commun 2022 | yes |
| 40592930 | PMC12217393 | Nat Commun 2025 | yes |
| 41124242 | PMC12542943 | Sci Adv 2025 | yes |

PMID 41124242's author list is `Lu X, Zhu Y, Wei C, **Cheng L**, …, Mao HQ` — **this paper's own
first author** co-authored a 2025 paper describing the same six-helper panel as SOPG. (The sentences
are near-identical across the three, so they may be one boilerplate copied forward; I say so rather
than count them as three independent witnesses.)

*The DOPG side is weaker than D-2 says.* D-2 lists three artefacts — SI Table 2's drawing, SI
Table 3's `Double_bonds = 2`, and the deposited CSV column. The Methods identify these as **one
authoring act**: *"the calculated log of the partition coefficient … was calculated using **ChemDraw
(v22.2.0, PerkinElmer)**"*, i.e. the helper descriptor block was computed by hand from drawn
structures. The drawing, the printed table and the deposited column are one ChemDraw file with three
outputs, not three witnesses.

*What I did settle.* The brief's caution — that `Double_bonds` might be an **ionizable**-lipid
descriptor — is wrong, and I can close it with primary evidence. SI Table 3 is titled
**"Helper lipid chemical parameterization"**, its columns are the six helper names, and its row
reads *"Double_bonds | Total # of carbon-carbon double bonds **in the lipid tails** | 2 0 2 0 0 2"*.
A 5/5 positive control confirms the scoping: my own RDKit C=C counts for the five uncontested
helpers reproduce the column exactly (14PA 0, DSPC 0, DOPE 2, DOTAP 2, DDAB 0), `Total_Carbon_Tails`
is 28 for 14PA and 36 for the C18 helpers, `cTPSA` is 0 and `Hbond_D/A` are 0 for DDAB, and
`P_charged_centers` is 0 for the two anionic helpers. MC3 has four C=C and appears nowhere. So
`Double_bonds = 2` for 18PG is a statement about the helper, and it means 18:1/18:1.

*A caution D-2 also misses.* **The project has already changed this once.** The corpus originally
stored 18PG as **DSPG** (18:0/18:0, C42H82O10P⁻, 777.5651, zero C=C). `WT_2026`'s curation rewrote
360 rows to **SOPG** on the authority of the same Nat Commun sentence
(`provenance/CHANGELOG.md:269`, `reagent_registry.py:9-10`, error class 48). The field has now had
three answers — 0, 1 and 2 double bonds — and no external datum has ever been brought to bear.

*Blast radius, re-measured.* The corpus holds one 18PG structure, `ZEFGRNLJASLRBZ-GXFQOPFSSA-M`
(SOPG, 1 C=C), on **1,474** rows: `WT_2026` 1,114, `YZ_2022` 180, `YZ_2024` 180. D-2's "360
already-merged rows" is right for the same-library subset; note that `WT_2026` is also a Mao-lab
paper (Toh WH is a co-author of *this* paper) and its own Methods say only *"ALC-0315, DSPC, 18PG
(sodium salt), DOPC, DOPE, and DOTAP were from Avanti Polar Lipids"* — same shorthand, no name, no
catalogue number. Calling it "a different purchase" is an assumption, so the exposure is plausibly
1,474 rows rather than 360.

*One way D-2's framing is too alarmed.* The `AC_2025` analogy does not carry. There, two molecules
both present in the data swapped their measurements. Here 18PG is a **helper**, constant within each
180-row block, and the block label never moves. The failure mode is one wrong helper structure on
972 (and possibly 1,474) rows — no measurement can land on the wrong molecule.

**Severity:** changes no value; bears directly on a structure decision over 972 new + 1,474 existing
rows. **Confidence: high on the facts; medium on the weighting.** My own read is that the evidence
now leans SOPG slightly more than D-2 implies, and that this should be stated in D-2 rather than
left for the next reader to discover.

### E-3 — `02_corpus_overlap.json` reports the overlap for a structure this folder does not emit, and misses a live DOPG collision

**File:** `curation/maohq_39375194/02_corpus_overlap.json`.

The entry labelled **`"18PG as DRAWN in this paper (HL)"`** carries
`full_inchikey: DSNRWDQKZIEDDB-SQYFZQSCSA-M` and reports 81 `FULL_key_hits` in `ZT_2026`.

That key is **not what the paper draws** (the drawing has no head-glycerol wedge) and **not what
`06_rows.csv` emits** (`DSNRWDQKZIEDDB-GCMPNPAFSA-M`). It is the registry's DOPG *base*, before the
head stereocentre was cleared — the step `03_structures_qc.json::drawn_18PG_derivation` documents
correctly (`specified_stereocentres` 2 → 1). So:

- the reported 81-row `ZT_2026` overlap belongs to a structure this curation does not emit;
- the **actually emitted** 18PG has **zero** corpus rows — the overlap file never says so.

The fact this conceals matters for D-2. The corpus already holds DOPG, under the name **`DOPG`**, in
`ZT_2026` — 81 rows, `…COP(=O)([O-])OC[C@@H](O)CO…`, head stereocentre **defined**. Merging as drawn
puts **two different DOPG stereo-variants in the corpus under two names**, which `C17` cannot see
because it keys on the reagent name. Neither `DECISIONS.md` nor `PROVENANCE.md` mentions `ZT_2026`'s
DOPG at all, although it constrains what a correct D-2 repair looks like: if DOPG wins, the right
move is almost certainly to adopt `ZT_2026`'s existing structure rather than introduce a third
variant.

**Severity:** metadata and decision support. No emitted value or structure is wrong, and the
derivation is correctly recorded elsewhere in the folder. **Confidence: high.**

### E-4 — `08_contradictions.csv` C-6 states the wrong row count

`C-6.what_we_emitted` reads *"Model_type = HEK293T on all **936** HEK rows"*. The file has
**1,152** (1,080 training + 72 validation). **Severity: cosmetic. Confidence: high.**

### E-5 — The plasmid-supplier contradiction is internal to Nat Commun 2022, not between the papers

`DECISIONS.md` D-3's table and `08_contradictions.csv` C-5 present ACS Nano = Aldevron against
Nat Commun 2022 = Translate. PMID 35879315's **Materials** section reads *"**All pDNA was purchased
from Aldevron.**"*; *"Luc pDNA (firefly mLuc, Translate)"* is a second, inconsistent statement in its
own LNP-synthesis paragraph. The two papers **agree** on Aldevron. C-5 uses this as "further
evidence that this paper's Methods paragraph was not written against the HepG2 arm" — that inference
does not survive, because the disagreement is inside the other paper.

**Severity:** cosmetic; no corpus field stores a supplier. **Confidence: high.**

### E-6 — 748 emitted rows sit at or below the paper's own stated blank, and no row says so (F-001 class)

This is the finding I would act on first.

The paper states its own background twice:

> *"Ln transformed RLU values were then floored to the average value of blank wells (around e^1.5
> RLU) to minimize noise at low and near zero transfection measurements."* — Methods, Data
> Preprocessing
>
> *"Transfection efficiency measured by the relative luminescence units (RLU) obtained with this LNP
> library was sufficiently broad (**e^1.5 – e^12** RLU/well) for most cell types."* — Results

Of the 5,112 emitted rows, **748 (14.6%)** carry a value **below ln 1.5**:

| experiment | cell | below 1.5 | of |
|---|---|---|---|
| LC_2024 | PC3 | **441** | 1,080 (40.8%) |
| LC_2024 | B16F10 | 200 | 1,080 |
| LC_2024 | ARPE-19 | 68 | 720 |
| LC_2024 | Neuro2a | 27 | 720 |
| LC_2024 | HEK293T | 2 | 1,080 |
| LC_2024_2 | HepG2 / ARPE-19 / PC3 | 4 / 4 / 2 | 72 ea |

`value_is_bound` is **empty on every one of them** (0/748 populated). The schema's definition of
`upper` — *"the truth is at or below it"* — is exactly the epistemic status of a reading at or below
blank. Their z-scores span −2.019 to −1.074, i.e. they occupy the bottom of the standardised range
as if they were measurements.

The counts were **computed**. They sit in
`07_qc.json::stage01_file_checks.scale_is_ln_not_minmax.n_below_RLU_floor_1.5`, and
`08_contradictions.csv` C-8 restates them ("1,282 of 5,760 below 1.5"). They then died at the emit
step. That is the F-001 shape verbatim: *"A flag that lives only in a per-folder QC file is not a
flag, it is a memo."*

To be fair to C-8: it asks a *different* question — whether the deposited file was floored and
min-max scaled — and retires it correctly (it was not; the authors' `utilities.py` applies both at
load time, and the deposit is upstream of both, which I confirmed from the value distribution: the
minima are 0.000–1.253, all below the floor, and `exp()` of the smallest distinct values lands on
half-integers, the fingerprint of ln(mean of small replicate counts)). Keeping the unfloored values
is the better choice. The gap is the absence of a **row-level** flag, not a scale error.

**Downstream consequence, concretely.** `LC_2024`'s PC3 arm is 40.8% background. Any
within-experiment rank statistic computed on that arm — which is precisely what Result 6 and the
`CLAUDE.md` §8 within-lab cross-context argument use — is ranking noise over two fifths of its
range. The curation's own D-6 table reports PC3 vs HepG2 ρ = 0.447 with no caveat attached.

**Severity:** changes no value or structure, but mislabels 14.6% of the largest staged block as
measurement. **Confidence: high on the numbers; medium on whether `upper` is the right vocabulary
item** (the schema's examples are colormap floors, and a per-row decision may be better captured by
`confidence` plus a note).

### E-7 — The C10 posture is circular for one of the seven anchors

`PROVENANCE.md` claims *"**7 of 7 anchors are exact matches on the FULL InChIKey**"*. For DOTAP that
match is guaranteed by construction. `make_reference.py:320-322` and the stored note both say the
anchor was composed *"with the stereocentre left unspecified **to match the registry convention**"* —
so in the one case where the paper's drawing and the emitted structure disagree on stereo (D-4, the
bold wedge at the propane C2), the anchor was aligned to the **emission**, not to the drawing. C10
therefore cannot detect the disagreement it exists to detect, and the mechanism is the intent of
`anchor_may_be_a_copy` even though the code path is genuinely different.

For 18PG the anchor honestly follows the drawing — but then the match is between two artefacts
produced from one reading, so C10 passing is not evidence the structure is right.

Effective external coverage is **5 of 7** (MC3, DDAB, DOPE, DSPC, 14PA — where drawing and registry
agree), not 7 of 7. The disagreements are not lost: both are recorded in `08_contradictions.csv` and
`DECISIONS.md`. It is the headline claim that overstates.

I should say what *is* strong here: the AST-based independence gate (`forbidden_imports_found: []`,
`forbidden_calls_found: []`), the genuinely opposite head-first traversal, the CIP-selected rather
than copied backbone tag, and the three seeded defects that each preserve molecular formula and
exact mass to six decimals while changing the InChIKey (C-10's self-inflicted MC3 fragment-direction
bug is a real catch, and the right kind). This folder's anchor machinery is the best I have seen in
the project; the overclaim is one line of prose.

**Severity:** verification integrity; no data error. **Confidence: high on the mechanism.**

### E-8 — No seed covers error class 22

The 13 file seeds cover the one-hot group, the factorial, SI Table 1's levels, error class 40 (two
parameterisations, row-wise), the helper descriptor block, measurement coverage, min-max detection
and the validation workbook. **None seeds a z taken over the wrong block.** That is the class behind
F-016 and the `WT_2026` repair, and this folder makes the block choice twice — two `Experiment_ID`s,
and cell types pooled rather than separated. The reasoning is recorded well in `07_qc.json::zscore`
(including the correct argument that a per-block z of ln equals the per-block z of log10), but it is
argued, not gated. A seed that re-standardises per `(experiment × Model_type)` and requires a
non-zero exit would close it.

**Severity:** coverage gap. **Confidence: high.**

---

## HYPOTHESES I TESTED AND DISCARDED

Recording these because a reviewer's own guesses are claims too, and three of mine were wrong.

- **"`Double_bonds` is an ionizable-lipid descriptor."** No — SI Table 3 is titled "Helper lipid
  chemical parameterization" and the row reads "in the lipid tails"; 5/5 positive control on the
  uncontested helpers. The curation's reading was right.
- **"`LC_2024_2` is a model-selected cherry-pick, so it is not an independent screen."** No — it is a
  complete 6 × 3 × 2 × 2 factorial at N:P 10 (SI Table 6, "Composition formulation parameters used
  in de novo validation"), confirmed row by row from the workbook.
- **"Separate z deletes a real library-level offset."** No — the per-cell deltas range from +4.07 to
  −0.95, so there is no scalar to preserve.
- **"The 72-formulation validation set is undocumented in the paper."** No — the PMC body never says
  "72", but SI Table 6 and Supplementary Figure 3 document it. I nearly reported this before
  checking the SI.
- **"`03_structures_qc.json` undercounts MC3 at 3,081 / 5 experiments."** No — that is the in-vitro
  scoping and it is exact.
- **The first composition join I ran returned 0 overlap everywhere** because my key included a column
  that is `NaN` in `YZ_2022`. The sibling failure the brief warns about, reproduced by me on the
  first attempt. Reported so the number in C3 is read as the *second* answer.

---

## UNVERIFIABLE

- **Which PG the lab actually bought.** No Avanti catalogue number appears in any Mao-lab paper I
  searched — 39375194 (article + SI), 35879315, 40592930, 41124242, 39967850, and the `WT_2026`
  FALCON preprint. A Europe PMC full-text search on `840475` (DOPG), `840503` (SOPG) and `840465`
  (DSPG) returns 58 / 14 / 0 papers, none from this group. Settleable only by the lab's order records
  or an author query, exactly as D-2 says.
- **cLogP as a discriminator.** It would work in principle — the Methods name ChemDraw v22.2.0, and
  DOPG and SOPG differ by one C=C (≈0.2–0.3 logP) — and their printed 18PG value is 8.8188. I have no
  ChemDraw, and RDKit cannot substitute: Crippen reproduces none of their values (DSPC published 6.12
  vs Crippen 12.17; DOTAP 12.515 vs 9.23). Worth one ChemDraw run by whoever adjudicates D-2; it
  would settle it in a minute. Their `cTPSA` is close to RDKit TPSA for the five certain helpers
  (within 0.4) but is insensitive to C=C, so it cannot discriminate.
- **The protocol for the 432 `LC_2024_2` rows.** SI Table 6 and Supplementary Figure 3 document the
  library and its results, but neither states its dose, timepoint or seeding. The emitted
  `Dose_ug_nucleicacid = 0.1` and the 48 h in `curation_note` are inherited from the training-screen
  Methods on those rows. Low severity, but the `curation_note` is identical on all 5,112 rows and
  does not distinguish a stated protocol from an inherited one.
- **D-3 for the HepG2 arm** — not resolvable from the two papers, as the folder says.
- **D-4 (the DOTAP stereocentre)** — I did not attempt it.

---

## FOR THE MERGE DECISION

- If merged, `rows_with_a_caveat` goes **1,092 → 3,036** (1,944 new: 972 `18PG` at
  `confidence = medium`, 972 `DOTAP` at `structure_assumed = 1`). `audit_baseline.json` gates that
  metric in both directions.
- Four new `Model_type` values enter the corpus (`ARPE-19`, `B16F10`, `Neuro2a`, `PC3`). Nothing in
  the pipeline checks cell-line spelling against a controlled list; there is no collision risk here
  (no existing B16/PC3/N2a entry), but it is a new-vocabulary event.
- D-5's consequence, verified arithmetically: `LC_2024`'s z is taken over five of the screen's six
  arms, so `YZ_2022` and `LC_2024` — one physical screen — sit on two scales. Including the HepG2 arm
  would shift `LC_2024`'s centre by about 0.15 z. The folder states this; I confirm it is real.
- A flagging asymmetry worth a sentence, not an error: `18PG`'s 972 rows carry
  `structure_assumed = 0` while `DOTAP`'s carry `1`. The implied rule is "assumed = we departed from
  what the paper drew", which is self-consistent but unstated — and it means a consumer filtering on
  `structure_assumed == 1` to find risky structures **misses the 972 rows carrying the project's most
  contested helper assignment**.

---

## WHAT I DID NOT GET TO

- I did not read `stage01_invert.py` / `stage04_values.py` line by line; I re-derived their outputs
  instead, per the brief.
- I did not verify that `renders/` was extracted at the stated 220–336 dpi — I extracted my own
  images at native resolution and read those.
- I did not check the Cholesterol and DMG-PEG2000 drawings against SI Table 2 (both anchors are
  deliberately empty); I only confirmed both emitted structures are byte-identical to the corpus.
  C-7's DMG-PEG2000 "bracket says 44, the molecule has 45" reasoning is plausible and the InChIKey
  argument is checkable, but I did not re-read the bracket.
- I did not audit the 3,456 excluded physicochemical rows (Size / Zeta / PDI) beyond confirming that
  `05_value_counts.json`'s breakdown sums to 10,368 and that 5,112 + 5,256 = 10,368.
- I did not attempt D-4, and I did not try to obtain ChemDraw to close the cLogP test on D-2.
