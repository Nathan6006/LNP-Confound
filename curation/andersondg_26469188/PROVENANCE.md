# KK_2015 — provenance and anchor posture

**Kauffman KJ, Dorkin JR, Yang JH, Heartlein MW, DeRosa F, Mir FF, Fenton OS,
Anderson DG.** "Optimization of Lipid Nanoparticle Formulations for mRNA Delivery
in Vivo with Fractional Factorial and Definitive Screening Designs."
*Nano Letters* **15**(11):7300–7306, 2015. PMID 26469188.
DOI 10.1021/acs.nanolett.5b02497. First published 2015-10-20.

**`Experiment_ID` = `KK_2015`**, from the **first author** (Kauffman, Kevin J.).
The folder is named `andersondg_26469188` after the corresponding author; a PI
name is not the corpus key. Checked against all **54** existing `Experiment_ID`s
by PMID rather than by year token: no collision, and PMID 26469188 appears
nowhere in the corpus.

**This paper is not a no-original-data record.** `abstractText` is 1,585
characters and the SI carries Methods and Results. `pubType` is `null` in the
Europe PMC record, which is why `pubType` is not the discriminator.

---

## 1. What was held, and why it belongs to this PMID

| | |
|---|---|
| file | `source/nl5b02497_si_001.pdf` |
| sha256 | `47b65dbac03f2bcacf039324dc5efa8a074052803706faf564096d48ee6ae654` |
| bytes | 2,221,826 |
| magic | PDF document, version 1.5, 14 pages |

Three independent attestations, not a file count:

1. the filename stem `nl5b02497` is the suffix of DOI `10.1021/acs.nanolett.5b02497`;
2. SI page 1 reproduces the Europe PMC title and **all eight authors in order**,
   with affiliations and the corresponding address;
3. every one of the 14 pages carries the publisher watermark path
   `pubs.acs.org/nalefd/article-supplement/1257894/pdf/nl5b02497_si_001/`.

The stronger check used by a sibling — byte-identity against the file inside
Europe PMC's own `supplementaryFiles` ZIP — **is not available here**:
`pmcid = null`, `inPMC = N`, `hasSuppl = N`, so no such ZIP exists.

**F-065 guard.** No file was copied into `source/` by this curation. The
directory held exactly one file on arrival and holds exactly one now; it is
asserted by **name and sha256** in `07_qc.json` C0, not by a count or a size.
No glob copy out of the shared session scratchpad was performed at any point.

**The main article is not held and was not fetched.** Europe PMC offers only
`availabilityCode = S` ("Subscription required"). `isOpenAccess = N`,
`inPMC = N`, and `authMan = epmcAuthMan = nihAuthMan = N`, so there is no PMC
author manuscript. `elink dbfrom=pubmed db=pmc` returns **only**
`pubmed_pmc_refs` (citing articles) and **no `pubmed_pmc` self-link** — i.e. no
deposit of any kind. No paywall or bot-check bypass was attempted. Consequently
Figures 1–4, the Methods, the mouse strain, the PEG-lipid identity and the EPO
species are all unavailable.

**Deposit / accession sweep — explicit null result.** 0 hits for every one of 13
patterns (`GSE[0-9]+`, `SRP[0-9]+`, `PRJNA[0-9]+`, `\bSRA\b`, `accession`,
`deposit`, `github`, `zenodo`, `dryad`, `figshare`, `osf.io`,
`data availability`, `source data`) across **18,676 characters** searched:
17,091 of SI text plus 1,585 of abstract. **0 characters of main text**, which is
the one gap in the sweep and is stated rather than papered over. The only URL in
the SI is the NIST Engineering Statistics Handbook, cited as background for
mixed-level fractional factorial designs.

**No author repository.** `find_repo.py` was run for **both** the first author
and the PI; the author-scoped repo queries (`Kauffman lipid nanoparticle`,
`Kauffman LNP`, `Anderson lipid nanoparticle`, `Anderson LNP`) each returned
**n = 0**, and a GitHub code search for the literal `"26469188"` returned 20
hits, none belonging to this paper. The sweep is not silently broken: among its
generic hits it returned `MaoResearchGroup/falcon-lnp-optimization`, which is
the real author repository of the Toh 2026 FALCON paper already curated as
`WT_2026` — a query that can surface a known-good LNP repository is a query that
is running.

---

## 2. Licence — fails §8

Off the **Europe PMC core record**, not the landing page:

| | |
|---|---|
| article | `license = null`, `isOpenAccess = N`, `availabilityCode = S`. **All Rights Reserved, © 2015 American Chemical Society.** |
| author manuscript | none (`authMan = epmcAuthMan = nihAuthMan = N`) |
| repository | none exists |
| funders recorded | Shire; Skolkovo Institute of Science and Technology / MIT — **no NIH grant**, consistent with there being no public-access deposit |
| CC0/CC-BY-compatible | **NO** |

`CLAUDE.md` §8 requires every constituent source to be CC0/CC-BY-compatible
before corpus release. This one is not. See **D-4**.

---

## 3. Value provenance — tier 4, zero pixels

Every value in this curation is a **text character inside the PDF's content
stream**. No page was rendered, no raster was sampled, no axis was calibrated,
no bar top was detected. `curation/lib`'s `axis`, `bars` and `heatmap` modules
are therefore **not imported**, and error classes 20, 23, 27, 28, 32, 33, 35 and
43 **cannot arise here**. This is the outcome `PIPELINE_V2.md` §3 names as the
best available for a paper with native text tables.

| locus | supplies | tier |
|---|---|---|
| SI Table S1 (p2) | serum EPO ng/mL ± SD, n=3, for 36 LNPs + PBS; the full recipe of all 39 formulations; EE / size / PDI | 4 |
| SI p10 DSD table | a **second** parameterisation of A-01…A-14 | 4 |
| SI p11 FFD table | a **second** parameterisation of B-15…B-32 | 4 |
| SI Table S2 (p3) | printed regression statistics, used as external checks | 4 |
| SI Table S3 (p4) | luciferase total flux, 2 formulations | 4 |
| SI Table S4 (p5) | serum FVII %, 2 formulations × 3 doses + PBS | 4 |

`triage.py`'s `table_pages` text heuristic reported only pages 2 and 11. It was
**wrong in the permissive direction**: `pdfplumber.find_tables()` succeeded on
**six** pages (2, 3, 4, 5, 10, 11), and all six were read and verified by hand.
`ruled_tables` correctly reported 0 pages missed.

**One place where geometry was unavoidable, and it is still not pixels.**
Table S4's reading order is not its visual order: the dose cell is vertically
centred, so each formulation's first dose line (`0.1 mg/kg: …`) sorts **above**
its own code label. A line-order parser recovered **6 of 8** value rows and
dropped the two 0.1 mg/kg points silently. Doses are therefore assigned by
**native text y-coordinate** to the nearest phospholipid anchor: anchor
separation 28.08 pt, row pitch ≈9.1 pt, worst assignment 9.40 pt to the correct
anchor against 18.68 pt to the wrong one — a **1.99× margin**. The parser does
not rely on that margin: it asserts each formulation receives exactly the dose
set {0.01, 0.03, 0.1} and raises otherwise. Seed `table_s4_geometry` removes an
anchor token and confirms it refuses rather than falling back.

---

## 4. ANCHOR POSTURE — recorded by hand

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **not run**, as instructed. This section is the
by-hand record.

### 4.1 This paper draws nothing

There is **not one chemical structure drawing** in the 14 SI pages. Per-page
object counts: the eight text pages (1–5, 10, 11, 14) hold **0 image objects**;
the six image pages are Fig S1 (regression diagnostics), Fig S2 (EE-vs-weight-
ratio and size-vs-PEG trends), Fig S3 (a methodology flowchart), Fig S4 (TEM
micrographs) and two pages of JMP screenshots. The main text, which might depict
C12-200, is closed-access and not held.

### 4.2 C12-200 has NO ANCHOR, and that is the load-bearing gap

`reference_structures.csv` carries `C12-200` with an **empty
`reference_SMILES`**, on purpose, so the absence is visible rather than implicit.

The paper gives **no drawing, no molecular formula, no HRMS, no molecular
weight and no synthesis**. Nor is one derivable: the C12-200:mRNA weight ratio
and the molar percentages are independent recipe parameters, so no molecular
mass follows from them. The emitted structure is carried over from the live
corpus — `AC_2025` (498 rows) plus `JW_2024` (20 rows, under the name `other`),
InChIKey `RVHYPUORVDKRTM-UHFFFAOYSA-N`, C70H145N5O5 — via
`corpus_components.smi()`, whose guard confirms the name is unambiguous.
All 36 emitted rows carry `confidence = medium`, `structure_assumed = 1` and a
`curation_note` saying exactly this.

### 4.3 The excipients the paper *names* are anchored instead

The Table S1/S3/S4 footnotes spell the systematic names out:
`DS = 1,2-distearyol-sn-glycero-` (the SI's own misspelling of *distearoyl*),
`DO = 1,2-dioleoyl-sn-glycero-` (Δ9-*cis*), `PC = 3-phosphocholine`,
`PE = 3-phosphoethanolamine`. `make_reference.py` re-derives **DSPC, DOPE and
DOPC** from those names.

**What this anchor is:** a **NAME-CONCORDANCE** check. All three legs of the
comparison — my build, `provenance/reagent_registry.py`, and 24,856 corpus `HL`
rows — are derivations from the *same* systematic name by different hands. It
agrees 3/3 on **full InChIKey**, maximum exact-mass delta **0.00000000 Da**.

**What this anchor is NOT:** an external structural attestation of *this
paper's* chemistry. It would catch a transcription error in any leg — the
briefing records a sibling whose hand-built DSPC came out **a ketone bearing an
ether** at the right mass to six decimals — but it **cannot catch a name that is
wrong for this paper**. It is therefore **not counted as a C10 pass.**

**Cholesterol is not anchored and not hand-written.** The paper gives only the
abbreviation `Chol.` in a column header. `CHL_SMILES` is populated from the
audited registry, and **that registry pin is explicitly not counted as an
anchor** — the same posture a sibling took rather than claiming it. A sibling's
hand-written cholesterol came out a diastereomer with 4 of 8 centres inverted at
identical formula, identical exact mass to six decimals and an identical
skeleton InChIKey block.

**The PEG-lipid is not anchored because it is never named** (§D-3).

### 4.4 How independence was verified — by AST, at run time

Byte-identity is **not** the signal and testing for it is the trap (F-037): a
correct re-derivation canonicalises byte-identically to a correct stored
structure, and my built DSPC/DOPE/DOPC do exactly that. So the **route** is what
is tested.

`make_reference.ast_guard()` parses **its own module** at run time and inspects
**`Call`, `Import` and `ImportFrom` nodes only**, so the ban list cannot match
itself — a substring scan over the source would flag the ban list inside its own
docstring, which is how one sibling's guard failed on its own module. Measured:
**8 import nodes and 168 call nodes examined, 0 findings**, against a ban list
of `corpus_components`, `reagent_registry`, `reagents`, `pandas`, the stage
modules, and the calls `smi`, `load`, `read_csv`.

A **template-provenance guard** additionally requires that no string literal of
≥25 characters in the module be a substring of any emitted structure:
**54 literals examined, 0 hits.**

**Both guards were made to fail on purpose, in a child process**
(`make_reference.py --seed-copied-template`): a seeded copy that imports
`reagent_registry` and carries a 40-character literal lifted from an emitted
SMILES produces `ast_guard` findings `[('import', 'reagent_registry', 51)]` and
**2** template hits. Both fired.

### 4.5 The route, and what it asserts

Head-first `assemble.attach` with atom-mapped attachment points only — phosphate
head onto the glycerol sn-3 oxygen, then the two acyl carbons — which is the
**opposite** order from writing an acyl-first lipid SMILES. No string
concatenation, no `RWMol`, no `RemoveAtom`, no "atom index 0 of the fragment".

The two chain lengths differ across the anchor set (18:0 saturated for DSPC,
18:1 Δ9-*cis* for DOPE/DOPC), and the head groups differ (quaternary vs primary
ammonium), so the anchors are not three copies of one shape.

* **Stereo assertions are NOT vacuous.** Each lipid has exactly **1** real
  stereocentre (glycerol C2), asserted **(R)** — derived from the stereospecific-
  numbering convention, not copied from any stored structure. Seed
  `make_reference_stereo` flips the backbone tag and the assertion fires,
  reporting **(S)**.
* **Double-bond geometry is asserted explicitly**, never left bare: 2 C=C, 2 *Z*,
  **0 unassigned** for the dioleoyl lipids; 0 C=C for DSPC.
* **The gate is proven to test connectivity, not mass.** Seed
  `make_reference_regioisomer` builds 1,3-diacyl-glycero-2-phospho DSPC:
  **identical formula C44H88NO8P, exact-mass delta 0.000000 Da**, different
  InChIKey. The route **refuses** it, via `assemble.check_stereo_preserved` —
  because C2 stops being stereogenic once both CH₂ arms become identical, which
  is the correct chemistry and the correct refusal.

### 4.6 MEASURED SENSITIVITY LIMIT — stated honestly

`make_reference`'s own assertions **cannot see double-bond position.** Seed
`three_way_delta6` moves the oleoyl double bond from Δ9 to Δ6: the C2 = (R)
assertion, the C=C count, the *Z* count and the unassigned-geometry check **all
still pass**, at identical formula and **exact-mass delta 0.000000 Da**. Only
the three-way InChIKey comparison catches it — and it does
(`KLJBLQFSKZACGJ-…` vs `SNKAWJBJQDLSFF-…`). That is the real limit of this
curation's structural gates, and it is why the three-way check exists rather
than being treated as a formality.

---

## 5. Token-pattern test — `C12-200`

`C12-200` is `C<tail length>-<amine catalogue index>`: `C12` is the tail (a
2-hydroxydodecyl arm from 1,2-epoxydodecane) and `200` is the index of the
polyamine core in this lab's amine catalogue. **`KK_2015` writes the two tokens
in the opposite order from `KL_2010`'s `200C12`.**

**Which token is which is proven structurally, not by substring.** Within a
homologous series that holds the amine fixed and walks the tail length,
`dC/d(tail)` **equals the number of arms**. Measured:

| series | experiment | N | `dC/d(tail)` | arms implied |
|---|---|---|---|---|
| `A200_C10 … A200_C14` | `KW_2014` | 5 | 5.0, 5.0, 5.0, 5.0 (constant) | **5** |
| `200C8 … 200C14` | `KL_2010` | 3 | 1.0, 1.0, 1.0 (constant) | **1** |

So `200` names the core (it is held fixed while the C-number varies) and the
C-number names the tail. The numbering is **not** contiguous — it is a sparse
catalogue index, not a local 1..N index (F-055).

**Naive-substring damage, measured:** `IL_name contains "200"` returns
**526 rows across 24 names and 24 distinct structures**, spanning at least three
incompatible chemotypes. Never used; every join here is on **InChIKey**.

---

## 6. Corpus overlap

**Positive control first** (F-070). cKK-E12: name join **216** rows → InChIKey
join **241** rows, recovering all three names `CKK-E12` / `other` / `MD-1`
across `AC_2025` (216), `JW_2024` (19) and `BL_2023` (6 in vivo). This
reproduces the documented figures exactly, so a zero from the same machinery
means something. A second control, `DLin-MC3-DMA`, returns both of its spellings.

**C12-200 is already heavily represented — and only in vitro.**

| | |
|---|---|
| name join | **498** rows (`AC_2025`) |
| **InChIKey join** | **518** rows — 20 more, under the name `other` in `JW_2024`, the same shape as the cKK-E12 control |
| in vitro | 518 |
| **in vivo** | **0** |

**No new in-vivo shared lipid.** The in-vivo pool's shared structures are
confirmed to be exactly **two**: `DLin-MC3-DMA` (under two spellings) and
`SM-102`. C12-200 appears in **0** in-vivo experiments today, so admitting
`KK_2015` takes it from 0 to 1 — and a structure needs ≥2 in-vivo experiments to
be shared. The in-vivo shared set would stay at exactly MC3 and SM-102. What
`KK_2015` *would* create is the **first in-vivo appearance of a 518-row
in-vitro structure**, a vitro↔vivo bridge, which is a different and smaller
thing and is reported as such.

**Substructure sweep — labelled a SUPERSET.** Positive control: C12-200 itself
matches at every tier, asserted in code so the sweep raises if it does not.

| tier | meaning | structures | rows |
|---|---|---|---|
| loose | any amine bearing a 2-hydroxyalkyl arm (the epoxide-ring-opening chemotype) | 1,535 / 19,919 | 2,386 |
| strict | the amine-200 pentaamine core **and** a hydroxyalkyl arm | 1 | 518 |
| exact | full InChIKey identity | 1 | 518 |

The loose tier is a **chemotype superset, not identity**. Only the exact tier is
identity.

**Helper lipids.** `DSPC` 7,038 corpus rows / 1 SMILES; `DOPE` 17,737 / 1;
`DOPC` 81 / 1 — all three three-way concordant with my builds and the registry.
**`DSPE` has 0 corpus rows**, and is also never measured in this paper: all
three DSPE designs failed to formulate.

**F-061 is reported, not repaired.** `data/` is untouched. See **D-9**.

---

## 7. Benchmarks — does the benchmark top the screen?

**No.** The benchmark is `Original`, the C12-200 formulation previously optimized
for siRNA (50 / 10 / 38.5 / 1.5 mol%, DSPC, 5:1 weight ratio) — the field-standard
starting point this paper sets out to beat.

| | EPO (ng/mL) | rank of 36 |
|---|---|---|
| `Original` (benchmark) | 962 | **20** |
| `C-35` (declared optimum) | 7,065 | **4** |
| `C-36` | 7,548 | 1 |
| `B-26` | 7,485 | 2 |
| `C-37` | 7,268 | 3 |

**19 of 36** formulations beat the benchmark. And the formulation the paper calls
"the optimized formulation" — C-35, the one carried into Tables S3/S4 and into
the headline claim — is **not** the screen maximum. That is not an error: C-35 is
the Library B *regression* optimum and Library C is the confirmation ladder run
around it, so "optimized" means model-optimal. It is recorded because a winner
selected into a second experiment can otherwise look tautologically best, and
here the arithmetic says otherwise.

The abstract's **"7-fold"** reproduces exactly on C-35: 7065 / 962 = **7.344×**.
On the screen maximum it would be 7.846×.

---

## 8. What this contributes, plainly

**One ionizable lipid.** Thirty-six individually-measured in-vivo arms, three
distinct helper lipids, eight weight-ratio levels — and **one** ionizable lipid,
C12-200, on every single row. For a corpus keyed on ionizable-lipid chemistry
the contribution is **1**, however many rows exist. The row count was not padded:
the 2 luciferase and 6 FVII values are recorded as admissible in `05_values.csv`
and deliberately **not** emitted (D-8).
