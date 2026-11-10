# PROVENANCE — goulart_41138890

Curated 2026-11-10. Read-only outside this folder; nothing in `data/` was
touched and no corpus-level check was run (`provenance/structure_check.py`,
`check_anchor_independence.py` and `check_anchor_coverage.py` write tracked
files other agents read, so the anchor posture is recorded here by hand
instead, as the run's brief requires).

## 1. Identity

| | |
|---|---|
| PMID | **41138890** |
| DOI | 10.1016/j.ijpharm.2025.126278 |
| PII | S0378517325011159 |
| Title | *High-Throughput in vivo and in vitro screening of lipid nanoparticles for nucleic acid delivery to the brain* |
| Journal | Int J Pharm 686:126278, first published 2025-10-23 |
| First author | **Marco Túllio Rodrigues Alves** (Crossref `given`/`family`, `sequence: first`) |
| Corresponding | Pedro Pires Goulart Guimarães, UFMG — `ppiresgo@reitoria.ufmg.br` |
| Proposed `Experiment_ID` | **`MR_2025`** |

`MR_2025` derives from the **first author**, not the PI the folder is named
for. Checked against all **54** corpus `Experiment_ID`s **by PMID**: the ID is
free and PMID 41138890 appears on no corpus row.

The no-original-data discriminator is clear: `abstractText` is **1,740
characters** with Methods and Results, so this is an ordinary research article.

## 2. Licence — fails CLAUDE.md §8

Read off the **Europe PMC REST core record**, not the publisher landing page:

```
isOpenAccess = N      license = null      subsetName = null
inEPMC = N            inPMC = N           hasSuppl = N
fullTextUrlList = [ { availability: "Subscription required" } ]
```

Crossref lists only the Elsevier **TDM user licence** and STM article-sharing
policies 017 / 037 / 012 / 029 / 004. **No Creative Commons licence of any
kind.** There is no repository, so no second licence to check.

## 3. Acquisition

| what | route | verdict |
|---|---|---|
| SI `mmc1.docx` | `ars.els-cdn.com/content/image/1-s2.0-S0378517325011159-mmc1.docx`, anonymous | 7,825,773 B, magic `504b0304`, sha256 `b4f0839…a977ab` |
| `mmc2`–`mmc5` | same route × `.pdf/.xlsx/.docx/.zip/.csv/.txt` | **30 probes, 29 × HTTP 404** with a 189-byte XML body. Only `mmc1.docx` exists. |
| Figures 1–6 | `…-gr<N>_lrg.jpg` | 6 × HTTP 200 at 2,247–2,647 × 3,440–3,603 px; `gr7`/`gr8` 404 |
| Article PDF | — | paywalled; **not** obtained, no bypass attempted |
| PMC deposit | `elink dbfrom=pubmed db=pmc` | only `pubmed_pmc_refs` (2 citing PMCIDs). **No `pubmed_pmc` self-link → no deposit**, so no PMCID and no `supplementaryFiles` route. |
| Author repo | `find_repo.py` ×2 (first author and PI), each with `--doi` and `--pmid` | the DOI/PMID code hits are a digit-string coincidence in an HPC cache-benchmark repo. **None.** |
| Accessions | regex over **1,002,901 characters** (docx paragraphs 6,382 + EMF text 19,724 + raw `word/document.xml` 976,795) for `GSE\d+`/SRP/PRJNA/SRX/SRR/SRA/accession/deposit/github/zenodo/dryad/figshare/osf.io/mendeley | **0 hits** |

**F-065 posture.** Figures were copied into `source/figs/` **by explicit name,
one `cp` per file, never a glob**, and every sha256 was recomputed after the
copy. The SI is attested against content only this paper can produce: its
`word/document.xml` carries this title, the 13-author list and
`ppiresgo@reitoria.ufmg.br`.

## 4. The SI is a VECTOR container, and a tag count misreads it

`word/document.xml` holds exactly **one** `<w:tbl>` — Table S4, 156 rows of
pairwise statistics. Tables **S1, S2 and S3 are EMF images**. A `<w:tbl>` count
therefore reports one table for a document holding four.

But EMF is a Windows **metafile**, not a bitmap. Walking the record stream:

| image | table | `EMR_EXTTEXTOUTW` records |
|---|---|---|
| `image15.emf` | Table S1 | **694** |
| `image16.emf` | Table S2 | 1,571 |
| `image17.emf` | Table S3 | **534** |

One Unicode glyph per record with an exact device-space reference point. So
every digit in `02_tableS1.csv` / `02_tableS3.csv` is a **transcription**, not
a pixel measurement — tier-4 precision reached from a tier-2 container. The
only inference is the grouping, and its threshold is measured: the widest
within-cell glyph advance is 147 (S1) / 146 (S3) device units and the narrowest
between-column gap is 387 / 306, ratios **2.633** and **2.096**.

> Getting the EMF record numbering wrong decoded 41 `EXTCREATEFONTINDIRECTW`
> records as text and printed 10.5 MB of binary. Above `0x51` the numbering is
> `STRETCHDIBITS=81, EXTCREATEFONTINDIRECTW=82, EXTTEXTOUTA=83,
> EXTTEXTOUTW=84, POLYBEZIER16=85`.

## 5. Anchor posture (C10 / C16 / C18, recorded by hand)

`reference_structures.csv` has the exact 7-column header and **5 rows**:
**1 with a structure, 4 empty on purpose**, 1 carrying an `expected_mismatch`.

**C12-200 — a real anchor, exact full-InChIKey match.**
`make_reference.py` assembles it from the Figure 1A drawing alone:
piperazine core written from its drawn connectivity, then the five N–H sites
alkylated one at a time through `assemble.attach`, each target chosen by
SMARTS and never by atom index, with the tail written as `"C" * 10` so no long
literal exists in the file at all. Result
`RVHYPUORVDKRTM-UHFFFAOYSA-N`, `C70H145N5O5`, exact mass **1136.124573** — the
corpus structure on **498** `AC_2025` rows, whose stored `Molecular.Weight` of
1136.96 agrees. Every tail is **printed as the text "C10H21"** in the figure,
so chain length is read rather than counted; what the anchor really tests is
the **connectivity** — which nitrogen carries how many arms, and on which side
of the piperazine the three-arm branch sits.

**Independence is executed, not asserted.** Two guards, both made to fail on
purpose in a child process with a clean control that comes back clean:

| guard | how | seeded | fires |
|---|---|---|---|
| AST independence | walks `make_reference.py`'s own AST, inspecting **`Import`/`ImportFrom`/`Call`/`Attribute` nodes only**, so the ban list cannot match itself | `import corpus_components` injected into a temp copy | **yes** |
| template provenance | no string literal ≥ 25 chars may be a substring of any structure the anchor is compared against | a 40-char slice of the corpus SMILES injected | **yes** |

**The gate tests connectivity, not composition.** Seeding a Markovnikov
epoxide opening — `N-CH(C10H21)-CH2OH` instead of `N-CH2-CH(OH)-C10H21` —
gives formula `C70H145N5O5` **identical**, exact mass **Δ 0.000000 Da**,
skeleton InChIKey block **different**, full InChIKey different. Only the
InChIKey sees it.

**Where there is NO anchor, and the guards are therefore VACUOUS, not passing:**

- **SM-102, ALC-0315** — not drawn anywhere; the paper gives the names only, in
  Fig 4G and Fig S8. Structures come from `provenance/reagent_registry.py`,
  which is a **name-derived pin, not an anchor**. Empty `reference_SMILES` on
  purpose.
- **Cholesterol, DOPE, DSPC, DOTAP** — drawn in Fig 1A, but pinned from the
  same registry, and comparing a registry pin against itself is **declared not
  independent** rather than counted. Not hand-written: an agent's hand-written
  cholesterol came out a diastereomer with 4 of 8 centres inverted at identical
  formula, identical exact mass and an identical skeleton InChIKey block.
- **C14PEG2k** — drawn in full but carries `expected_mismatch =
  head_group_class_differs_from_corpus_name`; see D-4.

## 6. What was measured, with numbers

See `07_qc.json`. The three that carry the most weight:

- **An independent tie between the main text and the SI.** Fig 2G's three
  legends print a pKa for nine formulations; Table S3 prints one for eleven.
  Different files, different renderers, different extraction paths, and
  different decimal separators (`7.47` vs `7,47`). **9/9 exact, 0 disagree.**
  One misread digit in the EMF transcription would break it.
- **A round number the authors DREW and never printed.** The `0` tick must land
  on the drawn axis baseline: measured offset **0.101 px = 1.14 axis units**.
  It earns its place because it is the *only* check that can see a uniform
  origin error — the `ticks_uniform_shift` seed moves every tick by +4 px, so
  `calibrate` and the LOO gate both still pass perfectly at 0.279 px while this
  offset moves to 4.101 px.
- **The factorial, proved rather than claimed.** 3 helper lipids × 6 molar
  contents = 18 = the row count of Table S1; the ionizable-lipid molar ratio is
  **35.0 on all 29 SI rows** and the PEG **2.5** on all 29; mol% sums exactly
  100.0 on **18/18** (S1) and **11/11** (S3); Table S3 is a composition-exact
  subset of Table S1 on **11/11** with 0 disagreements.

## 7. Reproducibility

`stage09_seeds.py` runs **8 seeds, 8 credited**, one child process per seed,
every verdict read from `returncode` and never from stdout or a pipeline.
Every seeded run writes to its **own temp directory** (F-069); afterwards both
tracked artifacts are regenerated into a fresh directory and compared by
sha256 — **byte-identical across two clean runs**.

Run order:

```
stage02_si_tables.py → stage03_crosscheck.py → stage04_fig4g.py
→ make_reference.py → stage05_anchor_check.py → stage06_emit.py
→ stage09_seeds.py → stage07_qc.py
```

## 8. Merge status

**BLOCKED, deliberately.** `DECISIONS.md` carries **five `STATUS: OPEN`**
entries (D-1, D-2, D-3, D-4, D-6) and three settled ones. The governing one is
D-1: the paper contributes **three** ionizable lipids and **zero** of them are
new, and all three are the field's universal controls.
