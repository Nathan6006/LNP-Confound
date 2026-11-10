# pabst_37602460 -- provenance

**Status: EXTRACTED 2026-11-10. NOTHING MERGED. `DECISIONS.md` D-1 is OPEN.**

| | |
|---|---|
| Paper | Jennings, Ascheric, Semeraro, Lohner, Malanovic, Pabst. *Combinatorial Screening of Cationic Lipidoids Reveals How Molecular Conformation Affects Membrane-Targeting Antimicrobial Activity.* ACS Appl. Mater. Interfaces 2023 |
| DOI | 10.1021/acsami.3c05481 |
| PMID | 37602460 (PMC10472336) -- indexed; no unindexed claim is made anywhere in this folder |
| Licence | CC-BY 4.0, read from the PMC record's `<permissions>` element |
| Experiment_ID | **none assigned.** No `Experiment_ID` is reserved, because whether this paper belongs in the corpus at all is `DECISIONS.md` D-1 |
| Rows emitted into `data/` | **0** |
| Rows previewed in the 79-column quality schema | 104 (`07_emit_preview.csv`) |
| Value source | see `00_value_sources.json` |

## The one thing to read first

**This is not a transfection paper.** The brief that commissioned the curation
gave the title as "...Reveals How Molecular *Structure* Affects..." with the
rest elided. The published title ends "...How Molecular **Conformation** Affects
**Membrane-Targeting Antimicrobial Activity**", and the measured quantity is
growth inhibition of *Bacillus subtilis* in Luria Broth. There is no
nucleic-acid cargo, no helper lipid, no cholesterol, no PEG-lipid, no
nanoparticle and no mammalian cell line anywhere in the assay. The chemistry is
squarely in family -- polyamine + acrylate lipidoids, and eleven of the 104
members are already in `data/data_vitro.csv` -- but the endpoint is not.

## Where the values come from

Every value in this paper is text. Nothing was digitised.

| what | how many | tier | locus |
|---|---|---|---|
| antimicrobial point score | 104 | `4_pdf_text_table` | SI Table S3, page S10, by pdfplumber word coordinate; confirmed cell-by-cell against a 600 dpi render |
| MIC / uM | 27 | `4_pdf_text_table` | main-text Table 1, from the PMC JATS `<table>` markup |
| therapeutic index | 27 | `4_pdf_text_table` | main-text Table 2, same route |
| HC50 / uM | 27 | derived | TI x MIC, which lands on the paper's own discrete hemolysis ladder 27/27 |

Tier 1 does not exist (ACS publishes no Source Data; the article's only
supplement is the one SI PDF). Tier 3 does not exist (`find_repo.py` returned
`repo_hits: []` on 2026-11-10 and the paper names no repository). Two panels are
recorded as `5_pixels` and neither supplies a value: Figure 1A and SI Figure S1
are read for **structures**, because the paper states no SMILES, InChI, MOL
block, CAS number or chemical name for any building block anywhere.

## The paper's own cross-checks

Fifteen checks, every one of them with its own seeded defect confirmed to move
its metric. Numbers in `07_qc.json`.

- tail acrylate MW vs SI Table S4: **13/13** within 0.05 Da, and the three
  constitutional-isomer MW ties the naming scheme predicts (6B/8, 10B'/12,
  13B/14) are all present in the printed table
- tail name token vs the length measured on the arm's own graph: **13/13**
- headgroup N/MW vs SI Table S5: **6/8** (the two misses are each other's
  values -- see below)
- tails per headgroup vs the article's prose: **8/8**
- SI Table S3 score == 6 **iff** main-text Table 1 MIC <= 25 uM: **27/27**,
  across two independently produced tables, an identity the authors never state
- all eight bounded therapeutic indices equal **400/MIC** exactly -- a constant
  that appears nowhere in the paper
- nine printed statistics reproduced: **9/9**
- three C10 anchors, built from explicit atom/bond tables with no SMILES string
  anywhere in the generator, all matching the enumeration by **full InChIKey**

**Passing all of them is necessary and not sufficient.** `QL_2026` reproduced
every printed hit rate on bars up to 28% too high. The load-bearing evidence
here is deliberately not a recomputed marginal: it is the 27/27 cross-table
identity, the unprinted 400, and the eleven members whose parent structure is
already in the corpus under `KW_2014`, `AA_2008` and `NC_2024` with a matching
full InChIKey.

## What the paper gets wrong

Nine entries in `08_contradictions.csv`. Three matter:

1. **SI Table S3's "Tail total" column contradicts its own cells** on 11 of 13
   rows. The cells are kept; they are corroborated by 7 of 8 printed column
   totals and by the 27/27 cross-table identity.
2. **SI Table S5's charge densities for 3NB and 4NC are transposed** -- each
   printed value is the other structure's computed value to four decimals. Four
   independent lines agree: the arithmetic, the paper's own stated trend (whose
   sign is wrong with the printed values and right with the computed ones), SI
   Figure S1 panel (E) drawing the cyclam spacers explicitly, and the corpus
   match to `KW_2014`'s `A138_C10`/`A138_C12`.
3. **SI Figure S1 panel (A) draws the 3N headgroup as a nine-membered triaza
   ring** carrying five arms, which is impossible -- a ring nitrogen has one
   N-H, so that core admits at most three. Figure 1A's linear diethylenetriamine
   is used. Recorded because connectivity is precisely what formula and mass
   cannot see.

## Corpus overlap -- a merge-level fact, reported and not acted on

Eleven of the 104 **parent** (unmethylated) structures already exist in
`data/data_vitro.csv` by full InChIKey, and one of those also in
`data/data_vivo.csv`:

| member | already in the corpus as |
|---|---|
| 2N14, 2N16 | `KW_2014` `A157_C10`, `A157_C12` |
| 3N13B | `NC_2024` `A113_Oi10` |
| 3N14, 3N16 | `KW_2014` `A114_C10`/`A114_C12` **and** `NC_2024` `A113_O10`/`A113_O12` |
| 3NB14, 3NB16 | `KW_2014` `A110_C10`, `A110_C12` (3NB14 also in vivo) |
| 4N14, 4N16 | `AA_2008` `98_O10`/`98_O12` **and** `KW_2014` `A98_C10`/`A98_C12` |
| 4NC14, 4NC16 | `KW_2014` `A138_C10`, `A138_C12` |

Zero of the 104 **methylated** cations match anything in the corpus, by full or
skeleton InChIKey. Zero `IL_name` strings collide. Two things worth noting for
whoever adjudicates D-2: the overlap is with the parent, and the assay was run
on the quaternised salt; and `3N14` resolves to **two** corpus amine labels
(`A114` in `KW_2014`, `A113` in `NC_2024`), which is a pre-existing corpus
condition, not something this paper introduces.

## What could not be verified

- **The degree of methylation.** Figure 1A draws every nitrogen methylated and
  positively charged, and the whole library was dosed in that form, but the
  paper reports no NMR integral for the N-CH3, no elemental analysis and no
  titration. SI Table S1's `IN/IA` ratios quantify tail substitution only.
- **Three tail isomer conventions** -- isodecyl (13B), tricyclodecanyl (9C),
  isobornyl (8CB) -- are commercial isomer mixtures; 24 of 104 members carry
  `structure_assumed = 1` and `confidence = medium` on the row for this reason.
- **Why the printed Tail total column is wrong.** One residual is 7, which no
  single dropped score (0-6) can explain, so it is not a missing ninth column.
