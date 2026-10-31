# Curation QC protocol

Written before any extraction, 2026-10-27. The point of this file is that the pipeline is
**file-based and resumable**: every stage writes its output to disk, so work survives a
context reset and can be audited without re-reading the PDFs.

## Jobs

| dir | paper | role | ground truth |
|---|---|---|---|
`liu_iphos` | Liu et al., *Nat Mater* 2021, iPhos, PMID 33542471 | **validation** | `SL_2021` (LNPDB, 572 ILs) **and** `Liu_iphos` (LION, 572 ILs) |
`amino_lipid_fg_2026` | Gong et al., *Nat Mater* 2026, PMID 41922838 | **validation** | `FG_2026` (960 ILs, continuous values) |
`siegwart_42634137` | Tian et al., *Adv Mater* 2026, PMID 42634137 | **new curation** | none |

Validation runs are **blind**: extract from the PDF only, write the result, and only then
load the curated rows to score. Never look at ground truth first — a scored comparison is
worthless if the answer was visible while extracting.

## Stages (each writes a file; re-runnable)

0. `00_pages/` — `pdftotext -layout` per page, plus `pdftoppm -r 200` PNGs of pages that
   hold structures or data figures. Cheap, done once.
1. `01_scheme.json` — the reaction scheme, read from the paper: building-block lists
   (heads, tails, linkers) with names, the reaction, and **degree of substitution**
   (how many tails attach per head). This is the single highest-leverage thing to get
   right: an error here is systematic across every lipid in the source.
2. `02_blocks.csv` — each building block with a SMILES, plus how it was obtained
   (`name` = from a chemical name in text, highest confidence; `drawn` = read from a
   structure diagram; `catalog` = vendor/CAS lookup).
3. `03_enumerated.csv` — the combinatorial product set, built by RDKit from 01+02.
4. `04_values.csv` — per-lipid activity, with `extraction_method` of `table` (SI table),
   `heatmap_bin` (discrete colour class) or `digitized` (measured from a plot), plus
   page/figure provenance.
5. `05_formulations.csv` — molar ratios, helper/chol/PEG identities, IL:mRNA ratios, and
   the mol→weight arithmetic with the MWs used recorded inline.
6. `06_merged.csv` — one row per LNP in the 67-column LNPDB schema plus provenance
   columns (stripped on final merge).
7. `07_qc.json` — every check below, pass/fail per row.

## Checks (all automated; `07_qc.json`)

Ordered by how much they actually catch.

- **C1 exact mass vs reported HRMS.** The strongest check. One CH2 = 14.0157 Da, so a
  ±1 carbon error — my worst failure mode on long chains — is unmissable at ppm
  tolerance. Any lipid whose HRMS is reported in the SI gets independently confirmed.
- **C2 benchmark lipids as positive controls.** These papers formulate MC3 / SM-102 /
  ALC-0315 / 306Oi10 as comparators, and the corpus already holds their canonical
  SMILES. Reproducing a benchmark exactly validates my reading of *that paper's* drawing
  conventions before its novel structures are trusted.
- **C3 enumeration count == stated library size.** If the paper says 572 and I produce
  560 or 1144, the degree of substitution is wrong. Systematic, so high value.
- **C4 no duplicate canonical SMILES within a source.** Duplicates mean regiochemistry
  or symmetry was mishandled.
- **C5 RDKit parse + valence sanity** on every SMILES; N count > 0 for an ionizable lipid.
- **C6 molar ratios sum to ~100**, and recomputed MW matches any stated MW.
- **C7 render-and-compare.** RDKit draws my SMILES; compare against the paper's diagram
  side by side. Catches gross topology errors. Will NOT reliably catch ±1 carbon — that
  is C1's job, and C7 must never be used as a substitute for it.
- **C8 overlap against the 19,827-lipid union** of all of the user's layers
  (`new_data/LNPDB_vitro_del_processed.csv`, `new_data/lnpcd_tox_processed.csv`,
  `confound/data/*.csv`). Catches accidental re-curation of something already held. For
  the two validation papers a HIGH overlap is the expected, desired result.
- **C9 value distribution sanity**: bin counts for heatmaps, range checks, and for a
  discretized source the number of distinct values must equal the number of legend bins.

## Scoring the validation runs

Against each ground truth, report:
- canonical-SMILES exact-match rate (recall and precision vs the curated set)
- Spearman and MAE on the activity value, over the SMILES intersection
- for binned sources, the confusion matrix across bins

**The reference bar is not 1.0.** LNPDB and LION independently curated iPhos and agree at
only **Spearman 0.897** on identical lipids from the same figure. That is the human-vs-human
ceiling, so my score must be read against 0.897, not against perfection. Report both.

## Rules

- Low-confidence rows are **held back**, not emitted silently. `confidence` and the failing
  check go in the output so the user's manual review is targeted.
- Every value carries `source_pdf`, `page`, `figure_or_table`, `extraction_method`.
- No value is ever invented to fill a gap. Missing stays missing.
- Where per-lipid values exist only inside a heatmap at a resolution that cannot be read
  reliably, say so and stop rather than produce plausible numbers.
