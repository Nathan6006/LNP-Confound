# PROVENANCE — Sago et al., PNAS 2018, PMID 30275336

## Paper identity, verified from the Europe PMC REST `core` record

One call: `search?query=EXT_ID:30275336 AND SRC:MED&resultType=core`.

| field | value |
|---|---|
| full published title | High-throughput in vivo screen of functional mRNA delivery identifies nanoparticles for endothelial cell gene editing. |
| journal | Proc Natl Acad Sci U S A **115**(42), 2018 |
| PMCID / DOI | PMC6196543 / 10.1073/pnas.1811276115 |
| `abstractText` length | **957** |
| Methods / Results present | yes / yes |
| verdict | ORIGINAL RESEARCH ARTICLE |

The no-original-data discriminator is `abstractText` length 0 **plus** no
Methods/Results — not `pubType`. Both tests pass here, so the paper carries
original data.

**THE ENDPOINT.** Pooled, competitive **normalised DNA-barcode delivery**
(`Normalized DNA Delivery (%)`) read by Illumina deep sequencing of 56-nt DNA
barcodes co-formulated with Cre mRNA and recovered from FACS-sorted tdTomato+
cells of Ai14 Cre-reporter mice. Secondary, individually measured endpoints for
the two lead LNPs: ICAM-2 protein silencing (% MFI, flow cytometry), ICAM-2
indel rate (%, TIDE, 2 loci), and % tdTomato+ cells.

## Licences

| artifact | licence | CC0/CC-BY compatible? |
|---|---|---|
| article + SI Appendix | **CC BY-NC-ND 4.0** (PMC JATS `<permissions>`, verbatim) | **NO** — NC and ND both fail |
| author repository | none exists | n/a |

`find_repo.py` on the PMID and on Dahlman repo-name queries returned 0 relevant
hits. The Methods name "a custom Python-based tool" and "an R script" and
deposit neither; no GEO/SRA accession appears anywhere. See `CLAUDE.md` §8 and
`DECISIONS.md` D-5.

## Acquisition

* `source/pnas.1811276115.sapp.pdf` — pre-existing in the folder,
  sha256 `32d3e9a6…4e9e8`, 27,953,612 bytes, 26 pages.
* `source/PMC6196543.jats.xml` — `pmc_fetch.py`, after NCBI `elink` returned a
  real `pubmed_pmc` **self**-link (`pubmed_pmc -> 6196543`), not merely
  `pubmed_pmc_refs`.
* Figures 1–4 as JPEG rasters, for scope confirmation only, never as a value
  source. **One magic-bytes catch:** the first `fig04.jpg` fetch returned a ZIP
  (the `supplementaryFiles` response), not a JPEG. It was detected by its magic
  bytes and re-fetched. A size floor would have passed it.

## Value tiers

No digitiser was written. **Tier 5 (pixels) was not used for any value.**
Everything in `05_values.csv` is a printed numeral from a PDF text cell or from
body prose; the formulation tables are tier 4 native `pdfplumber` grids. See
`00_value_sources.json`, which was written **first**.

## C10 ANCHOR POSTURE — recorded by hand

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **NOT run**: they write tracked files other
agents read. The posture is therefore stated here and is reproducible from
`make_reference.py`, which is self-contained.

**Ionizable-lipid anchors: 8 rows, every `reference_SMILES` EMPTY, on purpose.**
The paper draws no complete ionizable-lipid structure — only a single
repeat-unit cartoon with the backbone truncated to two stubs. `CLAUDE.md` §5:
store an empty `reference_SMILES` with a note so the *absence* of an anchor is
visible. There is also nothing for an anchor to test, since this curation emits
zero structures.

**Excipient anchors: 2, built and checked.** Cholesterol and DOPE are drawn in
full in SI Fig S2b. Both were written out by hand from the drawing and compared
by **FULL InChIKey** against what the corpus holds:

| anchor | full InChIKey | matches corpus |
|---|---|---|
| Cholesterol | `HVYWMOMLDIMFJA-DPAQBDIFSA-N` | yes |
| DOPE | `MWRBNPKJOOWZPW-NYVOMTAGSA-N` | yes |

**Independence is verified by PROVENANCE, not by byte-identity.**
`make_reference.py` walks **its own AST at run time** and asserts it imports no
assembler and calls none of `RWMol` / `AddBond` / `RemoveAtom` / `CombineMols` /
`esterify` / `attach` / `ReplaceSubstructs` / `EditableMol`. It walks the AST
rather than substring-scanning the file, because a substring scan would flag
this module's own docstring, which names every one of those functions — that is
how a sibling's first attempt failed.

Byte-identity is **not** the signal (F-037). Both anchors canonicalise to a
string byte-identical to the corpus's, which is the *expected* outcome for a
correct re-derivation and proves nothing on its own. The cholesterol SMILES was
written in a different atom order (`C[C@H](CCCC(C)C)[C@H]1…`) and RDKit
canonicalised it onto the corpus string.

**The gate was proved to fire, with a seed of identical formula and identical
exact mass.** `seed_epimer` inverts C3 to give epicholesterol — F-013 exactly,
the defect that put 3α-cholesterol on 537 `ZT_2026` rows:

| | control | `seed_epimer` |
|---|---|---|
| formula | C27H46O | C27H46O — **unchanged** |
| exact mass | 386.3549 | 386.3549 — **unchanged** |
| InChIKey skeleton block | `HVYWMOMLDIMFJA` | `HVYWMOMLDIMFJA` — **unchanged** |
| **full InChIKey** | `…-DPAQBDIFSA-N` | `…-VEIPTCAHSA-N` — **FIRES** |

`seed_banned_call` flips the AST gate to `CLEAN = false`.

## Seeded checks — the agreement count must MOVE

| check | control | seed | agreement moved? |
|---|---|---|---|
| readout reconstruction | 6/6, 18/18, 3/3, 18/18, residual 0.0 | `wrong_raw_count` | yes → 5/6, 15/18, 15/18 |
| | | `drop_input_normalisation` | yes → step3 1/18 |
| | | `no_renormalisation` | yes → step3 0/18, residual 213.37 |
| chain length | 14 carbons | `drop_a_real_bond` | yes → 13 |
| | | `inject_a_bond_sized_glyph` | yes → 15 |
| | | `inject_a_small_glyph` (negative control) | **no**, correctly |
| anchor | full InChIKey match | `seed_epimer` | yes → mismatch, formula/mass unchanged |
| AST gate | CLEAN | `seed_banned_call` | yes → not clean |
| emit | callable confidence accepted | literal `"high"` | raises `EmitError` |
| | | caveat with empty note | raises `EmitError` |
| | | out-of-vocabulary `"probably"` | raises `EmitError` |

**Two no-op seeds were written first and caught.** The chain check's original
`drop_a_bond` removed the short carbonyl→oxygen stub rather than a bond, and
`inject_a_glyph` injected something the size filter is *supposed* to exclude.
Both printed a clean `14`. They were replaced with seeds that move the count,
and the small-glyph seed was kept as a negative control.

## A refusal that was respected

`curation/lib/chains.vector_bonds` in its default stroked-only mode returns
**zero** segments on SI p.7 and `vector_chain_length` then raises — correctly.
This page is a second measured instance of the class that module documents from
`andersondg_41203968`: **every bond is a FILLED polygon with `width = None`**
(603 of page 7's 915 paths). The refusal was not silently widened. The chain was
instead measured from the PDF's own path bounding boxes under an explicit
uniformity gate (bond quads 5.92 × 3.87 pt, height spread 0.06 pt), which is
not geometry inference, and the whole arrangement was seeded.

## What was NOT verified

* **Compound 92's macrocycle.** Drawn as four ring NH with three propylene and
  one ethylene bridge; the ring size was not confirmed by a second method, so
  it is left unassigned rather than named.
* **The `R = C10/C12/C13` tail labels.** The alkyls are text labels, not drawn
  chains, so there is nothing to count. Whether `C13` includes the two epoxide
  carbons is not stated anywhere in the paper.
* **Degree of substitution** on every compound — the central reason nothing is
  emittable.
* **Which aggregation SI Fig S4a's `% TdTom+ Spleen ECs` reports** (n = 4 or
  n = 34). See `08_contradictions.csv` C-7.
* **Fig 4B/4E/4F bar values.** Not digitised: the structures are undeterminable,
  so a measured bar would key onto nothing.
