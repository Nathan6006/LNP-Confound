# siegwart_40644303 -- provenance

**Status: EXTRACTED 2026-11-10. NOT merged.** `DECISIONS.md` carries 10 `OPEN`
entries and the merge gate refuses the folder while it does; that is correct,
not a blocker to clear.

| | |
|---|---|
| Paper | Lee SM, Sun Y, Chatterjee S, Xiong H, Cheng Q, Wang X, Siegwart DJ. *Structure-Activity Relationship of Ionizable Lipids for siRNA and mRNA Lipid Nanoparticle **Design***. ACS Biomater. Sci. Eng. 2025, 11, 4844-4852 |
| DOI | 10.1021/acsbiomaterials.5c00463 |
| PMID | **40644303** (indexed; PMCID PMC12869450, NIHMSID NIHMS2133731) |
| Experiment_ID | `SL_2025` (in vitro, and the Figure 2b in-vivo cohort) + `SL_2025_2` (the Figure 3b in-vivo cohort) -- the split is PROVISIONAL, see D-1 |
| Rows emittable | **505** = 465 in vitro + 40 in vivo (of 42; 2 unreadable, D-3) |
| Structures | **489** = the 465-member screen + 24 new in-vivo variants |
| Value source | see `00_value_sources.json` |
| Licence | **not CC0/CC-BY-compatible** on the evidence we have. D-7. |

## Where the values come from

Everything that is a delivery endpoint in this paper is `5_pixels`. That is not
a shrug, and `00_value_sources.json` gives the checkable reason: the JATS XML
has `table-wrap` count **0** (the article contains no marked-up table at all),
the single supplementary item is the SI PDF we hold, `find_repo.py` returns 0
candidate repositories for the DOI/PMID/title/author, and `pymupdf` over all 28
SI pages returns no luminescence or radiance number anywhere -- the SI's only
numbers are diameters, PDIs, pK_a, E.E. percentages and NMR shifts.

The one tier-4 artefact, SI Figure S3's text table, is used for everything it
carries and is what makes the digitisation checkable: because we hold every
lipid's pK_a exactly, the unlabelled markers in Figure 5a can be assigned to
named lipids by their x coordinate and their y read independently of the bars.

## The paper's own cross-checks

Ten, in `07_qc.json` Q7: 9 pass, 1 fails (the "100 to 200 nm" claim is 40/42).
**Passing all of them would be necessary and not sufficient.** Every one is a
ranking or a comparison, so a uniform multiplicative error passes all ten --
`07_qc.json` seeds exactly that error and records that Q7 does not notice it.
The load-bearing checks are geometric:

* **Q2** -- the 465 cells' assignment, tested against `KZ_2016`: Spearman
  **0.899**, collapsing to <= 0.60 under every one-row and one-column shift.
* **Q4** -- a drawn rule nobody printed: the x-axis baseline must land on the
  zero the y ticks fit. Measured **0.41 px** (Fig 2b) and **-0.20 px** (Fig 3b).
* **Q5** -- two unrelated panels on the same 42 numbers: median
  scatter/bar = **1.0124**, sd 0.038, over 32 lipids whose marker lands within
  **0.016 pK_a units** of its independently known value.
* **Q8** -- the SI's own 1H NMR integrals, summed: **19 of 25** peak lists sum
  exactly to the structure's hydrogen count, 3A4-SC5 at 86 = 86.

## Contradictions

Ten, in `08_contradictions.csv`. Four matter: the Figure 1c Venn cannot be
reproduced from the Figure 1a heat map under the caption's own rule; SC3 is
drawn unbranched here and is neopentyl in the corpus; 3A2 is drawn as a
piperazine here and has two nitrogens in the corpus; and Figures 2b and 3b are
~55x apart with nothing linking them.

## C10 anchor

`make_reference.py` BUILDS it. Scheme 1's drawn product **3A4-SC5** is the only
complete lipid this paper draws, and the second row of
`reference_structures.csv` records that absence with an empty
`reference_SMILES` rather than leaving it invisible: C10 covers 1 of 489
structures here.

The generator writes an explicit **V2000 atom/bond table** and parses it with
`Chem.MolFromMolBlock`. No SMILES fragment, no string concatenation of
chemistry, no `RWMol`, no `AddBond`, no `RemoveAtom` -- while the enumerator
uses `curation/lib/assemble.attach`, i.e. `RWMol` + `AddBond` on SMILES
fragments. The two share RDKit's sanitiser and nothing else.

**And the obvious "independent" route would have been a copy.** The
enumerator's canonical SMILES for 3A4-SC5 already starts at the pentyl tail, so
writing the anchor as a SMILES read tail-first from the drawing would have
produced a byte-identical string -- the F-037 shape, in a folder where it was
deliberately being avoided. An atom-and-bond table cannot collide that way.
C10 reports `exact_match`; C16 reports all templates derived.
