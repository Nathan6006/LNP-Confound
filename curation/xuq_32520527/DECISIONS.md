# DECISIONS — XQ_2020 (Li et al., Bioconjugate Chem. 2020, PMID 32520527)

The pipeline detects that a source is silent; it does not get to decide what the
source meant. Everything below is queued, not resolved. `provenance/build_dataset.py`
refuses a folder with an `OPEN` entry, which is the intended state of this folder.

---

### D-1. Experiment_ID, and whether the three formulations are one experiment or three

The paper screens **eleven ionizable lipidoids** in **three formulations** each:
`R-O16CBA` (lipidoid alone), `-F1` (lipidoid/cholesterol/DOPE 4/1/1 by weight)
and `-F2` (4/1/4). All 33 bars are one flow-cytometry panel, one cell line, one
cargo, one dose, one timepoint, read on one y axis.

Proposed, not applied:

- **one** `Experiment_ID`, with the formulation carried in the row's
  cholesterol/DOPE columns. Precedent: `ZT_2026`, where a 486-formulation DoE on
  a single lipid is one experiment.
- `Experiment_ID` = **`YL_2020`**, by the corpus convention of
  {first-author initial}{surname initial}_{year} — first author **Y**amin **L**i.
  The folder is named for the corresponding author (Qiaobing **Xu**), which is
  the acquisition tracker's convention, not the corpus's. No `YL_*` exists in
  the live corpus (54 experiments checked 2026-11-10), but **several other Xu-lab
  papers are being curated in this same run** and a collision has to be settled
  across folders, not inside one.

Splitting into three would create three `Experiment_ID`s that share all eleven
lipids, which would feed Stage 5's shared-lipid concordance with pairs that are
the *same experiment* — the opposite of what that statistic measures.

STATUS: OPEN

---

### D-2. The acetal stereochemistry, and a stereo guard that cannot see it

The cyclic benzylidene acetal is the paper's entire point: it is what makes the
lipidoid pH-responsive. Its 1,3-dioxane ring has a **C2/C5 cis–trans
relationship** — the C2 aryl group sits cis or trans to the C5 methyl. The paper
draws no wedge anywhere, states no ratio, and reports no NOE. The emitted SMILES
is therefore flat.

Two things make this worth a decision rather than a footnote:

1. **Neither carbon is a CIP stereocentre on its own** (each sees two
   constitutionally identical ring branches), so this is a *para*-stereocentre
   pair, and `assemble.stereocentres()` — which is what
   `assemble.attach`'s `check_stereo_preserved` relies on — returns `{}` for the
   tail. `Chem.FindPotentialStereo` returns two unspecified centres.
   **`curation/lib`'s stereo guard is blind to the only stereochemistry in this
   molecule**, so its passing is not evidence here. Measured in
   `01_blocks_qc.json → tail`.
2. Acetalisation under thermodynamic control usually favours one diastereomer,
   so "flat" may be *under*-specific rather than neutral — but the paper gives
   nothing to pick with, and the run's rule is flag, do not guess.

Options: (a) emit flat with `structure_assumed = 0` and the note as written —
what stage 05 currently does; (b) emit flat with `structure_assumed = 1` on all
33 rows, on the view that omitting a configuration the chemistry turns on is
itself an assumption; (c) hold the paper out of the corpus.

The ring **regiochemistry** — which two of the three CH2OH groups close the ring
— is *not* part of this decision. It was a real question and it is settled: the
three CH2OH groups of 1,1,1-tris(hydroxymethyl)ethane have identical canonical
rank, so all three pairings give one constitution (`07_qc.json` P10b).

STATUS: OPEN

---

### D-3. Amine 77's unspecified pyrrolidine stereocentre

Head **77** is 2-(2-aminoethyl)-1-methylpyrrolidine. Ring C2 is a genuine
tetrahedral stereocentre. Figure 1C draws it flat; the Experimental Section says
only "commercially available amine-containing head groups" and names no supplier
or stereochemistry for it. The commercial reagent is normally sold as the
racemate, but the paper does not say so and this curation will not say it for
them.

Stage 05 currently emits 77's three rows with `structure_assumed = 1`,
`confidence = medium` and a note naming the feature. That is the only difference
between 77 and the other ten, and it is **3 of the 33 rows**. Confirm or change.

STATUS: OPEN

---

### D-4. The licence is not CC0/CC-BY

`CLAUDE.md` §8 requires every constituent source to be CC0/CC-BY-compatible
before corpus release. PMC's permissions statement for this record reads, in
full: *"This file is available for text mining. It may also be used consistent
with the principles of fair use under the copyright law."* That is the NIH
Public Access Policy deposit of an author manuscript, not an open licence.

Every value in this curation is derived from a figure in that deposit. Whether
derived numerical values fall inside "text mining" is a question for the release
check, not for this folder.

STATUS: OPEN

---

### D-5. The figure disagrees with the prose on two of 33 values

`08_contradictions.csv` row C-1. The RESULTS text says 76-, 77-, 113-, 306- and
400-O16CBA-F1 "induced delivery efficacies of 78-83%". Figure 5 puts 76-F1 at
**84.93%** and 77-F1 at **88.51%**, 6 px and 17 px above the 83% line — far
outside any plausible read error — and it puts 80-F1, which the sentence does
not mention at all, at 85.25%, above three formulations the sentence does list.

Stage 05 emitted the figure. Every other printed statistic in the paper
reproduces from the figure exactly (two counts, two maxima, one range at both
endpoints, seven point values, 33 pairwise signs), so the evidence is that the
sentence is a loose summary. But this is the one place a number we emit
contradicts a number the authors printed, and `PIPELINE_V2.md` §3 says anything
that would change a published number is queued, not resolved.

STATUS: OPEN

---

### D-6. Tier-5 raster as the sole value source, at 1050 px

There is no table, no Source Data file and no repository; the 33 numbers exist
only as bar heights in a 1050 x 439 px JPEG (`00_value_sources.json` records why
tiers 1-4 fail, each with a checkable sentence). The measured quantisation floor
is **0.33 percentage points** — 3.07 px per percentage point — which is an order
of magnitude below the s.d. the authors draw on every bar, so the limiting error
is their replicate scatter and not our pixel grid. The author-manuscript PDF,
which might carry the figure larger, is behind a Google reCAPTCHA that this run
may not bypass.

Decide whether a sole-source raster at this resolution is admissible, given that
`NG_2025` and `WT_2026` set a precedent of tier-1/tier-3 sources with zero
digitisation defects while the four digitised papers contributed eight findings.

STATUS: OPEN
