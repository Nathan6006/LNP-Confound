# Change log

Every change from the original LNPDB (vitro) to `data/data_vitro.csv`. The first
section is the log kept during the ML work; the second records what was reversed
on 2026-10-28 and why. Row counts are as verified by `build_dataset.py`.

## 1. Changes from the ML work

### Structure corrections

| change | detail | verified in final |
|---|---|---|
Canonicalized all SMILES | RDKit canonical form | 24,524/24,524 canonical |
`Liu_iphos` (SL_2021) | `IL_SMILES` was the charged form; switched to neutral | formal charges {0: 494, +1: 78}, no anionic forms |
`Lee_Unsat` / `SL_2020` E/Z | LNPDB stored `/C=C\` (Z) for the 3T, 4T, 2T and SC6 tail series; should be `/C=C/` (E). 27 lipids | 69 trans motifs present; V had 0 |
`Lee_Unsat` / `SL_2020` duplicates | 28 pairs where the T variant was given the C variant's SMILES (56 rows), e.g. `2A9-4T` = `2A9-4C` | 28 of 29 resolved; `4A3-SC6`/`4A3-SC8` still collide |
`LX_2024` (Xue) and `SL_2020` (Lee) | SMILES mistakes fixed against previous data as ground truth | carried through from P |

### Additions

`Li_3CR` (`BL_2023_3CR`, 480) · `COMET_LANCE` (`AC_2025`, 1,978; dose imputed by
mean) · `FG_2026` amino-acid lipid library (960) · plus the Bowen Li lab
additions `MO_2026` (1,152), `YX_2026` (1,920), `YX_2026_SDL` (1,781).

### Columns added

A `Cargo` column; one-hot encodings for `HL_name`, `Model_type`, `Cargo_type`;
and five RDKit descriptors: `Num_carbon_in_tail`, `num_protonatable_nitrogens`,
`num_permanent_cationic_N`, `formal_net_charge`, `num_unsaturated_cc_bonds`.

The five descriptors are kept and recomputed. The one-hot columns are **not
stored** — see README.md for why and how to regenerate them.

### Deletions

| deletion | rows | status now |
|---|---|---|
`ZL_2022` (floor −0.41, skew 4.61) and `YX_2023` (floor −0.33, skew 3.3): z-scored without taking log of RLU, so not comparable | 392 | **kept deleted** |
OOD cell lines: BMDC (36) + BMDM (12), a small tail inside `LM_2019` which is mostly HeLa; HBEC_ALI (29) inside `JW_2024`, mostly A549 | 77 | **kept deleted** |
`YX_2024` (AGILE) rows run in RAW264.7 | 1,200 | **RESTORED** |
`XH_2024`, near-duplicate SMILES with no formulation variation (median pairwise Tanimoto 0.947) | 30 | **RESTORED** |
`ZC_2023`, 48% of data in the top label so no relevancy label assignable | 100 | **RESTORED** |

## 2. Reversals, 2026-10-28

Three deletions were reversed. The two unlogged-z-score removals and the 77 OOD
cell-line rows stay deleted, which is defensible: the first two are on a
different scale from the rest of the corpus, and the OOD rows are small tails of
cell lines that appear nowhere else.

The three restorations were motivated by their reasons being **modelling
decisions, not data-quality ones** — RAW264.7 "introduced noise", the XH_2024
lipids were too similar to each other, ZC_2023 would not take a relevancy label.
None of those is a statement that the measurements are wrong, so they belong in a
corpus whose purpose is to study cross-laboratory confounding.

Mechanics:

- **`YX_2024` RAW264.7 (1,200).** The whole experiment is taken from the curation
  file, where it is z-scored over all 2,400 rows. P's 1,200-row values are
  bit-identical to V's (slope 1.000000, residual 7×10⁻¹⁵), so P never re-z-scored
  after the drop and no rescaling was needed. P's 20 corrected structures are
  carried onto the restored rows via a canonical-SMILES map.
- **`XH_2024` (30)** and **`ZC_2023` (100)** exist only in the curation file and
  are taken from it directly, then canonicalized. Neither was touched by any of
  the structure corrections above, so nothing is lost by their not having passed
  through P.

Restoring these changes the category sets for `Model_type` and `HL_name`, which
is the direct reason the one-hot columns are no longer stored.

## 3. Normalization repairs, 2026-10-28

Eleven experiments were off the corpus's population-sd z-score convention —
five z-scored with the sample sd, six short by a couple of rows' worth of
variance (rows evidently dropped after normalization). All experiments already
centred at zero were rescaled to exact population sd 1, an affine map that
changes no ordering. All 52 experiments across both files now sit at mean 0,
population sd 1 to within 1e-9.

`SB_2024` and `SY_2025` (vivo, 60 rows) were NOT rescaled: their
`Experiment_method` is `luminescence_relative_to_Spikevax`, i.e. benchmark
ratios floored at 0, not z-scores. Flagged for a decision — they are the vivo
analogue of `ZL_2022` / `YX_2023`.

## 4. Structure patches, 2026-10-28

| experiment | lipid(s) | was | now | evidence |
|---|---|---|---|---|
`SL_2020` | `4A3-SC6` | SC8 octyl structure | hexyl tails, MW 1354.99 | series pattern: 4 tails × 2 CH₂ = 112.22 below SC8 |
`KZ_2016` | 42 × `*-SC3` | SC2 isopentyl structure | neopentyl, formula unchanged | source metadata records SC2 and SC3 at equal MolWt but tail chain 4 vs 3 |

48 further duplicate SMILES remain across `KW_2014`, `BL_2023`, `AA_2008` and
`JW_2024`. These are head-group collisions repeating across whole tail series,
not tail-length typos, so resolving them needs the source publications. Left
untouched and flagged — see README.md.

## 5. Not a change-log item

328 `KL_2010` rows carry no `Experiment_value` upstream and are excluded on that
basis alone. Confirmed by going back to `Original_LNPDB.csv`, where those rows are
literally `NA`; the curation file only re-encoded `NA` as empty. `KL_2010` is Love
et al., *PNAS* 2010 (PMID 20080679): LNPDB catalogues all 378 synthesized lipids,
but the paper quantifies HeLa siRNA silencing for only 50. The experiment keeps
its 50 valued rows.

`BroadPharm_2024` is the same condition taken to the limit — 269 rows in
`Original_LNPDB.csv`, all `NA`, and absent from the curated corpus entirely.

## 6. Curated papers merged, 2026-10-28

| Experiment_ID | paper | vitro | vivo | validation |
|---|---|---|---|---|
`ZT_2026` | Tian et al., Adv Mater 2026, PMID 42634137 | 486 | 51 | 7 independent checks, see `curation/siegwart_42634137/PROVENANCE.md` |
`LL_2026` | Lin et al., JACS 2026, PMID 41689530 | 840 | 0 | HRMS 4/4 exact; 48/48 printed hit rates reproduced |
`LL_2026_2` | same paper, Fig 2a initial screen | 147 | 0 | same chemistry; adds 7 amine head groups |

Both arrive already z-scored within their own experiment blocks, so
`renormalize_ddof` is not re-run over them; only the 5 extra descriptors are
computed and rows are routed by `Model`. Index and LNP_ID are reassigned above
the running maximum, so a rebuild can never collide.

### `LL_2026` — Lin et al., aromatic ionizable lipids

840 lipids = 5 amine precursors (8A–12A) x 7 aromatic analogs (1B–7B) x 12
hydrophobic tails x {1 equiv, 2 equiv of aromatic analog}, screened in IGROV1
cells with Fluc mRNA. Figure 2e is a 4-bin discrete heatmap (`<1e4` / `1e4–1e5`
/ `1e5–1e6` / `>1e6`), so cells are nearest-neighbour matched to legend
swatches: max colour distance **2.24** on a 0–255 scale, **zero** ambiguous
cells.

Structures were enumerated from the reagents named outright in the SI Materials
list, with epoxide ring-opening and Michael addition decided per block from its
own functional group rather than a lookup (7B is a methacrylate and gets an
alpha-methyl arm accordingly).

Validation:
- **HRMS 4/4 exact**, to ±0.0001 Da, against the SI's `[M+H]+` values for
  10A2-3B-EP10, 9A2-3B-EP10, 9A2-4B-EP10 and 10A2-3B-EP12. This confirms the
  reaction model and the arm-count logic (2 aromatic + 4 tails on a 6-N–H amine).
- **48/48 printed hit rates** in Figs 2f/2g/2h reproduced from the digitized
  heatmap, max deviation **1.4 percentage points**, with tick calibration landing
  independently on 20%/tick in all three panels.
- Aromatic-ring counts partition exactly as the design predicts: 300/300/120/120
  for 1/2/3/6 rings per lipid.

Flagged: all 840 carry `positional_ambiguity` (the paper does not say which
nitrogens take the aromatic arms, so the iPhos convention is applied), and the
140 rows using Ac8-1/Ac10-1 are medium confidence because the SI names them
*iso*-octyl and *iso*-decyl acrylate — commercial isomer mixtures, so a single
SMILES is a drawing convention.

#### Regiochemistry: corrected, and the ambiguity is genuine

The aromatic arms go **one per nitrogen**, primary before secondary. This is not
a convention pulled from the iPhos work — the SI's "precise synthesis" of
10A2-3B-EP10 settles it: spermine's primary amines are masked as
trifluoroacetamides, the secondaries are Boc-protected, the **primaries** are
then freed and reacted with **2.0 eq of 3B** (one aromatic each) before
deprotection and capping with EP10. An earlier pass put both aromatic arms on a
single primary nitrogen. The molecular formula is identical either way, so the
HRMS check could not distinguish them; the connectivity was simply wrong and is
now fixed.

The `positional_ambiguity` flag nonetheless stays on all 840, and that is
correct rather than cautious. The SI's General Methods describe the library
itself as made without protection or purification — "one equivalent of the amine
analogue was reacted with one or two equivalents of aromatic epoxides or
acrylates ... the final crude products were obtained by removing the solvent
under vacuum." The 840 library members are therefore genuine positional-isomer
mixtures. Only the four representatives carried through the protection route are
single regiodefined compounds.

#### `LL_2026_2` — the 147-lipid initial screen (Fig 2a)

Collected as a separate experiment. 7 aromatics (1B–7B) x EP10 x 21 amine/equiv
columns = 147. It differs from the 840 in formulation (15/20/25/2 at IL:mRNA
20:1, the initial formulation, versus 15/15/25/2 at 10:1), which is why it gets
its own Experiment_ID rather than being pooled.

Worth having: it uses **all twelve** amine precursors, so it contributes seven
head groups (1A–7A) that appear nowhere else in the corpus. 70 of its 147
structures recur in the 840 under the other formulation, and the two screens
agree on bin for 52 of those 70 (74%) — reported as information, not as a
pass/fail check, since the formulation genuinely differs.

Its legend has three bins, not four: nothing fell below 1e4 RLU. Bins are
indexed 2–4 on the same RLU scale as Fig 2e so the two screens stay comparable,
rather than 1–3 which would silently shift every value by a decade.
Colour match: max distance 2.29, zero ambiguous cells.

Not extracted: the mechanism panels (ITC, Laurdan, FRET, pKa, DLS) are not
library data.

## 7. Publication-metadata corrections, 2026-10-29

Surfaced by an independent manifest audit, each verified against NCBI before
being applied. Implemented as `METADATA_PATCHES` in `build_dataset.py`, so a
rebuild reapplies them.

| experiment | field | was | now | evidence |
|---|---|---|---|---|
`SB_2024` | PMID | `29786478` | `38779592` | `29786478` is *"Multisociety Consensus ... Endovascular Therapy of Acute Ischemic Stroke"*, Int J Stroke 2018. **Wrong upstream** — the same value is in `Original_LNPDB.csv`. The experiment's own DOI resolves to 38779592 (Bae SH, *Bioact Mater* 2024), matching the `SB`/`2024` convention. |
`XH_2025` | DOI | `10.1101/2025.02.25.640222v1` | `…640222` | bioRxiv version suffix; the registered DOI has none, so `[doi]` lookups on the stored value failed |
`JW_2024` | link | `…s41587-024-02490-y` | `…s41565-023-01548-3` | PMID and DOI agree with each other and with NCBI (Jiang AY, *Nat Nanotechnol* 2024); only the link was wrong |
`MO_2026` | link | `molea.csv (user-supplied HTS…)` | DOI resolver | the field was not a URL |

Also corrected in the docs: the vivo corpus has **18** experiments, not 17. The
17 was the pre-merge count and was not updated when `ZT_2026`'s 51 vivo rows
were added.

`provenance/build_manifest.py` → `data/experiment_manifest.csv`: one row per
Experiment_ID per corpus (58 rows), with PMID, DOI, first author/journal/year
from NCBI, cell line, readout, library size and a `curated_by_us` flag. Intended
for pasting into a chat session to check whether a paper is already curated.

## 8. Structure corrections after independent verification, 2026-10-29

An adversarial verifier subagent audited `curation/siegwart_42634137/` and found
**three real structural errors**. I confirmed all three against Fig 1b/Fig 2a and
RDKit before changing anything.

| # | error | impact |
|---|---|---|
E1 | Cholesterol esterified at a terminal methyl of the **isooctyl side chain** instead of **C3**. `CHOL.replace("(O)","",1)` deleted the hydroxyl, and SMILES concatenation bonded the acyl to the string's first atom. | all 36 synthesized structures wrong |
E2 | The deleted hydroxyl left C3 as a two-coordinate carbon with a **radical electron**, and lost the 3β stereocentre. Formula/MolWt each 1 H low. | all 36 wrong; 537 corpus rows |
E3 | Quaternization site chosen by `max(atom index)`, not chemistry. For C6 it picked the terminal N bearing two β-ester arms — the least nucleophilic — while **Fig 2a plainly shows a piperazine ring N**. | C5/C6 wrong; C1–C4 unjustified |

**Why every existing check passed.** Mass and formula cannot see connectivity,
and E1 left the oxygen count unchanged — only the hydrogen count moved, by one.
The structures parsed, canonicalised, and had self-consistent descriptors. This
is the concrete argument for CLAUDE.md §6.

Fixed: the C3 oxygen is kept and esterified; the quaternization site prefers a
ring N, then fewest β-ester arms, then atom index. C5/C6 now match Fig 2a; the
24 C1–C4 rows carry `quaternization_site_assumed=1` and confidence `low`, since
the paper draws no structure for them. The same cholesterol bug existed in a
second copy inside `curation/build_new_curation.py` (the 15 commercial-core
conjugates) and was fixed there too. Both call sites now hard-fail on a radical.

C6O2B2: **C117H218N5O12+, MW 1887.05** (was C117H217N5O12+, 1886.05).

`provenance/audit.py` check C1 now tests for **radical electrons** — the gap
that let this ship.

Two further findings triaged and **not** errors: the Methods-4.8 dose/replicate
contradiction **is** recorded (in `curation/new_curation_qc.json`, which the
verifier's grep scope missed), and the corpus `molratio` columns were always
correct percentages — only the intermediate `04b` file mislabelled ratios as
`_molpct`, now renamed with genuine percentages added alongside.

## 9. Pending

Nothing. `curation/new_curation.csv` (Tian et al. 2026, 537 rows: 486 in vitro +
51 in vivo) was listed here as unmerged; it has been in `build_dataset.py`'s
`MERGE` since 2026-10-28 and lands as `ZT_2026`. Every curation folder that is
meant to be in the corpus is in `MERGE`; the two deliberate exclusions
(`xd_2026_lipidai`'s 485-row LipidDB compilation, and the four benchmark LNPs of
`ql_2026`) are documented at their entries in that list, not here.

## 10. 2026-10-29 — two component structures corrected, FALCON merged

Both errors were found while curating `wht_2026_falcon` (Toh et al., FALCON,
bioRxiv 2026) and both were confirmed before any file was touched. They are
applied by the new `COMPONENT_PATCHES` table in `build_dataset.py`, which does
whole-value replacement keyed by `(Experiment_ID, column, component name)` —
`PATCHES` could not express them, being substring edits on `IL_SMILES` only.

### 18PG was the wrong lipid — `YZ_2022`, `YZ_2024`, 360 rows

The corpus stored `18PG` as a symmetric 18:0/18:0 lipid (DSPG). Avanti's "18PG"
is not symmetric. `YZ_2022` (Nat Commun 2022, PMC9310361) defines it in its own
Materials section:

> the anionic 1,2-dimyristoyl-sn-glycero-3-phosphate (14PA) and
> **1-stearoyl-2-oleoyl**-sn-glycero-3-phospho-(1'-rac-glycerol) (18PG)

That is SOPG: sn-1 stearoyl 18:0, sn-2 oleoyl 18:1 cis-9.
C42H82O10P⁻ (777.5651) → **C42H80O10P⁻ (775.5495)**.

`YZ_2024` is not open access, but it is the same lab (Hai-Quan Mao, JHU), and
the two experiments share **all 360 formulations**, the same single ionizable
lipid and the same six helper lipids — one platform, one reagent. FALCON, also
Mao lab, cites the same shorthand and vendor. `14PA` and `DDAB` from the same
panel were checked and are correct.

### `AC_2025` had MC3 and KC2 mislabelled — 419 rows (REVISED, see below)

By definition DLin-MC3-DMA is (6Z,9Z,28Z,31Z)-heptatriaconta-6,9,28,31-tetraen-
19-yl 4-(dimethylamino)butanoate, an **acyclic ester**; DLin-KC2-DMA is
2,2-dilinoleyl-4-(2-dimethylaminoethyl)-**1,3-dioxolane**, a ring with no ester.
`AC_2025` had the dioxolane under the name `DLin-MC3-DMA` (203 rows) and the
ester under `KC2` (216 rows).

**The first fix was backwards and was replaced the same day.** It rewrote the
STRUCTURES to match the names. But which half is wrong decides which measured
value belongs to which molecule, and the upstream source settles it:
`github.com/alvinchangw/COMET`, master table
`experiments/data_json/in_house_lnp_data_overall_new_full_without_pbae_NPratios_updated_09222023_npratios.json`
(3,119 entries), stores each LNP as component SMILES plus labels and carries
**no lipid names at all**. The names entered the corpus downstream, from
`LNP_ML/data_files/COMET_LANCE/convert_json_to_csv.py`, whose
`IL_SMILES_TO_NAME` dict is annotated *"verified via RDKit MW"* — and MW cannot
separate two molecules that are both C43H79NO2 at 641.6111 Da. That dict maps
the dioxolane to `DLin-MC3-DMA` and the ester to `KC2`, the wrong way round.

So the structures carry the values and only the labels were wrong. The fix is a
**name patch keyed by canonical SMILES** (`NAME_PATCHES`), not a structure
patch. Verified against the pre-patch snapshot: **0 structures moved, 0 values
moved, 419 names relabelled**, and the (structure, value) multiset is identical
to the original corpus. `DLin-MC3-DMA` now names the ester (216 rows) and `KC2`
the dioxolane (203 rows) — the counts swap, which is the visible sign the fix
went the other way.

**Why nothing caught it.** Identical formula, identical exact mass. Mass,
formula, monotonicity and descriptor self-consistency are all blind to it. Worth
stating separately: a check that had only verified "name matches structure"
would have been satisfied by either repair, including the wrong one. The source
of truth had to be the upstream data, not internal consistency.

### Effect on `data/`

360 `HL_SMILES` rewritten (18PG) and 419 `IL_name` values relabelled
(`AC_2025`). No structure and no `Experiment_value` moved anywhere in the
corpus. Verified by a row-level diff against the pre-patch snapshot.

### FALCON merged — `WT_2026`, 2,544 rows

`curation/wht_2026_falcon/` holds the full curation: **2,544 measurements over
1,218 formulations**, read from the authors' own tables in
`github.com/MaoResearchGroup/falcon-lnp-optimization` (commit `cfe987f`), so
nothing is digitized from a figure.

**All 2,544 rows merge**, on seven cell lines — Ramos 959, HepG2 513, THP-1
409, DC2.4 180, NIH-3T3 180, C2C12 180, INS-1E 123 — under one `Experiment_ID`,
`WT_2026`: a source is a paper, and different cell lines from one paper are one
source. Five cell lines are new to the corpus (`Ramos`, `THP-1`, `INS-1E`,
`NIH-3T3`, `C2C12`), which is new information rather than a reason to discard
measurements; an earlier version of this entry held them back and that was
dropped as indefensible.

The fourteen contributing blocks were normalized by the authors against their
own per-file controls, and the repo warns values "cannot be cross compared
across files", so each block keeps the z-score computed inside its own
normalization unit and those are concatenated rather than re-pooled. The result
still meets the convention exactly: standardized blocks concatenate to mean 0
and population sd 1 by construction, and the merged block measures
1.15e-17 / 1.000000000000.

`data_vitro.csv`: **25,997 → 28,541 rows, 40 → 41 experiments.** Audit C1–C9,
C11 clean. Two known gaps, both recorded rather than papered over:

* `WT_2026` is the only experiment with a blank `Publication_PMID` — it is a
  preprint. Left blank rather than invented.
* C4 reports **one exact duplicate row**. It is in the authors' table: F1109
  (iteration 1) and F1120 (iteration 2) are the same composition with a
  byte-identical `nLnLE_RAMOS` (2.483377939) but different `nLnLE_HEPG2` (3.237
  vs 4.414), so the Ramos value appears carried forward between iterations.
  Kept — it is their published data, not a read error.

### C11 — named commercial reagents (new Tier-1 check)

`provenance/reagent_registry.py`, wired into `audit.py`. Pins 13 catalogue
reagents to structures written from their chemical names and checks every row in
the corpus by InChIKey skeleton. C10 cannot cover this class: it compares lipids
a curation folder ENUMERATED against a structure a paper DREW, and catalogue
reagents arrive already named, from no figure.

Verified to have teeth: run against the pre-patch snapshot it flags exactly the
779 affected rows (18PG in `YZ_2022`/`YZ_2024`, MC3 and KC2 in `AC_2025`) and
nothing else. Against current `data/`, all 13 are clean.

PEG-lipids are deliberately absent — the block length is a distribution, so no
single SMILES would be truthful. Add an entry whenever a new catalogue reagent
enters the corpus.

## 11. 2026-10-29 — `XD_2026` and `QL_2026` merged

Recorded here late: §10 was written while FALCON was the last merge, and these
two landed after it. Both are in `build_dataset.py`'s `MERGE`.

### `XD_2026` — Duan et al., Acta Pharm Sin B 2026 (PMID 42592089), 12 in-vivo rows

The paper's **own** screen only: E1–E9 plus ALC-0315, SM-102 and DLin-MC3-DMA,
BALB/c IV, whole-body IVIS luciferase at 6 h (Fig 5). `data_vivo.csv`:
**1,768 → 1,780 rows, 18 → 19 experiments.**

The 485-row LipidDB in the same repo is **deliberately not merged** — it is a
compilation across the studies that paper cites, so giving it one
`Experiment_ID` would fabricate a laboratory. The curation folder's own
`06_lipiddb_qc.json` records **56** distinct reference strings after
normalisation and **50** letters-only, against the 47 the paper's §2.1.1
claims - close enough that the discrepancy is plausibly a few merged or split
citations, not a 6x miscount. An earlier changelog entry here asserted "299
distinct publication identifiers"; that figure came from a reviewer's naive
`doi:`-only cut, was never checked against the folder's own recorded numbers,
and is WRONG. The exact count depends on the de-duplication rule (nine
defensible rules give answers from 50 to 116), so no single number is asserted;
the folder's 56/50 and the paper's 47 are all recorded instead. It is extracted standalone in the curation folder
as a second, independent literature-pooled corpus.

Structures were traced geometrically from Fig 5A (skeletonise → trace →
Douglas-Peucker → count vertices), never counted by eye. The benchmark control
found an error in the paper: Fig 5A draws ALC-0315 with a 3-carbon hydroxyalkyl
head where the compound has 4, while SM-102 and MC3 in the same panel reproduce
the corpus exactly. The three benchmarks are therefore emitted with corpus
SMILES, E1–E9 carry `confidence = medium` (no name, formula, HRMS or SI — no
external anchor of any kind), and C10 records the ALC-0315 override as
`expected_mismatch_paper_drawing_is_wrong`.

### `QL_2026` — Liu et al., Adv Healthc Mater 2026 (PMID 41145968), 344 in-vitro rows

344 thiosuccinimide lipids in HeLa, from a bispecific-linker combinatorial
library: 8 amine cores × 8 NHS-ester/maleimide linkers × 5 linear tails (320)
plus a 24-member six-tail library on the T3A core. `data_vitro.csv`:
**28,541 → 28,885 rows, 41 → 42 experiments.**

Values: 293 from the Fig 2b heatmap, 27 from the Fig 2c bar chart (which lifts
the heatmap's 6e4 colorbar ceiling for the strong hits), 24 from Fig 3c. The
four benchmark LNPs read off the same panels (MC3, SM102, ALC-0315, Lipid 5) are
**not emitted** — they are the calibration that ties the two libraries together,
not library members.

**One pooled z, not two.** Libraries 1 and 2 were measured in separate runs, and
their four shared benchmarks give a consistent scalar offset of **1.1116×**
(sd 0.0234), which is applied to library 2 before a single pooled z over all 344
rows. An earlier version standardized the two libraries separately; that is
wrong here and the error is instructive — separate standardization forces both
block means to 0, which deleted the paper's central finding (the six-tail
library's geometric mean is 154,630 RLU against 4,128 for the main library, ~37×)
and inverted the rank of 36% of cross-library pairs, while every validation
check still passed. Recorded as error class 21.

Flags carried into the corpus: 3 rows censored (above the Fig 2c axis),
113 floored (at the heatmap's lower bound), and 1 conflict — `T3A-L8-C14`, which
Fig 2c plots among the hits while its own Fig 2b cell reads 7.9e3, a 7.3×
disagreement between two figures of one paper. The heatmap value is carried, and
the conflict is flagged rather than reconciled.

C10 anchor: `T3A-L5-C10`, assembled from the Fig 2a building blocks by a
different code path than the enumerator — exact match. 21 validation checks in
`curation/ql_2026/validate.py` pass, including all three of Fig 3d's printed hit
rates, Fig 3e's eight per-linker rates, and the Results sentence "17 lipids …
greater than 1e5".

Its `Publication_PMID` was initially left blank as "not indexed"; PMID 41145968
was in fact already indexed with a matching DOI and first author. Corrected in
`curation/ql_2026/build_rows.py`. `WT_2026` (a preprint) is now the only
experiment with a blank PMID.

---

## 12. 2026-10-30 — round-1 corpus review: repairs applied

Context: `provenance/review/REVIEW_PROGRAM.md`. Round 1 ran eight read-only
reviewers (five per-paper, three corpus-level) over the five curations merged
since 2026-10-28, logged 21 findings in
`provenance/review/round_2026-10-30/FINDINGS.md`, and sent each repair to a
separate validation agent before it was applied. Two validators changed the
plan materially; both changes are noted below.

Corpus after: **`data_vitro.csv` 28,885 × 79, 42 experiments; `data_vivo.csv`
1,780 × 79.** Row counts unchanged. The schema went from 72 to 79 columns.

### Structures changed

**All 537 `ZT_2026` rows: epicholesterol → cholesterol** (F-013). Every row of
the experiment stored the 3α epimer. `enumerate.py::_cholesteryl_ester` bonded
the acid's oxygen to the sterol's C3 and then called `RWMol.RemoveAtom` on
cholesterol's own hydroxyl oxygen, which drops that oxygen out of C3's
neighbour list **without correcting C3's chiral tag** — so the stored tag came
to mean the opposite of what it meant. Reproduced from scratch: feeding genuine
cholesterol (InChIKey verified on input) through that code path gives C3 = `R`;
the same ester written as a plain SMILES string gives `S`.

Nothing could see it. The epimers share molecular formula `C117H218N5O12+`,
exact mass 1885.6597, the InChIKey **skeleton** block, and all 17 descriptor
columns bit-for-bit. C10 — the one check built for this — was comparing a value
to itself: `reference_structures.csv` had the same sha256 as `03_enumerated.csv`
(223 bytes), while its own note claimed independent fragment assembly "through
the C3 oxygen (3-beta)". C11 was blind twice over: it compares skeleton blocks
only, by design, and it inspects the `CHL_*` columns, which are empty on all 537
rows because this paper's claim is cholesterol-free LNPs. The sterol lives
inside `IL_SMILES`.

Fixed by doing the condensation the way the chemistry does it — the **acid**
loses its OH, the **alcohol's** oxygen is retained as the ester oxygen — so no
atom adjacent to a stereocentre is ever deleted. Verified against an
independent plain-string oracle that uses no `RWMol` at all: exact match for
all three linkers, and the **full** CIP label vector matches, not just C3.

The validator settled the right-half question against the paper rather than by
assumption: Fig 1b (p.3) draws B1/B2/B3 with the C3–O bond as a bold wedge
cofacial with the C18/C19 angular methyls; Fig 2a shows the same relation in a
different orientation; the ¹H C3-H sits at 4.71–4.58 (the broad axial *tt* of a
3β-acyloxy sterol, where a 3α ester would be near 5.0 and far narrower); and the
Fig S7 ¹³C spectrum has an isolated peak at 73.9 ppm where a 3α ester would be
near 70.9. No epi, allo or 3α sterol appears anywhere in the paper or SI. So
rewriting the structures destroys nothing real.

**It also found a live second copy the plan had missed**: `build_new_curation.py`
carried its own `quaternize` with the identical `RemoveAtom` construction,
affecting **15 more structures** (`MC3-CB*`, `SM-102-CB*`, `OF-2-CB*`, `TT3-CB*`,
`Lipid 14-CB*`). Fixing `enumerate.py` alone would have left all 15 wrong after
the fix was announced. That second copy is now **deleted** — the file imports
the one implementation rather than restating it.

**`MC3-CB1/2/3`: built on the KC2 dioxolane under an MC3 name** (F-014).
Measured 5 rings and 1 dioxolane against corpus `DLin-MC3-DMA`'s 0 and 0; the
dioxolane belongs to DLin-KC2-DMA. The name is right and the structure was
wrong — the opposite of the `AC_2025` repair, settled by Fig 1b (p.3), which
draws "Dlin-MC3-DMA (MC3)" as the acyclic ester with no ring of any kind.

The cause is worse than "a stale lookup". `load_corpus_components` did
`comp.setdefault(name, smiles)` — first row in file order wins — over a corpus
written in alphabetical experiment order. **292 names in that table map to more
than one structure.** `DLin-MC3-DMA` was *ambiguous*, not simply wrong: 2,190
rows of the correct ester plus 203 dioxolane rows from `AC_2025`, and `AC_2025`
sorts first. Regenerating after the `AC_2025` name repair fixed these three rows
**by luck of row order**, which is why the lookup now **raises** on an ambiguous
name instead of picking one. Every name this curation uses is unambiguous today;
it will fail loudly if that changes.

**`Lipid 14-CB1/2/3`: alkylated at the hindered central nitrogen** (F-015), now
the terminal dimethylamino. The duplicated rule still used `max(atom index)`.
The validator's finding here changed the fix: **porting the repaired rule
verbatim would not have worked.** `beta_ester_arms` walked only `N-C-C-C(=O)`,
so it could not see Lipid 14's `N-C-C-O-C(=O)` β-**acyloxy** arm; both candidate
nitrogens scored 0 and the rule fell through to the atom-index tiebreak —
reintroducing the SMILES-ordering artifact it was written to remove. Rewriting
the same molecule tail-first made the two rules swap answers. β-acyloxy is now
counted alongside β-ester, with a steric tiebreak below it; the choice is
invariant across four spellings of each parent, and C5/C6 still resolve to the
piperazine N drawn in Fig 2a.

The site itself remains **unverifiable** — the paper draws no `-C` conjugate for
any of the five commercial cores, Fig 2a draws only C6O2B2, and the SI NMR
covers no conjugate. So all 15 carry `structure_assumed = 1` with the reason on
the row, rather than a silent structure.

No `Experiment_value` changed in any of the above.

### Values changed

**`WT_2026`: 2,544 rows re-standardised, 14 blocks → 7** (F-003/F-016). The
curation z-scored each (source file × cell type) column separately, on the
authors' warning that values "cannot be cross compared across files". That
warning is about *across files*; **within** a file the authors cross-compare
cell types themselves — every multi-output file carries an `nLnLE_Selectivity`
column, defined in the Methods as on-target minus off-target Ln(nLE), and the
curation's own solver already verified that identity to 1e-15. Standardising
each column separately forces that difference to **identically zero**.

Measured on file 5: raw mean HepG2 − Ramos = **+4.0342 ln (a factor of 56.5)**
with 173 of 180 formulations HepG2 > Ramos, stored as exactly **0.0000**, and
Spearman(raw difference, stored difference) = **−0.517** — anti-correlated, so a
formulation that was genuinely more HepG2-selective looked more Ramos-selective
in the corpus.

Now one z-score per source file. The restored E5 offset is **+1.4422**, which is
4.0342 / 2.7972 (the file's pooled sd) exactly, and 173 of 180 signs are back.
The old code asserted that each block sat at mean 0 / sd 1 — true by
construction for any set of standardised blocks, including the wrong one. The
new check is real: with one (μ, σ) per file, `z_a − z_b = Selectivity / σ`, so
the stored values must reproduce the authors' printed column up to one constant
per file. **Max residual 1.8e-15**, asserted, and it fails if the normalisation
unit is ever changed back.

### Schema: 72 → 79 columns

Seven columns, defined once in **`provenance/schema_quality.json`** and read from
there by the builder, by every curation that emits a flag, and by the audit. The
list is not restated in code: the reason these columns exist is that per-curation
flags were restated in six places, drifted, and died before the merge.

`Experiment_value_within_model`, `value_normalization_unit`, `confidence`,
`positional_ambiguity`, `structure_assumed`, `value_is_bound`, `curation_note`.

**F-001 — 1,152 rows carried a convention and said nothing about it.**
`CLAUDE.md` claimed the project flagged such rows "in the data
(`positional_ambiguity`, `*_assumed`, `confidence`)". None of those columns
existed, in `data/` or in any file `MERGE` reads — so the loss was one stage
*earlier* than the merge and no builder change could have preserved it. All
1,152 now carry their flag: `LL_2026` 840, `LL_2026_2` 147, `QL_2026` 117
(113 floored + 3 censored + 1 cross-figure conflict), `ZT_2026` 39, `XD_2026` 9.

The count moved twice, both times upward, and both corrections are mine to own:
I first recorded ≥266 (counting only flags enumerable from QC JSONs, missing
`LL_2026`'s 987 `positional_ambiguity` rows), then adopted CORPUS-B's 1,137 as
"the better number" — which undercounts `ZT_2026` by 15, because it took only
the 24 lipids with an assumed quaternization site and not the 15 commercial-core
conjugates whose cholesterol conjugation site is equally inferred. Measured
against the emitted corpus: **1,152**.
§11 above says "flags carried into the corpus"; that was not true when written
and is true now.

**F-011 — 24 `ZT_2026` rows had their disclosed `low` confidence replaced by
`high`.** The mechanism was not the join I had diagnosed: `build_new_curation.py`
read `lipid_name → IL_SMILES` and then wrote the literal `"high"` for every
lipid it found, discarding the `confidence` and `quaternization_site_assumed`
columns sitting beside it. The 24 whose record reads verbatim *"ASSUMED: no ring
N and the paper draws no structure"* were labelled the project's highest
confidence. They now read `low`, with that sentence as the row's
`curation_note`. This had to precede F-001: propagating `confidence` first would
have published `high` on 24 assumed structures, which is worse than publishing
nothing, because the only surviving provenance would assert the opposite of the
truth. (The separate `06_merged.csv` — an intermediate, not a merge input — does
have a joined `confidence` that is a pure function of the O-tail; it is read only
for `lipid_name → IL_SMILES` and is left untouched.)

**F-016 — `Experiment_value_within_model`, and units declared on the row.**
`YX_2024` (2,400) and `LZ_2025` (96) are also per-cell-type standardised, and
unlike `WT_2026` they **cannot be repaired**: both arrive already z-scored in the
upstream P/V files with no raw signal anywhere in the repository, so the deleted
offset is not recoverable from what we hold. Rather than leave a reader to
infer that, every row states its `value_normalization_unit` — `experiment`,
`source_file`, or `experiment_x_model_type` — and the two unrecoverable
experiments carry the reason in `curation_note`.

`Experiment_value_within_model` is computed for every row on the z-score scale,
by re-standardising `Experiment_value` inside each (experiment × `Model_type`)
block. That is exact, not an approximation: a per-block z of an affine function
of the log signal is the per-block z of the signal. The reverse is not true,
which is the whole of F-016.

### The gate

**F-012 — the Tier-1 audit could not fail.** `audit.py` discarded `main()`'s
return value; C1–C9 printed and the exit code came only from the C10/C11
delegates. Eight simultaneous seeded violations exited **0**. Of 43 mutations, 30
appeared in the report and 5 failed a command. `build_dataset.py` never imported
`sys` at all.

The obvious fix is wrong on its own: failing on anything printed would fail
forever on conditions we accept on purpose — 48 duplicate SMILES groups, 212
unassigned stereogenic C=C, 2,544 blank PMIDs on a preprint, 3 `SB_2024` rows
whose ratios do not sum to 100. That is probably why it was left print-only. So:

- **`provenance/audit_baseline.json`** — 62 gated metrics, each with the reason
  it is accepted. Departure in **either** direction fails. A decrease fails
  deliberately: a known condition vanishing without anyone editing the baseline
  is how a check that stopped looking announces itself, and that failure mode
  reads as a clean report.
- **`provenance/test_checks.py`** — seeded-mutation suite. One field, one row, a
  copy of the corpus, and a required non-zero exit, plus a control run that must
  still pass. Its first run found two coverage gaps: gating only
  `blank_Publication_PMID` let a blanked `IL_name` through, and **C8 checked 7
  of the 14 published descriptor columns**, so a `LogP` of 42 passed. Both
  fixed, both mutations kept.
- **`curation/siegwart_42634137/test_enumerate.py`** — unit tests at the
  function, not 537 rows later: the string oracle, the full CIP vector, and the
  ordering-invariance of the quaternization rule.

### New checks

| | what |
|---|---|
| **C12** | error class 22 — flags an experiment whose per-cell-type blocks each sit independently at mean 0 / population sd 1 while the row does not declare that unit. Found by hand four times and never by a check. Specificity tested synthetically, including the worst case (true offset 0 and equal sds), because the corpus had **no** correctly-pooled multi-cell-type experiment to serve as a negative control — all 3 were affected. `WT_2026` is now that negative control. |
| **C13** | the quality vocabulary, from `schema_quality.json`; also counts rows flagged without a note. |
| **C14** | cleaves every ester inside `IL_SMILES` **corpus-wide** and pins the liberated sterol's **full** InChIKey against `reagent_registry.json`. Corpus-wide on purpose: scoping a check to the experiment where the bug was found is how the same bug survives the next curation. |
| **C10** | a stereo-only difference now **fails** when the reference's stereochemistry came from a pinned authority rather than a figure; and an anchor byte-identical to the enumeration it checks, in a folder with no generator script, is reported as a possible copy. |

### Correction to F-016's own scope

F-016 recorded "only these three; every other multi-cell-type experiment is
pooled correctly." That is wrong, and the truth is worse: the corpus contains
exactly **three** experiments with more than one `Model_type`, and **all three**
were per-cell-type standardised. 3 of 3, not 3 of N.

### Open

- **F-021, new, found by C14 on its first run** in an experiment no reviewer was
  assigned: `ZH_2023`'s `A4I18-1R2C70` (2 vivo rows) carries a steroid ester
  whose sterol cleaves to **C26H44O**, skeleton `DZKYBYMQFBHTTM` — not
  cholesterol (C27H46O, `HVYWMOMLDIMFJA`). It is missing the C19 angular methyl
  and carries no stereochemistry at all. Gated at its current count so it cannot
  spread; resolving it needs the `ZH_2023` source paper.
- **F-017/F-018/F-019/F-020**, the digitization findings, are under
  re-measurement and are not applied.
- **F-008/F-009** (`LL_2026` Fig 2a sampling window; a dropped enantiopure
  stereocentre on 141 rows) are reproduced but not applied. F-009 should stay a
  flag, not a structure change: the paper reports no optical rotation or chiral
  HPLC of the products.
- **F-010** — C2 skips 737 of 972 series because the tail token is not a bare
  integer, including all 140 `LL_2026` series, so the largest curation
  contributes zero checked series while the line reads "0 violations". A
  coverage floor now stops that number falling further; raising it is open work.
  The `LL_2026` arm decomposition it also reported missing is now populated.
- **7 C10 anchors** are byte-identical to their enumeration with no generator
  script: 6 in `ng_2025` (another session's in-flight work — not inspected) and
  `lin_aromatic_2026`'s `10A2-3B-EP10`.

### F-006/F-007 applied — `LL_2026`'s hit-rate validation, and the bias it was hiding

`digitize_fig2fgh.py` computed `slope`, the pixels-per-percent that converts
every bar, as a through-origin least-squares fit of the measured bar heights
onto the Fig 2e hit rates **that this panel exists to check**. The tick spacing
was measured and used only to print a diagnostic. The file had no `assert`.

Now re-anchored: a least-squares line through the detected tick centres gives
pixels per tick, the interval is asserted to be 20 pp, the tick fit residual is
asserted ≤2.5 px (unequal spacing would mean these are glyph centroids, class
20), and the agreement with Fig 2e is **asserted** rather than printed.

**Re-anchoring it exposed a real error the fit had been hiding.** With the
tick-only scale, all 48 deltas came out **negative**, −0.4 to −1.6 pp — a
systematic bias, not noise. The colour masks select only fully saturated bar
pixels, so every top edge read a few pixels low, and on ~400 px bars a 4 px loss
is about 1%. The fitted scale was 5.046 px/pp against a tick-implied 5.100, a
1.07% difference, which cancelled the ~1.1% height shortfall almost exactly.

So the fitted calibration was not merely blind to a scale error — **it was
correcting a real edge-detection error with a fudge factor derived from the
answer**, and that is why the agreement looked so good. This is the sharpest
instance of error class 27 in the project.

Fixed properly: bar tops now interpolate the half-intensity point of the
anti-aliasing ramp, uniformly for every bar. Max |delta| **1.60 → 1.00 pp**,
mean −0.87 → **−0.52 pp**. A one-sided residual of −0.52 pp remains and is
**recorded rather than tuned away** — its probable cause is the bottom edge (a
bar sits on the axis rule and the mask's last saturated row is 1–2 px above the
rule's top, worth ~0.28 pp on panel 2h). It is left open deliberately: these 48
numbers are a validation of the Fig 2e read and no corpus value is computed from
them, so the residual costs accuracy in a check rather than in the data — and
tuning it away is precisely what the old fitted scale did.

**F-007**: `min_area=800` dropped EP6's `>1e6` bar (≈50×9 = 184 px² as the mask
measures it), after which `m6` fell through to the literal `0.0` — so a
measurement *failure* was recorded as a measurement *of zero*, and that
fabricated zero became the set's single largest "disagreement". The floor is now
120 (the next-smallest component on that baseline is 3,676 px²), a miss records
`None` instead of a number, and `missed` is asserted empty. EP6 now reads 0.7%
against Fig 2e's 1.4%.

The false claim is corrected at all three sites it was repeated:
`11_LL_2026_qc.json` (via `build_rows.py`), `lit_review/build_data.py`, and
`CLAUDE.md` §5.

### F-020 applied — `XD_2026`'s 12 bar tops

`stage04_fig5c.py` declared the bar top at the first row where >60% of the inner
bar width is fill. The three replicate discs straddle the top edge and occlude
the **middle** of the width, so the rule fired inside the bar.

Reproduced independently, and the mechanism is pinned rather than inferred: at
the true top row the inner fill fraction is **0.21–0.57** on the nine affected
bars and **0.68–0.97** on the three that are not, so the 0.60 threshold is
exactly the discriminator — the rule *is* the error. Re-applying it reproduces
the stored tops on **12 of 12 at exactly 0 px**.

Now read from the **flank columns**, which a centred disc cannot reach, taking
the modal first-fill row over the 5 leftmost and 5 rightmost columns and
asserting the two flanks agree within 1 px. A global `min` over all columns
would not have worked: on E1 the anti-aliased halo of a disc floating above the
bar makes columns 23–25 report a top 8 px too high.

Pixel corrections (600 dpi): ALC-0315 +13, SM-102 +7, MC3 +7, E1 +18, E2 +2,
E3 +9, E4 +21, E5 +15, E6 +18, E7/E8/E9 unchanged. Mean bias **+0.051 decades**,
max **+0.117** (E4). All 12 `Experiment_value`s move — including the three whose
pixel top does not, because the block mean and sd shift under them, which is the
demonstration that the bias does not cancel.

The axis calibration is **sound** and was re-verified: 179.960 px/decade at
600 dpi (89.98 native), major-tick residuals ≤4e-4 decades, checked against 31
unlabelled minor ticks at rms 0.53 px. So this is not class 20, and there is one
error here rather than two compounding.

**The paper's own cross-check could not have caught it.** All 12 lipids still
fall in their stated Fig 5D A/B/C band, and both printed Fig 5D numbers (75%
overall, class B 1/3) still reproduce — they reproduced equally well on the wrong
values. Worth noting: correcting E3 narrows its margin to the B/A boundary from
0.132 to **0.082 decades (≈15 px)**, the thinnest in the panel. The
`stage07_emit.py` docstring's claim that two agreeing figures check the
digitisation overstates what that check can do.

### F-010 applied — C2's coverage, 235 → 634 series

The corpus line read *"235 verified-length series checked, 0 violations (737
series skipped)"*. 76% of series were skipped, and **all 140 `LL_2026` series
were among them**, so the largest curation in the corpus contributed zero
checked series while reporting a clean pass.

The cause was `NAMING_IS_LENGTH = {"KW_2014"}` — a hand-maintained allowlist of
one experiment. Everything else was skipped even when its mass steps were whole
multiples of CH₂.

Replaced by a **per-series** test derived from the data, which generalises
KW_2014's own justification (1,271 of its 1,273 names agreed with the computed
longest chain): a series' name token is treated as a chain length only if it
equals the measured `Num_carbon_in_tail` for **every** member of that series.
Self-verifying, nothing to maintain.

It is finer-grained than an experiment can be: in `LL_2026` the `EP` series pass
— EP6→EP18 steps are exactly 4 CH₂ per index unit, 4 arms × 2 carbons — while
the `Ac8-1`/`Ac10-1` series correctly do not, because there the trailing 1 is an
isomer marker, not a length.

**Coverage: 235 → 634 series in vitro across 11 experiments (`KW_2014` 234,
`LL_2026` 140, `AA_2008` 94, `QL_2026` 64, `LX_2024_3` 36, `NC_2024` 20,
`YY_2023` 17, `JM_2016` 12, `ZC_2023` 12, `SW_2024` 3, `XH_2024` 2), and 12 → 22
in vivo. Zero violations.** The coverage floors in `audit_baseline.json` were
raised to match, so the number cannot quietly fall back.

338 series are still skipped and that is the correct answer for each: `KS_2024`
`C24` is two C4 tails (token 24, measured 6), `KZ_2016` `SC1/SC2/SC3` are tail
identifiers (isobutyl/isopentyl/neopentyl), `YX_2024` `C1/C2/C3` are component
indices. Worth stating because it was tested: all three pass the whole-CH₂ step
gate, so **that gate alone is not sufficient** — relying on it contributed 380
spurious violations.

### F-019 applied — `QL_2026`'s colorbar was anchored on a border edge

A second instance of **error class 20**, on a colorbar rather than an axis.

`digitize_fig2b.py` mapped values with `VMAX + (VMIN-VMAX)*(y-BAR_TOP)/(BAR_BOT-BAR_TOP)`,
putting 0 RLU at `BAR_BOT = 1973` — the **outer edge of the bar's black
border**. Verified directly from the pixel profile: the border runs y = 1962-1973
(dark from 1962, saturated 1964-1971), so its outer edge is about 5.5 px below
the ramp's true zero. Stretching 0-6e4 over 870 px instead of ~860 biases
**every** Fig 2b value high — about +1% at the 2e4 hit threshold and **+38% at
the floor**, where a third of library 1 sits.

The LUT *sampling* window and the *value* anchors are now separate constants;
conflating them was the bug. The value anchors come from the colorbar's own 4e4
and 2e4 tick notches on the native embedded JPEG, `Y_VMAX, Y_VMIN = 1107.4,
1967.6`.

**Cross-checked by an independent route before applying,** because this changes
all 320 library-1 values: the bar's bottom border *centre* measures 1967.5
against the tick-derived 1967.23, agreeing to 0.27 px, and a border-centre
anchor pair (1104.5, 1967.5) gives a white-point value of 903.8 RLU against
913.7 from the ticks. **The two corrected routes differ by 9.9 RLU while the old
anchor was +362 RLU out**, so the correction is robust and the residual anchor
uncertainty is ~3% of it. The top anchor is the weaker of the two and is
documented as such: the bar's top border abuts a near-black dark-blue ramp, so
its inner edge cannot be read from colour — which is exactly why the ticks are
the authority.

**The floor now reads 948.6 RLU instead of 1310.3.** That matters beyond
accuracy: the 48 white cells have a feasible band of [0, ~982] RLU, and the
stored 1310.3 was **above the entire band** (1.43× its upper bound). It
contradicted the `value_is_bound = upper` flag this round had just added to those
very rows. The stored value now sits at the top of the band rather than outside
it.

Also recorded on those rows: the ramp reaches pure white at ~948 RLU with no
gradient below; adjacent LUT steps near the floor differ by 1.7-3.5 RGB units
while each spans 349-1,045 RLU, so the limit is the colormap's own quantisation
and not image noise; and **56% of the slice carries one of only six values, all
from the bottom six LUT steps** (344 rows hold 83 distinct `Experiment_value`s,
the three largest ties being 65, 57 and 48 rows).

The paper's own between-library conclusion is unaffected and in fact strengthens
under the correction.

### F-017 applied — and it refutes what F-017 claimed

`digitize_fig2c.py`'s comment read *"y-axis label centroids: 1e6 2342, 1e5 2748,
1e4 3156"*, whose gaps of 406 and 408 are unequal. I took that as proof the code
was anchored on glyph centroids rather than ticks and marked the finding
`confirmed`.

**The code was right and the comment was wrong.** Those three integers are
neither the label centroids (2344.3 / 2750.8 / 3158.9) nor the tick positions —
they are the tick positions *rounded to integers*, and the rounding is the whole
of the apparent inequality. Measured tick centroids are 2342.78 / 2749.54 /
3156.75, a linear fit gives **406.99 px/decade** with residuals ±0.037 px, so
`PX_PER_DECADE = 407.0` was correct to 0.01 px. Only `Y_1E4` moves, 3156.0 →
**3156.75**, worth +0.42% in value and **rms Δz 0.0006** across the slice.

So F-017 is a class-21 defect (provenance text contradicting the code) and not
the class-20 defect it was logged as. An earlier reviewer's proposed 406.0 /
y=3153.0 is **refuted** by a round number the paper never printed: Fig 2c's red
dotted reference line sits at y = 2872.21, which these constants put at 50,020
RLU — 0.04% off a round 5e4, against 49,804 stored and 49,159 under that
proposal.

Fig 3c's axis was re-derived as a control and the earlier repair there is
**confirmed exact** (150.913 px/decade and y(1e3) = 2911.50 against stored
151.0 / 2911.5), so the sibling panel needs no change.

### F-018 confirmed, NOT applied

The mechanism and the count reproduce — the `w//8` trim plus a 0.85 coverage
threshold reads replicate dots instead of the bar top on **15 of 60 bars** across
Figs 2c and 3c, by up to +28% in value — and the stored cross-library offset of
**1.1116 is an artifact of it**: ALC-0315 is dot-contaminated in both panels
(+24 px in 2c, +16 px in 3c) and Lipid 5 in 3c (+8 px). The honest reading of
the re-measurement is that there is **no measurable offset** between the two
screens (1.0070 ± 0.0247).

**Not applied, because I could not validate a replacement rule.** A side-stroke
rule written for hollow bars failed on the first bar: Fig 2c's bars are 12 px
wide and solid. Two approaches are already ruled out — a stricter coverage
threshold (merged dots reach 91-97% contiguous coverage) and pure flank reading
(wrong on anti-aliased outer stroke columns). The original rule stays in place
with a comment recording that it is known biased, by how much, and where the
per-bar measurements live. A detector that has not been validated bar by bar is
worse than a characterised bias, because the bias stops being characterised.
Returned to the re-measurement agent for a patch plus a per-bar verification
table.

One argument is **retired**: F-018 cited two pairs of byte-identical bar values
as a bug signature. Both pairs are real and both are spurious reads, but exact
ties are not by themselves evidence — at 151 px/decade on a 2 px row pitch Fig 3c
offers only ~75 distinguishable values per decade, and the corrected set still
contains one exact tie.

### The gate, verified end to end

`provenance/test_checks.py`: **29 of 29 seeded mutations produce a non-zero
exit, and the unmutated control still passes.** That includes a C14 mutation that
re-creates F-013 by inverting the steroid C3 tag with RDKit — an earlier version
of that mutation was a string substitution that silently matched nothing, so the
suite now asserts that each mutation actually changed the row. A mutation that
cannot land makes a working check look broken, which is the worst direction for
this suite to fail in.

### F-018 applied — and my reason for reverting it was itself a measurement error

Applied after a second pass produced a patch validated bar by bar. Post-apply
check as specified: `12_QL_2026_qc.json` `cross_figure_offset` reads
**1.0071, sd 0.0247**, matching the tested patch exactly.

**The geometry.** All 60 bars in Figs 2c and 3c are **hollow** — in a band above
the baseline each shows exactly two fully-coloured column groups (the side
strokes) with 57-76% white between them. Fig 2c bars are 66-72 px wide, Fig 3c
46-50 px. **None is solid.**

I had reverted the first attempt on the grounds that "Fig 2c's bars are 12 px
wide and solid". That was wrong, and wrong in a way worth recording: 12 px is the
**left side stroke** of `110-L1-C10` (cols 557-568); the bar is cols 557-627, 70
px wide, 69% white inside. I got it because I reimplemented the `coloured` mask
instead of reproducing the script's — the script computes `sat = img.max(2) -
img.min(2)` and I used a normalised saturation, which fragmented the mask so only
the stroke registered. With the script's own mask the segmentation gives 32 bars
at 70-72 px, which I have now confirmed directly.

So reverting was the right call for the wrong reason. The first attempt actually
matched an independent measurement to ≤2 px on **56 of 60** bars; it needed two
specific fixes, not a new geometry model:

- **a luminance mask, not a blue mask.** Fig 2c's red dotted 5e4 reference line
  crosses the top rule of `DLin-MC3-DMA` and `SM102`; on a blue-only mask the
  rule is punched through there, and `DLin-MC3-DMA` returned no top at all.
- **strictly contiguous runs (gap = 0).** Bridging even a 2 px gap lets a merged
  replicate-dot cluster span both stroke centres, which is exactly what it did
  on `T3A-L2-2C10` and `ALC-0315`.

**The rule now:** first row carrying one strictly contiguous run of ink spanning
≥72% of the bar width AND containing both side-stroke centres, sustained 2 rows
(one native pixel of the 300 dpi source). A dot cluster cannot satisfy it — it
cannot reach both strokes without a gap — while the bar's own top rule can,
because it is the thing that joins them.

**A trap that the mask change creates, and that I verified before applying:**
`digitize_fig3c.py`'s search start was `int(Y_1E3 - 4*PX_PER_DECADE) + 4 = 2311`,
and the panel box's top rule occupies **y = 2301-2314** — I confirmed 14
consecutive rows at >50% ink across the panel. That was harmless on a blue mask
because the box is black, but under a luminance mask **every Fig 3c bar would
return 2311**. Switching from hue to luminance makes black furniture visible, so
the search window has to be tightened rather than inherited: `TOP_SEARCH = 2320`
in both files, plus a guard assert that the first scanned row is not already a
full-width rule.

**Verification:** 58 of 60 bars agree with an independent hand measurement to
≤2 px. Two differ by 4 px, both the threshold-sensitive bars from a sweep of the
one free parameter (`INK_THR = 172`, mid-plateau: over [166, 178] no bar moves
more than 4 px; it breaks below ~160 and above ~180 ALC-0315 in Fig 3c jumps
back onto the dot answer). `T3A-L2-C8` is immaterial at 1.0% in value;
**`T3A-L2-2C10` is 6.2%**, because Fig 3c has only 151 px/decade, and its honest
value is y = 2521 ± 2 (374k-449k RLU) — flagged as not exact. Either way the
previous 2505 was wrong by 14-18 px. Rule thicknesses come out 6-14 px (2c) and
2-10 px (3c), consistent with a stroked rectangle outline, which is itself a
check that the detector is finding rules rather than dots. The CSVs now carry
`rule_thickness_px` and `frac_used`; the fallback coverage fractions never fire,
and `frac_used` exists so that if one ever does the bar surfaces as needing
review instead of being silently rescued.

**What the correction does to the cross-library offset — the substantive change.**
The stored offset was never a property of the two screens. ALC-0315 was
dot-contaminated in **both** panels (+24 px in 2c, +16 px in 3c) and Lipid 5 in 3c
(+8 px), and those are three of the four benchmarks the offset is measured from.
So the offset was mostly measuring our own detector error:

| | before | after |
|---|---|---|
| cross-library offset | 1.1069 ± 0.0233 | **1.0071 ± 0.0247** |
| mean z, library 1 / 2 | −0.1566 / +2.0885 | −0.1583 / +2.1106 |
| geometric-mean RLU, lib 1 / 2 | 3,569 / 154,630 | 3,556 / 144,926 |
| geometric-mean ratio | 43.3× | **40.8×** |

`rms Δz 0.0188, max |Δz| 0.1419`, 31 rows move >0.02 z and 12 move >0.05 z.

**The two screens agree to within ±2.5%, which is the spread of the four
benchmarks themselves — the offset is not distinguishable from 1.0.** The
division by it is kept rather than dropped, because 1.0071 is the honest estimate
and hard-coding 1.0 would assert an agreement the figures do not quite show; the
sd is reported beside it so a reader can see it covers 1.0. `build_rows.py`'s
prose, which asserted "It is a SCALAR ~1.1x", is rewritten — that claim was
false the moment this landed.

For scale: **F-019 remains the largest of the three effects.** It alone moved
library 1's geometric mean from 4,128 to 3,569 RLU and the ratio from 37.5× to
43.3×.

Also corrected: `digitize_fig2c.py`'s QC `anchors` string still quoted the stale
pre-F-017 "label centroids 1e6=2342…", and `06_fig3c_qc.json`'s
`resolution_note` described a borderline case at 701,135 RLU from a superseded
calibration that no longer exists. Both rewritten. The replacement note says
explicitly that the printed Fig 3d/3e rates **reproduced on the old values too**,
including on 15 bars up to 28% high, so they are not evidence the new rule is
right.

**One recommendation from the re-measurement was withdrawn by its own author**
and not applied: that `digitize_fig2b.py`'s LUT window ran 6 px into the border
skirt. `np.arange(BAR_TOP + 18, BAR_BOT - 17)` ends at 1955, exactly where the
ramp ends, so there was nothing to fix.
