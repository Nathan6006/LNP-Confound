# Dataset provenance

Everything a reviewer needs to trace `data/data_vitro.csv` and `data/data_vivo.csv`
back to their sources.

| | |
|---|---|
Built | 2026-10-28 |
Builder | `build_dataset.py` (rerun it to reproduce byte-for-byte) |
Machine check | `dataset_validation.json` |
Change log | `CHANGELOG.md` |
Literature sweep behind the candidate list | `../lit_review/provenance/` |
Per-paper figure extraction | `../curation/<paper>/` |

## The two files

| file | rows | cols | contents |
|---|---|---|---|
`data/data_vitro.csv` | 25,997 | 72 | 40 experiments, all in vitro |
`data/data_vivo.csv` | 1,768 | 72 | 18 experiments, all in vivo |
`data/Original_LNPDB.csv` | 19,797 | 66 | untouched upstream reference |

`data_vitro.csv` = original LNPDB (vitro) + the ML-work change log + the Bowen Li
lab additions, with three deletions reversed, plus two curated papers
(`ZT_2026`, 486 vitro / 51 vivo; `LL_2026`, 840 vitro; `LL_2026_2`, 147 vitro) — see CHANGELOG.md §6. `data_vivo.csv` is the vivo corpus
re-canonicalized with descriptors recomputed; the ML change log never touched it,
and there are no Bowen Li vivo additions.

## Lineage, and which source wins where

Two upstream files disagree, so the builder is explicit about precedence.

**P** — `LNP_ML/new_data/LNPDB_vitro_del_processed.csv` (23,194 × 83). The
post-change-log dataset. **Authoritative for `IL_SMILES` and `Experiment_value`**,
because it carries the manual corrections. It lacks 21 curation columns.

**V** — `data/LNPDB_vitro.csv` (19,028 × 67). The curation schema.
**Authoritative for those 21 columns**, and the only source for the restored rows.

Row assembly:

```
every P row except YX_2024                21,994
all of YX_2024, taken from V               2,400   (restores 1,200 RAW264.7)
ZC_2023 from V                               100   (restored)
XH_2024 from V                                30   (restored)
                                          --------
                                          24,524
```

Taking YX_2024 wholly from V was **checked, not assumed**: V's z-score over all
2,400 rows is bit-identical to P's over its 1,200 (affine fit slope 1.000000,
intercept 0.000000, max residual 7×10⁻¹⁵), so P never re-z-scored after dropping
RAW264.7 and the restored experiment is already correctly scaled. P had corrected
20 YX_2024 structures; those corrections are carried onto the restored RAW rows
through a V-canonical → P-canonical table, which is sound because RAW's SMILES set
is a strict subset of HeLa's.

## Columns

The 67 curation columns, then 5 RDKit descriptors from the change log:
`Num_carbon_in_tail`, `num_protonatable_nitrogens`, `num_permanent_cationic_N`,
`formal_net_charge`, `num_unsaturated_cc_bonds`.

All five definitions were **reverse-engineered and verified against the change-log
file** rather than guessed — each reproduces its column on 1,200/1,200 unique
structures. Two were not obvious:

- `num_protonatable_nitrogens` counts neutral non-amide amines **plus
  pyridine-type aromatic N**. The aromatic clause matters: imidazole-headed
  lipids carry one basic N and one non-basic pyrrole NH, and omitting it
  mis-scored every one of them.
- `Num_carbon_in_tail` is the **longest path through acyclic carbons**, not a
  parse of the lipid name (it merely coincides with the name suffix in the
  lipidoid libraries).

**One-hot columns are deliberately not stored.** The change log one-hot encoded
`HL_name`, `Model_type` and `Cargo_type`, but restoring ZC_2023, XH_2024 and the
RAW264.7 rows changes those category sets, and the Siegwart paper adds four more
helper lipids. Any one-hot baked in now would be stale on the next append.
Generate at model time:

```python
import pandas as pd
df = pd.read_csv("data/data_vitro.csv")
X = pd.get_dummies(df, columns=["HL_name", "Model_type", "Cargo_type"])
```

`rel` and `Experiment_weight` are likewise not stored — they are modelling
artefacts, and `rel` is undefined for the restored ZC_2023 anyway.

## `ExactMolWt` was renamed

The column was called `van.der.Waals.Molecular.Volume` and does not hold a van
der Waals volume — it holds the monoisotopic exact mass. Verified on
**40,472 of 40,542 rows** across all three upstream files (12,846 unique lipids).
The 70 exceptions are confined to `JW_2024` and `LX_2024` and are stale rows:
their `Molecular.Weight` disagrees with their SMILES too, and the
stored-exact / stored-average ratio is 0.99928 ± 0.00002, the monoisotopic
ratio. Renamed throughout; `.bak_rename_20261028-214051` backups were kept.

## Every descriptor was recomputed

Descriptor columns are computed fresh from the canonical `IL_SMILES`, not carried
over. **676 rows** had descriptors disagreeing with their own structure in the
curation file — the 70 stale rows above plus drift wherever P corrected a
structure without V being updated. An independent recheck of 4,000 random rows in
the finished file found 0 mismatches.

## Validation performed

All pass; see `dataset_validation.json` for the full machine output.

| check | result |
|---|---|
Row count | 24,524 = 21,994 + 2,400 + 100 + 30 |
Restorations | ZC_2023 100/100, XH_2024 30/30, YX_2024 RAW264.7 1,200/1,200 |
Retained deletions | ZL_2022 + YX_2023 absent; BMDC/BMDM/HBEC_ALI absent |
z-score per experiment | all 37 at mean 0, sd 1 — none off spec |
`IL_SMILES` parse | 24,524/24,524 parse; all canonical |
iPhos (SL_2021) neutral | formal charges {0: 494, +1: 78} — no anionic forms |
Lee_Unsat (SL_2020) E/Z | 69 trans motifs present (V had 0) |
T/C duplicate SMILES | 28 collisions in V resolved; 1 remains (see below) |
Descriptors | 676 corrected; 0/4,000 mismatches on recheck |

## Normalization convention was not uniform — now it is

`Experiment_value` is a per-experiment z-score, but eleven experiments were not
on the corpus's population-sd convention. Five (`SS_2018`, `LL_2025`, `SX_2025`,
`BL_2023_3CR`, `KL_2010`) were z-scored with the **sample** sd (ddof=1); six
others (`SP_2020`, `LZ_2025`, `YZ_2022`, `YZ_2024`, `AC_2025`, `YX_2024`) sit
slightly under, consistent with rows being dropped *after* normalization.
Deviations ran from 0.01% to 3.4%.

Every experiment already centred at zero was rescaled to exact population sd 1.
That is a pure affine map **within** an experiment, so it preserves every
ordering and relative spacing — a per-experiment z-score carries no more
information than that — while making the convention uniform. All 52 experiments
across both files now sit at mean 0, population sd 1 to within 1e-9.

## `SB_2024` and `SY_2025` are not z-scores at all

Flagged, **not** fixed and **not** dropped — this is your call. Both use
`Experiment_method = luminescence_relative_to_Spikevax`: ratios to a benchmark,
floored at 0, with mean ≈ 0.44 and sd ≈ 0.4. Forcing them to mean 0 / sd 1 would
disguise the problem rather than surface it.

They are the vivo analogue of the `ZL_2022` / `YX_2023` problem that got those
two dropped from the vitro corpus for being on a non-comparable scale, so the
same decision probably applies. 60 rows between them.

## Structure patches

Two collisions were traced to their source and fixed. In both cases the stored
SMILES was a sibling's, which is always a transcription bug: two names that
differ in tail length cannot share a structure.

**`SL_2020` / `4A3-SC6`** (1 lipid, 1 row). Held the SC8 octyl structure. The
SL_2020 tail series is 6-carbon — `2C`/`3C`/`4C`/`5N` name the double-bond
position — so `SC6` is the saturated hexyl control. Every other SC6/SC8 pair in
SL_2020 is already correct and differs by exactly 2 CH₂ per tail, which pins the
expected mass: 4A3 has 4 tails, so SC6 = 1467.21 − 112.22 = **1354.99**. The
patched structure hits 1354.99 exactly.

**`KZ_2016` / `*-SC3`** (42 lipids, 42 rows). Every SC3 held the SC2 isopentyl
structure. `LNP_ML/data_files/Zhou_dendrimer/individual_metadata.csv` resolves
what SC3 should be: `6A3-SC2` and `6A3-SC3` share MolWt 1876.692 but are recorded
with **different tail chain lengths, 4 and 3**. So SC3 is a more branched isomer
of the same C₅ formula. With SC1 = isobutyl (C4, chain 3) and SC2 = isopentyl
(C5, chain 4), **SC3 = neopentyl** (C5, chain 3). Verified: formula identical to
SC2 on all 42, structure now distinct on all 42, and the longest-chain descriptor
drops 4 → 3 on the 33 lipids whose head group does not itself contain a longer
chain (the other 9 split identically for SC2, so the pattern is consistent).

## The 328 `KL_2010` rows: answered

They are **`NA` in `Original_LNPDB.csv` itself** — genuinely absent upstream, not
lost in this pipeline. The curation file merely re-encoded `NA` as an empty
string. `KL_2010` is Love et al., *PNAS* 2010 (PMID 20080679, the C12-200 paper):
LNPDB catalogues all 378 synthesized lipids but only 50 carry a quantified HeLa
siRNA value, so the remaining 328 are structures without a measured endpoint in
this assay. The experiment keeps its 50 valued rows.

While checking this, one other experiment turned up with the same condition:
**`BroadPharm_2024` has 269 rows in `Original_LNPDB.csv`, every one `NA`.** It is
absent from the curated corpus entirely, which is correct — there is nothing to
model — but it is recorded here so the 19,797 → 19,028 gap is accounted for.

## Duplicate SMILES that remain — flagged, not patched

After the two patches, 48 collisions remain, and they are a different shape: not
tail-length errors but **head-group collisions**, where two library IDs resolve to
the same structure across an entire tail series.

| experiment | colliding head pairs | rows |
|---|---|---|
`KW_2014` | A165/A69, A191/A192, A163/A1, A103/A167, A164/A31 | 25 |
`BL_2023` | A55/A57, A24/A25 | 11 |
`AA_2008` | 21/22 | 10 |
`JW_2024` | M3/M81, M80/M95 | 2 |

Because each pair repeats across every tail length in its series, these are
consistent assignments rather than one-off typos — which means they are either
genuine (the same amine entered the library under two IDs) or a systematic
head-group mix-up. Deciding which needs the source publications, so they are
**left untouched and flagged** rather than guessed at. Note also that several
`BL_2023` `IL_name` values carry a trailing space.

## Upstream reference verified against a fresh download

A fresh copy of LNPDB was pulled from the official site on 2026-10-28 and is
**byte-identical** to the reference already in the repo
(MD5 `32b5187361fb1164622a9f6a04471d40`), archived as
`superseded/LNPDB_official_download_20261028.csv`. `data/Original_LNPDB.csv` is
that same file with one change: the `van.der.Waals.Molecular.Volume` → `ExactMolWt`
rename.

This settles two things. The `NA` values are in the official release, not
introduced here — so the 328 `KL_2010` and 269 `BroadPharm_2024` blanks are not
recoverable from upstream. And the known upstream defects (the SC6/SC8
collision, the mislabelled volume column) are genuinely upstream.

## Independent audit

`audit.py` → `audit_findings.json`. Written to look for error classes the
builder does not know about, so a clean build and a clean audit are separate
statements. Nine checks, both files.

| check | vitro | vivo |
|---|---|---|
C1 SMILES parse / canonical / no disconnected fragments | 0 / 0 / 0 | 0 / 0 / 0 |
C2 mass monotonicity in verified-length series | 0 violations | 0 violations |
C3 duplicate SMILES within an experiment | 48 name-groups (39 with conflicting values) | 12 (12 conflicting) |
C4 exact duplicate rows | **0** | **0** |
C5 formulation ratios | 0 off | 3 (explained below) |
C6 non-finite or extreme values | 0 | 0 |
C7 field consistency / blank required fields | 0 | 0 |
C8 descriptors vs stored SMILES | 0 of 24,524 | 0 of 1,717 |
C9 unassigned stereogenic C=C | 212 rows | 2 rows |

**Four apparent problems turned out to be my own checks being wrong**, and are
recorded because the corrected scoping is what makes the clean result meaningful:

- *Mass monotonicity* first reported 2,076 violations. The numeric token in a
  lipid name means a chain length in only some libraries. KS_2024's `C24` is
  **two C4 tails**, KZ_2016's `SC1/SC2/SC3` are tail **identifiers** (isobutyl /
  isopentyl / neopentyl), YX_2024's `C1/C2/C3` are component **indices**, and
  BL_2023's `T1/T2/T3` are tail identifiers where mass legitimately falls as
  unsaturation rises. The check now runs only where the convention is verified
  against source metadata (KW_2014: 1,180 of 1,182 names agree; the 2 exceptions
  carry a dodecyl on the *head*, which outranks a shorter tail). 0 violations.
- *Contradictory replicates* first reported 1,001 groups. 746 of them simply
  differed in `IL_to_nucleicacid_massratio`, a real experimental condition the
  key omitted. Keyed on all 68 descriptive columns: 333 true replicate groups
  (median spread 0.38 z-units, mostly the AC_2025 formulation sweep) and
  **zero exact duplicate rows**.
- *Stereochemistry* first reported 334 rows. SL_2020's 21 are terpenoid tails
  (citronellyl / farnesyl / nerolidyl) whose trisubstituted alkene carries two
  methyls on one carbon — there is no E/Z there to lose. A further 240 in
  BL_2023_3CR are C=N imines. Using `FindPotentialStereoBonds` and separating
  C=N leaves 212 genuine cases.
- *`4A3_C10`-style name checks* flagged 3,917 rows before the same naming issue
  was accounted for; the residue was 2 rows, both explained.

## Open findings from the audit

**C9 — 212 rows with an unassigned stereogenic C=C.** FG_2026 (120) and SS_2018
(2) carry no C=C stereo anywhere, so that is a source convention rather than an
inconsistency. The sharper cases are the **26 lipids that mix assigned and
unassigned C=C inside one molecule** (YX_2026_SDL 23, NC_2024 3). Inspected: in
every case the unassigned bond sits in a *different structural motif* from the
assigned ones — a butene linker in NC_2024's head group, an α,β-unsaturated
amide in YX_2026_SDL — so there is no correctly-assigned sibling to copy from.
Left unassigned rather than guessed. Worth knowing if you use stereo-aware
fingerprints.

**C3 — duplicate SMILES, now partly explained.** The colliding KW_2014 heads
(A165/A69, A191/A192, A163/A1, A103/A167, A164/A31) are a *known* upstream
condition: `LNP_ML/data_files/Whitehead_siRNA/` keeps only one of each pair and
its note says "removed duplicate smiles (for each unique pair, removed all but
the first time it happened)". So another curator hit the same thing
independently and chose to deduplicate. LNPDB retains both IDs with different
measured values. Whether the structures or the IDs are wrong cannot be settled
from here, so both are kept and flagged.

**C5 — 3 vivo rows off 100%.** Two are `SB_2024` SORT formulations where
6,6′-trehalose dioleate is added at 50 mol% *on top of* the four-component base,
giving 125 — that is the SORT convention, not an error. The third sums to 101.0
with no fifth component, a rounding artefact.

## One thing left open, deliberately
