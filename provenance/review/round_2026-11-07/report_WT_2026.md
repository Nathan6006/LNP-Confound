# Round-2 adversarial review — `WT_2026` (Toh, Cheng et al., FALCON, bioRxiv 2026)

*Saved by the orchestrator from the reviewer's returned text: the harness blocked
the agent from writing it. Read-only review; nothing in the repo was modified by
the agent. Findings are CLAIMS — the adjudication of each is in `FINDINGS.md`.*

## 0. Headline — F-004 recommendation

**The percent block is authoritative. The design block is the one out of
alignment. Do NOT reconstruct `N:P` through the permutation — BLANK it on the
129 URS rows.**

| conclusion | confidence |
|---|---|
| the four molar percents on the 129 URS rows are correct as stored | **certain** |
| `IL_to_nucleicacid_chargeratio` is not reliably the N/P of the formulation on its row | **certain** |
| the correct N/P is recoverable by the permutation | **rejected** — three mutually inconsistent candidates |

The halves separate because the design block's *shape* parameters are provably
permuted (they derive other rows' percents, injectively) while `NP_ratio` sits in
the same column block with **no observable saying whether it moved with them**.
On the one URS formulation the paper tabulates, the three candidates are 7.199
(stored), 16.766 (reconstructed), 6.0 (Suppl Table 5).

## 1. Confirmed correct

**F-003's repair holds**, re-derived from `data/` rather than from the curation.
All 2,544 rows invert to `nLnLE_*` at max residual 8.9e-16–2.7e-15 per block; all
2,544 four-molar-ratio sets match the source to 0.00e+00; per-`Model_type`
mean/sd are *not* 0/1 (C2C12 −0.333/0.913, DC2.4 +0.441/0.926, HepG2
+0.471/0.996, NIH-3T3 −0.359/0.951, Ramos −0.290/0.883, THP-1 +0.200/0.999;
INS-1E exactly 0/1 because E3 has one cell type). Selectivity identity holds:
file 5 mean z(HepG2) − mean z(Ramos) = +1.4422 × sd 2.7972 = **+4.0342 ln**
(56.5×), `max |(z_a − z_b)·sd − printed Selectivity| = 1.78e-15`, Spearman(raw,
stored) = **1.000000**, **0** of 180 sign inversions (round 1 measured 92).
**F-003 can be closed.**

**Design→percent conversion, six of seven files clean** at 0.15 pp: URS 0/129;
E0 273/273; E1 197/197 (median 0.040, max 0.110 = source rounding); E2 136/136;
E3 123/123; E4 180/180; E5 180/180.

**URS cross-match:** 126 design rows match exactly one percent row, **0** match
their own, 3 match none, mapping **injective**. N/P error median **6.856**, mean
**7.703**, max **21.535**, 81/126 off by >5 on a 4.01–25.54 range. Every F-004
number reproduces to three decimals.

**Structures:** the seven component names each map to exactly one SMILES
corpus-wide. Both C10 anchors are byte-**different** from the corpus strings
(different atom ordering) and InChIKey-identical — genuinely independent, not the
F-013 failure mode. MC3 is the linear 4-(dimethylamino)butanoate ester, not the
KC2 dioxolane, so the 641.6111 Da trap is not live here.

**Metadata** re-read from Methods: doses 150/25/5/10 ng/well by cell line and
10/2.5 for Fig 6 (= experiment 2 per Suppl Table 9) reproduce the stored µg
values including the ng→µg conversion; FLuc mRNA TriLink CleanCap 24 h; mixing is
MANTIS + Opentrons Flex (`handmixed` is the right bucket in the corpus
vocabulary); aqueous buffer genuinely not stated anywhere in main text or the
62-page SI. Molar ratio sums 99.900–100.100, 0 rows outside ±0.11. Suppl Table 9
library sizes 81/108/36/108/108 still reproduce. All 15 `F####` rows of Suppl
Table 5 reproduce both N/P and the four percents to the printed decimal.

## 2. Errors found

### E1 — F-004: percent block authoritative; `N:P` unrecoverable, not reconstructible
129 rows (`Formulation_ID LIKE 'WT_2026_U%'`), column
`IL_to_nucleicacid_chargeratio`. **Severity:** changes a value (a formulation
descriptor, not an endpoint). **Confidence:** certain on direction; certain that
reconstruction is unsafe.

**(a) `PEG_PEG+Chol` in the URS design block is strictly monotone decreasing over
all 129 rows, zero violations** — while nothing else in the file is sorted
(`NP_ratio` 64/64 asc/desc violations, `IL+HL` 64/64, `HL_IL+HL` 67/61, the four
percent columns ~60/68, `nLnLE_RAMOS` 62/66). A sorted column cannot arise from
uniform random sampling (p = 1/129!). No other file has a sorted design column.
This is the signature of a sort on a contiguous column selection —
`NP_ratio, IL+HL, HL_IL+HL, PEG_PEG+Chol` are columns F:I with the key last. The
alternative (percents moved) would require a sorted percent column; the strongest
|Spearman(percent, row index)| is 0.812, not 1.000. **Identifies which block
moved with no appeal to any external source.**

**(b) The paper ties ID to response.** Main text p. 10: *"F0162 achieved 2-fold
greater nLE than the top unguided search composition (U0036)"*; Fig 2i plots
U0036 ≈ 195 vs F0162 ≈ 390. In the CSV the `FALCON_ID = U0036` row carries the
**file maximum** `nLnLE_RAMOS` = 8.3756.

**(c) Suppl Table 5 ties ID to the percent block.** It is on **page 8** of
`media-1.pdf`, **not p. 59** as round 1 recorded, and was extracted with
`pdfplumber` word coordinates rather than reading order. Its `U0036` row prints
IL% 50.8 / HL% 2.0 / Chol% 43.0 / PEG% 4.2 (sum 100.0); the CSV **percent block**
row 36 is 50.8371/1.97161/43.0412/4.15005 — exact to the printed decimal. The
**design block** row 36 would give 12.72/20.71/62.31/4.26. (b)+(c) ⇒ percent
block ↔ response, so the curation's molar percents are right. Corroborated by
Suppl Table 7, which reports U0036's DLS (103.3 ± 2.2 nm, PDI 0.13, EE 99.8%) —
it was physically re-made at that composition and dosed into mice (Fig 2j).

**(d) A second independent tie.** URS percent row 1 (`U0001`) is byte-identical
to E0's `F0162` percents (38.687134/14.016387/46.780892/0.515587). In clean E0
that composition carries design params (52.7035, 26.5948, 1.0901) and N/P 25.056,
and Suppl Table 5's `F0162` row independently prints 52.7/26.6/1.1 and N/P 25.1.
The URS design block on that row holds (26.609, 3.998, 8.9516) and N/P 24.885.

**(e) The obvious statistical test is confounded — do not use it.** A naive
"which alignment predicts the response better" gives *both* alignments signal (RF
5-fold Spearman 0.318 vs 0.363, null 0.00 ± 0.11), because the response
correlates with row index (ρ = −0.370) and the design block's `PEG_PEG+Chol` is
*perfectly* anti-correlated with row index by (a), so the design alignment
inherits the index artifact. A cleaner external-model transfer (train on E0's 273
MC3/DSPC Ramos rows) favours the percent alignment: RF ρ = +0.194 / GBM +0.293
vs −0.008 / +0.008. **Corroboration only**; the case rests on (a)–(d).

**Why `N:P` must be blanked, not reconstructed.**
**(f)** Three mutually inconsistent candidates for U0036: CSV `NP_ratio` row 36 =
**7.1985**; reconstruction through the injective permutation (design row 3) =
**16.7659**; Suppl Table 5 = **6.0**. No URS row has `NP_ratio` within 0.02 of
6.0 (closest 6.0461); no row in any of the seven files does.
**(g)** The injective match uses only the three *shape* parameters. Whether
`NP_ratio` moved with them is an assumption about whether the Excel selection was
F:I or G:I, and nothing in the data distinguishes them — if G:I, the stored
values are already correct.
**(h)** Neither alignment carries signal: Spearman(stored N/P, `nLnLE_RAMOS`) =
**+0.042** on the 126 matched rows, reconstructed **−0.041**, 95th percentile of
|ρ| under random pairing (2,000 draws) = 0.167.
**(i) Suppl Table 5 is itself block-pasted, so it cannot rescue the
reconstruction.** The `U0036` row's eight numeric cells read
`[6.0, 50.8, 2.0, 43.0, 4.2, 6.0, 50.8, 2.0]` — a period-5 repeat; the last three
should be 52.8/3.7/8.8 and instead repeat the row's own N/P, IL%, HL%. Confirmed
by word coordinates (x = 349/408/488, the columns every other row uses), so not
an extraction artifact. At the table's bottom the last three columns slip by one
row: `F3017` prints 50.0/31.4/0.6 = **F3118's** true triple, and `F3118` prints
61.3/37.2/4.8, traced to **F3123** in file 3 — a formulation the table never
lists. IDs, N/P and percents are correct on those rows; only the derived block
slipped.

**Proposed repair:** on the 129 `WT_2026_U*` rows blank
`IL_to_nucleicacid_chargeratio`, leave the four molar ratios untouched, and write
a `curation_note`. Do **not** blank the percents and do **not** reconstruct N/P.
**How to verify it is the right half:** afterwards Suppl Table 5's `U0036`
percents must still equal the stored molratios on `WT_2026_U0036`, and that row's
`Experiment_value` must still invert to 8.3756 (the file max the main text calls
the top URS composition). Neither involves N/P.
**Error-class candidate:** *a source stating the same quantity in two
parameterisations where only one is row-aligned.* Automatable wherever both are
given — derive one from the other and require row-wise agreement; a cross-row
injective match with zero self-matches is the positive signature.

### E2 — `Publication_PMID` blank on all 2,544 rows; indexed since 2026-06-29
**Severity:** metadata. **Confidence:** certain. PMID **42367940**, bioRxiv,
`doi 10.64898/2026.06.13.731774`, entrez 2026/06/29, authors Toh WH, Cheng L, …
The DOI matches `Publication_DOI` exactly. The record was live **four months
before** this curation (2026-10-29) and before round 1. `audit_baseline.json`
accepted 2,544 blanks with the reason *"was not in PubMed as of 2026-10-30"* —
**that reason is factually false**, so the gate was accepting a real defect on a
wrong justification. Verbatim recurrence of the `QL_2026` lesson already written
into CLAUDE.md §3, on the same §3 row that records the previous instance.

### E3 — URS percent rows 1–3 are not uniform-random samples (NEW)
`WT_2026_U0001/U0002/U0003`. **Severity:** metadata, possibly a value
assignment. **Confidence:** certain on the observation, speculative on the cause.
The injective cross-match leaves exactly three unmatched rows on each side, and
on the percent side they are **rows 1, 2, 3** (1 in 357,760 under random
placement). `U0001` percents = `F0162`'s byte-identical to six decimals;
`U0002` = `F0192`'s byte-identical; `U0003` = 49.743590/10.256410/38.461538/
1.538462 = exactly 1940:400:1500:60 over 39, giving `(IL+HL)% = 60.0000` exactly
and `PEG/(PEG+Chol)% = 3.8462` = 1/26 — a round *designed* composition, close to
but not the Onpattro reference. The file name ends `_Ctrl_MC3_DSPC`. Their
`nLnLE_RAMOS` (4.776, 6.982, 6.782) sit above the file mean 4.671, so they are
not inert in any URS-vs-FALCON comparison, and they are the three rows for which
N/P is not even nominally recoverable.

### E4 — `Experiment_value_within_model` silently pools five incompatible normalisation units (NEW)
**Severity:** changes a value for any consumer of that column. **Confidence:**
certain. `WT_2026` is the corpus's only `value_normalization_unit = source_file`
experiment, because `dataset_info.txt` forbids comparing across the seven files.
But `Experiment_value_within_model` is defined within (`Experiment_ID` ×
`Model_type`), and one `Model_type` spans several source files — `Ramos` alone
spans URS + E0 + E1 + E4 + E5 = 959 rows:

| file | n | `Experiment_value` mean/sd | `..._within_model` mean/sd |
|---|---|---|---|
| URS | 129 | +0.000 / 1.000 | **+0.329** / 1.133 |
| E0 | 273 | −0.395 / 0.887 | **−0.119** / 1.005 |
| E1 | 197 | −0.434 / 0.710 | **−0.164** / 0.805 |
| E4 | 180 | +0.251 / 0.953 | **+0.613** / 1.080 |
| E5 | 180 | −0.721 / 0.450 | **−0.488** / 0.510 |

`Experiment_value`'s spread across files is expected and *declared*.
`Experiment_value_within_model` manufactures a new normalisation unit straddling
five scales the source says are incomparable, and **no column declares it** —
`value_normalization_unit` describes `Experiment_value` only. A consumer who
picks the within-model column to remove an assay offset gets exactly the
cross-file mixing the `source_file` flag exists to prevent.

### E5 — quality-flag columns blank on all 2,544 rows, where at least three conventions are in force
**Severity:** metadata. (1) **`18PG` (1,114 rows) is an inference from a
different paper** — FALCON says only *"18PG (sodium salt) … from Avanti"* and
never expands it; the corpus's SOPG comes from `YZ_2022`'s definition.
Corroboration the curation did not cite: the repo's own
`falcon_engine/run_mantis_formatter.py` prompts for a helper lipid from
`DOTAP, DSPC, 18PG, DOPE, DDAB, 14PA, 18MP`, the exact shorthand set of that
literature. Still an inference. (2) **Experiment 2's N/P is not a pure
ionizable-lipid N/P** (272 rows) — Suppl Table 3: *"numerically adjusted to
account for the number of nitrogens from both SM102 and DSPC"*; `PROVENANCE.md`
records it, the rows do not. (3) `Mixing_method = handmixed` is a vocabulary
mapping for robotic bulk dispensing.

### E6 — `PROVENANCE.md` still describes the pre-F-003 curation (error class 21)
Seven stale assertions, each contradicted by the folder's own regenerated
artifacts: "14 experiment blocks" (×5), "2,544 × 72 cols" (×2), and a corpus
total of "25,997 → 28,541 rows, 40 → 41 experiments".

### E7 — minor: `validate.py`'s N/P check cannot see E1
`04_validation.json` records `PASS N/P range 3.33-37.12`. A range check passes
identically on 129 misassigned values, because the permutation preserves the
multiset. Noted so the pass is not mistaken for coverage.

## 3. Unverifiable
1. The true N/P of the 129 URS formulations — three candidates, nothing selects
   one. Would be settled by the authors' MANTIS formulation sheet for the URS run
   (the repo ships only a blank template and a demo sheet filled with E0).
2. Whether `NP_ratio` was inside the sorted column selection. Not observable.
3. Why URS percent rows 1–3 are comparator compositions.
4. Whether `18PG` here is SOPG — strong circumstantial support, not proof. An
   Avanti catalogue number would settle it.
5. DMG-PEG2000's polydispersity stored as one SMILES (corpus-wide convention).
6. Suppl Table 5's `U0036` N/P of 6.0 — paste artifact or an as-made value.

## 4/5. Severity and confidence
| # | changes | confidence |
|---|---|---|
| E1 | a value (129 rows); no endpoint moves | certain on direction and on non-reconstructibility |
| E2 | metadata (2,544 rows + a false baseline `why`) | certain |
| E3 | metadata, possibly a value assignment (3 rows) | certain on observation, speculative on cause |
| E4 | a value for consumers of `Experiment_value_within_model` | certain |
| E5 | metadata (flagging posture) | certain on blanks, probable per flag |
| E6 | cosmetic for data, material for the next reviewer | certain |
| E7 | cosmetic | certain |

**No error found changes a structure, and none changes an `Experiment_value`.**
All 2,544 endpoint values and all 2,544 × 4 molar ratios round-trip exactly.

## 6. Not covered
`Experiment_value_within_model` for other experiments; `curation/ng_2025/` (not
read, per instruction). The reviewer did not re-run `build_rows.py`,
`validate.py`, `build_dataset.py` or `audit.py` — every number was derived from
the source CSVs, the PDFs, NCBI, or `data/data_vitro.csv` by its own scripts.

## 7. Corrections to round 1's record
- F-004's anchor is Supplementary Table 5 on **page 8** of `media-1.pdf`, not
  p. 59. *(Orchestrator's note: p. 59 was my error, not round 1's — I passed the
  agent a line number from a `pdftotext` dump as if it were a page number.)*
- Suppl Table 5 does contain a `U0036` row and it **does** establish the percent
  block's alignment — but its N/P cell (6.0) matches neither candidate and the
  same row's last three cells are demonstrably corrupt, which is why the repair
  is "blank", not "reconstruct".
