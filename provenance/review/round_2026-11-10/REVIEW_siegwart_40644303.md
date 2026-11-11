# Wave-1 adversarial review — `curation/siegwart_40644303/`

Reviewer: read-only Wave-1 agent. Date 2026-11-10.
Target: Lee SM et al., *ACS Biomater. Sci. Eng.* 2025, 11, 4844–4852, PMID **40644303**.
Scope: 505 emitted rows, 489 enumerated structures, the claimed 360-lipid corpus overlap,
the row-count discrepancy the packet reports, and the folder's own open decisions.

**Headline.** I could not break a single structure or a single value in this folder.
489/489 structures reproduce by full InChIKey against an independent enumerator I wrote,
and 465/465 heat-map cells reproduce against an independent pixel read with my own grid
detection and my own clustering. What I did break is **the corpus**: I resolved this
folder's open question D-4 against `data/data_vitro.csv`, and confirmed D-6 against it.
Both are live errors in `KZ_2016` today, on 78 rows, and this curation is right on both.

The project-level finding is in §2: the 360-lipid overlap is real chemistry but it is
**100% within the Siegwart laboratory**, 97.6% of it one prior library, and it is
**cross-cargo, not cross-laboratory**. Merging it as a Stage-5 pair would move Result 3's
headline by a confound.

---

## 1. CONFIRMED CORRECT

### 1.1 All 489 structures — independent re-derivation, full InChIKey

I wrote my own enumerator (`indep2.py`) using **pure f-string SMILES assembly**: no
`RWMol`, no `AddBond`, no isotope tags, no `curation/lib`, no import of any file in the
folder except the finished `03_enumerated.csv` for comparison. Core templates were typed by
hand from my own reading of Figure 1b and Scheme 1; the arm was written from the two-step
chemistry in the Methods, not copied from theirs.

| | |
|---|---|
| Full-InChIKey agreement | **489 / 489** |
| Formula disagreements | 0 |

My first pass scored 480/489. The 9 misses were cores `3A1`/`3A2`/`4A1` (ring cores) ×
tails `SN4`/`SN5`/`SN6` (ring tails), identical formula in every case. **The defect was
mine**: my template reused SMILES ring-closure digit `1` for the tail ring while the core's
`1` was still open, so RDKit paired the arm's digit with the core's. This is error class 18
(the piperazine/cholesterol ring-closure collision) reproduced live in a reviewer's own
check. Renumbering the three tail rings to `2` took it to 489/489. Recording it because the
near-miss is the evidence that the comparison had teeth: a wrong structure in their file
would have surfaced the same way.

Confidence: **high**. Caveat in §3.1 — my building blocks and theirs both come from reading
Figure 1b, so this tests the *enumerator*, not the blocks.

### 1.2 The arm's direction, read off Scheme 1

The brief's central warning ("a SMILES fragment has a DIRECTION that mass cannot see") is
the live hazard here, since one directional error multiplies across 489 structures. I
rendered Scheme 1's drawn product at 6× and read the thioether terminus directly
(`s1arm.png`): the drawing is unambiguous —

```
-O-CH2CH2-O-C(=O)-CH(Me)-CH2-S-CH2CH2CH2CH2-CH3
                    α       β
```

**Me on the α carbon, S on the β carbon.** That is exactly `ARM =
"[CH3:1]CC(=O)OCCOC(=O)C(C)[CH3:2]"` with map 2 bonding sulfur. The swapped assignment has
the identical formula and exact mass, so this is the only evidence that distinguishes it,
and the folder cites the right evidence (the written aza-/sulfa-Michael sequence plus the
drawing) rather than a mass. `stage03_enumerate.py`'s header also discloses that an earlier
version of this file bonded the arm through a phantom sulfur and built 465 **disulfides**,
caught only by corpus InChIKey comparison — a disclosure that is itself a good sign.

Confidence: **high**.

### 1.3 The factorial actually multiplies out

The brief asks that a claimed factorial be proved, not asserted.

| claim | check | result |
|---|---|---|
| 465 = 15 cores × 31 tails | set of emitted `IL_name` vs the full Cartesian product | **identical set, every cell exactly once, 0 missing, 0 extra** |
| 42 = 7 cores × 6 tails (in vivo) | same | 42-cell grid, **40 emitted**, missing exactly `3A4V-SC5` and `5A4-SC10` — the two documented in D-3 |

The 15 cores are independently corroborated by the Methods' chromatography list (line 309 of
the fulltext names exactly those 15), and the 31 tails by the Figure 1a column labels, which
I read myself. The paper's own "465" and "42 lipids in total, including the parent lipids"
are reproduced with no residual.

Confidence: **high**.

### 1.4 The 505-row count, and vitro/vivo disjointness

The packet is wrong here; see §2.4. Measured:

| file | rows |
|---|---|
| `09_rows_vitro.csv` | 465 |
| `09_rows_vivo.csv` | 40 |
| **total** | **505** — matches the queue and `09_emit_qc.json` |

Disjointness: overlap on `(Experiment_ID, IL_name, Model_target, Model_type)` is **0**.
17 `IL_name`s appear in both files, correctly — they are the SC5–SC10 parents measured both
in vitro (`Model_target = in_vitro`) and in liver (`Model_target = liver`), i.e. 18 grid
members minus the one dropped. All 505 rows' `IL_SMILES` match the enumeration by InChIKey
(0 mismatches). **The folder's row count is not ambiguous.**

Confidence: **high**.

### 1.5 The heat map — 465/465, from my own grid

I re-read Figure 1a from the raster without using any of their geometry:

- my own band detection (contiguous rows with >600 coloured pixels) → exactly **30 bands**,
  15 siRNA + 15 mRNA, alternating;
- my own column detection → exactly **31 column blocks** in the panel (x 177.5 → 1358.5);
- my own 3-means clustering of sampled patches, seeded by luminance extremes, with no
  reference to their swatches.

| | mine | theirs |
|---|---|---|
| column origin / pitch (px) | 177.5 / 39.367 | 177.381 / 39.3843 |
| siRNA row origin (px) | 135.0 | 134.583 |
| mRNA row origin (px) | 143.5 | 143.729 |
| mRNA band counts | 68 / 72 / 325 | 68 / 72 / 325 |
| siRNA band counts | 35 / 153 / 277 | 35 / 153 / 277 |
| **band agreement** | **465/465 mRNA, 465/465 siRNA** | |

My cluster centroids match theirs to <1 RGB unit. Their swatches are taken from the Figure 1c
legend rectangles — an anchor external to the heat map itself — and cross-checked against the
derived clusters, which is the right posture.

**External ground truth.** I reproduced their Q2 independently: the siRNA stripe against
`KZ_2016`'s stored values gives Spearman **0.8990** over all 465 (they report 0.899), and
every one-step grid shift collapses it (row ±1 → 0.528/0.547, col ±1 → 0.600/0.592,
diagonal → 0.399). The grid assignment is validated against data we did not read from this
figure. See §2.5 for a caveat about what the shift floor means.

Confidence: **high**.

### 1.6 The z-score convention

`Experiment_value` is a z of the **bin ordinal**, per CLAUDE.md §3's rule for a discretized
source. Verified arithmetically rather than taken on trust: with counts 325/72/68 at ordinals
1/2/3, mean = 1.447312, population sd = 0.734643, so ordinal 1 → **−0.608885**, which is the
stored value to all nine decimals. Blocks: vitro `SL_2025` n=465 mean −2.4e−10 sd 1.000000;
vivo `SL_2025` n=17 and `SL_2025_2` n=23, both exact. (Per the brief: mean 0 / sd 1 is a
tautology — what I checked is that the quantity z-scored is the ordinal, not a midpoint.)

Confidence: **high**.

### 1.7 Every printed statistic, including the one that fails

I verified the paper's printed numbers against the primary sources myself, not by rerunning
`stage08_qc.py`.

| printed claim | my independent check | verdict |
|---|---|---|
| "465 lipids" | 15 × 31, grid complete (§1.3) | PASS |
| "42 lipids in total, including the parent lipids" | 7 × 6 (§1.3) | PASS |
| formulation "15:15:30:3 (23.81%: 23.81%: 47.62%: 4.76%)" | emitted `IL/HL/CHL/PEG_molratio` = 23.81/23.81/47.62/4.76 | PASS |
| "25 ng of Luc mRNA per well" | `Dose_ug_nucleicacid` = 0.025 | PASS |
| "E.E. … greater than 80% across all formulations" | SI text table, min **86.00** (3A4V2-SC9) | PASS |
| "pKas of the 3A4 family ranged between 5.9 and 7" | SI, 5.99 (3A4-SC10) – 6.91 (3A4V-SC10) over 18 | PASS |
| "the size … ranged from 100 to 200 nm" | SI, **5A5-SC8 = 248 nm, 5A5-SC10 = 203 nm**; 40/42 | **FAIL — the paper's claim is wrong, and the folder says so** |
| "3A4 … SC6 tail presented the highest liver Luc expression" | my own bar-top read, §1.8 | PASS |
| "highest Luc activity observed in 4A3-SC8 and SC9" | their ranking puts SC9 > SC8 at the top of the 4A3 series | PASS |

I re-extracted SI Table S3 from the PDF text myself. The PDF streams it **column-by-column**,
so the parser had to interleave four separate blocks against a 42-name list — a classic place
for an off-by-one. I checked all 42 rows of the 3A4-family block (diameter, PDI, pKa, E.E.)
and all 24 diameters of the second block against the raw stream: `02_si_table_S3.csv` is
**exact**, including the alignment.

The one FAIL is recorded in `08_contradictions.csv` and is **not excused** — the folder notes
the word "generally" might cover it and reports 40/42 anyway. That is the posture the brief
asks for, and it is the opposite of the past pass that excused the one disagreeing hit rate.

Confidence: **high**.

### 1.8 In-vivo bar digitisation

I independently detected the Figure 2b bars by colour and measured their tops with a
full-width coverage rule:

| x centre (mine) | my top row | theirs |
|---|---|---|
| 273 | 1081 | 1081 (`3A4-SC5`) |
| **332** | **908** | **908** (`3A4-SC6`) |
| 388 | 994 | 994 |
| 444 | 1051 | 1051 |
| 499 | 1050 | 1050 |
| 555 | 1069 | 1069 |

6/6 exact. `3A4-SC6` has the minimum top row, i.e. the tallest bar, confirming the printed
claim. Converting with their tick fit (zero at 1096 px, 1e6 at 841 px): (1096−908)/255 × 1e6
= **737,255**, matching their stored 7.37e5.

Two things done right that I specifically went looking for:

- **The scale comes from axis furniture only.** `px_per_unit` is fitted on the six y ticks
  (residual 0.1147 px), *not* on the quantity being checked. The drawn baseline is measured
  (0.41 px from the fitted zero on Fig 2b, −0.20 px on Fig 3b) and explicitly **reported and
  not used** — this is the "round number the authors drew but never printed" test, used as a
  check rather than as an anchor. That is the correct side of error class 20.
- **`5A4-SC10` is refused, not rescued.** `bars.bar_top` declines it because the left 9 px
  flank has no covered row, which is the F-018/F-020 replicate-dot signature. Lowering the
  threshold to make it pass was explicitly not done.

Confidence: **high** on the 6 bars I re-measured; **medium** on the other 34, which I did not
re-measure (see §4).

### 1.9 The C10 anchor

Scheme 1's drawn product `3A4-SC5` is the only complete lipid in the paper, and the second
row of `reference_structures.csv` records the absence of any other with an **empty
`reference_SMILES`** and a note, rather than leaving the gap invisible. C10 therefore covers
**1 of 489** structures, and the file says so.

On independence: the anchor string is **byte-identical** to the enumeration. Per F-037 that
is *not* the signal and testing for it is the trap — a correct re-derivation canonicalises to
the same string. Judging by provenance instead, as the brief directs: `make_reference.py`
writes an explicit **V2000 atom/bond table** and parses it with `Chem.MolFromMolBlock`. I
read the whole file: there is no SMILES fragment of chemistry, no `RWMol`, no `AddBond`, no
`RemoveAtom`. The enumerator uses `curation/lib/assemble.attach`, i.e. `RWMol` + `AddBond` on
SMILES fragments. These are genuinely different code paths sharing only RDKit's sanitiser.

The file also notices, and documents, that the obvious "independent" route here would have
been a copy: the enumerator's canonical SMILES for this lipid already starts at the pentyl
tail, so writing the anchor tail-first from the drawing would have produced the F-037 shape
in the folder that was trying to avoid it.

My own third derivation of `3A4-SC5` agrees with both: `IRFDGVDUNLJDIH-UHFFFAOYSA-N`,
C47H86N2O13S3.

Confidence: **high** that the anchor is sound and independently built. But see §3.1 — one
anchor over 489 structures is the folder's largest residual risk, and it is disclosed.

### 1.10 Scheme 1 corroborates every in-vivo core, with printed spacings

Scheme 1's "Amine Core" box prints a carbon-spacing annotation under each core that the
folder does not cite as a check but which independently confirms all seven:

`3A4 (2-3)`, `3A4V (3-2)`, `3A4V2 (2-2)`, `4A3 (3)`, `4A5 (2)`, `5A4 (2)`, `5A5 (3)`.

Every one matches `01_blocks.csv` exactly (e.g. `3A4V` = `OCCC[NH]CC[NH2]`, a 3-carbon
hydroxyalkyl and a 2-carbon diamine arm = "(3-2)"). This is a free printed cross-check and it
passes 7/7. I also confirmed the drawn `SC5`–`SC10` are linear thiols of exactly 5–10
carbons, which bears on §2.1.

Confidence: **high**.

### 1.11 Flags, metadata and hygiene

- `rows_with_a_caveat` recomputed from the emitted rows = **180**, matching `09_emit_qc.json`
  exactly. Decomposition verified: 135 `positional_ambiguity` (4 partial cores × 31 vitro =
  124, plus 11 vivo) + 68 `value_is_bound = lower` − 23 overlap = 180.
- `confidence`: 370 high / 135 medium, and **`confidence != high` ⟺ `positional_ambiguity = 1`**
  on all 505 rows. No hard-coded literal (F-011); the value is derived from whether the core
  is partially substituted.
- **0 rows with a blank `curation_note`**, and 0 `medium`-confidence rows without one — the
  metric F-024 added, clean on first contact.
- `value_is_bound = lower` on exactly the 68 cells in the open top band (`> log2`); the
  closed bands carry none. Consistent with D-5's stated reading.
- `Publication_PMID` = 40644303 on every row, indexed — no repeat of F-027.
- Descriptor columns deliberately blank so `build_dataset.py` recomputes them, with the
  reason recorded in `09_emit_qc.json`. Correct.
- No `Experiment_ID` collision: `SL_2025` and `SL_2025_2` are both unused in a 54-experiment
  corpus.
- `DECISIONS.md` carries 10 OPEN entries and the merge gate refuses the folder while it does.
  **That is the correct state and I am not recommending it be cleared.**

Confidence: **high**.

---

## 2. ERRORS FOUND

None of the first three are errors in this curation. Two are errors in `data/data_vitro.csv`
that this curation detected and correctly queued rather than silently "fixed", and that I
have now **resolved** against primary sources the curation said it could not see. The fourth
is a packet error. The fifth is a project-level hazard.

---

### 2.1 **`KZ_2016`'s 42 `*-SC3` lipids are wrong in the corpus. The 2026-10-28 patch went the wrong way, and it was inferred from a stale derived field.** — *which half is wrong: the STRUCTURE, not the name*

**Severity: HIGH — 42 corpus rows carry a wrong tail. Confidence: HIGH.**

`provenance/build_dataset.py` line 272 patches `KZ_2016 *-SC3` from isopentyl to
**neopentyl**. `provenance/README.md` §180 gives the reasoning:

> `6A3-SC2` and `6A3-SC3` share MolWt 1876.692 but are recorded with **different tail chain
> lengths, 4 and 3**. So SC3 is a more branched isomer … **SC3 = neopentyl** (C5, chain 3).

That inference is drawn from `LNP_ML/data_files/Zhou_dendrimer/individual_metadata.csv`, a
**derived descriptor file**. The primary file it was derived from,
`main_data.csv`, is row-aligned with it and holds the SMILES. I read the primary:

```
3A1-SC3  metadata MolWt 979.419, chain 3
         main_data.csv SMILES → CCCSCC(C)C(=O)... (n-PROPYL)
         recomputed MolWt 895.257, C42H74N2O12S3
```

The metadata's MolWt **does not match its own row's SMILES**. I checked all 1,502 rows:

| | |
|---|---|
| metadata MolWt agrees with its own row's SMILES | **1,460** |
| disagrees | **42** |
| of which tail = SC3 | **42 / 42** |

The stale rows are *exactly* the SC3 rows, and `note.txt` in that folder says why:
**"fixed 42 smiles (SC3 tails were being represented as SC2)"**. The SMILES were corrected to
n-propyl and the descriptor file was never regenerated. **The MolWt equality the patch relied
on is the residue of the very bug the upstream curator had already fixed.** Its companion
field `Num_carbon_in_tail = 3` does not discriminate either: n-propanethiol `CCCS` has a
longest chain of 3, and so does neopentyl.

By substructure over all 42 upstream SC3 rows: **42/42 n-propyl, 0 neopentyl, 0 isopentyl.**
Positive controls on the same file behave: SC2 42/42 isopentyl, SC5 42/42 not-isopentyl.

Four independent lines now say n-propyl:

1. **Upstream primary `main_data.csv`** — 42/42 n-propyl.
2. **The corpus's own `IL_tail1_SMILES` = `CCCS`** on all 42 rows. This field was *not*
   touched by the patch (the patch rewrites `CC(C)CCS`→`CC(C)(C)CS`; `CCCS` is neither), and
   it is correct on SC1/SC2/SC4/SC5 — so it is a reliable, untouched, independent field.
3. **This paper's Figure 1b**, which I read at 5×: SC3 is an unbranched 3-carbon zigzag with
   one fewer vertex than SC4, nothing like the quaternary carbon neopentyl would require.
4. **The `SCn` naming rule**, confirmed from Scheme 1 where SC5–SC10 are drawn as linear
   C5–C10 thiols.

Against it: one stale MolWt.

Lineage, for the repair note: `Original_LNPDB.csv` holds **isopentyl** on all 42 (the
original bug); `LNP_ML/.../main_data.csv` holds the upstream **fix to n-propyl**; our build
patched the LNPDB copy to **neopentyl**.

**Correct value: `CCCS`, n-propanethiol, on all 42 `KZ_2016 *-SC3` rows.** This folder's 15
SC3 structures are right. **This resolves `DECISIONS.md` D-4 in favour of candidate (a)**,
and the repair belongs in `build_dataset.py`/`README.md`, not in this folder.

This is also a textbook recurrence of the class CLAUDE.md §9 (2026-10-30) names: *"Five
findings this round were conclusions drawn from a secondary artifact rather than the primary
one … all five were wrong."* The primary was one directory away and row-aligned.

---

### 2.2 **`KZ_2016`'s 36 `3A2` lipids contradict their own `IL_head_SMILES`, and the paper confirms the head field.** — *which half is wrong: the `IL_SMILES`, not the head name or the head SMILES*

**Severity: HIGH — 36 corpus rows carry a wrong core. Confidence: HIGH.**

Corpus `KZ_2016 3A2-SC6`:

- `IL_head_name` = `3A2`
- `IL_head_SMILES` = `NCCN1CCNCC1` — 1-(2-aminoethyl)piperazine, **3 nitrogens**
- `IL_SMILES` = `...CCN(arm)CCC1CCN(arm)CC1` — **2 nitrogens**, C52H94N2O12S3

The ring nitrogen has become a carbon and the linker has lengthened; the molecule is a
piperidine analogue, not the piperazine the head field names.

I built a corpus-wide test with its own controls: for every `KZ_2016` row, does
`N(IL_SMILES) = N(IL_head_SMILES) + arms × N(IL_tail1_SMILES)`, with the arm count taken
from the first digit of the core name (which Figure 1a's "# of Arms" bracket defines)?

| | |
|---|---|
| cores tested | 30 |
| **cores inconsistent** | **1 — `3A2`, 36/36 rows** |
| cores consistent (controls) | 29 (`1A1`…`6A3`) |

The test fires on exactly one core and passes on 29, so it is not a loose test.

Figure 1b independently confirms the head field: I rendered it at 4× and the 3A2 drawing is a
six-ring carrying `N` on the attachment side and `NH` on the far side, with a pendant `NH2`
— three nitrogens, which is also the only way the core supports its stated 3 arms. The
formulae differ (ours C51H93N3O12S3 vs corpus C52H94N2O12S3), so this is not a convention:
one of the two is a wrong structure.

**Correct value: `NCCN1CCNCC1` as the core, i.e. this curation's 31 `3A2` structures.**
This **confirms `DECISIONS.md` D-6 in favour of candidate (a)**.

Unlike §2.1 this one is **upstream, not ours**: `Original_LNPDB.csv` and
`main_data.csv` both already carry the 2-nitrogen structure, with
`Original_LNPDB.csv`'s own head field reading `C(CN1CCNCC1)N` — i.e. the upstream file has
contradicted itself since before our build. No check in the project looks at `IL_head_SMILES`.

---

### 2.3 **219 corpus rows carry an unparseable `IL_head_SMILES`, across 7 experiments, and nothing checks it**

**Severity: MEDIUM (no `IL_SMILES` is affected; it is a latent-redundancy loss). Confidence: HIGH.**

Found while building the §2.2 test — two cores dropped out because their head SMILES would
not parse. Scanning all 12 SMILES-bearing columns across both corpus files:

| column × experiment | n | example |
|---|---|---|
| `IL_head_SMILES` × `KZ_2016` | 110 | `'CN(CCCN)CCCN_'` (`4A3-SC2`) |
| `IL_head_SMILES` × `NC_2024` | 39 | `'C(CN1CCNCC1)N_'` |
| `IL_head_SMILES` × `YR_2024` | 22 | `'C1CCN(C1)CCN_'` |
| `IL_head_SMILES` × `YY_2023` | 21 | `'CN1CCN(CCCN)CC1_'` |
| `IL_head_SMILES` × `SL_2021` | 13 | `'CCCCN(CCCC)CCCC_'` |
| `IL_head_SMILES` × `ZC_2023` | 9 | `'C1=CC=C2C(=C1)C(=CN2)CCN_'` |
| `IL_head_SMILES` × `KW_2014` | 5 | `'C[N@@]1(C)CCCNCC[N@@](C)(C)CCCNCC1'` (over-valent N) |
| **total** | **219** | |

`IL_SMILES` itself is clean (0 unparseable of 30,998) — which is exactly why this survived:
**C1/C8 test `IL_SMILES` and nothing tests the component columns.** Most are a trailing `_`,
the same shape as F-022's trailing space (error class 37), and `KW_2014`'s is a different
defect (a quaternary nitrogen written without a charge).

This matters more than cosmetics: the component columns are precisely what caught both §2.1
and §2.2. A one-line parse check over all `*_SMILES` columns, gated in `audit_baseline.json`
at the 219 currently accepted, would be cheap and would have teeth. **Not caused by this
curation** — reported as a cross-paper Wave-1 finding.

---

### 2.4 **The packet's row-count "DISCREPANCY" is the packet's error**

**Severity: LOW (it would have sent a reviewer after a non-existent problem). Confidence: HIGH.**

`packet_siegwart_40644303.md` line 5: *"Rows read from `09_rows_vitro.csv` — 465 rows **(queue
says 505 — DISCREPANCY)**"*. The queue is right: 465 + 40 = **505**, stated in
`09_emit_qc.json`, in `PROVENANCE.md`, and confirmed by my own count. The packet lists
`09_rows_vivo.csv` three lines later and still reports a discrepancy, so `brief.py` is
counting rows from one file while comparing against a total that spans both. Any folder with
a vivo split will trip this.

Two smaller packet defects:

- The DECISIONS listing emits a spurious first entry **"D-n `<one line>`"** — it has parsed
  the format template in `DECISIONS.md` as an eleventh decision. The file has 10.
- The packet's `curation_note` breakdown shows only the top 5 of the distinct strings
  (they sum to 442 of 465), presented without an ellipsis, so it reads as exhaustive.

---

### 2.5 **The 360-lipid overlap is real, but it is 100% within one laboratory and it is cross-cargo. Admitting it to Stage 5 would move Result 3 by a confound.**

**Severity: HIGH for the project (not an error in this folder — D-10 raises it; this is the
measurement D-10 asks for). Confidence: HIGH on the counts, MEDIUM-HIGH on the interpretation.**

**The overlap is genuine.** Verified by full InChIKey against `data/data_vitro.csv` and
`data/data_vivo.csv`, with working positive controls (DSPC 7,038 rows; cholesterol 30,412;
and — see §4 — an MC3 control that I had to fix after my own first SMILES returned 0).

| corpus experiment | distinct structures shared | paper | lab |
|---|---|---|---|
| `KZ_2016` | **360** | Zhou K, PNAS 2016 (PMID 26729861) | Siegwart, UTSW |
| `SL_2020` | 8 | **Lee SM**, Angew Chem 2021 | Siegwart |
| `ZC_2023` | 1 | Chen Z, JACS 2023 | Siegwart |

**Every shared structure is shared with the Siegwart laboratory.** `SL_2020` has the same
first author as this paper. `SL_2021` (Liu S, Nat Mater 2021, iPhos) and `ZT_2026`
(`curation/siegwart_42634137`) are also Siegwart; `ZT_2026` shares no structure (different
chemistry), so there is no duplicate-merge risk there. This is **F-068 in its sharpest form**:
`Experiment_ID` is a screen key, and five screen keys here are one laboratory.

**It is not independent rediscovery.** The paper says so in its own words (fulltext line 192):
*"we revisited our previously reported library of dendrimer-like ionizable amino lipids,
initially studied for siRNA"*, and line 327: *"465 of the ionizable lipids from the initial
siRNA screen¹¹"*. `KZ_2016` is 42 cores × 36 tails = 1,512; this is a **15 × 31 sub-grid of
that same grid**. Same lab, same library, same compounds — the cargo changed.

**What it would do to Result 3, measured.** Stage 5 currently has 231 shared lipids over 13
pairs and a pooled n=236 at ρ=0.083, with the largest single pair (`AA_2008`↔`KW_2014`) at
135 lipids. One merge would add a single pair of **360**, i.e. ~61% of all shared lipids in
the corpus. I measured what that pair would contribute, using my own heat-map read against
`KZ_2016`'s stored values:

| slice | Spearman | n |
|---|---|---|
| **pooled, all 465** | **0.5407** | 465 |
| SC (alkyl) tails only | 0.4135 | 225 |
| SO (hydroxy) tails only | −0.1241 | 120 |
| SN (amino) tails only | undefined — **all 120 in one band** | 120 |
| median **within-core** (hold core, vary tail) | 0.6584 | 15 cores |
| median **within-tail** (hold tail, vary core) | **0.0829** | 15 of 31 defined |

The pooled 0.54 is almost entirely the **tail-class marginal**: alkyl tails work for both
cargoes, amino and hydroxy tails work for neither (all 120 SN lipids sit in the lowest mRNA
band, carrying zero information). Hold the tail fixed and vary the head — the comparison that
is actually about transferable structure–activity learning — and it is **0.083**, which is
indistinguishable from the project's own cross-laboratory pooled figure.

So admitting this pair would roughly triple Result 3's headline on the strength of a
marginal, in a pair that is same-lab and cross-cargo. That is the in-vivo
"benchmark-only pair" hazard of CLAUDE.md §8, at 360 lipids instead of 2. **D-10 is the right
place for this and it must not be resolved by whoever merges the folder.**

A second caveat on the Q2 shift control, worth recording: the shifted-grid correlations floor
at **0.53–0.60**, not at zero, because the row/column marginals survive a shift. Q2's
"collapses to ≤0.60" is still a clear discrimination from 0.899, but the floor is marginal
structure, not chance, and the write-up should say so.

**Finally, the overlap count depends on §2.1/§2.2.** Today 360 of 465 match and 105 are
"new" — and those 105 are *exactly* the three disagreements, with no residual:

```
3A2 × 31 tails        = 31
5A2 × 31 tails        = 31   (D-8 convention)
6A1 × 31 tails        = 31   (D-8 convention)
SC3 × 12 other cores  = 12
                       ---
                        105
```

If §2.1 and §2.2 are repaired in `KZ_2016`, the overlap becomes **406**; the remaining 62 turn
on D-8, which is probably permanently unresolvable. Whoever reports "360" should report which
corpus state it belongs to.

---

## 3. UNVERIFIABLE

### 3.1 488 of 489 structures have no external anchor — and my check does not close that

**Severity: MEDIUM. Confidence: HIGH that this is the residual risk.**

My 489/489 agreement (§1.1) tests the **enumerator**, not the **building blocks**. My core
and tail templates and theirs both come from reading Figure 1b, so a misread of the figure
would be common-mode and would reproduce identically in both. What is independently anchored:

- `3A4-SC5` against Scheme 1's drawn product (C10, 1 structure);
- the 7 in-vivo cores against Scheme 1's printed `(2-3)`/`(3-2)`/… spacings (§1.10);
- `SC5`–`SC10` against Scheme 1's drawn thiols;
- `3A2` and `SC3` against the corpus (§2.1, §2.2 — these now *contradict* the corpus, in this
  folder's favour).

What is **not** anchored anywhere: the 8 `SN` tails, the 8 `SO` tails, `SC1`/`SC2`,
`SC11`–`SC18`, and the cores `3A1`, `3A3`, `3A5`, `4A1`, `4A4`, `5A1`, `5A2`/`6A1`,
`5A3`/`6A2`. I read each of them off Figure 1b and they are consistent with `01_blocks.csv`,
but that is one more pair of eyes on the same raster, not an independent source. The SI's NMR
covers only the 42 in-vivo lipids (SC5–SC10), so it cannot help.

The folder states this honestly ("C10 therefore covers 1 of them") rather than implying
coverage it does not have.

### 3.2 D-4 is resolved; D-6 is confirmed; D-8 is probably permanently open

- **D-4** — resolved by §2.1. n-propyl. Requires a corpus repair, not a folder change.
- **D-6** — confirmed by §2.2. Piperazine. Same.
- **D-8** (`5A2`, `6A1` on tetraethylenepentamine, 62 structures) — I could not settle it and
  I do not think it is settleable. A one-pot aza-Michael on a polyamine gives a statistical
  mixture; the paper states only the arm count. I verified the folder's claim that the other
  three partial cores (`5A3`, `5A4`, `5A5`) are symmetry-unique and agree with the corpus:
  in my new-structure breakdown they appear only via their `SC3` member, so they match
  `KZ_2016` on every other tail. The folder's handling — a stated convention, applied once,
  `positional_ambiguity = 1` and `confidence = medium` on all 135 affected rows, with a note
  naming the convention — is the right answer. **The flag is the deliverable, not a structure.**

### 3.3 D-2's evidence requirement is narrower than the folder states — and the answer is in the fulltext it already holds

**Severity: LOW. Confidence: HIGH on the quotation, MEDIUM on the conclusion.**

D-2 records *"EVIDENCE NEEDED: PMID 26729861's own methods … We did not fetch that paper."*
For the question D-2 actually asks — *should the siRNA half be emitted?* — the paper's own
Methods already answer it (fulltext line 355):

> "We analyzed the data **from the prior siRNA screen** and qualitatively classified the
> ionizable amino lipids as poor, intermediate, or high-performing carriers of siRNA
> (0–25%, 25–70%, and 70–100% knockdown, respectively). This categorization … was **also**
> conducted with the Luc mRNA expression **from this IGROV-1 screen** (log-transformed
> luminescence values of 0–1, 1–2, and >2, respectively) to generate the heatmap."

The siRNA stripe is explicitly a **re-binning of `KZ_2016`'s data**, and the sentence
attributes the IGROV-1 screen to the mRNA half only. Emitting it would put one lab's same
measurements under a second `Experiment_ID` at coarser resolution — the "fabricate a lab"
failure that kept `XD_2026`'s LipidDB out. **The folder's choice (a) is right and is
better-evidenced than D-2 claims.** What the fetch *would* settle is the separate question of
whether `KZ_2016`'s `Model_type = HeLa` is correct — candidate (c) — and that remains open.

This also strengthens D-5: the same sentence presents "0–1" for mRNA in exactly the parallel
construction as "0–25%" for siRNA, where the bands demonstrably partition a closed range. The
folder's reading (a), closed band, no bound flag, is the literal one and the better-supported
one.

### 3.4 D-1, D-3, D-7, D-9 — I add nothing

- **D-1** (one in-vivo experiment or two): I confirmed the ~55× gap (Fig 2b max 7.37e5 at an
  axis topping 1e6; Fig 3b max 4.09e7 at an axis topping 6e7), confirmed the panels share no
  lipid, and found nothing in the fulltext or SI linking the scales. **One thing to add:** the
  provisional split is perfectly confounded with chemistry — `SL_2025` is the 3A4/3A4V/3A4V2
  family (17 rows), `SL_2025_2` is 4A3/4A5/5A4/5A5 (23 rows). A source classifier would
  separate them at ceiling from the core alone, which makes this a within-paper instance of
  the exact confound the manuscript is about. Worth a line wherever D-1 is resolved.
- **D-3**: I confirmed `3A4V-SC5` has no bar drawn and `5A4-SC10`'s refusal is the F-018
  signature. Dropping both is defensible; inventing a bound is not.
- **D-7** (licence): outside my competence and not specific to this paper.
- **D-9** (a three-level ordinal as `Experiment_value`): a project decision. I note the
  measured consequence — 325 of 465 in one band and all 120 `SN` lipids tied — which is
  relevant evidence for it, and which bears on D-10 (§2.5).

### 3.5 What a mass could and could not settle

For the record, since the brief asks: **HRMS would settle D-4** (propyl vs neopentyl differ
by C2H4 per arm, 84 Da over 3 arms) and **would settle D-6** (C51H93N3 vs C52H94N2). It
would **not** settle D-8, the arm placement on TEPA, nor the α/β regiochemistry of §1.2 —
both are mass-invariant. The SI contains no HRMS and no mass spectrometry of any kind; its
only numbers are diameters, PDIs, pKa, E.E. and NMR shifts, which I confirmed by reading all
28 pages.

---

## 4. CORRECTIONS TO MY OWN WORK

Per the standing instruction that reviewer findings are claims:

- **My first MC3 positive control returned 0 rows.** I had written MC3 with a propanoate and
  the wrong chain, giving `XUTAXVLZUHLQHU`, which is in no corpus row. Had I not re-derived
  it from the chemical name I would have reported a clean zero from a control that could not
  match a molecule known to be present — the exact failure the brief warns about. The
  overlap counts in §2.5 were taken with DSPC (7,038) and cholesterol (30,412) controls
  passing, and re-verified.
- **My first independent enumeration scored 480/489 and the defect was mine**, not theirs
  (§1.1, SMILES ring-closure digit collision).
- My cholesterol control counts 30,412 rows against the brief's 30,397; I matched on any cell
  containing "holesterol" rather than a named column, so mine is the looser count. Not a
  discrepancy in the data.

---

## 5. SEVERITY AND CONFIDENCE SUMMARY

| # | finding | where | severity | confidence |
|---|---|---|---|---|
| 2.1 | `KZ_2016` SC3 = neopentyl should be n-propyl, 42 rows; patch inferred from a stale derived field | **corpus** (`build_dataset.py`, `README.md`) | HIGH — changes 42 structures | HIGH |
| 2.2 | `KZ_2016` `3A2` `IL_SMILES` has 2 N, head field and paper say 3, 36 rows | **corpus** (upstream LNPDB) | HIGH — changes 36 structures | HIGH |
| 2.3 | 219 unparseable `IL_head_SMILES` across 7 experiments, unchecked | **corpus** | MEDIUM — latent | HIGH |
| 2.4 | packet reports a row-count discrepancy that does not exist; 2 smaller packet defects | **packet** (`brief.py`) | LOW | HIGH |
| 2.5 | 360-lipid overlap is 100% same-lab and cross-cargo; pooled ρ 0.54 is marginal-driven (within-tail 0.083) | **project / D-10** | HIGH for Result 3 | HIGH counts, MED-HIGH interpretation |
| 3.1 | C10 covers 1 of 489; 488 blocks rest on one reading of Fig 1b | folder, disclosed | MEDIUM | HIGH |
| 3.3 | D-2's stated evidence gap is already closed by the paper's own Methods | folder | LOW | HIGH |
| 3.4 | the provisional in-vivo split is confounded with core family | folder / D-1 | LOW-MEDIUM | HIGH |

**Changes values?** No — not one `Experiment_value` in this folder is affected by anything I
found. **Changes structures?** Not in this folder; 78 rows in `KZ_2016`. **Cosmetic?** 2.4
only.

### Recommendations

1. **Do not clear `DECISIONS.md` to merge.** The gate is doing its job. D-10 in particular is
   a project decision with a measured cost (§2.5).
2. **Repair `KZ_2016` separately**, under its own finding numbers, with the baseline and
   `CHANGELOG.md` entries in the same commit — and seed the defect to confirm any new check
   fires before trusting a clean run.
3. **Add a parse check over all `*_SMILES` columns** (§2.3), gated at the 219 currently
   accepted so a decrease fails too.
4. **Fix `brief.py`'s multi-file row count and its `DECISIONS.md` template parse** (§2.4).
5. **Consider populating `IL_head_SMILES` / `IL_tail1_SMILES` / `HL_SMILES` / `CHL_SMILES` /
   `PEG_SMILES` on the emitted rows.** All five are known exactly here and all five are
   currently blank. The component columns are the only reason §2.1 and §2.2 were findable in
   `KZ_2016`; leaving them blank forgoes that redundancy for 505 rows.

---

## 6. WHAT I DID NOT GET TO

- **I did not re-measure 34 of the 40 in-vivo bars.** I re-derived the Fig 2b axis and the six
  3A4-series bar tops (6/6 exact); the 3A4V, 3A4V2 and all of Fig 3b I checked only through
  their reported calibration residuals and refusal logic, not pixel by pixel.
- **I did not independently reproduce Q5**, the Figure 5a scatter-vs-bar cross-check (median
  ratio 1.0124 over 32 lipids). It is the strongest internal check in the folder and the one
  that would catch a uniform scale error, so it deserves a reviewer's own re-derivation.
- **I did not re-derive the SI's NMR hydrogen counts** (Q8, 19 of 25 peak lists summing
  exactly). I confirmed the SI contains the peak lists and covers only SC5–SC10.
- **I did not fetch PMID 26729861.** Its reagent list would settle D-4 a fifth time and its
  Methods would settle whether `KZ_2016`'s `Model_type = HeLa` is right (D-2 candidate (c)).
  §2.1 stands without it.
- **I did not attempt the Figure 1c Venn.** The folder reports it as irreproducible under the
  caption's own rule and does not use it; I accepted that rather than re-deriving it, so that
  contradiction is unconfirmed by me.
- **I did not examine `stage08_qc.py` (34 kB) line by line**, nor `curation/lib/assemble.py`.
  My posture was to re-derive outputs independently rather than audit their code, so a defect
  that happens to be common to their enumerator and my templates would have escaped me — but
  §1.1's 489/489 is across genuinely different implementations.
- **I did not check `lit_review`'s `ADDED` registry** or whether this PMID is tracked there.
