# Wave-1 adversarial review — the four folders that "claim rows but have no mergeable file"

`linyx_41851281` · `palczewski_41654520` · `xuq_32520527` · `sd_2026`
Read-only. Reviewer: Wave-1 agent, 2026-11-10. Nothing outside this file was written.
Not run, per brief: `structure_check.py`, `check_anchor_independence.py`, `check_anchor_coverage.py`.

---

## 0. HEADLINE — the premise of the assignment is wrong for three of the four

**The packets' "NO SCHEMA-VALID ROW FILE FOUND" is a FALSE NEGATIVE in three of
four cases, and the structures are already joined onto the value rows.**

`provenance/review/brief_staged.py::emitted_file()` accepts a file only if its
header contains the literal string `IL_SMILES`. Its own docstring claims it
selects "by SCHEMA, not filename" — but one literal column name *is* a filename
by another route. Measured:

| folder | SMILES column actually present on `05_values.csv` | rows with a structure | packet verdict |
|---|---|---|---|
| `linyx_41851281` | **`smiles`** (+ `inchikey`) | **138 / 140** | "no row file" — WRONG |
| `palczewski_41654520` | **`canonical_SMILES`** (+ `candidate_SMILES`) | **59 / 64** | "no row file" — WRONG |
| `sd_2026` | none (correctly — see §4) | 0 / 18 | "no row file" — right, wrong reason |
| `xuq_32520527` | none; structures in `03_enumerated.csv` | 0 / 33, joinable 33/33 | "no row file" — right |

So the task brief's claim that for `linyx`/`palczewski` "the structures sit in
`03_enumerated.csv`, never joined" is **false**. For `linyx` I verified the join
is not merely present but *consistent*: 134 of 140 value rows carry an `ml_id`,
all 134 resolve into `03_enumerated.csv`, all 116 enumerated ids are used, and
the `smiles` string on the value row is **byte-identical to the enumerated one on
all 134 rows (0 mismatches)**. For `palczewski`, `05_values.csv`'s
`canonical_SMILES` is byte-identical to `03_enumerated.csv`'s on **64/64** rows.

**This is error class 37's shape again** (a check that canonicalises or hard-codes
its own input surrenders the ability to see non-canonical input) crossed with
F-038 (`.get(k, 0)` after a column rename). Here it is worse than F-038 because
the miss is reported as an *affirmative finding* — "that is itself a finding: this
folder cannot be merged as it stands" — so a reviewer who trusted the packet
would have written off ~200 verified, tier-1 rows.

Corpus-wide measurement of the same defect, over all 46 `05_values.csv` files in
`curation/`: the SMILES column is spelled `IL_SMILES`, `smiles`, `SMILES`,
`canonical_SMILES`, `il_smiles`, `emitted_SMILES`, `product_SMILES` or
`parent_SMILES`, and the row file is named `06_rows.csv`, `06_emitted_rows.csv`,
`05_values.csv`, `05_values_flagged.csv`, `07_<EID>_rows.csv`,
`07_emit_preview.csv`, `03_enumerated.csv`, `09_rows_vitro.csv`,
`10_LL_2026_rows.csv` or `11_QL_2026_rows.csv`. See §6 for the check that should
exist.

---

## 1. `linyx_41851281` — Lin lab, Nat Biomed Eng 2026, PMID 41851281

### VERDICT: **CLERICAL.** There is no join to make; the join is already done and verified.

### The join key
`ml_id` (1–121), unique on both sides: **116 unique in `03_enumerated.csv`, 116
unique among the 134 value rows that carry one, 0 orphans in either direction.**
Not a name-token join — `ml_id` is a numeric primary key the authors publish in
`Lipid_ML_ID_Correspondence_Table.xlsx` (MOESM7) and reprint in Supplementary
Fig. 1b, and the curation enumerated it twice from those two independent loci.
The remaining 6 value rows are 4 benchmarks (joined by name to the reagent
registry) and 2 `external_ML_test` rows with no structure.

### What I re-derived independently (not re-ran)

- **Values: 112/112 verified against the publisher Source Data.** I read
  `xl/worksheets/*.xml` of `41551_2026_1640_MOESM4_ESM.xlsx` sheet
  `Supplementary Fig. 10` directly (my own XML reader, not their script):
  113 label rows = 3 benchmarks + library ids 1–110. Every replicate and mean in
  `05_values.csv` agrees. The 13 apparent mismatches are **6-significant-figure
  storage in the CSV** (max relative deviation 4e-6); `log10_mean` is computed
  from the full-precision value and agrees to <1e-4 on all of them. The only
  label the curation renamed is the workbook's `MC3` → `DLin-MC3-DMA`, which is
  correct and avoids the brief's own documented false zero (`MC-3` matches 0
  corpus rows for a lipid present 3,098 times).
- **Structures: 116/116 self-consistent and externally anchored.** RDKit
  recomputation from the stored SMILES: 0 unparseable, 0 non-canonical, 0 formula
  disagreements, 0 exact-mass disagreements, 0 InChIKey disagreements, **116
  distinct full InChIKeys.**
- **The HRMS anchor is real, not tautological.** I grepped the SI (MOESM1 p.41 ff)
  and the strings *are* there verbatim: "Lipid 1 … Mass calculated m/z 729.6271",
  "Lipid 2 … 727.6591", "Lipid 3 … 742.6588". The authors print the **neutral**
  monoisotopic mass while labelling it m/z; 115 of 116 reproduce at ±0.0001 Da.
- **The one deviation is already adjudicated.** ml_id 47 (`H8-OT1(3)`) is
  +0.0070 Da off its printed 1119.9811. `08_contradictions.csv` C-03 calls it a
  digit transposition (9881→9811) and pins head H8 independently off lipid 48
  (0.0000 Da). I reproduced both halves. Correct.
- **Head identity cross-confirmed from a *different paper*.** `linyx`'s heads are
  local labels H1–H14. Comparing against `xuq_32520527`'s independent read of the
  same vendor amine catalogue from a 2020 Bioconjugate Chem paper:
  `H2 = CN(C)CCCN = catalogue 80`, `H4 = NCCN(C)CCN = 113`,
  `H8 = CN(C)CCCNCCCN = 400`, `H9 = NCCCN(C)CCCN = 306` — all four byte-identical
  canonical SMILES, derived from two different figures in two different papers by
  two different agents. That is a stronger structural attestation than reading a
  raster twice, and it is five independent spot-checks in one.
- **15 structures are already in `data/`, and the curation found them.** I ran my
  own full-InChIKey sweep with the brief's working controls, and the controls pass
  by asserting named molecules: `HL_name` contains `DSPC` → **7,038** rows;
  `DLin-MC3-DMA` by InChIKey `NRLNQCOGCKAESA-KWXKLSQISA-N` → **3,098** across two
  spellings (3,087 `DLin-MC3-DMA` + 11 `Dlin-MC3-DMA`); cholesterol by
  `CHL_SMILES` InChIKey → **30,397**. All three reproduce the brief's numbers
  exactly. *(My first MC3 control returned 0 because I hand-wrote the SMILES
  without `/C=C\` stereo — F-070's amendment, on me. Retracted and redone off
  `provenance/reagent_registry.py`.)*
  My sweep finds **15 of 116** already present, under `AA_2008`, `KW_2014` and
  `NC_2024`; `07_qc.json → corpus_overlap` lists **exactly the same 15**, and
  `DECISIONS.md` D-7 queues the question. **I initially logged "the overlap was
  never measured" as a finding and am retracting it** — the result exists, it is
  just in `07_qc.json` rather than a file named `02_corpus_overlap.json`, so a
  filename-keyed scan (mine, and the packet's) misses it. Same class as §0.

### ERRORS FOUND

**L-1 (substantive, affects the project's headline measurement). The proposed
`Experiment_ID` splits one lab into two sources.** `DECISIONS.md` line 6:
*"Proposed `Experiment_ID`: **YL_2026** (in vitro) / **YL_2026_2** (in vivo)"*.
The corpus convention is the opposite: in-vitro and in-vivo rows of one screen
live in two **files** under **one** `Experiment_ID`. Measured in `data/`:
**9 of 54 Experiment_IDs appear in both pools** — `BL_2023`, `BL_2024`,
`JW_2024`, `KW_2014`, `KZ_2016`, `XH_2024`, `XH_2024_2`, `ZC_2023`, `ZT_2026`.
`_2` is reserved for a second *screen* (`LL_2026_2`, `LX_2024_2`). Emitting
`YL_2026_2` for the vivo arm would add a spurious source to S, and would
manufacture exactly the source-separability signal Results 1/4 measure.
*Severity: changes no value or structure; changes S and the fold structure of
Stages 3, 5, 6 and 7.*

**L-2 (substantive-procedural). The `Experiment_ID` is derived from the
CORRESPONDING author, not the first author, and the first author is not
recoverable from the assets in hand.** `YL_2026` = Yao-Xin Lin, named in
`00_value_sources.json` as `corresponding`. `EXTRACTION_BRIEF.md` §1 is explicit:
*"`Experiment_ID` derives from the FIRST AUTHOR, and a folder named for the PI is
wrong as a key."* The folder *is* named for the PI (`linyx_`). I looked for the
first author in everything on disk: MOESM1 p.1 is a title page with no author
list; MOESM3 (Peer Review File) names only "Corresponding Author: Professor
Yao-Xin Lin"; the main article PDF **was never acquired** and
`00_value_sources.json` says so. So the ID cannot be set correctly from this
folder as it stands. `YL_2026` is free in the corpus, as is any alternative.
*Severity: metadata; but it is the join key for every cross-source analysis, and
the brief records five curations that correctly overrode a PI-named folder.*

**L-3 (minor). `Experiment_value_within_model` and `value_normalization_unit` are
present as columns and empty on all 140 rows; no z was computed.** Correct
posture for a staged folder (nothing is merged), but it means the "140 emitted
rows" in the queue are not corpus rows yet.

**L-4 (minor). Two rows are not emittable at all.** `Lipid T1`/`Lipid T2`
(`external_ML_test`) carry a measured value with **no SMILES** and an explicit
note that connectivity was not resolved. 138, not 140.

### Not an error, but it should be in front of the user

**Three of the six lipids the paper presents as AI designs are molecules already
in this corpus.** P6 = ml_id 121 = `H1-OT7(2)` = `AA_2008`'s `28_O10` *and*
`KW_2014`'s `A28_C10`; P5 = ml_id 120 = `H2-OT5(2)` = `KW_2014`'s `A80_C14`;
P4 = ml_id 119 = `H9-OT7(4)` = `KW_2014`'s `A306_C10` *and* `NC_2024`'s
`A306_O10`. And `08_contradictions.csv` C-11 already records, from the peer-review
file, that P1 "was present in the first set of 60 lipids screened … before
computational modeling was performed". Taken together: 4 of 6 "predictions" are
not novel. The curation emits only measured values and no model output, which is
right. For the project this is **Result 3 material of the best available kind** —
4 lipids re-measured by a 2026 lab against 2008/2014/2024 labs, with chemistry
held constant.

### What it would take to make it mergeable
1. A user ruling on L-1 (one ID across both pools) and L-2 (first author).
2. One z of log10 per pool (`value_normalization_unit = experiment`), 119 vitro
   and 19 vivo rows.
3. Rename `smiles` → `IL_SMILES` and emit `06_rows.csv` in the live 79-column
   header. Drop the 2 structureless rows.
4. Decide D-6: the 3 benchmarks are formulated differently from the library
   (DSPC/ALC-0159 vs DOPE, different molar ratios). Flag or drop, don't pool
   silently.
5. Record the licence. `00_acquisition_raw.json` holds no licence string for this
   paper; nothing on disk establishes §8 status. **This is an acquisition gap, not
   a decision.**

### CONFIDENCE
**High** on the values (tier-1 workbook, re-read from the XML), **high** on the
structures (116/116 recomputed + 115/116 printed-mass anchors + 15 full-InChIKey
identities with three independently curated corpus experiments + 4 heads
corroborated by a second paper), **high** on L-1 (measured: 9 corpus IDs span
both pools), **high** on L-2. **Not verified:** the regiochemistry of the 3- and
4-arm H4/H9 members (one drawing, Supplementary Fig. 12, is the only external
evidence; mass is identical either way — the Lin-2026 class), and the alkene
positions in T1/T2 (the curation says so itself).

---

## 2. `palczewski_41654520` — Zhang et al., Nat Commun 2026, PMID 41654520, CC BY 4.0

### VERDICT: **NOT a join problem. ADMISSIBILITY-BLOCKED, and the block is real.**
The join is complete (64/64). The question queued as D-4 is whether this paper
belongs in the corpus at all, and on the brief's own §3 line I think the honest
answer is "not without a user ruling, and probably not".

### The join key
`compound` (`CBB1`…`CBB22`, `CBBZ0`…`CBBZ10`, `CBBox11`, `CBB18`, `CBB19`).
`canonical_SMILES` is **already on the value row** and byte-identical to
`03_enumerated.csv` on **64/64**. `CBB18`/`CBB19` resolve to *two* enumerated rows
each (`_KLVPM` / `_KLPVM`) and correctly carry a blank `canonical_SMILES` with
both options in `candidate_SMILES`. So the key is unique where it must be and
deliberately non-unique where the paper is ambiguous. **No join hazard at all.**

### What I re-derived independently

- **Values: 64/64 verified against the Source Data workbook.** Reading
  `41467_2026_69077_MOESM6_ESM.xlsx` sheets `Figure 2` and `Figure 4` from the XML:
  Fig 2c = 22 library + `LF`/`C`/`NE`; Fig 2d = 12 + 3; Fig 4c = 22 + 3;
  Fig 4d = 8 + 3. 22+12+22+8 = **64**, exactly the emitted count, and the three
  controls are excluded in every panel. Every `value_mean`, `value_sd`,
  `n_replicates` and `;`-joined replicate list reproduces; the 39 apparent
  mismatches are **4-decimal rounding of `value_mean`** (and sample-vs-population
  sd on three rows at the 4th decimal). No value is wrong.
- **Structures: HRMS anchors verified verbatim in the SI, 31/33 exact.** I grepped
  `41467_2026_69077_MOESM1_ESM.pdf` pp. 18–50 and found the printed strings:
  "HRMS calculated for C117H174N9O9S2 3+, 637.7620; found 637.7618 (error 0.3 ppm)",
  "C119H178N9O9S2 3+, 647.1057", "C95H150N9O7S2 3+, 531.3706",
  "CBB12 … C85H130N9O7S2+, 1450.9373", "CBB17 … C65H82F6N9O7S3+, 1278.5678",
  "CBB3 … C105H144N9O7S2+, 1707.0625". The three non-exact cases are all
  **the paper contradicting itself**, each flagged: CBB7's printed m/z is 0.3345
  high = exactly 1.0035 Da over charge 3 (D-3, with an independent homolog-spacing
  argument); CBB12's printed formula has H130 where the structure needs H128;
  CBB17's printed formula has S3 where the structure needs S2. `check_C` (our mass
  vs our formula) is 0.0000/−0.0001 on all 33. No tolerance was widened.
- **Chemotype, measured.** CBB1 C117H171N9O9S2, 1910.26 Da, 2 tertiary amines,
  2 sulfonamides, 14 rings; CBB2 the quaternised analogue at charge +2; CBB22 =
  Coomassie Brilliant Blue G-250 free acid, 831.30 Da, **0 tertiary amines**.
  SMARTS controls assert their premise: the tertiary-amine SMARTS returns 1 on
  N,N-dimethylbutylamine, the quaternary SMARTS returns 1 on tetramethylammonium.
  These are triarylmethane-dye bis-sulfonamide amphiphiles at 830–2,000 Da.

### ERRORS FOUND

**P-1 (substantive, and it decides the folder). There is no nucleic-acid cargo
anywhere in this paper.** Measured off the emitted rows: cargo is
`Cre recombinase (protein)` on 34 rows and `ABE8e ribonucleoprotein` on 30.
Readouts are `RFP conversion (%)` (Cre→loxP colour switch) and `eGFP/mCherry (%)`
(base-editing colour switch). Both readouts are individually admissible under
`EXTRACTION_BRIEF.md` §3 — Cre→reporter recombination and editing efficiency are
both listed — but the *delivered species* is a protein or an RNP, not mRNA/siRNA,
and the particle is a dye-amphiphile complex formed in OptiMEM, not an LNP.
`DECISIONS.md` D-4 states this plainly and refuses to decide. I agree with the
refusal and would add one measurement to it: **the corpus-emittable quantity here
is 30 distinct structures carrying a value** (not the queue's 34; 35 distinct
InChIKeys are enumerated, 30 of them have a measurement), and **none of those 30
is an ionizable lipid in the corpus's sense** — 5 of the 6 I inspected carry 0–2
tertiary amines on a sulfonated dye core at 1.4–2.0 kDa.
*Severity: decides whether 64 rows and 30 structures enter the corpus. Changes no
value and no structure.*

**P-2 (substantive-procedural). The proposed `Experiment_ID` is the
corresponding author's, and the first author is on disk.** D-4 proposes
`KP_2026` (Palczewski, corresponding). I read the JATS
`<contrib-group>` in the folder's own `PMC12992809.jats.xml`: the first author is
**Jianye Zhang** (equal-contrib, aff 1). `DECISIONS.md`'s own header already says
"Zhang et al." So the correct key is `JZ_2026` (free in the corpus), not
`KP_2026`. Same class as L-2, but here the evidence to fix it is already in the
folder. *Severity: metadata; join key.*

**P-3 (minor, but it is the §8 point). The CC BY 4.0 licence buys nothing unless
P-1 is ruled in.** The task notes this is one of only four §8-clean sources in
the run. That is true and I verified the licence string is taken off the PMC JATS
`<license>` block, per brief §11 — not the publisher landing page. But the rows
are only valuable if protein/RNP delivery is admitted, and the four blocks must
then stay four blocks: `08_contradictions.csv` row 9 measures the **same compound
differing 0.55× and 1.86×** between the flow-cytometry and microscopy panels of
one cargo. Pooling Fig2c with Fig2d would be error class 22 in the other
direction.

**P-4 (minor). 20 of 64 rows are a different measurement modality.** Fig2d (12)
and Fig4d (8) are `fluorescence microscopy only (Keyence BZ-X810), NOT flow
cytometry`. Flagged on the row, which is right, but it means the "22-compound
full-library panel at n=4" that D-4 calls the strongest case is 44 of the 64 rows,
not 64.

**P-5 (minor). No `corpus_admissible` column.** Brief §3: *"Put
`corpus_admissible` yes/no with a reason on every row of `05_values.csv`, counts
summing."* Absent. (So is `06_rows.csv`, correctly — 0 admissible rows pending
D-4 means a header-only file, which the brief asks for and the folder does not
have.) See §6: **28 of 46 `05_values.csv` files in `curation/` lack this column.**

### What it would take to make it mergeable
1. **A user ruling on D-4** — this is the whole blocker. If "no", the folder is a
   well-argued zero and should be recorded as one; the structural and value work
   is sound and worth keeping as provenance.
2. If "yes": fix the ID to `JZ_2026`, emit **four** `Experiment_ID`s (one per
   assay block, never pooled across instrument), one z of log10 per block, 59
   structure-bearing rows, `Model_type` carrying cargo + instrument.
3. Settle D-1 (`KLVPM` vs `KLPVM` — mass-blind, 2 compounds), D-2 (`CBBox11` is a
   two-isomer mixture — correctly emitted with no SMILES), D-5 (`CBB22` free acid
   vs mono-sodium; the SI's own "854 mg, 1.0 mmol" weighing says sodium salt to
   0.005%).

### CONFIDENCE
**High** on the values (64/64 from tier 1, controls excluded correctly),
**high** on the structures (31/33 printed-HRMS exact, 3 self-contradictions each
pinned by a second argument), **high** on P-1 and P-2 (both read off primary
objects — the emitted cargo strings and the JATS author list). **Not verified:**
I did not re-derive the 35 large SMILES atom-by-atom against the SI's IUPAC
names; formula + exact mass agreement is blind to connectivity on molecules this
size, and on a 117-carbon bis-sulfonamide that is a real residual risk. The
folder has **no C10 anchor with a structure** — `reference_structures.csv` is
1,403 bytes for 36 compounds — so against these the guards are vacuous, not
passing (brief §9).

---

## 3. `xuq_32520527` — Li Y et al., Bioconjug Chem 2020, PMID 32520527

### VERDICT: **CLERICAL.** The join exists, is bijective, and is not a name-token join.

### The join key
`lipidoid` (`75-O16CBA` … `400-O16CBA`). **11 unique in `05_values.csv`, 11 rows
in `03_enumerated.csv`, 11↔11, 0 orphans.** `03_enumerated.csv` already uses the
column name `IL_SMILES`. The 33 value rows are 11 lipidoids × 3 formulations
(bare `4/0/0`, `-F1` `4/1/1`, `-F2` `4/1/4` lipidoid:chol:DOPE by weight), so the
key is 3:1 on the value side by design, with `formulation_suffix` carrying the
difference. This is the *right* multiplicity, not an ambiguity.

**It is also the one folder where the brief's token trap was a live risk and was
avoided.** The heads are the shared Xu/Tufts amine catalogue — `75`, `76`, `77`,
`78`, `80`, `81`, `82`, `93`, `113`, `306`, `400` — and the brief's measured
damage figures apply directly (`80` → 94 rows / 85 structures; bare `1` → 20,837).
`04_corpus_overlap.json` records that **163 corpus `IL_name`s begin with one of
these eleven numbers** and joins on full InChIKey instead. I confirmed the result:
**0 of 11 match any corpus structure** by full InChIKey, with my own controls
asserting named molecules (DSPC 7,038 / MC3 3,098 across two spellings /
cholesterol 30,397 — all three reproducing the brief). The real overlap is at the
head group only (4 heads shared with `QL_2026`, 160 corpus rows built on them,
4/4 byte-identical canonical SMILES from two different evidence routes), and the
folder says so and does not act on it.

### What I re-derived independently

- **Structures: 11/11.** Formula, exact mass, canonical form and InChIKey all
  recomputed: 0 discrepancies, 11 distinct InChIKeys. Our [M+H]+ agrees with the
  **m/z the authors print in Figure 2F** on 11/11, worst |Δ| **0.0048 Da**.
  Michael mass balance closes independently: head exact mass + n_arms × 422.2305
  (the O16CBA monomer) reproduces every neutral mass to **+0.0001 Da**, which
  pins both the arm count and the tail inventory without any pixel.
  `77-O16CBA` and `82-O16CBA` are constitutional isomers at one formula and one
  printed m/z, and the two stored SMILES differ — the folder tests this (P10e).
- **Digitisation: independently confirmed within ~1 pp, and the discrepancy is
  mine.** I measured the drawn baseline rule in `nihms-1799810-f0005.jpg` myself:
  dark horizontal runs at rows 348–351, centre **349.5 px**. The curation's
  four-tick fit (90/60/30/0, max residual 0.113 px, worst LOO 0.334 px against a
  1.5 px gate) predicts zero at **350.004 px** — 0.5 px from a rule that took no
  part in the fit. Reading four bars with my own crude solidity walk and a
  calibration anchored only on that baseline and the paper's *printed* "29%" for
  113-O16CBA, I get 93.55 / 2.13 / 67.99 / 83.40 against stored 92.41 / 1.60 /
  67.35 / 82.32 — high by 0.5–1.1 pp throughout, consistent with my including the
  ~2 px bar outline and using a 0.7%-coarser scale. The sign is uniform, which is
  the signature of my own systematic offset, not of their error. Their value
  stands.
- **The paper's own prose validates 7 of 33 and I spot-checked two:**
  "113-O16CBA and 306-O16CBA induced 29% and 44% GFP+" vs stored 28.94 / 43.59.
- **The z is on-convention.** `Experiment_value_z_of_log10` over all 33 rows:
  population sd (ddof=0), mean −6.1e-8, sd 1.00000005, max residual against my
  recomputation **4.9e-7**. One z per experiment, pooling the three formulations
  — correct per brief §11 and *not* error class 22.
- **30 checks in `07_qc.json`, 29 gating, 1 non-gating failure (P6b: 3/5 inside a
  printed band), plus the unprinted-round-number class done properly** — four
  unlabelled minor ticks predicted by a fit that excluded them, worst deviation
  0.091 pp.

### ERRORS FOUND

**X-1 (minor but real: two names for one experiment inside one folder).**
`PROVENANCE.md` line 1 reads "PROVENANCE — XQ_2020 / proposed `YL_2020`";
`07_qc.json`'s `_what` reads "stage 07 QC for **XQ_2020**". Both are free in the
corpus. `YL_2020` is the correct one by the brief's first-author rule (first
author Li Y, folder named for the PI Xu Q) — so the folder got the *rule* right
and then left the PI-derived string in two artifacts. A merge driven off either
string would work; a merge driven off both would create two sources from one
screen. *Severity: metadata; but it is the §0 failure mode applied to an
`Experiment_ID`.*

**X-2 (minor). `value_pct_GFP_positive` and `log10_value` are mutually
inconsistent at the 3rd–4th decimal.** The percent is stored rounded to 2 dp
while `log10_value` is computed from the unrounded digitisation: e.g.
`93-O16CBA-F2` stores 1.60% and 0.203048, but log10(1.60) = 0.204120 (the
unrounded read is 1.5960). Harmless for the corpus, which uses the log column —
but a downstream reader recomputing log10 from the percent will not reproduce the
z, and nothing in the folder says which column is primary.

**X-3 (minor, a stated-precision gap). No sensitivity limit is quoted in the
units the corpus uses.** `05_values_qc.json` gives `px_per_percent = 3.0721`,
i.e. **1 px = 0.3256 pp**. On the log scale that is **0.0015 decades at 92.4%
and 0.0836 decades at 1.60%** — a 56× spread across the panel. The two smallest
values (`93-O16CBA-F1` 2.25%, `93-O16CBA-F2` 1.60%) sit 2–3 px above the baseline
and are also the two most extreme z-scores (−2.24, −2.53), so they carry the most
leverage on the experiment's scale and have the worst precision in the set. Brief
§10 asks for the measured sensitivity limit; the folder reports px-space residuals
and never converts. *Not a wrong value — an unstated uncertainty on the two rows
where it is largest.*

**X-4 (minor). The "150 dpi" is the file's own JFIF density field, and the folder
does not say so.** I read it: `jfif_density (150, 150)`, `jfif_unit 1`. So it is
*derivable* (brief §7's bar is met) rather than guessed — but a JFIF header is a
claim by whoever exported the JPEG, not a measured `px / rect_pt × 72`, and
`00_value_sources.json` presents it as a plain fact.

**X-5 (minor, queue not folder). `run_queue.json` lists `blockers: []` for this
folder**, while `DECISIONS.md` carries **6 OPEN decisions** and the licence is
explicitly not CC. The queue's own `blockers` field is the richest diagnosis in
the repo for the other three folders (§5); here it is empty and misleading.

### Admissibility
**Admissible.** `%GFP-positive` by flow cytometry, GFP mRNA, HeLa, 24 h — brief
§3 lists `%GFP+` explicitly. Not uptake, not viability, not pre-translation
cargo. Delivery.

### Licence
**FAILS §8**, and the folder says so with the evidence: `00_acquisition.json`
quotes the JATS permissions verbatim — *"available for text mining … consistent
with the principles of fair use"* — and names it as the NIHMS author manuscript,
not CC. This is one of the brief's three "NIH author manuscript, not CC" cases.

### What it would take to make it mergeable
1. Pick one `Experiment_ID` (`YL_2020` by the rule) and purge the other string.
2. Left-join `03_enumerated.csv` on `lipidoid` → 33 rows × `IL_SMILES`;
   rename `Experiment_value_z_of_log10` → `Experiment_value`; fill
   `Experiment_value_within_model` (identical, one `Model_type`); emit
   `06_rows.csv` on the live 79-column header.
3. Carry the three formulations in `Model_target`/formulation columns, **not** in
   separate standardisation.
4. Add X-3's per-row precision to `curation_note` on the two sub-3% rows, or at
   minimum to `07_qc.json`.
5. Resolve the 6 OPEN decisions and the §8 licence question.

### CONFIDENCE
**High** on the join (bijective, verified both directions), **high** on the
structures (11/11 recomputed + 11/11 against the paper's printed m/z + exact
Michael mass balance + an independent reagent-weighing ladder at worst 1.63%),
**medium-high** on the values (tier 5 raster, but the baseline lands 0.5 px from
an unseen prediction, 7 printed prose values reproduce within 0.65 pp, and my own
independent read agrees within 1.1 pp with a uniform sign explained by my method).
**Not verified:** the x-axis group labels, which are read off a 1050 px raster by
eye — the folder says this outright and constrains it with 11 prose sign
statements, but no check reads the glyphs; and the benzylidene acetal's cis/trans
relationship, which the paper never draws (stored flat, and the folder records
that `assemble.stereocentres()` returns `{}` for this tail while
`FindPotentialStereo` finds two, i.e. **the library's stereo guard is blind here
and its passing is not evidence** — an honest statement I could not fault).

---

## 4. `sd_2026` — Dong et al., Nat Biotechnol 2026, PMID 42806111

### VERDICT: **SUBSTANTIVE, and the 384/18 gap is a FINDING, not a defect.**
The 384 cannot be joined to the 18 by any key, because the paper does not publish
the mapping. The curation diagnosed this precisely and I reproduced the diagnosis.

### There is no join key, and here is why — re-derived

The 18 values are for **12 named items** across 3 experiments: `ALC-0315`,
`LP-01` (benchmarks) and `SC-1…SC-5`, `LC-1…LC-5`. The 384 enumerated rows are
keyed `H{1..12}-L{1..4}-Ta{1..4}-Tb{1..2}` — a complete 12×4×4×2 factorial,
which I confirmed closes exactly (384 rows, 384 distinct InChIKeys, 0
unparseable, 0 formula/mass/canonicalisation disagreements, 32 members per
headgroup and 96 per linker). **No column in either file names the other's key.**

The published potencies for the 384 exist only in Fig. 1c/1h (circos) and
Fig. 1d/1i (Sankey), and both are unreadable for two independent reasons the
folder states and I accept: (a) the heat-map legend is a continuous ramp
annotated `Max`/`Min` with no tick, no number and no unit, and the Sankey's
"RLU (log2)" axis is a bare arrow with zero ticks — one dashed ALC-0315 line
cannot fit a scale and `axis.calibrate` would be right to refuse; (b) the rings
are a **clustering's leaf order**, which the caption says and which is published
nowhere, so assigning columns would silently permute 384 values over 384
structures. **This is the brief §2 shape** — a design supporting 384 per-LNP
values that attributes 1.

For the 10 **named** leads, which do have values, the obstruction is different
and I measured it from the primary object. Grepping the SI for the printed
`(formula; MW:)` entries and counting how many of the 384 realise each formula:

| lead | printed formula | enumerated members with that formula |
|---|---|---|
| SC-Lipid-1 | C48H93N3O2 | 3 |
| SC-Lipid-2 | C50H98N4O4 | 4 |
| SC-Lipid-3 | C50H97N3O5 | 3 |
| SC-Lipid-4 | C48H93N3O5 | 2 |
| SC-Lipid-5 | C48H95N3O4 | 4 |
| LC-Lipid-1 | C50H97N3O4 | 11 |
| LC-Lipid-2 | C50H99N3O4 | 7 |
| LC-Lipid-3 | C54H105N3O4 | 7 |
| LC-Lipid-4 | C52H101N3O4 | 12 |
| LC-Lipid-5 | C54H105N3O4 | 7 |

**My counts are identical to `07_qc.json → P7.members_per_formula` on all ten.**
Mass cannot choose, and the folder proves it rather than asserting it: swapping
LC-1's headgroup H4→H3 leaves the molecular formula **identical** (C50H97N3O4)
while the InChIKey skeleton moves `GFWACCDGFPQKTE` → `PIFFKHUIMUSBJX`. The one
member that *is* pinned is **LC-1 = H4-L4-Ta1-Tb1**, by Supplementary Table 5's
four CAS numbers, and I verified the two numbers that pin corroborate it:
SI p.20 prints "LC-Lipid-1(C50H97N3O4; MW: 804.34)" and the enumerated
`H4-L4-Ta1-Tb1` is C50H97N3O4 / 804.34 / 803.7479 exactly. The NMR closes too:
the printed integrals sum **2+1+4+4+2+3+4+7+52+9+9 = 97 = the formula's H count
exactly**, with `1.02 (s, 9H)` the pivaldehyde *tert*-butyl (L4) and
`0.87 (t, 9H)` three terminal methyls.

**So the emittable content of this folder is 18 value rows, of which at most 6
can carry a structure:** `LC-1` × 3 (pinned) and `ALC-0315` × 3 (reagent
registry). `LP-01` × 3 has no structure anywhere in the folder, and
`SC-1…5`/`LC-2…5` × 9 have none either — `reference_structures.csv` correctly
holds nine rows with an empty `reference_SMILES` so the absence is visible.
**Net new ionizable-lipid contribution: 1 structure (LC-1).** The folder knows
this: `00_value_sources.json` says "exactly ONE ionizable lipid the corpus does
not already hold (LC-1; LP-01 and ALC-0315 are benchmarks)".

### ERRORS FOUND

**S-1 (substantive). D-6 asserts the corpus already holds `LP-01`, and it does
not.** D-6: *"two of which are the benchmarks ALC-0315 and LP-01 **that the
corpus already holds from other sources**."* I searched `data/` with controls
that assert their premise (`HL_name` contains DSPC → 7,038; `IL_name` ==
`DLin-MC3-DMA` → 3,087): **0 rows carry an `IL_name` matching `LP-?0?1`, and 0
`IL_name`s begin with `LP` at all.** `ALC-0315` by InChIKey is present 313 times
(`AC_2025` 216 vitro, `LZ_2025` 96 vitro, `XD_2026` 1 **whole_body vivo**), so
that half of the claim is right. The `LP-01` half is a claim about corpus
contents made with no corpus query — **this folder has no `02_corpus_overlap.json`
and no overlap section in `07_qc.json`**, which brief §11 lists as deliverable 2.
It cannot be settled by name alone and the folder stores no SMILES for `LP-01`,
so the only honest statement available is "not determinable from this folder".
*Severity: changes no value; it is a false premise under a merge decision, and it
is the round-1/round-2 pattern of concluding from a secondary artifact.*

**S-2 (substantive, and it is the good news this folder buries). `ALC-0315`
whole-body in vivo would be the corpus's THIRD cross-experiment in-vivo
structure.** The brief states the in-vivo pool has exactly two
(`DLin-MC3-DMA` 6 sources, `SM-102` 4). I measured `ALC-0315` in `data_vivo.csv`:
**1 row, `XD_2026`, `Model_target = whole_body`.** `SD_2026`'s Fig. 1n row is
`ALC-0315`, IV, **whole body** — so emitting it takes ALC-0315 from 1 in-vivo
source to 2 and, unlike SM-102's four ear/whole-body sources, it is
**organ-matched** with the existing one. The brief says such a lipid is "rare and
load-bearing; report it prominently". Neither `DECISIONS.md` nor `07_qc.json`
reports it; D-6 instead lumps ALC-0315 in as a reason for caution. Both readings
are true and the folder gives only one.

**S-3 (new, not in `08_contradictions.csv`). `LC-Lipid-3` and `LC-Lipid-5` are
printed with the SAME formula and the SAME MW** — C54H105N3O4, 860.45 (SI pp. 26
and 32) — so they are claimed as distinct constitutional isomers that mass can
never separate, and they compete for **the same pool of 7** candidate members.
That is a *stronger* joint constraint than two independent 7-way choices (any
future mapping must place two names inside one 7-member pool, and the pool must
therefore contain ≥2 valid members — a real consistency test), and it is also a
second reason to distrust the LC-5 entry: C-4 already records that LC-5 is the
only one of the ten whose printed ¹H integrals do not close (2+2+10+21+54+9+9 =
**107** against a formula H count of **105**), while LC-3's close exactly
(2+1+13+17+54+9+9 = **105**). `08_contradictions.csv` C-4 records the integration
but not the formula collision. *Severity: no value changes; it tightens D-2 and
adds a usable test.*

**S-4 (minor, retracted finding, recorded because the analogy was persuasive).** I
first measured "the SI prints two different formulas for LC-Lipid-1 (C50H97N3O4
on p.20, C51H99N3O5 on p.58)" and had it half-written. It is **my** regex:
`((?:SC|LC)-Lipid-\d)` matched the suffix of `N-LC-Lipid-1`, which is the *second*
library (D-4, Supplementary Fig. 68, conventional stepwise chemistry). Re-run with
a `(?<![A-Za-z-])` boundary, the two libraries separate cleanly and there is no
conflict. The lesson is the brief's: a pattern that can match a superstring of its
target is a control that does not assert its premise.

**S-5 (minor). `Experiment_value_within_model` is an empty column in BOTH
`03_enumerated.csv` and `05_values.csv`.** Its presence in `03_enumerated.csv` —
a structure file with no values — is what made the packet generator rank that
file as row-like; it is dead schema and should not be there.

**S-6 (the open decision I'd put in front of the user first). Two of the three
experiments are n = 3 with 2 benchmarks.** `SD_2026_2` (intrathecal, brain) and
`SD_2026_3` (intratracheal, lung) each hold `ALC-0315`, `LP-01`, `LC-1`. A z over
n = 3, population sd, is determined almost entirely by the choice of three, and
two of the three are the field's universal positive controls. I verified the z
arithmetic — population sd, max residual 5e-5, mean 0 / sd 1 on all three blocks
— so the computation is right and the *quantity* is the problem. `CLAUDE.md` §8's
open decision on benchmark-only pairs lands directly on this, and D-6 says so.

### Admissibility
**Admissible.** Whole-animal bioluminescence radiance of an ABE-NanoLuc mRNA
reporter at 6 h — a translated readout, downstream of delivery. Not
biodistribution, not pre-translation cargo (F-071). The SAXS phase fractions in
Supplementary Table 4 are correctly *not* emitted as values and used only as a
printed-statistic check.

### Licence
**FAILS §8**, recorded with the method (D-5): no `creativecommons.org/licenses/`
link on the landing page; "Springer Nature … holds exclusive rights".

### What it would take to make it mergeable
1. **Withdraw or substantiate S-1** (run a structural `LP-01` query, or state
   "not determinable" — it has no SMILES in the folder).
2. Decide S-6 (emit the two n=3 blocks, drop them, or exclude them from Stage 5's
   pooled concordance) and D-7 (`whole_body` vs `liver` for the 12 Fig. 1n rows).
3. Emit the 6 structure-bearing rows (`LC-1` × 3 + `ALC-0315` × 3 from the
   registry) and record `LP-01` × 3 and the 9 unmapped leads as values with no
   structure — or not at all. **The corpus has no structure-less row type**, which
   is D-1's candidate (a) and is a schema question, not a curation one.
4. Record S-2 prominently: this folder's single most valuable contribution to
   Result 3 is an organ-matched second in-vivo source for `ALC-0315`.
5. The 384 stay in the folder as provenance. **They are not mergeable and should
   not be made so.** D-1's candidate (d) — write to the corresponding author for
   the screening table and the circos leaf order — is the only route that turns
   384 structures into 384 values, and it is worth the email.

### CONFIDENCE
**High** that the 384/18 gap is substantive and correctly diagnosed (I reproduced
the formula-multiplicity table exactly, the factorial closure, the LC-1 pin, and
the NMR H-count closure). **High** on the values (vector read: tick positions are
stroked PDF segments and tick values are text spans, 5 ticks at residual 0.0000–
0.0001 px, Fig. 1p's 10⁵ tick held out and predicted at 0.0 error, 55/60
replicate markers straddling every drawn mean). **High** on S-1 by name search,
**medium** on S-1 overall — a name search cannot prove absence of a structure, and
I could not test `LP-01` structurally because no `LP-01` SMILES exists on disk and
I will not write a proprietary benchmark from memory (the brief records an agent
whose SM-102, MC3 *and* ALC-0315 were all three wrong). **Not verified:** the nine
unmapped leads' structures (the Fig. 1g/1l structure strips are ~270 px per
structure — below the ±1 chain-counting floor, as D-2 says), and Supplementary
Fig. 53's four extra organs, which I did not open.

---

## 5. Is this failure shape likely in the other 23 folders? — YES, and it already is

I measured it rather than guessed. Over all 46 `05_values.csv` files and every
`*.csv` in every `curation/` folder:

**(a) The row file has no agreed name.** It is `06_rows.csv` (14 folders),
`06_emitted_rows.csv`, `05_values.csv` (7 folders), `05_values_flagged.csv`,
`07_<EID>_rows.csv`, `07_emit_preview.csv`, `03_enumerated.csv`,
`09_rows_vitro.csv`/`09_rows_vivo.csv`, `10_LL_2026_rows.csv`,
`12_LL_2026_2_rows.csv`, `11_QL_2026_rows.csv`, `01_WT_2026_rows.csv`,
`05_WT_2026_merge.csv`, `08_NG_2025_merge.csv`. `brief_staged.py`'s own docstring
records that "a scan keyed on one filename found 14 of 27 and reported a lipid
count wrong by a factor of fourteen" — and its replacement then keyed on one
*column* name and got three false negatives out of four in my sample.

**(b) The SMILES column has no agreed name either:** `IL_SMILES`, `smiles`,
`SMILES`, `canonical_SMILES`, `il_smiles`, `emitted_SMILES`, `product_SMILES`,
`parent_SMILES`, `methylated_cation_SMILES`, `block_SMILES`, `head_SMILES`.
This is the actual cause of §0.

**(c) Two folders have MORE THAN ONE schema-valid file and nothing says which the
merge would take.** `mitchellmj_38982196` has `IL_SMILES` + an
`Experiment_value*` column in **both** `03_enumerated.csv` and `05_values.csv`;
`lin_aromatic_2026` has two (`10_…` and `12_…`, legitimately two screens);
`wht_2026_falcon` has two. `brief_staged.py` resolves ties by `max(n_rows)`,
which is a silent guess. A merge gate that picks a file by row count is one
rename away from merging an intermediate.

**(d) 28 of 46 `05_values.csv` files lack the `corpus_admissible` column** that
`EXTRACTION_BRIEF.md` §3 requires on every row with counts summing — including all
four folders reviewed here. That column is the only machine-readable record of the
F-071 / delivery-only judgement, so for 28 folders the admissibility decision
exists only in prose.

**(e) 10 folders carry a `06_rows.csv` with ZERO rows** (`dahlmanje_30216729`,
`30275336`, `30394729`, `33758781`, `maohq_38583366`, `mitchellmj_41603124`,
`songx_39664425`, `xuq_23651337`, `xuq_24574196`, `zhangn_42569578`). That is the
brief's correct posture for a well-argued zero (header-only so the zero is
visible) and it is working. The folders that produced *rows* are the ones that
skipped the file.

### The check that should exist and does not

**`C19` (next free class — `ERROR_CLASSES.md` currently runs to 48, so take 49+
and confirm unused per the 2026-11-09 renumbering rule): a per-folder SCHEMA
DECLARATION, and a gate that reads it.**

Concretely:

1. **Each curation folder emits `00_schema.json`** naming (i) the one file the
   merge is to take, (ii) which of its columns is `IL_SMILES`, (iii) which is
   `Experiment_value`, (iv) its row count, (v) its `Experiment_ID`(s). One file,
   five fields, written by the curation that knows the answer.
2. **`audit.py` gains a check** that for every folder either (a) `00_schema.json`
   exists and its named file/columns/row-count resolve exactly, or (b) the folder
   declares zero rows. Discovery by column-name sniffing is then a *fallback that
   warns*, never the primary route.
3. **The gate must be made to fail on purpose** before it is trusted — brief §10,
   and the rule that cost three incidents in round 2. Seed it four ways: rename
   `IL_SMILES` → `smiles` (must still pass, via the declaration); point the
   declaration at a file that does not exist (must fail); change the declared row
   count by one (must fail); put two schema-valid files in a folder with no
   declaration (must fail, not silently pick the larger). **Credit a seed only
   when the agreement count moves**, not when something merely raises.
4. **A second, cheaper check with teeth today:** assert that
   `run_queue.json`'s `rows_emitted` equals the row count of the declared file for
   every `status: extracted` folder. Today the queue's numbers (140 / 64 / 33 /
   18) are right for all four of these folders while the packet says the files do
   not exist — so the two sources of truth already disagree and nothing notices.

**The deeper lesson, in the form this project uses:** the packet generator was
written to replace a filename-keyed scan that had failed, and it replaced it with
a column-name-keyed scan that failed the same way. **A discovery rule derived from
a convention nobody enforces is not a check — it is the convention's second
victim.** The fix is not a better sniffer; it is making each folder *declare* what
it produced, so the gate compares a claim to a file instead of guessing at both.

One more, narrower: **`brief_staged.py`'s "NO SCHEMA-VALID ROW FILE FOUND … That
is itself a finding" should not state a conclusion.** It should say "no file
matched the discovery rule; the rule is `IL_SMILES` + `Experiment_value*`; here
are the candidate headers." The current wording converted a tool limitation into
an assertion about the data, and three of my four assignments were that assertion.

---

## 6. Summary table

| folder | verdict | join key | key unique? | rows claimed | rows with a structure | admissible? | licence |
|---|---|---|---|---|---|---|---|
| `linyx_41851281` | **clerical** | `ml_id` | 116↔116, 0 orphans | 140 | **138** | yes (luminescence, vitro+vivo) | **not recorded anywhere** |
| `palczewski_41654520` | **admissibility-blocked** | `compound` (already joined) | 64/64 consistent | 64 | 59 (30 distinct structures) | **OPEN** — protein/RNP cargo, no mRNA | CC BY 4.0 ✓ |
| `xuq_32520527` | **clerical** | `lipidoid` | 11↔11, bijective | 33 | 33 joinable | yes (%GFP+) | NIHMS fair-use ✗ |
| `sd_2026` | **substantive** | **none exists** | — | 18 | **6** (LC-1 ×3, ALC-0315 ×3) | yes (translated reporter) | All-rights-reserved ✗ |

255 claimed → **236 structure-bearing and joinable**, of which **177** are also
delivery-admissible today without a user ruling (138 linyx + 33 xuq + 6 sd_2026).
Two of the four need a small script; one needs a scope ruling; one needs an email
to the authors.

---

## 7. What I did NOT get to

- **`palczewski`: I did not re-derive the 35 large SMILES atom-by-atom against
  the SI's IUPAC names.** Formula and exact mass are blind to connectivity, and
  on a C117 bis-sulfonamide dye that is the largest residual risk in this review.
  The folder also has **no C10 anchor with a structure**, so brief §9's guards are
  vacuous there, not passing. If this paper is ruled in, that is the first job.
- **`sd_2026`: I did not open Supplementary Fig. 53** (four more organs, 48 rows,
  D-3), **Supplementary Fig. 68** (the second library, D-4), or Supplementary
  Table 5 itself — I corroborated the LC-1 CAS pin only through its printed
  formula, MW and NMR integral sum, not by reading the four CAS numbers.
- **`sd_2026`: S-1 is unresolved structurally.** I could not test whether `LP-01`
  is in the corpus under another name, because no `LP-01` SMILES exists in the
  folder and writing a proprietary benchmark from memory is the exact failure the
  brief documents.
- **`linyx`: I did not verify the 19 in-vivo values** against MOESM9 `Fig. 3`
  W2:X43 / MOESM10 `Fig. 4`, nor the claim that Fig. 3i and Fig. 4c in-vivo are
  one experiment (which the curation argues from ALC-0315's identical triplicate).
  I verified only the 112 in-vitro MOESM4 values. I also did not independently
  re-derive the 3-/4-arm regiochemistry of H4/H9 — mass is identical either way
  and one drawing is the only evidence, which is the Lin-2026 class and the
  largest open structural risk in that folder.
- **`xuq`: I did not read the Figure 5 x-axis group labels myself.** That is the
  folder's own stated blind spot and it is the one place where a single
  misreading would permute all 33 values across 11 structures while every other
  check still passed.
- **I did not look at the other 23 folders' contents** beyond the header/column
  census in §5. The census says the shape is systemic; it does not say which of
  the 23 is actually mis-joined.
- **I did not run `test_checks.py`, `audit.py` or `check_docs.py`** — out of scope
  for a read-only Wave-1 pass, and three of the four folders are unmerged so
  C1–C16 have nothing to see.
