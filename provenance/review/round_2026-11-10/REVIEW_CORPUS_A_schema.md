# REVIEW — Wave 2 / CORPUS-A (schema, vocabulary, normalisation, identity)

Read-only adversarial review across **all 29 staged folders at once**, plus the
live corpus where a staged folder's claim is only checkable against it.
Everything below is a measurement made in this session with `python3 -I`.
Nothing was edited; this file is the only thing written.

Scope read first: `provenance/EXTRACTION_BRIEF.md`,
`provenance/schema_quality.json`, `provenance/ERROR_CLASSES.md`,
`curation/lib/emit.py`, `provenance/build_dataset.py` (merge path),
`provenance/audit.py` C12, `provenance/reagent_registry.py`,
`provenance/run_queue.json`, all 29 `packet_*.md`.

**Headline.** The 29 folders are in better shape on the things round 1 and 2
were about (flags, vocabulary, emit discipline) than I expected, and worse on
three things nobody owns: **a staged folder's emitted file has no schema at all
in 17 of 26 cases and nothing in the repository can inspect it before it is
merged**; **`SL_2025` is a second inter-folder `Experiment_ID` collision, and
the free-ID check every folder runs is structurally incapable of finding it**;
and **`audit.py` C12 partitions only on `Model_type`, so the corpus carries two
undeclared per-ORGAN z-scores today — `ZH_2023` 325 rows and `BL_2023` 53 rows,
378 in total, both declaring `value_normalization_unit = experiment`.**
The last one is in `data/`, not in a staging folder.

The safety net holds: **all 29 folders have a `DECISIONS.md` with ≥5 `STATUS:
OPEN` entries (220 in total), so `build_dataset.py`'s gate would refuse every
one of them today.** Every schema finding below is therefore "corrupts or
crashes the merge when the decisions are resolved", not "already in `data/`" —
except the C12 finding, which is live.

**Added after the coordinator's mid-task correction, verified from the files:**
the brief's claim that four folders are structureless was **false on two of
four** (`linyx_41851281` has 138/140 structures under `smiles`, 101 new to the
corpus; `palczewski_41654520` 59/64 under `canonical_SMILES`, 30 new), and my
first draft of §1.1–§1.2 **reproduced the brief's rule, got the same four
folders, and reported the agreement as confirmation**. §1.2 is corrected and
the class is **N-16 / error class 51**: a discovery rule standing in for a
declaration, which has now failed three times in the same place (filename →
literal column name → tuple of column names). The durable fix is **C31**, a
per-folder `00_schema.json` the producer writes and the consumer resolves
exactly or refuses. One more measurement fell out: `brief_staged.py`'s tiebreak
for `mitchellmj_38982196` is **non-deterministic** — six runs split 3/3 — and
the committed packet briefed its reviewer on an enumeration intermediate
(**N-17**).

---

## 1. CONFIRMED (your three numbered claims, measured)

### 1.1 Emitted-row files have no consistent name — CONFIRMED, and the
### schema-based replacement rule you proposed is itself degenerate

Identifying by schema (`IL_SMILES` present AND a column starting
`Experiment_value`) finds **26 files across 25 of 29 folders**, under **6
different names**:

| name | folders |
|---|---|
| `06_rows.csv` | 13 |
| `05_values.csv` | 8 |
| `09_rows_vitro.csv` + `09_rows_vivo.csv` | 1 (siegwart, two files) |
| `06_emitted_rows.csv` | 1 |
| `05_values_flagged.csv` | 1 |
| `07_QG_2026_rows.csv` | 1 |

**But the rule does not do what it says.** `Experiment_value_within_model` is
one of the seven quality columns and `curation/lib/emit.py` writes it for every
folder, so "a column starting `Experiment_value`" is satisfied by *every* file
that carries the flag block. The rule reduces to **"has `IL_SMILES`"**. Proof:
it matches `curation/mitchellmj_38982196/03_enumerated.csv`, a 17-column
*structures* file with no values at all, because that folder put the quality
columns on its enumeration output. And **9 of the 26 matched files have no bare
`Experiment_value` column at all** (`andersondg_41203968`, `chamolstad_41691524`,
`khanof_36652629`, `li_39883839`, `lius_40500862`, `mckinlaycj_39105751`,
`mitchellmj_38982196`, `xuq_32438474`, `zhangz_40433897` — 818 rows).

The rule that actually discriminates is **`IL_SMILES` AND `Experiment_ID` AND
`Model` AND a bare `Experiment_value`** — which only **16 of 26** files satisfy.

### 1.2 Four folders with values and no `IL_SMILES` — **REFUTED on 2 of 4,
### re-measured independently.** The 255 figure was the generator's defect.

The coordinator's mid-task correction is **correct, and I verified it from the
files rather than taking it** (parser control asserted first: the registry DSPC
SMILES keys to `NRJAVPSFFCBXDT-HUESYALOSA-N`, so a `None` below means
unparseable, not a broken parser):

| folder | structure column | rows | parseable | blank | distinct | already in corpus | **NEW** |
|---|---|---|---|---|---|---|---|
| `linyx_41851281` | **`smiles`** | 140 | **138** | 2 | 119 | 18 | **101** |
| `palczewski_41654520` | **`canonical_SMILES`** | 64 | **59** | 5 | 30 | 0 | **30** |
| `xuq_32520527` | none in `05_values.csv` | 33 | — | — | — | — | — |
| `sd_2026` | none anywhere | 18 | 0 | — | 0 | — | 0 |

For `xuq_32520527` I tested the join rather than assuming it: `05_values.csv`
has **11 distinct `lipidoid`** values and `03_enumerated.csv` carries
`IL_SMILES` keyed on both `lipidoid` and `IL_name`, each covering **11/11**
(`head_name` covers 0/11 and is the wrong key). So all **33 rows are
recoverable by a bijective join**, and the recovery is checkable.

**So the honest figures are:** rows with no ionizable-lipid structure anywhere
= **18** (`sd_2026` only), not 255. Rows recoverable only by a declared join =
**33**. Rows whose structure column is present but blank on the row = **7**
(`linyx` 2, `palczewski` 5) plus the 70 below. And **131 structures new to the
corpus (101 + 30) would have been written off** by the claim as stated.

My own §1.1 and §1.2 as first drafted inherited this: I ran the discovery rule
I was given, got the same four folders, and **reported the agreement as
confirmation.** That is the error the round-1 ledger names five times — a
conclusion drawn from a secondary artifact rather than the primary one — and
here the secondary artifact was *the reviewer's own brief*, which is class 50's
shape (ground truth taken from the project's own narrative) one level worse
than class 50 because it was being used to scope a review rather than score
one. **Corrected above; see N-16 for the class.**

**A second form of the shape does exist**, and this part stands: the column is
present and **blank on the row**.

- `levkinpa_29558113/06_rows.csv` — **40 of 45 rows** have a blank `IL_SMILES`.
  They are its `emittable = False` rows (head groups A3/A4/A8/A13/A14 that the
  SI never defines). Honest, and disclosed in a bespoke column.
- `xuq_32438474/05_values.csv` — **30 of 120 rows** blank. These are
  `is_library_member = 0` controls: `Lpf2k` ×10, `Control` ×10, `GFP-Cre` ×6,
  `mRNA` ×4. Correctly structureless.

Neither is a defect. The finding is that **a column-presence test reports both
folders as "has structures" and a column-name test reports `linyx` and
`palczewski` as having none** — the two tests fail in opposite directions on
the same question. Corrected total staged rows with no recoverable
ionizable-lipid structure: **18 + 7 blanks = 25 of 7,296 (0.34%)**, with 70
more legitimately structureless (40 undefined head groups, 30 controls) and 33
recoverable by join.

And the mirror-image shape, which nothing names: **892 of 7,296 staged rows
(12.2%) carry no `Experiment_value`** — 818 in files without the column, plus
**74 rows where the column exists and is blank on every row**
(`levkinpa_23909648` 29, `levkinpa_29558113` 45). Both of those folders'
`06_rows.csv` — the file whose *name* means emitted — contains zero
corpus-scale values.

### 1.3 The four queue/file row-count disagreements — your explanation holds
### for two, is the wrong predicate for one, and is simply wrong for siegwart

| folder | queue | file | explanation |
|---|---|---|---|
| `siegwart_40644303` | 505 | 465 | **NOT admissibility.** `09_rows_vitro.csv` 465 + `09_rows_vivo.csv` 40 = **505 exactly**. The scan read one of two files. |
| `lius_40500862` | 4 | 23 | **Your explanation holds.** `05_values.csv`, no `06_rows.csv`; `emitted_as_corpus_candidate` = True 4 / False 19. |
| `xuq_32438474` | 90 | 120 | **Right shape, different predicate.** The 30 excluded rows are `is_library_member = 0` (controls), not an admissibility judgement. |
| `levkinpa_29558113` | 5 | 45 | **Your explanation is wrong about which file.** `05_values.csv` is *also* 45 rows, and the 45-row file is `06_rows.csv` itself. The 5 are marked by a bespoke `emittable` column (True 5 / False 40) that exists in no other folder and in no schema. So this folder writes 40 self-declared non-emittable rows into the emitted-rows deliverable. |

### 1.4 F-062 — nobody "corrected" the misspelling in an emitted row

`protein_adundance_normalized` is used as `Experiment_method` by
`andersondg_26469188`, `dahlmanje_41183194`, `dixonjb_41269090`,
`namjh_39264000`; **no staged folder emits `protein_abundance_normalized`.**
The gate is safe.

Latent hazard, correctly scoped: `curation/mitchellmj_41603124/05_values.csv`
and `curation/dahlmanje_30394729/05_values.csv` carry a **third** vocabulary in
a column called `readout_class`, spelled *correctly* —
`protein_abundance` (2 rows) and `protein_abundance_knockdown` (7 rows). Both
folders are `extracted_zero_admissible`, so nothing is emitted and no gate is
threatened. But if anyone ever maps `readout_class` → `Experiment_method` the
766-row category splits in two.

### 1.5 The seven quality columns — F-001's posture is fully carried, and
### `emit.emit_rows` is doing real work

All 26 emitted files carry **all 7** quality columns. Across **7,296 staged
rows**: **zero** out-of-vocabulary values in any of the five enum columns, and
**zero** rows with a caveat and an empty `curation_note`. `rows_with_a_caveat`
by the F-024 definition = **3,079**. Every one of the 29 folders imports
`emit.emit_rows`. The two defects C13 and `emit.py` exist for do not occur.

`value_is_bound` is used **consistently** across the 6 folders that set it:
`upper` is always an assay/colormap/axis floor, `lower` is always an open top
band or a bar above an axis break. I checked the `curation_note` on one
instance per folder; no folder has the senses inverted. (I provisionally
flagged `levkinpa_23909648`'s single `upper` row as having a note that never
mentions the bound — **that was my error**: the note's last sentence says
"printed as 0 with an SD of 0, i.e. a non-detect: the true transfection
efficiency is at or below the assay's detection limit". Retracted.)

---

## 2. NEW FINDINGS

Ordered by what they cost if nobody acts.

### N-1 (SEVERE, live in `data/`) — `audit.py` C12 partitions only on
### `Model_type`, and the corpus has two undeclared per-ORGAN z-scores

C12's docstring is right about the mechanism and even states the tautology
("pooling blocks each at mean 0 / sd 1 gives mean 0 / sd 1"). Its
implementation partitions on `r["Model_type"]` only, and bails on
`len(usable) < 2`. Every in-vivo experiment in the corpus has a **single**
`Model_type`, so **C12 checks zero in-vivo experiments.**

I re-ran C12's exact logic (same 5e-3 / 0.05 pooled gate, same 1e-3 / 5e-3
per-part gate, same `n >= 3` usable rule) with the partition key changed to
`Model_target`:

```
VITRO  Model_type  : checked 3  flagged 2  undeclared 0        <- what the project reads today
VITRO  Model_target: checked 0  flagged 0  undeclared 0
VIVO   Model_type  : checked 0  flagged 0  undeclared 0
VIVO   Model_target: checked 3  flagged 3  undeclared 2  rows 378
    BL_2023   rows=53  blocks={lung_epithelium:49, muscle:4}  declared=['experiment']                 <-- UNDECLARED
    NG_2025_2 rows=23  blocks={liver:12, muscle:11}           declared=['experiment_x_model_type']
    ZH_2023   rows=325 blocks={liver:163, spleen:162}         declared=['experiment']                 <-- UNDECLARED
```

Exact residuals, `data/data_vivo.csv`:

- **`ZH_2023`, 325 rows.** liver n=163 mean −4.294e-11 sd **0.999999999975**;
  spleen n=162 mean +3.704e-11 sd **1.000000000005**. Declared `experiment` on
  all 325 rows. The liver–spleen offset is deleted and the row says it is not.
  `Experiment_value_within_model` is identical to `Experiment_value` to
  **exactly 0.000e+00** (single `Model_type = Mouse_B6`), so the derived column
  offers no recovery either — **error class 41 in its pure form**.
  Independent corroboration from the project's own record: CLAUDE.md §3 states
  that `ZH_2023`'s trailing-space padding "tracked the organ exactly (every
  liver row padded, every spleen row not), which is the fingerprint of two
  source tables concatenated upstream." Two source tables, one z each.
- **`BL_2023`, 53 rows.** lung_epithelium n=49 sd **1.000000000067**; muscle
  n=4 sd **1.0 exactly**. Declared `experiment`. Here the partition is
  confounded with `Route_of_administration` (intratracheal / intramuscular)
  and dose (5.0 / 10.0 µg), which argues these are two experiments rather than
  that the declaration is right. Either way `experiment` is false.

**`NG_2025_2` is the built-in positive control and it passes:** identical
signature, declares a per-part unit, and the generalised check correctly does
*not* flag it. So the generalisation is not a detector that fires on
everything — it flags 3 of 3 in-vivo candidates and exonerates the one that
declared itself.

Two further consequences:

- `NG_2025_2`'s z is per **`Model_target`** while its only declarable option
  was `experiment_x_model_type`. The schema enum (`experiment`, `source_file`,
  `experiment_x_model_type`, `""`) **cannot express the unit this corpus
  actually used**, so the one honest declaration had to borrow the wrong name.
- C12 never even reaches `NG_2025_2` (single `Model_type` → `continue`), so
  that declaration has never been confirmed by any check.

Severity: **changes values.** 378 corpus rows whose between-organ offset is
gone and whose provenance column asserts it is not. This is F-016 on a second
axis, found by generalising F-016's own check.

Confidence: **very high.** Two independent gates (pooled-scale gate, per-part
gate) plus a passing control plus the documented two-table upstream structure.
What I did **not** establish: whether the per-organ z is correct for these two
papers' data (it may be unrecoverable, as for `WT_2026`). I am reporting that
the declaration is false, not which half to repair — per §6.

### N-2 (SEVERE at merge) — nothing in the repository can inspect a staged
### folder's emitted file, and the merge fills every absent column silently

`build_dataset.py::merge_curations` does:

```python
row = {c: r.get(c, "") for c in cols}
...
(vitro if r["Model"] == "in_vitro" else vivo).append(row)
```

Measured against the 79-column corpus header:

| emitted files | column count |
|---|---|
| 9 files | exactly 79 |
| 3 files | 74 (missing exactly the 5 RDKit descriptors the builder recomputes — harmless) |
| **14 files** | **46, 46, 40, 40, 39, 38, 38, 35, 30, 29, 28, 26, 26, 18** |

`zhangz_40433897/06_rows.csv` has **18** columns; **66 of the 79 corpus columns
are absent** and would merge as empty strings with no error. `r.get(c, "")` is
a silent default on a 61-column gap — error class 28 at schema scale.

Two distinct failure modes if a folder is added to `MERGE` today:

- **Loud (acceptable):** 9 folders have no `Experiment_ID` column and 10 have
  no `Model` column; both are indexed with `r["..."]`, so the merge raises
  `KeyError`.
- **SILENT (not acceptable):** `lij_42130331/07_QG_2026_rows.csv` **has** a
  `Model` column and it is **blank on all 300 rows** (its D-1 records that the
  paper never states in vitro vs in vivo). `r["Model"] == "in_vitro"` is False,
  so the `else` branch sends **all 300 rows into `data_vivo.csv`** — 16.6% of
  the 1,803-row in-vivo pool, from a paper that does not say it is in vivo.
  Nothing warns.

The structural point: `audit.py` runs on `data/`; `build_dataset.py` only sees
files named in `MERGE`. **29 staged folders, 7,296 rows, 220 open decisions,
and no check in `provenance/` reads any of them.** The brief's §11 instruction
to "read the 79-column header **live** from `data_vitro.csv`" was followed by
12 of 26 files and nothing noticed.

Severity: **changes values / crashes the build.** Confidence: **high** (read
the merge source; the routing is one expression).

### N-3 (SEVERE for the manuscript) — `siegwart_40644303` would put a
### 360-structure SAME-LAB pair into Result 3, larger than the entire existing
### shared-lipid pool

Measured by full InChIKey over `09_rows_vitro.csv` (465 rows, 465 distinct
structures, 0 unparseable):

```
overlap with KZ_2016 : 360 structures   (PMID 26729861)
overlap with SL_2020 :   8              (PMID 33305471)
overlap with ZC_2023 :   1
```

Both overlapping Siegwart-lab papers, confirmed from
`provenance/esummary_cache.json`, not inferred from initials:

- `KZ_2016` = PMID 26729861, first author **Zhou K**, **last author Siegwart DJ**
- `SL_2020` = PMID 33305471, first author **Lee SM**, **last author Siegwart DJ**
- the folder is `siegwart_40644303`, i.e. a third Siegwart-lab paper

Against the **current** in-vitro corpus, computed the same way:

```
19 shared-structure source-pairs; 321 shared-structure instances in total
largest pairs: AA_2008<->KW_2014 135 | LL_2026<->LL_2026_2 70 | KW_2014<->NC_2024 55
               JC_2023<->YX_2024 30 | KZ_2016<->SL_2020 12 | ... (14 more, each <=4)
```

So merging this folder adds **one pair of 360** — **2.67× the current largest
pair and larger than all 19 current pairs combined.** Stage 5's primary metric
is the concatenation of all cross-source-pair points, so a single **same-lab**
pair would be **>50%** of the pooled set. The frozen Result 3 (ρ = 0.083,
n=236) already flags that its top pair is 57% of the points; this would make
the top pair same-lab and bigger than everything else.

F-068 recorded that `Experiment_ID` is a screen key and not a lab key, and that
**1 of 22** frozen Stage 5 pairs is same-lab with n=1 and therefore no ρ. This
finding is that F-068 is about to go from n=1 to n=360.

**Credit where due:** the folder's `DECISIONS.md` **D-10** asks exactly this
question ("This experiment shares all 465 of its lipids with KZ_2016. Is that a
shared-lipid pair, or a duplicate?"), `STATUS: OPEN`. What does not exist
anywhere is the quantification above, and there is no lab key in the schema to
express the answer with.

Severity: **changes a headline manuscript number.** Confidence: **very high**
on the measurement; the *interpretation* (that same-lab concordance is higher
and would drag the pooled ρ up) is an inference, supported by the project's own
within-lab datapoints (molea 0.30, split-half ceiling 0.400) but not measured
here — I cannot compute it because the values are not yet on one scale.

### N-4 (SEVERE, needs adjudication before merge) — 105 of `siegwart_40644303`'s
### 488 lipid names carry a DIFFERENT structure from the same name in `KZ_2016`,
### and 60 of them at identical molecular formula

Of 488 distinct (`IL_name`, InChIKey) pairs in the two siegwart row files:

```
name present in KZ_2016 and structure AGREES : 368
name present in KZ_2016 and structure DIFFERS: 105     <- all 105 against KZ_2016, no other experiment
name absent from the corpus                  :  23
siegwart's structure found under a DIFFERENT KZ_2016 name:  0   (so it is NOT a relabelling/permutation)
of the 105, SAME molecular formula: 60   (invisible to C1-C9 and C11)
```

The 105 decompose into **exactly two orthogonal axes**, which is the
diagnostic part no per-paper reviewer could see:

| axis | disagreement |
|---|---|
| three amine cores `3A2`, `5A2`, `6A1` | **all 31 tails each** (93 names) |
| the `SC3` tail | **all 15 amine cores** (15 names; `3A2/5A2/6A1-SC3` double-counted) |

and the 12 other amine cores agree on **30 of 31** tails each (the exception
being their `SC3`).

Formula evidence:

```
3A1-SC3   siegwart C42H74N2O12S3    corpus C48H86N2O12S3     (-6 C, 2-arm core)
5A2-SC3   siegwart C68H123N5O20S5   corpus C78H143N5O20S5   (-10 C, 5-arm core)
6A1-SC3   siegwart C80H143N5O24S6   corpus C92H167N5O24S6   (-12 C, 6-arm core)
5A2-SC2   siegwart C78H143N5O20S5   corpus C78H143N5O20S5    IDENTICAL formula, different InChIKey
6A1-SC4   siegwart C86H155N5O24S6   corpus C86H155N5O24S6    IDENTICAL formula, different InChIKey
3A2-SC1   siegwart C45H81N3O12S3    corpus C46H82N2O12S3     3 N vs 2 N
```

So: the `SC3` axis is a **per-arm chain-length** disagreement (2–3 C per arm);
the `3A2` axis is a **different amine** (3 N vs 2 N); the `5A2`/`6A1` axis is a
**pure connectivity** disagreement at identical formula and identical exact
mass. I also tested the obvious off-by-one hypothesis: siegwart's `<A>-SC3`
equals the corpus's `<A>-SC2` in **0 of 15** cases, and the reverse in **0 of
15**. It is not a shift.

The corpus side is internally consistent with its own changelog — the corpus's
`3A1-SC2` and `3A1-SC3` share formula C48H86 with different InChIKeys, which is
exactly the documented 2026-10-28 patch ("42 `KZ_2016` `*-SC3` lipids held the
SC2 isopentyl structure; source metadata pins neopentyl"). Siegwart's `SC3` is
*six carbons lighter than both*, so it is a third candidate, not the pre-patch
one.

**The folder discloses all three axes**, as `DECISIONS.md` D-4 ("SC3: n-propyl
(this paper's Figure 1b) or neopentyl (corpus KZ_2016)?"), D-6 ("3A2:
piperazine (Figure 1b) or the corpus's two-nitrogen core?") and D-8 ("The
partial-substitution convention on 5A2 and 6A1"), all `STATUS: OPEN`. So this
is honest, gated work.

**What is missing, and is the cross-folder half:**

1. The exposure is **two-sided and the corpus side is unowned.**
   `data/experiment_manifest.csv` records `KZ_2016` as `curated_by_us = no`,
   `curation_folder` empty. It contributes **1,512 vitro + 26 vivo rows** to
   the *frozen* corpus and to the frozen Stage 5 results. If siegwart's
   reading is right on any of the three axes, those rows are wrong and no
   reviewer, per-paper or otherwise, is assigned to them.
2. The decomposition above (105 / 488, three axes, 60 formula-identical) is
   recorded nowhere. D-4/D-6/D-8 each state their ambiguity qualitatively and
   quote one formula.
3. This is the `AC_2025` question at 105× scale, and the §6 rule applies
   verbatim: a name and a structure disagree, and **rewriting the structure
   keeps the value with the name while rewriting the name keeps the value with
   the structure.** I am explicitly **not** adjudicating.

Severity: **changes structures, on both sides.** Confidence: **very high** on
the measurement and the decomposition; **zero** on which side is right.

### N-5 (SEVERE at merge, scale) — 5,112 of the 7,296 staged rows (70%) carry
### ONE ionizable lipid, and it is `DLin-MC3-DMA`

`maohq_39375194/06_rows.csv`: 5,112 rows, **1 distinct `IL_SMILES`**,
`IL_name = DLin-MC3-DMA` on all 5,112, split into two Experiment_IDs
(`LC_2024` 4,680 over 5 cell lines, `LC_2024_2` 432 over 6).

Full staged chemistry census (distinct ionizable lipids per folder, by full
InChIKey):

```
lij_42130331 300 | siegwart 488 | andersondg_41203968 120 | mitchellmj_38982196 100
whiteheadka_38888972 49 | namjh 36 | lius 13 | peer 12 | whiteheadka_30637934 11
khanof 10 | xuq_32438474 9 | li_39883839 9 | levkinpa_29558113 6 | mckinlaycj 5
mitchellmj_41555650 5 | zhangz 3 | levkinpa_23909648 3 | goulart 3
dahlmanje_41183194 2 | dixonjb 2
ONE lipid: andersondg_26469188 (36 rows) | chamolstad (77) | dahlmanje_36701517 (37)
           xuq_32662132 (13) | maohq_39375194 (5,112)
ZERO:      linyx (140) | palczewski (64) | sd_2026 (18) | xuq_32520527 (33)
```

The brief §2 already says "a corpus keyed on ionizable-lipid chemistry gains
that many **however many rows exist**. Do not pad the row count." The
cross-folder measurement is that the run's own headline row count is **70% one
molecule**, and the project's manuscript consequences are specific:

- MC3 is already **3,101 rows** in the corpus (`IL_name` containing "MC3", and
  the same 3,101 by InChIKey over those 4 keys). Adding 5,112 more would make
  **8,213 MC3 rows across ~8 "sources"**, and the two new sources would be
  **indistinguishable from each other and from every existing MC3 row** on
  ionizable-lipid features alone. Result 1 measures exactly that separability.
- `LC_2024` and `LC_2024_2` share their single lipid, so they would form a
  same-lab "cross-source pair" with n=1 lipid in Stage 5.
- Chance baseline moves 1/43 → 1/45 on 5,112 rows of one molecule.
- §8's open decision ("a source whose only shared lipids are benchmarks that
  top every screen can push cross-lab concordance toward +1 for reasons that
  have nothing to do with transferable SAR") applies here not to a source's
  *only* shared lipid but to the **entire source, twice over**.

The folder is honest about what it is: `corpus_admissible` reasons record that
1,080 of its rows were excluded because they **already join `YZ_2022`
one-to-one at Spearman 1.000000**, and that 3,456 more are physicochemical.
Nothing is hidden. The finding is that `rows_emitted = 5112` enters the queue's
headline with no chemistry qualifier attached.

Severity: **changes manuscript numbers if merged; changes no structure.**
Confidence: **very high** (one `IL_SMILES`, counted).

### N-6 (SEVERE) — `SL_2025` is a second inter-folder `Experiment_ID`
### collision, and the free-ID check cannot see collisions by construction

Enumerating every ID-shaped string from every emitted file, every
`00_value_sources.json` and every `07_qc.json` across all of `curation/`:

```
proposed by MORE THAN ONE folder:
  CS_2018   ['dahlmanje_30216729', 'dahlmanje_30394729']      <- yours
  SL_2025   ['lius_40500862',      'siegwart_40644303']       <- NEW
```

(`KW_2014`, `KZ_2016`, `NC_2024`, `QL_2026` also appear in two folders each,
but in `02_corpus_overlap.json`-style reference contexts, not as proposals.
I checked those and they are references.)

`SL_2025` is real on both sides:

- `siegwart_40644303` **emits** it: `09_rows_vitro.csv` 465 rows `SL_2025`,
  `09_rows_vivo.csv` 17 `SL_2025` + 23 `SL_2025_2`, PMID 40644303.
- `lius_40500862/00_value_sources.json` sets
  `"proposed_experiment_id": "SL_2025"`, and `DECISIONS.md` repeats it. PMID
  40500862.

**The mechanism is the finding, and it is shared with CS_2018.**
`lius_40500862` ran a free-ID check, recorded a *retraction* of its own earlier
collision suspicion in `08_contradictions.csv`, and declared `SL_2025` "FREE,
**verified by enumerating every Experiment_ID in both corpus files**". That
check is correct and complete against `data/`'s 54 IDs — and **structurally
incapable of seeing another staged folder**. So the retraction was retired by a
measurement that could not observe the thing it exonerated. The same shape
produced CS_2018.

For completeness, the 29 staged folders propose these IDs, none of which
collides with `data/`: `AR_2023`, `DG_2025`, `EH_2026`, `EM_2025`, `GN_2025`,
`KG_2025`, `KK_2015`, `KW_2019`, `KW_2024`, `LC_2024(_2)`, `LL_2013`,
`MR_2025`, `MW_2024`, `PL_2018`, `QG_2026`, `SB_2025(_2)`, `SD_2026(_2,_3)`,
`SL_2025(_2)`, `XZ_2020`, `YL_2020(_2,_3)`, `YL_2026`, `YX_2025`.

Severity: **two papers would become one source, or one would silently
overwrite the other's ID.** Confidence: **very high.**

### N-7 (MODERATE, and it disabled Wave 1's ability to find N-6) —
### 13 of 29 packets report `proposed Experiment_ID: none` while the folder
### carries one

`brief_staged.py` line 90 reads `qrow.get('proposed_experiment_id')` from
`run_queue.json`, and 16 of 29 queue entries never set it. Measured against
each folder's own emitted file and `00_value_sources.json`:

```
packets saying "none" while the folder carries an ID: 13
  andersondg_41203968 DG_2025   | dahlmanje_36701517 AR_2023
  levkinpa_29558113   PL_2018   | lij_42130331      QG_2026
  linyx_41851281      YL_2026   | lius_40500862     SL_2025
  mitchellmj_38982196 MW_2024   | namjh_39264000    SB_2025, SB_2025_2
  sd_2026  SD_2026/_2/_3        | siegwart_40644303 SL_2025, SL_2025_2
  whiteheadka_38888972 KW_2024  | xuq_32438474      YL_2020/_2/_3
  xuq_32662132        XZ_2020
```

**Both halves of the `SL_2025` collision are in that list.** The field a Wave-1
reviewer would use to check for a collision is blank for exactly the pair that
collides. "If a packet disagrees with the repository, the packet is wrong" —
it does, 13 times, and the right source is the emitted file, not the queue.

Severity: **a reviewer cannot find a collision from the packet.** Confidence:
**very high** (read the generator line, measured both sides).

### N-8 (MODERATE) — `peer_39912611` has `Cargo` and `Cargo_type` transposed on
### all 144 rows

```
corpus Cargo      ∈ {mRNA, pDNA, siRNA}
corpus Cargo_type ∈ {FLuc, FVII, GFP, RLuc, base_editor, hEPO, peptide_barcode}

peer_39912611  Cargo      = 'FLuc' × 144
peer_39912611  Cargo_type = 'mRNA' × 144
```

Each value is a legal member of **the other** column's vocabulary. Nothing in
`peer_39912611/DECISIONS.md` or `08_contradictions.csv` mentions it, so it is
not a declared convention. This is error class 47 (two names swapped) applied
to two columns instead of two reagents, and it is invisible to any per-column
vocabulary check — both columns are "new terms" and both are individually
plausible.

Severity: **cosmetic in the corpus, breaks any `Cargo`-stratified analysis.**
Confidence: **very high.**

### N-9 (MODERATE) — vocabulary drift: two sentinel strings for "not stated",
### three spellings of one organ, a readout synonym, a cargo synonym, and four
### columns holding another column's values

Full tally of every distinct value in `Experiment_method`, `Cargo_type`,
`Model_type`, `Model_target`, `Experiment_batching`, `Model`,
`Route_of_administration`, `Cargo` across all 26 emitted files, against both
corpus files.

**`Model_target`** — corpus has 9 values; staged folders add 17:

- **`'NOT STATED IN THE SI'`** (`levkinpa_23909648`, in `Model_type` *and*
  `Model_target`) and **`'UNSTATED'`** (`lius_40500862`). CLAUDE.md §3 is
  explicit: *"A blank is 'not stated', not 'clean.'"* Two folders invented two
  different sentinel strings for the one state the schema expresses as a blank.
  With `pd.get_dummies` at model time (§3's stated plan) this becomes a
  one-hot feature named `Model_type_NOT STATED IN THE SI`.
- **`luminescence_normalized` in `Model_target`** — `lius_40500862` and
  `mitchellmj_38982196`. That is an `Experiment_method` value.
- **`eGFP plasmid DNA` in `Model_target`** — `levkinpa_29558113`. A cargo.
- **`spleen_over_liver_ratio` in `Model_target`** — `lius_40500862`. §8 already
  rules that a selectivity contrast is not a poolable readout.
- **`lung` (mitchellmj_41555650) vs `lungs` (peer_39912611) vs the corpus's
  `lung_epithelium`** — three spellings of one organ across the staged set.
- Genuinely new and legitimate but undeclared: `heart`, `kidney` (two folders
  each), `serum`, and seven lymphatic sites from `dixonjb_41269090` including
  `lymph_node_popliteal_contralateral_nondraining`.

Worth noting on the last: `dixonjb_41269090` **excluded** 56 FACS-subset rows
with the reason *"The corpus `Model_target` vocabulary is organ-level … and has
no term for a FACS-sorted subset inside an organ"*, and then emitted seven
anatomical sub-sites into that same column. Defensible (site ≠ cell type), but
6 of its 28 emitted rows are `_nondraining`/`contralateral` **negative
controls** sitting inside the z block at means −1.53, −1.65, −1.55. A negative
control inside the standardisation shifts every other row's z.

**`Model_type`** — corpus has 22; staged add 14. Most are honest new cell lines
(`ARPE-19`, `B16F10`, `Neuro2a`, `PC3`, `4T1-Luc`, `MC38-Luc`,
`primary_human_CD8_T_lymphocyte`, `primary_rat_cortical_neuron`). Two hold
another column's values: **`in_vitro`** (`xuq_32438474`) and **`in_vivo`**
(`mitchellmj_41555650`), which are `Model`/`Model_target` values, not a cell
line or strain. Plus `C57BL/6 mouse` (`lius_40500862`) against the corpus's
`Mouse_B6` — a free-text restatement of an existing category.

**`Experiment_method`** — corpus has 6; staged add 3, two of which are
synonyms rather than new kinds:

- **`radiance_normalized`** (`peer_39912611`) — IVIS radiance *is* the
  corpus's `luminescence_normalized`. A synonym split.
- **`siLuc_knockdown_normalized`** (`zhangz_40433897`) — the corpus already
  has `LRP6_knockdown_normalized`, so this convention generates **one term per
  silenced gene** and no two knockdown experiments will ever share a method
  value.
- `GFP_fluorescence_normalized` (`goulart_41138890`) — genuinely new in kind.

**`Cargo_type`** — corpus has 7; staged add 5. **`EPO`**
(`andersondg_26469188`) against the corpus's existing **`hEPO`** is a synonym
split of the same cargo. `Cre`, `Thy1.1`, `aVHH` are genuinely new.
`mRNA` is N-8's transposition.

**`Experiment_batching`**: only `individual`, used by every folder that sets
it. No drift. `barcoded` unused. `Model`: `in_vitro`/`in_vivo` only, no drift.
`Route_of_administration`: one addition, `intravenous_retro_orbital`
(`peer_39912611`) — a genuine refinement, undeclared.

Severity: **silently splits categories at model time.** Confidence: **very
high** (exhaustive tally).

### N-10 (MODERATE) — the `_N` suffix has become an escape hatch, with five
### different unstated criteria, and each use adds a "source" to Result 1/4

§3 defines `_2` as "a second screen from one paper". Five staged folders use
it, and the criterion differs every time:

| folder | split | what distinguishes the blocks |
|---|---|---|
| `namjh_39264000` | SB_2025 (72) / SB_2025_2 (36) | ear/intradermal/RLuc/5 µg vs serum/intramuscular/hEPO/10 µg — **genuinely a different assay** |
| `maohq_39375194` | LC_2024 (4,680) / LC_2024_2 (432) | training library vs de-novo validation library, disjoint composition levels — defensible; 5 of 6 cell lines shared; 1 of 1 lipid shared |
| `xuq_32438474` | YL_2020 (48) / _2 (36) / _3 (36) | Fig 3B (4 cell lines, %GFP+) / S24C (HeLa-DsRed, %GFP+) / S24F (HeLa-DsRed, %DsRed+) — **_2 and _3 are the SAME cohort, two readouts** |
| `sd_2026` | SD_2026 (12) / _2 (3) / _3 (3) | **identical `readout` on all three**; no field in the file distinguishes them |
| `siegwart_40644303` | SL_2025 (17 vivo) / SL_2025_2 (23 vivo) | **identical on every metadata field I can read**: in_vivo, Mouse_B6, liver, luminescence_normalized, mRNA, FLuc, intravenous, 5 µg |

The last two are the problem. A split with no declared distinction, separately
z-scored, is **a per-cohort z wearing a different column** — exactly what the
brief names as error class 22 ("a per-cell, per-organ or **per-cohort** z"),
implemented as two Experiment_IDs instead of two normalisation units. And
because it is an ID and not a normalisation unit, **C12 can never see it**: two
IDs are two experiments, and each is dutifully mean 0 / sd 1.

For this project specifically, every `_N` adds an environment to LOSO (Result
4), a class to the source classifier (Result 1), and moves 1/S. `siegwart`'s
two in-vivo blocks would be a same-lab, same-organ, same-assay "source pair"
whose separability is trivially perfect.

To be fair to siegwart: its **D-1** is *"Are Figure 2b and Figure 3b one
in-vivo experiment or two?"*, `STATUS: OPEN`. So the question is asked. It is
not asked by `sd_2026`, and `xuq_32438474`'s `_3` is explicitly labelled in its
own `readout` field as *"a DOWNSTREAM functional readout of the same delivery
event, not a second measurement of it"* — which argues against giving it its
own Experiment_ID at all.

Severity: **changes S, changes LOSO folds, changes Result 1.** Confidence:
**high** on the measurement, **medium** on `sd_2026` (3 rows each; I could not
find the distinguishing field, which may exist outside `05_values.csv`).

### N-11 (MODERATE) — `18PG` on 972 staged rows contradicts the registry and
### 1,474 existing corpus rows, and the registry's authority is one sentence

Cross-folder reagent check over `IL/HL/CHL/PEG/fifthcomponent` name→SMILES
pairs in all 26 emitted files, by full InChIKey, against
`provenance/reagent_registry.py` (15 names) and against `data/`.

Positive control asserted first and passing: registry `Cholesterol` =
`HVYWMOMLDIMFJA-DPAQBDIFSA-N`, `DSPC` = `NRJAVPSFFCBXDT-HUESYALOSA-N`,
`18PG` = `ZEFGRNLJASLRBZ-GXFQOPFSSA-M`. (My first attempt loaded
`reagent_registry.**json**`, which is the check's *output*, not the registry —
the control correctly failed on "no cholesterol in registry" and I switched to
the `.py`. Recording it because a non-zero entry count would have looked fine.)

Results:

- **a staged NAME carrying more than one structure across folders: 0.**
- **a staged name disagreeing with the registry: 1** —
  **`18PG` in `maohq_39375194`, 972 rows**, stored as DOPG
  (`DSNRWDQKZIEDDB-GCMPNPAFSA-M`) where the registry and **1,474 corpus rows**
  in `WT_2026`/`YZ_2022`/`YZ_2024` carry SOPG (`ZEFGRNLJASLRBZ-GXFQOPFSSA-M`).
  **C11 would fail on merge**, correctly.
- a staged name whose structure is absent from the corpus's set for that name:
  106 — 1 is this `18PG`, 105 are N-4's siegwart names.

The folder discloses it fully (`08_contradictions.csv` C-2, `DECISIONS.md` D-2,
`confidence = medium` + the conflict in `curation_note` on all 972 rows). But
the cross-folder observation is about **the registry**, not the folder: the
registry pins `18PG` to SOPG on the authority of **one sentence in one paper's
Materials section** (`YZ_2022`, PMC9310361), and this paper — *from the same
lab, describing the same 1,080-formulation library* — says DOPG in **three
independent places**: SI Table 2 draws both tails with a cis C=C, SI Table 3
prints `Double_bonds = 2`, and `Raw_Data/1080_LNP_OHE.csv` carries
`Double_bonds = 2` on all 180 `18PG` rows. SOPG is C42H80O10P⁻ and DOPG is
C42H78O10P⁻, so mass *can* see this one. **1,474 already-merged corpus rows
depend on which sentence wins**, and that is not the staging folder's decision
to make. §6's rule applies: go to the source and find which the measurement
was recorded against.

Severity: **changes excipient structures on up to 2,446 rows.** Confidence:
**very high** on the conflict, **zero** on which side is right.

### N-12 (MINOR, but it is class 29 surviving its own guard) —
### `emit.emit_rows`'s literal check is defeated by a constant-returning lambda,
### and `lij_42130331` flags 300 rows `structure_assumed = 0` whose own note
### says two regiochemistries were assumed

`emit.py` refuses `confidence="high"` and says so at length — it makes error
class 29 "unwritable rather than discouraged". It cannot refuse
`lambda i, r: "0"`, which is the same thing with a lambda around it.

`curation/lij_42130331/stage07_emit.py:89`:

```python
structure_assumed=lambda n, rec: "0",
positional_ambiguity=lambda n, rec: "0",
value_is_bound=lambda n, rec: "",
value_normalization_unit=lambda n, rec: "experiment",
```

while the `note` callable (also effectively constant — **2 distinct notes over
300 distinct structures**) reads:

> "…the paper publishes no HRMS, NMR table, SMILES or structure file, so the
> **thiolactone ring-opening regiochemistry and the thiol-epoxide
> regiochemistry are read from the drawing and have no independent check**."

`schema_quality.json` defines `structure_assumed = 1` as "where a specific
structural feature was assumed rather than read from the paper … a quaternization
site the paper draws no structure for, an **ester connectivity inferred from a
single NMR shift**". Two regiochemistries read from a drawing with no
independent check is that condition. So **300 rows assert `structure_assumed =
0` and their own `curation_note` describes two assumed features.** The row still
counts as a caveat (`confidence = medium`), so `emit.py`'s caveat-needs-a-note
rule fires and passes; the *specific* flag is inverted.

For context, three other folders set `structure_assumed = 1` on 100% of their
rows for conditions of the same kind (`andersondg_26469188` 36/36,
`chamolstad_41691524` 77/77, `dahlmanje_36701517` 37/37,
`whiteheadka_30637934` 11/11, `zhangz_40433897` 75/75). So the convention is
not in doubt; `lij_42130331` is the outlier.

Severity: **changes a flag, not a structure or a value.** Confidence: **high**
on the inversion (read the schema, the lambda and the note); **medium** on
whether a reviewer who read the paper would agree the thiol-epoxide
regiochemistry is "a specific feature" rather than library-wide uncertainty.

### N-16 (SEVERE, and it is the class behind §1.2) — the structure column has
### no agreed name either, and a discovery rule keyed on a convention nobody
### enforces fails the same way a filename-keyed one did

I agree with the coordinator's framing and the sibling's phrasing: *a discovery
rule derived from a convention nobody enforces is not a check, it is the
convention's second victim.* Measured, because the claim deserves the same
treatment as the data.

**Across the 29 staged folders' CSV headers: 30 distinct SMILES-like column
names.** Three of them are used as a folder's **primary emitted structure
column**:

| primary structure column | folders |
|---|---|
| `IL_SMILES` | 26 |
| `smiles` | 1 (`linyx_41851281`) |
| `canonical_SMILES` | 1 (`palczewski_41654520`) |
| none in the values file | 2 (`xuq_32520527` — join-recoverable; `sd_2026` — genuinely none) |

Widening to **all** of `curation/`, the coordinator's 8-name list checks out —
all four names I could not find in the staged set exist as real CSV headers in
folders outside it, which is why my 30-name tally missed them:

```
il_smiles       curation/dahlmanje_30216729/03_enumerated.csv, dahlmanje_30394729/03_enumerated.csv
emitted_SMILES  curation/mitchellmj_41603124/03_enumerated.csv
product_SMILES  curation/dahlmanje_29489381/03_enumerated.csv
parent_SMILES   curation/pabst_37602460/{03_enumerated,05_values,07_emit_preview}.csv
```

`pabst_37602460/05_values.csv` matters: that is a **values** file whose
structure column is `parent_SMILES`, so the discovery rule fails there too, not
only on enumeration intermediates. The other four are `03_enumerated.csv` only.

The remaining 22 names are legitimate auxiliary columns (`IL_head_SMILES`,
`alcohol_SMILES`, `arm_ester_SMILES`, `SMILES_with_attachment_map`,
`free_amine_SMILES`, `block_SMILES`, `adduct_fragment_SMILES`,
`tail_acrylate_SMILES`, …). They are not a defect, but they are the reason a
sniffer cannot be made safe by adding candidates: `whiteheadka_38888972` alone
carries **19** SMILES-like column names in one folder, so any tuple of
candidate names large enough to catch `canonical_SMILES` is also large enough
to catch `candidate_SMILES` — which sits beside it in `palczewski_41654520` and
is, by its own name, *not* the chosen structure.

**Error class 51** (see the numbering note below), stated generally:

> **A discovery rule standing in for a declaration.** Three successive scans
> for "the file a curation emits" keyed on a filename (found 14 of 27,
> under-reported new lipids by 14×), then on a literal column name (reported 4
> folders as structureless; 2 were false and 131 new structures would have been
> written off), then on a tuple of column names (still a guess, and it must
> admit `candidate_SMILES` to admit `canonical_SMILES`). Each replacement
> failed in the same place, because the thing being discovered is a
> **convention**, and no scan can distinguish a convention that was followed
> from one that was not. The fix is never a better sniffer: the producer must
> **declare** what it produced, and the consumer must resolve the declaration
> exactly or refuse.

Severity: **the discovery rule scopes every review and every merge.** Its
first failure mis-stated a lipid count by 14×; its second would have discarded
131 structures; both printed a clean, affirmative result. Confidence: **very
high** (measured in both directions).

### N-17 (SEVERE) — `packet_mitchellmj_38982196.md` describes an enumeration
### intermediate, chosen by a NON-DETERMINISTIC tiebreak

Two files in that folder satisfy the discovery rule, and — the part worth
measuring — **both have exactly 100 rows**, which is also the folder's
`rows_emitted`. `brief_staged.py:105`:

```python
exact = [c for c in cands if c[1] == want]
pick  = exact[0] if exact else max(cands, key=lambda t: t[1])
```

So `exact` has **two** entries and the pick is `exact[0]` — the first element
of a list built by iterating a `set()` union of two globs over path strings.
Python randomises `str` hashing per process, so this is not "a silent guess one
rename away" (the coordinator's reading, which assumed the `max` branch); it is
**non-deterministic between runs of the same unchanged repository**. Six
separate `python3 -I` processes:

```
run 1  05_values.csv      run 4  03_enumerated.csv
run 2  05_values.csv      run 5  05_values.csv
run 3  03_enumerated.csv  run 6  03_enumerated.csv
```

3 and 3. And the **committed** packet landed on the wrong one:

```
packet_mitchellmj_38982196.md:5
## Rows read from `curation/mitchellmj_38982196/03_enumerated.csv` — 100 rows
```

`03_enumerated.csv` is a 17-column structures file with **no `Experiment_ID`,
no `Model`, no `IL_name` and no bare `Experiment_value`**. Its Wave-1 reviewer
was briefed on an intermediate as "the rows", by coin flip, and a regeneration
of the packets would silently disagree with the committed one.

This is error class 30's mechanism (*"resolved by file order"*) in a new place,
and it is why I flagged in §1.1 that the discovery rule matches that file at
all. Note the folder is also the one that put the seven quality columns on its
enumeration output, which is what makes the file match.

Severity: **a reviewer is briefed on the wrong file; packet output is not
reproducible.** Confidence: **very high** (read the source, ran it six times,
read the committed packet).

### N-18 (SEVERE for the manuscript, and it is the cost of §1.2's false claim) —
### `linyx_41851281` carries 15 NOVEL CROSS-LAB shared lipids, the thing Result 3
### is short of, and the brief's claim would have discarded them

Having recovered `linyx_41851281`'s structure column, I re-ran the corpus
overlap that §5 had reported clean. `linyx_41851281` has **18 of its 119
distinct structures already in the corpus**:

```
AA_2008  9 | KW_2014 10 | NC_2024  2 | AC_2025  3 | XD_2026  3 | SB_2024 2 | WT_2026 2
SP_2020 1 | SS_2018 1 | YZ_2022 1 | YZ_2024 1 | ZC_2023 1 | ZH_2023 1 | BL_2023 1
LZ_2025 1 | SY_2025 1 | NG_2025_2 1
```

**3 are the universal benchmarks** (`DLin-MC3-DMA`, `ALC-0315`, `SM-102`),
which is the §8 hazard and carries no information. **The other 15 are novel
lipidoids shared with `AA_2008`, `KW_2014` and `NC_2024` under completely
different names:**

| `linyx` name | `AA_2008` | `KW_2014` | `NC_2024` |
|---|---|---|---|
| `H2-OT7(2)` | `80_O10` | `A80_C10` | |
| `H3-OT5(3)` | `109_O14` | `A109_C14` | |
| `H3-OT7(3)` | `109_O10` | `A109_C10` | |
| `H3-OT3(3)` | `109_O18` | | |
| `H3-NT3(3)` | `109_N18` | | |
| `H4-OT5(4)` | | `A113_C14` | |
| `H4-OT6(4)` | | `A113_C12` | |
| `H4-OT7(4)` | | `A113_C10` | |
| `H9-OT5(4)` | | `A306_C14` | `A306_O14` |

This matters more than anything else I found on the positive side:

- `AA_2008`↔`KW_2014` is the **current largest shared-structure pair in the
  corpus (135)** and, per the frozen Stage 5, **57% of Result 3's pooled
  points**. A third, independent source landing on that same catalogue is a
  genuine **cross-laboratory** addition — the opposite sign from N-3, which
  would add a 360-structure *same-lab* pair.
- These lipids are **joinable only on structure**. Three labs give one molecule
  three names (`80_O10` / `A80_C10` / `H2-OT7(2)`), which is F-063 and the
  brief's §5 token-pattern warning in live form: a bare `80` substring join
  returns 94 rows across 85 structures, and `H2-OT7(2)` shares no token with
  either corpus name. **Only the InChIKey finds these 15.**
- My §5 statement that "no cross-folder duplicate structure is a novel
  non-benchmark lipid shared between two different labs" was **false**, and it
  was false *because* I accepted the brief's premise and excluded this folder.
  The brief's error was not merely conservative: it would have written off the
  single most Result-3-relevant folder in the run.

`palczewski_41654520`'s 30 structures are **0** already in the corpus, so it
adds no overlap either way.

Severity: **changes Result 3's available cross-lab evidence, favourably.**
Confidence: **very high** on the 18/15 split and the name mappings (full
InChIKey, parser control asserted). What I did **not** do: check whether the
`linyx` values are on a comparable scale — the folder has **no bare
`Experiment_value` column**, so this is a structure overlap, not yet a
concordance pair.

### N-13 (MINOR) — `corpus_admissible` is mandated by the brief and exists
### under 8 names, 5 value vocabularies, and not at all in 15 of 29 staged
### folders (27 of 46 repo-wide)

Brief §3: *"Put `corpus_admissible` yes/no **with a reason on every row** of
`05_values.csv`, counts summing."* Measured in every folder's `05_values.csv`:

| form | folders |
|---|---|
| `corpus_admissible` yes/no **+ a reason column** | 7 (`andersondg_26469188`, `dahlmanje_41183194`, `dixonjb_41269090`, `goulart_41138890`, `maohq_39375194`, `mitchellmj_41555650`, `zhangz_40433897`) |
| `corpus_admissible` yes/no, **no reason column at all** | 1 (`dahlmanje_36701517`, 109 rows) |
| `emitted` yes/no | 2 (`siegwart_40644303` 505/467, `andersondg_26469188` also) |
| `emitted_to_06_rows` + `reason_if_not_emitted` | 1 (`whiteheadka_30637934`) |
| `emitted_as_corpus_candidate` True/False | 1 (`lius_40500862`) |
| `is_library_member` `1`/`0` | 1 (`xuq_32438474`) |
| `is_library_member` `yes`/`no (benchmark)` | 1 (`sd_2026`) |
| `emittable` True/False (in `06_rows.csv`) | 1 (`levkinpa_29558113`) |
| `refused` (free text) | 1 (`li_39883839`) |
| **nothing of any name** | **15** (`andersondg_41203968`, `chamolstad_41691524`, `khanof_36652629`, `levkinpa_23909648`, `lij_42130331`, `linyx_41851281`, `mckinlaycj_39105751`, `mitchellmj_38982196`, `namjh_39264000`, `palczewski_41654520`, `peer_39912611`, `whiteheadka_38888972`, `xuq_32520527`, `xuq_32662132`, and `levkinpa_29558113`'s `05_values.csv`) |

Reason columns come under three names: `admissible_reason`,
`admissibility_reason`, `reason_if_not_emitted`. Value vocabularies: `yes/no`,
`True/False`, `1/0`, `yes`/`no (benchmark)`, free text.

**Repo-wide count, since the coordinator raised it:** `curation/` holds **46**
`05_values.csv` files, of which **27 lack `corpus_admissible`** — not 28. The
27 are the 14 staged ones above plus `dahlmanje_29489381`, `kooh_42690843`,
`liu_41162324`, `pabst_37602460`, `siegwartdj_42150080`, `songx_39664425`, and
`levkinpa_29558113`/`siegwart_40644303`/`whiteheadka_30637934`/
`whiteheadka_38888972`/`xuq_32438474`/`xuq_32520527`/`xuq_32662132` whose
judgement lives under one of the other seven column names. For those folders
the delivery-only judgement exists **only in prose**, exactly as the
coordinator says.

Where it exists, the quality is genuinely high — `dixonjb_41269090`,
`maohq_39375194`, `mitchellmj_41555650` and `zhangz_40433897` write reasons
that would survive a methods review. The finding is that **no scan can count
admissible rows across the run**, which is why the three row-count
disagreements in §1.3 each needed a different predicate to explain.

Severity: **the run has no machine-readable admissibility ledger.**
Confidence: **very high.**

### N-14 (MINOR) — the run's own row count has three values and no definition
### in code

```
run_queue.json  "rows_staged"                     : 7,045
sum of run_queue.json's own "rows_emitted" entries : 7,207
sum of the 26 schema-identified emitted files      : 7,296
```

`rows_staged` appears in no `.py` file in the repository (I grepped
`provenance/tools/queue_update.py`, `provenance/*.py`,
`provenance/review/*.py`), so it is hand-maintained. The 7,207 → 7,296 delta is
exactly the three §1.3 folders (+40 levkinpa_29558113, +19 lius, +30
xuq_32438474); the 7,045 → 7,207 delta of **162** I could not attribute.

This is error class 39 — *"a quantity with a prose definition and no
implementation gets re-implemented, slightly differently, by every pass that
touches it"* — applied to the headline number of the run that is producing the
corpus. It is exactly the F-024/F-026 shape, which C15 was built to prevent for
`CLAUDE.md` and which does not cover `run_queue.json`.

Severity: **cosmetic now; it is the number that will be quoted in Methods.**
Confidence: **very high.**

### N-15 (MINOR) — `dahlmanje_41183194/06_rows.csv` has a duplicated header

`Experiment_value_within_model` appears **twice** in the header line (35 header
fields, 34 distinct). `csv.DictReader` silently keeps the **last** one; `pandas`
would rename it `Experiment_value_within_model.1`. Ten rows, no value
consequence today, but a duplicated column name in a file destined for a
79-column merge is a silent data-loss path and no check looks for it.

Confidence: **very high** (counted the header tokens).

---

## 3. The checks that should exist and do not

Numbered from C19 (C1–C18 taken). Each states what it asserts, where it lives,
and which finding above it would have caught. Per §6, **none of these is done
until it has been made to fail on purpose.**

### C19 — the z-score's partition key, generalised beyond `Model_type`
**Where:** generalise `audit.py::c12_zscore_unit`.
**Asserts:** for every experiment on the z scale, and for **each** of
`Model_type`, `Model_target`, `Route_of_administration`,
`Dose_ug_nucleicacid`, `Experiment_method`, `Cargo_type`: if every part with
n≥3 is independently mean 0 (<1e-3) and sd 1 (±5e-3), the experiment's
`value_normalization_unit` must declare that partition.
**Would have caught:** N-1 — `ZH_2023` 325 rows and `BL_2023` 53 rows, both
declaring `experiment`.
**Positive control it already has:** `NG_2025_2` (per-organ, declared) must be
flagged-and-exonerated, and `LZ_2025`/`YX_2024` (per-cell-type, declared) must
stay exonerated. Make it fail by re-standardising one organ of a correctly
pooled experiment in a copy.
**Also needed:** add `experiment_x_model_target` (and a generic
`experiment_x_<column>`) to `schema_quality.json`'s enum, because the one honest
per-organ declaration in the corpus had to borrow the wrong name.

### C20 — a pre-merge schema gate on every staged curation
**Where:** new `provenance/check_staged.py`, run over every folder in
`run_queue.json` whose status starts `extracted`, **not** only over `MERGE`.
**Asserts, per folder:** exactly one file carries the emitted-row schema; that
file's header is a **superset** of the 79-column corpus header minus the 5
recomputed descriptors; no duplicated header token; `Experiment_ID`, `Model`,
`Experiment_value`, `IL_name`, `IL_SMILES` present **and non-blank on every
row**; `Model` ∈ {`in_vitro`,`in_vivo`}; every other column either in the
corpus header or listed in a declared-extras file.
**Would have caught:** N-2 (14 files off-schema, 66 missing columns in the
worst, 300 rows silently routed to vivo), §1.1 (6 filenames, 9 files with no
bare `Experiment_value`), §1.2 (892 valueless rows), N-15 (duplicate header).
**Make it fail:** blank one `Model` cell, delete one required column, duplicate
one header token — three seeds, three distinct messages.

### C21 — `build_dataset.py` must refuse a partial row instead of defaulting
**Where:** `merge_curations`, replacing `row = {c: r.get(c, "") for c in cols}`.
**Asserts:** every column the curation does not supply is on an explicit
per-folder allow-list of "the builder computes this" (the 5 descriptors,
`Index`, `LNP_ID`); anything else raises, naming the folder and the columns.
And `Model` is read with `row["Model"]` into an explicit three-way branch whose
third arm raises.
**Would have caught:** N-2's silent half. This is the schema-scale instance of
error class 28 (*"never let a detector miss fall through to a numeric
literal"*) — here it falls through to `""`.

### C22 — `Experiment_ID` uniqueness across `data/` **and** every staged folder
**Where:** new check, or `check_staged.py`. One function, called by every
curation instead of each one re-enumerating `data/`.
**Asserts:** the multiset of proposed IDs over {`data/`'s 54} ∪ {every staged
folder's emitted + `00_value_sources.json`} has no duplicate; a `_N` suffix is
accepted only when the two blocks share a `Publication_PMID`.
**Would have caught:** N-6 (`SL_2025`, two papers), your `CS_2018`, and F-075's
`YZ_2024`. **This is the one check that would have made `lius_40500862`'s own
retraction correct** — it ran the right query against the wrong universe.

### C23 — a derived `_N` split must declare its criterion
**Where:** `check_staged.py`, plus a required `split_criterion` field per
`_N` block in `07_qc.json`.
**Asserts:** for any two Experiment_IDs differing only by a `_N` suffix, at
least one metadata column differs on every row, and the differing column(s) are
named in the folder's QC file. A split with **zero** differing columns raises.
**Would have caught:** N-10 — `SL_2025`/`SL_2025_2` (identical on all 8
metadata columns I could read) and `sd_2026`'s three blocks (identical
`readout`). This closes the hole C19 cannot reach: a per-cohort z expressed as
two IDs is two experiments to C19 and each is honestly mean 0 / sd 1.

### C24 — a controlled vocabulary with a declaration file
**Where:** new `provenance/schema_vocabulary.json`, the `schema_quality.json`
pattern applied to the six categorical columns; enforced by `check_staged.py`
and by `audit.py` over `data/`.
**Asserts:** every value of `Experiment_method`, `Cargo`, `Cargo_type`,
`Model`, `Model_type`, `Model_target`, `Route_of_administration`,
`Experiment_batching` is either in the file or in the folder's own
`new_terms` list with one sentence of justification. Plus two cheap
cross-column rules that catch the cases a per-column check cannot:
**(a)** no value of column X may be a member of column Y's vocabulary
(catches N-8's `Cargo`/`Cargo_type` transposition and N-9's four
`Experiment_method`-value-in-`Model_target` cases);
**(b)** a new term within edit distance 2, or differing only by a trailing `s`,
of an existing term must be justified explicitly (catches `lung`/`lungs`,
`EPO`/`hEPO`, and would have caught `protein_adundance`→`protein_abundance`
from the other direction).
And the sentinel rule: **"not stated" is a blank**, so a value matching
`/^(NOT |UN)?STATED/i` raises.
**Would have caught:** N-8, N-9 in full.

### C25 — one `corpus_admissible` column, one vocabulary, counts summing
**Where:** `check_staged.py`; add the column to the scaffold in
`provenance/new_curation.py` so a folder cannot invent its own.
**Asserts:** `05_values.csv` has `corpus_admissible ∈ {yes,no}` and a non-empty
`admissible_reason` on **every** row; `count(yes)` equals the emitted-row file's
row count; `count(yes)+count(no)` equals the values file's row count.
**Would have caught:** N-13 (8 names, 5 vocabularies, 15 folders with none),
and it makes the §1.3 row-count reconciliation mechanical instead of a
per-folder puzzle.

### C26 — `rows_staged` must be computed, never written
**Where:** `provenance/tools/queue_update.py`.
**Asserts:** `rows_staged` is derived from the emitted files at write time, and
each entry's `rows_emitted` equals its emitted file's admissible row count.
**Would have caught:** N-14 (7,045 / 7,207 / 7,296). This is C15's principle
(*the expected value is transcribed by hand on purpose, so the pair either
agrees or the build fails*) applied to the machine-state file.

### C27 — the reverse of C11: one structure, one name — and the corpus as its
### own anchor
**Where:** `audit.py`, over `data/` ∪ every staged folder.
**Asserts:** (a) for every `IL_name` appearing in more than one experiment, the
skeleton InChIKeys agree, or the pair is on an explicit "generic combinatorial
label" allow-list; (b) for every reagent *name* appearing in two experiments,
full InChIKeys agree modulo a small alias table (`OF-2`/`OF-02`).
**Would have caught:** N-4's 105 siegwart↔`KZ_2016` disagreements — which no
per-paper reviewer can see, because `KZ_2016` has `curated_by_us = no` and no
curation folder — and N-11's `18PG`. This is F-063 and error class 46, which
`ERROR_CLASSES.md` already marks *"not yet automated … otherwise mechanical"*
and which my run shows is now carrying 105 live instances, 60 of them at
identical molecular formula.
Note the required escape hatch: CLAUDE.md §3 records **288 `IL_name`s that
legitimately carry different structures in different experiments**, so (a)
must be an allow-list check, not a blanket equality. The signal that
distinguishes N-4 from the benign case is that **all 105 collide against one
experiment, from the same lab, on the same library** — so the check should
report *per-experiment-pair* collision counts and gate on a pair, not on a name.

### C28 — the chemistry qualifier on a row count
**Where:** `check_staged.py` and the queue.
**Asserts:** every `rows_emitted` is accompanied by
`distinct_ionizable_lipids`, computed by full InChIKey, and a merge that would
take any single structure above a stated share of its pool prints that share.
**Would have caught:** N-5 — 5,112 of 7,296 staged rows (70%) being one
molecule, and that molecule being the corpus's most-shared benchmark. The brief
§2 states this rule in prose; nothing measures it.

### C29 — `emit.py`'s literal guard must see through a constant callable
**Where:** `curation/lib/emit.py::_resolve`.
**Asserts:** for each flag spec, evaluate it over **all** items; if the result
set has size 1 while the item set has size > 1 **and** the records differ on any
field, raise unless the caller passes an explicit
`uniform_because="<one sentence>"` that is written into `curation_note`.
**Would have caught:** N-12 — `structure_assumed=lambda n, rec: "0"` on 300
rows whose own note describes two assumed regiochemistries. A blanket flag can
be right; it should have to say so. Make it fail by wrapping `"high"` in a
lambda — the exact construct F-011/F-025/F-046 took three rounds to kill.

### C30 — a per-flag note requirement
**Where:** `emit.py`.
**Asserts:** when more than one caveat flag is set on a row, the
`curation_note` must contain a token associated with each (`bound`/`at or
below`/`censored` for `value_is_bound`; `assume`/`ASSUMED` for
`structure_assumed`; `isomer`/`regio`/`position`/`mixture` for
`positional_ambiguity`).
**Caught nothing in this run** — I tested it by hand across all 6 folders that
set `value_is_bound` and every note explains its bound. Proposing it anyway and
saying plainly that it is **unearned today**: one note currently satisfies up
to four flags, which is the shape F-025 had (a note that blamed the wrong
component while the real reason sat in `03_enumerated.csv`).

### C31 — a per-folder `00_schema.json` DECLARATION, and a gate that resolves
### it exactly or refuses
**This is the one that matters most, and it replaces C20's sniffing half.**
I agree with the coordinator: the durable fix is a declaration, not a better
sniffer, for the reason N-16 states — no scan can tell a followed convention
from an unfollowed one.

**Where:** `curation/<folder>/00_schema.json`, written by
`provenance/new_curation.py`'s scaffold so a folder cannot exist without one,
and resolved by `provenance/check_staged.py` (C20) and by
`build_dataset.py::merge_curations` before it reads a single row.

**Shape** (every field required; no defaults, no inference):

```json
{
  "emitted_file": "05_values.csv",        // exactly one, or null
  "structure_column": "canonical_SMILES", // the ONE chosen structure
  "value_column": "Experiment_value",
  "experiment_id_column": "Experiment_ID",
  "model_column": "Model",
  "rows_emitted": 59,
  "experiment_ids": ["XX_2025"],
  "structure_join": null,                 // or {"from","key","to","to_key","to_column"}
  "zero_because": null                    // required iff emitted_file is null
}
```

**Asserts, in this order, raising on the first failure and naming the folder:**

1. `00_schema.json` exists and every field is present. A missing field is not a
   blank.
2. `emitted_file` resolves to exactly one existing file. Not a glob, not a
   list, not a tiebreak — **this is what kills N-17**, because a declaration
   cannot tie with itself.
3. `structure_column`, `value_column`, `experiment_id_column` and
   `model_column` all exist in that file's header, and **the header contains no
   duplicate token** (N-15).
4. `rows_emitted` equals the file's row count **and** the folder's
   `run_queue.json` entry (N-14), and equals `count(corpus_admissible == "yes")`
   once C25 lands.
5. `structure_column` is non-blank and RDKit-parseable on every row, **or** the
   row is covered by `structure_join`, which must be **bijective and total** on
   the declared key — the `xuq_32520527` case, where 11/11 `lipidoid` values
   resolve and `head_name` resolves 0/11, so the join is checkable rather than
   asserted.
6. `experiment_ids` equals the distinct set in the file, and feeds C22's
   cross-folder uniqueness check. **This is what makes C22 possible at all**:
   N-6's `SL_2025` collision is only visible once both folders declare.
7. `emitted_file: null` requires a non-empty `zero_because`, so **a declared
   zero is distinguishable from a discovery failure.** That distinction is the
   whole point: today the two are the same output, which is how `linyx` and
   `palczewski` were reported as structureless.

**What it would have caught:** N-16 and §1.2 in full (`smiles`,
`canonical_SMILES`, the `xuq_32520527` join, `sd_2026`'s honest zero — four
different answers that the current rule flattens into one wrong one); N-17 (no
tiebreak exists); §1.1 (6 filenames, 9 files with no bare `Experiment_value`,
`mitchellmj_38982196`'s two matches); N-2's silent half; N-14; N-15; and it is
the precondition for C22 and C25.

**Error class 51.** I checked the number rather than taking it: the coordinator
said `ERROR_CLASSES.md` "runs to 48", and it **runs to 50** — class 49 is the
per-folder re-implemented hazard (`curation/lib/`) and class 50 is ground truth
transcribed from the project's own narrative. I enumerated every table row
number (1–50, with 18 and 19 each appearing three times in the pre-renumbering
rows) and grepped every `.md` in `provenance/`, `provenance/review/` and
`curation/*/` for `class N` citations: the highest cited anywhere is **50**.
So **51 is the next free number and is unused.** Taking 49 as instructed would
have collided with a class added on 2026-11-09 — which is precisely the failure
the file's own numbering note documents and tells the next author to prevent by
taking `max + 1` *and confirming it is unused*.

**The four ways it must be made to fail on purpose** before anyone trusts it
(per §6: a new check is not done until it has been made to fail, and my first
two attempts at an analogous guard would have reported clean while seeing
nothing). Each in **its own child process**, verdict read from `returncode`,
never from stdout or a pipeline:

1. **Wrong name, right data.** In a *copy*, rename `canonical_SMILES` to
   `candidate_SMILES` in `palczewski_41654520` without touching
   `00_schema.json`. Must fail at assertion 3. This is the seed for N-16: the
   declaration and the file disagree.
2. **Right name, wrong file.** Point `emitted_file` at
   `mitchellmj_38982196/03_enumerated.csv` — which satisfies the old discovery
   rule and has the declared row count. Must fail at assertion 3 on the absent
   `Experiment_ID`/`Model`/`Experiment_value`. This is the seed for N-17, and
   it is the one that proves the gate reads the *header* and not just the row
   count: **a seed that only changes `rows_emitted` would pass assertion 4 and
   stop there**, which is the "seed exits non-zero while moving nothing" shape.
3. **A silent zero.** Set `emitted_file: null` and leave `zero_because` empty
   on `sd_2026`. Must fail at assertion 7. Then set a real `zero_because` and
   confirm it **passes** — because the check that cannot accept an honest zero
   will be disabled by whoever hits it first.
4. **A plausible but non-total join.** On `xuq_32520527`, change
   `structure_join.to_key` from `lipidoid` to `head_name`. Coverage goes 11/11
   → **0/11** and it must fail at assertion 5. The companion seed is the subtle
   one: drop a single row from `03_enumerated.csv` so coverage is 10/11 — a
   gate that only tests "the join returned something" passes that, and the
   agreement count must **move** for the seed to be credited.

Plus the mandatory clean control: all 29 folders with correct declarations must
come back **exit 0**, and the suite must write only to a temp `OUTDIR` —
`00_schema.json` is a tracked file in the folder being mutated, so F-069's
hazard (a seeded copy writing its output into the real curation folder, with a
non-zero exit concealing the write) applies directly. Assert every
`00_schema.json` is byte-identical across two clean runs **after** the seed
suite.

---

## 4. UNVERIFIABLE

- **Which side of N-4 is right.** Needs PMID 26729861's reagent list (which
  thiol "SC3" names), or HRMS/NMR for any `SC3`, `3A2`, `5A2` or `6A1` lipid.
  Both candidate structures parse, canonicalise and carry self-consistent
  descriptors, and 60 of the 105 share molecular formula and exact mass, so no
  composition check can decide it. Siegwart's D-4 names exactly this evidence
  requirement.
- **Which side of N-11 is right** (`18PG` = SOPG or DOPG). Two papers from one
  lab, one sentence of prose against three tabulated fields. Mass can see it;
  authority cannot be settled from the repository.
- **Whether `ZH_2023`'s and `BL_2023`'s per-organ z is recoverable.** I can
  prove the declaration is false; I cannot tell whether the raw per-organ
  signal exists anywhere. Both are `curated_by_us = no` with no curation folder,
  so there is no QC file to consult. `WT_2026`'s precedent (F-003) is that it
  sometimes is.
- **The actual effect of N-3 on Result 3's ρ.** Needs the frozen pipeline run
  on a corpus that includes `SL_2025`, which §2 and §8 forbid doing casually.
  I report the point counts, not the correlation.
- **`sd_2026`'s three-way split criterion.** Nothing in `05_values.csv`
  distinguishes `SD_2026`, `_2` and `_3`; it may be stated in a file I did not
  read (that folder has no `IL_SMILES`, so it is another agent's).
- **Whether `lij_42130331`'s 300 rows are in vitro or in vivo.** Its own D-1
  says the paper does not state it and the main text is paywalled. I verified
  the merge consequence, not the answer.
- **The 162-row gap** between `run_queue.json`'s `rows_staged` (7,045) and the
  sum of its own `rows_emitted` (7,207). I could not attribute it to any entry.

---

## 5. CONFIDENCE, per finding

| # | finding | severity | confidence |
|---|---|---|---|
| N-1 | C12 blind to the organ axis; 378 undeclared per-organ z rows in `data/` | changes values | **very high** on the defect; nil on the repair |
| N-2 | no pre-merge gate; 14 off-schema files; 300 rows silently routed to vivo | changes values / crashes | **high** |
| N-3 | 360-structure same-lab pair would dominate Result 3 | changes a headline | **very high** measurement, **medium** inference |
| N-4 | 105 name↔structure disagreements vs `KZ_2016`, 60 formula-identical | changes structures, both sides | **very high** measurement, nil adjudication |
| N-5 | 5,112 of 7,296 staged rows are one lipid (MC3) | changes manuscript numbers | **very high** |
| N-6 | `SL_2025` inter-folder collision; the free-ID check cannot see it | two papers become one source | **very high** |
| N-7 | 13 packets say `proposed Experiment_ID: none` | disables Wave 1 | **very high** |
| N-8 | `peer_39912611` `Cargo`/`Cargo_type` transposed, 144 rows | breaks stratification | **very high** |
| N-9 | vocabulary drift: 2 sentinels, `lung`/`lungs`, `EPO`/`hEPO`, 4 misplaced columns | splits categories | **very high** |
| N-10 | `_N` as an escape hatch; 2 splits with no declared distinction | changes S and LOSO | **high** (medium for `sd_2026`) |
| N-11 | `18PG` 972 staged rows vs registry + 1,474 corpus rows | changes excipients | **very high** conflict, nil adjudication |
| N-12 | `lij_42130331` `structure_assumed = 0` on 300 rows its note contradicts | changes a flag | **high** |
| N-16 | 30 structure-column names; a discovery rule in place of a declaration (**class 51**) | scopes every review and merge | **very high** |
| N-18 | `linyx_41851281` carries 15 novel **cross-lab** shared lipids (AA_2008/KW_2014/NC_2024) | changes Result 3's evidence, favourably | **very high** |
| N-17 | packet pick non-deterministic; committed packet names an intermediate | reviewer briefed on wrong file | **very high** |
| N-13 | `corpus_admissible` under 8 names, absent in 15 of 29 (27 of 46) | no ledger | **very high** |
| N-14 | three values for `rows_staged` | cosmetic now | **very high** |
| N-15 | duplicated header token in `dahlmanje_41183194` | silent data loss path | **very high** |

**Things I checked and found clean**, stated so the next pass does not re-pay
for them:

- 0 unparseable `IL_SMILES` across all 7,296 staged rows.
- 0 out-of-vocabulary values in the 5 quality enum columns, 7,296 rows.
- 0 rows with a caveat and an empty `curation_note`, 7,296 rows.
- All 7 quality columns present in all 26 emitted files; all 29 folders import
  `emit.emit_rows`.
- `value_is_bound` senses consistent across all 6 folders that use it.
- 0 staged names carrying more than one structure across folders.
- 0 staged folders with an undeclared within-`Model` per-part z (the C19
  generalisation run over all 16 staged folders that have a bare
  `Experiment_value` returns empty — so the detector that found 2 in `data/`
  finds 0 here).
- `Experiment_batching` and `Model` show no drift at all.
- The F-062 misspelling is reused, never "corrected", in every emitted row.
- All 29 folders carry a `DECISIONS.md` with ≥5 `STATUS: OPEN` entries (220
  total) — nothing can merge today.
- Cross-folder duplicate structures: **8 InChIKeys** emitted by more than one
  folder, all 8 already in the corpus. 5 are universal benchmarks (SM-102 in 5
  folders, MC3 in 4, C12-200 in 3, cKK-E12 in 3). The other 3 are
  `whiteheadka_30637934` ↔ `whiteheadka_38888972` (`306O10`, `306O11`,
  `306O12`, plus `306Oi10` = 4 lipids), a genuine **same-lab cross-paper**
  shared set, under two naming conventions (`A306_O10` vs `306O10`) — F-063's
  one-structure-many-names, benign but worth the alias table C27 needs.
  ~~No cross-folder duplicate structure is a novel non-benchmark lipid shared
  between two different labs.~~ **STRUCK — this was wrong, and wrong because the
  scan excluded `linyx_41851281` on the brief's false premise. See N-18: that
  folder carries 15 novel cross-LAB shared lipids with `AA_2008`/`KW_2014`/
  `NC_2024`.** The corrected statement is that the *26 files my scan covered*
  add no new cross-lab chemistry, and the two files it wrongly excluded add 15.
  This tally must be re-run over all 29 folders with the corrected structure
  columns before anyone quotes it.

**Errors of mine, recorded per §6:**

1. I flagged `levkinpa_23909648`'s `value_is_bound` note as not explaining its
   bound. It does, in its final sentence. Retracted.
2. My first reagent cross-check loaded `provenance/reagent_registry.json` — the
   check's **output** — instead of `reagent_registry.py`. It returned 1 "entry"
   and would have reported 0 disagreements. Only the asserted cholesterol
   positive control caught it. Recording it because that is the F-070 shape
   exactly, and a count of 1 looked like a parse quirk rather than the wrong
   file.
3. **The worst one.** §1.1 and §1.2 as first written ran the discovery rule the
   brief handed me, reproduced its four folders exactly, and wrote that up as
   **"CONFIRMED EXACTLY"** — including "no fifth folder is missing the column",
   which is true and irrelevant, because two of the four were not missing a
   structure at all. I did notice and report that the rule was *degenerate*
   (it reduces to "has `IL_SMILES`", §1.1), and then failed to draw the
   obvious consequence: a rule keyed on the literal name `IL_SMILES` cannot see
   a folder that chose a different name. I had the premise and did not take the
   step. ~200 verified rows and 131 corpus-new structures were inside that gap.
   Corrected only because the coordinator sent the correction; I would not have
   found it from the direction I was working.

**Claims of the coordinator's that I checked and that did NOT survive**, per
§6's "treat findings as claims", in both directions:

- *"`ERROR_CLASSES.md` runs to 48"* — it runs to **50**. Taking the instructed
  49 would have collided with a class added 2026-11-09. I took **51** after
  enumerating all table-row numbers and grepping every citation.
- *"28 of 46 `05_values.csv` files lack `corpus_admissible`"* — it is **27**.
- *"my generator breaks the tie by `max(n_rows)`, a silent guess one rename
  away"* — the `max` branch is never reached for `mitchellmj_38982196`, because
  `rows_emitted` equals **both** candidates' row counts, so `exact` has two
  entries and `exact[0]` resolves by `set` iteration order. It is worse than a
  guess: it is **not reproducible between runs** (N-17).
- *"8 structure-column names"* — correct repo-wide, and I initially could not
  find four of them because my tally was scoped to the 29 staged folders. All
  four are real headers in non-staged folders, one of them
  (`pabst_37602460/05_values.csv`, `parent_SMILES`) in a **values** file rather
  than an enumeration intermediate. Both scopes are now reported.

---

## 6. What I did NOT get to

- **Formulation columns.** I never tallied `HL_name`/`CHL_name`/`PEG_name`
  vocabularies, molar-percent sums (C5's tautology-vs-real distinction per the
  brief §10), `IL_to_nucleicacid_*` ratios, or `Dose_ug_nucleicacid` units
  across folders. `mitchellmj_41555650` has a `molratio_sum` column suggesting
  it checked its own; nobody checked them against each other.
- **`08_contradictions.csv` across folders.** I read two files' worth. A
  cross-folder read would show whether the same contradiction is recorded twice
  with different verdicts, and whether any folder's retracted suspicion was
  retired by a measurement as scoped-wrong as `lius_40500862`'s was (N-6).
  That is a targeted, cheap sweep and I would run it next.
- **The 220 open decisions.** I counted them and parsed their STATUS lines with
  `build_dataset.py`'s own parser; I did not read them for cross-folder
  conflict (two folders deciding the same question opposite ways). Siegwart's
  D-4/D-6/D-8/D-10 are the obvious place to start since N-3 and N-4 both land
  there.
- **`reference_structures.csv` / C10 posture.** Out of scope by instruction (I
  did not run `structure_check.py`, `check_anchor_independence.py` or
  `check_anchor_coverage.py`, and did not inspect the anchors).
- **The four "no-`IL_SMILES`" folders' own content** — another agent owns
  `linyx_41851281`, `palczewski_41654520`, `xuq_32520527`, `sd_2026`. I
  re-measured their *structures* (§1.2) because the brief's claim about them was
  the thing under test, and found the second form of the shape (70 blank-cell
  rows in two more folders). I did **not** look at their values, their
  normalisation or their flags — and note that because those four were excluded
  from my emitted-file set on the brief's premise, **they are absent from every
  other tally in this report**: the vocabulary tally (N-9), the quality-column
  audit, the z-partition test (N-1's staged half), the chemistry census (N-5)
  and the cross-folder InChIKey overlap (§5) all exclude up to 255 rows and, now
  that `linyx`'s 119 and `palczewski`'s 30 structures are known, **149 distinct
  ionizable lipids**. Those five checks should be re-run over all 29 folders
  with the corrected structure columns. `linyx_41851281` in particular shows
  **18 structures already in the corpus**, which my cross-folder overlap section
  could not see and which may add to N-3's shared-lipid accounting.
- **§5's "0 staged folders with an undeclared per-part z"** — I can now close
  this rather than leave it open: **neither `linyx_41851281` nor
  `palczewski_41654520` has a bare `Experiment_value` column**, so there is no
  z to partition and nothing was missed. The clean result stands, and the
  reason it stands is that those folders have no z at all — which puts them in
  §1.2's 892-valueless-rows bucket instead.
- **I did re-run the corpus-overlap tally** over the two recovered folders
  rather than deferring it, because it bore on a clean result I had already
  stated. That produced **N-18**, which is the most consequential thing in this
  report on the positive side. The vocabulary tally (N-9), the quality-column
  audit and the chemistry census (N-5) are still scoped to 26 of 28 files and
  should be re-run; I judged them lower-yield than the overlap because neither
  folder can contribute an `Experiment_method`/`Model_target` value to a merge
  without a value column first.
- **Whether the staged folders' new `Model_type` values collide with each other
  semantically** (`4T1-Luc` vs a plain `4T1`; `MC38-Luc` vs the corpus's
  `MC38`). The corpus has `MC38`; `zhangz_40433897` adds `MC38-Luc`. That is
  very likely the same cell line luciferase-transfected, i.e. a fifth instance
  of N-9's synonym-split class, and C24's edit-distance rule would catch it —
  but I did not verify the biology and so did not count it as a finding.
- **No seeded failure test on any proposed check.** C19–C30 are specifications,
  not checks. Per §6 and the brief §10, **none of them is done until it has
  been made to fail on purpose with the right message, with a clean control
  run** — and my own N-1 result was produced by a hand replication of C12's
  logic in a separate process, not by a check that exists.
