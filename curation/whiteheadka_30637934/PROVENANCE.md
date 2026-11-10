# PROVENANCE -- whiteheadka_30637934 -> `KW_2019`

Hajj, Ball, Deluty, Singh, Strelkova, Knapp, Whitehead. *Branched-Tail Lipid
Nanoparticles Potently Deliver mRNA In Vivo due to Enhanced Ionization at
Endosomal pH.* **Small** 15(6), 2019. PMID **30637934**, DOI
**10.1002/smll.201805097**.

**Title and endpoint verified FIRST**, before any chemistry, from the Europe PMC
REST `core` record (`SRC:MED EXT_ID:30637934`, fetched 2026-11-10). The full
published title is the one above -- nothing was truncated -- and the endpoint is
**mRNA delivery / transfection**: in vitro firefly-luciferase RLU in HeLa, and
in vivo bioluminescence and serum EPO. There is no antimicrobial, hemolysis or
cytotoxicity endpoint anywhere in the record.

## What is in this folder

| file | what |
|---|---|
| `00_value_sources.json` | the intake gate, **written first** |
| `00_acquisition.json` | hashes, and the two provenance gaps stated as gaps |
| `01_blocks.csv`, `03_enumerated.csv` | 12 blocks, 11 structures |
| `01_mass_identity.json` | the abstract's own unprinted mass claim, tested |
| `02_corpus_overlap.json` | overlap three ways + the near-miss sweep |
| `02_fig_s6_panelA.csv`, `02_fig_s6_panelB.csv`, `02_fig_s6_geometry.json` | the values and every geometric number behind them |
| `04_table_s1.csv`, `04_table_s2.csv`, `04_tables_qc.json` | tier-4 SI tables and the pKa cross-checks |
| `05_values.csv`, `05_values_qc.json` | all 26 located values, 11 emitted, 15 with a per-row reason |
| `06_rows.csv`, `06_rows_qc.json` | the 11 corpus rows in the 79-column schema |
| `07_qc.json`, `07_seeded_checks.json` | the rollup, and every check clean **and** seeded |
| `08_contradictions.csv` | 10 rows, including 2 of my own errors with what retired them |
| `DECISIONS.md` | 5 `OPEN` entries; the merge is blocked on them by design |
| `make_reference.py`, `reference_structures.csv` | the C10 anchors and their gate |

Run order: `stage01_enumerate` -> `stage02_fig_s6` -> `stage03_corpus_overlap`
-> `stage04_tables` -> `stage05_values` -> `stage06_rows` ->
`stage07_seeded_checks` -> `stage08_qc`. Each with `python3 -I`.

## Lineage of every emitted number

**Structures** (11): built by `curation/lib/assemble` only -- `esterify` for the
ester bond (the acid loses its OH, the alcohol keeps its oxygen) and `attach` on
**atom-mapped** attachment atoms for all four N-C bonds. Reagents read from the
SI Methods p.2 by chemical name: *3,3-diamino-N-methyldipropylamine* (306) plus
an alkyl acrylate per tail. No `RemoveAtom`, no index-0 bonding, no string
surgery. 4 Michael arms asserted on every product; 0 radical atoms.

**Values** (11): **tier 4b**, the PDF's own path coordinates. No page was
rendered to a raster and no pixel was read anywhere in this curation. Bar tops
are the 8.51 pt bar **outline** rules; error-bar caps are 4.19 pt and pattern
tiles 2.80 pt, so bar / whisker / marker separate by **exact width** rather than
by a coverage threshold. The y scale comes from 6 detected tick marks through
`curation/lib/axis.calibrate` (residual 0.0008 pt, worst leave-one-out 0.0015
pt, against 2.5 / 1.5 pt tolerances), and the outline's legs start exactly on
the y=0 tick, so no origin is supplied separately.

**Which bar is which lipidoid**: from the 14 x-axis tick marks (pitch spread
0.061 pt; every bar centre on a tick to 0.059 pt). A seeded one-slot rotation
scores 0/11 on the absolute nearest-tick test.

## Anchor posture -- RECORDED BY HAND

`provenance/structure_check.py`, `check_anchor_independence.py` and
`check_anchor_coverage.py` were **NOT run**, by instruction: they write tracked
files other agents read. So the posture is stated here.

- **The paper draws no structure at all.** All 17 SI pages were enumerated:
  no reaction scheme, no structure panel, no drawn lipidoid. The main article,
  which might contain one, is behind a subscription.
- Therefore the anchor is derived from the **published chemical NAMES** on SI
  p.2, which is external ground truth of a different kind from a drawing.
  `reference_structures.csv` holds **three** anchors -- `306O10` (linear,
  mid-series), `306Oi10` (the one branched tail, the paper's headline compound)
  and `306O18` (the longest) -- with the header exactly
  `lipid_name,reference_SMILES,source_figure,read_on,read_by,expected_mismatch,note`.
  Three rather than one because a single anchor at one chain length cannot
  detect an off-by-one that scales with n.
- **Independence is established by PROVENANCE, not by byte difference.** F-037
  is the precedent and the trap: a correct hand-derivation is byte-identical to
  a correct enumeration once canonicalised, so testing for byte difference
  rewards a wrong anchor. What makes `make_reference.py` independent is
  checkable by reading it, and is also **executed**:
  `assert_no_assembler()` walks the module's own AST at run time and confirms it
  imports no assembler module and calls none of `RWMol`, `AddBond`,
  `RemoveAtom`, `CombineMols`, `EditableMol`, `ReactionFromSmarts`,
  `RunReactants`, `esterify`, `attach`. Its clean output is recorded in
  `07_qc.json`. The anchors are source string literals typed head-first, the
  opposite traversal from `stage01`'s arm-first build.
- **The first version of `assert_no_assembler` could not pass**: it
  substring-scanned a file that contains its own blocklist and reported all
  eight forbidden tokens on a clean run. AST inspection replaced it. A check
  that cannot pass is as useless as one that cannot fail.
- **The gate was proven by seeding.** Clean: 3/3 anchors match by FULL
  InChIKey, exit **0**. Seeded `swap_tail_length` -> exit **1**; seeded
  `move_an_ester` -> exit **1**; seeded `drop_an_arm` -> exit **1**. Exit codes
  read from `subprocess.run(...).returncode`, never through a pipe. The
  `move_an_ester` seed has the **identical molecular formula and identical exact
  mass to every decimal place**, so the gate demonstrably tests connectivity and
  not mass.
- `C18`'s copy-posture question, had it been run, would be answered: these
  anchors ARE expected to be canonically identical to the enumeration, and that
  is the correct outcome, not a copy. The folder has a generator script and that
  script's independence is machine-checked.

## What was NOT done, and why

- **Nothing was merged.** `data/` is untouched; `06_rows.csv` is not registered
  in `build_dataset.py`'s `MERGE` list; five `DECISIONS.md` entries are `OPEN`,
  which blocks the merge by design.
- **No file outside this folder was created or modified.** No `curation/lib`
  edit, no shared Playwright browser, none of the three tracked-file-writing
  checks.
- **No `curation/lib` threshold was loosened.** `reagents.smi` refused
  `DMPE-PEG2000` (two structures in the corpus) and the refusal was respected
  and reported, not worked around; this paper prints `C14-PEG2000`, which
  resolves uniquely. `axis.calibrate` refused a seeded 6 pt tick error outright.
- **Figures S1, S3D, S4 and S5 were not digitised.** S1 is Cy5-mRNA
  biodistribution, which is label distribution/uptake and is excluded from the
  corpus by definition. S3D, S4 and S5 are in vivo but resolve only 2-3
  lipidoids each, S4's bars are **six-organ stacked totals** (measured, not
  assumed -- `08_contradictions.csv` X-4) and S5 is a different cargo and
  readout class. Queued as `D-2`.
- **Figure S6B was not emitted.** It is a second plate that disagrees with
  Figure S6A by 0.50x-1.28x and rank-inverts the four shared lipidoids
  (Spearman rho = -0.20) at the same nominal dose. Its 12 values are recorded
  in `05_values.csv`. Neither panel was averaged into the other and neither was
  rescaled by the other. Queued as `D-3`.

## Residual risk, stated plainly

Three things cannot be settled from what we hold, and all three are in
`DECISIONS.md`: which excipient takes 16 mol% and which 46.5 (`D-1`); the exact
branch of isodecyl acrylate (`D-4`, which `CLAUDE.md` already lists as
unverifiable); and the fact that every emitted value was measured on **crude**
1:4 reaction product while the stored structure is the nominal tetra-substituted
species (`D-5`). The last is flagged `structure_assumed = 1` on all 11 rows with
the reason on the row.
