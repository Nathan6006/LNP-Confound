# Cross-laboratory prediction of lipid-nanoparticle transfection is structurally confound-limited

Data corpus and provenance machinery for the manuscript of that title.

**Corpus:** `data/data_vitro.csv` (28,885 rows × 79 columns, 42 experiments) and
`data/data_vivo.csv` (1,780 × 79), built from `data/Original_LNPDB.csv` plus six
curations of primary literature. Nothing in `data/` is hand-edited — it is
emitted by `provenance/build_dataset.py` and regenerable from source with
`run_all.py`.

## Layout

| path | what it is |
|---|---|
| `data/` | the emitted corpus and the upstream LNPDB download |
| `provenance/build_dataset.py` | the only writer of `data/data_*.csv` |
| `provenance/audit.py` | checks C1–C14; exits non-zero on any departure from `audit_baseline.json` |
| `provenance/test_checks.py` | 29 seeded mutations, each of which must make the audit fail |
| `provenance/schema_quality.json` | single source of truth for the 7 quality/provenance columns |
| `provenance/CHANGELOG.md` | every merge and repair, with its evidence |
| `provenance/ERROR_CLASSES.md` | the 36 named failure modes this project has hit |
| `provenance/review/` | corpus review rounds: program, findings, overlays |
| `curation/<paper>/` | one folder per paper: digitizers, enumerators, validators, QC JSON |
| `curation/SOURCES.md` | **which paper each folder is, and which PDFs to download** |
| `scripts/`, `figures/`, `results/`, `reports/` | analysis, figures, outputs |
| `CLAUDE.md` | working rules — read this before changing anything |

## Making a fresh clone runnable

The publisher PDFs, the page renders made from them, and two vendored upstream
repos are excluded from version control. `curation/SOURCES.md` lists every one
with its DOI and the exact filename expected, so the only manual step is
downloading the papers you have access to and dropping them in.

```sh
pip install -r environment/requirements.txt
python3 provenance/build_dataset.py     # rebuild data/ from source
python3 provenance/audit.py             # must print AUDIT PASS
python3 provenance/test_checks.py       # 29/29; ~15 min, pre-merge gate
```

## State

Round-1 corpus review is closed: 21 findings, 16 applied, **5 open** and listed
at the end of `provenance/review/round_2026-10-30/FINDINGS.md`. Round 2 is due
after the next 3–5 merges, or before any frozen analysis run.
