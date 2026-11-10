#!/bin/sh
# Every stage, in dependency order. python3 -I throughout: the sources in
# source/ are untrusted input and an interpreter must not pick code up from
# beside them.
set -e
cd "$(dirname "$0")"
python3 -I stage01_blocks.py          # 01_blocks.csv       blocks from NAMES
python3 -I stage02_fig2a.py           # 06_tns.csv          TNS, not emitted
python3 -I stage03_enumerate.py       # 03_enumerated.csv   all 580 structures
python3 -I stage04_fig2b.py           # 04_fig2b_bars.csv   the 49 values
python3 -I stage05_fig2c.py           # 07_fig2c_points.csv cross-panel check
python3 -I stage06_rows.py            # 05_values.csv       the corpus rows
python3 -I stage07_fig3c_pie.py       # 08_pies.csv         printed-stat check
python3 -I make_reference.py          # reference_structures.csv  the C10 anchor
python3 -I 09_seeded_checks.py        # 09_seeded_checks.json     teeth
python3 -I stage08_crosschecks.py     # 07_qc.json          every check
