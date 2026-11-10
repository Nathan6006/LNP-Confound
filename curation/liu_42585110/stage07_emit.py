#!/usr/bin/env python3
"""emit rows in the 79-column schema

GENERATED STUB -- replace this docstring with what this stage actually does,
and say where each number came from.

Run from this folder:   python3 -I stage07_emit.py

Imports from curation/lib, which holds ONE audited implementation per hazard
class. Do not copy a geometry helper in here: emit already hold the audited versions.
"""
import json
import os
import sys

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
import emit  # noqa: F401


# ---------------------------------------------------------------------------
# HAZARD: THE SESSION SCRATCHPAD IS SHARED BY EVERY AGENT IN THE RUN.
#
# Measured 2026-11-10 (`mitchellmj_41603124`): a glob copy out of the scratchpad
# put **17 figure files belonging to OTHER PAPERS** into that curation's
# `source/figs/`. It was caught only because the agent listed the directory
# immediately after the copy and re-hashed the survivors against its own
# download record (6/6 identical). **A file-count or a size check would have
# passed**, and a digitiser pointed at a foreign figure produces values that are
# internally consistent, geometrically valid and attributed to the wrong paper.
#
# That is the worst failure shape this project has: not a wrong number, but a
# right number on the wrong molecule, which is the `AC_2025` class (419
# measurements, F-014).
#
# RULES:
#   * never `cp <scratchpad>/*` into a curation folder -- copy named files only;
#   * after any copy into `source/`, list the directory and assert every file is
#     one you fetched, by name AND by hash against your own download record;
#   * prefer fetching straight into `<folder>/source/` so the scratchpad is
#     never in the path at all.
# ---------------------------------------------------------------------------

def main():
    raise SystemExit(
        "stage07_emit.py is still the generated stub. Write it, or delete it and "
        "record in DECISIONS.md why this stage does not apply to this paper. "
        "An empty stage that silently returns is how three ZT_2026 marginal "
        "checks scored a literal 0 for twelve days while the summary beside "
        "them said ALL FOUR EXACT (F-038).")


if __name__ == "__main__":
    main()
