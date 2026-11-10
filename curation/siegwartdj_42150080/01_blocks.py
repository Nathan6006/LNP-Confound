#!/usr/bin/env python3
"""Stage 01 - establish the library size. Run: python3 -I 01_blocks.py

The task brief said "Library size not recorded - establishing it is your first
job", and the answer is trap 1 from the brief: the names on the Fig S2/S3 axes
are FORMULATIONS, not ionizable lipids.

The names are read SYMBOLICALLY out of the PDF's own text layer (page 5 carries
Fig S2A, S2B and S3, and each repeats the same category axis), never typed by
hand and never read off pixels. The script asserts that all three axes carry the
same 16 names in the same order; if they ever stop agreeing it raises instead of
emitting a plausible list.
"""
import csv
import os
import re
import sys

import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "pnas.2528477123.sapp.pdf")
PAGE = 4  # 0-based; SI p.5, which carries Fig S2A, S2B and S3

# A formulation name is "Base LNP" or "18:1 <HELPER> <molpct>".
NAME_RE = re.compile(r"^(Base LNP|18:1 (TAP|DAP|PA) (10|20|30|40|50))$")

# Avanti catalogue identities of the three helper lipids. These are NOT read
# from the paper (the SI never spells them out); they are the catalogue meaning
# of the shorthand, recorded so a reviewer can object to the mapping itself.
HELPER = {
    "TAP": ("DOTAP", "1,2-dioleoyl-3-trimethylammonium-propane",
            "permanently cationic quaternary ammonium"),
    "DAP": ("DODAP", "1,2-dioleoyl-3-dimethylammonium-propane",
            "ionizable tertiary amine"),
    "PA":  ("DOPA",  "1,2-dioleoyl-sn-glycero-3-phosphate",
            "anionic phosphate"),
}


def axis_runs(lines):
    """Every maximal run of consecutive lines that are all formulation names."""
    runs, cur = [], []
    for ln in lines:
        if NAME_RE.match(ln):
            cur.append(ln)
        else:
            if len(cur) > 1:
                runs.append(cur)
            cur = []
    if len(cur) > 1:
        runs.append(cur)
    return runs


def main():
    doc = pymupdf.open(PDF)
    lines = [ln.strip() for ln in doc[PAGE].get_text().splitlines()]
    runs = axis_runs(lines)

    if len(runs) != 3:
        sys.exit(f"01_blocks: expected 3 category axes on SI p.{PAGE+1} "
                 f"(Fig S2A, S2B, S3), found {len(runs)}. Refusing to emit.")
    if runs[0] != runs[1] or runs[1] != runs[2]:
        sys.exit("01_blocks: the three category axes do NOT carry identical "
                 "name lists. Refusing to emit a merged list.")

    names = runs[0]
    if len(names) != len(set(names)):
        sys.exit("01_blocks: duplicate name within one axis. Refusing.")

    rows = []
    for n in names:
        if n == "Base LNP":
            rows.append(dict(
                formulation=n, helper_lipid_shorthand="", helper_common_name="",
                helper_systematic_name="", helper_charge_class="",
                helper_mol_pct="", ionizable_lipid="6A1-SC8 (inferred, see note)",
                n_distinct_ionizable_structures=1,
                source="SI Fig S2A/S2B/S3 category axis, PDF text layer",
                note="the formulation every other row is built on; its full "
                     "composition is defined in the main text, which is not held"))
        else:
            _, helper, pct = NAME_RE.match(n).groups()
            common, systematic, charge = HELPER[helper]
            rows.append(dict(
                formulation=n, helper_lipid_shorthand="18:1 " + helper,
                helper_common_name=common, helper_systematic_name=systematic,
                helper_charge_class=charge, helper_mol_pct=pct,
                ionizable_lipid="6A1-SC8 (inferred, see note)",
                n_distinct_ionizable_structures=1,
                source="SI Fig S2A/S2B/S3 category axis, PDF text layer",
                note="Base LNP + " + common + " at " + pct + " mol%; the "
                     "ionizable lipid does not change across this series"))

    out = os.path.join(HERE, "01_blocks.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"formulations            {len(rows)}")
    print(f"helper lipids distinct  {len(set(r['helper_common_name'] for r in rows) - {''})}")
    print(f"helper mol%% levels      {sorted(set(r['helper_mol_pct'] for r in rows) - {''}, key=int)}")
    print(f"distinct ionizable structures visible in the SI   1")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
