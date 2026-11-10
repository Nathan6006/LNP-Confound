#!/usr/bin/env python3
"""Stage 03 -- enumerate the 18 distinct ionizable structures, and check each
one against the SI's own [M+H]+ molecular formula.

    python3 -I stage03_enumerate.py

THE RULE, AND WHERE IT COMES FROM. The SI (p. S4) says amine heads were mixed
with acrylate tails at 1/2.5 molar ratio for 48 h at 70 C, and calls the
products "fully substituted". So EVERY AMINE N-H IS CONSUMED: a primary amine
takes two tails, a secondary amine one. The single-tailed series was made at
1/1 and takes exactly one.

That rule is not asserted, it is CHECKED. The SI's formula tables fix the tail
count independently of it:
    75-OCholB   C76H131N4O8S4  -> S4  -> 2 OCholB tails (2 S each)
    304-OCholB  C114H196N7O12S6 -> S6 -> 3 tails, and 304 has exactly 3 N-H.
A formula still cannot see connectivity or regiochemistry (error class 14), so
the C10 anchor in make_reference.py is what tests topology.

CARBAMATE NITROGEN IS NOT AN AMINE. The OCholB and OChol tails each carry an
N-H inside a carbamate. A naive "react every N-H" loop would alkylate it. The
SMARTS below excludes any N bonded to a carbonyl; the formula check is what
would catch the mistake (an extra tail changes C, H, O and S at once), but the
exclusion is written where a reviewer reads it rather than left to the check.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
# python3 -I drops the script's own directory from sys.path.
sys.path.insert(0, HERE)

from rdkit import Chem                                          # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula          # noqa: E402

import assemble                                                 # noqa: E402
from stage01_blocks import HEADS, build_tails                   # noqa: E402

# An amine N-H. Excludes amide/carbamate N (N bonded to C=O) and sulfonamide.
AMINE_NH = Chem.MolFromSmarts("[NX3;H1,H2;!$(N[CX3]=[OX1]);!$(N[SX4]=[OX1])]")

# The SI's own tables. Transcribed by hand from the PDF text layer, which is
# why they are an EXTERNAL check and not a restatement of the enumeration.
SI_FORMULA = {
    # Figure S13, p. S25
    "75-OCholB": "C76H131N4O8S4", "76-OCholB": "C76H129N4O8S4",
    "77-OCholB": "C77H131N4O8S4",
    "75-O16B": "C40H81N2O4S4", "76-O16B": "C40H79N2O4S4",
    "77-O16B": "C41H81N2O4S4",
    "75-OChol": "C80H139N4O8", "76-OChol": "C80H137N4O8",
    "77-OChol": "C81H139N4O8",
    # Figure S16, p. S28
    "78-OCholB": "C74H127N4O8S4", "80-OCholB": "C75H129N4O8S4",
    "81-OCholB": "C77H133N4O8S4", "87-OCholB": "C77H133N4O10S4",
    "90-OCholB": "C76H129N4O9S4", "304-OCholB": "C114H196N7O12S6",
    # Figure S22E, p. S34
    "75-OCholB-1": "C41H74N3O4S2", "76-OCholB-1": "C41H72N3O4S2",
    "77-OCholB-1": "C42H74N3O4S2",
}
SI_CAL_MZ = {
    "75-OCholB": 1357.13, "76-OCholB": 1355.12, "77-OCholB": 1369.15,
    "75-O16B": 781.51, "76-O16B": 779.49, "77-O16B": 793.51,
    "75-OChol": 1284.06, "76-OChol": 1282.04, "77-OChol": 1296.06,
    "78-OCholB": 1329.08, "80-OCholB": 1343.11, "81-OCholB": 1371.16,
    "87-OCholB": 1403.16, "90-OCholB": 1371.12, "304-OCholB": 2049.21,
    "75-OCholB-1": 736.51, "76-OCholB-1": 734.50, "77-OCholB-1": 748.51,
}
SI_OBS_MZ = {
    "75-OCholB": 1357.17, "76-OCholB": 1356.84, "77-OCholB": 1369.87,
    "75-O16B": 781.33, "76-O16B": 779.73, "77-O16B": 793.63,
    "75-OChol": 1284.09, "76-OChol": 1282.36, "77-OChol": 1295.91,
    "78-OCholB": 1328.50, "80-OCholB": 1342.96, "81-OCholB": 1370.97,
    "87-OCholB": 1403.01, "90-OCholB": 1370.86, "304-OCholB": 2049.22,
    "75-OCholB-1": 736.55, "76-OCholB-1": 734.64, "77-OCholB-1": 748.73,
}
AVG = {"C": 12.011, "H": 1.008, "N": 14.007, "O": 15.999, "S": 32.06}


def avg_mass_plus_h(formula_mh):
    """Average mass of a [M+H]+ formula string. The SI's 'Cal.' column is an
    AVERAGE mass, not monoisotopic -- established by reproducing 1357.13 for
    C76H131N4O8S4 (monoisotopic would be 1355.885)."""
    import re
    tot = 0.0
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula_mh):
        if not el:
            continue
        tot += AVG[el] * (int(n) if n else 1)
    return tot


def plus_h(formula):
    """Neutral formula -> [M+H]+ formula string, in Hill order."""
    import re
    counts = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if el:
            counts[el] = counts.get(el, 0) + (int(n) if n else 1)
    counts["H"] = counts.get("H", 0) + 1
    order = ["C", "H"] + sorted(k for k in counts if k not in ("C", "H"))
    return "".join(f"{e}{counts[e]}" if counts[e] != 1 else e
                   for e in order if e in counts)


def _tag_one_amine(smiles, mapno=1):
    """Map ONE amine N-H, chosen by canonical rank, not by atom index.

    'Never bond to atom index 0 of the fragment' (CLAUDE.md s6) -- canonical
    ranking is a property of the molecule, an index is a property of how the
    string happened to be written.
    """
    m = Chem.MolFromSmiles(smiles)
    hits = [h[0] for h in m.GetSubstructMatches(AMINE_NH)]
    if not hits:
        return None
    ranks = list(Chem.CanonicalRankAtoms(m, breakTies=True))
    idx = min(hits, key=lambda i: ranks[i])
    m.GetAtomWithIdx(idx).SetAtomMapNum(mapno)
    return Chem.MolToSmiles(m)


def michael(head_smiles, handle, n_tails=None, label="", max_tails=6):
    """Aza-Michael additions until every amine N-H is consumed (n_tails=None),
    or exactly n_tails of them.

    `max_tails` is a RUNAWAY BOUND, not a chemical statement. If the carbamate
    exclusion in AMINE_NH is ever lost, each OCholB tail added brings its own
    carbamate N-H and the loop never terminates -- the failure would be a hang,
    which reads as an environment problem rather than as a defect. Stopping at
    `max_tails` turns it into a formula mismatch, which is a number that moves.
    Every real member here takes 1, 2 or 3 tails; main() asserts that."""
    cur = head_smiles
    made = 0
    while True:
        if n_tails is not None and made >= n_tails:
            break
        if made >= max_tails:
            break
        tagged = _tag_one_amine(cur)
        if tagged is None:
            break
        cur = assemble.attach(tagged, handle, a_map=1, b_map=2,
                              label=f"{label}#{made + 1}")
        made += 1
    return cur, made


def main():
    tails, sterol = build_tails()
    heads = {num: smi for num, smi, _f, _n in HEADS}
    members = []
    # fully substituted
    for num in [h[0] for h in HEADS]:
        for tail in ("OCholB",) if num not in ("75", "76", "77") else ("OCholB", "O16B", "OChol"):
            name = f"{num}-{tail}"
            smi, n = michael(heads[num], tails[tail]["handle"], label=name)
            members.append((name, num, tail, n, smi, "full"))
    # single-tailed, 1/1 molar ratio (SI p. S4; Figure S22)
    for num in ("75", "76", "77"):
        name = f"{num}-OCholB-1"
        smi, n = michael(heads[num], tails["OCholB"]["handle"], n_tails=1,
                         label=name)
        members.append((name, num, "OCholB", n, smi, "single_tail"))

    rows, bad = [], []
    counts = sorted({n for _nm, _h, _t, n, _s, _x in members})
    if not set(counts) <= {1, 2, 3}:
        raise SystemExit(f"tail counts {counts} outside the expected {{1,2,3}}; "
                         "the SI's S4/S6 sulfur counts fix 2 tails for the "
                         "eight primary-amine heads and 3 for 304.")
    for name, head, tail, n, smi, series in members:
        m = Chem.MolFromSmiles(smi)
        f = CalcMolFormula(m)
        mh = plus_h(f)
        want = SI_FORMULA[name]
        ok = (mh == want)
        if not ok:
            bad.append((name, mh, want))
        rows.append({
            "IL_name": name, "head": head, "tail": tail, "n_tails": n,
            "series": series,
            "IL_SMILES": smi,
            "formula_neutral": f,
            "formula_MH_computed": mh,
            "formula_MH_SI": want,
            "formula_agrees": int(ok),
            "cal_mz_SI": SI_CAL_MZ[name],
            "cal_mz_recomputed": round(avg_mass_plus_h(want), 2),
            "obs_mz_SI": SI_OBS_MZ[name],
            "inchikey": Chem.MolToInchiKey(m),
            "skeleton": Chem.MolToInchiKey(m).split("-")[0],
        })
    out = os.path.join(HERE, "03_enumerated.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}: {len(rows)} distinct ionizable structures")
    print(f"  formula agrees with SI: {sum(r['formula_agrees'] for r in rows)}/{len(rows)}")
    if bad:
        print("  MISMATCHES:", bad)
        return 1
    keys = {r["inchikey"] for r in rows}
    print(f"  distinct full InChIKeys: {len(keys)}  (rows {len(rows)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
