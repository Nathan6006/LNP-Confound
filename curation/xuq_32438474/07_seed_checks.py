#!/usr/bin/env python3
"""Seed each structural check on purpose and require the AGREEMENT COUNT to MOVE.

    python3 -I 07_seed_checks.py

CLAUDE.md s6: "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON
PURPOSE." Three curations in this run shipped seeds that were no-ops -- one
shadowed by a default argument, one that swapped two identical monomers, one
that merely raised an exception somewhere unrelated. So every seed here does
two things:

  1. asserts the mutated input is actually DIFFERENT from the real one, and
  2. reports the formula-agreement count, which must move off 18/18.

Seed C is deliberately a seed that the formula check CANNOT catch. It is here
to make that blindness a measured number rather than a sentence.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)

from rdkit import Chem                                          # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula          # noqa: E402

import stage01_blocks as S1                                     # noqa: E402
import stage03_enumerate as S3                                  # noqa: E402

EPICHOLESTEROL = ("CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@H](O)CC[C@]"
                  "4(C)[C@H]3CC[C@]12C")
PIPERAZINE_304 = "CNCCN1CCN(CCNC)CC1"


def enumerate_all(heads=None, amine_smarts=None, sterol=None):
    """Re-run the enumeration with one thing swapped. Returns
    (agreeing, total, rows)."""
    saved_smarts = S3.AMINE_NH
    saved_reg = S1.reagents.registered
    if amine_smarts is not None:
        S3.AMINE_NH = Chem.MolFromSmarts(amine_smarts)
    if sterol is not None:
        S1.reagents.registered = lambda name: (sterol if name == "Cholesterol"
                                               else saved_reg(name))
    try:
        tails, _ = S1.build_tails()
        h = dict(heads or {n: s for n, s, _f, _x in S1.HEADS})
        agree = 0
        rows = []
        names = [n for n, _s, _f, _x in S1.HEADS]
        for num in names:
            tl = ("OCholB", "O16B", "OChol") if num in ("75", "76", "77") else ("OCholB",)
            for t in tl:
                nm = f"{num}-{t}"
                smi, _n = S3.michael(h[num], tails[t]["handle"], label=nm)
                rows.append((nm, smi))
        for num in ("75", "76", "77"):
            nm = f"{num}-OCholB-1"
            smi, _n = S3.michael(h[num], tails["OCholB"]["handle"], n_tails=1,
                                 label=nm)
            rows.append((nm, smi))
        out = []
        for nm, smi in rows:
            m = Chem.MolFromSmiles(smi)
            mh = S3.plus_h(CalcMolFormula(m))
            ok = mh == S3.SI_FORMULA[nm]
            agree += ok
            out.append((nm, mh, S3.SI_FORMULA[nm], ok,
                        Chem.MolToInchiKey(m)))
        return agree, len(out), out
    finally:
        S3.AMINE_NH = saved_smarts
        S1.reagents.registered = saved_reg


def main():
    res = {}
    base_a, base_n, base_rows = enumerate_all()
    res["control"] = {"agree": base_a, "total": base_n,
                      "note": "unmutated run"}
    base_keys = {r[0]: r[4] for r in base_rows}

    # ---- seed A: amine 304 as the piperazine isomer the drawing suggests ----
    real304 = dict((n, s) for n, s, _f, _x in S1.HEADS)["304"]
    assert Chem.MolToSmiles(Chem.MolFromSmiles(real304)) != \
        Chem.MolToSmiles(Chem.MolFromSmiles(PIPERAZINE_304)), \
        "seed A is a no-op: the two 304 candidates canonicalise identically."
    heads = {n: s for n, s, _f, _x in S1.HEADS}
    heads["304"] = PIPERAZINE_304
    a, n, rows = enumerate_all(heads=heads)
    res["seedA_304_as_piperazine"] = {
        "agree": a, "total": n, "moved": a != base_a,
        "what_changed": [r[0] for r in rows if not r[3]],
        "note": "1,4-bis(2-(methylamino)ethyl)piperazine C10H24N4 instead of "
                "tris(2-(methylamino)ethyl)amine C9H24N4. This is the real "
                "ambiguity in Figure 2A: the drawing's two upper arms close "
                "visually into what looks like a hexagon.",
    }

    # ---- seed B: let the carbamate N-H react ----
    a, n, rows = enumerate_all(amine_smarts="[NX3;H1,H2]")
    res["seedB_carbamate_N_reacts"] = {
        "agree": a, "total": n, "moved": a != base_a,
        "what_changed": [r[0] for r in rows if not r[3]][:6],
        "note": "drops the !$(N[CX3]=[OX1]) exclusion, so every OCholB/OChol "
                "carbamate N-H takes a tail as well.",
    }

    # ---- seed C: epicholesterol. THE FORMULA CHECK CANNOT SEE THIS. ----
    assert EPICHOLESTEROL != S1.reagents.registered("Cholesterol"), \
        "seed C is a no-op: epicholesterol string equals the registry string."
    a, n, rows = enumerate_all(sterol=EPICHOLESTEROL)
    changed_keys = [nm for nm, _mh, _w, _ok, k in rows
                    if base_keys.get(nm) and k != base_keys[nm]]
    same_skeleton = [nm for nm, _mh, _w, _ok, k in rows
                     if base_keys.get(nm)
                     and k.split("-")[0] == base_keys[nm].split("-")[0]]
    res["seedC_epicholesterol"] = {
        "agree": a, "total": n, "moved": a != base_a,
        "full_inchikeys_changed": len(changed_keys),
        "skeleton_inchikeys_unchanged": len(same_skeleton),
        "note": "F-013's exact mechanism, as a measured blindness. The "
                "formula check does NOT move (18/18 still agree) and the "
                "InChIKey SKELETON block does not change on any affected "
                "structure. Only the FULL InChIKey moves. This is why "
                "03_enumerated.csv stores the full key and why the C10 anchor "
                "in make_reference.py is required.",
    }

    verdict = {
        "seedA_moved": res["seedA_304_as_piperazine"]["moved"],
        "seedB_moved": res["seedB_carbamate_N_reacts"]["moved"],
        "seedC_moved_formula": res["seedC_epicholesterol"]["moved"],
        "seedC_moved_full_inchikey":
            res["seedC_epicholesterol"]["full_inchikeys_changed"] > 0,
    }
    res["_verdict"] = verdict
    with open(os.path.join(HERE, "07_seed_checks.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))
    ok = (verdict["seedA_moved"] and verdict["seedB_moved"]
          and not verdict["seedC_moved_formula"]
          and verdict["seedC_moved_full_inchikey"])
    print("SEEDS OK" if ok else "SEEDS NOT AS EXPECTED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
