#!/usr/bin/env python3
"""Stage 03 -- enumerate what is enumerable, and record what is not.

ENUMERATED: the two tail reagents, built through curation/lib/assemble with
atom-mapped attachment points. No string concatenation anywhere.

NOT ENUMERATED: all 16 lipidoids. A lipidoid here is (amine head + n x oleyl
acrylamide) and the SI names, draws and formulates no amine. See
01_blocks_qc.json -> amine_heads.

STEREO POSTURE, STATED SO A CLEAN LINE IS NOT MISREAD AS EVIDENCE.
Neither product has a single stereogenic CARBON. assemble.check_stereo_
preserved compares CIP code counts, so on these molecules it compares 0 to 0:
IT IS VACUOUS, NOT PASSING. The one piece of stereochemistry that does exist
is the oleyl C9-C10 double bond's Z geometry, which CIP-centre counting cannot
see at all, so this stage asserts it explicitly and 07_qc.json seeds its loss.

Run:  python3 -I stage03_enumerate.py
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import assemble                                           # noqa: E402
from rdkit import Chem, RDLogger                           # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt              # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula     # noqa: E402

RDLogger.DisableLog("rdApp.*")

# --- fragments, each with an atom-mapped attachment atom -------------------
# The acryloyl group's NITROGEN is the attachment atom (map 1): this is an
# acrylAMIDE. Map 2 is C1 of the alkyl chain -- the end that bonds to N, which
# is what fixes the chain's DIRECTION and therefore puts the alkene at C9
# counting FROM NITROGEN.
ACRYLAMIDE_N = "C=CC(=O)[NH:1]"
TAILS = {
    # (9Z)-octadec-9-en-1-yl: mapped CH2 = C1, then 7 CH2 = C2..C8, then the
    # cis alkene C9=C10, then 8 carbons C11..C18.
    "oleyl_acrylamide": "[CH2:2]CCCCCCC/C=C\\CCCCCCCC",
    # octadecyl: mapped CH2 = C1, then 17 more carbons = C2..C18.
    "octadecyl_acrylamide": "[CH2:2]CCCCCCCCCCCCCCCCC",
}

EXPECTED = {
    # Independent of the build: computed from the chemical name's atom count by
    # hand and written here, so a build that silently loses or gains an atom
    # fails rather than being believed.
    "oleyl_acrylamide": {"formula": "C21H39NO", "chain_C": 18,
                         "n_C_C_double_bonds_in_chain": 1},
    "octadecyl_acrylamide": {"formula": "C21H41NO", "chain_C": 18,
                             "n_C_C_double_bonds_in_chain": 0},
}


def chain_carbons(frag_smiles):
    """Carbon count of the alkyl fragment, measured from the fragment."""
    m = Chem.MolFromSmiles(frag_smiles)
    if m is None:
        raise ValueError(f"unparseable fragment {frag_smiles!r}")
    return sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C")


def alkene_position_from_nitrogen(product_smiles):
    """Which chain carbon, counting from N, carries the C=C. Returns None if
    the chain has no C=C.

    This is the quantity mass CANNOT see. Reversing the chain past the alkene
    gives an identical formula and an identical exact mass; only the graph
    distance from nitrogen separates them.
    """
    m = Chem.MolFromSmiles(product_smiles)
    n_idx = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N"]
    if len(n_idx) != 1:
        raise ValueError("expected exactly one nitrogen")
    n = n_idx[0]
    # the chain is the N-substituent that is not the carbonyl carbon
    nbrs = [a for a in m.GetAtomWithIdx(n).GetNeighbors()
            if a.GetSymbol() == "C"
            and not any(b.GetBondType() == Chem.BondType.DOUBLE
                        and b.GetOtherAtom(a).GetSymbol() == "O"
                        for b in a.GetBonds())]
    if len(nbrs) != 1:
        raise ValueError("expected exactly one alkyl substituent on N")
    # walk the chain
    prev, cur, pos = n, nbrs[0].GetIdx(), 1
    while True:
        at = m.GetAtomWithIdx(cur)
        for b in at.GetBonds():
            if b.GetBondType() == Chem.BondType.DOUBLE and \
                    b.GetOtherAtom(at).GetSymbol() == "C":
                other = b.GetOtherAtom(at).GetIdx()
                if other != prev:
                    return pos      # cur is the lower-numbered alkene carbon
        nxt = [a.GetIdx() for a in at.GetNeighbors()
               if a.GetIdx() != prev and a.GetSymbol() == "C"]
        if not nxt:
            return None
        if len(nxt) > 1:
            raise ValueError("chain branches; not a linear alkyl")
        prev, cur, pos = cur, nxt[0], pos + 1


def double_bond_stereo(product_smiles):
    """The chain alkene's assigned geometry, or None."""
    m = Chem.MolFromSmiles(product_smiles)
    Chem.AssignStereochemistry(m, cleanIt=True, force=True)
    out = []
    for b in m.GetBonds():
        if b.GetBondType() != Chem.BondType.DOUBLE:
            continue
        if b.GetBeginAtom().GetSymbol() != "C" or \
                b.GetEndAtom().GetSymbol() != "C":
            continue
        # skip the acryloyl vinyl -- it is terminal, so it has no E/Z
        if b.GetBeginAtom().GetTotalNumHs() == 2 or \
                b.GetEndAtom().GetTotalNumHs() == 2:
            continue
        out.append(str(b.GetStereo()))
    return out[0] if len(out) == 1 else (None if not out else out)


def build(name, tail_smiles):
    prod = assemble.attach(ACRYLAMIDE_N, tail_smiles,
                           a_map=1, b_map=2, bond="single", label=name)
    m = Chem.MolFromSmiles(prod)
    if m is None:
        raise ValueError(f"{name}: product does not parse")
    assemble.check_no_radicals(m, label=name)
    # VACUOUS on these molecules -- recorded so, not credited.
    cip = assemble.check_stereo_preserved(
        prod, [ACRYLAMIDE_N, tail_smiles], label=name)
    return {
        "name": name,
        "SMILES": prod,
        "formula": CalcMolFormula(m),
        "exact_mass": round(ExactMolWt(m), 4),
        "inchikey": Chem.MolToInchiKey(m),
        "chain_C": chain_carbons(tail_smiles),
        "alkene_pos_from_N": alkene_position_from_nitrogen(prod),
        "alkene_stereo": double_bond_stereo(prod),
        "n_cip_centres": len(assemble.stereocentres(m)),
        "cip_check": cip,
    }


# ---------------------------------------------------------------------------
# Seeds. Each is applied in its own process by stage07_qc.py, which reads the
# exit code from returncode. The important one is seed_chain_direction: it has
# an IDENTICAL molecular formula and an IDENTICAL exact mass, so it is a
# defect only a connectivity test can see.
# ---------------------------------------------------------------------------
SEEDS = {
    # oleyl written from the METHYL end: alkene lands at C10 from nitrogen.
    # C21H39NO, 321.3032 -- identical on both counts.
    "seed_chain_direction": ("oleyl_acrylamide",
                             "[CH2:2]CCCCCCCC/C=C\\CCCCCCC"),
    # the Z geometry dropped. Formula and exact mass again identical.
    "seed_alkene_stereo": ("oleyl_acrylamide",
                           "[CH2:2]CCCCCCCC=CCCCCCCCC"),
    # one carbon short: this one the FORMULA can see, which is the point of
    # having it -- it separates "mass sees it" from "only the graph sees it".
    "seed_chain_length": ("oleyl_acrylamide",
                          "[CH2:2]CCCCCC/C=C\\CCCCCCCC"),
}


def main(seed=None):
    built, fails = {}, []
    tails = dict(TAILS)
    if seed is not None:
        if seed not in SEEDS:
            raise SystemExit(f"unknown seed {seed!r}")
        k, v = SEEDS[seed]
        tails[k] = v
    for name, tail in tails.items():
        b = build(name, tail)
        built[name] = b
        exp = EXPECTED[name]
        if b["formula"] != exp["formula"]:
            fails.append(f"{name}: formula {b['formula']} != {exp['formula']}")
        if b["chain_C"] != exp["chain_C"]:
            fails.append(f"{name}: chain C {b['chain_C']} != {exp['chain_C']}")
        nd = 0 if b["alkene_pos_from_N"] is None else 1
        if nd != exp["n_C_C_double_bonds_in_chain"]:
            fails.append(f"{name}: chain C=C count {nd} != "
                         f"{exp['n_C_C_double_bonds_in_chain']}")
        if name == "oleyl_acrylamide":
            if b["alkene_pos_from_N"] != 9:
                fails.append(f"{name}: alkene at C{b['alkene_pos_from_N']} "
                             "from N, the name says C9")
            if b["alkene_stereo"] != "STEREOZ":
                fails.append(f"{name}: alkene stereo {b['alkene_stereo']!r}, "
                             "the name 'oleyl' says Z (cis)")
        if b["n_cip_centres"] != 0:
            fails.append(f"{name}: {b['n_cip_centres']} CIP centres, expected 0")

    # ---- write 03_enumerated.csv -----------------------------------------
    # The corpus column names are IL_SMILES / IL_name. ZERO LIPIDOIDS are
    # enumerable, so no row carries an IL_SMILES. The two tail reagents are
    # written with their own columns and IL_SMILES deliberately BLANK, so a
    # reader cannot mistake a reagent for a library member.
    cols = ["IL_name", "IL_SMILES", "role", "tail_SMILES", "formula",
            "exact_mass", "inchikey", "chain_C", "alkene_pos_from_N",
            "alkene_stereo", "n_cip_centres", "confidence", "note"]
    rows = []
    conf = {"oleyl_acrylamide": "high", "octadecyl_acrylamide": "medium"}
    notes = {
        "oleyl_acrylamide":
            "Tail REAGENT, not a library member. Structure derived from the "
            "chemical name printed in the abstract ('oleyl acrylamide'); the "
            "paper draws nothing. IL_SMILES is blank because the amine head "
            "this reagent was added to is not disclosed anywhere in the SI.",
        "octadecyl_acrylamide":
            "Tail REAGENT, not a library member. ASSUMED FEATURE: the "
            "acrylamide linkage, carried over by parallel construction -- the "
            "abstract names only 'n-octadecyl tails' and never prints "
            "'octadecyl acrylamide'. IL_SMILES blank for the same reason as "
            "the oleyl reagent.",
    }
    for name, b in built.items():
        rows.append({
            "IL_name": "", "IL_SMILES": "",
            "role": f"tail_reagent:{name}",
            "tail_SMILES": b["SMILES"], "formula": b["formula"],
            "exact_mass": b["exact_mass"], "inchikey": b["inchikey"],
            "chain_C": b["chain_C"],
            "alkene_pos_from_N": ("" if b["alkene_pos_from_N"] is None
                                  else b["alkene_pos_from_N"]),
            "alkene_stereo": b["alkene_stereo"] or "",
            "n_cip_centres": b["n_cip_centres"],
            "confidence": conf[name], "note": notes[name],
        })
    if seed is not None:
        for name, b in built.items():
            print(f"[seed={seed}] {name}: {b['formula']} {b['exact_mass']} "
                  f"chainC={b['chain_C']} alkene@C{b['alkene_pos_from_N']} "
                  f"{b['alkene_stereo']}")
        print(f"n_agree={3 - len(fails)} n_total=3")
        if fails:
            print("FAILED:", *fails, sep="\n  - ")
            return 1
        return 0
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    qc = {
        "lipidoids_enumerated": 0,
        "lipidoids_that_exist": 16,
        "why_zero": "the amine heads are not disclosed; see 01_blocks_qc.json",
        "tail_reagents_enumerated": len(rows),
        "builds": {k: {kk: vv for kk, vv in v.items() if kk != "cip_check"}
                   for k, v in built.items()},
        "stereo_posture": {
            "cip_centres_in_either_product": 0,
            "assemble_check_stereo_preserved": "RUN AND PASSED, BUT VACUOUS. "
                "It compares CIP code counts and both sides are 0. A clean "
                "line here is evidence of nothing; it is recorded only so the "
                "absence is explicit rather than silent.",
            "the_stereochemistry_that_does_exist": "the oleyl C9=C10 Z "
                "(cis) double bond, which CIP-centre counting cannot see. "
                "Asserted explicitly by double_bond_stereo(); seeded in "
                "07_qc.json (seed_alkene_stereo).",
            "assemble_parity_preserved": "NOT CALLED. It compares spatial "
                "parity of assigned stereoCENTRES, of which there are none.",
        },
        "direction_posture": {
            "what_mass_cannot_see": "the chain's direction. Writing the oleyl "
                "chain from the methyl end puts the alkene at C10 from "
                "nitrogen with formula C21H39NO and exact mass 321.3032 -- "
                "IDENTICAL on both counts. alkene_pos_from_N() measures the "
                "graph distance from nitrogen, which is the only quantity "
                "that separates them; seeded in 07_qc.json "
                "(seed_chain_direction).",
        },
        "checks_failed": fails,
    }
    with open(os.path.join(HERE, "03_enumerated_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=2)

    for name, b in built.items():
        print(f"{name}: {b['SMILES']}")
        print(f"   {b['formula']} {b['exact_mass']} {b['inchikey']} "
              f"chainC={b['chain_C']} alkene@C{b['alkene_pos_from_N']}-from-N "
              f"{b['alkene_stereo']} cip={b['n_cip_centres']}")
    print(f"lipidoids enumerated: 0 of 16 (amine heads not disclosed)")
    print(f"n_agree={3 - len(fails)} n_total=3")
    if fails:
        print("FAILED:")
        for f in fails:
            print("  -", f)
        return 1
    return 0


if __name__ == "__main__":
    _s = None
    if "--seed" in sys.argv:
        _s = sys.argv[sys.argv.index("--seed") + 1]
    sys.exit(main(seed=_s))
