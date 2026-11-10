#!/usr/bin/env python3
"""The C10 comparison, run locally.

provenance/structure_check.py, check_anchor_independence.py and
check_anchor_coverage.py write tracked files other agents read, so NONE of them
is run. This is the local equivalent; its output is 04_c10_anchor_check.json.

THE GATE IS PROVED TO TEST CONNECTIVITY, not formula and not mass: each anchor
is also compared against a seeded defect with IDENTICAL molecular formula and
IDENTICAL exact mass, and the comparison must fail on the seed while passing on
the real build.
"""
import json, os, sys
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, Descriptors
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..")
sys.path.insert(0, os.path.join(HERE, ".."))
import make_reference as mr                                  # noqa: E402

OUT = os.environ.get("STAGE05_OUT",
                     os.path.join(HERE, "..", "04_c10_anchor_check.json"))


def info(smi):
    m = Chem.MolFromSmiles(smi)
    k = Chem.MolToInchiKey(m)
    return {"smiles": Chem.MolToSmiles(m),
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "exact_mass": round(Descriptors.ExactMolWt(m), 6),
            "inchikey": k, "skeleton": k.split("-")[0]}


# --- seeds: same formula, same exact mass, different connectivity ---------
def seed_mc3_regioisomer():
    """Move the ester-bearing carbon from C19 to C18 of the same C37
    tetraene: the arms become 17 and 19 carbons instead of 18 and 18. Same
    37 chain carbons, same four cis olefins, same C43H79NO2, same
    641.611081 Da -- only the attachment point moves. This is the exact shape
    of the MC3 defect this project has used before as a mass-blind seed."""
    def arm(a, b):
        return "C" * a + "/C=C\\C/C=C\\" + "C" * b + "C"
    head = "CN(C)" + "C" * 3 + "C(=O)O"
    return head + "C(" + arm(7, 4) + ")" + arm(9, 4)


def seed_peg_reversed_ester():
    """Reverse one myristoyl ester: C(=O)O -> OC(=O) on the sn-2 arm."""
    myr = "C(=O)" + "C" * 13
    rev = "C" * 13 + "C(=O)"            # same atoms, opposite direction
    eo = "OCC"
    return "C" + eo * 45 + "OC" + "C(O" + myr + ")" + "CO" + rev


def main():
    rows, qc = mr.build()
    built = {r["lipid_name"]: r["reference_SMILES"]
             for r in rows if r["reference_SMILES"]}

    vit = pd.read_csv(os.path.join(ROOT, "data", "data_vitro.csv"),
                      low_memory=False)
    z = vit[vit.Experiment_ID == "YZ_2022"]
    corpus = {"DLin-MC3-DMA": z.IL_SMILES.dropna().unique()[0],
              "DMG-PEG2000": z.PEG_SMILES.dropna().unique()[0],
              "Cholesterol": z.CHL_SMILES.dropna().unique()[0]}

    results = []
    for name, smi in built.items():
        a, b = info(smi), info(corpus[name])
        results.append({
            "lipid_name": name, "anchor": a, "corpus": b,
            "full_inchikey_match": a["inchikey"] == b["inchikey"],
            "skeleton_match": a["skeleton"] == b["skeleton"],
            "byte_identical_canonical": a["smiles"] == b["smiles"],
            "_note": "Byte identity is EXPECTED and is NOT the signal (F-037): "
                     "a correct re-derivation is byte-identical to a correct "
                     "original once canonicalised. Independence is established "
                     "by make_reference.py's two run-time guards, both of "
                     "which were made to fail on purpose."})

    seeds = []
    for name, fn in (("DLin-MC3-DMA", seed_mc3_regioisomer),
                     ("DMG-PEG2000", seed_peg_reversed_ester)):
        real, bad = info(built[name]), info(fn())
        seeds.append({
            "anchor": name, "defect_class":
            "regiochemistry" if "mc3" in fn.__name__ else "ester direction",
            "formula_identical": real["formula"] == bad["formula"],
            "exact_mass_identical": real["exact_mass"] == bad["exact_mass"],
            "delta_mass_Da": round(bad["exact_mass"] - real["exact_mass"], 6),
            "clean_inchikey": real["inchikey"], "seeded_inchikey":
            bad["inchikey"],
            "skeleton_block_also_differs":
            real["skeleton"] != bad["skeleton"],
            "credited": (real["formula"] == bad["formula"]
                         and real["exact_mass"] == bad["exact_mass"]
                         and real["inchikey"] != bad["inchikey"])})

    out = {"_what": "C10 for maohq_38583366, run locally.",
           "anchors_built": len(built),
           "anchors_empty_on_purpose": [e[0] for e in mr.EMPTY],
           "results": results,
           "exact_matches_on_full_inchikey":
               sum(1 for r in results if r["full_inchikey_match"]),
           "connectivity_seeds": seeds,
           "seeds_credited": sum(1 for s in seeds if s["credited"]),
           "cholesterol_posture":
               "NO ANCHOR. The registry/corpus string is NOT counted as one -- "
               "comparing the corpus to the corpus is comparing a value to "
               "itself. Against an empty anchor the guards are VACUOUS, not "
               "passing. What Fig 1 settles qualitatively is in PROVENANCE.md.",
           "helper_lipid_posture":
               "NO ANCHOR for DOTAP/DDAB/DOPE/DSPC/14PA/18PG: this paper names "
               "them in a box and draws none of them."}
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
    if out["exact_matches_on_full_inchikey"] != len(built) \
            or out["seeds_credited"] != len(seeds):
        raise SystemExit("stage05: C10 FAILED")


if __name__ == "__main__":
    main()
