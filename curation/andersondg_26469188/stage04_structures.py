#!/usr/bin/env python3
"""KK_2015 -- three-way structure concordance, and the arm-count derivation.

TWO JOBS.

1. THREE-WAY HELPER CHECK. For each measured phospholipid, compare
   (a) make_reference.py's independent head-first build from the SI's own
       systematic name,
   (b) provenance/reagent_registry.py's structure written from the chemical name,
   (c) the live corpus's HL_SMILES.
   Compared by FULL InChIKey. This is the only gate in this curation that can
   see a double-bond POSITION error: make_reference's own assertions cannot
   (see its --seed-regioisomer sensitivity note), so --seed-delta6 below seeds
   exactly that defect against THIS check.

2. ARM-COUNT DERIVATION for the amine-200 token claim. Within a homologous
   series that holds the amine fixed and walks the tail length, the carbon
   count rises by (number of arms) per unit of tail length. So the series
   itself COUNTS THE ARMS, from formulas alone, with no structure reading.
   This is what makes the F-061 statement quantitative instead of an opinion.
"""
import argparse
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)   # python3 -I does NOT prepend the script dir
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "provenance"))

from rdkit import Chem, RDLogger                                  # noqa: E402
from rdkit.Chem import rdMolDescriptors as rdm                    # noqa: E402

import make_reference                                             # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUTDIR = os.environ.get("KK2015_OUTDIR", HERE)
DATA = os.path.join(HERE, "..", "..", "data")

MEASURED_HELPERS = ("DSPC", "DOPE", "DOPC")


def registry_entry(name):
    """Read the registry's SMILES WITHOUT importing it into make_reference's
    process -- this module is allowed to, make_reference.py is not."""
    import reagent_registry
    return reagent_registry.REGISTRY[name]


def corpus_helper(name):
    smis = set()
    rows = 0
    for f, in (("data_vitro.csv",), ("data_vivo.csv",)):
        with open(os.path.join(DATA, f), newline="") as fh:
            for r in csv.DictReader(fh):
                if (r.get("HL_name") or "").strip() == name and r.get("HL_SMILES"):
                    smis.add(r["HL_SMILES"])
                    rows += 1
    return smis, rows


def ik(s):
    m = Chem.MolFromSmiles(s)
    return Chem.MolToInchiKey(m), rdm.CalcMolFormula(m), rdm.CalcExactMolWt(m)


def three_way(override=None):
    out = {}
    for name in MEASURED_HELPERS:
        built = override[name] if override and name in override \
            else make_reference.build(name)[0]
        reg = registry_entry(name)[0]
        cs, nrows = corpus_helper(name)
        b_ik, b_f, b_m = ik(built)
        r_ik, r_f, r_m = ik(reg)
        if len(cs) != 1:
            raise AssertionError("%s: corpus holds %d distinct SMILES" % (name, len(cs)))
        c_smi = next(iter(cs))
        c_ik, c_f, c_m = ik(c_smi)
        out[name] = {
            "built_inchikey": b_ik, "registry_inchikey": r_ik,
            "corpus_inchikey": c_ik,
            "all_three_agree": b_ik == r_ik == c_ik,
            "built_eq_registry": b_ik == r_ik,
            "built_eq_corpus": b_ik == c_ik,
            "formula": b_f,
            "formula_all_agree": b_f == r_f == c_f,
            "exact_mass": round(b_m, 6),
            "exact_mass_max_delta": round(max(abs(b_m - r_m), abs(b_m - c_m)), 8),
            "corpus_rows": nrows,
            "built_is_byte_identical_to_registry": built == reg,
            "byte_identity_note": "Byte identity is EXPECTED and is NOT the "
                                  "signal (F-037): a correct re-derivation "
                                  "canonicalises to the same string as a correct "
                                  "stored structure. What matters is that the "
                                  "build ROUTE shares nothing with the registry, "
                                  "which make_reference.ast_guard tests at run "
                                  "time on its own AST.",
        }
    return out


def arm_count():
    """Count arms per amine core from formulas alone, across homologous series."""
    rows = []
    for f, pool in (("data_vitro.csv", "vitro"), ("data_vivo.csv", "vivo")):
        with open(os.path.join(DATA, f), newline="") as fh:
            for r in csv.DictReader(fh):
                if r.get("IL_name") and r.get("IL_SMILES"):
                    rows.append((r["Experiment_ID"], r["IL_name"], r["IL_SMILES"]))
    series = {
        "KL_2010 '200C<n>'": ("KL_2010", re.compile(r"^200C(\d+)$")),
        "KW_2014 'A200_C<n>'": ("KW_2014", re.compile(r"^A200_C(\d+)$")),
    }
    out = {}
    for label, (exp, pat) in series.items():
        pts = {}
        for e, n, s in rows:
            if e != exp:
                continue
            m = pat.match(n)
            if m:
                mol = Chem.MolFromSmiles(s)
                pts[int(m.group(1))] = (
                    sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C"),
                    sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N"),
                    rdm.CalcMolFormula(mol), n)
        ks = sorted(pts)
        slopes = []
        for a, b in zip(ks, ks[1:]):
            slopes.append((pts[b][0] - pts[a][0]) / (b - a))
        out[label] = {
            "members": [{"tail_token": k, "name": pts[k][3], "n_C": pts[k][0],
                         "n_N": pts[k][1], "formula": pts[k][2]} for k in ks],
            "dC_per_unit_tail": slopes,
            "slope_constant": len(set(slopes)) == 1,
            "ARMS_IMPLIED": slopes[0] if slopes and len(set(slopes)) == 1 else None,
            "n_nitrogen": sorted({pts[k][1] for k in ks}),
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed-delta6", action="store_true",
                    help="Move the oleoyl double bond from delta-9 to delta-6. "
                         "IDENTICAL formula, IDENTICAL exact mass, every one of "
                         "make_reference's own assertions still passes. Proves "
                         "the three-way InChIKey check is what catches it.")
    args = ap.parse_args()

    if args.seed_delta6:
        good = make_reference.build("DOPC")[0]
        saved = dict(make_reference.ACYL)
        try:
            make_reference.ACYL["O"] = (
                "CCCCCCCCCCC/C=C\\CCCC[C:%d](=O)",
                "(6Z)-octadec-6-enoyl -- SEEDED DEFECT", 18, 1)
            bad, bprov = make_reference.build("DOPC")
        finally:
            make_reference.ACYL.clear()
            make_reference.ACYL.update(saved)
        g_ik, g_f, g_m = ik(good)
        b_ik, b_f, b_m = ik(bad)
        print("SEED delta6: delta-9 oleoyl -> delta-6 oleoyl in DOPC")
        print("  make_reference's own assertions all PASSED on the defect:")
        print("    C2 CIP = %s, C=C = %d, Z = %d, unassigned = %d"
              % (bprov["glycerol_C2_CIP"], bprov["n_CC_double_bonds"],
                 bprov["n_Z"], bprov["n_unassigned_CC"]))
        print("  formula identical : %s  (%s vs %s)" % (g_f == b_f, g_f, b_f))
        print("  exact-mass delta  : %.6f Da" % abs(g_m - b_m))
        tw = three_way(override={"DOPC": bad})
        fired = not tw["DOPC"]["all_three_agree"]
        print("  three_way DOPC    : built=%s registry=%s corpus=%s"
              % (tw["DOPC"]["built_inchikey"], tw["DOPC"]["registry_inchikey"],
                 tw["DOPC"]["corpus_inchikey"]))
        print("  three_way agrees  : %s" % tw["DOPC"]["all_three_agree"])
        ok = (g_f == b_f) and abs(g_m - b_m) < 1e-6 and fired
        print("SEED delta6: %s" % (
            "THREE-WAY CHECK FIRED on a defect with identical formula and "
            "identical exact mass (delta 0.000000 Da) that every internal "
            "assertion passed" if ok else "CHECK STAYED SILENT -- no teeth"))
        return 0 if ok else 1

    out = {"three_way_helper_concordance": three_way(),
           "amine_200_arm_count": arm_count()}
    out["three_way_helper_concordance"]["_independence_posture"] = (
        "This is a NAME-CONCORDANCE check, NOT an external attestation of this "
        "paper's chemistry. The paper draws no structure; it prints a systematic "
        "name, and all three legs of this comparison are derivations from that "
        "same name by different hands. It catches a transcription error in any "
        "leg -- and the briefing records a sibling whose hand-built DSPC came out "
        "a ketone bearing an ether at the right mass to six decimals -- but it "
        "cannot catch a name that is wrong FOR THIS PAPER. Counted as such in "
        "PROVENANCE.md and NOT counted as a C10 pass.")
    out["amine_200_arm_count"]["_what_this_proves"] = (
        "Within a homologous series holding the amine fixed and walking the tail "
        "length, dC/d(tail) EQUALS the number of arms. KW_2014's A200_C10..C14 "
        "gives 5 arms on a 5-nitrogen core; KL_2010's 200C8..C14 gives 1 arm on a "
        "3-nitrogen core. Both series are INTERNALLY self-consistent, which is "
        "precisely why F-061 is a real conflict and not a typo: two curations "
        "each coherently encode a different molecule for the same catalogue "
        "amine. Reported; data/ untouched.")
    with open(os.path.join(OUTDIR, "_structures.json"), "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True, default=str)

    for n, d in out["three_way_helper_concordance"].items():
        if n.startswith("_"):
            continue
        print("%-5s %-12s mass %.6f | 3-way agree=%s | max mass delta %.8f | "
              "corpus rows %d" % (n, d["formula"], d["exact_mass"],
                                  d["all_three_agree"], d["exact_mass_max_delta"],
                                  d["corpus_rows"]))
        print("      %s" % d["built_inchikey"])
    for lab, d in out["amine_200_arm_count"].items():
        if lab.startswith("_"):
            continue
        print("%-22s N=%s  dC/dtail=%s  constant=%s  ARMS=%s"
              % (lab, d["n_nitrogen"], d["dC_per_unit_tail"],
                 d["slope_constant"], d["ARMS_IMPLIED"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
