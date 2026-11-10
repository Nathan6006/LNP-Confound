#!/usr/bin/env python3
"""Stage 03 - what structure this paper actually determines, and what it does not.

The library is <amine>D<maleate>: a one-step aza-Michael mono-addition of an
amine N-H onto a dialkyl maleate, neat, 80 C, 24 h (Supplementary Figure 2).
20 amines x 6 maleates = 120 members, plus 13 catalogue reagents = 133.

THE PAPER DRAWS 3 OF THE 120. Supplementary Figure 2 draws the amine, the
maleate and the product for 17D9.2, 18D4 and 20D8i, and prints a mass spectrum
(Supp Figs 4-6) and a 1H NMR (Supp Figs 7-9) for those same three. Every other
member's amine or maleate is defined only in reference 34 (Xu, Han, Xu,
Ricciardi et al., "Plug-and-play assembly of biodegradable ionizable lipids...",
bioRxiv Mar 2025), which is a DIFFERENT PAPER and is not in source/.

THE ARITHMETIC THIS FILE RUNS. For a product = amine + k x maleate, the printed
[M+H]+/[M+2H]2+/[M+3H]3+ pin the TOTAL CARBON COUNT of the two ester chains:

    a + b = (M - amine - k*116.072) / (14.027 * k)

That is an integer if and only if the read is right, so the fractional part is
the residual. It is validated on 18D4, where the answer is independently known
(the maleate is drawn as plainly n-butyl on both sides and the token is "4").

WHAT IT CANNOT DO: mass cannot see connectivity. a+b = 19 is satisfied by
11+8 and by 10+9 alike, and by every branch position of each. So a confirmed
carbon count is NOT a determined structure, and this file says so per compound
rather than emitting one.
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, os.path.join(HERE, ".."))

from rdkit import Chem, RDLogger                          # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors      # noqa: E402
RDLogger.DisableLog("rdApp.*")

# Average atomic masses, for the arithmetic above.
MALEATE_CORE = 114.0557        # C4H2O4, the diester core without the two R
CH2 = 14.0266
H = 1.008

# m/z as PRINTED in the SI mass spectra. Read off the rasters at 2x; each peak
# carries its own ion label in the figure, so the charge state is the paper's
# statement, not our inference.
MS = {
    "17D9.2": {"fig": "Supplementary Figure 4 (SI p.18)",
               "peaks": {"[M+2H]2+": 648.123, "[M+H]+": 1294.693,
                         "[M+Na]+": 1316.125}},
    "18D4":   {"fig": "Supplementary Figure 5 (SI p.19)",
               "peaks": {"[M+2H]2+": 451.589, "[M+H]+": 901.829,
                         "[M+Na]+": 923.681}},
    "20D8i":  {"fig": "Supplementary Figure 6 (SI p.20)",
               "peaks": {"[M+3H]3+": 627.269, "[M+2H]2+": 940.195}},
}
NA, PROTON = 22.9898, 1.008


def neutral_mass(label, mz):
    if label == "[M+H]+":
        return mz - PROTON
    if label == "[M+Na]+":
        return mz - NA
    if label.startswith("[M+") and label.endswith("H]2+"):
        return 2 * mz - 2 * PROTON
    if label.startswith("[M+") and label.endswith("H]3+"):
        return 3 * mz - 3 * PROTON
    raise ValueError(label)


def avg_mw(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise SystemExit("unparseable: " + smiles)
    return Descriptors.MolWt(m), rdMolDescriptors.CalcMolFormula(m)


def chain_carbons(M, amine_mw, k):
    """Total carbons across the TWO ester chains of each maleate arm."""
    per_arm = (M - amine_mw) / k
    return (per_arm - MALEATE_CORE - 2 * H) / CH2


def main():
    out, rep = [], {}

    # --- the amines the paper DRAWS -------------------------------------
    amine17 = "NCCN(CCN)CCN"                       # TREN, Supp Fig 2
    amine18 = "NCCNCCN1CCN(CCN)CC1"                # Supp Fig 2
    a17_mw, a17_f = avg_mw(amine17)
    a18_mw, a18_f = avg_mw(amine18)
    rep["amines_drawn"] = {
        "17": {"smiles": amine17, "formula": a17_f, "avg_mw": round(a17_mw, 3),
               "name": "tris(2-aminoethyl)amine (TREN)"},
        "18": {"smiles": amine18, "formula": a18_f, "avg_mw": round(a18_mw, 3),
               "name": "N-(2-aminoethyl)-2-[4-(2-aminoethyl)piperazin-1-yl]"
                       "ethan-1-amine"},
        "20": {"smiles": None,
               "why": "a branched polyamidoamine with multiple -CH2CH2C(=O)NH-"
                      "CH2CH2NH2 arms. Drawn in Supp Fig 2, but the core and "
                      "the arm count were NOT resolved at the 137 dpi this SI "
                      "offers, and no reading is recorded rather than a guess."},
    }

    # --- the arithmetic, per compound -----------------------------------
    arith = {}
    for name, amine_mw, k in [("18D4", a18_mw, 3), ("17D9.2", a17_mw, 3)]:
        rows = []
        for lab, mz in MS[name]["peaks"].items():
            M = neutral_mass(lab, mz)
            ab = chain_carbons(M, amine_mw, k)
            rows.append({"ion": lab, "mz_printed": mz,
                         "neutral_M": round(M, 3),
                         "a_plus_b": round(ab, 3),
                         "residual_from_integer": round(ab - round(ab), 3)})
        Ms = [r["neutral_M"] for r in rows]
        arith[name] = {
            "figure": MS[name]["fig"], "arms_k": k,
            "amine_avg_mw": round(amine_mw, 3),
            "per_ion": rows,
            "instrument_internal_spread_Da": round(max(Ms) - min(Ms), 3),
            "a_plus_b_mean": round(sum(r["a_plus_b"] for r in rows) / len(rows), 3),
        }
    # 20D8i: the amine is unread, so only M is reported
    Ms = [neutral_mass(l, m) for l, m in MS["20D8i"]["peaks"].items()]
    arith["20D8i"] = {"figure": MS["20D8i"]["fig"],
                      "neutral_M_per_ion": [round(x, 3) for x in Ms],
                      "instrument_internal_spread_Da": round(max(Ms) - min(Ms), 3),
                      "a_plus_b": None,
                      "why": "the amine's formula is unread, so the arithmetic "
                             "has two unknowns and is not solvable."}
    rep["ms_arithmetic"] = arith

    # --- 18D4, the one compound this paper fully determines -------------
    # Built by plain SMILES composition from the two drawn reagents: three
    # aspartate-diester arms, one on each of the amine's three N-H nitrogens
    # (2 primary + 1 secondary), the two piperazine tertiary N untouched.
    # "Favor mono-addition" is the scheme's own words and the drawing shows
    # exactly three arms and six butyl esters.
    ARM = "C(CC(=O)OCCCC)C(=O)OCCCC"
    d18 = "C(CC(=O)OCCCC)(C(=O)OCCCC)N(CCN1CCN(CC1)CCN%s)CCN%s" % (ARM, ARM)
    mw, formula = avg_mw(d18)
    can = Chem.MolToSmiles(Chem.MolFromSmiles(d18))
    rep["18D4"] = {
        "smiles": can, "formula": formula, "avg_mw": round(mw, 3),
        "inchikey": Chem.MolToInchiKey(Chem.MolFromSmiles(can)),
        "predicted_MH": round(mw + PROTON, 3),
        "printed_MH": MS["18D4"]["peaks"]["[M+H]+"],
        "delta_Da": round(MS["18D4"]["peaks"]["[M+H]+"] - (mw + PROTON), 3),
        "nearest_alternative_Da": round(2 * CH2, 3),
        "note": "one CH2 on each of two chains is the nearest alternative; the "
                "discrimination margin is that over the delta.",
    }

    # --- the per-member determinability table ---------------------------
    blocks = list(csv.DictReader(open(os.path.join(HERE, "01_blocks.csv"))))
    MEASURED = {"2D8": 9, "7D4": 37, "17D8i": 94, "17D9.2": 96, "18D4": 97,
                "20D8i": 112, "SM-102": 125}
    DRAWN_AMINE = {"17", "18", "20"}
    DRAWN_MALEATE = {"9.2", "4", "8i"}
    import re
    pat = re.compile(r"^(\d+)D(4|6\.2|8|8i|18|9\.2)$")
    for b in blocks:
        n = b["IL_name"]
        m = pat.match(n)
        if m:
            am, ma = m.group(1), m.group(2)
            amine_ok = am in DRAWN_AMINE
            mal_ok = ma in DRAWN_MALEATE
            if n in ("17D9.2", "18D4", "20D8i"):
                det = "drawn_in_full"
            elif amine_ok and mal_ok:
                det = "composable_from_two_other_drawings"
            else:
                det = "NOT_determinable"
            why = ("amine %s %s drawn; maleate %s %s drawn"
                   % (am, "is" if amine_ok else "is NOT",
                      ma, "is" if mal_ok else "is NOT"))
        else:
            det = "catalogue_reagent"
            why = "a named catalogue lipid; structure from the reagent registry"
        out.append({"Index": b["Index"], "IL_name": n,
                    "class": "nDm_library" if m else "catalogue",
                    "determinability": det, "why": why,
                    "individually_measured_in_Fig3": MEASURED.get(n, ""),
                    "emitted_SMILES": can if n == "18D4" else "",
                    })

    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)

    from collections import Counter
    rep["determinability_counts"] = dict(Counter(r["determinability"] for r in out))
    rep["among_the_7_individually_measured"] = {
        r["IL_name"]: r["determinability"] for r in out
        if r["individually_measured_in_Fig3"] != ""}
    with open(os.path.join(HERE, "03_structures.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    print(json.dumps(rep, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
