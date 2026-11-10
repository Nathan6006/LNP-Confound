#!/usr/bin/env python3
"""Enumerate the An-CE-Cm lipids named in the SI, and mass-check them 3 ways.

THE LINKAGE IS A CARBONATE, AND THAT IS THE POINT OF THE PAPER.
A carbonate  R-O-C(=O)-O-R'  is isomeric with an ester bearing one extra
oxygen somewhere else, and carbonate / carbamate / ester are invisible to a
molecular formula. Three independent pieces of evidence pin it here:

  1. THE SYNTHESIS (SI p.3). 4-nitrophenyl chloroformate + alcohol gives the
     alkyl 4-nitrophenyl CARBONATE; the hydroxylamine core's hydroxyls then
     displace 4-nitrophenolate. Both oxygens of the O-CO-O come from the two
     alcohols. A carbamate would need the core's NITROGEN to attack, and the
     core's nitrogens are tertiary amines with no N-H.
  2. 13C (Fig S12, SI p.24). delta 155.16 ppm. A dialkyl carbonate carbonyl is
     ~155; an aliphatic ESTER carbonyl is ~173. The spectrum has NO peak near
     173. 79.32 ppm is the secondary O-CH of the tridecan-7-yl group, pushed
     downfield by the carbonate oxygen.
  3. HRMS. Each carbonate differs from its named carboxylate counterpart by
     exactly one oxygen per tail: A5-CE-C7-6 689.5674 - A5-C7-6 657.5776 =
     31.9898 = 2 x 15.9949, on 2 tails. Checked below.

REGIOCHEMISTRY IS READ, NOT ASSUMED. For a 2-tail lipid the two tails could
sit on one nitrogen or one on each. Figure S12's caption settles it: "2 peaks
at around 53 ppm indicated that the hydrophobic tails were linked to both
amine group sides", and the spectrum lists 52.93 and 52.67. One tail per
nitrogen. This is the Lin-regiochemistry class (identical formula either way,
HRMS blind), and the paper supplies the discriminating evidence.

THE NAME TOKENS ARE NOT CHAIN LENGTHS, they are a carbinol convention:
"Cm-n" is a SECONDARY alcohol whose carbinol carbon is position m of an
(m+n)-carbon chain, so its two arms are (m-1) and n carbons. Pinned by three
HRMS anchors (C5-4 -> C9 alcohol, C7-6 -> C13, C7-7 -> C14) AND, independently
of mass, by the SI naming its carboxylate counterpart's acid: A5-C7-6 is made
from 2-HEXYLOCTANOIC acid, whose alpha carbon bears two hexyls -- exactly the
two hexyl arms of tridecan-7-ol. "Cm" with no second token is a PRIMARY
alcohol of m carbons.

Assembly goes through curation/lib/assemble.attach with atom-mapped attachment
atoms. No canonical SMILES is ever rewritten by string surgery: that is how a
sibling curation silently mislocated every ester while both HRMS anchors still
passed at 0.3-0.5 ppm.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                "..", "lib")))
from rdkit import Chem                                          # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt                   # noqa: E402

import assemble                                                 # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

PROTON = 1.00727646688          # mass of H+ (H atom minus one electron)

# The core, named in full in the SI: N,N,N',N'-tetrakis(2-hydroxyethyl)-
# ethylenediamine. Map 1 and 2 are the two arms of one nitrogen, 3 and 4 the
# two arms of the other, so "one tail per nitrogen" is expressible.
CORE = "[OH:1]CCN(CC[OH:2])CCN(CC[OH:3])CC[OH:4]"
CORE_PLAIN = "OCCN(CCO)CCN(CCO)CCO"

# Which oxygens carry a tail, for each tail count. Four arms are equivalent by
# symmetry in the free tetraol, so 3 tails has only one topology (2 + 1) and
# 4 tails only one. For 2 tails the choice is REAL and is fixed by Fig S12.
SITES = {2: (1, 3), 3: (1, 2, 3), 4: (1, 2, 3, 4)}

# alcohol token -> (arm_a_carbons, arm_b_carbons) for a secondary carbinol,
# or an int for a primary alcohol of that many carbons.
TOKENS = {
    "C6": 6, "C8": 8, "C9": 9, "C10": 10,
    "C5-4": (4, 4), "C6-5": (5, 5), "C7-6": (6, 6), "C7-7": (6, 7),
    "C8-7": (7, 7), "C8-8": (7, 8), "C10-6": (9, 6),
}


def carbonate_tail(token):
    """[C:9](=O)O-R, the acyl half of a carbonate. R comes from TOKENS."""
    spec = TOKENS[token]
    if isinstance(spec, int):
        r = "C" * spec                       # primary: O-CH2(CH2)n-1CH3
    else:
        a, b = spec
        r = f"C({'C' * a}){'C' * b}"         # secondary carbinol, two arms
    return f"[C:9](=O)O{r}"


def ester_tail(token):
    """[C:9](=O)-R', the acyl half of a CARBOXYLATE ester -- one oxygen fewer.

    A5-C10 comes from hendecanoic (undecanoic) acid, so the acyl is C10H21-CO;
    A5-C7-6 from 2-hexyloctanoic acid, so the alpha carbon bears two hexyls.
    Both named verbatim in the SI, p.5.
    """
    spec = TOKENS[token]
    if isinstance(spec, int):
        return "[C:9](=O)" + "C" * spec      # acid has spec+1 carbons
    a, b = spec
    return f"[C:9](=O)C({'C' * a}){'C' * b}"


def build(token, n_tails, kind):
    tail = carbonate_tail(token) if kind == "carbonate" else ester_tail(token)
    smi = CORE
    for site in SITES[n_tails]:
        smi = assemble.attach(smi, tail, a_map=site, b_map=9,
                              label=f"{kind}-{token}-site{site}")
    # Any unreacted arm still carries its map number; a map is not a chemical
    # statement and must not reach the corpus.
    m = Chem.MolFromSmiles(smi)
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    out = Chem.MolToSmiles(m)
    assemble.check_no_radicals(out, label=f"{kind}-{token}")
    return out


# name, alcohol token, n_tails, kind, HRMS calculated [M+H]+ as PRINTED in the
# SI, HRMS found as printed, the SI page it is printed on, and whether the tail
# count is anchored by that HRMS or inferred from the series.
LIBRARY = [
    # ---- carbonates with an HRMS anchor --------------------------------
    ("A5-CE-C8",   "C8",   3, "carbonate", 705.5260, 705.5300, 18, "hrms"),
    ("A5-CE-C9",   "C9",   3, "carbonate", 747.5729, 747.5760, 19, "hrms"),
    ("A5-CE-C10",  "C10",  3, "carbonate", 789.6199, 789.6233, 20, "hrms"),
    ("A5-CE-C5-4", "C5-4", 2, "carbonate", 577.4422, 577.4437, 21, "hrms"),
    ("A5-CE-C7-6", "C7-6", 2, "carbonate", 689.5674, 689.5709, 22, "hrms"),
    ("A5-CE-C7-7", "C7-7", 2, "carbonate", 717.5987, 717.6024, 23, "hrms"),
    # ---- the two extra-tail controls, both with an HRMS anchor ---------
    ("A5-CE-C9-4tails",   "C9",   4, "carbonate", 917.7036, 917.7086, 26, "hrms"),
    ("A5-CE-C7-6-3tails", "C7-6", 3, "carbonate", 915.7607, 915.7665, 27, "hrms"),
    # ---- the two carboxylate counterparts, both with an HRMS anchor ----
    ("A5-C10",  "C10",  3, "ester", 741.6351, 741.6393, 31, "hrms"),
    ("A5-C7-6", "C7-6", 2, "ester", 657.5776, 657.5792, 32, "hrms"),
    # ---- named in Figure S5, NO HRMS, tail count inferred from the series
    ("A5-CE-C6",    "C6",    3, "carbonate", None, None, 17, "inferred_linear"),
    ("A5-CE-C6-5",  "C6-5",  2, "carbonate", None, None, 17, "inferred_branched"),
    ("A5-CE-C8-7",  "C8-7",  2, "carbonate", None, None, 17, "inferred_branched"),
    ("A5-CE-C8-8",  "C8-8",  2, "carbonate", None, None, 17, "inferred_branched"),
    ("A5-CE-C10-6", "C10-6", 2, "carbonate", None, None, 17, "inferred_branched"),
]


def main():
    qc = {"core_smiles": CORE_PLAIN, "n_built": 0,
          "hrms": {"n_compared": 0, "ours_vs_their_calc": [],
                   "their_calc_vs_their_found": [],
                   "ours_vs_their_found": []},
          "linkage_evidence": {
              "13C_carbonyl_ppm": 155.16,
              "13C_source": "Figure S12, SI p.24 (A5-CE-C7-6, CDCl3)",
              "dialkyl_carbonate_expected_ppm": "~155",
              "aliphatic_ester_expected_ppm": "~173",
              "peak_near_173_present": False,
              "13C_OCH_ppm": 79.32,
              "regiochemistry_peaks_ppm": [52.93, 52.67],
              "regiochemistry_claim_in_caption":
                  "2 peaks at around 53 ppm indicated that the hydrophobic "
                  "tails were linked to both amine group sides"},
          }
    rows, blocks = [], []
    for (name, tok, nt, kind, calc, found, page, anchor) in LIBRARY:
        smi = build(tok, nt, kind)
        m = Chem.MolFromSmiles(smi)
        ours_mh = ExactMolWt(m) + PROTON
        formula = Chem.rdMolDescriptors.CalcMolFormula(m)
        inchikey = Chem.MolToInchiKey(m)
        n_unspec = len(Chem.FindMolChiralCenters(
            m, includeUnassigned=True, useLegacyImplementation=False))

        rec = {"IL_name": name, "IL_SMILES": smi, "formula": formula,
               "InChIKey": inchikey, "linkage": kind, "n_tails": nt,
               "alcohol_token": tok, "si_page": page,
               "tail_count_anchor": anchor,
               "ours_MH_plus": round(ours_mh, 4),
               "si_calc_MH_plus": calc, "si_found_MH_plus": found,
               "unassigned_stereocentres": n_unspec}
        if calc is not None:
            d1 = (ours_mh - calc) / calc * 1e6
            d2 = (found - calc) / calc * 1e6
            d3 = (ours_mh - found) / found * 1e6
            rec["ppm_ours_vs_si_calc"] = round(d1, 3)
            rec["ppm_si_found_vs_si_calc"] = round(d2, 3)
            rec["ppm_ours_vs_si_found"] = round(d3, 3)
            qc["hrms"]["n_compared"] += 1
            qc["hrms"]["ours_vs_their_calc"].append(round(d1, 3))
            qc["hrms"]["their_calc_vs_their_found"].append(round(d2, 3))
            qc["hrms"]["ours_vs_their_found"].append(round(d3, 3))
        rows.append(rec)
        blocks.append({"IL_name": name, "core": "A5", "core_name":
                       "N,N,N',N'-tetrakis(2-hydroxyethyl)ethylenediamine",
                       "linker": "carbonate (O-CO-O)" if kind == "carbonate"
                                 else "carboxylate ester (O-CO)",
                       "alcohol_token": tok,
                       "alcohol_or_acid": _reagent(tok, kind),
                       "n_tails": nt,
                       "free_hydroxyls": 4 - nt,
                       "tails_per_nitrogen": "1 + 1" if nt == 2 else
                                             ("2 + 1" if nt == 3 else "2 + 2")})
        qc["n_built"] += 1

    for k in ("ours_vs_their_calc", "their_calc_vs_their_found",
              "ours_vs_their_found"):
        v = qc["hrms"][k]
        qc["hrms"][k + "_max_abs_ppm"] = round(max(abs(x) for x in v), 3)

    # The carbonate-minus-carboxylate oxygen count, as a 4th mass comparison.
    by = {r["IL_name"]: r for r in rows}
    qc["carbonate_minus_carboxylate_Da"] = {}
    for carb, est, nt in [("A5-CE-C10", "A5-C10", 3),
                          ("A5-CE-C7-6", "A5-C7-6", 2)]:
        d = by[carb]["si_calc_MH_plus"] - by[est]["si_calc_MH_plus"]
        qc["carbonate_minus_carboxylate_Da"][f"{carb} - {est}"] = {
            "observed_Da": round(d, 4), "n_tails": nt,
            "expected_n_O_x_15.9949": round(nt * 15.9949146, 4),
            "residual_mDa": round((d - nt * 15.9949146) * 1000, 3)}

    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(blocks[0].keys()))
        w.writeheader()
        for b in blocks:
            w.writerow(b)
    with open(os.path.join(HERE, "02_structure_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)

    for r in rows:
        print(f"{r['IL_name']:<20} {r['formula']:<14} "
              f"{r['ours_MH_plus']:>10.4f}  "
              f"{('%+.2f ppm' % r['ppm_ours_vs_si_calc']) if r['si_calc_MH_plus'] else '   no HRMS':>12}"
              f"  {r['InChIKey']}")
    print(f"\n{qc['n_built']} structures; HRMS compared on "
          f"{qc['hrms']['n_compared']}")
    for k in ("ours_vs_their_calc", "their_calc_vs_their_found",
              "ours_vs_their_found"):
        print(f"  max |{k}| = {qc['hrms'][k + '_max_abs_ppm']} ppm")
    print("  carbonate-carboxylate oxygen check:",
          json.dumps(qc["carbonate_minus_carboxylate_Da"]))


def _reagent(tok, kind):
    spec = TOKENS[tok]
    if kind == "ester":
        if isinstance(spec, int):
            return {10: "hendecanoic (undecanoic) acid, SI p.5"}.get(
                spec, f"C{spec + 1} n-alkanoic acid")
        return "2-hexyloctanoic acid, SI p.5"
    if isinstance(spec, int):
        return {8: "1-octanol, SI p.3", 9: "1-nonanol", 10: "1-decanol",
                6: "1-hexanol"}[spec]
    a, b = spec
    n = a + b + 1
    names = {9: "5-nonanol, SI p.4", 11: "undecan-6-ol",
             13: "tridecan-7-ol", 14: "tetradecan-7-ol",
             15: "pentadecan-8-ol", 16: "hexadecan-8-ol"}
    if tok == "C10-6":
        return "hexadecan-10-ol"
    return names.get(n, f"C{n} secondary alcohol")


if __name__ == "__main__":
    main()
