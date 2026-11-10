#!/usr/bin/env python3
"""ZG_2020 (PMID 33758781) -- enumerate the paper's lipids and check them.

Run:  python3 -I stage01_enumerate.py

WHAT THIS PAPER ACTUALLY VARIES. All 109 LNPs share ONE ionizable lipid,
cKK-E15. The library axis is the HELPER PHOSPHOLIPID (10 levels), crossed with
sterol (3) and molar ratio (4). So this script's main job is not to produce
corpus rows -- there are none -- but to pin the ten phospholipids hard enough
that the regiochemistry finding below is defensible.

EVERY STRUCTURAL INPUT IS TEXT. The ten phospholipids are built from the
PRINTED IUPAC name of intermediate 4 plus the PRINTED fatty acid used in each
Steglich esterification (Appendix S1 pp.2-8 / JATS Methods 3.2), and checked
against the PRINTED HRMS calcd mass. No pixel is read anywhere in this file.

THE ADAMANTYL GROUP SITS ON THE PRIMARY (C3) HYDROXYL, FOR ALL TEN. That is
what the route says -- "dibutyltin oxide-mediated selective PRIMARY alcohol
acylation" gives intermediate 4, whose own printed name puts the adamantylacetoxy
at position 3 and the free OH at position 2; the fatty acid then goes onto that
remaining secondary OH. Two of the ten product names (A-10, A-11) print the two
acyl groups the other way round. Mass cannot see the difference -- the formula
and the exact mass are identical either way -- so this is resolved from the
synthesis, exactly as the Lin regiochemistry was. See 08_contradictions.csv.
"""
import csv, json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import assemble                                            # noqa: E402
from rdkit import Chem                                     # noqa: E402
from rdkit.Chem import rdCIPLabeler, rdMolDescriptors      # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt              # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# Outputs go to OUTDIR, which is HERE for a normal run. The seed harnesses
# redirect it to a temp directory: a MUTATED copy must never overwrite the
# real curation outputs. It did, and the determinism check caught it --
# reference_structures.csv was left holding seed G3's REGIOISOMER anchors.
OUTDIR = HERE

# --- masses used for adduct arithmetic (CODATA/IUPAC, written out so the
#     convention is visible rather than buried in a helper) -----------------
M_NA = 22.9897692820          # 23Na
M_H = 1.00782503207           # 1H
M_E = 0.00054857990907        # electron

# --- intermediate 4, transcribed from its PRINTED IUPAC name ---------------
# "(R)-3-(2-((3r,5r,7r)-Adamantan-1-yl)acetoxy)-2-hydroxypropyl
#  (2-(trimethylammonio)ethyl) phosphate"
#   propyl C1 -> phosphate ester O ; C2 -> free OH ; C3 -> adamantylacetoxy
#
# WHICH ENANTIOMER, DERIVED RATHER THAN COPIED. At C2 the substituents are
# OH, C1H2-O-P, C3H2-O-C(=O)R and H. C1 and C3 tie at (O,H,H), so the branches
# decide: C1's oxygen carries (P,C) and C3's carries (C,C), and P(15) > C(6).
# So the order is OH > C1 > C3 > H. The same order holds in the starting
# material sn-glycero-3-phosphocholine (there C3 is CH2-OH, whose oxygen
# carries (H,C), still below (P,C)) and in every product A-N (where C3's
# oxygen carries (C,C)). The ranking is therefore IDENTICAL at all three
# stages, so no CIP relabel can occur along this route and the printed (R)
# of GPC, of 4 and of all ten products must be the same spatial arrangement.
# We pick the tag RDKit's modern CIP labeller calls R, and assert it below.
COMPOUND_4 = ("[C@@H](COP(=O)([O-])OCC[N+](C)(C)C)(O)"
              "COC(=O)CC12CC3CC(CC(C3)C1)C2")
COMPOUND_4_HRMS = ("C20H36NNaO7P", "[M+Na]+", 456.2122, 456.2118)

# --- the ten products: name -> (acid name as printed, acid SMILES, printed
#     HRMS formula, printed adduct, printed calcd, printed found) -----------
PRODUCTS = [
    ("A-9",      "decanoic acid",      "OC(=O)CCCCCCCCC",
     "C30H54NNaO8P", "[M+Na]+", 610.3479, 610.3468),
    ("A-10",     "undecanoic acid",    "OC(=O)CCCCCCCCCC",
     "C31H56NNaO8P", "[M+Na]+", 624.3636, 624.3624),
    ("A-11",     "lauric acid",        "OC(=O)CCCCCCCCCCC",
     "C32H58NNaO8P", "[M+Na]+", 638.3792, 638.3780),
    ("A-12",     "tridecanoic acid",   "OC(=O)CCCCCCCCCCCC",
     "C33H60NNaO8P", "[M+Na]+", 652.3949, 652.3939),
    ("A-13",     "myristic acid",      "OC(=O)CCCCCCCCCCCCC",
     "C34H63NO8P",   "[M+H]+",  644.4286, 644.4278),
    ("A-14",     "pentadecanoic acid", "OC(=O)CCCCCCCCCCCCCC",
     "C35H64NNaO8P", "[M+Na]+", 680.4262, 680.4253),
    ("A-15",     "palmitic acid",      "OC(=O)CCCCCCCCCCCCCCC",
     "C36H66NNaO8P", "[M+Na]+", 694.4418, 694.4417),
    ("A-16",     "heptadecanoic acid", "OC(=O)CCCCCCCCCCCCCCCC",
     "C37H68NO8P",   "[M+Na]+", 708.4575, 708.4565),
    ("A-17",     "stearic acid",       "OC(=O)CCCCCCCCCCCCCCCCC",
     "C38H70NNaO8P", "[M+Na]+", 722.4731, 722.4723),
    ("A-17-2Z",  "linoleic acid",      r"OC(=O)CCCCCCC/C=C\C/C=C\CCCCC",
     "C38H66NNaO8P", "[M+Na]+", 718.4418, 718.4399),
]


def adduct_mass(neutral_mol, adduct):
    m = ExactMolWt(neutral_mol)
    if adduct == "[M+Na]+":
        return m + M_NA - M_E
    if adduct == "[M+H]+":
        return m + M_H - M_E
    raise ValueError(adduct)


def formula_with_adduct(neutral_mol, adduct):
    """Hill formula of the CATION, so it can be compared to what is printed."""
    rw = Chem.RWMol(neutral_mol)
    if adduct == "[M+Na]+":
        a = Chem.Atom(11)
        a.SetFormalCharge(1)
        rw.AddAtom(a)
    elif adduct == "[M+H]+":
        # protonate the phosphate [O-]; keeps the Hill count honest
        for at in rw.GetAtoms():
            if at.GetSymbol() == "O" and at.GetFormalCharge() == -1:
                at.SetFormalCharge(0)
                at.SetNumExplicitHs(1)
                break
        else:
            raise ValueError("no anionic O to protonate")
    m = rw.GetMol()
    Chem.SanitizeMol(m)
    return rdMolDescriptors.CalcMolFormula(m)


def hill_strip(f):
    """'C30H54NNaO8P+' -> 'C30H54NNaO8P' (drop a trailing charge token)."""
    return f.rstrip("+-").replace("+", "").replace("-", "")


def main():
    out = {"compound_4": {}, "products": [], "checks": {}}

    c4 = Chem.MolFromSmiles(COMPOUND_4)
    if c4 is None:
        raise SystemExit("compound 4 did not parse")
    rdCIPLabeler.AssignCIPLabels(c4)          # modern labeller, not legacy
    cip = [a.GetPropsAsDict().get("_CIPCode")
           for a in c4.GetAtoms() if a.HasProp("_CIPCode")]
    # The paper prints (R) for 4 and for all ten products. Assert it rather
    # than report it: a silently-flipped centre is the F-013 shape. The metric
    # is PRINTED FIRST so a seed aimed at this gate reaches a metric line and
    # can be credited.
    print(f"METRIC compound4_cip={','.join(c or '?' for c in cip)}")
    if cip != ["R"]:
        raise SystemExit(f"compound 4: expected exactly one R centre, got {cip}")
    f_c4, a_c4, calc_c4, found_c4 = COMPOUND_4_HRMS
    got_f = hill_strip(formula_with_adduct(c4, a_c4))
    got_m = adduct_mass(c4, a_c4)
    out["compound_4"] = {
        "smiles": Chem.MolToSmiles(c4),
        "neutral_formula": rdMolDescriptors.CalcMolFormula(c4),
        "cip_codes": cip,
        "printed_formula": f_c4, "our_adduct_formula": got_f,
        "formula_match": got_f == f_c4,
        "printed_calcd": calc_c4, "our_calcd": round(got_m, 4),
        "delta_mDa": round((got_m - calc_c4) * 1000, 2),
        "printed_found": found_c4,
    }

    acid_patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    n_formula_ok = n_mass_ok = 0
    for name, acid_name, acid_smi, pf, adduct, calcd, found in PRODUCTS:
        acid = Chem.MolFromSmiles(acid_smi)
        if acid is None:
            raise SystemExit(f"{name}: acid did not parse")
        # chain length implied by the acid, measured not assumed
        n_c_acid = sum(1 for a in acid.GetAtoms() if a.GetSymbol() == "C")
        # ONE audited assembler; the acid loses its OH, compound 4 keeps its O
        prod_smi = assemble.esterify(acid_smi, COMPOUND_4, label=name)
        prod = Chem.MolFromSmiles(prod_smi)
        rdCIPLabeler.AssignCIPLabels(prod)
        pcip = sorted(a.GetPropsAsDict().get("_CIPCode")
                      for a in prod.GetAtoms() if a.HasProp("_CIPCode"))
        print(f"METRIC cip_{name}={','.join(c or '?' for c in pcip)}")
        if pcip != ["R"]:
            raise SystemExit(f"{name}: expected one R centre, got {pcip}")
        our_f = hill_strip(formula_with_adduct(prod, adduct))
        our_m = adduct_mass(prod, adduct)
        fmt_ok = our_f == pf
        mass_ok = abs(our_m - calcd) < 0.0006          # 0.6 mDa, ~1 ppm here
        n_formula_ok += fmt_ok
        n_mass_ok += mass_ok
        out["products"].append({
            "name": name,
            "acid_printed": acid_name,
            "acid_carbons": n_c_acid,
            "name_token": name.split("-")[1],
            "token_vs_acyl_carbons": f"token {name.split('-')[1]} vs acyl "
                                     f"{n_c_acid} carbons (acid INCLUDES the "
                                     f"carbonyl C; token counts the alkyl "
                                     f"carbons beyond it)",
            "smiles": prod_smi,
            "neutral_formula": rdMolDescriptors.CalcMolFormula(prod),
            "cip_codes": pcip,
            "printed_adduct": adduct,
            "printed_formula": pf, "our_adduct_formula": our_f,
            "formula_match": fmt_ok,
            "printed_calcd": calcd, "our_calcd": round(our_m, 4),
            "delta_mDa": round((our_m - calcd) * 1000, 2),
            "printed_found": found,
            "printed_found_minus_calcd_mDa": round((found - calcd) * 1000, 2),
        })

    out["checks"] = {
        "n_products": len(PRODUCTS),
        "formula_matches": n_formula_ok,
        "mass_matches_within_0.6_mDa": n_mass_ok,
        "adduct_conventions_used": sorted({p[4] for p in PRODUCTS}),
        "all_products_single_R_centre": True,
        "_adamantyl_stereo_is_VACUOUS": (
            "The printed names decorate the cage as (3r,5r,7r)- or (3R,5R,7R)-"
            "adamantan-1-yl. RDKit finds ONE stereocentre in every product, "
            "the glycerol C2. Adamantane's bridgeheads are not stereogenic, so "
            "those descriptors constrain nothing and the stereo assertion "
            "above says nothing about them -- it is vacuous there, not "
            "passing. The paper itself writes them lower-case 'r' in some "
            "entries and upper-case 'R' in others, which is the same point."),
    }

    with open(os.path.join(OUTDIR, "03_enumerated.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    with open(os.path.join(OUTDIR, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "role", "SMILES", "neutral_formula",
                    "exact_mass_neutral", "acyl_donor_printed",
                    "acyl_carbons", "name_token", "cip", "printed_adduct",
                    "printed_formula", "printed_calcd", "our_calcd",
                    "delta_mDa", "printed_found", "ionizable_lipid",
                    "confidence", "structure_assumed", "curation_note"])
        from rdkit.Chem.Descriptors import ExactMolWt as _emw
        for p_ in out["products"]:
            w.writerow([
                p_["name"], "helper_phospholipid (NOT the ionizable lipid)",
                p_["smiles"], p_["neutral_formula"],
                round(_emw(Chem.MolFromSmiles(p_["smiles"])), 6),
                p_["acid_printed"], p_["acid_carbons"], p_["name_token"],
                ",".join(p_["cip_codes"]), p_["printed_adduct"],
                p_["printed_formula"], p_["printed_calcd"], p_["our_calcd"],
                p_["delta_mDa"], p_["printed_found"], "cKK-E15",
                "high", 0,
                "Structure from the printed IUPAC name of intermediate 4 plus "
                "the printed acyl donor; confirmed against the SI HRMS calcd "
                "and against Fig 1c's independent neutral-formula table. "
                "Adamantylacetyl placed on the PRIMARY (C3) hydroxyl per the "
                "synthesis route; A-10 and A-11 are NAMED the other way round "
                "in the SI -- mass cannot see the difference (0.000000 Da). "
                "See 08_contradictions.csv."])
        # the one ionizable lipid, and the excipients, for completeness
        w.writerow(["cKK-E15", "IONIZABLE LIPID (the only one in the paper)",
                    "", "", "", "", "", "", "", "", "", "", "", "", "",
                    "cKK-E15", "low", 1,
                    "DRAWN in Fig 2a but the tail length is not readable: "
                    "each tail is abbreviated to a CnH2n+1 label whose "
                    "subscripts are 2-3 px tall in the only figure "
                    "resolution available. SMILES deliberately left EMPTY "
                    "rather than guessed. The topology IS readable and "
                    "matches corpus cKK-E12. See DECISIONS.md D-2."])

    with open(os.path.join(OUTDIR, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["block_role", "block_name", "block_SMILES", "source"])
        w.writerow(["backbone_intermediate", "compound 4", COMPOUND_4,
                    "Appendix S1 p.2, printed IUPAC name"])
        for name, acid_name, acid_smi, *_ in PRODUCTS:
            w.writerow(["acyl_donor", acid_name, acid_smi,
                        f"Appendix S1, synthesis of {name} (printed reagent)"])

    # ---- THE GATES. Computing a metric and printing it is not a check. ----
    # Exactly one printed formula may differ, and it is named here with the
    # reason, so a NEW mismatch still fails the run.
    KNOWN_FORMULA_TYPO = {
        "A-16": ("printed 'C37H68NO8P' carries no Na although its own adduct "
                 "is [M+Na]+; our C37H68NNaO8P+ reproduces the printed calcd "
                 "708.4575 to -0.02 mDa, so the FORMULA is the typo and the "
                 "MASS is right. Deciding which half is wrong from the number, "
                 "not from preference (the AC_2025 lesson).")}
    bad_f = [p["name"] for p in out["products"]
             if not p["formula_match"] and p["name"] not in KNOWN_FORMULA_TYPO]
    bad_m = [p["name"] for p in out["products"]
             if abs(p["delta_mDa"]) >= 0.6]
    print(f"METRIC formula_mismatches_unexplained={len(bad_f)} {bad_f}")
    print(f"METRIC mass_failures_over_0.6mDa={len(bad_m)} {bad_m}")
    print(f"METRIC max_abs_delta_mDa="
          f"{max(abs(p['delta_mDa']) for p in out['products']):.3f}")
    out["checks"]["known_formula_typo"] = KNOWN_FORMULA_TYPO
    out["checks"]["formula_mismatches_unexplained"] = bad_f
    out["checks"]["mass_failures_over_0.6_mDa"] = bad_m
    out["checks"]["max_abs_delta_mDa"] = max(abs(p["delta_mDa"])
                                             for p in out["products"])
    with open(os.path.join(OUTDIR, "03_enumerated.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    if bad_f or bad_m:
        raise SystemExit(
            f"GATE FAILED: formula mismatches {bad_f}, mass failures {bad_m}")

    print(json.dumps(out["checks"], indent=1))
    print("\ncompound 4:", json.dumps(out["compound_4"], indent=1))
    for p in out["products"]:
        print(f"{p['name']:9s} {p['neutral_formula']:14s} "
              f"printed {p['printed_formula']:14s} {p['printed_adduct']:8s} "
              f"fmt={'OK ' if p['formula_match'] else 'MISMATCH'} "
              f"calcd {p['printed_calcd']:9.4f} ours {p['our_calcd']:9.4f} "
              f"d={p['delta_mDa']:+7.2f} mDa  CIP={p['cip_codes']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
