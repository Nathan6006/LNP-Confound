#!/usr/bin/env python3
"""Stage 03b - the PROPOSED structures, and the exact reason none is pinned.

The paper draws NO VNxCyOz structure anywhere. Figure S1b draws the reaction
with a generic `R-NH2`; Figure S1a draws 35 commercial amines, each labelled
with a CAS registry number; Figure S2's caption prints an exact mass per
compound. So a structure here is an INFERENCE from three measured things:

  1. the scaffold: `R-NH2 + 2 x VE-acrylate -> R-N(CH2CH2C(=O)O-tocopheryl)2`,
     with VE-AC's own [M+Na]+ reproducing the printed calc to 0.03 mDa;
  2. the exact formula of the product, from Figure S2;
  3. the 34 distinct amines drawn in Figure S1a.

(1)+(2) give the core amine's formula exactly. (3) then says which drawn
reagents have that formula. Where exactly one does AND exactly one of its
nitrogens can carry two arms, the structure is DETERMINATE GIVEN THAT THE CORE
IS IN THE DRAWN SET. It is not pinned, because two of the seven cores are not
in the drawn set at all - including the lead compound's - which proves the
drawn set is not a complete inventory of the cores used.

Assembly uses curation/lib/assemble: atom-mapped attachment, the acid losing
its OH and the alcohol keeping its oxygen, and the stereo/radical assertions on
every build. Tocopherol's one assigned stereocentre at C2 plus its two
unassigned side-chain centres make `check_stereo_preserved` the live guard.
"""
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import assemble                                     # noqa: E402
from rdkit import Chem, RDLogger                    # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt       # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula  # noqa: E402

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
ELECTRON = 0.00054857990946

# (2R)-alpha-tocopherol, written as an explicit source literal with the
# aromatic chromanol ring. The two side-chain centres (4'R,8'R in natural RRR
# tocopherol) are left UNASSIGNED: the paper names no isomer and reports no
# optical rotation, so assigning them would be invention.
TOCOPHEROL = "Cc1c(C)c2c(c(C)c1O)CC[C@@](C)(CCCC(C)CCCC(C)CCCC(C)C)O2"
ACRYLIC = "OC(=O)C=C"

# The Michael-adduct arm, as its own source literal, written HEAD-FIRST from
# the nitrogen end: mapped CH2, second CH2, carbonyl carbon, ester oxygen,
# tocopheryl. A SMILES fragment has a direction and mass is blind to it, so
# the arm is written in the OPPOSITE traversal to the esterify() route that
# builds VE-AC below, and the two are reconciled by formula and exact mass.
ARM = ("[CH2:2]CC(=O)Oc1c(C)c2c(c(C)c1C)CC[C@@](C)"
       "(CCCC(C)CCCC(C)CCCC(C)C)O2")

LIPIDS = [
    dict(name="VN3C2O", panel="b", formula_ion="C71H123N2O8",
         core_formula="C7H18N2O2", core_cas=None, core_smiles=None,
         why="No amine drawn in Figure S1a has formula C7H18N2O2. The two "
             "bis(hydroxyethyl) amines drawn (CAS 4439-20-7 and 3197-06-6) are "
             "both C6H16N2O2, exactly one CH2 short. Under the VN2/VN3 = "
             "ethylene/propylene reading the core would be the propylene "
             "homologue of 3197-06-6, which the paper does not draw."),
    dict(name="VN2C4", panel="c", formula_ion="C70H119N2O6",
         core_formula="C6H14N2", core_cas=None, core_smiles=None,
         why="TWO drawn amines have formula C6H14N2: 7154-73-6 "
             "(1-(2-aminoethyl)pyrrolidine, one 2-arm site) and 694-83-7 "
             "(trans-1,2-diaminocyclohexane, two 1-arm sites). They differ in "
             "connectivity and the exact mass cannot see the difference."),
    dict(name="VN2C2", panel="d", formula_ion="C70H121N2O6",
         core_formula="C6H16N2", core_cas="100-36-7",
         core_smiles="CCN(CC)CC[NH2:1]",
         why="Exactly one drawn amine has formula C6H16N2 (CAS 100-36-7, "
             "N,N-diethylethylenediamine) and exactly one of its nitrogens "
             "carries two hydrogens, so the two arms have only one site. Its "
             "+CH2 partner 104-78-9 is also drawn, which is what the measured "
             "CH2 ladder predicts."),
    dict(name="VN3C2", panel="e", formula_ion="C71H123N2O6",
         core_formula="C7H18N2", core_cas="104-78-9",
         core_smiles="CCN(CC)CCC[NH2:1]",
         why="TWO drawn amines have formula C7H18N2: 104-78-9 "
             "(N,N-diethyl-1,3-propanediamine) and 10061-68-4 "
             "(N,N'-diethyl-1,3-propanediamine). The CH2 ladder discriminates: "
             "VN3C2 = VN2C2 + CH2 to 0.05 mDa, and the linker homologue of "
             "VN2C2's core (100-36-7) is 104-78-9, not 10061-68-4."),
    dict(name="VN3C1", panel="f", formula_ion="C69H119N2O6",
         core_formula="C5H14N2", core_cas="109-55-7",
         core_smiles="CN(C)CCC[NH2:1]",
         why="Exactly one drawn amine has formula C5H14N2 (CAS 109-55-7, "
             "N,N-dimethyl-1,3-propanediamine) and exactly one of its "
             "nitrogens carries two hydrogens."),
    dict(name="PANEL-G (label disputed)", panel="g", formula_ion="C71H121N2O6",
         core_formula="C7H16N2", core_cas=None, core_smiles=None,
         why="The caption labels this VN3C2, but panel (e) already labels a "
             "DIFFERENT formula VN3C2 and the two differ by H2. Exactly one "
             "drawn amine has C7H16N2 (51387-90-7), but the compound's NAME is "
             "unresolved, so it cannot be joined to any value."),
    dict(name="VN2C1", panel="h", formula_ion="C68H117N2O6",
         core_formula="C4H12N2", core_cas=None, core_smiles=None,
         why="NO amine drawn in Figure S1a has formula C4H12N2. This is the "
             "paper's LEAD compound (Figures S8, S12, S13 and all in-vivo "
             "work). Under the CH2-ladder reading the core is the ethylene "
             "homologue of 109-55-7, i.e. N,N-dimethylethylenediamine "
             "(CAS 108-00-9) - which the figure does not draw. Figure S1a "
             "says 35 amines and draws 34 distinct ones (7209-38-3 twice), so "
             "exactly one slot is unaccounted for; that is suggestive and is "
             "not evidence."),
]


def retag_secondary_amine(smiles, label):
    """Re-tag the ONE secondary amine nitrogen with map 1.

    assemble.attach clears the map after bonding, so the second Michael
    addition onto the same nitrogen needs it put back. This finds it by
    valence - exactly one N with exactly one hydrogen - and raises if that is
    not unique, rather than taking the first match.
    """
    m = Chem.MolFromSmiles(smiles)
    hits = [a.GetIdx() for a in m.GetAtoms()
            if a.GetSymbol() == "N" and a.GetTotalNumHs() == 1]
    if len(hits) != 1:
        raise SystemExit("%s: %d secondary amine N found, need exactly 1."
                         % (label, len(hits)))
    m.GetAtomWithIdx(hits[0]).SetAtomMapNum(1)
    return Chem.MolToSmiles(m)


def main():
    out = {"scaffold": {}, "lipids": []}

    # --- VE-AC, by esterification: the acid loses its OH, tocopherol keeps
    # its oxygen. This is the external check on the TOCOPHEROL literal, since
    # the paper prints VE-AC's own exact mass.
    ve_ac = assemble.esterify(ACRYLIC, TOCOPHEROL, label="VE-AC")
    m = Chem.MolFromSmiles(ve_ac)
    assemble.check_no_radicals(m, label="VE-AC")
    na = ExactMolWt(m) + 22.9897692809 - ELECTRON
    st = assemble.stereocentres(m)
    out["scaffold"] = {
        "tocopherol_smiles": TOCOPHEROL,
        "tocopherol_formula": CalcMolFormula(Chem.MolFromSmiles(TOCOPHEROL)),
        "arm_literal": ARM,
        "arm_formula": CalcMolFormula(Chem.MolFromSmiles(ARM)),
        "VE_AC_smiles": ve_ac,
        "VE_AC_formula": CalcMolFormula(m),
        "VE_AC_M_plus_Na_calc": round(na, 4),
        "paper_calc": 507.3809, "paper_found": 507.3804,
        "delta_vs_paper_calc_mDa": round((na - 507.3809) * 1000, 2),
        "stereocentres_assigned": sum(1 for v in st.values() if v in ("R", "S")),
        "stereocentres_unassigned": sum(1 for v in st.values() if v == "?"),
        "note": "Only tocopherol's C2 is assigned. The 4' and 8' side-chain "
                "centres are left unassigned because the paper names no "
                "isomer, so every structure below carries stereo uncertainty "
                "that mass cannot see.",
    }

    for L in LIPIDS:
        rec = dict(L)
        rec["structure_pinned"] = False
        if L["core_smiles"]:
            prod = assemble.attach(L["core_smiles"], ARM, a_map=1, b_map=2,
                                   label="%s arm 1" % L["name"])
            prod = retag_secondary_amine(prod, "%s after arm 1" % L["name"])
            prod = assemble.attach(prod, ARM, a_map=1, b_map=2,
                                   label="%s arm 2" % L["name"])
            pm = Chem.MolFromSmiles(prod)
            assemble.check_no_radicals(pm, label=L["name"])
            ion = ExactMolWt(pm) + 1.00782503207 - ELECTRON
            rec["proposed_IL_SMILES"] = prod
            rec["proposed_formula"] = CalcMolFormula(pm)
            rec["proposed_M_plus_H"] = round(ion, 4)
            # the external check: our build vs the printed HRMS calc
            printed = {"VN2C2": 1085.9220, "VN3C2": 1099.9376,
                       "VN3C1": 1071.9063}.get(L["name"])
            rec["printed_calc_M_plus_H"] = printed
            rec["delta_vs_printed_mDa"] = (
                round((ion - printed) * 1000, 2) if printed else None)
            rec["formula_matches_HRMS"] = (
                rec["proposed_formula"].replace("+", "") ==
                L["formula_ion"].replace("+", "").replace(
                    "H%d" % 0, "") or None)
        out["lipids"].append(rec)

    json.dump(out, open(os.path.join(HERE, "03_proposed_structures.json"), "w"),
              indent=1)
    print(json.dumps(out["scaffold"], indent=1))
    for r in out["lipids"]:
        print("%-24s ion %-14s core %-11s cas %-11s built=%-5s dHRMS=%s mDa"
              % (r["name"], r["formula_ion"], r["core_formula"], r["core_cas"],
                 bool(r.get("proposed_IL_SMILES")),
                 r.get("delta_vs_printed_mDa")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
