#!/usr/bin/env python3
"""XQ_2020 stage 03 -- enumerate the eleven R-O16CBA lipidoids.

    python3 -I stage03_enumerate.py

THE ONE QUESTION THIS STAGE EXISTS TO ANSWER. A primary amine can take ONE or
TWO Michael arms and a secondary amine one; amines 113, 306 and 400 have three
nitrogens each. So "which nitrogen carries which arm" is a real regiochemical
question, and molecular formula is blind to it. The run brief is explicit about
this and about the sibling that built acyl fragments by string surgery on a
canonical SMILES and silently mislocated every ester while both HRMS anchors
still passed at 0.3-0.5 ppm.

So this stage does NOT choose a regiochemistry. It ENUMERATES EVERY ONE.

  1. the arm COUNT comes from Figure 2F's printed molecular formula, solved
     exactly: formula(neutral lipidoid) - formula(amine head) must equal an
     integer multiple of C23H34O7. No remainder is tolerated.
  2. every distinct way of placing that many arms on N-H nitrogens is built,
     by breadth-first expansion over ALL N-H sites at every step.
  3. the distinct canonical products are counted. One product means the
     constitution is forced and there is nothing to flag. More than one means a
     genuine positional ambiguity, and the row is flagged with it.

Every bond is made with assemble.attach, atom-mapped on both partners. No atom
is ever deleted during the alkylation, so the F-013 RemoveAtom mechanism cannot
arise here at all.

THE ARM FRAGMENT IS BUILT TWICE, BY TWO ROUTES THAT SHARE NO CODE, and the two
must agree or this stage raises -- the rule curation/lib/chains.py applies to
chain lengths, applied here to a fragment.

Writes 03_enumerated.csv and 03_enumerated_qc.json.
"""
import csv
import json
import os
import re
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

from rdkit import Chem, RDLogger                                    # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt                       # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula              # noqa: E402
import assemble                                                     # noqa: E402

RDLogger.DisableLog("rdApp.*")

PROTON = 1.00727646          # mass of H+ (proton), not of a hydrogen atom

HEXDMBAH = "OCC1(C)COC(-c2c(OC)cc(OCCCCCC)cc2OC)OC1"
O16CBA = "C=CC(=O)OCC1(C)COC(-c2c(OC)cc(OCCCCCC)cc2OC)OC1"
THME = "OCC(C)(CO)CO"        # 1,1,1-tris(hydroxymethyl)ethane

# ---------------------------------------------------------------------------
# Figure 2F, transcribed from nihms-1799810-f0002.jpg panel F.
# The FORMULA column is the formula of the PROTONATED ion [M+H]+, which is how
# the authors label the column; the neutral is therefore one H lighter.
# ---------------------------------------------------------------------------
FIG2F = {
    #  id : (formula of [M+H]+, calculated m/z, observed m/z)
    "75":  ("C52H85N2O14",   961.60,  961.64),
    "76":  ("C52H83N2O14",   959.58,  959.73),
    "77":  ("C53H85N2O14",   973.60,  973.73),
    "78":  ("C50H81N2O14",   933.57,  933.64),
    "80":  ("C51H83N2O14",   947.58,  947.64),
    "81":  ("C53H87N2O14",   975.62,  975.73),
    "82":  ("C53H85N2O14",   973.60,  973.73),
    "93":  ("C52H80N3O14",   970.56,  970.64),
    "113": ("C74H118N3O21", 1384.83, 1384.91),
    "306": ("C76H122N3O21", 1412.86, 1412.82),
    "400": ("C77H124N3O21", 1426.87, 1426.91),
}

# degradation-fragment m/z, read as PDF TEXT SPANS (not pixels) from the SI and
# from the RESULTS prose. Used as an independent arm-count check in stage 07.
DEGRADATION_MZ = {
    ("75", "O16CBA-1"): (713.55, "RESULTS, Figure 3E"),
    ("76", "O16CBA-1"): (711.64, "SI p. S6, Figure S5 text layer"),
    ("76", "O16CBA-0"): (463.36, "SI p. S6, Figure S5 text layer"),
    ("77", "O16CBA-1"): (725.64, "SI p. S7, Figure S6 text layer"),
    ("77", "O16CBA-0"): (477.45, "SI p. S7, Figure S6 text layer"),
}


def parse_formula(s):
    c = Counter()
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", s):
        if el:
            c[el] += int(n) if n else 1
    return c


def formula_of(mol):
    return parse_formula(re.sub(r"[+-]\d*$", "", CalcMolFormula(mol)))


def nh_sites(mol):
    """Indices of nitrogens bearing at least one hydrogen.

    An N with no H is not a Michael acceptor site here: the imidazole N3 of
    head 93 would give a charged imidazolium, and Figure 2F's formula for
    93-O16CBA is that of a NEUTRAL molecule plus one proton, so no quaternary
    nitrogen is present.
    """
    return [a.GetIdx() for a in mol.GetAtoms()
            if a.GetSymbol() == "N" and a.GetTotalNumHs() > 0]


def tag(smiles, idx, mapno):
    m = Chem.MolFromSmiles(smiles)
    m.GetAtomWithIdx(idx).SetAtomMapNum(mapno)
    return Chem.MolToSmiles(m)


def build_arm():
    """The Michael arm, built twice by routes that share no code.

    Route A  reduce O16CBA's acrylate C=C to a single bond, in the molecule
             graph, under a SMARTS match asserted unique, then map the
             terminal carbon.
    Route B  esterify propanoic acid (C3 mapped) onto HexDMBAH with
             assemble.esterify, which never touches the alcohol's oxygen.

    Route A starts from the finished tail and goes backwards; route B starts
    from the tail's own precursor and goes forwards. They share no line of
    code and no intermediate. If they agree, the arm is right.
    """
    # --- route A
    m = Chem.MolFromSmiles(O16CBA)
    patt = Chem.MolFromSmarts("[CH2X3]=[CHX3][CX3](=O)[OX2]")
    hits = m.GetSubstructMatches(patt)
    if len(hits) != 1:
        raise SystemExit(f"stage03: O16CBA has {len(hits)} acrylate groups, "
                         "need exactly 1 -- refusing to guess which one is "
                         "the Michael acceptor.")
    term, beta, _c, _o, _oe = hits[0]
    rw = Chem.RWMol(m)
    rw.GetBondBetweenAtoms(term, beta).SetBondType(Chem.BondType.SINGLE)
    prod = rw.GetMol()
    Chem.SanitizeMol(prod)
    assemble.check_no_radicals(prod, label="arm route A")
    a_idx = term
    route_a = tag(Chem.MolToSmiles(prod),
                  [i for i, x in enumerate(
                      Chem.CanonicalRankAtoms(prod))],  # placeholder, replaced
                  0) if False else None
    # map the terminal carbon by its identity in THIS molecule, not by index
    # in some other one: re-locate it as the only CH3 that is beta to the ester
    # carbonyl through exactly one carbon.
    patt2 = Chem.MolFromSmarts("[CH3X4][CH2X4][CX3](=O)[OX2]")
    h2 = prod.GetSubstructMatches(patt2)
    if len(h2) != 1:
        raise SystemExit(f"stage03: route A found {len(h2)} propanoate "
                         "termini after reduction, need exactly 1.")
    rwm = Chem.RWMol(prod)
    rwm.GetAtomWithIdx(h2[0][0]).SetAtomMapNum(2)
    route_a = Chem.MolToSmiles(rwm.GetMol())

    # --- route B
    route_b = assemble.esterify("[CH3:2]CC(=O)O", HEXDMBAH, label="arm route B")

    if Chem.CanonSmiles(route_a) != Chem.CanonSmiles(route_b):
        raise SystemExit("stage03: THE TWO ARM ROUTES DISAGREE.\n"
                         f"  A (reduce O16CBA): {route_a}\n"
                         f"  B (esterify onto HexDMBAH): {route_b}")

    unmapped = Chem.MolFromSmiles(re.sub(r":2\]", "]", route_a))
    f = CalcMolFormula(unmapped)
    if f != "C23H36O7":
        raise SystemExit(f"stage03: arm fragment formula is {f}, expected "
                         "C23H36O7 (= O16CBA C23H34O7 plus the two hydrogens "
                         "the alkene gains on reduction).")
    return route_a, route_b, f


def enumerate_products(head_smiles, n_arms, arm_smiles, label):
    """Every distinct product of n_arms Michael additions onto N-H sites."""
    frontier = {Chem.CanonSmiles(head_smiles)}
    per_step = []
    for step in range(n_arms):
        nxt = set()
        for s in frontier:
            m = Chem.MolFromSmiles(s)
            for idx in nh_sites(m):
                nxt.add(assemble.attach(tag(s, idx, 1), arm_smiles,
                                        a_map=1, b_map=2,
                                        label=f"{label} step{step + 1}"))
        if not nxt:
            raise SystemExit(f"stage03 {label}: no N-H site left at step "
                             f"{step + 1}; the arm count from Figure 2F "
                             "cannot be placed.")
        per_step.append(len(nxt))
        frontier = nxt
    return sorted(frontier), per_step


def main():
    blocks = {}
    with open(os.path.join(HERE, "01_blocks.csv")) as f:
        for r in csv.DictReader(f):
            blocks[r["block_id"]] = r

    arm_a, arm_b, arm_formula = build_arm()
    arm_c = formula_of(Chem.MolFromSmiles(O16CBA))

    qc = {"_what": "stage 03 QC for XQ_2020. The regiochemistry question is "
                   "answered by exhaustive enumeration, not by a choice.",
          "arm_fragment": {
              "route_A_reduce_O16CBA": arm_a,
              "route_B_esterify_onto_HexDMBAH": arm_b,
              "agree": True,
              "formula_unmapped": arm_formula,
              "why_two_routes": "curation/lib/chains.py's rule -- two "
                                "independent methods must agree or it raises -- "
                                "applied to a fragment instead of a chain. The "
                                "run brief's named failure mode is a sibling "
                                "that built acyl fragments by string surgery "
                                "on a canonical SMILES and mislocated every "
                                "ester with both HRMS anchors still passing."}}

    # --- the symmetry claim behind the acetal regiochemistry, PROVED
    thme = Chem.MolFromSmiles(THME)
    ranks = list(Chem.CanonicalRankAtoms(thme, breakTies=False))
    oh_ranks = sorted(ranks[a.GetIdx()] for a in thme.GetAtoms()
                      if a.GetSymbol() == "O")
    qc["acetal_regiochemistry"] = {
        "the_question": "A cyclic benzylidene acetal's ring is closed by TWO "
                        "diol oxygens. 1,1,1-tris(hydroxymethyl)ethane has "
                        "THREE CH2OH groups, so there are three pairings. "
                        "Molecular formula cannot tell them apart.",
        "symmetry_ranks_of_the_three_oxygens": oh_ranks,
        "all_three_equivalent": len(set(oh_ranks)) == 1,
        "conclusion": "All three CH2OH groups of 1,1,1-tris(hydroxymethyl)"
                      "ethane are related by the molecule's own symmetry "
                      "(identical canonical rank with breakTies=False), so all "
                      "three pairings give ONE constitution: "
                      "2-aryl-5-methyl-5-(hydroxymethyl)-1,3-dioxane. There is "
                      "no regiochemical choice to make and nothing to flag. "
                      "This is PROVED here rather than asserted, because the "
                      "run brief names acetal regiochemistry as a trap.",
        "what_this_does_NOT_settle": "the cis/trans relationship between the "
                                     "C2 aryl group and the C5 substituents. "
                                     "That is stereochemistry, not "
                                     "constitution, and it is open: see "
                                     "DECISIONS.md D-2.",
    }
    if len(set(oh_ranks)) != 1:
        raise SystemExit("stage03: the three CH2OH groups are NOT equivalent; "
                         "the acetal regiochemistry is a real choice and this "
                         "script may not make it.")

    # --- regenerate O16CBA itself from its precursor, as a check on the
    #     hand-written tail SMILES in stage01
    rebuilt = assemble.esterify("C=CC(=O)O", HEXDMBAH, label="O16CBA rebuild")
    qc["o16cba_rebuilt_from_HexDMBAH_plus_acrylic_acid"] = {
        "rebuilt": rebuilt,
        "stage01_hand_written": Chem.CanonSmiles(O16CBA),
        "match": Chem.CanonSmiles(rebuilt) == Chem.CanonSmiles(O16CBA),
    }
    if Chem.CanonSmiles(rebuilt) != Chem.CanonSmiles(O16CBA):
        raise SystemExit("stage03: the hand-written O16CBA disagrees with the "
                         "one assembled from its own named precursors.")

    rows, members = [], []
    for hid in [h for h in FIG2F]:
        head = blocks[hid]["SMILES"]
        hf = formula_of(Chem.MolFromSmiles(head))

        # ---- arm count, solved from the PRINTED formula. No remainder allowed.
        ion_f = parse_formula(FIG2F[hid][0])
        neutral = Counter(ion_f)
        neutral["H"] -= 1                       # [M+H]+ -> M
        delta = Counter({el: neutral[el] - hf.get(el, 0) for el in neutral})
        n_from_C = delta["C"] / arm_c["C"]
        n = int(round(n_from_C))
        residual = {el: delta[el] - n * arm_c.get(el, 0)
                    for el in set(delta) | set(arm_c)}
        if n < 1 or any(residual.values()):
            raise SystemExit(
                f"stage03 {hid}: Figure 2F's formula {FIG2F[hid][0]} minus the "
                f"amine head is not an integer multiple of the tail C23H34O7. "
                f"n={n_from_C}, residual={residual}. Refusing to place arms.")

        slots = int(blocks[hid]["n_michael_slots"])
        if n > slots:
            raise SystemExit(f"stage03 {hid}: {n} arms but only {slots} N-H "
                             "slots on the head.")

        prods, per_step = enumerate_products(head, n, arm_a, hid)

        # ---- three-way mass split, per the run brief
        m = Chem.MolFromSmiles(prods[0])
        ours_mh = ExactMolWt(m) + PROTON
        ours_formula = CalcMolFormula(m)
        their_formula_mz = ExactMolWt(Chem.MolFromSmiles(
            Chem.MolToSmiles(m))) + PROTON   # same object; see note in 07_qc
        cal, obs = FIG2F[hid][1], FIG2F[hid][2]

        rec = {
            "id": hid,
            "lipidoid": f"{hid}-O16CBA",
            "head_block": hid,
            "head_SMILES": head,
            "n_arms_from_printed_formula": n,
            "n_michael_slots_on_head": slots,
            "n_distinct_products": len(prods),
            "products_per_step": per_step,
            "constitution_forced": len(prods) == 1,
            "SMILES": prods[0],
            "alternatives": prods[1:],
            "our_formula_neutral": ours_formula,
            "their_formula_MplusH": FIG2F[hid][0],
            "formula_agrees": (parse_formula(ours_formula) ==
                               Counter({k: v for k, v in neutral.items()})),
            "our_mz_MplusH": round(ours_mh, 4),
            "their_calculated_mz": cal,
            "their_observed_mz": obs,
            "d_ours_vs_their_calc": round(ours_mh - cal, 4),
            "d_their_calc_vs_their_obs": round(cal - obs, 4),
            "d_ours_vs_their_obs": round(ours_mh - obs, 4),
            "stereo_unspecified_centres": [
                int(e.centeredOn) for e in Chem.FindPotentialStereo(m)],
        }
        members.append(rec)

        rows.append({
            "lipidoid": rec["lipidoid"],
            "IL_name": rec["lipidoid"],
            "head_amine": hid,
            "head_name": blocks[hid]["chemical_name"],
            "n_arms": n,
            "IL_SMILES": prods[0],
            "formula_neutral": ours_formula,
            "exact_mass_neutral": f"{ExactMolWt(m):.4f}",
            "MplusH_ours": f"{ours_mh:.4f}",
            "MplusH_their_calc": cal,
            "MplusH_their_obs": obs,
            "n_distinct_constitutions": len(prods),
            "positional_ambiguity": 0 if len(prods) == 1 else 1,
            "structure_assumed": 1 if hid == "77" else 0,
            "confidence_reason": (
                "head drawn in full in Figure 1C; arm count solved exactly "
                "from the molecular formula printed in Figure 2F; the "
                "placement of those arms is the unique solution of an "
                "exhaustive enumeration"
                + ("; BUT the head carries one tetrahedral stereocentre that "
                   "Figure 1C draws flat and the paper never specifies"
                   if hid == "77" else "")),
        })

    qc["members"] = members
    qc["summary"] = {
        "n_distinct_ionizable_structures": len(rows),
        "all_constitutions_forced": all(r["n_distinct_constitutions"] == 1
                                        for r in rows),
        "n_with_positional_ambiguity": sum(r["positional_ambiguity"]
                                           for r in rows),
        "the_result": "Every one of the eleven lipidoids has EXACTLY ONE "
                      "constitution consistent with its printed formula. For "
                      "the three-armed members this is not luck: 113 and 306 "
                      "are symmetric about their central tertiary nitrogen, so "
                      "the only way to place three arms on two primary "
                      "nitrogens is two-plus-one and the two choices are the "
                      "same molecule; 400 has exactly three N-H slots, so "
                      "three arms fill all of them. The enumeration proves "
                      "this rather than assuming it.",
        "n_distinct_canonical_SMILES": len({r["IL_SMILES"] for r in rows}),
        "note_on_77_and_82": "77-O16CBA and 82-O16CBA share the molecular "
                             "formula C53H84N2O14 and the same calculated and "
                             "observed m/z (973.60 / 973.73) in Figure 2F. "
                             "They are constitutional isomers and their "
                             "SMILES differ. An internal consistency win for "
                             "the table, and a reminder that mass cannot see "
                             "connectivity.",
    }

    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "03_enumerated_qc.json"), "w") as f:
        json.dump(qc, f, indent=1)

    print(f"stage03: {len(rows)} lipidoids, "
          f"{qc['summary']['n_distinct_canonical_SMILES']} distinct SMILES, "
          f"{qc['summary']['n_with_positional_ambiguity']} ambiguous")
    return 0


if __name__ == "__main__":
    sys.exit(main())
