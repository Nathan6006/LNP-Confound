#!/usr/bin/env python3
"""Stage 03 -- the 384-member SD_2026 library, enumerated by the Ugi 4CR.

THE REACTION IS THE PAPER'S OWN CLAIM, NOT OURS. Results, first line: "We
generated a library of 384 ionizable lipids using a four-component Ugi reaction
(4CR), combining 12 headgroups, four tails, two subtails and four linkers."
Supplementary Fig. 1 repeats the arithmetic: 12 x 4 x 2 x 4 = 384.

WHY A REACTION SMARTS AND NOT A CHAIN OF attach() CALLS. The Ugi expels one
water, and the atom it expels is the ALDEHYDE's oxygen -- the acid keeps both
of its own (one becomes the acyl C=O on the nitrogen, the other becomes the
C=O on the isocyanide-derived carbon, by the Mumm rearrangement). curation/lib/
assemble exposes `attach` (no atom leaves) and `esterify` (acid + alcohol), and
neither expresses that. Writing the deletion by hand is precisely the F-013
hazard, so the condensation is instead stated once, declaratively, as a mapped
reaction SMARTS that a chemist can read in one line -- and EVERY product is
then put through the library's own assertions:

    assemble.check_no_radicals      -- the Siegwart radical-carbon class
    assemble.check_stereo_preserved -- F-013; H8 and H11 carry a ring
                                       stereocentre that must survive intact

`check_stereo_preserved` permits a centre to be CREATED, which matters here:
L2, L3 and L4 make the former aldehyde carbon a new stereocentre in 288 of the
384 products. The Ugi is not stereoselective and the paper reports no
resolution, so those are emitted without stereochemistry, and that is recorded
on the row rather than left implicit.

WHAT CHECKS THIS. Mass cannot see connectivity, and a formula check is
necessary and not sufficient -- so it is not the only check here:
  * the SI prints a molecular formula AND an MW for each of the ten lipids the
    paper carried forward (SC-lipid-1..5, LC-lipid-1..5). Every one of them
    must be realised by some member of this enumeration, exactly.
  * Supplementary Table 5 names LC-1's four reagents by CAS number. That fixes
    ONE member completely -- H4 + L4 + Ta1 + Tb1 -- so its formula, its MW and
    its structure are all pinned by something other than a drawing.
  * the drawn structure of LC-1 in Fig. 1l is the C10 anchor, built by
    make_reference.py through a different code path (no reaction SMARTS).
"""
import csv
import json
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import assemble                                                  # noqa: E402
import emit                                                      # noqa: E402
from errors import StructureError                                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# The Ugi four-component condensation, mapped atom by atom.
#   aldehyde  R2-CHO        its oxygen is the one that leaves
#   amine     R1-NH2        becomes the tertiary amide nitrogen
#   acid      R3-COOH       BOTH oxygens survive (Mumm rearrangement)
#   isocyanide R4-NC        its carbon becomes the second amide carbonyl
UGI = ("[C;H1,H2;$(C=O):1]=[O;X1]"
       ".[N;H2;X3:2]"
       ".[C;X3:3](=[O;X1:4])[O;X2;H1:5]"
       ".[C-;X1:6]#[N+;X2:7]"
       ">>[C:3](=[O:4])[N;H0:2][C:1][C+0:6](=[O:5])[N+0;H1:7]")

# What the SI prints for the ten lipids the paper carried forward. Page numbers
# are the SI PDF's own, 1-based. These are transcribed by hand ON PURPOSE: the
# point of the check is that a number we did not compute agrees with one we did.
SI_NAMED = {
    "SC-lipid-1": ("C48H93N3O2",  744.29, 5),
    "SC-lipid-2": ("C50H98N4O4",  819.36, 8),
    "SC-lipid-3": ("C50H97N3O5",  820.34, 11),
    "SC-lipid-4": ("C48H93N3O5",  792.29, 14),
    "SC-lipid-5": ("C48H95N3O4",  778.31, 17),
    "LC-lipid-1": ("C50H97N3O4",  804.34, 20),
    "LC-lipid-2": ("C50H99N3O4",  806.36, 23),
    "LC-lipid-3": ("C54H105N3O4", 860.45, 26),
    "LC-lipid-4": ("C52H101N3O4", 832.40, 29),
    "LC-lipid-5": ("C54H105N3O4", 860.45, 32),
}

# Supplementary Table 5 names LC-1's four reagents by CAS. This is the one
# member of the 384 whose block assignment is stated rather than inferred.
LC1_BLOCKS = ("H4", "L4", "Ta1", "Tb1")


def blocks():
    out = {}
    with open(os.path.join(HERE, "01_blocks.csv")) as fh:
        for r in csv.DictReader(fh):
            out.setdefault(r["block_class"], {})[r["block_id"]] = r
    return out


def ugi(aldehyde, amine, acid, isocyanide, label=""):
    """One Ugi product, with the library's assertions run on it."""
    rxn = AllChem.ReactionFromSmarts(UGI)
    mols = tuple(Chem.MolFromSmiles(s)
                 for s in (aldehyde, amine, acid, isocyanide))
    prods = rxn.RunReactants(mols)
    smis = set()
    for p in prods:
        m = p[0]
        try:
            Chem.SanitizeMol(m)
        except Exception as e:
            raise StructureError(f"stage03 {label}: sanitize failed: {e}")
        smis.add(Chem.MolToSmiles(m))
    if len(smis) != 1:
        raise StructureError(
            f"stage03 {label}: the Ugi template matched {len(prods)} ways and "
            f"gave {len(smis)} DISTINCT products: {sorted(smis)}. A block with "
            "two reactive sites of the same kind makes the product ambiguous; "
            "that is a positional_ambiguity to flag, not a set to pick from.")
    prod = smis.pop()
    assemble.check_no_radicals(prod, label=label)
    stereo = assemble.check_stereo_preserved(
        prod, [aldehyde, amine, acid, isocyanide], label=label)
    return prod, stereo


def main():
    b = blocks()
    H = b["headgroup_amine"]
    L = b["linker_aldehyde"]
    TA = b["tail_a_acid"]
    TB = b["tail_b_isocyanide"]

    qc = {"_what": "Stage 03: Ugi 4CR enumeration of the SD_2026 library.",
          "reaction_smarts": UGI,
          "block_counts": {"headgroups": len(H), "linkers": len(L),
                           "tails_a": len(TA), "tails_b": len(TB)},
          "product_of_counts": len(H) * len(L) * len(TA) * len(TB)}

    rows, seen = [], {}
    created_centres = 0
    for h in sorted(H, key=lambda k: int(k[1:])):
        for ta in sorted(TA, key=lambda k: int(k[2:])):
            for tb in sorted(TB, key=lambda k: int(k[2:])):
                for l in sorted(L, key=lambda k: int(k[1:])):
                    name = f"{h}-{l}-{ta}-{tb}"
                    smi, st = ugi(L[l]["SMILES"], H[h]["SMILES"],
                                  TA[ta]["SMILES"], TB[tb]["SMILES"],
                                  label=name)
                    m = Chem.MolFromSmiles(smi)
                    if st.get("centres_created"):
                        created_centres += 1
                    rows.append({
                        "lipid_id": name, "headgroup": h, "linker": l,
                        "tail_a": ta, "tail_b": tb,
                        "canonical_SMILES": smi,
                        "formula": rdMolDescriptors.CalcMolFormula(m),
                        "MolWt": f"{Descriptors.MolWt(m):.2f}",
                        "ExactMolWt": f"{Descriptors.ExactMolWt(m):.4f}",
                        "unassigned_stereocentres": str(len(
                            Chem.FindMolChiralCenters(
                                m, includeUnassigned=True,
                                useLegacyImplementation=False))),
                        "stereo_fingerprint": json.dumps(st["product"],
                                                         sort_keys=True)})
                    seen.setdefault(smi, []).append(name)

    qc["n_enumerated"] = len(rows)
    qc["n_distinct_canonical_smiles"] = len(seen)
    qc["duplicate_structures"] = {k: v for k, v in seen.items() if len(v) > 1}
    qc["n_products_with_a_new_stereocentre"] = created_centres

    # --- reconcile against the stated 384 ----------------------------------
    stated = 384
    qc["stated_library_size"] = {
        "abstract": 384, "results_first_line": 384,
        "supplementary_fig_1_caption": "12 x 4 x 2 x 4 = 384",
        "introduction": "more than 300",
        "reconciliation":
            f"{len(rows)} enumerated against {stated} stated, from "
            f"{len(H)} x {len(L)} x {len(TA)} x {len(TB)}. The Introduction's "
            "'a combinatorial library of more than 300 ionizable lipids' is "
            "the paper's own looser restatement of the same number and is "
            "recorded in 08_contradictions.csv, not treated as a second count."}
    if len(rows) != stated:
        raise StructureError(
            f"stage03: enumerated {len(rows)}, the paper states {stated}.")

    # --- the SI's ten printed formulas must each be realised ---------------
    by_formula = {}
    for r in rows:
        by_formula.setdefault(r["formula"], []).append(r["lipid_id"])
    named = {}
    for nm, (formula, mw, page) in SI_NAMED.items():
        hits = by_formula.get(formula, [])
        ref = Chem.MolFromSmiles("C")   # placeholder, replaced below
        named[nm] = {"si_formula": formula, "si_mw": mw, "si_page": page,
                     "n_library_members_with_this_formula": len(hits),
                     "members": hits}
        if not hits:
            raise StructureError(
                f"stage03: the SI prints {formula} for {nm} (SI p. {page}) and "
                "NO member of the enumerated library has that formula. A block "
                "is wrong.")
        got = float(rows[[r["lipid_id"] for r in rows].index(hits[0])]["MolWt"])
        named[nm]["our_mw_for_that_formula"] = round(got, 2)
        named[nm]["mw_delta"] = round(got - mw, 3)
        if abs(got - mw) > 0.06:
            raise StructureError(
                f"stage03: {nm} formula {formula} gives MW {got:.2f}; the SI "
                f"prints {mw}.")
    qc["si_named_lipids"] = named
    qc["si_named_lipids_all_realised"] = True
    qc["si_formula_check_is_necessary_not_sufficient"] = (
        "Ten printed formulas all land on a library member, and the MW we "
        "compute from each agrees with the printed MW to <= 0.06. That tests "
        "the BLOCK INVENTORY and the reaction's atom balance. It does NOT test "
        "connectivity -- several members share a formula, and mass is blind to "
        "which nitrogen an arm sits on. Connectivity is tested by the Fig. 1l "
        "anchor in reference_structures.csv (C10).")

    # --- the one member Supplementary Table 5 pins outright -----------------
    lc1_id = "{}-{}-{}-{}".format(LC1_BLOCKS[0], LC1_BLOCKS[1],
                                  LC1_BLOCKS[2], LC1_BLOCKS[3])
    lc1 = [r for r in rows if r["lipid_id"] == lc1_id]
    if not lc1:
        raise StructureError(f"stage03: {lc1_id} is not in the enumeration.")
    lc1 = lc1[0]
    si_f, si_mw, _ = SI_NAMED["LC-lipid-1"]
    qc["table5_pinned_member"] = {
        "lipid_id": lc1_id,
        "reagents_named_by_cas_in_table5": {
            "amine": "3-(pyrrolidin-1-yl)propan-1-amine, CAS 23159-07-1 -> H4",
            "aldehyde": "pivaldehyde, CAS 630-19-3 -> L4",
            "acid": "myristic acid, CAS 544-63-8 -> Ta1",
            "isocyanide": "6-isocyanohexyl 2-hexylundecanoate, MW 379.63 -> Tb1"},
        "our_formula": lc1["formula"], "si_formula": si_f,
        "our_mw": lc1["MolWt"], "si_mw": si_mw,
        "formula_matches": lc1["formula"] == si_f,
        "mw_delta": round(float(lc1["MolWt"]) - si_mw, 3),
        "canonical_SMILES": lc1["canonical_SMILES"]}
    if lc1["formula"] != si_f:
        raise StructureError(
            f"stage03: Supplementary Table 5 builds LC-1 from {LC1_BLOCKS}, "
            f"which this enumeration gives as {lc1['formula']}, but the SI "
            f"prints {si_f} for LC-lipid-1.")

    # --- flags, through the ONE writer -------------------------------------
    recs = {r["lipid_id"]: r for r in rows}
    ring_stereo_heads = {"H8", "H11"}

    def confidence(lid, rec):
        # Derived from the record, never a literal. The pinned member is the
        # only one whose block assignment is stated rather than inferred; it is
        # still not `high`, because its BLOCKS come from a schematic figure.
        return "medium" if lid == lc1_id else "low"

    def note(lid, rec):
        bits = []
        if lid == lc1_id:
            bits.append("the one library member whose four reagents "
                        "Supplementary Table 5 names by CAS number")
        else:
            bits.append("block assignment is the combinatorial enumeration the "
                        "paper describes; this specific member is not named, "
                        "drawn or measured anywhere in the paper")
        if rec["linker"] != "L1":
            bits.append("the Ugi creates a stereocentre at the linker-derived "
                        "carbon and the reaction is not stereoselective, so "
                        "this is stored without stereochemistry")
        if rec["headgroup"] in ring_stereo_heads:
            bits.append(f"{rec['headgroup']} is drawn with no wedge at its "
                        "ring stereocentre")
        return "; ".join(bits) + "."

    flags = emit.emit_rows(
        [r["lipid_id"] for r in rows], recs,
        confidence=confidence, note=note,
        structure_assumed=lambda lid, rec: "0" if lid == lc1_id else "1",
        positional_ambiguity=lambda lid, rec: "0",
        value_normalization_unit=lambda lid, rec: "",
        label="sd_2026 enumeration")
    for r, f in zip(rows, flags):
        r.update(f)
    qc["rows_with_a_caveat"] = emit.caveat_count(flags)

    cols = (["lipid_id", "headgroup", "linker", "tail_a", "tail_b",
             "canonical_SMILES", "formula", "MolWt", "ExactMolWt",
             "unassigned_stereocentres", "stereo_fingerprint"]
            + emit.columns())
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    with open(os.path.join(HERE, "03_enumerated_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)

    print(f"03_enumerated.csv: {len(rows)} members, "
          f"{len(seen)} distinct structures, "
          f"{len(qc['duplicate_structures'])} duplicate groups")
    print("LC-1 (H4-L4-Ta1-Tb1):", lc1["formula"], lc1["MolWt"],
          "vs SI", si_f, si_mw)
    print("SI named formulas realised:",
          {k: v["n_library_members_with_this_formula"] for k, v in named.items()})


if __name__ == "__main__":
    main()
