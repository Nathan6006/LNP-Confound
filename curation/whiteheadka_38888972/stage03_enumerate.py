#!/usr/bin/env python3
"""Enumerate the 580-member lipidoid library -> 03_enumerated.csv.

THE CHEMISTRY, AND THE THING MASS CANNOT SEE. Each lipidoid is the Michael
adduct of ONE amine head and ONE acrylate tail. The acrylate adds across an
N-H bond, so THE NUMBER OF TAILS PER MOLECULE IS THE HEAD'S N-H COUNT -- not a
constant, and not something a molecular formula or an exact mass can confirm
after the fact, because formula and mass are *consequences* of the count rather
than evidence for it. PIPELINE_V2: "Mass cannot see connectivity. Read the
synthesis."

WHAT THE SYNTHESIS SAYS. Methods, "Lipidoid library synthesis": "Each of the
ionizable lipidoids was synthesized by combining an alkyl-amine head with
alkyl-acrylate tails. Stoichiometric ratio of 1:4 was used. Reagents were
combined in a scintillation vial kept at 90 C for 3 days while stirring." And
Figure 1A DRAWS the reaction for 306 + 4 Oi10 -> 306Oi10, with all four tails
on the two primary amines of 306 -- a head whose N-H count is exactly 4.

So the paper pins the rule on a head where 4 equivalents and 4 N-H coincide,
and says nothing about the 14 heads where they do not. Measured here:

    N-H = 2   60, 80, 81, 602                              4 heads
    N-H = 3   123, 205, 303, 304, 601                      5 heads
    N-H = 4   113, 211, 305, 306, 313, 347, 400, 401,
              402, 500, 501, 503, 509, 514, 516           15 heads
    N-H = 5   174, 200                                     2 heads
    N-H = 6   110, 122                                     2 heads
    N-H = 7   111                                          1 head

FULL SUBSTITUTION IS WHAT THE CORPUS ALREADY HOLDS for this exact chemistry,
independently curated: AA_2008 ('87_O15' -> both N-H of head 87 esterified,
its two hydroxyls untouched), KL_2010 ('113C12' -> 4 tails on the 4 N-H of
CN(CCN)CCN) and KW_2014 ('A125_C10' -> 2 tails on a primary amine). So full
substitution is the enumerated structure here too.

BUT ON FIVE HEADS THAT IS NOT WHAT 4 EQUIVALENTS CAN DELIVER. 110, 111, 122,
174 and 200 have 5-7 N-H and only 4 equivalents of acrylate, so the product is
necessarily an under-substituted mixture and the stored structure is the
limiting case. Those 100 lipidoids carry structure_assumed = 1 and
confidence = low, and the question of what to store is DECISIONS.md D-2. It is
not resolved here.

NOTHING IS BUILT BY STRING CONCATENATION. Every bond is made by
curation/lib assemble.attach, which resolves the attachment atom by ATOM MAP
(never by index), hands the hydrogen count back to RDKit instead of
decrementing it, and asserts on every build that no radical appeared and that
the product's stereocentres are the fragments' sum.

Run: python3 -I stage03_enumerate.py
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

from rdkit import Chem                                          # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors            # noqa: E402

import assemble                                                 # noqa: E402
import emit                                                     # noqa: E402

# --- the three caveat populations, each stated as a condition on the RECORD --
# A head whose drawing in Figure 1B contradicts its Table S1 name. The name
# wins (these are catalogue reagents bought under it), but which half to
# believe is DECISIONS.md D-1 and is not settled here.
HEADS_NAME_VS_DRAWING = {"113", "200", "303"}
# A tail whose stored structure contains a feature neither the name nor the
# drawing pins.
TAILS_ASSUMED = {
    "Oi10": "isodecyl acrylate is sold as a branched-isomer mixture; stored "
            "as the 8-methylnonyl isomer drawn in Figure 1A/1C",
    "X8": "the name gives no E/Z for the 2-ene while every other X name states "
          "cis or trans; E is taken from the Figure 1C drawing, and C7/C11 "
          "are left unassigned",
}
RANK = {"high": 2, "medium": 1, "low": 0}


def load_blocks():
    heads, tails = {}, {}
    with open(os.path.join(HERE, "01_blocks.csv")) as fh:
        for r in csv.DictReader(fh):
            (heads if r["block_type"] == "amine_head" else tails)[r["code"]] = r
    if len(heads) != 29 or len(tails) != 20:
        raise SystemExit(f"blocks: {len(heads)} heads, {len(tails)} tails")
    return heads, tails


# A carbinol carbon flanked by two aminomethyl arms: the C2 of head 401,
# 1,3-diaminopropan-2-ol. See `prochiral_pivot` below.
PROCHIRAL_PIVOT = Chem.MolFromSmarts("[CX4;H1]([OX2H1])([CH2][NX3])[CH2][NX3]")


def prochiral_pivot(head_mol, label):
    """The one PROCHIRAL carbon that partial N-substitution desymmetrises.

    STATUS NOTE, 2026-11-10. `curation/lib/assemble.py` WAS EDITED MID-SESSION
    (mtime 00:35, by concurrent work on `andersondg_41203968`) to relax
    `check_stereo_preserved` so that it now refuses only a LOSS or INVERSION of
    an ASSIGNED centre, not the CREATION of a new unassignable one. Under that
    newer rule head 401 no longer refuses and this marker is not needed. It is
    kept because it is correct under either rule and because removing it would
    leave this stage silently dependent on which version of the library is on
    disk. VERIFIED: all 580 structures are byte-identical with the marker and
    without it (580 compared, 0 differences). Everything below describes the
    rule that was in force when this stage was written.

    WHY THIS EXISTS, AND WHY IT IS NOT A WEAKENED CHECK.
    `assemble.attach` asserts that a product's stereocentres are its fragments'
    sum. On head 401 (1,3-diaminopropan-2-ol) that assertion REFUSES, correctly
    and for a chemically real reason: C2 carries OH, H and two CH2-N arms, so
    it is prochiral in the free amine and in the tetra-substituted product, but
    a POTENTIAL stereocentre in between, when the two arms carry different
    numbers of tails. Measured: of the 29 heads, exactly one refuses, at tail 1
    (the centre appearing) and again at tail 4 (it disappearing).

    That is not the F-013 mechanism. F-013 is a centre going from ASSIGNED to
    something else through RWMol.RemoveAtom; here nothing is deleted, no R or S
    changes, and the count of unassignable centres moves transiently.

    The handling is a 13C MARKER on one of the two aminomethyl carbons, set
    before the build and stripped after. An isotope makes the two arms
    distinguishable to RDKit's canonical ranking, so the pivot is an
    unassignable ('?') centre in the head AND in every intermediate AND in the
    tetra-substituted product -- the count never moves and every intermediate
    comparison balances on its own terms. It invents no stereochemistry: '?'
    means "could be a centre, is not assigned", which is exactly true
    throughout. Nothing R or S is ever created.

    Then the marker is removed and `assemble.check_stereo_preserved` runs on
    the FINAL product against the REAL head and the REAL fragments. That
    assertion is the one that matters, it is not relaxed for anyone, and it
    runs on all 580 molecules whether or not they have a pivot.

    Returns (pivot_atom_idx, marker_atom_idx), or None. Raises if the head
    already carries an assigned stereocentre, because then a marker could
    change a real CIP label instead of only breaking a tie.
    """
    hits = head_mol.GetSubstructMatches(PROCHIRAL_PIVOT)
    if not hits:
        return None
    if len(hits) != 1:
        raise SystemExit(f"{label}: {len(hits)} prochiral pivots; refusing.")
    assigned = [c for c in assemble.stereocentres(head_mol).values()
                if c != "?"]
    if assigned:
        raise SystemExit(
            f"{label}: head carries assigned stereocentre(s) {assigned}; a "
            "marker could change a real CIP label. Refusing.")
    return hits[0][0], hits[0][2]


def build(head_smiles, adduct_smiles, n_nh, label):
    """Add one adduct fragment across every N-H of the head.

    The fragment's attachment carbon carries atom map 2 (set in
    stage01_blocks.py by esterifying a map-tagged propanoic acid). The head's
    reacting nitrogen is tagged with map 1 on the Mol, so attach never has to
    guess an index.
    """
    mol = Chem.MolFromSmiles(head_smiles)
    if mol is None:
        raise SystemExit(f"{label}: unparseable head")
    pv = prochiral_pivot(mol, label)
    if pv is not None:
        pivot, marker = pv
        if mol.GetAtomWithIdx(marker).GetIsotope():
            raise SystemExit(f"{label}: marker atom already isotopic.")
        mol.GetAtomWithIdx(marker).SetIsotope(13)
    added = 0
    while True:
        target = None
        for a in mol.GetAtoms():
            if a.GetSymbol() == "N" and a.GetTotalNumHs() > 0:
                target = a.GetIdx()
                break
        if target is None:
            break
        for a in mol.GetAtoms():
            a.SetAtomMapNum(0)
        mol.GetAtomWithIdx(target).SetAtomMapNum(1)
        smi = assemble.attach(mol, adduct_smiles, a_map=1, b_map=2,
                              label=f"{label} tail {added + 1}")
        mol = Chem.MolFromSmiles(smi)
        added += 1
        if added > 12:
            raise SystemExit(f"{label}: runaway substitution")
    if added != n_nh:
        raise SystemExit(
            f"{label}: added {added} tails but the head has {n_nh} N-H. The "
            "substitution count IS the structural claim; refusing.")
    for a in mol.GetAtoms():
        a.SetAtomMapNum(0)
    if pv is not None:
        # Strip the marker. Exactly one atom may carry it, and it must still be
        # the methylene next to the carbinol carbon -- asserted, not assumed,
        # because attach rebuilds the molecule from SMILES on every step.
        marked = [a for a in mol.GetAtoms() if a.GetIsotope() == 13]
        if len(marked) != 1:
            raise SystemExit(f"{label}: {len(marked)} marker atoms at the end.")
        a = marked[0]
        if not (a.GetSymbol() == "C" and a.GetTotalNumHs() == 2):
            raise SystemExit(f"{label}: marker moved off its methylene.")
        a.SetIsotope(0)
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
    # THE ASSERTION THAT MATTERS, and it runs on every molecule, pivot or not:
    # the finished lipidoid's stereocentres must be the real head's plus
    # n_nh copies of the real tail fragment. This is not relaxed for 401.
    assemble.check_stereo_preserved(
        mol, [head_smiles] + [adduct_smiles] * n_nh, label=label)
    assemble.check_no_radicals(mol, label=label)
    return Chem.MolToSmiles(mol)


def main():
    heads, tails = load_blocks()
    order_h = [r["code"] for r in csv.DictReader(
        open(os.path.join(HERE, "01_blocks.csv")))
        if r["block_type"] == "amine_head"]
    order_t = [r["code"] for r in csv.DictReader(
        open(os.path.join(HERE, "01_blocks.csv")))
        if r["block_type"] == "acrylate_tail"]

    records, items = {}, []
    for h in order_h:
        H = heads[h]
        n_nh = int(H["n_NH_bonds"])
        for t in order_t:
            T = tails[t]
            name = f"{h}{t}"
            smi = build(H["SMILES"], T["adduct_fragment_SMILES"], n_nh, name)
            m = Chem.MolFromSmiles(smi)

            # --- post-build arithmetic, on the product, not on a comment ----
            n_ester = len(m.GetSubstructMatches(
                Chem.MolFromSmarts("[CX3](=O)[OX2][CX4]")))
            n_nh_left = sum(a.GetTotalNumHs() for a in m.GetAtoms()
                            if a.GetSymbol() == "N")
            if n_ester != n_nh:
                raise SystemExit(f"{name}: {n_ester} esters for {n_nh} tails")
            if n_nh_left != 0:
                raise SystemExit(f"{name}: {n_nh_left} N-H survived")

            reasons, conf = [], "high"
            if n_nh > 4:
                reasons.append(
                    f"head {h} has {n_nh} N-H but the synthesis used 4 "
                    "equivalents of acrylate, so full substitution is a "
                    "limiting case, not a product")
                conf = "low"
            if h in HEADS_NAME_VS_DRAWING:
                reasons.append(
                    f"Figure 1B draws head {h} differently from its Table S1 "
                    "name; the catalogue name is stored (08_contradictions.csv)")
                conf = "low"
            if t in TAILS_ASSUMED:
                reasons.append(TAILS_ASSUMED[t])
                if RANK[conf] > RANK["medium"]:
                    conf = "medium"
            note = ""
            if reasons:
                note = ("Whitehead 2024 lipidoid " + name + ": "
                        + "; ".join(reasons) + ".")
            elif h == "516":
                note = ("The but-2-ene bridge of head 516 is drawn with a WAVY "
                        "bond in Figure 1B and its Table S1 name gives no E/Z, "
                        "so the geometry is stored unassigned, as declared.")

            rec = {
                "IL_name": name, "head": h, "tail": t,
                "head_name": H["name_verbatim"], "tail_name": T["name_verbatim"],
                "head_SMILES": H["SMILES"],
                "tail_acrylate_SMILES": T["SMILES"],
                "n_tails": n_nh, "n_NH_in_head": n_nh,
                "IL_SMILES": smi,
                "formula": rdMolDescriptors.CalcMolFormula(m),
                "exact_mass": round(Descriptors.ExactMolWt(m), 4),
                "mol_weight": round(Descriptors.MolWt(m), 3),
                "n_ester": n_ester,
                "inchikey": Chem.MolToInchiKey(m),
                "confidence_basis": conf,
                "note_basis": note,
                "reasons": " | ".join(reasons),
            }
            records[name] = rec
            items.append(name)

    # --- flags on the row, through the one writer that refuses a literal ----
    flags = emit.emit_rows(
        items, records,
        confidence=lambda i, r: r["confidence_basis"],
        note=lambda i, r: r["note_basis"],
        structure_assumed=lambda i, r: "1" if r["reasons"] else "0",
        positional_ambiguity=lambda i, r: "0",
        label="whiteheadka_38888972 enumeration")

    cols = ["IL_name", "head", "tail", "head_name", "tail_name", "head_SMILES",
            "tail_acrylate_SMILES", "n_NH_in_head", "n_tails", "n_ester",
            "IL_SMILES", "formula", "exact_mass", "mol_weight", "inchikey",
            "confidence", "structure_assumed", "positional_ambiguity",
            "curation_note", "reasons"]
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for name, fl in zip(items, flags):
            row = dict(records[name])
            row.update(fl)
            w.writerow(row)

    keys = [records[i]["inchikey"] for i in items]
    qc = {
        "n_enumerated": len(items),
        "stated_library_size": 580,
        "difference": len(items) - 580,
        "n_unique_inchikey": len(set(keys)),
        "n_unique_canonical_smiles": len({records[i]["IL_SMILES"]
                                          for i in items}),
        "tails_per_head": {h: int(heads[h]["n_NH_bonds"]) for h in order_h},
        "rows_with_a_caveat": emit.caveat_count(flags),
        "rows_confidence_low": sum(1 for f in flags if f["confidence"] == "low"),
        "rows_confidence_medium": sum(1 for f in flags
                                      if f["confidence"] == "medium"),
        "rows_confidence_high": sum(1 for f in flags
                                    if f["confidence"] == "high"),
    }
    if qc["n_unique_inchikey"] != len(items):
        qc["WARNING_duplicate_structures"] = True
    with open(os.path.join(HERE, "03_enumerated_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1, sort_keys=True)
    print(json.dumps(qc, indent=1, sort_keys=True))


if __name__ == "__main__":
    main()
