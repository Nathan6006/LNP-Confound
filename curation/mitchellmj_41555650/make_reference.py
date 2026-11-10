#!/usr/bin/env python3
"""Build the C10 reference anchors for this curation -- INDEPENDENTLY.

INDEPENDENCE, and why it is argued this way.

Byte-identity of the emitted anchor is NOT the signal, and testing for it is a
trap (F-037): a correct re-derivation is byte-identical to a correct enumeration
once canonicalised, which is exactly why C10's own `anchor_may_be_a_copy` missed
a copied template. Independence here is a property of the ROUTE, and it is
CHECKED AT RUN TIME by walking this module's own AST (`audit_independence`
below) to confirm that:

  * it imports no assembler (`assemble`, `curation.lib.assemble`), and
  * it calls none of RWMol / AddBond / RemoveAtom / CombineMols / esterify /
    attach / AddAtom / ReplaceSubstructs.

The AST is walked rather than the file text substring-scanned, because a
substring scan flags the words in this very docstring -- a sibling's first
attempt did exactly that and reported its own prose as a violation.

THE ROUTE IS DELIBERATELY OPPOSITE TO THE ENUMERATOR'S.
  03_enumerate.py : starts from the AMINE CORE and attaches a pre-built arm
                    fragment to each N-H in turn, HEAD -> TAIL, with
                    assemble.attach doing every bond.
  this file       : writes the whole molecule as ONE SMILES string read off the
                    SI's systematic name, traversed TAIL -> HEAD, starting at
                    the far methyl of a lipid tail. It shares no fragment
                    string, no core string and no bond-forming code with the
                    enumerator. Nothing is concatenated: each string is a single
                    literal, written once, by hand.

A SMILES FRAGMENT HAS A DIRECTION MASS CANNOT SEE, so writing the chain
backwards is a real test of the enumerator's arm orientation and not a
formality: an arm attached through C2 instead of C1 has the same formula and
the same exact mass and is a different molecule.

TWO ANCHORS, DIFFERENT CHAIN LENGTHS AND DIFFERENT CORES, so a single
systematic slip cannot pass both:
  C10-306  C10 tails, 4 arms, acyclic triamine core
  C14-200  C14 tails, 5 arms, piperazine-containing pentamine core
"""
import ast, csv, json, os, sys
from rdkit import Chem
from rdkit.Chem import inchi
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

HERE = os.path.dirname(os.path.abspath(__file__))

FORBIDDEN_CALLS = {"RWMol", "AddBond", "RemoveAtom", "CombineMols",
                   "esterify", "attach", "AddAtom", "ReplaceSubstructs",
                   "ReactionFromSmarts", "RunReactants"}
FORBIDDEN_IMPORTS = {"assemble", "curation.lib.assemble", "curation.lib"}

# ---------------------------------------------------------------------------
# Hand-written, read off the SI systematic name, traversed TAIL -> HEAD.
#
# C10-306, SI p.3:
#   "1,1',1'',1'''-(((methylazanediyl)bis(propane-3,1-diyl))bis(azanetriyl))
#    tetrakis(decan-2-ol)"
#   -> CH3-N(-CH2CH2CH2-N(CH2CH(OH)C8H17)2)2
#   The string below starts at C10 of one decan-2-ol and walks inward to C2
#   (which carries the OH), then C1, then the nitrogen.
#
# C14-200, SI p.2:
#   "1,1'-((2-(4-(2-((2-(bis(2-hydroxytetradecyl)amino)ethyl)
#    (2-hydroxytetradecyl)amino)ethyl)piperazin-1-yl)ethyl)azanediyl)
#    bis(tetradecan-2-ol)"
#   -> five 2-hydroxytetradecyl arms on a
#      H2N-CH2CH2-[piperazine]-CH2CH2-NH-CH2CH2-NH2 skeleton.
# ---------------------------------------------------------------------------
ANCHORS = {
    "C10-306": dict(
        smiles=("CCCCCCCCC(O)CN(CC(O)CCCCCCCC)CCCN(C)CCCN"
                "(CC(O)CCCCCCCC)CC(O)CCCCCCCC"),
        source_figure="SI p.3, systematic name of C10-306 (no drawn structure)",
        note=("Written tail-first from the far methyl of one decan-2-ol. "
              "Four arms, methyl-bis(3-aminopropyl)amine core.")),
    "C14-200": dict(
        smiles=("CCCCCCCCCCCCC(O)CN(CC(O)CCCCCCCCCCCC)CCN1CCN(CC1)CCN"
                "(CC(O)CCCCCCCCCCCC)CCN(CC(O)CCCCCCCCCCCC)"
                "CC(O)CCCCCCCCCCCC"),
        source_figure="SI p.2, systematic name of C14-200 (no drawn structure)",
        note=("Written tail-first from the far methyl of one tetradecan-2-ol. "
              "Five arms; the piperazine is written as a ring closure in the "
              "middle of the traversal, where the enumerator reaches it last.")),
}

# A seeded defect with IDENTICAL FORMULA AND IDENTICAL EXACT MASS: the hydroxyl
# of ONE arm moved from C2 to C3. Mass cannot see it; connectivity can.
SEED = ("CCCCCCCCC(O)CN(CC(O)CCCCCCCC)CCCN(C)CCCN"
        "(CC(O)CCCCCCCC)CCC(O)CCCCCCC")


def audit_independence(path):
    """Walk THIS module's AST. Returns (ok, findings)."""
    tree = ast.parse(open(path).read())
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name.split(".")[0] in FORBIDDEN_IMPORTS or a.name in FORBIDDEN_IMPORTS:
                    bad.append(f"import {a.name}")
        elif isinstance(node, ast.ImportFrom):
            mod = node.module or ""
            if mod in FORBIDDEN_IMPORTS or mod.split(".")[0] in FORBIDDEN_IMPORTS:
                bad.append(f"from {mod} import ...")
        elif isinstance(node, ast.Call):
            f = node.func
            name = (f.attr if isinstance(f, ast.Attribute)
                    else f.id if isinstance(f, ast.Name) else None)
            if name in FORBIDDEN_CALLS:
                bad.append(f"call {name}() at line {node.lineno}")
    return (not bad), bad


def canon(s, label):
    m = Chem.MolFromSmiles(s)
    if m is None:
        raise SystemExit(f"{label}: anchor SMILES does not parse")
    if any(a.GetNumRadicalElectrons() for a in m.GetAtoms()):
        raise SystemExit(f"{label}: anchor carries a radical")
    return m


def main():
    ok, bad = audit_independence(os.path.abspath(__file__))
    enum = {r["IL_name"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}

    report = {
        "independence_checked_how":
            "ast.walk over this module at run time; see audit_independence()",
        "independent": ok, "findings": bad,
        "forbidden_calls": sorted(FORBIDDEN_CALLS),
        "forbidden_imports": sorted(FORBIDDEN_IMPORTS),
        "anchors": {},
    }
    out = []
    for name, d in ANCHORS.items():
        m = canon(d["smiles"], name)
        can = Chem.MolToSmiles(m)
        key = inchi.MolToInchiKey(m)
        e = enum[name]
        agree_skel = key.split("-")[0] == e["InChIKey_skeleton"]
        agree_full = key == e["InChIKey"]
        report["anchors"][name] = dict(
            anchor_formula=CalcMolFormula(m),
            anchor_exact_mass=round(ExactMolWt(m), 5),
            enumeration_formula=e["formula"],
            anchor_InChIKey=key, enumeration_InChIKey=e["InChIKey"],
            skeleton_match=agree_skel, full_match=agree_full,
            anchor_written_SMILES=d["smiles"],
            anchor_canonical=can,
            enumeration_SMILES=e["IL_SMILES"],
            written_string_is_byte_identical_to_enumeration=(d["smiles"] == e["IL_SMILES"]),
            canonical_is_byte_identical_to_enumeration=(can == e["IL_SMILES"]),
            byte_identity_note=(
                "Byte-identity of the CANONICAL form is EXPECTED when both are "
                "right and is NOT evidence of copying (F-037). The written "
                "string above is what was authored; compare THAT."))
        out.append(dict(
            lipid_name=name, reference_SMILES=can,
            source_figure=d["source_figure"],
            read_on="2026-11-10", read_by="curation agent (mitchellmj_41555650)",
            expected_mismatch="", note=d["note"]))

    with open(os.path.join(HERE, "reference_structures.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        for r in out:
            w.writerow([r["lipid_name"], r["reference_SMILES"], r["source_figure"],
                        r["read_on"], r["read_by"], r["expected_mismatch"],
                        r["note"]])
        # The remaining 13 lipids have NO anchor. Recorded as EMPTY rows so the
        # absence is visible rather than silent (CLAUDE.md section 5).
        for nm in sorted(set(enum) - set(ANCHORS)):
            w.writerow([nm, "", "no structure is drawn anywhere in the SI; the "
                        "main article was not acquired (ACS Cloudflare)",
                        "2026-11-10", "curation agent (mitchellmj_41555650)", "",
                        "NO ANCHOR. Structure rests on the SI systematic name "
                        "plus the stoichiometry check P1. Two of the fifteen "
                        "carry an independent anchor; this is not one."])

    json.dump(report, open(os.path.join(HERE, "reference_report.json"), "w"), indent=1)
    for n, a in report["anchors"].items():
        print(f"{n}: skeleton_match={a['skeleton_match']} full_match={a['full_match']} "
              f"{a['anchor_formula']} vs {a['enumeration_formula']}")
    print("independent:", ok, bad)
    bad_n = [n for n, a in report["anchors"].items() if not a["full_match"]]
    if bad_n or not ok:
        print("ANCHOR GATE FAILED:", bad_n)
        sys.exit(1)
    print("anchor gate: PASS")


if __name__ == "__main__":
    main()
