#!/usr/bin/env python3
"""The C10 anchors for KW_2019 -- derived from the paper's CHEMICAL NAMES.

THE PAPER DRAWS NOTHING. All 17 SI pages were enumerated: there is no reaction
scheme, no structure panel and no drawn lipidoid anywhere, and the main article
is behind a subscription (00_value_sources.json). So the usual C10 anchor --
"the structure the paper DRAWS" -- does not exist for this paper, and storing
an empty reference_SMILES would make the gap visible but test nothing.

What the paper DOES state, in words, on SI p.2, is both reagents by chemical
name: the amine "3,3-Diamino-N-methyldipropylamine (306, Sigma Aldrich)" and
"alkyl-acrylate tails of varying lengths". A chemical name is external ground
truth of a different kind from a drawing and it is what pins connectivity here.
So each anchor below is a SOURCE STRING LITERAL, typed by hand from that name
plus the alkyl group the "O<n>" token denotes, written HEAD-FIRST in a single
traversal starting at the N-methyl nitrogen.

HOW INDEPENDENCE IS ESTABLISHED -- BY PROVENANCE, NOT BY BYTE DIFFERENCE.
F-037 is the precedent and it is a trap: a correct hand-derivation is
byte-identical to a correct enumeration once canonicalised, so testing for byte
difference rewards a wrong anchor and punishes a right one. What makes this
module independent is checkable by reading it:

  * it imports NO assembler -- no `assemble`, no `RWMol`, `AddBond`,
    `RemoveAtom`, `CombineMols`, `EditableMol` or `ReactionFromSmarts`;
  * it reads NOTHING from 03_enumerated.csv except, at comparison time, the
    strings it is checking AGAINST;
  * the anchors are literals in this file, typed from the chemical names, by a
    traversal order (core -> arms) opposite to stage01's (arm -> core);
  * `assert_no_assembler()` below enforces the first bullet by inspecting this
    module's own source and import table at run time, so the claim is
    EXECUTED, not asserted in a comment. CLAUDE.md: "'A different code path'
    has to be executed, not asserted."

THE GATE. `python3 -I make_reference.py` exits 0 only if every anchor matches
its enumerated counterpart by FULL InChIKey. `--seed <kind>` injects a defect
and must exit 1; 07_qc.json records that it does.
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# --- the anchors ------------------------------------------------------------
# Typed from: 3,3'-diamino-N-methyldipropylamine, both primary nitrogens fully
# Michael-added with two equivalents of the named alkyl acrylate, giving
# N-CH2CH2-C(=O)O-alkyl. Traversal: N-methyl core first, then each arm.
#
#   decyl        = n-C10H21            -> OCCCCCCCCCC
#   isodecyl     = 8-methylnonyl       -> OCCCCCCCC(C)C
#   octadecyl    = n-C18H37            -> OCCCCCCCCCCCCCCCCCC
#
# Three anchors on purpose: one linear mid-series, the one branched tail (the
# paper's headline compound), and the longest tail in the set. A single anchor
# at one chain length cannot detect an off-by-one that scales with n.
ANCHORS = {
    "306O10": (
        "CN(CCCN(CCC(=O)OCCCCCCCCCC)CCC(=O)OCCCCCCCCCC)"
        "CCCN(CCC(=O)OCCCCCCCCCC)CCC(=O)OCCCCCCCCCC",
        "SI p.2 Methods, reagent names: 3,3-Diamino-N-methyldipropylamine "
        "(306) + decyl acrylate; tetra-substituted at 1:4 amine:acrylate"),
    "306Oi10": (
        "CN(CCCN(CCC(=O)OCCCCCCCC(C)C)CCC(=O)OCCCCCCCC(C)C)"
        "CCCN(CCC(=O)OCCCCCCCC(C)C)CCC(=O)OCCCCCCCC(C)C",
        "SI p.2 Methods, reagent names: 306 + ISODECYL acrylate, written as "
        "8-methylnonyl; the specific branch is a convention, see D-4"),
    "306O18": (
        "CN(CCCN(CCC(=O)OCCCCCCCCCCCCCCCCCC)CCC(=O)OCCCCCCCCCCCCCCCCCC)"
        "CCCN(CCC(=O)OCCCCCCCCCCCCCCCCCC)CCC(=O)OCCCCCCCCCCCCCCCCCC",
        "SI p.2 Methods, reagent names: 306 + octadecyl acrylate"),
}

# Module names and RDKit editing entry points that an ASSEMBLER would use.
# These are the thing the anchor must not share with the enumerator.
FORBIDDEN_MODULES = ("assemble",)
FORBIDDEN_CALLS = ("RWMol", "AddBond", "RemoveAtom", "CombineMols",
                   "EditableMol", "ReactionFromSmarts", "RunReactants",
                   "esterify", "attach")


def assert_no_assembler():
    """Execute the independence claim instead of writing it in a comment.

    Inspected by AST, not by substring. A substring scan of a file that
    contains its own blocklist can only false-positive, and the first version
    of this function did exactly that: it reported all eight forbidden tokens
    "found in the module body" on a clean run, because the tuple listing them
    is in the module body. A check that cannot pass is as useless as one that
    cannot fail -- F-038's three marginal checks scored a literal 0 for twelve
    days while the summary beside them said "ALL FOUR EXACT". The AST sees the
    blocklist as string Constants and an assembler call as a Name or
    Attribute, which is the distinction that matters.
    """
    import ast
    src = open(os.path.abspath(__file__)).read()
    tree = ast.parse(src)
    imported, called = set(), set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                imported.add(a.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module.split(".")[0])
        elif isinstance(node, ast.Call):
            f = node.func
            if isinstance(f, ast.Name):
                called.add(f.id)
            elif isinstance(f, ast.Attribute):
                called.add(f.attr)
    bad_mod = sorted(imported & set(FORBIDDEN_MODULES))
    bad_call = sorted(called & set(FORBIDDEN_CALLS))
    live = sorted(m for m in sys.modules
                  if m.split(".")[0] in FORBIDDEN_MODULES)
    if bad_mod or bad_call or live:
        raise SystemExit(
            f"make_reference: assembler imports {bad_mod}, assembler calls "
            f"{bad_call}, loaded assembler modules {live}. The anchor must not "
            "share a code path with the enumerator; that is F-013 (a byte "
            "copy) and F-037 (the file written to fix F-013, which re-copied "
            "the other half).")
    return {"ast_imports": sorted(imported),
            "forbidden_modules_imported": bad_mod,
            "forbidden_assembler_calls_used": bad_call,
            "assembler_modules_loaded_at_runtime": live,
            "how": "ast.walk over this module's own source, not a substring scan"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", choices=["swap_tail_length", "move_an_ester",
                                       "drop_an_arm"], default=None,
                    help="inject a defect; the gate must then exit 1")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    prov = assert_no_assembler()
    from rdkit import Chem, RDLogger
    from rdkit.Chem import inchi
    RDLogger.DisableLog("rdApp.*")

    anchors = {k: list(v) for k, v in ANCHORS.items()}
    if args.seed == "swap_tail_length":
        # one CH2 removed from every decyl arm: C58 instead of C59
        anchors["306O10"][0] = anchors["306O10"][0].replace(
            "OCCCCCCCCCC", "OCCCCCCCCC")
    elif args.seed == "move_an_ester":
        # the ester oxygen moved one carbon along ONE arm. Same formula, same
        # exact mass to every decimal place -- invisible to any mass check.
        anchors["306O10"][0] = anchors["306O10"][0].replace(
            "CCC(=O)OCCCCCCCCCC", "CCC(=O)OCC(C)CCCCCCC", 1)
    elif args.seed == "drop_an_arm":
        anchors["306O18"][0] = anchors["306O18"][0].replace(
            "CCCN(CCC(=O)OCCCCCCCCCCCCCCCCCC)CCC(=O)OCCCCCCCCCCCCCCCCCC",
            "CCCNCCC(=O)OCCCCCCCCCCCCCCCCCC", 1)

    with open(os.path.join(HERE, "03_enumerated.csv"), newline="") as fh:
        enum = {r["paper_name"]: r["IL_SMILES"] for r in csv.DictReader(fh)}

    out, failures = [], []
    for name, (smi, src) in anchors.items():
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise SystemExit(f"{name}: anchor does not parse")
        a_key = inchi.MolToInchiKey(m)
        e_smi = enum[name]
        e_key = inchi.MolToInchiKey(Chem.MolFromSmiles(e_smi))
        ok = a_key == e_key
        skel_ok = a_key.split("-")[0] == e_key.split("-")[0]
        if not ok:
            failures.append((name, a_key, e_key))
        out.append({
            "lipid_name": name,
            "reference_SMILES": Chem.MolToSmiles(m),
            "source_figure": src,
            "read_on": "2026-11-10",
            "read_by": "curation agent, from the SI Methods reagent NAMES "
                       "(the paper draws no structure anywhere)",
            "expected_mismatch": "",
            "note": (f"Anchor typed as a source string literal, head-first "
                     f"traversal, no assembler imported (verified at run time "
                     f"by assert_no_assembler). full InChIKey "
                     f"{'MATCH' if ok else 'MISMATCH'}; skeleton "
                     f"{'match' if skel_ok else 'MISMATCH'}. Byte-identity to "
                     f"the enumeration is NOT the independence signal (F-037); "
                     f"provenance is."),
        })

    if not args.seed:
        with open(os.path.join(HERE, "reference_structures.csv"), "w",
                  newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=[
                "lipid_name", "reference_SMILES", "source_figure", "read_on",
                "read_by", "expected_mismatch", "note"])
            w.writeheader(); w.writerows(out)

    if not args.quiet:
        print(f"independence: {prov}")
        for r in out:
            print(f"  {r['lipid_name']:9s} "
                  f"{'MATCH' if r['lipid_name'] not in [f[0] for f in failures] else 'MISMATCH'}")
    if failures:
        for name, a, e in failures:
            print(f"ANCHOR MISMATCH {name}: anchor {a} != enumerated {e}",
                  file=sys.stderr)
        print(f"{len(failures)}/{len(out)} anchors disagree with the "
              f"enumeration.", file=sys.stderr)
        return 1
    print(f"{len(out)}/{len(out)} anchors match by FULL InChIKey.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
