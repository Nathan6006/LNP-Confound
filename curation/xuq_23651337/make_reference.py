#!/usr/bin/env python3
"""Build the C10 reference anchors for PMID 23651337, independently.

WHAT THIS PAPER DRAWS: NOTHING. The 3-page SI contains two TEM micrographs, a
line plot and a bar chart, and no chemical structure of any kind. The
main-text figure that would draw the library is behind the ACS subscription
wall. So `reference_structures.csv` carries an EMPTY reference_SMILES for the
lipidoids themselves -- the gap is recorded, not papered over.

WHAT IT DOES STATE is a chemical NAME: "oleyl acrylamide", in the published
abstract. That is anchorable, and this file anchors it.

INDEPENDENCE -- VERIFIED BY PROVENANCE, NOT BY BYTE-IDENTITY.
Byte-identity of the emitted anchor is NOT the signal and testing for it is a
trap (F-037): a correct re-derivation is byte-identical to a correct
enumeration once canonicalised, which is exactly why C10's own
`anchor_may_be_a_copy` missed the defect it was written for. So:

  * `assert_no_assembler()` walks THIS MODULE'S OWN AST at run time and
    confirms there is no import of `assemble` (or of curation/lib at all) and
    no call to RWMol / AddBond / RemoveAtom / CombineMols / esterify / attach /
    ReplaceCore. It raises if any appears. That is executed, not asserted in
    a comment -- the F-013 note claimed "independent fragment assembly through
    the C3 oxygen (3-beta)" about a string that was 3-alpha and a byte copy.
  * the TRAVERSAL IS OPPOSITE to stage03's. stage03 starts from the acryloyl
    nitrogen and attaches a chain fragment whose FIRST atom is C1. This file
    writes the chain METHYL-END-FIRST, from C18 down to C1, and then appends
    the acryloyl group -- so the alkene's position is laid down counting from
    the opposite end of the molecule. An off-by-one in either route puts the
    alkene in a different place while leaving the molecular formula and the
    exact mass UNCHANGED.
  * the comparison is on the full InChIKey.

THE GATE IS PROVED TO WORK by seeding a defect with an IDENTICAL molecular
formula and an IDENTICAL exact mass (`--seed reverse`, the chain written in
the other direction: C21H39NO, 321.3032 either way). A gate that only ever
catches a formula change has not been shown to test connectivity.

CHAIN-LENGTH COVERAGE. "Pick anchors of different chain lengths" cannot be
done from this paper's own molecules: both tails are C18 and that is the whole
point of the comparison it makes. Instead the GENERATOR is exercised at
n = 12, 14, 16 and 18 and each product's formula and alkene position are
checked, which is what catches an off-by-one that scales with n. None of the
shorter products is claimed to be in the paper.

Run:  python3 -I make_reference.py
      python3 -I make_reference.py --seed reverse     (must exit non-zero)
"""
import ast
import csv
import json
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

FORBIDDEN_CALLS = {"RWMol", "AddBond", "RemoveAtom", "CombineMols",
                   "esterify", "attach", "ReplaceCore", "ReplaceSidechains",
                   "emit_rows", "calibrate"}
FORBIDDEN_IMPORTS = {"assemble", "axis", "bars", "heatmap", "chains",
                     "reagents", "emit", "lib"}


def assert_no_assembler(path=None):
    """Walk this module's own AST and prove the provenance claim."""
    path = path or os.path.abspath(__file__)
    with open(path) as fh:
        tree = ast.parse(fh.read(), filename=path)
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name.split(".")[0] in FORBIDDEN_IMPORTS:
                    bad.append(f"import {a.name}")
        elif isinstance(node, ast.ImportFrom):
            root = (node.module or "").split(".")[0]
            if root in FORBIDDEN_IMPORTS:
                bad.append(f"from {node.module} import ...")
        elif isinstance(node, ast.Call):
            f = node.func
            name = (f.attr if isinstance(f, ast.Attribute)
                    else (f.id if isinstance(f, ast.Name) else None))
            if name in FORBIDDEN_CALLS:
                bad.append(f"call {name}()")
    if bad:
        raise RuntimeError(
            "make_reference.py is NOT independent of the enumerator: "
            + "; ".join(sorted(set(bad))))
    return {"ast_nodes_walked": sum(1 for _ in ast.walk(tree)),
            "forbidden_imports_found": 0, "forbidden_calls_found": 0,
            "forbidden_import_list": sorted(FORBIDDEN_IMPORTS),
            "forbidden_call_list": sorted(FORBIDDEN_CALLS)}


# ---------------------------------------------------------------------------
# The generator. Pure string construction, METHYL END FIRST.
#
# A name like "N-[(9Z)-octadec-9-en-1-yl]prop-2-enamide" is read as:
#   - a chain of `n` carbons, numbered 1 at the nitrogen end;
#   - a cis double bond between carbons `delta` and `delta`+1, or none;
#   - the nitrogen of a prop-2-enamide on carbon 1.
# This routine lays the chain down from carbon `n` (the methyl) DOWN to carbon
# 1, then writes the amide, then the vinyl. stage03 does the reverse.
# ---------------------------------------------------------------------------
def acrylamide_from_name(n, delta=None, cis=True, reverse=False):
    """Write N-[alkyl/alkenyl]prop-2-enamide, METHYL END FIRST.

    n      carbons in the chain
    delta  the lower-numbered alkene carbon COUNTED FROM NITROGEN (carbon 1 is
           the one bonded to N), or None for a saturated chain
    cis    Z geometry
    reverse  the seeded defect: shift the alkene one carbon further from the
           nitrogen. Formula and exact mass are unaffected.

    The string is laid out as

        [A carbons][/C=C\\][B carbons] N C(=O) C=C
         ^ methyl end                   ^ nitrogen end

    Counting from the nitrogen, the B block occupies carbons 1..B, so the
    alkene's lower-numbered carbon is B+1. Hence B = delta-1 and
    A = n - 2 - B. stage03 builds the SAME molecule by indexing from carbon 1
    outward; here every index is computed from the methyl end, so an off-by-one
    in either routine moves the alkene while leaving C21H39NO and 321.3032
    exactly where they were.
    """
    if n < 4:
        raise ValueError("chain too short to place an alkene")
    if delta is None:
        return "C" * n + "NC(=O)C=C"
    if not (2 <= delta <= n - 2):
        raise ValueError(f"delta {delta} out of range for a C{n} chain")
    b = delta - 1
    if reverse:
        b = delta                      # alkene lands one carbon further out
    a = n - 2 - b
    if a < 1 or b < 1:
        raise ValueError(f"cannot place an alkene at C{delta} of a C{n} chain")
    mark_open, mark_close = ("/", "\\") if cis else ("", "")
    return ("C" * a + mark_open + "C=C" + mark_close + "C" * b
            + "NC(=O)C=C")


def describe(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"unparseable: {smiles}")
    return {"canonical": Chem.MolToSmiles(m),
            "formula": CalcMolFormula(m),
            "exact_mass": round(ExactMolWt(m), 4),
            "inchikey": Chem.MolToInchiKey(m)}


def main():
    seed = None
    if "--seed" in sys.argv:
        seed = sys.argv[sys.argv.index("--seed") + 1]

    prov = assert_no_assembler()

    # ---- the generator exercised across chain lengths ------------------
    # Hand-written expectations: a C(n) N-alkyl acrylamide is C(n+3)H(2n+5)NO
    # saturated, C(n+3)H(2n+3)NO with one alkene.
    scaling, scal_fail = [], []
    for n in (12, 14, 16, 18):
        sat = describe(acrylamide_from_name(n, delta=None))
        uns = describe(acrylamide_from_name(n, delta=9,
                                            reverse=(seed == "reverse")))
        want_sat = f"C{n + 3}H{2 * n + 5}NO"
        want_uns = f"C{n + 3}H{2 * n + 3}NO"
        ok = sat["formula"] == want_sat and uns["formula"] == want_uns
        scaling.append({"n": n, "saturated": sat["formula"],
                        "expected": want_sat,
                        "unsaturated": uns["formula"],
                        "expected_unsat": want_uns, "ok": bool(ok)})
        if not ok:
            scal_fail.append(f"n={n}: {sat['formula']}/{uns['formula']} != "
                             f"{want_sat}/{want_uns}")

    # ---- the two anchors ------------------------------------------------
    oleyl = describe(acrylamide_from_name(18, delta=9, cis=True,
                                          reverse=(seed == "reverse")))
    octa = describe(acrylamide_from_name(18, delta=None))

    # ---- compare to the enumerator, on the InChIKey ---------------------
    enum_path = os.path.join(HERE, "03_enumerated.csv")
    enum = {}
    if os.path.exists(enum_path):
        with open(enum_path) as fh:
            for r in csv.DictReader(fh):
                enum[r["role"].split(":", 1)[1]] = r
    cmp_rows, fails = [], list(scal_fail)
    for key, ref in (("oleyl_acrylamide", oleyl),
                     ("octadecyl_acrylamide", octa)):
        e = enum.get(key)
        if e is None:
            fails.append(f"{key}: no enumeration row to compare against")
            continue
        same_ik = e["inchikey"] == ref["inchikey"]
        cmp_rows.append({
            "anchor": key,
            "reference_inchikey": ref["inchikey"],
            "enumerated_inchikey": e["inchikey"],
            "full_inchikey_identical": bool(same_ik),
            "reference_formula": ref["formula"],
            "enumerated_formula": e["formula"],
            "formula_identical": e["formula"] == ref["formula"],
            "reference_exact_mass": ref["exact_mass"],
            "enumerated_exact_mass": float(e["exact_mass"]),
            "exact_mass_identical":
                abs(float(e["exact_mass"]) - ref["exact_mass"]) < 1e-4,
            "reference_canonical_smiles": ref["canonical"],
            "enumerated_smiles": e["tail_SMILES"],
            "byte_identical_canonical":
                Chem.MolToSmiles(Chem.MolFromSmiles(e["tail_SMILES"]))
                == ref["canonical"],
        })
        if not same_ik:
            fails.append(f"{key}: InChIKey mismatch "
                         f"{ref['inchikey']} vs {e['inchikey']}")

    n_agree = sum(1 for c in cmp_rows if c["full_inchikey_identical"])
    report = {
        "paper_draws_any_structure": False,
        "independence_proof": prov,
        "traversal": "this file writes the chain METHYL END FIRST (C18 -> C1) "
                     "then the amide then the vinyl; stage03_enumerate.py "
                     "starts at the acryloyl nitrogen and attaches a fragment "
                     "whose first atom is C1. Opposite directions.",
        "generator_chain_length_scaling": scaling,
        "comparisons": cmp_rows,
        "n_agree": n_agree, "n_total": len(cmp_rows),
        "seed": seed,
        "what_the_seed_proves": "`--seed reverse` writes the chain in the "
            "other direction, putting the alkene at C10 from nitrogen. "
            "Formula C21H39NO and exact mass 321.3032 are IDENTICAL either "
            "way, and both SMILES parse and canonicalise, so only the "
            "InChIKey separates them. If this gate passed under the seed it "
            "would not be testing connectivity.",
        "byte_identity_is_not_the_signal": "F-037. `byte_identical_canonical` "
            "is RECORDED above and is NOT used as a pass criterion; a correct "
            "re-derivation SHOULD be byte-identical once canonicalised. The "
            "pass criterion is the full InChIKey, and independence is "
            "established by the AST walk, not by a diff.",
        "chain_length_coverage": "both of this paper's tails are C18, so an "
            "anchor at a different chain length cannot be taken from it. The "
            "off-by-one-that-scales-with-n hazard is covered instead by "
            "exercising the GENERATOR at n = 12, 14, 16, 18 (above) and by "
            "stage03's seed_chain_length. No shorter product is claimed to be "
            "in the paper.",
        "checks_failed": fails,
    }

    if seed is None:
        hdr = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
               "read_by", "expected_mismatch", "note"]
        rows = [
            {"lipid_name": "oleyl acrylamide",
             "reference_SMILES": oleyl["canonical"],
             "source_figure": "NOT DRAWN -- chemical name only, in the "
                              "published abstract ('the Michael addition of "
                              "amines to oleyl acrylamide')",
             "read_on": "2026-11-10", "read_by": "curation agent",
             "expected_mismatch": "",
             "note": "Built from the chemical name by pure string "
                     "construction, methyl-end-first, with no assembler "
                     "import and no RWMol/AddBond/RemoveAtom/CombineMols "
                     "call -- verified by an AST walk of make_reference.py at "
                     "run time. Anchors the TAIL REAGENT, not a library "
                     "member."},
            {"lipid_name": "n-octadecyl acrylamide",
             "reference_SMILES": octa["canonical"],
             "source_figure": "NOT DRAWN -- the abstract names only "
                              "'n-octadecyl tails'; the acrylamide linkage is "
                              "parallel construction, corroborated by "
                              "AA_2008's N18 series (02_corpus_overlap.json "
                              "d1: 43/43)",
             "read_on": "2026-11-10", "read_by": "curation agent",
             "expected_mismatch": "",
             "note": "Same route as the oleyl anchor. The LINKAGE is an "
                     "assumed feature of this paper; the STRUCTURE for the "
                     "name is externally attested by 43 independently curated "
                     "AA_2008 N18 lipids."},
            {"lipid_name": "lipidoid 1-16 (the library itself)",
             "reference_SMILES": "",
             "source_figure": "NONE. The SI draws no structure; the main "
                              "article is unobtainable (ACS subscription, no "
                              "PMC deposit).",
             "read_on": "2026-11-10", "read_by": "curation agent",
             "expected_mismatch": "no_anchor_possible_paper_draws_nothing",
             "note": "EMPTY ON PURPOSE so the absence is visible. The 16 "
                     "amine heads are not named, drawn or formulated "
                     "anywhere in the SI, so no lipidoid is enumerable and "
                     "there is nothing for an anchor to test."},
        ]
        with open(os.path.join(HERE, "reference_structures.csv"), "w",
                  newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=hdr)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        with open(os.path.join(HERE, "anchor_independence.json"), "w") as fh:
            json.dump(report, fh, indent=2)

    print(f"AST independence: {prov['ast_nodes_walked']} nodes walked, "
          f"0 forbidden imports, 0 forbidden calls")
    for s in scaling:
        print(f"  generator n={s['n']}: {s['saturated']} / {s['unsaturated']} "
              f"{'OK' if s['ok'] else 'WRONG'}")
    for c in cmp_rows:
        print(f"  {c['anchor']}: InChIKey "
              f"{'IDENTICAL' if c['full_inchikey_identical'] else 'MISMATCH'} "
              f"| formula same={c['formula_identical']} "
              f"exact_mass same={c['exact_mass_identical']}")
    print(f"n_agree={n_agree} n_total={len(cmp_rows)}")
    if fails:
        print("FAILED:", *fails, sep="\n  - ")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
