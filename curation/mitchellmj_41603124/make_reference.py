#!/usr/bin/env python3
"""The C10 anchor for `mitchellmj_41603124`, BUILT, not copied.

C10 only tests anything if the anchor came from a different code path than the
enumerator, and that has been claimed twice in this project and been false
twice (F-013, a byte copy; F-037, the file written to fix F-013 re-copied the
other half). So:

  * 03_enumerate.py composes 18D4 by PLAIN SMILES STRING INTERPOLATION,
    arm-first: it writes one arm literal, then splices it into a template.
  * this file builds the same molecule with curation/lib/assemble, HEAD-FIRST:
    it starts at the piperazine, walks outward to each nitrogen in turn, and
    makes every bond with `assemble.attach` on atom-mapped attachment atoms.
    There is no f-string of a fragment anywhere in the build.

`_assert_independent()` enforces that by WALKING THIS MODULE'S OWN AST at run
time -- not by scanning the source text, which would flag this very docstring.
It requires that `build()` contains no JoinedStr/BinOp-on-str/%-format and at
least one call to assemble.attach, and it is made to FAIL ON PURPOSE by
`--selftest`, which re-runs it against a deliberately string-built function.

THE GATE TESTS CONNECTIVITY, and `--selftest` proves it by seeding a REGIOISOMER
of 18D4 -- both arms of the secondary nitrogen moved onto one primary nitrogen
instead. Identical formula C46H85N5O12, identical average mass to 3 dp and
identical monoisotopic mass to 6 dp. Only the InChIKey sees it.

17D9.2 and 20D8i get EMPTY reference_SMILES rows, so the absence of an anchor
is visible rather than implied. The reason is on the row.
"""
import argparse, ast, csv, inspect, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "provenance"))

import assemble                                           # noqa: E402
from rdkit import Chem, RDLogger                          # noqa: E402
from rdkit.Chem import Descriptors, Draw                  # noqa: E402
RDLogger.DisableLog("rdApp.*")

OUT = os.path.join(HERE, "reference_structures.csv")

# The two reagents as Supplementary Figure 2 draws them, with the attachment
# atoms atom-mapped. Nothing here is taken from 03_enumerate.py.
PIPERAZINE_CORE = "C1CN(CC[NH:1])CCN1CC[NH:2]"     # amine 18, both primary N
SEC_N = "[NH:3]"                                    # its secondary N-H
ARM_A = "[CH:21](CC(=O)OCCCC)C(=O)OCCCC"   # the maleate after aza-Michael
ARM_B = "[CH:22](CC(=O)OCCCC)C(=O)OCCCC"
ARM_C = "[CH:23](CC(=O)OCCCC)C(=O)OCCCC"


def build():
    """18D4, head-first from the piperazine. One arm per N-H.

    Every fragment below is a LITERAL. There is no f-string, no `+` on a str
    and no %-format anywhere in this function; `_assert_independent` walks this
    function's AST and refuses if one appears.
    """
    amine = "C1CN(CC[NH2:1])CCN1CC[NH:2]CC[NH2:3]"
    prod = assemble.attach(amine, ARM_A, a_map=1, b_map=21,
                           bond="single", label="18D4 arm on primary N (ring side)")
    prod = assemble.attach(prod, ARM_B, a_map=2, b_map=22,
                           bond="single", label="18D4 arm on the secondary N")
    prod = assemble.attach(prod, ARM_C, a_map=3, b_map=23,
                           bond="single", label="18D4 arm on primary N (tail)")
    return assemble.canonical(prod)


def build_seeded():
    """A REGIOISOMER of 18D4, as a literal: the secondary N keeps its H and one
    primary N carries BOTH of its hydrogens' worth of arms.

    Same formula, same average mass, same monoisotopic mass to 6 dp. A
    different molecule. This is what proves the C10 gate reads connectivity
    and not composition.
    """
    return assemble.canonical(
        "CCCCOC(=O)CC(NCCN1CCN(CCNCCN(C(CC(=O)OCCCC)C(=O)OCCCC)"
        "C(CC(=O)OCCCC)C(=O)OCCCC)CC1)C(=O)OCCCC")


def _string_built():                      # the NEGATIVE control for the guard
    frag = "CCCC"
    return "O" + frag + "C(=O)C"


def _assert_independent(fn, must_have_attach=True):
    """Walk FN's own AST. Raise if it builds SMILES by string surgery."""
    src = inspect.getsource(fn)
    tree = ast.parse(src.lstrip())
    bad, calls = [], []
    for node in ast.walk(tree):
        if isinstance(node, ast.JoinedStr):
            bad.append("f-string")
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Mod)):
            bad.append("str concat/%% at line %d" % node.lineno)
        if isinstance(node, ast.Call):
            f = node.func
            name = (f.attr if isinstance(f, ast.Attribute)
                    else getattr(f, "id", ""))
            calls.append(name)
            if name in ("join", "format"):
                bad.append("str.%s" % name)
    if bad:
        raise AssertionError("NOT INDEPENDENT: %s builds SMILES by string "
                             "surgery: %s" % (fn.__name__, sorted(set(bad))))
    if must_have_attach and "attach" not in calls:
        raise AssertionError("NOT INDEPENDENT: %s never calls assemble.attach"
                             % fn.__name__)
    return sorted(set(calls))


def masses(smiles):
    m = Chem.MolFromSmiles(smiles)
    from rdkit.Chem.Descriptors import ExactMolWt
    return (Chem.rdMolDescriptors.CalcMolFormula(m), round(Descriptors.MolWt(m), 3),
            round(ExactMolWt(m), 6), Chem.MolToInchiKey(m))


def selftest():
    ok = True
    # 1. the guard must PASS on build()
    try:
        calls = _assert_independent(build)
        print("PASS  guard accepts build(); calls seen: %s" % calls)
    except AssertionError as e:
        print("FAIL  guard rejected build(): %s" % e); ok = False
    # 2. THE GUARD MUST FAIL ON PURPOSE
    try:
        _assert_independent(_string_built)
        print("FAIL  guard accepted a string-built function -- it sees nothing")
        ok = False
    except AssertionError as e:
        print("PASS  guard fires on a string-built function: %s" % str(e)[:90])
    # 3. the ANCHOR must separate a same-mass regioisomer
    a, b = build(), build_seeded()
    fa, wa, ea, ka = masses(a)
    fb, wb, eb, kb = masses(b)
    print("      anchor   %s  avg=%.3f  exact=%.6f  %s" % (fa, wa, ea, ka))
    print("      seeded   %s  avg=%.3f  exact=%.6f  %s" % (fb, wb, eb, kb))
    same_mass = (fa == fb and wa == wb and ea == eb)
    if same_mass and ka != kb:
        print("PASS  seeded regioisomer: formula, average mass AND exact mass "
              "identical (delta 0.000000 Da); InChIKey differs -> the gate "
              "tests CONNECTIVITY, not composition")
    elif not same_mass:
        print("FAIL  the seed changed the mass, so it does not prove the gate "
              "tests connectivity (formula %s vs %s)" % (fa, fb)); ok = False
    else:
        print("FAIL  the seeded regioisomer has the SAME InChIKey"); ok = False
    return 0 if ok else 1


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--selftest", action="store_true")
    a = p.parse_args()
    if a.selftest:
        return selftest()

    _assert_independent(build)
    smi = build()
    formula, avg, exact, key = masses(smi)
    print("18D4 anchor: %s  %s  avg=%.3f  %s" % (smi, formula, avg, key))

    rows = [
        ["18D4", smi, "Supplementary Figure 2 (SI p.16), product panel",
         "2026-11-10", "make_reference.py::build (assemble.attach, head-first "
                       "from the piperazine; no string surgery, enforced by "
                       "_assert_independent walking this module's AST)",
         "",
         "The one member this paper both draws and confirms. Printed "
         "[M+H]+ 901.829 (Supp Fig 5) against 901.217 for this structure, "
         "delta 0.612 Da; the nearest alternative (one CH2 on each of the two "
         "ester chains) is 28.053 Da away, a 45.8x discrimination margin. "
         "Three arms, one on each N-H; the two piperazine tertiary N are "
         "untouched, which is the scheme's own 'favor mono-addition'."],
        ["17D9.2", "", "Supplementary Figure 2 (SI p.16), product panel",
         "2026-11-10", "NOT BUILT -- deliberately empty", "",
         "DRAWN but NOT DETERMINABLE. The three printed ions pin the two ester "
         "chains at a+b = 18.993 carbons total (residual -0.007 from 19), on "
         "TREN with three arms. But mass cannot see connectivity: 11+8 and "
         "10+9 are isomeric, as is every branch position, and the SI raster is "
         "137 dpi. No structure is recorded rather than a guessed one."],
        ["20D8i", "", "Supplementary Figure 2 (SI p.16), product panel",
         "2026-11-10", "NOT BUILT -- deliberately empty", "",
         "DRAWN but NOT DETERMINABLE. M = 1878.58 +/- 0.41 Da from the printed "
         "[M+3H]3+ and [M+2H]2+. Amine 20 is a branched polyamidoamine whose "
         "core and arm count were not resolved at 137 dpi, so the arithmetic "
         "has two unknowns and is not solvable."],
    ]
    # The excipients and benchmarks are NAMED (Supp Table 3), never DRAWN, so
    # there is nothing in this paper to anchor against. Recorded as such.
    for n in ("DOPE", "Cholesterol", "DSPC", "SM-102", "DLin-MC3-DMA",
              "ALC-0315"):
        rows.append([n, "", "Supplementary Table 3 (SI p.11-13) -- NAMED only",
                     "2026-11-10", "not built: the paper draws no excipient "
                     "and no benchmark lipid anywhere", "",
                     "No anchor is possible from this paper: Supp Table 3 "
                     "gives the NAME and the mol%, and Supplementary Figure 2 "
                     "draws only the three nDm products. Recorded instead in "
                     "07_qc.json: all six names resolve identically in "
                     "provenance/reagent_registry.py and in the live corpus, "
                     "6/6 by FULL InChIKey."])

    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        w.writerows(rows)
    print("wrote %s (%d rows, %d with a SMILES)"
          % (OUT, len(rows), sum(1 for r in rows if r[1])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
