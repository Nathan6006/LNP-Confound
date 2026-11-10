#!/usr/bin/env python3
"""C10 anchors for maohq_38583366, built from CHEMICAL NAMES, head first.

WHAT IS ANCHORED AND WHAT IS NOT. This paper enumerates nothing -- it has one
ionizable lipid and six catalogue helpers -- but Fig 1 DRAWS three of its
components in full: DLin-MC3-DMA, cholesterol and DMG-PEG2000. Two of those are
anchored here. The third is not, deliberately:

  * CHOLESTEROL -- NO ANCHOR, empty reference_SMILES. Hand-writing an
    eight-stereocentre steroid in this file would be a worse anchor than none:
    one agent's hand-written cholesterol came out a diastereomer with 4 of 8
    centres inverted, at identical formula, identical exact mass to six
    decimals and an IDENTICAL SKELETON InChIKey block -- F-013 reproduced from
    scratch. The registry string is NOT counted as an anchor either, because
    comparing the corpus to the corpus is comparing a value to itself. What
    Fig 1 does settle is recorded as prose in PROVENANCE.md: the C3 oxygen is
    on a BOLD WEDGE (3-beta, the exact feature F-013 inverted on 537 rows),
    the Delta-5 olefin is drawn, and the side chain is iso-octyl.

  * The six HELPER lipids -- NO ANCHOR, empty reference_SMILES. They are named
    in a box in Fig 1 and nowhere drawn, so this paper offers no external
    ground truth for them. That absence is the point of storing the empty row.
    (It is also why this paper cannot settle F-057; see 02_corpus_overlap.json.)

INDEPENDENCE IS ESTABLISHED AT RUN TIME, TWO WAYS, AND BOTH ARE MADE TO FAIL.
  1. `ast_independence()` walks THIS module's own parsed AST and inspects only
     Import / ImportFrom / Call / Attribute nodes, so the ban list cannot match
     itself -- a substring scan flagged its own docstring once, and another
     agent's guard failed on its own module for exactly that reason.
  2. `template_provenance()` requires that no string literal of >= 25 chars in
     this module be a substring of any target structure. A correct
     re-derivation is byte-identical to a correct original once canonicalised
     (F-037), so byte-identity is NOT the signal and testing for it is a trap.
  --selftest seeds a defect into each guard in a CHILD process and requires the
  alarm, with a clean control that must come back clean.

TRAVERSAL. The corpus writes these lipids TAIL FIRST. This file composes them
HEAD FIRST -- from the dimethylamino nitrogen outward for MC3, and from the
methyl cap of the polyether inward for DMG-PEG2000 -- and builds every chain by
multiplying a one-atom token, so a chain length is a number in one place and
never a hand-typed run of C's.
"""
import ast, json, os, subprocess, sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, Descriptors
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))

BANNED_IMPORTS = ("reagent_registry", "reagents", "corpus_components",
                  "assemble", "structure_check", "emit")
BANNED_CALLS = ("RWMol", "EditableMol", "AddBond", "RemoveAtom", "CombineMols",
                "ReplaceSubstructs", "SetChiralTag", "esterify", "attach")

CH2 = "C"
ME = "C"
# selftest() injects a seeded defect by replacing the single line below. The
# token is unique in this file so the substitution cannot damage the selftest's
# own source, which is what happened on the first attempt.
# __SEED_POINT__


def alkyl(n):
    """An unbranched saturated chain of n carbons, written left to right."""
    return CH2 * n


def dlin_mc3_dma():
    """heptatriaconta-6,9,28,31-tetraen-19-yl 4-(dimethylamino)butanoate.

    Head first: N(Me)2 -> three methylenes -> the ester carbonyl -> the ester
    oxygen -> C19 of the heptatriacontatetraene, which carries the two
    identical C18 arms. Each arm, counted OUTWARD from C19, is
    8 x CH2, then cis-CH=CH, then CH2, then cis-CH=CH, then 4 x CH2, then CH3
    = 18 carbons with two cis olefins five and eight bonds from the far end.
    """
    arm = (alkyl(8) + "/C=C\\" + CH2 + "/C=C\\" + alkyl(4) + ME)
    head = "C" + "N(" + ME + ")" + alkyl(3) + "C(=O)O"
    return head + "C(" + arm + ")" + arm


def dmg_peg2000(n_eo=45):
    """1,2-dimyristoyl-rac-glycero-3-methoxypolyethylene glycol.

    Head first from the methyl cap: CH3-O-(CH2CH2-O)n -> the glycerol sn-3
    methylene -> C2 bearing a myristoyl ester -> the sn-1 methylene bearing the
    other myristoyl. Fig 1 brackets ONE ethyleneoxy unit and subscripts it 44,
    with a further CH2CH2-O drawn OUTSIDE the bracket before the methyl cap, so
    the DRAWN total is 45 -- the subscript alone would be a 44 Da error that
    only the InChIKey catches.
    """
    myristoyl = "C(=O)" + alkyl(13)
    eo = "O" + CH2 + CH2
    return (ME + (eo * n_eo) + "OC" + "C(O" + myristoyl + ")"
            + "CO" + myristoyl)


TARGETS = [
    {"lipid_name": "DLin-MC3-DMA", "builder": dlin_mc3_dma,
     "source_figure": "Fig 1, 'Lipid Nanoparticle Library Design' box",
     "note": "Built from the IUPAC name, head first. Fig 1 draws the head as "
             "(CH3)2N-CH2CH2CH2-C(=O)-O- and the ester oxygen on the chain "
             "carbon, both confirmed by eye at 3x native; the tail vertex "
             "count was NOT counted by eye -- it comes from the name."},
    {"lipid_name": "DMG-PEG2000", "builder": dmg_peg2000,
     "source_figure": "Fig 1, 'Lipid Nanoparticle Library Design' box",
     "note": "n_eo = 45: 44 inside the drawn bracket plus one drawn outside "
             "it. The block length is a distribution, so this anchor tests "
             "CONNECTIVITY and the drawn unit count, not the reagent."},
]
EMPTY = [("Cholesterol", "Fig 1", "NO ANCHOR ON PURPOSE -- a hand-written "
          "8-stereocentre sterol is a worse anchor than none (F-013 was "
          "reproduced from scratch that way), and the registry string is the "
          "corpus, so comparing to it is comparing a value to itself. Fig 1 "
          "draws the C3 oxygen on a BOLD WEDGE (3-beta), the Delta-5 olefin "
          "and an iso-octyl side chain; that reading is in PROVENANCE.md."),
         ("DOTAP", "not drawn", "Named in a box in Fig 1, never drawn."),
         ("DDAB", "not drawn", "Named in a box in Fig 1, never drawn."),
         ("DOPE", "not drawn", "Named in a box in Fig 1, never drawn."),
         ("DSPC", "not drawn", "Named in a box in Fig 1, never drawn."),
         ("14PA", "not drawn", "Named in a box in Fig 1, never drawn."),
         ("18PG", "not drawn", "Named in a box in Fig 1, never drawn. This is "
          "why this paper cannot settle F-057.")]


def ast_independence(path=None):
    """Walk this module's own AST. Import/ImportFrom/Call/Attribute only."""
    path = path or os.path.abspath(__file__)
    tree = ast.parse(open(path).read())
    imports, calls = [], []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if any(b in a.name for b in BANNED_IMPORTS):
                    imports.append(a.name)
        elif isinstance(node, ast.ImportFrom):
            if node.module and any(b in node.module for b in BANNED_IMPORTS):
                imports.append(node.module)
        elif isinstance(node, ast.Call):
            f = node.func
            name = getattr(f, "id", None) or getattr(f, "attr", None)
            if name in BANNED_CALLS:
                calls.append(name)
        elif isinstance(node, ast.Attribute):
            if node.attr in BANNED_CALLS:
                calls.append(node.attr)
    return {"forbidden_imports_found": imports, "forbidden_calls_found": calls,
            "independent": not imports and not calls}


def template_provenance(structures, path=None):
    """No long literal in this module may be a substring of a target."""
    path = path or os.path.abspath(__file__)
    tree = ast.parse(open(path).read())
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) \
                and len(node.value) >= 25:
            for s in structures:
                if node.value in s:
                    bad.append(node.value[:60])
    return {"long_literals_found_inside_a_target": bad, "clean": not bad}


def build():
    rows, qc = [], []
    for t in TARGETS:
        smi = t["builder"]()
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise SystemExit("make_reference: %s does not parse" % t["lipid_name"])
        rows.append({"lipid_name": t["lipid_name"],
                     "reference_SMILES": Chem.MolToSmiles(m),
                     "source_figure": t["source_figure"],
                     "read_on": "2026-11-10", "read_by": "curation agent",
                     "expected_mismatch": "", "note": t["note"]})
        qc.append({"lipid_name": t["lipid_name"],
                   "formula": rdMolDescriptors.CalcMolFormula(m),
                   "exact_mass": round(Descriptors.ExactMolWt(m), 6),
                   "inchikey": Chem.MolToInchiKey(m)})
    for name, fig, note in EMPTY:
        rows.append({"lipid_name": name, "reference_SMILES": "",
                     "source_figure": fig, "read_on": "2026-11-10",
                     "read_by": "curation agent", "expected_mismatch": "",
                     "note": note})
    return rows, qc


def selftest():
    """Seed a defect into each guard, in a child process, and require it."""
    import tempfile, shutil, re
    src = open(os.path.abspath(__file__)).read()
    tmp = tempfile.mkdtemp(prefix="maohq38_anchor_")
    results = []
    ctrl = os.path.join(tmp, "ctrl.py")
    open(ctrl, "w").write(src)
    r = subprocess.run([sys.executable, "-I", ctrl, "--guards"],
                       capture_output=True, text=True)
    results.append({"seed": "CONTROL", "returncode": r.returncode,
                    "verdict": "PASS" if r.returncode == 0
                    else "FAIL -- control not clean"})
    if r.returncode != 0:
        print(r.stdout, r.stderr)
        raise SystemExit("make_reference selftest: control failed")

    marker = "# __SEED" + "_POINT__"
    s1 = src.replace(marker, "def _never_called():\n"
                     "    from reagent_registry import REAGENTS\n")
    p1 = os.path.join(tmp, "seed_import.py"); open(p1, "w").write(s1)
    r1 = subprocess.run([sys.executable, "-I", p1, "--guards"],
                        capture_output=True, text=True)
    results.append({"seed": "forbidden_import(reagent_registry, in a never-called function so the AST sees it and nothing executes)",
                    "returncode": r1.returncode,
                    "guard_message": [l for l in (r1.stdout + r1.stderr)
                                      .splitlines() if "GUARD" in l][:2],
                    "verdict": "CREDITED" if r1.returncode != 0 and
                    "forbidden_imports_found" in (r1.stdout + r1.stderr)
                    else "NOT CREDITED"})

    target = Chem.MolToSmiles(Chem.MolFromSmiles(dlin_mc3_dma()))
    lit = target[:40]
    s2 = src.replace(marker, "_COPIED = %r" % lit)
    p2 = os.path.join(tmp, "seed_literal.py"); open(p2, "w").write(s2)
    r2 = subprocess.run([sys.executable, "-I", p2, "--guards"],
                        capture_output=True, text=True)
    results.append({"seed": "template_copy(40-char slice of the target)",
                    "returncode": r2.returncode,
                    "guard_message": [l for l in (r2.stdout + r2.stderr)
                                      .splitlines() if "GUARD" in l][:2],
                    "verdict": "CREDITED" if r2.returncode != 0 and
                    "long_literals_found_inside_a_target"
                    in (r2.stdout + r2.stderr) else "NOT CREDITED"})
    shutil.rmtree(tmp, ignore_errors=True)
    return results


def guards_or_die():
    rows, qc = build()
    structures = [r["reference_SMILES"] for r in rows if r["reference_SMILES"]]
    a = ast_independence()
    t = template_provenance(structures)
    print("GUARD ast_independence %s" % json.dumps(a))
    print("GUARD template_provenance %s" % json.dumps(t))
    if not a["independent"] or not t["clean"]:
        raise SystemExit("make_reference: INDEPENDENCE GUARD FAILED")
    return rows, qc, a, t


def main():
    if "--guards" in sys.argv:
        guards_or_die()
        return
    rows, qc, a, t = guards_or_die()
    import csv
    cols = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
            "read_by", "expected_mismatch", "note"]
    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    seeds = selftest() if "--selftest" in sys.argv else None
    out = {"guards": {"ast_independence": a, "template_provenance": t},
           "built": qc, "empty_on_purpose": [e[0] for e in EMPTY],
           "selftest": seeds}
    with open(os.path.join(HERE, "04_anchor_build.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
