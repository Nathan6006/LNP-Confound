#!/usr/bin/env python3
"""Build the C10 anchor(s) for this folder and write reference_structures.csv.

WHAT THIS PAPER GIVES AN ANCHOR FOR, AND WHAT IT DOES NOT.

main-article Figure 1A DRAWS every lipid component, cKK-E12 and MC3 included.
But the only copy of Figure 1 obtainable is 708 x 391 px for twelve panels
(source/pmc_oa/nn5c12080_0001.jpg): the cKK-E12 panel is 130 x 70 px. At that
size the ring carbonyls, the two NH, the two tertiary nitrogens and the
'C10H21' arm labels are all legible, but the number of methylenes in each
lysine side chain -- the one thing an anchor exists to test -- is not. Counting
them by eye is exactly what this project forbids, and a vertex trace needs
resolution this raster does not have. So cKK-E12 gets an EMPTY
reference_SMILES, with the absence recorded rather than papered over.

MC3 does get a real anchor, and it is honest about what it tests. It is NOT
derived from the drawing; it is derived from the IUPAC name that
provenance/reagent_registry.py stores as that reagent's authority,
   (6Z,9Z,28Z,31Z)-heptatriaconta-6,9,28,31-tetraen-19-yl
   4-(dimethylamino)butanoate
by walking the locants and composing the chain carbon by carbon. The registry
holds a literal SMILES; this file holds none. Agreement means two independent
routes from the same chemical name reach the same molecule.

INDEPENDENCE IS VERIFIED BY PROVENANCE, NOT BY BYTES (F-037). Byte-identity of
the emitted anchor is not the signal -- a correct re-derivation is
byte-identical to a correct enumeration once canonicalised. So `_ast_guard`
walks THIS MODULE'S OWN AST at run time and refuses if it imports an assembler
or calls RWMol / AddBond / RemoveAtom / CombineMols / esterify / attach. The
guard is itself made to fail on purpose by --prove-guard, which writes a
temporary copy with `import assemble` added and requires the guard to reject
it.

THE GATE IS SHOWN TO TEST CONNECTIVITY by --seed-regioisomer, which moves the
ester from C19 to C18. That is a regioisomer: identical molecular formula and
identical exact mass to every decimal place. Only the InChIKey sees it.

Run:  python3 -I make_reference.py                 writes reference_structures.csv
      python3 -I make_reference.py --prove-guard    the AST guard must reject a
                                                    seeded copy (exit 0 on success)
      python3 -I make_reference.py --seed-regioisomer
                                                    must exit 1 on an InChIKey
                                                    mismatch at identical mass
"""
import ast
import csv
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

from rdkit import Chem, RDLogger                      # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors  # noqa: E402

RDLogger.DisableLog("rdApp.*")

BANNED_CALLS = {"RWMol", "AddBond", "RemoveAtom", "CombineMols", "esterify",
                "attach", "ReplaceSubstructs", "EditableMol"}
BANNED_IMPORTS = {"assemble", "curation.lib.assemble", "build_new_curation"}

HEADER = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
          "read_by", "expected_mismatch", "note"]


def _ast_guard(path):
    """Refuse if this module could have borrowed the enumerator's machinery.

    Walks the AST, not the text. A substring scan over the file would flag this
    module's own docstring, which is how a sibling's first attempt failed.
    """
    tree = ast.parse(open(path).read())
    bad = []
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            for a in n.names:
                if a.name.split(".")[0] in BANNED_IMPORTS or a.name in BANNED_IMPORTS:
                    bad.append(f"import {a.name}")
        elif isinstance(n, ast.ImportFrom):
            mod = (n.module or "")
            if mod.split(".")[0] in BANNED_IMPORTS or mod in BANNED_IMPORTS:
                bad.append(f"from {mod} import ...")
        elif isinstance(n, ast.Call):
            f = n.func
            name = getattr(f, "attr", None) or getattr(f, "id", None)
            if name in BANNED_CALLS:
                bad.append(f"call {name}()")
    return bad


# --- the derivation: locants -> SMILES, composed carbon by carbon ----------
NAME = ("(6Z,9Z,28Z,31Z)-heptatriaconta-6,9,28,31-tetraen-19-yl "
        "4-(dimethylamino)butanoate")
CHAIN = 37                       # heptatriaconta-
ENE = [6, 9, 28, 31]             # -6,9,28,31-tetraen-, all Z
ACYL = "OC(=O)CCCN(C)C"          # 4-(dimethylamino)butanoate, as its ester O


def from_locants(ester_at):
    """Compose the SMILES by walking C1..C37 and reading the locants.

    No molecule object is edited anywhere: the result is a string, built one
    carbon at a time, and only then parsed.
    """
    dbl = {(i, i + 1) for i in ENE}
    starts = {a for a, _ in dbl}
    ends = {b for _, b in dbl}
    out = []
    for i in range(1, CHAIN + 1):
        if i > 1:
            a, b = i - 1, i
            if (a, b) in dbl:
                out.append("=")
            elif b in starts:
                out.append("/")      # the single bond BEFORE a Z alkene
            elif a in ends:
                out.append("\\")     # the single bond AFTER it
        out.append("C")
        if i == ester_at:
            out.append("(" + ACYL + ")")
    return "".join(out)


def main(argv):
    bad = _ast_guard(os.path.abspath(__file__))
    if bad:
        print("make_reference: AST guard REFUSED this module -- it is not "
              "independent of the enumerator: " + "; ".join(bad),
              file=sys.stderr)
        return 1

    if "--prove-guard" in argv:
        src = open(os.path.abspath(__file__)).read()
        seeded = src.replace("import ast\n", "import ast\nimport assemble\n", 1)
        assert seeded != src
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "seeded.py")
            open(p, "w").write(seeded)
            found = _ast_guard(p)
        ok = any("assemble" in f for f in found)
        print(f"AST guard on a copy with 'import assemble' added: "
              f"{'REJECTED' if ok else 'ACCEPTED'} ({found})")
        return 0 if ok else 1

    seed = "--seed-regioisomer" in argv
    smi = from_locants(18 if seed else 19)
    m = Chem.MolFromSmiles(smi)
    if m is None:
        print("make_reference: derived SMILES does not parse.", file=sys.stderr)
        return 1

    # the enumerated structure this anchor is checked against
    with open(os.path.join(HERE, "03_enumerated.csv")) as fh:
        enum = {r["formulation_id"]: r for r in csv.DictReader(fh)}
    tgt = Chem.MolFromSmiles(enum["MC3"]["IL_SMILES"])

    f_a = rdMolDescriptors.CalcMolFormula(m)
    f_b = rdMolDescriptors.CalcMolFormula(tgt)
    m_a = Descriptors.ExactMolWt(m)
    m_b = Descriptors.ExactMolWt(tgt)
    k_a, k_b = Chem.MolToInchiKey(m), Chem.MolToInchiKey(tgt)

    print(f"derived from the name, ester at C{18 if seed else 19}:")
    print(f"  SMILES   {smi}")
    print(f"  formula  {f_a}   (enumerated: {f_b})   same={f_a == f_b}")
    print(f"  mass     {m_a:.6f} (enumerated: {m_b:.6f})  "
          f"delta={abs(m_a - m_b):.6f}")
    print(f"  InChIKey {k_a}")
    print(f"           {k_b}   same={k_a == k_b}")
    print(f"  byte-identical to the registry string? "
          f"{smi == enum['MC3']['IL_SMILES']}  "
          f"(NOT the signal -- see F-037)")

    if seed:
        if k_a == k_b:
            print("SEED NOT CAUGHT: the C18 regioisomer matched.",
                  file=sys.stderr)
            return 1
        print("SEED CAUGHT: identical formula and identical exact mass, "
              "different InChIKey. The gate tests connectivity.")
        return 1

    if k_a != k_b:
        print("make_reference: the anchor does NOT match the enumerated "
              "structure. Refusing to write a passing row.", file=sys.stderr)
        return 1

    rows = [
        {"lipid_name": "DLin-MC3-DMA",
         "reference_SMILES": smi,
         "source_figure": "NOT a figure: the IUPAC name recorded as this "
                          "reagent's authority in provenance/reagent_registry.py. "
                          "Figure 1A of the paper does draw MC3, but the only "
                          "copy obtainable is 708x391 px for twelve panels.",
         "read_on": "2026-11-10",
         "read_by": "make_reference.py from_locants(): chain length, alkene "
                    "locants and ester locant walked out of the name and "
                    "composed carbon by carbon as a string; no molecule object "
                    "is edited (AST-guarded, guard proved by --prove-guard)",
         "expected_mismatch": "",
         "note": "Exact InChIKey match to the enumerated structure. What this "
                 "tests is the CORPUS/REGISTRY string against the chemical "
                 "name, by a second route; it does not test the paper's own "
                 "drawing, which is unreadable here. The gate is shown to see "
                 "connectivity by --seed-regioisomer (ester C19->C18: "
                 "identical formula, identical exact mass, different key)."},
        {"lipid_name": "cKK-E12",
         "reference_SMILES": "",
         "source_figure": "main-article Figure 1A, cKK-E12 panel "
                          "(source/pmc_oa/nn5c12080_0001.jpg)",
         "read_on": "2026-11-10",
         "read_by": "",
         "expected_mismatch": "",
         "note": "EMPTY ON PURPOSE. The paper DOES draw this lipid, but the "
                 "panel is 130 x 70 px inside a 708 x 391 px figure, the only "
                 "copy obtainable (Europe PMC fulltextRepo returned an HTML "
                 "error page and the PMC OA service 404'd; no bot check was "
                 "bypassed). At that size the 2,5-diketopiperazine with its "
                 "two ring carbonyls and two NH, the two tertiary nitrogens, "
                 "the four CH2-CH(OH)-R arms and their 'C10H21' labels are all "
                 "legible and CONSISTENT with the emitted structure -- but the "
                 "methylene count of each lysine side chain, the one thing an "
                 "anchor exists to test, is not resolvable, and counting it by "
                 "eye is what this project forbids. No anchor is better than a "
                 "guessed one."},
    ]
    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER)
        w.writeheader()
        w.writerows(rows)
    print(f"\nreference_structures.csv: {len(rows)} rows, "
          f"{sum(1 for r in rows if not r['reference_SMILES'])} with an empty "
          f"reference_SMILES (the absence is visible)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
