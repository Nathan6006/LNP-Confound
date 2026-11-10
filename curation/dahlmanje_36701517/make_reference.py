"""Build the C10 reference structures for Radmand 2023 (PMID 36701517).

WHAT THIS PAPER DRAWS, AND WHAT IT DOES NOT
-------------------------------------------
The ionizable lipid is cKK-E12 in all 144 formulations and the paper NEVER
DRAWS IT. It is named once, with a vendor and a catalogue number (Organix Inc.,
O-8744), in the SI Materials. So cKK-E12 gets a row with an EMPTY
reference_SMILES, which is the protocol's way of making the ABSENCE of an
anchor visible rather than inventing one.

What the paper does draw in full are the helper lipids, and it draws the three
that matter TWICE each -- once in the Supplementary Figure 2 catalogue and
again beside the LNP that uses it:

    18:0 DDAB   Fig S2a row C3   and Fig 3a  (Cat-LNP)
    18:0 PC     Fig S2b row N7   and Fig S15a (Neu-LNP, labelled DSPC there)
    18:1 PA     Fig S2c row A1   and Fig S15c (An-LNP)

INDEPENDENCE (F-013, F-037)
---------------------------
These strings are composed here from fragments named after what the figure
draws. They are NOT taken from provenance/reagent_registry.py, from
curation/corpus_components.py, or from any corpus file -- `_ast_guard()` below
walks THIS MODULE'S OWN AST at run time and refuses if it imports an assembler
or calls RWMol / AddBond / RemoveAtom / CombineMols / esterify / attach, and
`_no_external_source()` refuses if any registry or corpus path is referenced.
The guard is made to fail on purpose by `--seed-guard`.

WHY STRING COMPOSITION IS SAFE *HERE*, SPECIFICALLY
---------------------------------------------------
CLAUDE.md mandates fragment assembly because string concatenation once fused a
piperazine to cholesterol through a RING-CLOSURE DIGIT COLLISION. All three
molecules built here are ACYCLIC -- zero rings, so there is no ring-closure
digit anywhere in any fragment and that collision cannot occur. The assertions
below check that explicitly (ring count 0) rather than leaving it as a claim.

Usage:
    python3 -I make_reference.py                 # write reference_structures.csv
    python3 -I make_reference.py --seed-ester    # seed an isomeric defect
    python3 -I make_reference.py --seed-guard    # make the AST guard fire
"""
import ast
import csv
import os
import sys

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

FORBIDDEN_CALLS = {"RWMol", "AddBond", "RemoveAtom", "CombineMols",
                   "esterify", "attach", "EditableMol", "ReplaceSubstructs"}
FORBIDDEN_IMPORTS = {"assemble", "reagents", "corpus_components",
                     "reagent_registry", "structure_check"}
FORBIDDEN_STRINGS = ("reagent_registry", "corpus_components",
                     "data_vitro", "data_vivo", "03_enumerated")


def _ast_guard(path=None, _extra_src=""):
    """Walk THIS file's AST. Raise if it could have copied the enumeration."""
    src = _extra_src or open(path or os.path.abspath(__file__)).read()
    tree = ast.parse(src)
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name.split(".")[0] in FORBIDDEN_IMPORTS:
                    bad.append(f"import {a.name}")
        elif isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] in FORBIDDEN_IMPORTS:
                bad.append(f"from {node.module} import ...")
        elif isinstance(node, ast.Call):
            f = node.func
            nm = getattr(f, "attr", None) or getattr(f, "id", None)
            if nm in FORBIDDEN_CALLS:
                bad.append(f"call {nm}()")
    # The literals that could NAME an external structure source. Prose is
    # exempt -- this module's own docstrings have to be able to say which
    # files it is independent of, and the guard's own blocklist has to be able
    # to spell them. Only executable literals are scanned, which is where a
    # path to a registry or a corpus file would actually appear.
    exempt = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.ClassDef)):
            body = getattr(node, "body", [])
            if body and isinstance(body[0], ast.Expr) \
                    and isinstance(body[0].value, ast.Constant):
                exempt.add(id(body[0].value))
        if isinstance(node, ast.Assign):
            names = [t.id for t in node.targets if isinstance(t, ast.Name)]
            if any(n.startswith("FORBIDDEN") for n in names):
                for sub in ast.walk(node.value):
                    exempt.add(id(sub))
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) \
                and id(node) not in exempt:
            for s in FORBIDDEN_STRINGS:
                if s in node.value:
                    bad.append(f"executable string literal naming {s!r}")
    if bad:
        raise SystemExit(
            "make_reference: AST GUARD FIRED -- this module is not an "
            "independent code path.\n  " + "\n  ".join(sorted(set(bad))) +
            "\nA C10 anchor built by, or copied from, the enumerator compares "
            "a value to itself (F-013, F-037).")
    return {"forbidden_calls_found": 0, "forbidden_imports_found": 0,
            "nodes_walked": sum(1 for _ in ast.walk(tree))}


# --- fragments, named after what the figure DRAWS --------------------------
#
# A SMILES FRAGMENT HAS A DIRECTION THAT MASS CANNOT SEE. An acyl group spliced
# in BEFORE the ester oxygen must end with its carbonyl carbon; the same group
# spliced in AFTER the ester oxygen must start with it. Writing the wrong one
# gives a molecule with the IDENTICAL FORMULA and the IDENTICAL EXACT MASS that
# is a ketone bearing an ether, not a diester. That is not a hypothetical: the
# first version of this file did exactly that and reported C44H88NO8P at
# 789.624755 Da, which is right for DSPC, while the structure was wrong. It is
# the seeded defect below, and `_direction_checks` is what catches it.
C18_SATURATED_CHAIN = "CCCCCCCCCCCCCCCCCC"            # 18 carbons, no C=C
STEAROYL_LEFT = "CCCCCCCCCCCCCCCCCC(=O)"              # 18:0 acyl, carbonyl LAST
STEAROYL_RIGHT = "C(=O)CCCCCCCCCCCCCCCCC"             # 18:0 acyl, carbonyl FIRST
OLEOYL_LEFT = "CCCCCCCC/C=C\\CCCCCCCC(=O)"            # 18:1 cis-9, carbonyl LAST
OLEOYL_RIGHT = "C(=O)CCCCCCC/C=C\\CCCCCCCC"           # 18:1 cis-9, carbonyl FIRST
PHOSPHOCHOLINE_HEAD = "P(=O)([O-])OCC[N+](C)(C)C"
PHOSPHATE_HEAD = "P(=O)(O)[O-]"
# sn-1 acyl is prepended; this is the glycerol with its sn-3 head and sn-2 acyl
SN_GLYCERO_BACKBONE = "OC[C@H](CO{head})OC(=O){sn2_tail}"


def build(seed_ester=False):
    """Compose the three drawn helper lipids. No molecule here has a ring."""
    out = {}

    # 18:0 DDAB -- Fig S2a row C3 and Fig 3a. Two C18 saturated chains on a
    # permanently quaternary nitrogen bearing two methyls. Counterion Br-,
    # drawn but not part of the lipid cation; stored as the cation, which is
    # how every phospholipid in this corpus is stored.
    out["18:0 DDAB"] = (C18_SATURATED_CHAIN + "[N+](C)(C)"
                        + C18_SATURATED_CHAIN)

    # 18:0 PC (DSPC) -- Fig S2b row N7 and Fig S15a. Two stearoyl chains on
    # sn-1 and sn-2 of a glycerophosphocholine; the sn-2 H is drawn on a wedge.
    sn1 = STEAROYL_RIGHT if seed_ester else STEAROYL_LEFT
    out["18:0 PC (DSPC)"] = sn1 + SN_GLYCERO_BACKBONE.format(
        head=PHOSPHOCHOLINE_HEAD, sn2_tail="CCCCCCCCCCCCCCCCC")

    # 18:1 PA -- Fig S2c row A1 and Fig S15c. Two oleoyl (18:1 cis-9) chains
    # and a free phosphate drawn as P(=O)(OH)(O-) with Na+.
    out["18:1 PA"] = OLEOYL_LEFT + SN_GLYCERO_BACKBONE.format(
        head=PHOSPHATE_HEAD, sn2_tail="CCCCCCC/C=C\\CCCCCCCC")
    return out


# name -> (n ester groups, n ketones, n dialkyl ethers, n quaternary N)
EXPECTED_CONNECTIVITY = {
    "18:0 DDAB": (0, 0, 0, 1),
    "18:0 PC (DSPC)": (2, 0, 0, 1),
    "18:1 PA": (2, 0, 0, 0),
}
_ESTER = Chem.MolFromSmarts("[CX4][CX3](=O)[OX2][CX4]")
_KETONE = Chem.MolFromSmarts("[CX4][CX3](=O)[CX4]")
_ETHER = Chem.MolFromSmarts("[CX4][OX2][CX4]")
_QUAT_N = Chem.MolFromSmarts("[NX4+]")


def _direction_checks(name, m):
    """The assertions a formula and an exact mass cannot make."""
    got = (len(m.GetSubstructMatches(_ESTER)),
           len(m.GetSubstructMatches(_KETONE)),
           len(m.GetSubstructMatches(_ETHER)),
           len(m.GetSubstructMatches(_QUAT_N)))
    want = EXPECTED_CONNECTIVITY[name]
    if got != want:
        raise SystemExit(
            f"{name}: connectivity (esters, ketones, dialkyl ethers, "
            f"quaternary N) is {got}, expected {want}. A reversed acyl "
            "fragment produces exactly this signature at an unchanged formula "
            "and an unchanged exact mass.")
    return {"n_esters": got[0], "n_ketones": got[1],
            "n_dialkyl_ethers": got[2], "n_quaternary_N": got[3]}


def check(name, smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise SystemExit(f"{name}: does not parse: {smiles}")
    if any(a.GetNumRadicalElectrons() for a in m.GetAtoms()):
        raise SystemExit(f"{name}: radical electron present.")
    rings = rdMolDescriptors.CalcNumRings(m)
    if rings != 0:
        raise SystemExit(
            f"{name}: {rings} ring(s). Every molecule this module builds is "
            "acyclic; a ring here means a ring-closure digit was introduced "
            "by string composition, which is the exact failure the fragment-"
            "assembly rule exists to prevent.")
    conn = _direction_checks(name, m)
    stereo = Chem.FindMolChiralCenters(m, includeUnassigned=True,
                                       useLegacyImplementation=False)
    unassigned = [c for c in stereo if c[1] == "?"]
    db = [b for b in m.GetBonds()
          if b.GetBondType() == Chem.BondType.DOUBLE
          and b.GetBeginAtom().GetSymbol() == "C"
          and b.GetEndAtom().GetSymbol() == "C"]
    db_geom = [str(b.GetStereo()) for b in db]
    return {
        "connectivity": conn,
        "smiles_canonical": Chem.MolToSmiles(m),
        "formula": rdMolDescriptors.CalcMolFormula(m),
        "exact_mass": round(Descriptors.ExactMolWt(m), 6),
        "n_rings": rings,
        "n_stereocentres": len(stereo),
        "n_unassigned_stereocentres": len(unassigned),
        "stereo_assertion_vacuous": len(stereo) == 0,
        "cc_double_bonds": len(db),
        "cc_double_bond_geometry": db_geom,
        "inchikey": Chem.MolToInchiKey(m),
        "skeleton": Chem.MolToInchiKey(m).split("-")[0],
    }


ROWS = [
    # lipid_name, source_figure, note
    ("cKK-E12", "", "NOT DRAWN ANYWHERE IN THE PAPER. The ionizable lipid of "
     "all 144 formulations is named once, with a vendor and catalogue number "
     "(Organix Inc., O-8744), and no structure is shown in any main or "
     "supplementary figure. Empty reference_SMILES is deliberate: the absence "
     "of an anchor is the finding."),
    ("18:0 DDAB", "Fig S2a row C3; redrawn in Fig 3a beside Cat-LNP",
     "Drawn twice, identically: two C18 saturated chains on a quaternary "
     "nitrogen with two methyls, Br- counterion. Stored as the cation."),
    ("18:0 PC (DSPC)", "Fig S2b row N7; redrawn in Fig S15a beside Neu-LNP",
     "Drawn twice. Fig S2b labels it '18:0 PC' and Fig S15a labels the same "
     "structure 'DSPC'; the main text expands DSPC as 1,2-distearoyl-sn-"
     "glycero-3-phosphocholine, so the two labels agree."),
    ("18:1 PA", "Fig S2c row A1; redrawn in Fig S15c beside An-LNP",
     "Drawn twice. Two 18:1 cis-9 acyl chains and a free phosphate shown as "
     "P(=O)(OH)(O-) with Na+. Stored as the mono-anion without sodium."),
]


def main(argv):
    seed_ester = "--seed-ester" in argv
    if "--seed-guard" in argv:
        src = open(os.path.abspath(__file__)).read()
        src = "import assemble\n" + src
        _ast_guard(_extra_src=src)
        raise SystemExit("make_reference: --seed-guard did NOT fire. The AST "
                         "guard is not testing anything.")
    guard = _ast_guard()
    built = build(seed_ester=seed_ester)
    recs, refused = {}, {}
    for name, smiles in built.items():
        try:
            recs[name] = check(name, smiles)
        except SystemExit as exc:
            if not seed_ester:
                raise
            refused[name] = str(exc)
            from rdkit.Chem import rdMolDescriptors as _rmd
            mm = Chem.MolFromSmiles(smiles)
            print("REFUSED %s: %s" % (name, str(exc).splitlines()[0]))
            print("   seeded formula %s exact %.6f"
                  % (_rmd.CalcMolFormula(mm), Descriptors.ExactMolWt(mm)))
            continue
        print(f"{name:18s} {recs[name]['formula']:14s} "
              f"{recs[name]['exact_mass']:.6f}  rings={recs[name]['n_rings']} "
              f"stereo={recs[name]['n_stereocentres']} "
              f"C=C={recs[name]['cc_double_bond_geometry']}")
    print("AST guard:", guard)

    if seed_ester:
        clean = {n: check(n, s) for n, s in build().items()}
        n = "18:0 PC (DSPC)"
        print("\nSEEDED ESTER REVERSAL on", n)
        mm = Chem.MolFromSmiles(build(seed_ester=True)[n])
        sk = Chem.MolToInchiKey(mm).split("-")[0]
        print("  formula    clean %s   seeded %s   same=%s"
              % (clean[n]["formula"], rdMolDescriptors.CalcMolFormula(mm),
                 clean[n]["formula"] == rdMolDescriptors.CalcMolFormula(mm)))
        print("  exact mass clean %.6f seeded %.6f  delta %.9f"
              % (clean[n]["exact_mass"], Descriptors.ExactMolWt(mm),
                 abs(clean[n]["exact_mass"] - Descriptors.ExactMolWt(mm))))
        print("  skeleton   clean %s seeded %s   DIFFERENT=%s"
              % (clean[n]["skeleton"], sk, clean[n]["skeleton"] != sk))
        print("  gate refused the seeded build:", n in refused)
        sys.exit(1 if (n in refused and clean[n]["skeleton"] != sk) else 0)

    path = os.path.join(HERE, "reference_structures.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        for name, fig, note in ROWS:
            smiles = recs[name]["smiles_canonical"] if name in recs else ""
            w.writerow([name, smiles, fig, "2026-11-10",
                        "curation agent, dahlmanje_36701517", "", note])
    print("wrote", path)


if __name__ == "__main__":
    main(sys.argv[1:])
