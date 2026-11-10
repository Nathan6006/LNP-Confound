#!/usr/bin/env python3
"""C10 anchors for maohq_39375194 -- BUILT, not copied.

RUN:  python3 -I make_reference.py                 # write reference_structures.csv
      python3 -I make_reference.py --selftest      # AST + seeded-defect proof
      python3 -I make_reference.py --seed <name>   # emit one seeded variant

WHY THIS FILE LOOKS THE WAY IT DOES
-----------------------------------
C10 compares each structure we emit against the structure the paper DRAWS. It
tests nothing at all if the anchor came from the same place as the thing it
checks. That has been claimed twice in this project and been false twice:
F-013 (the anchor was a byte copy of the enumeration, same sha256, while its
note claimed independent fragment assembly) and F-037 (the file written to FIX
F-013 re-copied the other half -- 141 bytes of the enumeration's canonical
SMILES with the sterol swapped for `%`).

So independence here is verified BY PROVENANCE, at run time, by walking THIS
MODULE'S OWN AST:

  * no import of reagent_registry, curation/lib/reagents, corpus_components or
    assemble -- the sources stage03_structures.py uses;
  * no call to RWMol, EditableMol, AddBond, RemoveAtom, CombineMols,
    ReplaceSubstructs, esterify or attach.

A sibling's first attempt at this substring-scanned the file and flagged its
own docstring, so the scan is an `ast` walk over parsed nodes, never text.

BYTE-IDENTITY IS NOT THE SIGNAL, and testing for it is the trap F-037 names: a
correct re-derivation is byte-identical to a correct original once canonicalised.
These anchors ARE expected to canonicalise to the same strings. What makes them
evidence is that they were composed by a DIFFERENT TRAVERSAL from fragments
written out of chemical names and the drawn depictions:

  stage03 / the registry : TAIL-FIRST. The sn-1 acyl chain is written first and
                           the head group appears in a branch.
  this file              : HEAD-FIRST. Composition starts at the choline
                           nitrogen / ethanolamine nitrogen / head glycerol /
                           phosphate and walks outward to the acyl chains,
                           which are appended last.

WHAT THE STEREOCENTRE COMES FROM. The backbone C2 tag is NOT copied. Both tags
are generated and the one whose CIP code is **R** is selected. R is a
NAME-derived constant: in 1,2-diacyl-sn-glycero-3-phospho lipids the C2
substituent priorities are O-acyl > CH2-O-P > CH2-O-acyl > H for every head
group here (PC, PE, PG, PA all attach through phosphorus), and the natural sn
configuration puts that at R. The registry arrives at its tag from its own
written-from-name SMILES. Agreement is therefore a test, not an input.

WHAT IS NOT ANCHORED, AND WHY -- stated so the gap is visible rather than
absent: Cholesterol is recorded as a COPY of the registry string. Hand-writing
a 8-stereocentre steroid here would be a worse anchor than none, and
cholesterol already has a dedicated corpus gate (C14 cleaves every ester in the
corpus and pins the liberated sterol's FULL InChIKey). DMG-PEG2000 is stored
with an EMPTY reference_SMILES because no single SMILES is the reagent.
"""
import ast
import csv
import sys

from rdkit import Chem
from rdkit import RDLogger
from rdkit.Chem import Descriptors, rdCIPLabeler, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

# ---------------------------------------------------------------- fragments
# Each is written out of a CHEMICAL NAME and checked against the depiction in
# SI Table 2 (renders/p2_x6*.png, extracted at native 220-336 dpi).
ACYL = {
    # octadecanoyl; 1 carbonyl C + 17 chain C = 18
    "stearoyl": "C(=O)CCCCCCCCCCCCCCCCC",
    # tetradecanoyl; 1 + 13 = 14
    "myristoyl": "C(=O)CCCCCCCCCCCCC",
    # (9Z)-octadec-9-enoyl; 1 + C2-C8 (7) + C9=C10 (2) + C11-C18 (8) = 18
    "oleoyl": "C(=O)CCCCCCC/C=C\\CCCCCCCC",
    # 4-(dimethylamino)butanoyl
    "dmab": "C(=O)CCCN(C)C",
}
ALKYL = {"octadecyl": "C" * 18}

# heptatriaconta-6,9,28,31-tetraen-19-yl, split at the C19 carbinol so the two
# halves are appended on either side of the attachment carbon.
#
# A FRAGMENT HAS A DIRECTION THAT MASS CANNOT SEE, and this file was written
# with that error in it first. The left half was written C1->C18 and then
# appended as a branch off C19, which bonds its FIRST atom -- C1 -- to C19.
# The result had its unsaturation at Delta-6/Delta-9 measured from the ester
# instead of from C1: the correct molecular formula C43H79NO2, the correct
# exact mass 641.611081 to six decimals, and a different InChIKey
# (RBLCKFZMIVYNAG vs CLYUJYGYPALXTO). Only the anchor comparison saw it.
#
# The fix is a derivation, not a fit. From C19 toward C1 the chain runs
# C18..C11 (8), C10=C9, C8, C7=C6, C5..C1 (5); from C19 toward C37 it runs
# C20..C27 (8), C28=C29, C30, C31=C32, C33..C37 (5). The name is symmetric
# about C19, so BOTH branches are the same string, written outward from the
# attachment carbon. `X/C=C\X` is the cis pattern and is unchanged by
# reversing the traversal, so the Z descriptors carry over.
MC3_BRANCH = "C" * 8 + "/C=C\\" + "C" + "/C=C\\" + "CCCCC"
MC3_LEFT = MC3_BRANCH       # C18 -> C1, written outward from C19
MC3_RIGHT = MC3_BRANCH      # C20 -> C37, written outward from C19

HEADS = {
    # phosphocholine, written from the quaternary nitrogen outward
    "PC": "C[N+](C)(C)CCOP(=O)([O-])O",
    # phosphoethanolamine, from the ammonium nitrogen outward
    "PE": "[NH3+]CCOP(=O)([O-])O",
    # phospho-(1'-rac-glycerol): head glycerol first, NO stereocentre, which is
    # how SI Table 2 draws it (plain OH, no wedge) and how Avanti sells PG
    "PG_rac": "OCC(O)COP(=O)([O-])O",
    # phosphate monoester anion
    "PA": "OP(=O)([O-])O",
}

CIP_BACKBONE = "R"   # name-derived; see the module docstring


def _glycerol(head, acyl_sn1, acyl_sn2, tag):
    """HEAD-FIRST: <head> - C3H2 - C2(H)(O-acyl2) - C1H2 - O-acyl1."""
    return "%sC[C%sH](COC(=O)%s)O%s" % (
        head, tag, acyl_sn1.split("C(=O)", 1)[1], acyl_sn2)


def phospholipid(head, acyl_sn1, acyl_sn2, label):
    """Compose, then SELECT the tag by CIP code rather than copying one."""
    cands = []
    for tag in ("@", "@@"):
        smi = _glycerol(HEADS[head], ACYL[acyl_sn1], ACYL[acyl_sn2], tag)
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        rdCIPLabeler.AssignCIPLabels(m)
        codes = [a.GetPropsAsDict().get("_CIPCode")
                 for a in m.GetAtoms()
                 if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED]
        cands.append((tag, smi, codes))
    keep = [c for c in cands if c[2] == [CIP_BACKBONE]]
    if len(keep) != 1:
        raise SystemExit(
            "%s: expected exactly one tag giving CIP %s at the single "
            "backbone centre, got %s" % (label, CIP_BACKBONE,
                                         [(c[0], c[2]) for c in cands]))
    return Chem.MolToSmiles(Chem.MolFromSmiles(keep[0][1]))


def dotap():
    """1,2-dioleoyl-3-trimethylammonium-propane. Written from the quaternary
    nitrogen outward. The stereocentre is left UNSPECIFIED to match the
    reagent registry's convention for this reagent; SI Table 2 draws a wedge
    there, which is logged in 08_contradictions.csv rather than acted on."""
    o = ACYL["oleoyl"]
    return Chem.MolToSmiles(Chem.MolFromSmiles(
        "C[N+](C)(C)CC(O%s)CO%s" % (o, o)))


def ddab():
    """Dimethyldioctadecylammonium, from the nitrogen outward. No counterion,
    matching the registry; SI Table 2 draws the bromide."""
    a = ALKYL["octadecyl"]
    return Chem.MolToSmiles(Chem.MolFromSmiles(
        "C[N+](C)(%s)%s" % (a, a)))


def mc3(left=None, right=None):
    """(6Z,9Z,28Z,31Z)-heptatriaconta-6,9,28,31-tetraen-19-yl
    4-(dimethylamino)butanoate. Written HEAD-FIRST: the dimethylamino nitrogen
    is the first atom, then the butanoate, then the ester oxygen, then the C19
    carbinol with its two chain halves."""
    left = MC3_LEFT if left is None else left
    right = MC3_RIGHT if right is None else right
    smi = "CN(C)CCCC(=O)OC(%s)%s" % (left, right)
    return Chem.MolToSmiles(Chem.MolFromSmiles(smi))


def dmg_peg2000(n=45):
    """1,2-dimyristoyl-rac-glycero-3-methoxypolyethylene glycol 2000.

    n COUNTS ETHYLENE-OXIDE UNITS, AND THE BRACKET IS NOT n. SI Table 2 draws
    glycerol-CH2-O-[CH2CH2-O]44-CH2CH2-O-CH3: the subscript is 44, but one
    further CH2CH2-O is drawn EXPLICITLY outside the closing bracket, so the
    molecule has 45. Reading the subscript as the answer gives a structure one
    ethylene oxide short, with a different InChIKey and a 44 Da mass error.
    Measured both ways against the corpus: n=45 reproduces YZ_2022's stored
    PEG_SMILES exactly (NLDCRRFPGOSAHV-UHFFFAOYSA-N), n=44 does not
    (FHMUYYLWQDJUPK). The glycerol is drawn without a wedge, so it is built
    without a stereocentre."""
    m = ACYL["myristoyl"]
    return Chem.MolToSmiles(Chem.MolFromSmiles(
        "COCC" + "OCC" * (n - 1) + "OCC(O%s)CO%s" % (m, m)))


def build(seed=None):
    """The anchors. `seed` injects ONE defect with an identical molecular
    formula and an identical exact mass, to prove the gate can fire."""
    sn1 = {}
    if seed == "mc3_ester_at_C18":
        # move the attachment one carbon along a C37 chain: C18 ester instead
        # of C19. Same atoms, same formula, same exact mass, different molecule.
        out_mc3 = mc3(left=MC3_LEFT[:-1], right="C" + MC3_RIGHT)
    else:
        out_mc3 = mc3()

    if seed == "dspc_1_3_diacyl":
        # phosphocholine moved to C2, acyls to C1 and C3. Same formula and mass.
        s = ACYL["stearoyl"]
        dspc = Chem.MolToSmiles(Chem.MolFromSmiles(
            "C[N+](C)(C)CCOP(=O)([O-])OC(CO%s)CO%s" % (s, s)))
    else:
        dspc = phospholipid("PC", "stearoyl", "stearoyl", "DSPC")

    if seed == "pg_oleoyl_delta10":
        # the SAME fragment written with the double bond one bond further out:
        # a DIRECTION/position error that mass cannot see.
        ACYL["_o10"] = "C(=O)CCCCCCCC/C=C\\CCCCCCC"
        pg = phospholipid("PG_rac", "_o10", "_o10", "18PG")
        del ACYL["_o10"]
    else:
        pg = phospholipid("PG_rac", "oleoyl", "oleoyl", "18PG")

    sn1["DLin-MC3-DMA"] = out_mc3
    sn1["DSPC"] = dspc
    sn1["18PG"] = pg
    sn1["DOPE"] = phospholipid("PE", "oleoyl", "oleoyl", "DOPE")
    sn1["14PA"] = phospholipid("PA", "myristoyl", "myristoyl", "14PA")
    sn1["DOTAP"] = dotap()
    sn1["DDAB"] = ddab()
    return sn1


# ----------------------------------------------------------- independence
FORBIDDEN_IMPORTS = {"reagent_registry", "reagents", "corpus_components",
                     "assemble", "structure_check"}
FORBIDDEN_CALLS = {"RWMol", "EditableMol", "AddBond", "RemoveAtom",
                   "CombineMols", "ReplaceSubstructs", "esterify", "attach",
                   "SetChiralTag"}


def ast_independence(path=None):
    """Walk THIS MODULE'S parsed AST. Not a substring scan: a sibling's first
    version grepped the file and flagged its own docstring."""
    path = path or __file__
    tree = ast.parse(open(path).read())
    bad_imports, bad_calls = [], []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name.split(".")[0] in FORBIDDEN_IMPORTS:
                    bad_imports.append(a.name)
        elif isinstance(node, ast.ImportFrom):
            root = (node.module or "").split(".")[0]
            if root in FORBIDDEN_IMPORTS:
                bad_imports.append(node.module)
            for a in node.names:
                if a.name in FORBIDDEN_IMPORTS:
                    bad_imports.append(a.name)
        elif isinstance(node, ast.Call):
            fn = node.func
            name = (fn.attr if isinstance(fn, ast.Attribute)
                    else getattr(fn, "id", None))
            if name in FORBIDDEN_CALLS:
                bad_calls.append(name)
    return {"forbidden_imports_found": sorted(set(bad_imports)),
            "forbidden_calls_found": sorted(set(bad_calls)),
            "independent": not bad_imports and not bad_calls,
            "checked_file": path,
            "method": "ast.walk over parsed nodes, not a text scan"}


def keyinfo(smi):
    m = Chem.MolFromSmiles(smi)
    return {"inchikey": Chem.MolToInchiKey(m),
            "skeleton": Chem.MolToInchiKey(m).split("-")[0],
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "exact_mass": round(Descriptors.ExactMolWt(m), 6)}


def selftest():
    base = build()
    out = {"ast": ast_independence(), "seeds": []}
    for seed, target in (("mc3_ester_at_C18", "DLin-MC3-DMA"),
                         ("dspc_1_3_diacyl", "DSPC"),
                         ("pg_oleoyl_delta10", "18PG")):
        s = build(seed)
        a, b = keyinfo(base[target]), keyinfo(s[target])
        out["seeds"].append({
            "seed": seed, "anchor": target,
            "formula_identical": a["formula"] == b["formula"],
            "exact_mass_identical": abs(a["exact_mass"] - b["exact_mass"]) < 1e-9,
            "formula": a["formula"],
            "exact_mass": a["exact_mass"],
            "inchikey_clean": a["inchikey"], "inchikey_seeded": b["inchikey"],
            "skeleton_differs": a["skeleton"] != b["skeleton"],
            "GATE_FIRES": a["inchikey"] != b["inchikey"],
            "CREDITED": (a["formula"] == b["formula"]
                         and abs(a["exact_mass"] - b["exact_mass"]) < 1e-9
                         and a["inchikey"] != b["inchikey"]),
            "why_this_seed": {
                "mc3_ester_at_C18": "the ester moved one carbon along the C37 "
                                    "chain -- a regiochemistry error",
                "dspc_1_3_diacyl": "phosphocholine on C2 and acyls on C1/C3 "
                                   "-- a connectivity error",
                "pg_oleoyl_delta10": "the acyl fragment written with its C=C "
                                     "one bond further out -- a fragment "
                                     "DIRECTION/position error, the class "
                                     "that scored 0/8 then 1/8 then 5/8 "
                                     "elsewhere with a correct formula every "
                                     "time",
            }[seed],
        })
    out["n_seeds"] = len(out["seeds"])
    out["n_credited"] = sum(1 for s in out["seeds"] if s["CREDITED"])
    out["all_credited"] = out["n_credited"] == out["n_seeds"]
    return out


ROWS = [
    ("DLin-MC3-DMA", "SI Table 2 (p2, image xref 61)",
     "composed head-first from 4-(dimethylamino)butanoate and the two halves "
     "of heptatriaconta-6,9,28,31-tetraen-19-yl"),
    ("DOTAP", "SI Table 2 (p2, image xref 64)",
     "composed from the quaternary nitrogen outward; stereocentre left "
     "unspecified to match the registry convention, see 08_contradictions.csv"),
    ("DDAB", "SI Table 2 (p2, image xref 65)",
     "composed from the quaternary nitrogen outward; counterion omitted"),
    ("DOPE", "SI Table 2 (p2, image xref 66)",
     "composed head-first from phosphoethanolamine; C2 tag chosen by CIP = R"),
    ("DSPC", "SI Table 2 (p2, image xref 67)",
     "composed head-first from phosphocholine; C2 tag chosen by CIP = R"),
    ("18PG", "SI Table 2 (p2, image xref 68; head magnified in p2_x68_head.png)",
     "composed head-first from phospho-(1'-rac-glycerol) with TWO oleoyl "
     "chains, because the panel draws a cis C=C in both tails and SI Table 3 "
     "prints Double_bonds = 2; head glycerol left unspecified because the "
     "panel draws it without a wedge"),
    ("14PA", "SI Table 2 (p2, image xref 69)",
     "composed head-first from the phosphate monoester; counterion omitted"),
]


def main():
    if "--selftest" in sys.argv:
        import json
        print(json.dumps(selftest(), indent=1))
        t = selftest()
        return 0 if (t["ast"]["independent"] and t["all_credited"]) else 1
    if "--seed" in sys.argv:
        import json
        print(json.dumps(build(sys.argv[sys.argv.index("--seed") + 1]),
                         indent=1))
        return 0

    built = build()
    ind = ast_independence()
    if not ind["independent"]:
        raise SystemExit("make_reference: this module is not independent: %s"
                         % ind)
    with open("reference_structures.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        for name, fig, how in ROWS:
            w.writerow([name, built[name], fig, "2026-11-10",
                        "curation agent (maohq_39375194)", "",
                        "DERIVED, not copied: %s. Independence verified at run "
                        "time by walking this module's AST (no registry/"
                        "assembler import, no RWMol/AddBond/RemoveAtom/"
                        "CombineMols call); byte-identity with the emitted "
                        "structure is EXPECTED and is not the signal (F-037)."
                        % how])
        w.writerow(["Cholesterol", "", "SI Table 2 (p2, image xref 62)",
                    "2026-11-10", "curation agent (maohq_39375194)", "",
                    "NO INDEPENDENT ANCHOR. Left empty on purpose so the gap "
                    "is visible. Hand-writing an 8-stereocentre steroid here "
                    "would be a worse anchor than none, and cholesterol is "
                    "already covered corpus-wide by C14, which cleaves every "
                    "ester in IL_SMILES and pins the liberated sterol's FULL "
                    "InChIKey against the registry."])
        w.writerow(["DMG-PEG2000", "", "SI Table 2 (p2, image xref 63)",
                    "2026-11-10", "curation agent (maohq_39375194)", "",
                    "EMPTY ON PURPOSE. The panel draws n = 44, but the block "
                    "length is a distribution, so no single SMILES is the "
                    "reagent and the reagent registry excludes PEG lipids on "
                    "purpose (CLAUDE.md section 3). The DRAWN member is "
                    "dmg_peg2000(45) -- 45 ethylene oxides, not the 44 in the "
                    "bracket, because one more is drawn outside it -- and it "
                    "reproduces YZ_2022's stored PEG_SMILES exactly, which is "
                    "what 06_rows.csv emits for corpus consistency."])
    print("wrote reference_structures.csv (%d built + 2 deliberately empty)"
          % len(ROWS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
