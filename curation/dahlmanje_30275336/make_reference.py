"""Build this paper's C10 reference anchors, and PROVE the build is independent.

THE SITUATION HERE IS UNUSUAL AND IS RECORDED AS SUCH. This paper emits ZERO
ionizable-lipid structures (03_enumerated.csv), so there is no enumeration for
an ionizable-lipid anchor to test. Per CLAUDE.md section 5 the absence is made
VISIBLE rather than skipped: every one of the 8 compounds gets a row with an
EMPTY reference_SMILES and a note saying what the paper draws instead.

What the paper DOES draw in full is its excipients. Two of them are built here
from the drawing, by hand, and compared by FULL InChIKey against the structure
the corpus already holds. That is a real external check and it is the one that
matters most for this project: F-013 put epicholesterol (3-alpha) on 537 rows,
and the epimers share molecular formula, exact mass to four decimals, the
InChIKey SKELETON block and all 17 descriptor columns. Only the full InChIKey
sees it.

INDEPENDENCE IS VERIFIED BY PROVENANCE, NOT BY BYTE-IDENTITY. Byte-identity of
an emitted anchor is not the signal and testing for it is a trap (F-037): a
correct re-derivation is byte-identical to a correct enumeration once
canonicalised. So this module walks ITS OWN AST at run time and asserts that it
imports no assembler and calls no RWMol/AddBond/RemoveAtom/CombineMols/
esterify/attach. It walks the AST rather than substring-scanning the file,
because a sibling's substring scan flagged its own docstring -- as this
docstring, which names every one of those functions, would also do.
"""
import ast, json, sys, os
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

SEED = sys.argv[1] if len(sys.argv) > 1 else "none"
HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- AST gate
BANNED_CALLS = {"RWMol", "AddBond", "RemoveAtom", "CombineMols",
                "esterify", "attach", "ReplaceSubstructs", "EditableMol"}
BANNED_IMPORTS = {"assemble", "curation.lib.assemble", "enumerate"}

tree = ast.parse(open(os.path.join(HERE, "make_reference.py")).read())
calls, imports = set(), set()
for node in ast.walk(tree):
    if isinstance(node, ast.Call):
        f = node.func
        name = (f.id if isinstance(f, ast.Name)
                else f.attr if isinstance(f, ast.Attribute) else None)
        if name:
            calls.add(name)
    elif isinstance(node, ast.Import):
        for a in node.names:
            imports.add(a.name)
    elif isinstance(node, ast.ImportFrom):
        imports.add(node.module or "")

if SEED == "seed_banned_call":
    calls.add("RWMol")

bad_calls = sorted(calls & BANNED_CALLS)
bad_imports = sorted(imports & BANNED_IMPORTS)
ast_clean = not bad_calls and not bad_imports

# ------------------------------------------------- anchors read from the PDF
# Written by reading what SI Fig S2b (p.7) DRAWS, atom by atom, and typing the
# SMILES out. No corpus string was copied and no assembler was called; the AST
# gate above is the mechanical proof of the second half of that claim.
#
# Cholesterol as drawn: 3-BETA hydroxyl (bold wedge to HO), Delta-5 ring
# alkene, 10-beta and 13-beta angular methyls, 17-beta iso-octyl side chain,
# 8-beta/9-alpha/14-alpha ring-junction H's shown explicitly.
CHOL = "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
# DOPE as drawn: 1,2-di-(cis-9-octadecenoyl)-sn-glycero-3-phosphoethanolamine,
# zwitterion (phosphate O- and ammonium NH3+ both drawn).
DOPE = ("CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[NH3+])"
        "OC(=O)CCCCCCC/C=C\\CCCCCCCC")

if SEED == "seed_epimer":
    # IDENTICAL formula C27H46O and IDENTICAL exact mass; only C3 inverts.
    # This is precisely F-013, and it is invisible to every composition check.
    CHOL = CHOL.replace("C[C@@H](O)CC", "C[C@H](O)CC")

CORPUS = {   # what data/ holds, for the comparison only -- never as the source
  "Cholesterol": "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
  "DOPE": "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[NH3+])OC(=O)CCCCCCC/C=C\\CCCCCCCC",
}

def describe(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    ik = Chem.MolToInchiKey(m)
    return {"canonical": Chem.MolToSmiles(m), "inchikey": ik,
            "skeleton_block": ik.split("-")[0], "formula": CalcMolFormula(m),
            "exact_mass": round(ExactMolWt(m), 4)}

anchors = {}
for name, built in (("Cholesterol", CHOL), ("DOPE", DOPE)):
    mine, theirs = describe(built), describe(CORPUS[name])
    anchors[name] = {
      "built_here": mine, "corpus": theirs,
      "FULL_InChIKey_match": mine["inchikey"] == theirs["inchikey"],
      "skeleton_match": mine["skeleton_block"] == theirs["skeleton_block"],
      "formula_match": mine["formula"] == theirs["formula"],
      "exact_mass_match": mine["exact_mass"] == theirs["exact_mass"],
      "canonical_SMILES_byte_identical": mine["canonical"] == theirs["canonical"],
    }

# ------------------------------------------------------- the anchor CSV
HEADER = "lipid_name,reference_SMILES,source_figure,read_on,read_by,expected_mismatch,note"
COMPOUNDS = ["7C1", "73", "78", "80", "83", "88", "90", "92"]
lines = [HEADER]
for c in COMPOUNDS:
    fig = "SI Fig S1j (p.5)" if c == "92" else "SI Fig S2b (p.7)"
    lines.append(
        '%s,,%s,2026-11-10,curation-agent,no_structure_is_drawn,'
        '"The paper draws a single repeat-unit cartoon - one tertiary nitrogen '
        'bearing one arm with the backbone truncated to two stubs - not a '
        'complete molecule. Degree of substitution is never stated and the '
        'PEI-based compounds are polydisperse polymers. EMPTY ON PURPOSE so '
        'the absence of an anchor is visible (CLAUDE.md section 5)."' % (c, fig))
for name in ("Cholesterol", "DOPE"):
    a = anchors[name]
    lines.append(
        '%s,%s,SI Fig S2b (p.7),2026-11-10,curation-agent,%s,'
        '"Excipient drawn in full. Built here by hand from the drawing; '
        'compared by FULL InChIKey against the corpus. full_match=%s '
        'skeleton_match=%s."'
        % (name, a["built_here"]["canonical"],
           "none" if a["FULL_InChIKey_match"] else "STEREO_MISMATCH",
           a["FULL_InChIKey_match"], a["skeleton_match"]))
if SEED == "none":
    with open(os.path.join(HERE, "reference_structures.csv"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

print(json.dumps({
  "seed": SEED,
  "AST_independence_gate": {
    "method": "ast.walk of this module at run time, not a substring scan "
              "(a sibling's substring scan flagged its own docstring)",
    "calls_seen_n": len(calls), "banned_calls_found": bad_calls,
    "banned_imports_found": bad_imports, "CLEAN": ast_clean},
  "anchors": anchors,
  "anchor_rows_written": len(lines) - 1,
  "ionizable_anchors_with_a_SMILES": 0,
  "ionizable_anchors_empty_on_purpose": len(COMPOUNDS),
}, indent=1))
