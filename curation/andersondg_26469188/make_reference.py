#!/usr/bin/env python3
"""KK_2015 -- build the C10 reference anchors, independently of any enumerator.

WHAT THIS PAPER DRAWS: NOTHING. There is not one chemical structure in the 14
SI pages (pdfplumber reports 0 image objects on every text page and the 6 image
pages are plots, a flowchart, TEM micrographs and JMP screenshots). The main
text is closed-access and not held. So:

  * C12-200, the only ionizable lipid, HAS NO ANCHOR. Its reference_SMILES is
    EMPTY on purpose, so the absence is visible rather than implicit. The paper
    gives no drawing, no formula, no HRMS, no MW and no synthesis, and the
    weight ratio and molar percentages are independent recipe parameters from
    which no molecular mass is derivable.

  * The three MEASURED helper lipids are anchored instead, because the SI does
    spell their systematic names out in the Table S1/S3/S4 footnotes:
        DS = 1,2-distearyol-sn-glycero-  (sic: 'distearyol', the SI's spelling)
        DO = 1,2-dioleoyl-sn-glycero-    (delta-9 cis unsaturated)
        PC = 3-phosphocholine
        PE = 3-phosphoethanolamine

HOW INDEPENDENT IS THAT, HONESTLY? It is a NAME-CONCORDANCE anchor, not an
external attestation of this paper's chemistry. Deriving a structure from a
systematic name is what provenance/reagent_registry.py already does, so this
agrees with the registry by construction of the same name. It will catch a
transcription error in either place -- and the briefing records a sibling whose
hand-built DSPC came out A KETONE BEARING AN ETHER at the right mass to six
decimal places -- but it cannot catch a name that is wrong FOR THIS PAPER.
PROVENANCE.md says so in those words rather than counting this as a C10 pass.

The stereocentre is NOT copied from anywhere. It is derived: in stereospecific
numbering (sn), 1,2-diacyl-sn-glycero-3-phospholipids carry C2 as (R), and this
module ASSERTS R on the built product. If the assertion fails the build raises;
nothing is flipped until it matches.

Route: `assemble.attach` with atom-mapped attachment points only. No string
concatenation, no RWMol, no RemoveAtom, no 'atom index 0 of the fragment'. The
traversal is HEAD-FIRST (phosphate head onto the glycerol sn-3 oxygen, then the
two acyl carbons), which is the opposite order from writing an acyl-first
lipid SMILES.
"""
import argparse
import ast
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import assemble                                                   # noqa: E402
from rdkit import Chem, RDLogger                                  # noqa: E402
from rdkit.Chem import rdMolDescriptors as rdm                    # noqa: E402

RDLogger.DisableLog("rdApp.*")
OUTDIR = os.environ.get("KK2015_OUTDIR", HERE)

# --------------------------------------------------------------- fragments ---
# Mapped attachment atoms. Each map number is a DECLARED bond site, so the
# regiochemistry is stated, never inferred from an atom index.
#
# sn-glycerol backbone. [O:1] = sn-1, [O:2] = sn-2, [O:3] = sn-3.
# The chiral tag on C2 is a HYPOTHESIS here; the (R) assertion below decides it.
BACKBONE = "[OH:1]C[C@@H]([OH:2])C[OH:3]"

HEADS = {
    # 3-phosphocholine: phosphate + choline. [P:3] bonds to the sn-3 oxygen.
    "PC": ("[P:3](=O)([O-])OCC[N+](C)(C)C", "3-phosphocholine"),
    # 3-phosphoethanolamine: phosphate + ethanolamine, primary ammonium.
    "PE": ("[P:3](=O)([O-])OCC[NH3+]", "3-phosphoethanolamine"),
}

ACYL = {
    # stearoyl, 18:0. [C:4]/[C:5] is the carbonyl carbon that bonds to sn-1/sn-2.
    "S": ("CCCCCCCCCCCCCCCCC[C:%d](=O)", "octadecanoyl (stearoyl), 18:0", 18, 0),
    # oleoyl, 18:1 delta-9 cis. Geometry written explicitly, never left bare.
    "O": ("CCCCCCCC/C=C\\CCCCCCC[C:%d](=O)", "(9Z)-octadec-9-enoyl (oleoyl), 18:1 cis-9", 18, 1),
}

TARGETS = {"DSPC": ("S", "PC"), "DOPE": ("O", "PE"), "DOPC": ("O", "PC")}


def _remap(smi, frm, to):
    """Change an atom map number WITHOUT touching the structure."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise ValueError("unparseable fragment %r" % smi)
    n = 0
    for a in m.GetAtoms():
        if a.GetAtomMapNum() == frm:
            a.SetAtomMapNum(to)
            n += 1
    if n != 1:
        raise ValueError("expected exactly one atom mapped %d in %r, found %d"
                         % (frm, smi, n))
    return Chem.MolToSmiles(m)


def _glycerol_c2_cip(mol):
    """CIP code of the glycerol C2: the sp3 carbon bearing an ester O, a CH2
    carrying an ester O, and a CH2 carrying the phosphate O."""
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    patt = Chem.MolFromSmarts(
        "[CX4;H1](OC(=O)[#6])([CH2]OC(=O)[#6])[CH2]O[PX4]")
    hits = mol.GetSubstructMatches(patt)
    if len(hits) != 1:
        raise AssertionError("glycerol C2 not uniquely located (%d hits)" % len(hits))
    a = mol.GetAtomWithIdx(hits[0][0])
    if not a.HasProp("_CIPCode"):
        raise AssertionError("glycerol C2 carries no CIP code -- the build lost "
                             "its stereochemistry")
    return a.GetProp("_CIPCode")


def build(name, backbone=BACKBONE, maps=(3, 1, 2)):
    """Head-first assembly. Returns (smiles, provenance dict).

    `maps` = (which glycerol oxygen takes the PHOSPHATE, which takes the sn-1
    acyl, which takes the sn-2 acyl). The default (3, 1, 2) is the natural
    1,2-diacyl-sn-glycero-3-phospho arrangement. It is a parameter ONLY so the
    --seed-regioisomer harness can build a same-formula, same-exact-mass
    connectivity defect and prove this builder's gate sees connectivity.
    """
    acyl_key, head_key = TARGETS[name]
    head_smi, head_name = HEADS[head_key]
    acyl_tmpl, acyl_name, _nc, _nu = ACYL[acyl_key]

    m_p, m_a1, m_a2 = maps
    # 1. phosphate head onto the sn-3 oxygen  (O:3 x P:3 -> both maps cleared)
    inter = assemble.attach(backbone, _remap(head_smi, 3, m_p),
                            a_map=m_p, b_map=m_p, label="%s O-P" % name)
    # 2. acyl carbon onto the sn-1 oxygen
    inter = assemble.attach(inter, acyl_tmpl % 4, a_map=m_a1, b_map=4,
                            label="%s ester A" % name)
    # 3. acyl carbon onto the sn-2 oxygen
    prod = assemble.attach(inter, acyl_tmpl % 5, a_map=m_a2, b_map=5,
                           label="%s ester B" % name)

    mol = Chem.MolFromSmiles(prod)
    if mol is None:
        raise AssertionError("%s: product does not parse" % name)
    assemble.check_no_radicals(mol, label=name)

    # --- the derived stereo assertion. NOT copied from any stored structure. --
    cip = _glycerol_c2_cip(mol)
    if cip != "R":
        raise AssertionError(
            "%s: glycerol C2 came out (%s), but stereospecific numbering makes "
            "1,2-diacyl-sn-glycero-3-phospholipids (R). The backbone chiral tag "
            "is wrong. This assertion exists so the tag is decided by the "
            "nomenclature convention, not by flipping it until an InChIKey "
            "matches." % (name, cip))

    # --- double-bond geometry asserted explicitly, never assumed ------------
    n_db = sum(1 for b in mol.GetBonds()
               if b.GetBondType() == Chem.BondType.DOUBLE
               and b.GetBeginAtom().GetSymbol() == "C"
               and b.GetEndAtom().GetSymbol() == "C")
    n_z = sum(1 for b in mol.GetBonds()
              if b.GetStereo() == Chem.BondStereo.STEREOZ)
    n_unassigned = sum(1 for b in mol.GetBonds()
                       if b.GetBondType() == Chem.BondType.DOUBLE
                       and b.GetBeginAtom().GetSymbol() == "C"
                       and b.GetEndAtom().GetSymbol() == "C"
                       and b.GetStereo() == Chem.BondStereo.STEREONONE)
    expected_db = 2 if acyl_key == "O" else 0
    if n_db != expected_db:
        raise AssertionError("%s: %d C=C bonds, expected %d" % (name, n_db, expected_db))
    if n_unassigned:
        raise AssertionError("%s: %d C=C bonds have UNASSIGNED geometry. The SI "
                             "states delta-9 CIS, so a bare C=C is a lost fact, "
                             "not a neutral default." % (name, n_unassigned))
    if n_z != expected_db:
        raise AssertionError("%s: %d Z bonds, expected %d" % (name, n_z, expected_db))

    prov = {
        "head": head_name, "acyl": acyl_name,
        "n_stereocentres": len(assemble.stereocentres(mol)),
        "glycerol_C2_CIP": cip,
        "n_CC_double_bonds": n_db, "n_Z": n_z, "n_unassigned_CC": n_unassigned,
        "formula": rdm.CalcMolFormula(mol),
        "exact_mass": round(rdm.CalcExactMolWt(mol), 6),
        "inchikey": Chem.MolToInchiKey(mol),
    }
    return Chem.MolToSmiles(mol), prov


# ----------------------------------------------------------- the AST guard ---
# F-037: byte-identity of the EMITTED anchor is not the signal, and testing for
# it is a trap -- a correct re-derivation is byte-identical to a correct
# enumeration once canonicalised. So the guard inspects what this module may
# REACH, by walking its own AST at run time and looking at Call, Import and
# ImportFrom nodes ONLY. A substring scan over the source would match the ban
# list inside its own docstring and inside the ban list itself; one sibling's
# guard failed on its own module for exactly that reason.
BANNED_MODULES = {"corpus_components", "reagent_registry", "reagents",
                  "pandas", "stage01_extract", "stage03_corpus", "stage05_emit"}
BANNED_CALLS = {"smi", "load", "read_csv"}


def ast_guard(path=None, banned_modules=None, banned_calls=None):
    """Walk this module's AST. Returns (findings, counts). Never a bare bool."""
    path = path or os.path.abspath(__file__)
    bm = BANNED_MODULES if banned_modules is None else banned_modules
    bc = BANNED_CALLS if banned_calls is None else banned_calls
    tree = ast.parse(open(path).read())
    found, n_import, n_call = [], 0, 0
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for al in node.names:
                n_import += 1
                if al.name.split(".")[0] in bm:
                    found.append(("import", al.name, node.lineno))
        elif isinstance(node, ast.ImportFrom):
            n_import += 1
            root = (node.module or "").split(".")[0]
            if root in bm:
                found.append(("from-import", node.module, node.lineno))
            for al in node.names:
                if al.name in bc:
                    found.append(("from-import-name", al.name, node.lineno))
        elif isinstance(node, ast.Call):
            n_call += 1
            f = node.func
            nm = (f.attr if isinstance(f, ast.Attribute)
                  else f.id if isinstance(f, ast.Name) else None)
            if nm in bc:
                found.append(("call", nm, node.lineno))
    return found, {"import_nodes_examined": n_import,
                   "call_nodes_examined": n_call,
                   "banned_modules": sorted(bm), "banned_calls": sorted(bc)}


def template_guard(path=None, enumerated=None, minlen=25):
    """No string literal of >= minlen chars in this module may be a substring
    of any emitted structure. That is what a copied template looks like."""
    path = path or os.path.abspath(__file__)
    enumerated = enumerated or []
    lits = []
    for node in ast.walk(ast.parse(open(path).read())):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            if len(node.value) >= minlen:
                lits.append(node.value)
    hits = [(li[:60], s[:40]) for li in lits for s in enumerated if li in s]
    return hits, {"literals_examined": len(lits), "minlen": minlen,
                  "structures_examined": len(enumerated)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed-copied-template", action="store_true",
                    help="MAKE THE GUARDS FAIL ON PURPOSE. Imports a banned "
                         "module and defines a literal lifted from an emitted "
                         "structure, then reports whether both guards fire.")
    ap.add_argument("--seed-regioisomer", action="store_true",
                    help="Build 1,3-diacyl-sn-glycero-2-phospho instead. IDENTICAL "
                         "molecular formula and IDENTICAL exact mass to 6 dp, "
                         "different connectivity -- proves the gate is not a mass "
                         "check.")
    ap.add_argument("--seed-wrong-stereo", action="store_true",
                    help="Flip the backbone chiral tag and confirm the derived "
                         "(R) assertion raises.")
    args = ap.parse_args()

    if args.seed_regioisomer:
        good, _ = build("DSPC")
        m_good = Chem.MolFromSmiles(good)
        # The defect, written explicitly HERE IN THE HARNESS ONLY:
        # 1,3-distearoyl-glycero-2-phosphocholine. Same three fragments, same
        # three bonds, different attachment points -> a constitutional isomer.
        bad_smi = ("CCCCCCCCCCCCCCCCCC(=O)OCC(COC(=O)CCCCCCCCCCCCCCCCC)"
                   "OP(=O)([O-])OCC[N+](C)(C)C")
        m_bad = Chem.MolFromSmiles(bad_smi)
        same_f = rdm.CalcMolFormula(m_bad) == rdm.CalcMolFormula(m_good)
        dm = abs(rdm.CalcExactMolWt(m_bad) - rdm.CalcExactMolWt(m_good))
        diff_ik = Chem.MolToInchiKey(m_bad) != Chem.MolToInchiKey(m_good)
        print("SEED regioisomer: 1,3-diacyl-glycero-2-phospho vs 1,2-diacyl-3-phospho")
        print("  formula identical : %s  (%s vs %s)"
              % (same_f, rdm.CalcMolFormula(m_good), rdm.CalcMolFormula(m_bad)))
        print("  exact-mass delta  : %.6f Da" % dm)
        print("  InChIKey differs  : %s" % diff_ik)
        print("  good=%s" % Chem.MolToInchiKey(m_good))
        print("  bad =%s" % Chem.MolToInchiKey(m_bad))
        # now: does the REAL builder refuse to produce it?
        gate, why = None, None
        try:
            build("DSPC", maps=(2, 1, 3))
        except Exception as e:
            gate, why = type(e).__name__, str(e).split(".")[0]
        print("  builder refuses   : %s  (%s: %s)"
              % (gate is not None, gate, why))
        print("  WHICH GATE        : assemble.check_stereo_preserved, not a mass "
              "check. In the 1,3-diacyl-2-phospho isomer the glycerol C2 bears two "
              "IDENTICAL CH2-O-stearoyl arms, so it stops being stereogenic; the "
              "library sees the assigned centre disappear and refuses. That is the "
              "correct chemistry and the correct refusal.")
        ok = same_f and dm < 1e-6 and diff_ik and gate is not None
        print("SEED regioisomer: %s" % (
            "CONNECTIVITY DEFECT REFUSED AT IDENTICAL FORMULA AND IDENTICAL "
            "EXACT MASS (delta %.6f Da)" % dm if ok
            else "seed did not isolate connectivity"))
        print("  SENSITIVITY LIMIT : this builder's OWN assertions (C2=R, C=C "
              "count, Z count) CANNOT see double-bond POSITION. A delta-6 oleoyl "
              "passes every one of them at identical formula and identical exact "
              "mass. Only the three-way InChIKey comparison in stage04 catches "
              "that, and stage04 seeds exactly that defect.")
        return 0 if ok else 1

    if args.seed_wrong_stereo:
        bad = BACKBONE.replace("[C@@H]", "[C@H]")
        assert bad != BACKBONE, "seed is a no-op: the tag did not change"
        try:
            build("DSPC", backbone=bad)
        except AssertionError as e:
            print("SEED wrong-stereo: ASSERTION FIRED -- %s" % str(e).split(".")[0])
            return 0
        print("SEED wrong-stereo: NOT CAUGHT -- the (R) assertion has no teeth")
        return 1

    built = {}
    for name in sorted(TARGETS):
        smi, prov = build(name)
        built[name] = (smi, prov)
        print("%-5s %-14s %s  C2=%s  nZ=%d  %s"
              % (name, prov["formula"], prov["inchikey"], prov["glycerol_C2_CIP"],
                 prov["n_Z"], smi))

    enumerated = [s for s, _ in built.values()]

    if args.seed_copied_template:
        path = os.path.join(OUTDIR, "_seed_make_reference.py")
        stolen = enumerated[0][:40]
        src = (open(os.path.abspath(__file__)).read()
               .replace("import assemble  ", "import reagent_registry  # SEED\nimport assemble  ", 1))
        src += "\n_SEEDED_TEMPLATE = %r\n" % stolen
        open(path, "w").write(src)
        f1, c1 = ast_guard(path)
        f2, c2 = template_guard(path, enumerated)
        print("SEED copied-template: ast_guard findings=%r (examined %s)" % (f1, c1))
        print("SEED copied-template: template_guard hits=%d (examined %s)"
              % (len(f2), c2))
        os.remove(path)
        ok = bool(f1) and bool(f2)
        print("SEED copied-template: %s" % ("BOTH GUARDS FIRED" if ok
                                            else "A GUARD STAYED SILENT -- no teeth"))
        return 0 if ok else 1

    f1, c1 = ast_guard()
    f2, c2 = template_guard(None, enumerated)
    print("ast_guard      : findings=%r  %s" % (f1, c1))
    print("template_guard : hits=%r  %s" % (f2, c2))
    if f1 or f2:
        print("GUARD FAILED on the real module -- not writing an anchor.")
        return 1

    hdr = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
           "read_by", "expected_mismatch", "note"]
    rows = [
        ["C12-200", "", "none -- this paper draws no structure anywhere",
         "2026-11-10", "stage01_extract.py + manual read of all 14 SI pages", "",
         "NO ANCHOR EXISTS. The SI contains zero chemical structure drawings "
         "(0 image objects on every text page; the 6 image pages are plots, a "
         "flowchart, TEM micrographs and JMP screenshots) and prints no formula, "
         "no HRMS, no MW and no synthesis for C12-200. The main text is "
         "closed-access and not held. Empty on purpose so the gap is visible."],
    ]
    for name in ("DSPC", "DOPE", "DOPC"):
        smi, prov = built[name]
        rows.append([
            name, smi,
            "SI Table S1/S3/S4 footnote (systematic name, not a drawing)",
            "2026-11-10", "make_reference.py (head-first assemble.attach route)", "",
            "NAME-CONCORDANCE anchor, NOT an independent attestation of this "
            "paper's chemistry: the paper draws nothing, so this is the SI's own "
            "systematic name re-derived by a different hand. %s + %s; glycerol C2 "
            "asserted (R) from the sn convention, not copied; %d Z double bond(s) "
            "asserted explicitly; %s, exact mass %.6f."
            % (prov["acyl"], prov["head"], prov["n_Z"], prov["formula"],
               prov["exact_mass"])])
    rows.append([
        "DSPE", "", "SI Table S1 footnote (abbreviation only)", "2026-11-10",
        "stage01_extract.py", "",
        "No anchor and no structure: all three DSPE designs (A-07, A-12, A-14) "
        "failed to formulate, so no value exists. 'DSPE' is also an incomplete "
        "specification on its own."])
    rows.append([
        "Cholesterol", "", "SI Table S1 column header ('Chol.')", "2026-11-10",
        "stage01_extract.py", "",
        "Deliberately NOT hand-written. A sibling's hand-built cholesterol came "
        "out a diastereomer with 4 of 8 centres inverted at identical formula, "
        "identical exact mass to six decimals and an identical skeleton InChIKey "
        "block. This paper gives only the abbreviation 'Chol.', so no sterol is "
        "asserted here; the registry pin is also NOT counted as an anchor."])
    rows.append([
        "PEG-lipid", "", "none -- never named in the SI", "2026-11-10",
        "stage01_extract.py", "",
        "The PEG-lipid is never identified: the Table S1 column header is the "
        "bare word 'PEG'. Emitted as PEG_name='Unknown'. Guessing was refused "
        "because C14-PEG2000 is drawn in the literature as a "
        "phosphoethanolamine-carbamate-PEG and stored in this corpus as the "
        "glycerol ether."])

    path = os.path.join(OUTDIR, "reference_structures.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(hdr)
        w.writerows(rows)
    print("wrote %s  (%d rows, %d with a structure, %d empty on purpose)"
          % (path, len(rows), sum(1 for r in rows if r[1]),
             sum(1 for r in rows if not r[1])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
