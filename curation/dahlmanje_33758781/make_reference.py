#!/usr/bin/env python3
"""ZG_2020 -- build the C10 anchors by a route independent of the enumerator.

Run:  python3 -I make_reference.py   ->  reference_structures.csv

WHY THIS FILE LOOKS PARANOID. C10 only tests anything if the anchor was built
by a different code path, and that has been claimed twice in this project and
been false twice: F-013 (a byte copy with a note claiming fragment assembly)
and F-037 (the file written to FIX F-013, whose template was the enumeration's
canonical SMILES with the sterol swapped for '%' and one ring digit
renumbered). BYTE-IDENTITY OF THE EMITTED ANCHOR IS NOT THE SIGNAL, because a
correct re-derivation canonicalises to the same string as a correct
enumeration -- which is exactly why C10's own `anchor_may_be_a_copy` missed
F-037. So independence is tested on the GENERATOR, two ways:

  GUARD 1 (AST, at run time): this module's own syntax tree must contain no
    call to an RDKit editing primitive and no import of `assemble`. It walks
    Call/Attribute/Name/Import nodes ONLY -- never string constants -- because
    a guard whose ban list is itself matchable flags its own source (a sibling
    hit exactly that today).
  GUARD 2 (template provenance): no string literal of 25+ characters in this
    module may be a substring of any enumerated SMILES, or vice versa. That is
    F-037's defect stated directly, and it is what byte-comparison of the
    OUTPUT cannot see.

BOTH GUARDS ARE MADE TO FAIL ON PURPOSE by `--selftest`, which writes mutated
copies of this file into a temp directory and runs each in its own child
process. A guard that has only ever passed is not a guard.

HOW THE ANCHOR IS BUILT DIFFERENTLY. The enumerator starts from intermediate 4
and esterifies with an acid through `assemble.esterify` (RDKit RWMol: combine,
AddBond, RemoveAtom), and RDKit canonicalises choline-first. This file writes
the molecule as plain text in the OPPOSITE traversal -- adamantyl cage first,
then the acetyl ester, the glycerol, the fatty acyl, and the phosphocholine
last -- reading off the structure Figure 1a/1b draws. No RDKit editing call
appears anywhere.

THE STEREOCENTRE IS SELECTED, NOT GUESSED. `_pick_R` builds BOTH tags, labels
each with RDKit's modern CIP labeller, and returns the one labelled R --
R being the descriptor the paper PRINTS for intermediate 4 and for all ten
products. That is a deterministic choice from an external fact, not a human
flipping tokens until something matched (the trap that produced a cholesterol
diastereomer for a sibling today).
"""
import ast, csv, json, os, subprocess, sys

from rdkit import Chem
from rdkit.Chem import rdCIPLabeler
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))

# Outputs go to OUTDIR, which is HERE for a normal run. The seed harnesses
# redirect it to a temp directory: a MUTATED copy must never overwrite the
# real curation outputs. It did, and the determinism check caught it --
# reference_structures.csv was left holding seed G3's REGIOISOMER anchors.
OUTDIR = HERE
ME = os.path.abspath(__file__)

BANNED_CALLS = {
    "RWMol", "EditableMol", "AddBond", "RemoveAtom", "RemoveBond",
    "CombineMols", "ReplaceCore", "ReplaceSidechains", "DeleteSubstructs",
    "ReactionFromSmarts", "RunReactants", "esterify", "attach",
    "SetChiralTag", "SetAtomMapNum",
}
BANNED_IMPORTS = {"assemble"}
LITERAL_MIN = 25

# --- the anchors, written adamantyl-cage-first off Figure 1a/1b -----------
# {} is the glycerol C2 chiral tag; _pick_R fills it.
ADAMANTYL_FIRST = ("C12(CC3CC(C1)CC(C3)C2)CC(=O)OC[C{}H]({acyl})"
                   "COP(=O)([O-])OCC[N+](C)(C)C")
ANCHORS = {
    # name: (acyloxy group written acid-end-first, source panel, what it is)
    "A-9":  ("OC(=O)CCCCCCCCC",
             "Figure 1a/1b + Appendix S1 p.2",
             "decanoyl (C10) tail, Fig 1b tail ID 9"),
    "A-17": ("OC(=O)CCCCCCCCCCCCCCCCC",
             "Figure 1a/1b + Appendix S1 p.7",
             "stearoyl (C18) tail, Fig 1b tail ID 17"),
}

# Anchors we deliberately do NOT write. Each says what is missing, so the
# ABSENCE is visible rather than silent (CLAUDE.md section 6).
GAPS = {
    "cKK-E15": (
        "Figure 2a",
        "DRAWN, BUT NOT READABLE. Fig 2a draws the full connectivity -- a "
        "2,5-diketopiperazine, two (CH2)4 linkers to two tertiary nitrogens, "
        "each bearing two CH2-CH(OH)-R arms -- but abbreviates every tail to "
        "a 'CnH2n+1' label whose subscripts are 2-3 px tall in the only "
        "figure we can obtain (709 px web JPEG; the publisher's full-"
        "resolution figure is behind a browser we are not using). n is "
        "therefore not readable, the paper cites no source for cKK-E15 and "
        "mentions it exactly twice, so the tail length is an inference from "
        "the cKK-E12 naming convention and NOT an anchor. See DECISIONS.md "
        "D-2."),
    "beta-Sitosterol": (
        "Figure 2a",
        "NOT WRITTEN ON PURPOSE. It is cholesterol plus a C24 ethyl, and that "
        "ethyl CREATES a stereocentre (24R in natural sitosterol) that a "
        "709 px drawing cannot resolve. Hand-writing a sterol is how a "
        "sibling produced a 4-of-8-inverted diastereomer today with identical "
        "formula, identical exact mass and an identical skeleton InChIKey "
        "block. Not in provenance/reagent_registry.py either. Declared "
        "unpinned rather than guessed."),
    "20alpha-OH Cholesterol": (
        "Figure 2a",
        "NOT WRITTEN ON PURPOSE, same reason: the 20-hydroxyl's "
        "configuration is the whole point of the name and is not resolvable "
        "at this resolution. Not in the reagent registry."),
    "C14PEG2000": (
        "Figure 2a",
        "NOT ANCHORABLE IN PRINCIPLE. The PEG block length is a "
        "distribution, so no single SMILES is right -- which is why "
        "provenance/reagent_registry.py deliberately holds no PEG-lipid."),
}


# ---------------------------------------------------------------------------
def guard_ast(path=ME):
    """GUARD 1. Walk THIS module's AST; refuse any construction primitive."""
    tree = ast.parse(open(path).read())
    hits = []
    for node in ast.walk(tree):
        # Only call targets and imports. String constants are never inspected,
        # so the ban list above cannot match itself.
        if isinstance(node, ast.Call):
            f = node.func
            nm = (f.attr if isinstance(f, ast.Attribute)
                  else f.id if isinstance(f, ast.Name) else None)
            if nm in BANNED_CALLS:
                hits.append(("call", nm, getattr(node, "lineno", -1)))
        elif isinstance(node, ast.Import):
            for a in node.names:
                if a.name.split(".")[0] in BANNED_IMPORTS:
                    hits.append(("import", a.name, node.lineno))
        elif isinstance(node, ast.ImportFrom):
            if (node.module or "").split(".")[0] in BANNED_IMPORTS:
                hits.append(("import_from", node.module, node.lineno))
            for a in node.names:
                if a.name in BANNED_IMPORTS or a.name in BANNED_CALLS:
                    hits.append(("import_name", a.name, node.lineno))
    return hits


def guard_template_provenance(enumerated, path=ME):
    """GUARD 2. No long literal here may be a piece of the enumeration.

    This is F-037 stated directly. It works on the GENERATOR's text, which is
    the only place the defect is visible: the emitted anchor of a CORRECT
    re-derivation is byte-identical to a correct enumeration.
    """
    tree = ast.parse(open(path).read())
    hits = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            lit = node.value
            if len(lit) < LITERAL_MIN:
                continue
            for name, smi in enumerated.items():
                if lit in smi or smi in lit:
                    hits.append((name, lit[:60], node.lineno))
    return hits


def _pick_R(template, **fmt):
    """Return the enantiomer RDKit's modern CIP labeller calls R.

    R is what the paper PRINTS for intermediate 4 and all ten products. The
    choice is made from that external descriptor, mechanically.
    """
    out = []
    for tag in ("@", "@@"):
        smi = template.format(tag, **fmt)
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise SystemExit(f"anchor template did not parse: {smi}")
        rdCIPLabeler.AssignCIPLabels(m)
        cips = [a.GetPropsAsDict().get("_CIPCode")
                for a in m.GetAtoms() if a.HasProp("_CIPCode")]
        out.append((tag, cips, Chem.MolToSmiles(m)))
    rs = [o for o in out if o[1] == ["R"]]
    if len(rs) != 1:
        raise SystemExit(
            f"_pick_R: expected exactly one R enantiomer, got {out}")
    return rs[0][2], out


def build():
    built = {}
    for name, (acyl, _src, _what) in ANCHORS.items():
        smi, _both = _pick_R(ADAMANTYL_FIRST, acyl=acyl)
        built[name] = smi
    return built


def main(argv):
    enum = json.load(open(os.path.join(HERE, "03_enumerated.json")))
    enumerated = {p["name"]: p["smiles"] for p in enum["products"]}

    if "--selftest" in argv:
        return selftest(enumerated)

    g1 = guard_ast()
    g2 = guard_template_provenance(enumerated)
    print(f"METRIC guard_ast_violations={len(g1)} {g1}")
    print(f"METRIC guard_template_provenance_violations={len(g2)} {g2}")
    if g1 or g2:
        raise SystemExit(f"ANCHOR INDEPENDENCE GUARD FAILED: ast={g1} "
                         f"template={g2}")

    built = build()

    rows, agree = [], 0
    for name, (acyl, src, what) in ANCHORS.items():
        anchor = built[name]
        enum_smi = enumerated[name]
        ak = Chem.MolToInchiKey(Chem.MolFromSmiles(anchor))
        ek = Chem.MolToInchiKey(Chem.MolFromSmiles(enum_smi))
        same = ak == ek
        agree += same
        print(f"ANCHOR {name:8s} {what}")
        print(f"   anchor      {anchor}")
        print(f"   enumeration {enum_smi}")
        print(f"   InChIKey    anchor={ak}  enum={ek}  MATCH={same}")
        print(f"   byte-identical strings: {anchor == enum_smi}  "
              "(NOT the signal either way -- see the module docstring)")
        rows.append({
            "lipid_name": name, "reference_SMILES": anchor,
            "source_figure": src, "read_on": "2026-11-10",
            "read_by": "curation agent, from the Europe PMC figure JPEG and "
                       "the Appendix S1 text layer",
            "expected_mismatch": "",
            "note": f"{what}. Built adamantyl-cage-first as plain text, the "
                    f"OPPOSITE traversal to the enumerator's "
                    f"assemble.esterify route; guards 1 and 2 in "
                    f"make_reference.py are executed at run time and are "
                    f"made to fail on purpose by --selftest. InChIKey "
                    f"{'agrees with' if same else 'DISAGREES WITH'} the "
                    f"enumeration.",
        })

    rows.append({
        "lipid_name": "Cholesterol", "reference_SMILES":
            "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)"
            "[C@H]3CC[C@]12C",
        "source_figure": "Figure 2a", "read_on": "2026-11-10",
        "read_by": "provenance/reagent_registry.py, taken verbatim",
        "expected_mismatch": "",
        "note": "NOT AN INDEPENDENT DERIVATION, and says so. This string is "
                "the reagent registry's pin, copied verbatim so the "
                "provenance claim is mechanically true. The paper draws "
                "cholesterol in Fig 2a and the SKELETON is legible at 709 px "
                "(four fused rings, 3-OH, 5-ene, iso-octyl tail) and agrees; "
                "the EIGHT stereocentres are not resolvable at that "
                "resolution, so they are inherited from the registry and are "
                "not evidence from this paper.",
    })
    for name, (src, why) in GAPS.items():
        rows.append({"lipid_name": name, "reference_SMILES": "",
                     "source_figure": src, "read_on": "2026-11-10",
                     "read_by": "curation agent", "expected_mismatch": "",
                     "note": why})

    out = os.path.join(OUTDIR, "reference_structures.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "lipid_name", "reference_SMILES", "source_figure", "read_on",
            "read_by", "expected_mismatch", "note"])
        w.writeheader()
        w.writerows(rows)

    print(f"\nMETRIC anchors_derived={len(ANCHORS)} "
          f"inchikey_agreements={agree}/{len(ANCHORS)}")
    print(f"METRIC anchors_declared_absent={len(GAPS)} {sorted(GAPS)}")
    print(f"METRIC registry_pins_not_independent=1 ['Cholesterol']")
    if agree != len(ANCHORS):
        raise SystemExit("GATE FAILED: an anchor disagrees with the "
                         "enumeration by InChIKey")
    return 0


# ---------------------------------------------------------------------------
def selftest(enumerated):
    """Make both guards, and the connectivity gate, fail on purpose."""
    import shutil, tempfile
    tmp = tempfile.mkdtemp(prefix="zg2020_anchor_")
    # The copy lives in a temp directory, so HERE must be pinned to the real
    # folder or it dies on a missing 03_enumerated.json BEFORE reaching the
    # guard -- which the credit rule below correctly refuses to count. ME is
    # deliberately NOT pinned: both guards must read the MUTATED source.
    src = open(ME).read().replace(
        'HERE = os.path.dirname(os.path.abspath(__file__))',
        f'HERE = {HERE!r}', 1).replace('OUTDIR = HERE', f'OUTDIR = {tmp!r}', 1)
    cases = [
        ("G1_banned_call_injected",
         "Insert a call to Chem.RWMol -- an RDKit editing primitive -- into "
         "build(). GUARD 1 must see it in the AST.",
         "    built = {}\n",
         "    built = {}\n    _ = Chem.RWMol(Chem.MolFromSmiles('CC'))\n",
         "guard_ast_violations"),
        ("G2_template_copied_from_enumeration",
         "Replace the adamantyl-first template with the ENUMERATION'S OWN "
         "canonical SMILES for A-9 -- F-037's exact defect, a copied "
         "template. GUARD 2 must see it even though the emitted anchor would "
         "be byte-identical to a correct one.",
         'ADAMANTYL_FIRST = ("C12(CC3CC(C1)CC(C3)C2)CC(=O)OC[C{}H]({acyl})"\n'
         '                   "COP(=O)([O-])OCC[N+](C)(C)C")',
         'ADAMANTYL_FIRST = ("' + enumerated["A-9"].replace(
             "[C@H](OC(=O)CCCCCCCCC)", "[C{}H]({acyl})") + '")',
         "guard_template_provenance_violations"),
        ("G3_regioisomer_identical_mass",
         "Swap the adamantylacetyl and the fatty acyl on the anchor's "
         "glycerol. Formula IDENTICAL, exact mass IDENTICAL to 0.000000 Da. "
         "Only the InChIKey comparison can see it, so this proves the gate "
         "tests CONNECTIVITY and not composition.",
         'ADAMANTYL_FIRST = ("C12(CC3CC(C1)CC(C3)C2)CC(=O)OC[C{}H]({acyl})"\n'
         '                   "COP(=O)([O-])OCC[N+](C)(C)C")',
         'ADAMANTYL_FIRST = ("C12(CC3CC(C1)CC(C3)C2)CC(=O)O[C{}H](CO{acyl_b})"\n'
         '                   "COP(=O)([O-])OCC[N+](C)(C)C")',
         "anchors_derived"),
    ]
    bad = []
    for cid, why, find, repl, metric in cases:
        if find not in src:
            print(f"{cid:38s} SEED DID NOT APPLY. NOT CREDITED.")
            bad.append(cid)
            continue
        mut = src.replace(find, repl, 1)
        if cid == "G3_regioisomer_identical_mass":
            # the swapped template needs the acyl written the other way round
            mut = mut.replace('smi, _both = _pick_R(ADAMANTYL_FIRST, acyl=acyl)',
                              'smi, _both = _pick_R(ADAMANTYL_FIRST, '
                              'acyl_b=acyl.replace("OC(=O)", "C(=O)"))')
        if mut == src:
            print(f"{cid:38s} SEED IS A NO-OP. NOT CREDITED.")
            bad.append(cid)
            continue
        path = os.path.join(tmp, cid + ".py")
        open(path, "w").write(mut)
        r = subprocess.run([sys.executable, "-I", path],
                           capture_output=True, text=True, cwd=HERE,
                           timeout=600)
        out = r.stdout + r.stderr
        reached = f"METRIC {metric}=" in out
        fired = r.returncode != 0
        verdict = ("FIRED" if (fired and reached) else
                   "FIRED_but_metric_absent_NOT_CREDITED" if fired else
                   "DID_NOT_FIRE")
        if verdict != "FIRED":
            bad.append(cid)
        print(f"{cid:38s} rc={r.returncode} VERDICT={verdict}")
        for l in out.splitlines():
            if l.startswith("METRIC") or "GUARD FAILED" in l \
                    or "GATE FAILED" in l or "MATCH=" in l:
                print(f"      {l.strip()[:150]}")
    shutil.rmtree(tmp, ignore_errors=True)
    print(f"\n{len(cases)} guard seeds; {len(bad)} uncredited: {bad}")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
