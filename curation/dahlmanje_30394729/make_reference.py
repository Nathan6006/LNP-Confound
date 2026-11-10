#!/usr/bin/env python3
"""Build the C10 anchors for curation/dahlmanje_30394729, independently.

    python3 -I make_reference.py          # build, check, write reference_structures.csv
    python3 -I make_reference.py --selftest   # additionally run the seeded defects

THE IONIZABLE LIPID IS NEVER DRAWN. 7C1 is a polydisperse branched-PEI/epoxide
polymer; Fig 3A and Fig 5A draw ONE representative arm with the backbone nitrogen
replaced by a diethylamino group, and Fig S1B draws the synthesis, not the
product. So `7C1` gets a row with an EMPTY reference_SMILES, and the anchors are
the things the paper does draw in full: cholesterol (Fig 3A, Fig 5A), DSPC
(Fig 5A), and that representative arm (flagged as an expected mismatch, because
it is a cartoon of the lipid and not the lipid).

INDEPENDENCE. There is no enumerator in this folder -- 03_enumerated.csv carries
no SMILES, because no structure could be enumerated. What this file must not do
is copy the structures that 01_blocks.csv and 02_corpus_overlap.py take from
provenance/reagent_registry.py, which is what the comparison at the end tests
against. `_ast_guard()` walks THIS MODULE'S OWN AST at run time and refuses if

  * it imports any stage module of this folder, or
  * it opens any file this folder emits, or
  * any string literal in it is byte-identical to a registry SMILES.

A substring scan would flag this docstring; an AST walk does not. The guard is
made to fail on purpose by --selftest, which injects each violation into a temp
copy and requires the refusal.

ROUTES, chosen to share as little as possible with a registry lookup:
  cholesterol  re-traversed by hand from the drawn structure, starting at the
               C3 hydroxyl and going round the A ring the other way, with
               different ring-closure digits. Canonicalisation makes a correct
               re-traversal byte-identical to a correct original -- that is
               expected and is NOT the independence signal (F-037).
  DSPC         fragment assembly: sn-glycero-3-phosphocholine + two stearoyl
               fragments, bonded through atom-mapped attachment atoms by
               curation/lib/assemble.attach. No atom is deleted anywhere near
               the sn-2 stereocentre.
  7C1 arm      fragment assembly: diethylamine + the epoxide-opened C15 chain,
               attached N->CH2. The direction is pinned by the atom map, because
               the other regiochemistry (N on the substituted carbon) has the
               SAME FORMULA AND THE SAME EXACT MASS -- that is the seeded defect.
"""
import ast, csv, json, os, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))

from rdkit import Chem, RDLogger                               # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt                  # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula         # noqa: E402

import assemble                                                # noqa: E402
import reagents                                                # noqa: E402

RDLogger.DisableLog("rdApp.*")

EMITTED = ("01_blocks.csv", "03_enumerated.csv", "04_si_tables.json",
           "02_corpus_overlap.json", "05_values.csv", "06_rows.csv",
           "07_qc.json", "07_qc_checks.json")
STAGE_MODULES = ("01_si_tables", "02_corpus_overlap", "06_checks", "09_seeds",
                 "stage01", "stage02", "stage06")


# ---------------------------------------------------------------- the guard
def _ast_guard(path=None):
    """Walk this module's own AST. Returns the findings; raises on any."""
    path = path or os.path.abspath(__file__)
    with open(path) as f:
        tree = ast.parse(f.read(), filename=path)

    banned_literals = set()
    for name in ("Cholesterol", "DSPC"):
        s = reagents.registered(name)
        if s:
            banned_literals.add(s)

    # The guard's OWN ban lists are string literals naming those files, so a
    # naive walk flags this module for defining what it forbids -- which is what
    # it did on its first run. Exclude exactly those two assignments, by node
    # identity, and nothing else.
    exempt = set()
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id in ("EMITTED", "STAGE_MODULES",
                                                     "SELFTESTS")
                for t in node.targets):
            for sub in ast.walk(node):
                exempt.add(id(sub))

    bad = []
    for node in ast.walk(tree):
        if id(node) in exempt:
            continue
        if isinstance(node, ast.Import):
            for a in node.names:
                if a.name in STAGE_MODULES:
                    bad.append(f"imports stage module {a.name!r}")
        elif isinstance(node, ast.ImportFrom):
            if node.module in STAGE_MODULES:
                bad.append(f"imports from stage module {node.module!r}")
        elif isinstance(node, ast.Constant) and isinstance(node.value, str):
            v = node.value
            if v in banned_literals:
                bad.append("contains a string literal byte-identical to a "
                           "reagent_registry SMILES")
            if v in EMITTED:
                bad.append(f"names an emitted file {v!r}")
    if bad:
        raise RuntimeError("make_reference independence guard FAILED: "
                           + "; ".join(sorted(set(bad))))
    return {"checked": os.path.basename(path), "violations": 0,
            "banned_literals_checked": len(banned_literals),
            "stage_modules_checked": len(STAGE_MODULES),
            "emitted_files_checked": len(EMITTED)}


# ------------------------------------------------------------ the anchors
#: The hand re-traversal's stereochemistry, written from the drawn wedges and
#: then DISCARDED. Kept verbatim as evidence, never emitted. See
#: `cholesterol()` and 08_contradictions.csv.
_DISCARDED_CHOLESTEROL_STEREO = (
    "O[C@@H]9CC[C@]8(C)C(=CC[C@@H]7[C@@H]8CC[C@]6(C)[C@H]7CC"
    "[C@@H]6[C@@H](C)CCCC(C)C)C9")


def cholesterol():
    """Re-traversed from the drawn structure (Fig 3A / Fig 5A), C3-OH first.

    CONNECTIVITY ONLY, on purpose.

    The first attempt wrote all eight centres by hand from the drawn wedges,
    traversing C3-OH first and round the A ring the other way. It produced
    HVYWMOMLDIMFJA-OXMOYTRBSA-N: a DIASTEREOMER of cholesterol with four of the
    eight centres inverted (C3, C10, C13, C20), identical in molecular formula,
    identical in exact mass to six decimals, and identical in the InChIKey
    SKELETON block -- the F-013 failure class exactly, and nothing but the
    stereo block saw it. It is not ent-cholesterol either; the enantiomer is
    HVYWMOMLDIMFJA-HFRBDLMJSA-N.

    The honest repair is NOT to flip tokens until the key matches the registry:
    that is fitting the anchor to the answer it is supposed to test. So the
    anchor is emitted without stereochemistry, it tests CONNECTIVITY, and the
    gap is declared on the row (expected_mismatch = stereo_not_derived).
    """
    flat = "OC9CCC8(C)C(=CCC7C8CCC6(C)C7CCC6C(C)CCCC(C)C)C9"
    m = Chem.MolFromSmiles(flat)
    if m is None:
        raise RuntimeError("cholesterol re-traversal does not parse")
    assemble.check_no_radicals(m, label="cholesterol")
    n_unassigned = len([v for v in assemble.stereocentres(m).values() if v == "?"])
    if n_unassigned != 8:
        raise RuntimeError(f"cholesterol skeleton: expected 8 potential centres, "
                           f"found {n_unassigned}")
    return assemble.canonical(flat)


def dspc():
    """sn-glycero-3-phosphocholine + 2 x stearoyl, by atom-mapped attachment."""
    # Both glycerol oxygens are written as real hydroxyls ([OH:n], not [O:n]):
    # a bracket oxygen with no hydrogen is a RADICAL, and assemble refuses it --
    # which is how this build found its own first bug.
    gpc = "[OH:1]C[C@@H]([OH:3])COP(=O)([O-])OCC[N+](C)(C)C"
    stearoyl = "[CH:2](=O)CCCCCCCCCCCCCCCCC"
    one = assemble.attach(gpc, stearoyl, a_map=1, b_map=2, label="DSPC sn-1")
    # The sn-2 oxygen still carries map 3; attach clears only what it bonds.
    two = assemble.attach(one, stearoyl, a_map=3, b_map=2, label="DSPC sn-2")
    return assemble.canonical(two)


def c7_arm(markovnikov=False):
    """The arm Fig 3A / Fig 5A draw for 7C1: Et2N-CH2-CH(OH)-C13H27.

    `markovnikov=True` is the SEEDED DEFECT: the amine on the substituted
    carbon instead. Identical formula, identical exact mass, different
    connectivity -- the only thing that sees it is the InChIKey.
    """
    amine = "[N:1](CC)CC"
    if markovnikov:
        chain = "OC[CH:2]CCCCCCCCCCCCC"        # N lands on C2, OH on C1
    else:
        chain = "[CH2:2][CH](O)CCCCCCCCCCCCC"  # N lands on C1 (anti-Markovnikov)
    return assemble.canonical(assemble.attach(amine, chain, a_map=1, b_map=2,
                                              label="7C1 representative arm"))


def describe(smi):
    m = Chem.MolFromSmiles(smi)
    return {"smiles": smi, "formula": CalcMolFormula(m),
            "exact_mass": round(ExactMolWt(m), 6),
            "inchikey": reagents.full_key(smi),
            "skeleton": reagents.skeleton(smi),
            "n_assigned_stereocentres": len(
                [v for v in assemble.stereocentres(m).values() if v != "?"]),
            "n_unassigned_stereocentres": len(
                [v for v in assemble.stereocentres(m).values() if v == "?"])}


# ------------------------------------------------------------- the selftest
SELFTESTS = {
    "import_a_stage_module": ("import ast, csv, json, os",
                              "import ast, csv, json, os\nimport stage06"),
    "paste_a_registry_smiles": ("RDLogger.DisableLog(\"rdApp.*\")",
                                "RDLogger.DisableLog(\"rdApp.*\")\n"
                                "_PASTED = \"%s\""),
    "open_an_emitted_file": ("RDLogger.DisableLog(\"rdApp.*\")",
                             "RDLogger.DisableLog(\"rdApp.*\")\n"
                             "_F = \"04_si_tables.json\""),
}


def selftest():
    src = open(os.path.abspath(__file__)).read()
    out = {}
    for name, (old, new) in SELFTESTS.items():
        if name == "paste_a_registry_smiles":
            new = new % reagents.registered("Cholesterol")
        mutated = src.replace(old, new, 1)
        if mutated == src:
            out[name] = {"seeded": False, "why": "the anchor text was not found",
                         "guard_fired": False, "credited": False}
            continue
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "mutant.py")
            with open(p, "w") as f:
                f.write(mutated)
            # Run the guard against the MUTATED FILE, in a child process.
            drv = os.path.join(td, "drive.py")
            with open(drv, "w") as f:
                f.write(
                    "import importlib.util, sys, json\n"
                    f"sys.path.insert(0, {os.path.join(ROOT, 'curation', 'lib')!r})\n"
                    f"spec = importlib.util.spec_from_file_location('mr', {os.path.abspath(__file__)!r})\n"
                    "m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)\n"
                    "try:\n"
                    f"    r = m._ast_guard({p!r})\n"
                    "    print('METRIC guard_violations 0')\n"
                    "    sys.exit(0)\n"
                    "except RuntimeError as e:\n"
                    "    print('METRIC guard_violations 1')\n"
                    "    print('MESSAGE', e)\n"
                    "    sys.exit(3)\n")
            r = subprocess.run([sys.executable, "-I", drv],
                               capture_output=True, text=True)
            has_metric = "METRIC guard_violations" in r.stdout
            fired = "METRIC guard_violations 1" in r.stdout
            out[name] = {
                "seeded": True, "metric_line_present": has_metric,
                "guard_fired": fired, "returncode": r.returncode,
                "message": next((l[8:] for l in r.stdout.splitlines()
                                 if l.startswith("MESSAGE ")), None),
                "credited": bool(has_metric and fired and r.returncode == 3),
            }
    # And the control: the real file must pass.
    try:
        _ast_guard()
        out["control_unmutated"] = {"guard_fired": False, "credited": True}
    except RuntimeError as e:
        out["control_unmutated"] = {"guard_fired": True, "credited": False,
                                    "message": str(e)}

    # The structural seed: same formula, same exact mass, different connectivity.
    good, bad = describe(c7_arm()), describe(c7_arm(markovnikov=True))
    out["structural_seed_regiochemistry"] = {
        "correct": good, "seeded": bad,
        "same_formula": good["formula"] == bad["formula"],
        "same_exact_mass": good["exact_mass"] == bad["exact_mass"],
        "mass_difference_Da": round(abs(good["exact_mass"] - bad["exact_mass"]), 9),
        "same_inchikey": good["inchikey"] == bad["inchikey"],
        "same_skeleton_block": good["skeleton"] == bad["skeleton"],
        "credited": (good["formula"] == bad["formula"]
                     and good["exact_mass"] == bad["exact_mass"]
                     and good["inchikey"] != bad["inchikey"]),
    }
    return out


def main():
    guard = _ast_guard()

    chol = cholesterol()
    ds = dspc()
    arm = c7_arm()

    rows = [["lipid_name", "reference_SMILES", "source_figure", "read_on",
             "read_by", "expected_mismatch", "note"]]
    rows.append([
        "7C1", "", "Fig S1B (synthesis scheme); Fig 3A / Fig 5A (one arm)",
        "2026-11-10", "curation agent (claude)", "no_structure_is_drawn",
        "The ionizable lipid is NEVER drawn as a complete structure. It is "
        "branched PEI Mn~600 (1 eq) ring-opened with 14 eq of a C15 1,2-epoxide, "
        "i.e. a polydisperse statistically substituted polymer: no single SMILES "
        "is correct, so the anchor is deliberately EMPTY and the gap is visible."])
    rows.append([
        "7C1_representative_arm", arm, "Fig 3A, Fig 5A (highlighted inset)",
        "2026-11-10", "curation agent (claude)",
        "expected_mismatch_drawing_is_a_cartoon_of_the_polymer",
        "Et2N-CH2-CH(OH)-C13H27, the arm as drawn, with the PEI backbone "
        "nitrogen replaced by the drawing's diethylamino group. Built by "
        "atom-mapped attachment so the regiochemistry is pinned: the "
        "Markovnikov isomer has the identical formula and exact mass."])
    rows.append([
        "Cholesterol", chol, "Fig 3A, Fig 5A", "2026-11-10",
        "curation agent (claude)", "stereo_not_independently_derived",
        "CONNECTIVITY ONLY. Re-traversed from the drawn structure starting at "
        "the C3 hydroxyl. The hand-written stereo attempt gave "
        "HVYWMOMLDIMFJA-OXMOYTRBSA-N, a diastereomer with 4 of 8 centres "
        "inverted, identical in formula, exact mass and skeleton block; it was "
        "discarded rather than tuned against the registry."])
    rows.append([
        "DSPC", ds, "Fig 5A", "2026-11-10", "curation agent (claude)", "",
        "Built by fragment assembly from sn-glycero-3-phosphocholine and two "
        "stearoyl fragments through atom-mapped attachment; no atom deleted."])
    with open(os.path.join(HERE, "reference_structures.csv"), "w", newline="") as f:
        csv.writer(f).writerows(rows)

    report = {
        "_what": "C10 anchors for curation/dahlmanje_30394729, built independently "
                 "of 01_blocks.csv and of provenance/reagent_registry.py. "
                 "provenance/structure_check.py was NOT run from this folder (it "
                 "writes a tracked file other agents read); the comparison below "
                 "is done here, in-process, and recorded by hand in PROVENANCE.md.",
        "independence_guard": guard,
        "anchors": {},
    }
    for name, smi in (("Cholesterol", chol), ("DSPC", ds),
                      ("7C1_representative_arm", arm)):
        d = describe(smi)
        reg = reagents.registered(name)
        d["registry_smiles"] = reg
        d["byte_identical_to_registry"] = (reg == smi)
        d["inchikey_matches_registry"] = (
            reagents.full_key(reg) == d["inchikey"] if reg else None)
        report["anchors"][name] = d
    report["anchors"]["7C1"] = {"reference_SMILES": "", "why_empty":
                                "the paper draws no complete structure"}
    disc = _DISCARDED_CHOLESTEROL_STEREO
    report["discarded_cholesterol_stereo_attempt"] = {
        "smiles": disc,
        "inchikey": reagents.full_key(disc),
        "skeleton": reagents.skeleton(disc),
        "registry_inchikey": reagents.full_key(reagents.registered("Cholesterol")),
        "same_skeleton_block": (reagents.skeleton(disc)
                                == reagents.skeleton(reagents.registered("Cholesterol"))),
        "same_full_inchikey": (reagents.full_key(disc)
                               == reagents.full_key(reagents.registered("Cholesterol"))),
        "formula_and_exact_mass_identical": True,
        "verdict": "DISCARDED, not tuned. Recorded so the hazard is visible."}
    report["note_on_byte_identity"] = (
        "DSPC's assembled anchor canonicalises BYTE-IDENTICAL to the registry "
        "string. That is not evidence of copying and not evidence of "
        "independence either (F-037): a correct re-derivation is byte-identical "
        "to a correct original once canonicalised, which is exactly why C10's "
        "anchor_may_be_a_copy test missed F-037. Independence here is the AST "
        "guard and the build route, both executed above, plus the fact that the "
        "assembly route refused two of its own first attempts (a radical "
        "oxygen, and the cholesterol stereo).")

    if "--selftest" in sys.argv:
        report["selftest"] = selftest()
        n = sum(1 for v in report["selftest"].values()
                if isinstance(v, dict) and v.get("credited"))
        report["selftest_credited"] = f"{n}/{len(report['selftest'])}"

    with open(os.path.join(HERE, "10_anchor_report.json"), "w") as f:
        json.dump(report, f, indent=1)
    print(json.dumps(report, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
