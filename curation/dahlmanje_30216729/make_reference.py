#!/usr/bin/env python3
"""Build reference_structures.csv -- the C10 anchor, and an HONEST statement of
its posture.

    python3 -I make_reference.py [OUTDIR]
    python3 -I make_reference.py --selftest      # make the guards fail on purpose

THE SHORT VERSION. This paper's ionizable lipid is never drawn as a structure,
so there is no anchor for it; the two excipients it DOES draw in full
(cholesterol, DSPC) can only be written here from provenance/reagent_registry.py,
which is the same source the corpus itself uses. That is a REGISTRY PIN, NOT AN
INDEPENDENT ANCHOR, and this file declares it as such rather than counting it.
The guards below exist to make that declaration mechanical instead of a claim
in prose -- CLAUDE.md section 6 records C10's independence being asserted twice
and false twice (F-013, a byte copy; F-037, the file written to fix F-013
re-copying the other half).

WHY NO INDEPENDENT ANCHOR WAS BUILT, stated so a reviewer can check each leg:

  1. 7C1 -- the ionizable component. SI p14 draws a GENERIC SCHEME (an epoxide
     labelled R=C13 plus PEI600, 90 C / 72 h / EtOH), not a molecule. PEI600 is
     branched polydisperse polyethylenimine and the degree of substitution is
     not stated, so no atom-resolved structure exists to anchor. Empty
     reference_SMILES, per the rule that the ABSENCE of an anchor must be
     visible.

  2. Cholesterol -- drawn in full on SI p14. The brief's own measured warning
     applies: a hand-written sterol came out a diastereomer with 4 of 8 centres
     inverted, identical in formula, in exact mass to six decimals and in the
     skeleton InChIKey block. Writing one here would be that same gamble, and
     the correct response to a mismatch would be to discard the anchor, not to
     flip tokens until it matched. So the structure is taken from the registry
     and the row is marked NOT INDEPENDENT.

  3. DSPC -- drawn in full on SI p14. curation/lib/assemble.esterify would
     build it safely from stearic acid and sn-glycero-3-phosphocholine, and
     that is a genuinely different code path. But the headgroup fragment
     written from the chemical name, 'COP(=O)([O-])OCC[N+](C)(C)C', is a
     28-character SUBSTRING of the registry's DSPC string, so the
     template-provenance guard below correctly refuses to certify the result as
     independent. Writing the one stereocentre of sn-glycerol by hand is the
     same gamble as the sterol, one centre instead of eight. Marked NOT
     INDEPENDENT.

  4. The four PEG-lipids -- drawn in full, with explicit (OCH2CH2)n subscripts.
     Empty reference_SMILES: the block length is a DISTRIBUTION, which is
     exactly why CLAUDE.md section 3 keeps PEG-lipids out of the registry. What
     IS checkable about them is checked, in stage03_drawn_checks.py: the drawn
     repeat counts 45/45/45/67 against the table's nominal 2000/2000/2000/3000,
     ratio residual 0.74%.

  5. No bond trace was attempted. The structure band is VECTOR -- in the band
     above the library table, 391 filled paths carrying 1,598 segments against
     only 4 stroked paths, so bonds here are FILLED POLYGONS with width=None
     and any trace must pass include_filled=True. A full connectivity trace of
     cholesterol and DSPC was judged not worth it for a paper that emits zero
     rows, and that is recorded as a decision (DECISIONS.md D-4), not disguised
     as an impossibility.

Net C10 result for this folder: ZERO independent anchors, honestly zero.
"""
import ast
import csv
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
OUTDIR = os.path.abspath(sys.argv[-1]) if (
    len(sys.argv) > 1 and not sys.argv[-1].startswith("-")) else HERE
sys.path.insert(0, os.path.join(REPO, "provenance"))
import reagent_registry as rr   # noqa: E402  (the dependency the guard reports)

HEADER = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
          "read_by", "expected_mismatch", "note"]

READ_ON = "2026-11-10"
READ_BY = "curation agent, dahlmanje_30216729"

# Modules/functions whose use means the anchor is NOT independent of the
# project's own stored structures.
NON_INDEPENDENT_SOURCES = {
    "reagent_registry", "corpus_components", "build_dataset",
    "build_new_curation", "enumerate",
}
NON_INDEPENDENT_CALLS = {"smi", "load", "REGISTRY"}


def ast_provenance(path):
    """Walk a module's AST and report which structure sources it uses.

    INSPECTS Call AND Import NODES ONLY. That restriction is the point: a
    guard that greps the file's text matches its own ban list, its own
    docstring and its own comments. One sibling's guard failed on its own
    module for exactly that reason and another's substring scan flagged its own
    docstring, so this never looks at a string literal to decide provenance.
    """
    tree = ast.parse(open(path).read())
    imports, calls = set(), set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for a in node.names:
                imports.add(a.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                imports.add(node.module.split(".")[0])
        elif isinstance(node, ast.Call):
            f = node.func
            if isinstance(f, ast.Name):
                calls.add(f.id)
            elif isinstance(f, ast.Attribute):
                calls.add(f.attr)
    used = sorted(imports & NON_INDEPENDENT_SOURCES)
    # An Attribute LOAD of REGISTRY is not a Call, so catch it separately --
    # still structurally, never textually.
    attrs = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Attribute):
            attrs.add(node.attr)
    return {
        "module": os.path.basename(path),
        "imports_detected": sorted(imports),
        "non_independent_imports": used,
        "non_independent_calls": sorted((calls | attrs) & NON_INDEPENDENT_CALLS),
        "is_independent_of_project_structure_stores": not used,
        "nodes_inspected": "Import, ImportFrom, Call, Attribute only -- never "
                           "string literals, so the ban list cannot match "
                           "itself",
    }


def template_provenance(literals, enumerated):
    """No string literal of >=25 chars may be a substring of any structure this
    curation would emit, or of any registry structure. That is the F-037 test:
    byte-identity of the EMITTED anchor is not the signal, because a correct
    re-derivation canonicalises to the same bytes as a correct enumeration.
    What matters is whether the TEMPLATE was copied."""
    hits = []
    for lit in literals:
        if len(lit) < 25:
            continue
        for name, target in enumerated.items():
            if target and lit in target:
                hits.append({"literal": lit, "found_inside": name,
                             "literal_len": len(lit)})
    return {"n_literals_checked": sum(1 for l in literals if len(l) >= 25),
            "violations": hits, "n_violations": len(hits),
            "clean": not hits}


def module_literals(path):
    tree = ast.parse(open(path).read())
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            out.append(node.value)
    return out


def rows():
    chol = rr.REGISTRY["Cholesterol"][0]
    dspc = rr.REGISTRY["DSPC"][0]
    return [
        {"lipid_name": "7C1",
         "reference_SMILES": "",
         "source_figure": "SI p14, synthesis scheme (R=C13 epoxide + PEI600)",
         "read_on": READ_ON, "read_by": READ_BY,
         "expected_mismatch": "no_structure_is_drawn",
         "note": ("The ionizable component is drawn ONLY as a generic R-group "
                  "scheme. PEI600 is branched polydisperse polyethylenimine, "
                  "nominal MW 600, and the degree of epoxide substitution is "
                  "not stated; the paper gives no formula, no HRMS and no NMR. "
                  "Empty on purpose so the absence of an anchor is visible.")},
        {"lipid_name": "Cholesterol",
         "reference_SMILES": chol,
         "source_figure": "SI p14, drawn in full",
         "read_on": READ_ON, "read_by": READ_BY,
         "expected_mismatch": "anchor_is_a_registry_pin_NOT_independent",
         "note": ("Structure taken from provenance/reagent_registry.py, the "
                  "same source the corpus uses, so this compares a value to "
                  "itself and is NOT counted as a C10 result. Not hand-written: "
                  "a hand-written sterol has already come out a diastereomer "
                  "with 4 of 8 centres inverted, identical in formula, exact "
                  "mass and skeleton InChIKey block. The paper's drawing was "
                  "not traced (see D-4).")},
        {"lipid_name": "DSPC",
         "reference_SMILES": dspc,
         "source_figure": "SI p14, drawn in full",
         "read_on": READ_ON, "read_by": READ_BY,
         "expected_mismatch": "anchor_is_a_registry_pin_NOT_independent",
         "note": ("Structure taken from the registry. An assemble.esterify "
                  "build from stearic acid + sn-glycero-3-phosphocholine would "
                  "be a different code path, but the headgroup fragment "
                  "written from the chemical name is a 28-character substring "
                  "of the registry string, so the template-provenance guard "
                  "would refuse to certify it as independent. NOT counted as a "
                  "C10 result.")},
    ] + [
        {"lipid_name": nm,
         "reference_SMILES": "",
         "source_figure": f"SI p14, drawn in full with (OCH2CH2){n}OCH3",
         "read_on": READ_ON, "read_by": READ_BY,
         "expected_mismatch": "block_length_is_a_distribution",
         "note": (f"PEG-lipid drawn with an explicit repeat count of {n}. Empty "
                  f"reference_SMILES because the polyether block length is a "
                  f"distribution and no single SMILES is right -- the same "
                  f"reason CLAUDE.md section 3 keeps PEG-lipids out of the "
                  f"reagent registry. The drawn repeat count IS checked, "
                  f"against the table's nominal PEG MW, in "
                  f"stage03_drawn_checks.py.")}
        for nm, n in (("C14PEG2000", 45), ("C16PEG2000", 45),
                      ("C18PEG2000", 45), ("C14PEG3000", 67))
    ]


def selftest():
    """MAKE THE GUARDS FAIL ON PURPOSE, in child processes.

    A new check is not done until it has been made to fail (CLAUDE.md section
    6). Two defects are seeded:
      1. a module that imports reagent_registry while claiming independence --
         ast_provenance must report is_independent = False;
      2. a module holding a long literal copied out of a target structure --
         template_provenance must flag it;
      3. a control module that does neither -- both guards must come back clean,
         or the guards are simply always-on and prove nothing.
    """
    import tempfile
    tmp = tempfile.mkdtemp(prefix="anchor_guard_")
    chol = rr.REGISTRY["Cholesterol"][0]
    cases = {
        "seeded_non_independent_import":
            "import reagent_registry\nx = reagent_registry.REGISTRY['DSPC']\n",
        "seeded_copied_template":
            f"S = {chol[:40]!r}\n",
        "control_clean":
            "import math\nS = 'CCCCCCCCCCCC'\nprint(math.pi)\n",
    }
    child = os.path.join(tmp, "guard_child.py")
    with open(child, "w") as f:
        f.write(
            "import json, os, sys\n"
            f"sys.path.insert(0, {HERE!r})\n"
            f"sys.path.insert(0, {os.path.join(REPO, 'provenance')!r})\n"
            "import make_reference as mr\n"
            "p = os.environ['CASE_PATH']\n"
            "prov = mr.ast_provenance(p)\n"
            "tmpl = mr.template_provenance(mr.module_literals(p),\n"
            "                              {'Cholesterol': mr.rr.REGISTRY['Cholesterol'][0]})\n"
            "print('GUARD', json.dumps({'independent': prov['is_independent_of_project_structure_stores'],\n"
            "                           'template_clean': tmpl['clean']}))\n")
    results = {}
    for name, src in cases.items():
        p = os.path.join(tmp, name + ".py")
        with open(p, "w") as f:
            f.write(src)
        env = dict(os.environ, CASE_PATH=p, PYTHONDONTWRITEBYTECODE="1")
        r = subprocess.run([sys.executable, "-I", child], capture_output=True,
                           text=True, env=env, cwd=tmp)
        line = [l for l in r.stdout.splitlines() if l.startswith("GUARD ")]
        if not line:
            results[name] = {"credited": False,
                             "verdict": ("NOT CREDITED: no GUARD line -- the "
                                         "child died before reporting"),
                             "returncode": r.returncode,
                             "stderr_tail": r.stderr.strip()[-300:]}
            continue
        got = json.loads(line[-1][len("GUARD "):])
        if name == "seeded_non_independent_import":
            ok = got["independent"] is False
            why = "AST guard must see the import"
        elif name == "seeded_copied_template":
            ok = got["template_clean"] is False
            why = "template guard must see the copied literal"
        else:
            ok = got["independent"] is True and got["template_clean"] is True
            why = ("control must come back clean, or the guards are always-on "
                   "and certify nothing")
        results[name] = {"credited": ok, "observed": got, "expectation": why,
                         "returncode": r.returncode}
    import shutil
    shutil.rmtree(tmp, ignore_errors=True)
    out = {"n_cases": len(results),
           "n_credited": sum(1 for v in results.values() if v["credited"]),
           "all_credited": all(v["credited"] for v in results.values()),
           "cases": results}
    print(json.dumps(out, indent=1))
    return 0 if out["all_credited"] else 1


def main():
    if "--selftest" in sys.argv:
        return selftest()
    data = rows()
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, "reference_structures.csv"), "w",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADER)
        w.writeheader()
        w.writerows(data)

    targets = {r["lipid_name"]: r["reference_SMILES"] for r in data}
    targets["registry_Cholesterol"] = rr.REGISTRY["Cholesterol"][0]
    targets["registry_DSPC"] = rr.REGISTRY["DSPC"][0]
    qc = {
        "_what": ("C10 anchor posture for PMID 30216729, recorded by hand "
                  "because provenance/structure_check.py and "
                  "check_anchor_independence.py are out of scope for this "
                  "task."),
        "n_anchor_rows": len(data),
        "n_with_a_structure": sum(1 for r in data if r["reference_SMILES"]),
        "n_empty_on_purpose": sum(1 for r in data if not r["reference_SMILES"]),
        "n_INDEPENDENT_anchors": 0,
        "why_zero_independent": ("the only two structures the paper draws that "
                                 "could be written at all are written from "
                                 "provenance/reagent_registry.py, which is the "
                                 "corpus's own source, so each would compare a "
                                 "value to itself"),
        "ast_provenance_of_this_module": ast_provenance(__file__),
        "template_provenance_of_this_module":
            template_provenance(module_literals(__file__), targets),
        "c10_would_report": ("nothing, because no enumerated structure exists "
                             "to compare an anchor against: this paper "
                             "contributes one ionizable component and it has no "
                             "atom-resolved structure"),
    }
    with open(os.path.join(OUTDIR, "04_anchor_qc.json"), "w") as f:
        json.dump(qc, f, indent=1)
    print(json.dumps({k: v for k, v in qc.items() if not k.startswith("_")},
                     indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
