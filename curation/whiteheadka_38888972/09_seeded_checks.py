#!/usr/bin/env python3
"""Seed each of this folder's own assertions and confirm it FIRES.

CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON
PURPOSE." That rule cost three incidents -- a re-copied anchor that could only
pass, three ZT_2026 marginal checks that scored a literal 0 for twelve days
while the summary beside them said ALL FOUR EXACT, and two versions of C16 that
reported "all derived" while examining nothing. Every one of them printed a
clean line.

curation/lib's own invariants are covered by provenance/test_lib.py. What is
NOT covered there is the handful of assertions written in THIS folder, and the
one place where a lib refusal was handled rather than propagated (the prochiral
marker in stage03_enumerate.py). Those are seeded here.

Each case states the defect, runs it, and requires the named failure. A case
that does NOT raise is itself a failure of this script.

Run: python3 -I 09_seeded_checks.py      (non-zero exit if a check is toothless)
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)

from rdkit import Chem                                          # noqa: E402

import assemble                                                 # noqa: E402
import emit                                                     # noqa: E402
from errors import EmitError, StructureError                    # noqa: E402

import importlib.util                                           # noqa: E402

spec = importlib.util.spec_from_file_location(
    "s3", os.path.join(HERE, "stage03_enumerate.py"))
s3 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(s3)

RESULTS = []


def case(name, defect, want, fn):
    try:
        fn()
    except BaseException as e:                       # noqa: BLE001
        ok = want in type(e).__name__ or want in str(e)
        RESULTS.append({"case": name, "defect": defect, "fired": True,
                        "expected": want, "matched": bool(ok),
                        "message": str(e).split("\n")[0][:220]})
        return
    RESULTS.append({"case": name, "defect": defect, "fired": False,
                    "expected": want, "matched": False,
                    "message": "DID NOT RAISE -- this check has no teeth."})


def main():
    blocks = {r["code"]: r for r in csv.DictReader(
        open(os.path.join(HERE, "01_blocks.csv")))}
    head306 = blocks["306"]["SMILES"]
    adduct = blocks["O10"]["adduct_fragment_SMILES"]

    # 1. The final stereo assertion must still catch the failure mode it
    #    exists for. NOTE: curation/lib/assemble.py was relaxed mid-session
    #    (2026-11-10 00:35) so that CREATING an unassignable centre no longer
    #    refuses -- only LOSING or INVERTING an assigned one does. An earlier
    #    version of this case seeded a created centre and reported
    #    "DID NOT RAISE", which is how the library change was noticed at all.
    #    So the case now seeds a LOST assigned centre, which is the F-013
    #    shape and is what the assertion still guards.
    def lose_a_centre():
        smi = s3.build(head306, adduct, 4, "306O10")
        assemble.check_stereo_preserved(
            Chem.MolFromSmiles(smi),
            [head306] + [adduct] * 3 + ["C[C@H](N)CC"],
            label="306O10 seeded")
    case("assigned stereocentre lost",
         "the fragment list carries an assigned (S) centre that the product "
         "does not -- the F-013 shape",
         "StructureError", lose_a_centre)

    # 1b. POSITIVE CONTROL for the prochiral marker: head 401's finished
    #     lipidoid must be stereochemically identical to head + 4 tails, i.e.
    #     the 13C marker really is gone and the transient centre really has
    #     resolved. A seeded check that can only fire is half a check.
    def pivot_resolved():
        smi = s3.build(blocks["401"]["SMILES"], adduct, 4, "401O10")
        if "13C" in smi or "[13" in smi:
            raise AssertionError("401O10 still carries the 13C marker")
        fp = assemble.stereo_fingerprint(Chem.MolFromSmiles(smi))
        if fp != {"R": 0, "S": 0, "?": 0}:
            raise AssertionError(f"401O10 stereo fingerprint is {fp}, want 0")
        RESULTS.append({"case": "401 pivot resolves (positive control)",
                        "defect": "none -- this one must NOT raise",
                        "fired": False, "expected": "no exception",
                        "matched": True,
                        "message": "401O10 is marker-free and achiral, as the "
                                   "symmetric product must be"})
    try:
        pivot_resolved()
    except BaseException as e:                       # noqa: BLE001
        RESULTS.append({"case": "401 pivot resolves (positive control)",
                        "defect": "none -- this one must NOT raise",
                        "fired": True, "expected": "no exception",
                        "matched": False, "message": str(e)[:200]})

    # 2. The substitution count IS the structural claim. A head whose N-H count
    #    is misstated must not silently produce a differently-substituted
    #    molecule.
    case("wrong tail count",
         "head 306 (4 N-H) declared as 3 N-H",
         "added 4 tails but the head has 3",
         lambda: s3.build(head306, adduct, 3, "306O10 seeded"))

    # 3. A head left with an N-H must not pass as fully substituted.
    def undersubstituted():
        smi = s3.build(head306, adduct, 4, "306O10")
        m = Chem.MolFromSmiles(smi)
        n = len(m.GetSubstructMatches(
            Chem.MolFromSmarts("[CX3](=O)[OX2][CX4]")))
        if n != 3:
            raise AssertionError(
                f"ester-count assertion: {n} esters for a claimed 3 tails")
    case("ester count vs tail count",
         "the product's ester count is compared against a tail count of 3",
         "ester-count assertion", undersubstituted)

    # 4. emit.emit_rows must refuse a literal confidence. This is the shape
    #    F-011, F-025 and F-046 all took.
    case("literal confidence",
         'confidence="high" passed as a string instead of per-item',
         "EmitError",
         lambda: emit.emit_rows(["x"], {"x": {}}, confidence="high",
                                note=lambda i, r: "n", label="seeded"))

    # 5. A caveat with no note must be refused.
    case("caveat without a note",
         "structure_assumed=1 with an empty curation_note",
         "EmitError",
         lambda: emit.emit_rows(["x"], {"x": {}},
                                confidence=lambda i, r: "low",
                                note=lambda i, r: "",
                                structure_assumed=lambda i, r: "1",
                                label="seeded"))

    # 6. The Figure 2B label/tick count guard. A dropped label would shift
    #    every value one molecule to the left -- the AC_2025 shape.
    def drop_a_label():
        import stage04_fig2b as s4
        keep = list(s4.LABELS)
        s4.LABELS = keep[:-1]
        try:
            s4.main()
        finally:
            s4.LABELS = keep
    case("Figure 2B label dropped",
         "one x-axis label removed from the transcription",
         "x ticks but", drop_a_label)

    # 7. The anchor's own assertions. Removing one tail from the template must
    #    not produce a quietly-wrong anchor.
    def anchor_three_tails():
        import make_reference as mr
        t = mr.TEMPLATE
        mr.TEMPLATE = "{A}N({B})CCCN(C)CCCN{B}"   # one tail short
        try:
            mr.main()
        finally:
            mr.TEMPLATE = t
    case("anchor with 3 tails",
         "the C10 anchor template is built with 3 tails instead of the 4 "
         "Figure 1A draws",
         "esters, Figure 1A draws 4", anchor_three_tails)

    out = {"cases": RESULTS,
           "n": len(RESULTS),
           "n_fired": sum(1 for r in RESULTS if r["fired"]),
           "n_matched": sum(1 for r in RESULTS if r["matched"])}
    with open(os.path.join(HERE, "09_seeded_checks.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for r in RESULTS:
        print(("PASS " if r["matched"] else "FAIL ") + r["case"]
              + " -> " + r["message"][:150])
    print(f"{out['n_matched']}/{out['n']} seeded defects were caught")
    # Re-run the real stages so the folder is left with correct artifacts.
    import stage04_fig2b as s4
    s4.main()
    import make_reference as mr
    mr.main()
    return 0 if out["n_matched"] == out["n"] else 1


if __name__ == "__main__":
    sys.exit(main())
