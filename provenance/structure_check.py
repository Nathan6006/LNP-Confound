#!/usr/bin/env python3
"""Tier-1 check: enumerated structures vs the structure the PAPER DRAWS.

Why this exists. Every structural bug found so far survived mass, formula and
descriptor checks, because those validate what a molecule is MADE OF and never
how it is CONNECTED:

  * cholesterol esterified at a side-chain methyl instead of C3 - oxygen count
    unchanged, only one hydrogen moved
  * a radical carbon where a hydroxyl was deleted - parsed and canonicalised fine
  * an aromatic arm on the wrong nitrogen - molecular formula identical

Papers almost always draw their lead compound in full. Comparing one drawn
structure per paper against the enumerated one closes all three classes at once,
because it tests topology directly.

HOW IT WORKS. Each curation folder may carry `reference_structures.csv`:

    lipid_name, reference_SMILES, source_figure, read_on, read_by, note

`reference_SMILES` is read off the figure ONCE by a human or an agent with eyes.
After that this check is fully automatic and costs milliseconds, so the class can
never silently recur - the one-time reading is what buys permanent coverage.

Comparison is by InChIKey, reported in two parts:
  * SKELETON (first block) - connectivity. A mismatch here is a real structural
    error and fails the check.
  * STEREO (second block) - stereochemistry. Reported separately because a
    hand-typed steroid stereo descriptor is easy to get wrong in the REFERENCE,
    and that is not evidence the enumeration is wrong. Flagged, not failed.

An optional `expected_mismatch` column covers the case where the PAPER's drawing
is the thing that is wrong and we deliberately emit the correct structure. Those
rows get their own status instead of CONNECTIVITY_MISMATCH, so a documented
drawing error stays permanently visible in this JSON - rather than living only in
a prose note, which is exactly what ERROR_CLASSES says does not persist - without
being confusable with a fresh bug. Leave the column empty for every normal row.

The reference must be assembled by a DIFFERENT code path than the enumeration
(fragment assembly, or simply a different atom ordering). A reference string
copied from the enumerated one can only ever pass, and tests nothing. This is
NOT checked automatically and cannot be: once both sides are canonicalised, an
honestly re-derived reference is byte-identical to a copied one, so a string
comparison would fail exactly the curations that did it right. The `read_by`
and `note` columns are where the independent code path is recorded, and a
reviewer reads them.

A paper that draws no full structure gets a row with an empty reference_SMILES
and a note saying so, so the ABSENCE of an anchor is visible rather than silent.
"""
import csv
import glob
import json
import os
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
CURATION = os.path.join(HERE, "..", "curation")
OUT = os.path.join(HERE, "structure_check.json")


def key_parts(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, None, None
    k = Chem.MolToInchiKey(m)
    return k, k.split("-")[0], (k.split("-")[1] if "-" in k else "")


def enumerated_lookup(folder):
    """lipid_name -> SMILES from whatever enumeration file the folder has."""
    out = {}
    for name in ("03_enumerated.csv", "06_merged.csv", "13_enumerated_147.csv",
                 "09_enumerated_344.csv"):
        p = os.path.join(folder, name)
        if not os.path.exists(p):
            continue
        for r in csv.DictReader(open(p)):
            if r.get("lipid_name") and r.get("IL_SMILES"):
                out.setdefault(r["lipid_name"], r["IL_SMILES"])
    return out


def main():
    results, n_checked, n_fail, n_stereo, n_missing = [], 0, 0, 0, 0
    n_copy = 0
    n_expected = 0
    for ref_path in sorted(glob.glob(os.path.join(CURATION, "*", "reference_structures.csv"))):
        folder = os.path.dirname(ref_path)
        paper = os.path.basename(folder)
        enum = enumerated_lookup(folder)
        for r in csv.DictReader(open(ref_path)):
            name = r["lipid_name"]
            ref = (r.get("reference_SMILES") or "").strip()
            entry = {"paper": paper, "lipid_name": name,
                     "source_figure": r.get("source_figure", ""),
                     "read_on": r.get("read_on", ""), "note": r.get("note", "")}
            if not ref:
                entry["status"] = "NO_DRAWN_STRUCTURE"
                n_missing += 1
                results.append(entry)
                continue
            got = enum.get(name)
            if not got:
                entry["status"] = "ENUMERATED_STRUCTURE_NOT_FOUND"
                n_fail += 1
                results.append(entry)
                continue
            rk, rskel, rstereo = key_parts(ref)
            gk, gskel, gstereo = key_parts(got)
            entry.update({"reference_inchikey": rk, "enumerated_inchikey": gk,
                          "skeleton_match": rskel == gskel,
                          "stereo_match": rstereo == gstereo})
            n_checked += 1
            expected = (r.get("expected_mismatch") or "").strip()
            if rskel != gskel and expected:
                # The paper DRAWS the molecule wrong and we deliberately emit the
                # correct structure instead. Recorded as its own status so a known
                # drawing error stays visible without being indistinguishable from
                # a new bug; `expected_mismatch` says why, in the reference file.
                entry["status"] = "expected_mismatch_paper_drawing_is_wrong"
                entry["expected_mismatch"] = expected
                n_expected += 1
            elif rskel != gskel:
                entry["status"] = "CONNECTIVITY_MISMATCH"
                n_fail += 1
            elif rstereo != gstereo:
                # A stereo-only difference used to be NON-FAILING on the theory
                # that hand-typed stereo descriptors are the less reliable half.
                # F-013 is what that costs: all 537 ZT_2026 rows carried
                # epicholesterol, the skeleton block of the InChIKey was
                # identical, the formula and all 17 descriptors were identical,
                # and a stereo-only flag would have been the only signal there
                # was. It would have been printed and not acted on.
                #
                # So stereo now FAILS when the reference says it was taken from
                # a pinned authority rather than typed off a figure -- which is
                # exactly when the reference's stereo is the trustworthy half.
                # `stereo_is_authoritative` in the reference file marks that;
                # make_reference.py sets it for anchors built from
                # reagent_registry.json. Everything else stays advisory.
                authoritative = (str(r.get("stereo_is_authoritative", "")).strip()
                                 or ("reagent_registry.json" in entry["note"]))
                if authoritative:
                    entry["status"] = "STEREO_MISMATCH_AGAINST_AUTHORITY"
                    n_fail += 1
                else:
                    entry["status"] = "stereo_differs_connectivity_ok"
                    n_stereo += 1
            else:
                entry["status"] = "exact_match"

            # Can this anchor disagree at all? An anchor byte-identical to the
            # enumeration it checks compares a value to itself and can only
            # pass. That is how F-013 survived: reference_structures.csv's
            # C6O2B2 had the same sha256 as 03_enumerated.csv's, 223 bytes, while
            # its note claimed independent assembly. Identical CANONICAL strings
            # are the correct outcome once both are right, so this cannot be an
            # assertion -- but the stored reference being byte-equal AND already
            # in canonical form is the signature of a copy, and it is worth
            # saying out loud every run.
            if ref == got and ref == Chem.MolToSmiles(Chem.MolFromSmiles(ref)):
                # Suppressed only by something CHECKABLE: a generator script in
                # the folder that produced this row. A prose note claiming
                # independent derivation is not enough -- the note on the
                # anchor that hid F-013 made exactly that claim, and was false.
                gen = os.path.join(folder, "make_reference.py")
                if os.path.exists(gen):
                    entry["anchor_derivation"] = "make_reference.py"
                else:
                    entry["anchor_may_be_a_copy"] = (
                        "stored reference is byte-identical to the enumeration "
                        "and already canonical, and this folder has no "
                        "make_reference.py showing a separate derivation. An "
                        "anchor copied from the enumeration can only pass "
                        "(F-013). Byte-equality is not proof of a copy -- two "
                        "correct routes give the same canonical SMILES -- so "
                        "this is a prompt to confirm, not a failure.")
                    n_copy += 1
            results.append(entry)

    summary = {
        "what": "enumerated structures vs the structure the paper draws",
        "anchors_compared": n_checked,
        "connectivity_failures": n_fail,
        "stereo_only_differences": n_stereo,
        "expected_mismatches_paper_drawing_wrong": n_expected,
        "papers_with_no_drawn_structure": n_missing,
        "anchors_that_may_be_a_copy_of_the_enumeration": n_copy,
        "results": results,
    }
    json.dump(summary, open(OUT, "w"), indent=1)

    print(f"drawn-structure anchors compared: {n_checked}")
    for e in results:
        tag = {"exact_match": "OK  ",
               "stereo_differs_connectivity_ok": "STEREO",
               "CONNECTIVITY_MISMATCH": "FAIL",
               "STEREO_MISMATCH_AGAINST_AUTHORITY": "FAIL",
               "expected_mismatch_paper_drawing_is_wrong": "DRAW",
           "NO_DRAWN_STRUCTURE": "none",
               "ENUMERATED_STRUCTURE_NOT_FOUND": "FAIL"}[e["status"]]
        print(f"  {tag:7} {e['paper']:24} {e['lipid_name']:14} {e['source_figure']:8} {e['status']}")
        if e.get("anchor_may_be_a_copy"):
            print(f"  {'COPY?':7} {e['paper']:24} {e['lipid_name']:14} "
                  f"{e['anchor_may_be_a_copy'][:60]}")
    print(f"\nconnectivity failures: {n_fail}   stereo-only: {n_stereo}   "
          f"papers with no drawn structure: {n_missing}   "
          f"anchors that may be a copy: {n_copy}")
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
