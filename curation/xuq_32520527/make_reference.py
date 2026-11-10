#!/usr/bin/env python3
"""XQ_2020 -- BUILD the C10 reference anchors. Run, do not hand-edit the CSV.

    python3 -I make_reference.py

WHAT AN ANCHOR IS FOR. C10 compares each enumerated structure against the
structure the PAPER DRAWS. It only tests anything if the anchor was produced by
a code path that shares as little as possible with the enumerator. That claim
has been made twice in this project and been false twice: F-013 (a byte copy
whose note claimed independent fragment assembly) and F-037 (inside the file
written to fix F-013, where the amine half was still the enumeration's
canonical SMILES with the sterol swapped for '%').

WHAT THIS PAPER ACTUALLY DRAWS, stated plainly so the anchor is not oversold:

  Figure 1C draws all ELEVEN amine heads in full.
  Figure 1B draws the O16CBA tail in full.
  Figure 1B draws the LIPIDOID only as 'R-N[arm]2' -- R is a letter and the
  arm count is a bracket subscript.

So the paper draws no complete lipidoid anywhere. Every anchor below is a
COMPOSITION of two separately drawn fragments, and that is what it tests: it
tests that the composition is right, not that a drawn whole was copied
correctly. Saying so is the point of the `note` column.

HOW THIS IS INDEPENDENT OF stage03_enumerate.py, concretely:

  stage03   builds a molecule GRAPH. It derives the arm by reducing O16CBA's
            alkene under a SMARTS match, then joins fragments with
            assemble.attach, which is Chem.RWMol + AddBond on atom-mapped
            attachment atoms, once per arm, breadth-first over every N-H site.
  this file builds a STRING. Plain str.format substitution of a hand-written
            arm into a hand-written head template, with no RDKit call of any
            kind until the finished string is parsed for comparison. There is
            no RWMol, no AddBond, no atom map, no SMARTS and no import of
            curation/lib/assemble anywhere in this module.

  The heads are re-written here in a DIFFERENT TRAVERSAL from 01_blocks.csv --
  amine-first instead of substituent-first -- and 01_blocks.csv is never read.
  The aryl ring and the dioxane ring are likewise entered at different atoms
  than stage03 enters them.

WHY BYTE-DIFFERENCE IS CHECKED ON THE TEMPLATE AND NOT ON THE RESULT. Once
canonicalised, a CORRECT re-derivation is byte-identical to a CORRECT
enumeration. Testing the emitted anchor for byte-difference would therefore
flag correct work and miss the actual defect, which is exactly why C10's own
`anchor_may_be_a_copy` missed F-037. The test that means something is on the
INPUT: the template this file composes must not be the enumeration's canonical
SMILES, and must not share a long literal run with it.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

READ_ON = "2026-11-10"
READ_BY = "curation agent, from source/nihms-1799810-f0001.jpg at 3-4x"

# ---------------------------------------------------------------------------
# The arm, hand-written from Figure 1B, read ester-first.
# Figure 1B draws:  -N-CH2CH2-C(=O)-O-CH2-C(CH3)< ring >
# ring = 1,3-dioxane whose C2 carries 4-hexyloxy-2,6-dimethoxyphenyl.
# {a}/{b} are ring-closure digits, supplied per arm so several arms can coexist
# in one string.
# ---------------------------------------------------------------------------
ARM = ("CCC(=O)OCC{a}(C)COC(c{b}c(OC)cc(OCCCCCC)cc{b}OC)OC{a}")

# ---------------------------------------------------------------------------
# The eleven heads, hand-written from Figure 1C, entered at the REACTING
# nitrogen and written outward -- the opposite of 01_blocks.csv, which is
# entered at the far substituent. {A} marks an arm.
#
# The number of {A} slots is the arm count Figure 2F's molecular formula forces
# (see 03_enumerated_qc.json). For 113 and 306 the 2+1 split across the two
# primary nitrogens is the only possibility; for 400 the three N-H slots are
# all filled.
# ---------------------------------------------------------------------------
HEADS = {
    "75":  "N({A})({A})CCN(CC)CC",
    "76":  "N({A})({A})CCN1CCCC1",
    "77":  "N({A})({A})CCC1CCCN1C",
    "78":  "N({A})({A})CCN(C)C",
    "80":  "N({A})({A})CCCN(C)C",
    "81":  "N({A})({A})CCCN(CC)CC",
    "82":  "N({A})({A})CCN1CCCCC1",
    "93":  "N({A})({A})CCCn1ccnc1",
    "113": "N({A})({A})CCN(C)CCN{A}",
    "306": "N({A})({A})CCCN(C)CCCN{A}",
    "400": "N({A})({A})CCCN({A})CCCN(C)C",
}

NOTE = (
    "COMPOSED ANCHOR, not a drawn whole: the paper draws the amine head in "
    "full (Figure 1C) and the O16CBA tail in full (Figure 1B) but draws the "
    "lipidoid only as R-N[tail]2, so no complete structure is drawn anywhere. "
    "Built here by plain string substitution of a hand-written arm into a "
    "hand-written head template, in a different traversal from "
    "stage03_enumerate.py and with no RDKit assembly call; the arm count is "
    "the one Figure 2F's printed molecular formula forces. Stereochemistry is "
    "absent on purpose: the 1,3-dioxane's C2/C5 cis-trans relationship is "
    "never drawn (DECISIONS.md D-2)."
)
NOTE_77 = (
    "As the other ten, and additionally: amine 77's pyrrolidine C2 is a "
    "genuine stereocentre drawn flat in Figure 1C; the anchor is flat too, so "
    "it cannot test a configuration the paper never states."
)


def compose(head_id):
    ring = iter(range(1, 10))
    # ring digit 1 is already used by the pyrrolidine/piperidine/imidazole in
    # five of the heads, so arms start at 2.
    digits = [d for d in range(2, 10)]
    out = HEADS[head_id]
    i = 0
    while "{A}" in out:
        a, b = digits[i], digits[i + 1]
        out = out.replace("{A}", ARM.format(a=a, b=b), 1)
        i += 2
    return out


def main():
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")

    enum = {r["head_amine"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}

    rows, report = [], []
    for hid in HEADS:
        template = compose(hid)
        m = Chem.MolFromSmiles(template)
        if m is None:
            raise SystemExit(f"make_reference: {hid} template does not parse:\n"
                             f"  {template}")
        anchor = Chem.MolToSmiles(m)
        theirs = Chem.CanonSmiles(enum[hid]["IL_SMILES"])

        # --- the independence test, on the TEMPLATE (see module docstring)
        longest = 0
        for L in range(len(template), 0, -1):
            if any(template[s:s + L] in theirs
                   for s in range(0, len(template) - L + 1)):
                longest = L
                break
        report.append({
            "lipidoid": f"{hid}-O16CBA",
            "template_is_byte_identical_to_the_enumeration":
                template == theirs,
            "longest_shared_literal_run_chars": longest,
            "template_length_chars": len(template),
            "inchikey_anchor": Chem.MolToInchiKey(m),
            "inchikey_enumeration": Chem.MolToInchiKey(
                Chem.MolFromSmiles(theirs)),
            "exact_match": Chem.MolToInchiKey(m) == Chem.MolToInchiKey(
                Chem.MolFromSmiles(theirs)),
            "skeleton_match": Chem.MolToInchiKey(m).split("-")[0]
                == Chem.MolToInchiKey(Chem.MolFromSmiles(theirs)).split("-")[0],
        })
        rows.append({
            "lipid_name": f"{hid}-O16CBA",
            "reference_SMILES": anchor,
            "source_figure": "Figure 1B (O16CBA tail, drawn in full) + "
                             "Figure 1C (amine head, drawn in full); arm count "
                             "from Figure 2F's printed molecular formula",
            "read_on": READ_ON,
            "read_by": READ_BY,
            "expected_mismatch": "",
            "note": NOTE_77 if hid == "77" else NOTE,
        })

    copies = [r for r in report
              if r["template_is_byte_identical_to_the_enumeration"]]
    if copies:
        raise SystemExit(
            "make_reference: the composed TEMPLATE is byte-identical to the "
            f"enumeration for {[r['lipidoid'] for r in copies]}. That is "
            "F-013/F-037 exactly: an anchor built from the enumerator's own "
            "output compares a value to itself. Re-derive the template from "
            "the figure.")
    bad = [r for r in report if not r["exact_match"]]
    if bad:
        raise SystemExit(
            "make_reference: anchor and enumeration DISAGREE for "
            f"{[r['lipidoid'] for r in bad]}. One of the two is wrong; this "
            "is what the anchor exists to say. Do not 'fix' the anchor to "
            "match -- find out which half is wrong.")

    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=[
            "lipid_name", "reference_SMILES", "source_figure", "read_on",
            "read_by", "expected_mismatch", "note"])
        w.writeheader()
        w.writerows(rows)
    with open(os.path.join(HERE, "anchor_independence.json"), "w") as f:
        json.dump({
            "_what": "C16-style evidence, computed here because this run may "
                     "not execute provenance/check_anchor_independence.py "
                     "(it writes tracked files other agents are using). The "
                     "test is on the GENERATOR'S COMPOSED TEMPLATE, not on "
                     "the emitted anchor: a correct re-derivation is "
                     "byte-identical to a correct enumeration once "
                     "canonicalised, which is why C10's own "
                     "`anchor_may_be_a_copy` missed F-037.",
            "n_anchors": len(report),
            "n_byte_identical_templates": len(copies),
            "worst_shared_literal_run_chars": max(
                r["longest_shared_literal_run_chars"] for r in report),
            "n_exact_inchikey_matches": sum(1 for r in report
                                            if r["exact_match"]),
            "anchors": report,
        }, f, indent=1)
    print(f"make_reference: {len(rows)} anchors, "
          f"{sum(1 for r in report if r['exact_match'])}/{len(report)} exact "
          f"InChIKey match, worst shared literal run "
          f"{max(r['longest_shared_literal_run_chars'] for r in report)} chars "
          f"of {min(r['template_length_chars'] for r in report)}+")
    return 0


if __name__ == "__main__":
    sys.exit(main())
