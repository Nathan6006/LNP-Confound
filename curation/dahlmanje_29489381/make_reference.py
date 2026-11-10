#!/usr/bin/env python3
"""The C10 anchor for curation/dahlmanje_29489381.

THE PAPER DRAWS NO COMPLETE IONIZABLE LIPID.  Every lipomer scheme in the SI
(Fig S2a p7, Fig S3d p12, Fig S3f p13, Fig S5a p23) is a generic Markush: an
epoxide or (meth)acrylate carrying an unlabelled 'R', plus an amine that is
either drawn in full (the oligoamines) or abstracted to diethylamine with
'PEI_600' written beneath it, going to a product drawn with one bond formed
and the rest of the molecule left as 'R'.  So `reference_SMILES` is EMPTY on
all fifteen lipomer rows and the note says why.  PROTOCOL: store the empty row
so the ABSENCE of an anchor is visible rather than silent.

What CAN be anchored is the block set, and this file does it by a route that
shares nothing with the enumerator.  stage05_blocks.py holds a typed SMILES
literal per reagent.  Here each structure is GENERATED from the numbers in the
reagent's own chemical NAME -- segment lengths for the polyamines, ring size
and nitrogen positions for the macrocycle, chain length minus two for the
epoxides -- and the two are compared by full InChIKey.  No SMILES literal is
copied across, which is the thing F-013 and F-037 both turned out to have
done while claiming otherwise.
"""
import csv
import json
import os
import sys

from rdkit import Chem

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def polyamine(segments):
    """H2N-(CH2)a-NH-(CH2)b-NH-...-NH2 from the segment lengths the NAME gives.

    spermidine  = N-(3-AMINOPROPYL)BUTANE-1,4-DIAMINE  -> [3, 4]
    spermine    = N,N'-BIS(3-AMINOPROPYL)BUTANE-1,4-DIAMINE -> [3, 4, 3]
    triethylenetetramine = N,N'-BIS(2-AMINOETHYL)ETHANE-1,2-DIAMINE -> [2,2,2]
    """
    parts = ["N"]
    for n in segments:
        parts.append("C" * n)
        parts.append("N")
    return "".join(parts)


def azamacrocycle(ring_size, n_positions):
    """A saturated azamacrocycle from its ring size and nitrogen LOCANTS.

    1,4,8,12-tetraazacyclopentadecane -> ring_size 15, locants (1,4,8,12).
    The carbon runs between nitrogens are derived from the locants, never
    typed: gaps = (next locant - this locant - 1), wrapping at the ring.
    """
    locs = list(n_positions)
    gaps = []
    for i, a in enumerate(locs):
        b = locs[(i + 1) % len(locs)]
        gaps.append((b - a - 1) if i + 1 < len(locs)
                    else (ring_size + locs[0] - a - 1))
    out = "N1"
    for g in gaps[:-1]:
        out += "C" * g + "N"
    out += "C" * gaps[-1] + "1"
    return out


def epoxyalkane(total_carbons):
    """1,2-epoxy-n-alkane.  The oxirane consumes 2 of the n carbons, so the
    alkyl substituent is C(n-2) -- which is what the SI's own colour key
    confirms ('C12 Epoxy' <-> 'R1=C10', 'C15 Epoxy' <-> 'R2=C13')."""
    return "C" * (total_carbons - 2) + "C1CO1"


def alkyl_ester(acid_smiles, alkyl_carbons):
    return acid_smiles + "O" + "C" * alkyl_carbons


GENERATED = {
    "Spermidine": polyamine([3, 4]),
    "Spermine Free Base": polyamine([3, 4, 3]),
    "Triethylenetetramine": polyamine([2, 2, 2]),
    "1,4,8,12-Tetraazacyclopentadecane": azamacrocycle(15, (1, 4, 8, 12)),
    "C12 Epoxy": epoxyalkane(12),
    "C15 Epoxy": epoxyalkane(15),
    "Lauryl Acrylate": alkyl_ester("C=CC(=O)", 12),
    "Stearyl Methacrylate": alkyl_ester("CC(=C)C(=O)", 18),
}

HEADER = ["lipid_name", "reference_SMILES", "source_figure", "read_on",
          "read_by", "expected_mismatch", "note"]

LIPOMERS = ["7C1", "8", "48", "58", "68", "69", "78", "79", "88", "92",
            "100", "102", "104", "106", "108"]


def main():
    import stage05_blocks as B

    results, mismatches = [], []
    for name, gen in sorted(GENERATED.items()):
        typed = (B.AMINES.get(name) or B.TAILS.get(name))[0]
        mg, mt = Chem.MolFromSmiles(gen), Chem.MolFromSmiles(typed)
        kg, kt = Chem.MolToInchiKey(mg), Chem.MolToInchiKey(mt)
        same = kg == kt
        results.append({
            "block": name,
            "generated_from_the_name": gen,
            "generated_inchikey": kg,
            "enumerator_inchikey": kt,
            "match": same,
            "byte_identical_to_enumerator_string": gen == typed,
        })
        if not same:
            mismatches.append(name)

    rows = []
    for lid in LIPOMERS:
        am, tl = B.LIPOMERS[lid]
        rows.append({
            "lipid_name": "lipomer %s (%s + %s)" % (lid, am, tl),
            "reference_SMILES": "",
            "source_figure": ("SI Fig S2a p7 / Fig S3d p12 / Fig S3f p13 / "
                              "Fig S5a p23 (generic scheme only)"),
            "read_on": "2026-11-10",
            "read_by": "curation agent, dahlmanje_29489381",
            "expected_mismatch": "",
            "note": ("NO ANCHOR POSSIBLE: the paper draws no complete lipomer "
                     "anywhere. The scheme is a Markush with an unlabelled R "
                     "on the tail and, for the PEI lipomers, the amine itself "
                     "abstracted to diethylamine. The degree of substitution "
                     "is never stated. This curation emits 0 corpus rows, so "
                     "the anchor guards nothing yet; the row exists so the "
                     "absence is visible."),
        })

    with open(os.path.join(HERE, "reference_structures.csv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=HEADER)
        w.writeheader()
        for r in rows:
            w.writerow(r)

    rep = {
        "n_anchor_rows": len(rows),
        "n_with_a_reference_SMILES": 0,
        "why_empty": ("the paper draws only a generic Markush for every "
                      "lipomer; see this file's docstring"),
        "block_independent_rederivation": {
            "n_blocks_checked": len(results),
            "n_inchikey_matches": sum(1 for r in results if r["match"]),
            "n_mismatches": len(mismatches),
            "mismatched": mismatches,
            "n_byte_identical_strings": sum(
                1 for r in results if r["byte_identical_to_enumerator_string"]),
            "detail": results,
        },
    }
    with open(os.path.join(HERE, "06_reference_qc.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    for r in results:
        print("%-36s %s  gen=%s" % (r["block"],
                                    "MATCH " if r["match"] else "MISMATCH",
                                    r["generated_from_the_name"]))
    print("\n%d/%d blocks re-derived from the name agree by full InChIKey; "
          "%d generated strings are byte-identical to the enumerator's "
          "(after canonicalisation they must be -- that is not the signal)."
          % (rep["block_independent_rederivation"]["n_inchikey_matches"],
             len(results),
             rep["block_independent_rederivation"]["n_byte_identical_strings"]))
    return 1 if mismatches else 0


if __name__ == "__main__":
    sys.exit(main())
