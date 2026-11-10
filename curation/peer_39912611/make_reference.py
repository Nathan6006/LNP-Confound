#!/usr/bin/env python3
"""peer_39912611 -- build the C10 anchors by a route that cannot collide
with stage01_structures.py.

WHY THE ROUTE MATTERS MORE THAN THE STRING. F-013: a C10 anchor whose note
claimed independent fragment assembly was a byte copy of the enumeration, same
sha256, and C10 reported exact_match while comparing a value to itself, on 537
rows of epicholesterol. F-037: the file written to FIX F-013 re-copied the
other half -- the amine template WAS the enumeration's canonical SMILES with
the sterol replaced by '%' and one ring digit renumbered.

And byte identity is NOT the signal either. A correct re-derivation canonicalises
to the same string as a correct enumeration; that is what SHOULD happen. The
signal is whether the two routes share code.

WHAT THIS ROUTE IS. stage01 builds each lipid with RDKit: an atom-mapped acyl
fragment, assemble.esterify / assemble.attach, bonds added to an RWMol. This
file never imports curation/lib, never constructs an RWMol, never calls
esterify or attach, and never touches a canonical SMILES. It emits characters,
walking each lipid's PRINTED IUPAC NAME from the head group inward, one atom
token per locant. The two routes share no function and no constant:
  * stage01's arm comes from ARM_MAPPED, a single literal;
  * this file's arm is assembled here from the two name fragments
    "octane-8,1-diyl" (eight carbons, N at 8, O at 1) and "2-hexyldecanoate"
    (C1 carbonyl, C2 bearing a hexyl, C3..C10 the rest), counted out in the
    loop below rather than written down.
The agreement of two such routes is evidence; the agreement of one route with
itself is not.

RUN:  python3 -I make_reference.py
"""
import csv
import hashlib
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
READ_ON = "2026-11-10"
READ_BY = "curation agent, from the SI at 600 dpi (names: vector text)"

# --- "octane-8,1-diyl 2-hexyldecanoate", counted out from the two names ---
OCTANE_DIYL_N = 8          # octane-8,1-diyl: N on C8, O on C1
HEXYL_N = 6                # 2-HEXYL-
DECANOATE_N = 10           # DECANOATE: C1 is the carbonyl


def arm():
    """Emit the tail arm as characters. No fragment constant is read."""
    chain = "C" * OCTANE_DIYL_N            # C8 (to N) ... C1 (to O)
    acyl_c2_rest = "C" * (DECANOATE_N - 2)  # C3..C10
    hexyl = "C" * HEXYL_N
    return f"{chain}OC(=O)C({hexyl}){acyl_c2_rest}"


# One entry per lipid: the atom tokens of the third N-substituent, read off the
# PRINTED NAME from the core nitrogen outward. The name is quoted so a reviewer
# can check each token against it without opening another file.
HEAD_CHAIN = {
 "Lipid 32": (["C", "C", "C(=O)", "O", "C", "C", "C", "N(C)C"],
              "((3-(3-(Dimethylamino)propoxy)-3-oxopropyl)azanediyl)"
              "bis(octane-8,1-diyl) bis(2-hexyldecanoate)", "SI p.21"),
 "Lipid 33": (["C", "C", "O", "C(=O)", "O", "C", "C", "C", "N(C)C"],
              "11-(8-((2-hexyldecanoyl)oxy)octyl)-2-methyl-7-oxo-6,8-dioxa-"
              "2,11-diazanonadecan-19-yl 2-hexyldecanoate "
              "(walked N11->C10,C9,O8,C7(=O),O6,C5,C4,C3,N2+2-methyl)",
              "SI p.23"),
 "Lipid 34": (["C", "C", "N", "C(=O)", "C", "C", "C", "N(C)C"],
              "((2-(4-(dimethylamino)butanamido)ethyl)azanediyl)"
              "bis(octane-8,1-diyl) bis(2-hexyldecanoate)", "SI p.25"),
 "Lipid 35": (["C", "C", "N", "C(=O)", "N", "C", "C", "C", "N(C)C"],
              "11-(8-((2-hexyldecanoyl)oxy)octyl)-2-methyl-7-oxo-2,6,8,11-"
              "tetraazanonadecan-19-yl 2-hexyldecanoate "
              "(walked N11->C10,C9,N8,C7(=O),N6,C5,C4,C3,N2+2-methyl)",
              "SI p.27"),
 "Lipid 36": (["C", "C", "C(=O)", "N", "C", "C", "C", "N(C)C"],
              "((3-((3-(Dimethylamino)propyl)amino)-3-oxopropyl)azanediyl)"
              "bis(octane-8,1-diyl) bis(2-hexyldecanoate)", "SI p.29"),
 "Lipid 37": (["C", "C", "O", "C(=O)", "C", "C", "N1CCN(C)CC1"],
              "((2-((3-(4-Methylpiperazin-1-yl)propanoyl)oxy)ethyl)azanediyl)"
              "bis(octane-8,1-diyl) bis(2-hexyldecanoate)", "SI p.31"),
 "Lipid 38": (["C", "C", "C(=O)", "O", "C", "C", "N1CCN(C)CC1"],
              "((3-(2-(4-Methylpiperazin-1-yl)ethoxy)-3-oxopropyl)azanediyl)"
              "bis(octane-8,1-diyl) bis(2-hexyldecanoate)", "SI p.33"),
 "Lipid 39": (["C", "C", "O", "C(=O)", "O", "C", "C", "N1CCN(C)CC1"],
              "((2-(((2-(4-methylpiperazin-1-yl)ethoxy)carbonyl)oxy)ethyl)"
              "azanediyl)bis(octane-8,1-diyl) bis(2-hexyldecanoate)",
              "SI p.35"),
 "Lipid 40": (["C", "C", "N", "C(=O)", "C", "C", "N1CCN(C)CC1"],
              "((2-(3-(4-methylpiperazin-1-yl)propanamido)ethyl)azanediyl)"
              "bis(octane-8,1-diyl) bis(2-hexyldecanoate)", "SI p.37"),
 "Lipid 41": (["C", "C", "N", "C(=O)", "N", "C", "C", "N1CCN(C)CC1"],
              "((2-(3-(2-(4-methylpiperazin-1-yl)ethyl)ureido)ethyl)azanediyl)"
              "bis(octane-8,1-diyl) bis(2-hexyldecanoate)", "SI p.39"),
 "Lipid 42": (["C", "C", "C(=O)", "N", "C", "C", "N1CCN(C)CC1"],
              "((3-((2-(4-Methylpiperazin-1-yl)ethyl)amino)-3-oxopropyl)"
              "azanediyl)bis(octane-8,1-diyl) bis(2-hexyldecanoate)",
              "SI p.41"),
 # Lipid 14 has NO write-up and NO name in this SI -- the authors synthesised
 # it "according to our previously reported procedures" (their ref 41). Its
 # anchor is read off the DRAWING instead: Supp Fig 15 row 1 column 1, at
 # 600 dpi, three CH2 between the ester carbonyl and the dimethylamino N.
 "Lipid 14": (["C", "C", "O", "C(=O)", "C", "C", "C", "N(C)C"],
              "NO IUPAC NAME IN THIS SI. Read from the drawn product in "
              "Supplementary Figure 15, row 1 column 1, and from its drawn "
              "reagent 4-(dimethylamino)butyric acid.", "SI p.17 Supp Fig 15"),
}


def build(tokens):
    a = arm()
    return f"N({a})({a})" + "".join(tokens)


def main():
    rows = []
    for name in ("Lipid 14", "Lipid 32", "Lipid 33", "Lipid 34", "Lipid 35",
                 "Lipid 36", "Lipid 37", "Lipid 38", "Lipid 39", "Lipid 40",
                 "Lipid 41", "Lipid 42"):
        tokens, src, page = HEAD_CHAIN[name]
        smi = build(tokens)
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise SystemExit(f"{name}: anchor does not parse: {smi}")
        # Ring-closure sanity: the only ring any of these may contain is the
        # N-methylpiperazine. String composition once fused a piperazine to a
        # cholesterol through a ring-digit collision; this is the guard.
        nring = m.GetRingInfo().NumRings()
        want_ring = 1 if any("1" in t for t in tokens) else 0
        if nring != want_ring:
            raise SystemExit(f"{name}: {nring} rings, expected {want_ring} -- "
                             "ring-closure digit collision in the composed "
                             "string.")
        if any(at.GetNumRadicalElectrons() for at in m.GetAtoms()):
            raise SystemExit(f"{name}: radical electron in the anchor.")
        rows.append({
            "lipid_name": name,
            "reference_SMILES": Chem.MolToSmiles(m),
            "source_figure": page,
            "read_on": READ_ON,
            "read_by": READ_BY,
            "expected_mismatch": "",
            "note": (f"Anchor composed character-by-character from: {src} "
                     "No curation/lib import, no RWMol, no esterify/attach, "
                     "no rewriting of a canonical SMILES. Stereochemistry is "
                     "deliberately ABSENT: the two 2-hexyldecanoate C2 centres "
                     "are drawn flat everywhere in the paper and the reagent "
                     "is a racemate, so an anchor carrying a configuration "
                     "would assert something the source does not. "
                     f"formula {CalcMolFormula(m)}."),
        })
    # An explicitly EMPTY row, so the one real gap is visible rather than
    # inferred from its absence.
    rows.append({
        "lipid_name": "Lipid 14 (as published in ref 41)",
        "reference_SMILES": "",
        "source_figure": "",
        "read_on": READ_ON, "read_by": READ_BY, "expected_mismatch": "",
        "note": "NO ANCHOR FROM THIS PAPER'S OWN CHARACTERISATION. Lipid 14 "
                "has no synthesis write-up, no IUPAC name, no NMR and no "
                "ESI-MS in this SI; it is cited to the authors' earlier work. "
                "The Lipid 14 row above is anchored on the SCHEME DRAWING "
                "only. An independent corroboration exists outside this "
                "paper and is recorded in 08_contradictions.csv row 1: "
                "ZT_2026's Lipid 14-CB1/2/3 in data/data_vivo.csv carry this "
                "exact scaffold with the same C4 acid.",
    })

    out = os.path.join(HERE, "reference_structures.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        for r in rows:
            w.writerow([r["lipid_name"], r["reference_SMILES"],
                        r["source_figure"], r["read_on"], r["read_by"],
                        r["expected_mismatch"], r["note"]])

    # Report the comparison here too, since provenance/structure_check.py may
    # NOT be run in this session (it writes tracked files other agents use).
    enum = {}
    with open(os.path.join(HERE, "03_enumerated.csv")) as fh:
        for r in csv.DictReader(fh):
            enum[r["lipid_name"]] = r["SMILES"]
    nmatch = 0
    print(f"{'lipid':10s} {'anchor==enum':13s} {'same sha256':12s} inchikey")
    for r in rows[:-1]:
        a = r["reference_SMILES"]
        b = enum[r["lipid_name"]]
        same = Chem.MolToInchiKey(Chem.MolFromSmiles(a)) == \
            Chem.MolToInchiKey(Chem.MolFromSmiles(b))
        nmatch += same
        print(f"{r['lipid_name']:10s} {str(same):13s} "
              f"{str(hashlib.sha256(a.encode()).hexdigest() == hashlib.sha256(b.encode()).hexdigest()):12s} "
              f"{Chem.MolToInchiKey(Chem.MolFromSmiles(a))}")
    print(f"C10 (by hand): {nmatch}/12 anchors InChIKey-identical to the "
          "enumeration.")
    return 0 if nmatch == 12 else 1


if __name__ == "__main__":
    sys.exit(main())
