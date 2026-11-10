#!/usr/bin/env python3
"""The monomer pool. Writes 01_blocks.csv, and MEASURES two chain lengths.

Run from this folder:   python3 -I stage01_blocks.py

blocks.py holds the side-chain structures read off main-article Figure 1B. This
stage does the two things that make those readings checkable:

1. VECTOR CHAIN LENGTHS. SI page 20 (Figure S7) redraws Dod and Ehx in its own
   legend as vector ChemDraw output, so the PDF carries their bond coordinates.
   curation/lib chains reads those coordinates. This is independent of the
   three-letter code and independent of Table S4's arithmetic -- it is the only
   route of the three that touches geometry at all.

   Carbons = merged vertices - 1, because the chain's first vertex is the
   attachment point (where the squiggle meets the backbone nitrogen), not a
   carbon. That convention is DERIVED, not assumed: it is the only convention
   that gives the right answer on BOTH chains at once, and the two chains have
   different lengths, so a wrong convention could not satisfy both.

2. CARBON ARITHMETIC. SI Table S4 prints "Total Lipid Carbons" for 34 peptoids.
   Summing blocks.py's per-monomer carbon counts over each sequence must
   reproduce it. With 34 equations over 6 unknowns the system is massively
   over-determined, so this pins the whole lipid monomer set at once -- it is
   the check that actually decides that Hex/Oct/Dec/Dod are 6/8/10/12 rather
   than the name token saying so. CLAUDE.md section 5: name tokens are not
   chain lengths.

   curation/lib chains.check_name_token REPORTS the token-vs-measured
   comparison for each. It never sets a value.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)

import chains                      # noqa: E402
import pymupdf                     # noqa: E402
from rdkit import Chem, RDLogger   # noqa: E402

import blocks                      # noqa: E402
import seqtables                   # noqa: E402

RDLogger.DisableLog("rdApp.*")

SI = os.path.join(HERE, "source", "nn4c05513_si_001.pdf")

# Bounding boxes, in PDF points, of the two legend chains on SI page 20. Taken
# from the page's own word positions: the "=" sign at x~107 and the labels
# "Dod"/"(C12)" at x~192 and "Ehx"/"(C8)" at x~154 bracket each chain.
LEGEND = {
    "Dod": dict(page=20, bbox=(110, 70, 191, 90)),
    "Ehx": dict(page=20, bbox=(110, 88, 154, 115)),
}

# How many carbons the three-letter code ITSELF suggests, for the report-only
# comparison. None where the code implies nothing numeric.
NAME_TOKEN = {"Hex": 6, "Oct": 8, "Dec": 10, "Dod": 12, "Ole": None,
              "Ehx": None}


def vector_lengths():
    doc = pymupdf.open(SI)
    out = {}
    for code, spec in LEGEND.items():
        pg = doc[spec["page"] - 1]
        bonds = chains.vector_bonds(pg, spec["bbox"])
        v = chains.vector_chain_length(bonds)
        out[code] = {"stroked_segments": len(bonds), "merged_vertices": v,
                     "carbons": v - 1}
    doc.close()
    return out


def carbon_arithmetic():
    s4, _ = seqtables.table_s4()
    rows, agree = [], 0
    for pid, rec in sorted(s4.items()):
        lipids = [m for m in rec["monomers"]
                  if blocks.MONOMERS[m]["kind"] == "lipid"]
        n_c = sum(blocks.MONOMERS[m]["n_carbons"] for m in lipids)
        n_bu = sum(1 for m in lipids
                   if blocks.MONOMERS[m]["branched"]
                   or blocks.MONOMERS[m]["unsaturated"])
        ok = (n_c == rec["lipid_carbons"] and len(lipids) == rec["lipid_monomers"]
              and n_bu == rec["branched_unsaturated"])
        agree += ok
        rows.append({"peptoid": pid, "ours_lipid_monomers": len(lipids),
                     "paper_lipid_monomers": rec["lipid_monomers"],
                     "ours_lipid_carbons": n_c,
                     "paper_lipid_carbons": rec["lipid_carbons"],
                     "ours_branched_unsat": n_bu,
                     "paper_branched_unsat": rec["branched_unsaturated"],
                     "agree": ok})
    return rows, agree


def main():
    vec = vector_lengths()
    arith, agree = carbon_arithmetic()
    s9, _ = seqtables.table_s9()
    used = sorted({m for r in s9.values() for m in r["monomers"]})

    name_token_reports = []
    for code, n in NAME_TOKEN.items():
        if n is None:
            continue
        name_token_reports.append(
            chains.check_name_token(n, blocks.MONOMERS[code]["n_carbons"],
                                    label=code))

    out = os.path.join(HERE, "01_blocks.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["code", "kind", "sidechain_SMILES", "n_carbons",
                    "amine_order", "branched", "unsaturated",
                    "used_in_222_266", "resolved", "drawn_in_Figure_1B",
                    "vector_measured_carbons", "note"])
        for code in sorted(blocks.MONOMERS):
            m = blocks.MONOMERS[code]
            side = m["sidechain"]
            if side is None and m["amine_order"] == "secondary":
                side = f"N({m['arms'][0]})({m['arms'][1]})-  [on the backbone N]"
            w.writerow([
                code, m["kind"], side or "", m["n_carbons"], m["amine_order"],
                int(m["branched"]), int(m["unsaturated"]),
                "yes" if code in used else "no",
                "yes" if (m["sidechain"] or m.get("arms")) else "no",
                "yes" if m["drawn_as"] else "no",
                vec[code]["carbons"] if code in vec else "",
                m["note"]])
        for code in sorted(set(blocks.NOT_ENUMERATED_POOL)):
            w.writerow([code, "", "", "", "", "", "", "no", "no", "yes", "",
                        "appears only in SI Table S1's 221-member screening "
                        "library, which has no extractable value; not "
                        "transcribed. See blocks.py."])

    report = {"vector_chain_lengths": vec,
              "name_token_reports": name_token_reports,
              "carbon_arithmetic_agree": f"{agree}/{len(arith)}",
              "carbon_arithmetic_rows": arith,
              "monomers_used_in_222_266": used}
    with open(os.path.join(HERE, "01_blocks_qc.json"), "w") as fh:
        json.dump(report, fh, indent=1)

    print("vector chain lengths, from the PDF's own bond coordinates (SI p20):")
    for code, v in vec.items():
        print(f"  {code}: {v['stroked_segments']} stroked segments -> "
              f"{v['merged_vertices']} merged vertices -> "
              f"{v['carbons']} carbons "
              f"(blocks.py says {blocks.MONOMERS[code]['n_carbons']})")
    print(f"\nTable S4 carbon arithmetic: {agree}/{len(arith)} peptoids agree "
          "on all three of (lipid monomers, lipid carbons, branched/unsat)")
    print(f"monomers used by 222-266   : {len(used)} -> {', '.join(used)}")
    print(f"\nwrote {out}")
    bad = [r["peptoid"] for r in arith if not r["agree"]]
    if bad or any(v["carbons"] != blocks.MONOMERS[c]["n_carbons"]
                  for c, v in vec.items()):
        print(f"\nDISAGREEMENTS: {bad}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
