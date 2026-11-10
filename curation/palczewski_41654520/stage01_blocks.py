#!/usr/bin/env python3
"""01_blocks.csv -- the paper's OWN two-position decomposition.

    python3 -I stage01_blocks.py

Read off Supplementary Figure 1 (SI p2), which states it as a Markush table
"CBB Structure = (Lipid1, Lipid2, Linker)" over 12 lipid blocks and 3 linkers.

THIS IS NOT A COMBINATORIAL ENUMERATION and the file says so on its own rows.
12 x 12 x 3 would be 432; the paper made 21 of those cells and bought a 22nd
compound, and every one is individually named, purified and characterised. So
"stated vs enumerated" is 21 drawn cells vs 21 built, not a block product.
"""
import csv
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
import paper_library as P                                       # noqa: E402

def main():
    rows = []
    for code, desc in P.DRAWN_LIPIDS.items():
        used = sorted((c for c, t in P.DRAWN_TABLE.items()
                       if code in (t[0], t[1])),
                      key=lambda s: int(s[3:]))
        rows.append(dict(position="1_lipid", block_id=code,
                         block_description=desc,
                         structure_or_tail_id=P.LIPID_TO_TAIL[code],
                         smiles_fragment=P.ACYL.get(P.LIPID_TO_TAIL[code],
                                                    ("", ""))[0],
                         n_compounds=len(used), compounds=" ".join(used),
                         source="Supplementary Figure 1 'Lipids' box, SI p2"))
    for code, desc in P.DRAWN_LINKERS.items():
        used = sorted((c for c, t in P.DRAWN_TABLE.items() if t[2] == code),
                      key=lambda s: int(s[3:]))
        rows.append(dict(position="2_linker", block_id=code,
                         block_description=desc,
                         structure_or_tail_id=P.LINKER_TO_SCAFFOLD[code],
                         smiles_fragment="", n_compounds=len(used),
                         compounds=" ".join(used),
                         source="Supplementary Figure 1 'Linkers' box, SI p2"))
    rows.append(dict(position="0_scaffold", block_id="CBB-G250",
                     block_description="Coomassie Brilliant Blue G-250 "
                                       "triarylmethane dye, both benzenesulfonyl "
                                       "groups derivatised",
                     structure_or_tail_id="dye core",
                     smiles_fragment=P.dye(P.HEAD_ETHOXY, "{a}", "{b}"),
                     n_compounds=22, compounds="CBB1..CBB22",
                     source="Supplementary Figure 1 scheme, SI p2 (drawn in "
                            "full) + the IUPAC names on SI pp 18-50"))
    rows.append(dict(position="0_scaffold", block_id="CBBZ-G250-propargyl",
                     block_description="the same dye with 4-(prop-2-yn-1-"
                                       "yloxy)aniline in place of 4-ethoxy"
                                       "aniline, both sulfonic acids left free",
                     structure_or_tail_id="CBBZ dye core",
                     smiles_fragment=P.dye(P.HEAD_PROPYNYL, "O", "O"),
                     n_compounds=11, compounds="CBBZ0..CBBZ10",
                     source="Supplementary Figure 3 + SI p42 CBBZ0 synthesis"))
    for zn, (tail, why) in sorted(P.AZIDE_TAIL.items(),
                                  key=lambda kv: int(kv[0][1:])):
        rows.append(dict(position="1_lipid", block_id=zn + "-azide",
                         block_description=why,
                         structure_or_tail_id="triazole-" + zn,
                         smiles_fragment="[N-]=[N+]=N" + tail,
                         n_compounds=1, compounds="CBBZ" + zn[1:],
                         source="SI pp 31-42, intermediates I5-I24"))
    cols = ["position", "block_id", "block_description",
            "structure_or_tail_id", "smiles_fragment", "n_compounds",
            "compounds", "source"]
    with open(os.path.join(_HERE, "01_blocks.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("blocks: %d lipid, %d linker, %d scaffold; drawn cells %d"
          % (sum(r["position"] == "1_lipid" for r in rows),
             sum(r["position"] == "2_linker" for r in rows),
             sum(r["position"] == "0_scaffold" for r in rows),
             len(P.DRAWN_TABLE)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
