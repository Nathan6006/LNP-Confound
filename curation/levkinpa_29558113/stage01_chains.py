#!/usr/bin/env python3
"""Stage 01 -- measure every building block's chain length. NEVER by eye.

Scheme S1 (SI p.2) is VECTOR: pymupdf reports 358 stroked paths at width
0.776 pt on that page, so the drawing's own bond coordinates are available and
curation/lib/chains.py's primary (vector) route applies. CLAUDE.md section 5:
"Name tokens are not chain lengths" -- T18/T16/T14 and PY18/PY16/PY14/PY12 are
compared against the measurement by chains.check_name_token, which reports and
never corrects.

A SECOND, GEOMETRY-INDEPENDENT COUNT comes from the printed HRMS molecular
formula: for 2-(alkyldisulfanyl)pyridine the formula is C(5+n)H(2n+5)NS2, so
the carbon count fixes n with no pixels and no coordinates involved. The two
routes are reported separately and compared; this is the "two independent
methods must agree or it raises" contract, with mass standing in for the raster
route (which is unavailable and unnecessary -- the page is not flattened).

MASS CANNOT SEE CONNECTIVITY. The formula fixes how many carbons are in the
chain; it says nothing about whether the chain is linear or where the double
bond sits. Linearity and the oleoyl 9,10 position come from the drawn geometry
below, not from the mass.

Run:  python3 -I stage01_chains.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import pymupdf                                              # noqa: E402
import chains                                               # noqa: E402

SI = os.path.join(HERE, "source", "bc8b00007_si_001.pdf")
PAGE = 1                                                    # 0-based; SI p.2

# Panel bounds, in PDF points, read off the page's own text anchors (the word
# boxes for "PY18, 4 =", "T18, 2 =" etc. are at the y values below) and then
# widened to the full page width on the open side. These are this paper's own
# geometry and so live here, not in curation/lib.
THIOL_BANDS = {                       # panel b), the HS-CnH(2n+1) reagents
    "PY18": (340, 152, 612, 162),
    "PY16": (340, 163, 612, 173),
    "PY14": (340, 177, 612, 187),
    "PY12": (340, 190, 612, 200),
}
ACYL_BANDS = {                        # panel a), the acid chlorides Cl-CO-R
    "T18": (136, 154, 340, 170),
    "T16": (136, 173, 340, 189),
    "T14": (136, 191, 340, 206),
}
OLEOYL_BAND = (166, 85, 292, 97)     # compound 2: the oleoyl chain ONLY -- x starts
# right of the Cl-C1 bond and y stops above the carbonyl O, so the measured
# vertices are carbons and nothing else. Both endpoints of a segment must lie
# inside the box (chains.vector_bonds), so the Cl-C1 and C=O strokes drop out.

# HRMS, transcribed from the SI text (pdfplumber, tier 4). Stored as the
# CARBON COUNT of the printed formula; the H count is NOT used because the SI
# is inconsistent about whether it prints the neutral or the [M+H]+ hydrogen
# count (see 08_contradictions.csv, C-3).
HRMS_C = {"PY18": 23, "PY16": 21, "PY14": 19, "PY12": 17,   # C(5+n)
          "T18": 22, "T16": 20, "T14": 18, "T18U": 22}      # C(4+n)
HRMS_MASS = {"PY18": 396.2760, "PY16": 368.2447, "PY14": 340.2134,
             "PY12": 312.1821, "T18": 384.2938, "T16": 356.2625,
             "T14": 328.2312, "T18U": 382.2782}


def name_token(block):
    """The integer in the label. T18U -> 18, PY12 -> 12."""
    return int("".join(c for c in block if c.isdigit()))


def main():
    doc = pymupdf.open(SI)
    pg = doc[PAGE]
    out = {
        "_what": "Chain length per building block, measured from the PDF's own "
                 "bond coordinates (curation/lib/chains.py vector route) and "
                 "cross-checked against the printed HRMS carbon count.",
        "page": PAGE + 1,
        "stroked_paths_on_page": sum(
            1 for p in pg.get_drawings() if p.get("width")),
        "blocks": {},
        "name_token_checks": [],
    }

    for blk, bbox in THIOL_BANDS.items():
        segs = chains.vector_bonds(pg, bbox)
        v = chains.vector_chain_length(segs)
        # HS-C1...Cn drawn as n bonds (S-C1 plus n-1 C-C) and n+1 vertices,
        # because the S sits at the start of the first stroked segment while
        # "HS" itself is a text label.
        n_vec = len(segs)
        n_mass = HRMS_C[blk] - 5                 # C(5+n)H(2n+5)NS2
        res = chains.count_carbons(vector=n_vec, name_token=name_token(blk),
                                   label=blk)
        out["blocks"][blk] = {
            "kind": "alkanethiol -> S-S-alkyl arm",
            "bbox_pt": bbox, "n_stroked_segments": n_vec,
            "n_merged_vertices": v,
            "carbons_from_geometry": n_vec,
            "carbons_from_HRMS_formula": n_mass,
            "geometry_and_mass_agree": n_vec == n_mass,
            "HRMS_mz_printed": HRMS_MASS[blk],
            "count_carbons": res,
        }
        out["name_token_checks"].append(
            chains.check_name_token(name_token(blk), n_vec, label=blk))

    for blk, bbox in ACYL_BANDS.items():
        segs = chains.vector_bonds(pg, bbox)
        v_nd = chains.vector_chain_length(segs, drop_double=False)
        v = chains.vector_chain_length(segs)
        # Cl-C1(=O)-C2...Cn drawn as: 1 Cl-C bond + 1 C=O (two parallel
        # strokes, one removed by drop_double) + (n-1) C-C bonds.
        # vertices after dedup = Cl + O + n carbons.
        n_vec = v - 2
        n_mass = HRMS_C[blk] - 4                 # thiolactone C4H5NOS + acyl Cn
        res = chains.count_carbons(vector=n_vec, name_token=name_token(blk),
                                   label=blk)
        out["blocks"][blk] = {
            "kind": "acyl chloride -> amide tail",
            "bbox_pt": bbox, "n_stroked_segments": len(segs),
            "n_vertices_no_double_dedup": v_nd,
            "n_merged_vertices": v,
            "carbons_from_geometry": n_vec,
            "carbons_from_HRMS_formula": n_mass,
            "geometry_and_mass_agree": n_vec == n_mass,
            "HRMS_mz_printed": HRMS_MASS[blk],
            "count_carbons": res,
        }
        out["name_token_checks"].append(
            chains.check_name_token(name_token(blk), n_vec, label=blk))

    # --- T18U: the oleoyl chain of compound 2, and WHERE its double bond is ---
    segs = chains.vector_bonds(pg, OLEOYL_BAND)
    verts = sorted({(round(x, 1), round(y, 1))
                    for s in segs for x, y in ((s[0], s[1]), (s[2], s[3]))})
    # Collapse the double bond's inner stroke, then walk left to right.
    walk = []
    for x, y in verts:
        if walk and abs(x - walk[-1][0]) < 2.0 and abs(y - walk[-1][1]) < 2.5:
            continue
        walk.append((x, y))
    n_vec = len(walk)
    # WHERE the double bond sits. A cis (Z) bond in this drawing is the one
    # HORIZONTAL pair of strokes in an otherwise strictly alternating zig-zag;
    # its position is read by counting chain vertices to its left. Measured, not
    # assumed: "oleoyl" is a name and names are not measurements.
    horiz = [s for s in segs if abs(s[1] - s[3]) < 0.6]
    if len(horiz) != 2:
        raise SystemExit(
            f"stage01: expected exactly 2 horizontal strokes (the cis double "
            f"bond's outer and inner line) in the oleoyl chain, found "
            f"{len(horiz)}. Refusing to guess the unsaturation position.")
    xl = min(min(s[0], s[2]) for s in horiz)
    n_left = sum(1 for x, _ in walk if x < xl - 0.5)
    # n_left counts the carbons strictly LEFT of the bond; the bond itself joins
    # the next two, so C(n_left+1)=C(n_left+2). An off-by-one here would put the
    # unsaturation on the wrong carbon and nothing downstream would see it --
    # the formula is identical for every positional isomer.
    dbl = [n_left + 1, n_left + 2]
    if dbl != [9, 10]:
        print(f"stage01 WARNING: oleoyl double bond measured at C{dbl[0]}=C{dbl[1]}, "
              "not the C9=C10 of oleic acid. Do not silently accept this.")
    out["blocks"]["T18U"] = {
        "kind": "oleoyl chloride -> mono-unsaturated amide tail",
        "bbox_pt": OLEOYL_BAND,
        "n_stroked_segments": len(segs),
        "chain_vertices_left_to_right": n_vec,
        "carbons_from_geometry": n_vec,
        "carbons_from_HRMS_formula": HRMS_C["T18U"] - 4,
        "geometry_and_mass_agree": n_vec == HRMS_C["T18U"] - 4,
        "HRMS_mz_printed": HRMS_MASS["T18U"],
        "n_horizontal_strokes": len(horiz),
        "double_bond_between_carbons": dbl,
        "double_bond_geometry": "cis/Z -- drawn as a HORIZONTAL bond joining "
                                "two up-vertices of the zig-zag, which is the "
                                "standard ChemDraw cis depiction. Corroborated "
                                "by 13C 129.99/129.73 and 1H 5.36-5.27 (m, 2H).",
        "count_carbons": chains.count_carbons(
            vector=n_vec, name_token=18, label="T18U"),
    }
    out["name_token_checks"].append(
        chains.check_name_token(18, n_vec, label="T18U"))

    out["summary"] = {
        "all_geometry_mass_agree": all(
            b["geometry_and_mass_agree"] for b in out["blocks"].values()),
        "name_token_mismatches": [c["label"] for c in out["name_token_checks"]
                                  if not c["matches"]],
    }
    bad = [k for k, b in out["blocks"].items()
           if not b["geometry_and_mass_agree"]]
    if bad:
        raise SystemExit(
            f"stage01: geometry and HRMS disagree on {bad}. A chain length is "
            "an integer and two independent methods differing means one is "
            "wrong, not that the answer is between them "
            "(curation/lib/chains.count_carbons' contract). Refusing to emit.")
    with open(os.path.join(HERE, "02_chain_counts.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out["summary"], indent=1))
    for k, b in out["blocks"].items():
        print(f"{k:6s} geometry={b['carbons_from_geometry']:3d}  "
              f"HRMS={b['carbons_from_HRMS_formula']:3d}  "
              f"agree={b['geometry_and_mass_agree']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
