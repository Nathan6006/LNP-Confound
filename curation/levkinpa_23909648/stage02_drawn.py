#!/usr/bin/env python3
"""Stage 02 -- GEOMETRY. Count the bonds the authors DREW and compare the
drawn topology to the structures we built.

WHY THIS EXISTS. "Passing every printed statistic is necessary, not
sufficient" -- QL_2026 reproduced every printed hit rate on bars up to 28% too
high. Geometry has to be checked geometrically. Every OTHER confirmation in
this folder is a number the authors printed: an IUPAC name, an HRMS pair, a
formula, an exact mass, an NMR integral. All of those are mutually consistent
arithmetic over the same claimed formula, and all of them would still agree if
the picture showed something else.

WHY NOT chains.vector_bonds. It was tried first and returned 0 bonds on all
three pages, which is the library refusing correctly rather than inventing a
number: it reads STROKED paths, and ChemDraw exported these bonds as FILLED
OUTLINES. Measured on SI p.10: `page.get_drawings()` returns 24 objects, every
one of type 'f' (filled), each a thin closed quadrilateral -- two long 'l'
edges capped by two short 'c' curves -- whose rect is e.g. 0.6 x 14.2 pt. So
ONE FILLED OBJECT IS ONE BOND, and the count is read directly. This is a
paper-specific encoding, which PIPELINE_V2 section 3 puts in the curation
folder next to the figure, not in curation/lib.

The ~76 small images on each page are NOT bonds. They are the text glyphs for
the compound label and the Chemical Formula / Exact Mass / Molecular Weight
block, which is why none of that text appears in the PDF text layer.

THE INVARIANT. For a connected molecule, bonds = heavy atoms - 1 + rings.
Each compound here has exactly one ring (the pyrrolidine), so the drawn bond
count is an independent witness on BOTH the chain lengths and the number of
tails: a C12 tail mistaken for C14 moves it by 2, and a mono- mistaken for a
di-alkylamine moves it by 16. Nothing else in this folder would notice either.

Run:  python3 -I stage02_drawn.py
"""
import json
import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import chains                                                # noqa: E402
import pymupdf                                               # noqa: E402
from rdkit import Chem, RDLogger                             # noqa: E402

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "bc400158w_si_001.pdf")

DRAWN = [
    (6, "1C12-D", "CCCCCCCCCCCCN(CCCCCCCCCCCC)CCN1CCCC1"),
    (9, "1C16-S", "CCCCCCCCCCCCCCCCNCCN1CCCC1"),
    (12, "1C16-D", "CCCCCCCCCCCCCCCCN(CCCCCCCCCCCCCCCC)CCN1CCCC1"),
]
# Generous and FIXED, covering the whole upper half of the page. It is not
# fitted to the expected answer -- error class 27's signature is a fit whose
# x-vector is the quantity being validated.
BBOX = (60.0, 100.0, 560.0, 320.0)
VERTEX_TOL_PT = 2.6     # endpoints closer than this are the same atom


def bond_axis(d):
    """The two endpoints of a filled bond outline: its most distant points."""
    pts = []
    for it in d["items"]:
        for p in it[1:]:
            if hasattr(p, "x"):
                pts.append((p.x, p.y))
            elif hasattr(p, "x0"):
                pts.append((p.x0, p.y0))
                pts.append((p.x1, p.y1))
    best, pa, pb = -1.0, None, None
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            dd = math.dist(pts[i], pts[j])
            if dd > best:
                best, pa, pb = dd, pts[i], pts[j]
    return pa, pb, best


def main():
    doc = pymupdf.open(PDF)
    out = []
    for pi, label, smi in DRAWN:
        page = doc[pi]

        # BOTH library modes are RUN and RECORDED, so every claim about what
        # they return is executed rather than asserted.
        #
        # include_filled=False -> 0 segments: these bonds are not strokes.
        # include_filled=True  -> each filled bond outline yields its TWO long
        #   'l' edges, so the count must be EXACTLY 2 x the bond count. The
        #   library's docstring warns that filled mode also returns GLYPH
        #   OUTLINES as bonds, and that they must be separated by stroke
        #   thickness, aspect ratio or label proximity. NO SEPARATION IS
        #   NEEDED HERE, and that is measured rather than assumed: on these
        #   three pages every glyph ('N', 'H', the compound label, the whole
        #   Chemical Formula / Exact Mass / Molecular Weight block) is a
        #   bitmap IMAGE, not a path -- which is why none of that text is in
        #   the PDF text layer. The exact 2:1 ratio is itself the evidence
        #   that nothing else leaked in: a glyph outline would contribute an
        #   irregular number of segments and break the identity.
        stroked = chains.vector_bonds(page, BBOX)
        filled_segs = chains.vector_bonds(page, BBOX, include_filled=True)

        x0, y0, x1, y1 = BBOX
        filled = [d for d in page.get_drawings()
                  if d["type"] == "f"
                  and x0 <= d["rect"].x0 and d["rect"].x1 <= x1
                  and y0 <= d["rect"].y0 and d["rect"].y1 <= y1]

        axes = [bond_axis(d) for d in filled]
        lens = sorted(round(a[2], 2) for a in axes)

        # --- vertex graph, by merging coincident bond endpoints -----------
        verts = []
        deg = []
        for pa, pb, _ in axes:
            for p in (pa, pb):
                hit = None
                for k, v in enumerate(verts):
                    if math.dist(p, v) <= VERTEX_TOL_PT:
                        hit = k
                        break
                if hit is None:
                    verts.append(p)
                    deg.append(1)
                else:
                    deg[hit] += 1

        m = Chem.MolFromSmiles(smi)
        n_heavy = m.GetNumAtoms()
        n_rings = m.GetRingInfo().NumRings()
        n_bonds = m.GetNumBonds()
        our_deg = sorted((a.GetDegree() for a in m.GetAtoms()), reverse=True)

        # --- THE CLOSED IDENTITY ------------------------------------------
        # ChemDraw emits one path per heavy-atom bond plus one more per
        # double bond (the second line of the double). Kekulise so aromatic
        # bonds are counted as drawn, not as a bond order of 1.5. This
        # chemistry is fully saturated, so n_double is 0 and the identity
        # reduces to drawn == heavy bonds -- but the term is computed, not
        # dropped, so the check stays correct if it is ever reused on a paper
        # with an unsaturated tail.
        km = Chem.Mol(m)
        Chem.Kekulize(km, clearAromaticFlags=True)
        n_double = sum(1 for b in km.GetBonds()
                       if b.GetBondType() == Chem.BondType.DOUBLE)
        expected_paths = km.GetNumBonds() + n_double

        out.append(dict(
            page=pi + 1, label=label,
            drawn_filled_bond_objects=len(filled),
            chains_vector_bonds_stroked=len(stroked),
            chains_vector_bonds_filled_segments=len(filled_segs),
            filled_segments_per_bond=(
                round(len(filled_segs) / len(filled), 4) if filled else None),
            filled_segments_identity_holds=(
                len(filled_segs) == 2 * len(filled)),
            n_double_bonds_kekulised=n_double,
            expected_drawn_paths="heavy bonds + double bonds (kekulised)",
            expected_drawn_paths_value=expected_paths,
            closed_identity_holds=len(filled) == expected_paths,
            glyph_separation_mode="NONE REQUIRED -- measured, not assumed: "
                                  "all glyphs on these pages are bitmap "
                                  "images, not paths, so filled mode returns "
                                  "bond outlines only. The exact 2:1 "
                                  "segments-per-bond ratio is the evidence.",
            why_stroked_is_zero="the bonds are FILLED outlines, not strokes; "
                                "curation/lib/chains reads strokes and "
                                "correctly declines to report a number here",
            our_heavy_atoms=n_heavy, our_rings=n_rings,
            our_bond_count=n_bonds,
            invariant="bonds = heavy atoms - 1 + rings",
            invariant_check=n_heavy - 1 + n_rings == n_bonds,
            bond_count_agrees=len(filled) == n_bonds,
            bond_length_pt=dict(n=len(lens), min=lens[0] if lens else None,
                                median=lens[len(lens) // 2] if lens else None,
                                max=lens[-1] if lens else None),
            drawn_vertices=len(verts),
            drawn_degree_sequence=sorted(deg, reverse=True),
            our_degree_sequence=our_deg,
            degree_note="the drawn degree sequence is DESCRIPTIVE, not "
                        "asserted: bonds stop short of the 'N' and 'H' glyphs "
                        "(which are images, not text), so endpoints at a "
                        "heteroatom do not merge cleanly and that atom is "
                        "split into separate vertices. The BOND COUNT is the "
                        "assertion; the degree sequence is reported so a "
                        "reviewer sees exactly how far the topology was "
                        "actually checked.",
            n_images_on_page=len(page.get_image_info())))
    doc.close()

    n_ok = sum(1 for r in out if r["bond_count_agrees"]
               and r["closed_identity_holds"]
               and r["filled_segments_identity_holds"])
    res = dict(
        _what="Stage 02. Drawn-bond geometry for the three purified "
              "structures, from the PDF's own filled bond outlines.",
        bbox_used=BBOX, vertex_tol_pt=VERTEX_TOL_PT,
        n_agree=n_ok, n_total=len(out), structures=out)
    with open(os.path.join(HERE, "02_drawn_geometry.json"), "w") as fh:
        json.dump(res, fh, indent=1)

    print(f"OK  drawn-bond geometry: {n_ok}/{len(out)} agree")
    for r in out:
        print(f"    p{r['page']:<3} {r['label']:<7} drawn filled bonds "
              f"{r['drawn_filled_bond_objects']:>3}  vs our bond count "
              f"{r['our_bond_count']:>3}  "
              f"({r['our_heavy_atoms']} heavy atoms - 1 + {r['our_rings']} "
              f"ring)  {'AGREE' if r['bond_count_agrees'] else 'DISAGREE'}"
              f"; closed identity (heavy bonds + "
              f"{r['n_double_bonds_kekulised']} double) = "
              f"{r['expected_drawn_paths_value']} "
              f"{'OK' if r['closed_identity_holds'] else 'FAIL'}")
        print(f"            chains.vector_bonds: stroked="
              f"{r['chains_vector_bonds_stroked']}, "
              f"include_filled={r['chains_vector_bonds_filled_segments']} "
              f"= {r['filled_segments_per_bond']} segments/bond "
              f"{'OK' if r['filled_segments_identity_holds'] else 'FAIL'}; "
              f"bond length median {r['bond_length_pt']['median']} pt; "
              f"{r['drawn_vertices']} merged vertices")
    return 0 if n_ok == len(out) else 1


if __name__ == "__main__":
    sys.exit(main())
