#!/usr/bin/env python3
"""The building blocks: 19 amine cores, 31 thiol tails, 1 diester linker.

    python3 -I stage01_blocks.py

CHAIN LENGTHS ARE MEASURED, NOT COUNTED BY EYE. Error class 48. Figure 1b draws
all 31 tails; the 15 SC drawings segment cleanly as connected components, and
each is skeletonised, walked from one endpoint and simplified with
curation/lib/chains.raster_chain_length. For an unbranched thiol drawn as
HS-(zigzag)-Me, every drawn bond is one C-S or C-C bond and the terminal atoms
are TEXT, so the carbon count is the segment count. The name token is compared
through chains.check_name_token, which reports and never corrects: CLAUDE.md
records KS_2024's `C24` (two C4 tails) and this very library's SC1/SC2/SC3.

Only the PDF's own bond coordinates would be better (chains.vector_chain_length),
and we do not have them: the publisher PDF is subscription-only and the PMC
author-manuscript PDF sits behind a proof-of-work challenge, so the figures we
hold are the published rasters. 00_acquisition.json records that.

WHAT THE MEASUREMENT FINDS, AND WHY IT MATTERS HERE. SC1 and SC2 are BRANCHED
(a degree-3 skeleton node), so a path trace does not describe them and this
file refuses to produce a number for them. SC3 measures as THREE segments, the
same as a propyl chain and one fewer than SC4 -- so in this library SCn is NOT
n carbons for n < 4: SC1 is isobutyl (C4), SC2 isopentyl (C5), SC3 n-propyl
(C3), SC4 n-butyl (C4), and only from SC4 up does SCn mean n carbons.

THE SN AND SO TAILS ARE NOT CHAINS. They carry a tertiary amine, a ring, a
carboxylic acid or a diol, so a segment count is not their structure. They are
transcribed from the drawing, each row carries read_by = "drawn structure,
Fig 1b", and each is checked against the corpus's independently curated
KZ_2016 entry in stage03.
"""
import csv
import os
import sys

import numpy as np
from PIL import Image
from scipy import ndimage
from skimage.morphology import skeletonize

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import chains   # noqa: E402

FIG = os.path.join(HERE, "source", "figures", "fig1.jpg")
OUT_BLOCKS = os.path.join(HERE, "01_blocks.csv")
OUT_QC = os.path.join(HERE, "01_blocks_qc.json")

# Figure 1b's tail box, read off this figure and kept beside it.
TAIL_BOX = (412, 982, 528, 730)      # x0, x1, y0, y1
# y1 = 730 stops ABOVE the SN/SO drawings, which start at y ~ 738. The box is
# the SC block only: the segmenter must return exactly 15 components, and
# including the SN/SO rows returns 35 and would make the name mapping a guess.
SC_ORDER = ["SC1", "SC2", "SC3", "SC4", "SC5",
            "SC6", "SC7", "SC8", "SC9", "SC10",
            "SC11", "SC12", "SC14", "SC16", "SC18"]
SC_COLUMN_SPLIT = (520, 700)         # x boundaries between the three columns

# --- the amine cores, read from Figure 1b and Scheme 1 --------------------
# Isotope tags mark the nitrogens that carry arms. They survive canonicalisation
# and attach() does not touch them, which is how a single nitrogen can be found
# again for its SECOND arm after the first bond renumbered the molecule. They
# are stripped from the emitted structure.
CORES = {
    # name: (tagged SMILES, {tag: n_arms}, n_NH_total, source)
    "3A1":  ("[101NH2]CC1CC[102NH]CC1", {101: 2, 102: 1}, 3, "Fig 1b"),
    "3A2":  ("[101NH2]CCN1CC[102NH]CC1", {101: 2, 102: 1}, 3, "Fig 1b"),
    "3A3":  ("CN(C)CCC[101NH]CCC[102NH2]", {101: 1, 102: 2}, 3, "Fig 1b"),
    "3A4":  ("OCC[101NH]CCC[102NH2]", {101: 1, 102: 2}, 3, "Fig 1b, Scheme 1"),
    "3A5":  ("C[101NH]CCN(CC[102NH]C)CC[103NH]C", {101: 1, 102: 1, 103: 1}, 3,
             "Fig 1b"),
    "4A1":  ("[101NH2]CCCN1CCN(CCC[102NH2])CC1", {101: 2, 102: 2}, 4, "Fig 1b"),
    "4A3":  ("[101NH2]CCCN(C)CCC[102NH2]", {101: 2, 102: 2}, 4,
             "Fig 1b, Scheme 1"),
    "4A4":  ("[101NH2]CCOCCOCC[102NH2]", {101: 2, 102: 2}, 4, "Fig 1b"),
    "5A1":  ("[101NH2]CCC[102NH]CCC[103NH2]", {101: 2, 102: 1, 103: 2}, 5,
             "Fig 1b"),
    "5A2":  ("[101NH2]CC[102NH]CC[103NH]CC[104NH]CC[105NH2]",
             {101: 2, 102: 1, 103: 0, 104: 0, 105: 2}, 7, "Fig 1b"),
    "6A1":  ("[101NH2]CC[102NH]CC[103NH]CC[104NH]CC[105NH2]",
             {101: 2, 102: 1, 103: 0, 104: 1, 105: 2}, 7, "Fig 1b"),
    "5A3":  ("[101NH2]CCC[102NH]CCC[103NH]CCC[104NH2]",
             {101: 2, 102: 1, 103: 0, 104: 2}, 6, "Fig 1b"),
    "6A2":  ("[101NH2]CCC[102NH]CCC[103NH]CCC[104NH2]",
             {101: 2, 102: 1, 103: 1, 104: 2}, 6, "Fig 1b"),
    "5A4":  ("[101NH2]CCN(CC[102NH2])CC[103NH2]", {101: 1, 102: 2, 103: 2}, 6,
             "Fig 1b, Scheme 1"),
    "6A3":  ("[101NH2]CCN(CC[102NH2])CC[103NH2]", {101: 2, 102: 2, 103: 2}, 6,
             "Fig 1b"),
    # the four variants, drawn only in Scheme 1 / Figures 2b and 3a
    "3A4V": ("OCCC[101NH]CC[102NH2]", {101: 1, 102: 2}, 3, "Scheme 1, Fig 2b"),
    "3A4V2": ("OCC[101NH]CC[102NH2]", {101: 1, 102: 2}, 3, "Scheme 1, Fig 2b"),
    "4A5":  ("[101NH2]CCN(C)CC[102NH2]", {101: 2, 102: 2}, 4,
             "Scheme 1, Fig 3a"),
    "5A5":  ("[101NH2]CCCN(CCC[102NH2])CCC[103NH2]", {101: 1, 102: 2, 103: 2},
             6, "Scheme 1, Fig 3a"),
}

# --- the SN and SO tails, transcribed from Figure 1b ----------------------
# Not chains: a segment count is not their structure. Each is checked against
# the independently curated KZ_2016 corpus entry in stage03.
# CORRECTION, recorded because it is the point of the cross-check. SO2 and SO8
# were first transcribed one carbon SHORT (5 instead of 6). Both were caught by
# the InChIKey comparison against the corpus's independently curated KZ_2016
# entry in stage03, then confirmed by tracing the drawings: the SO components
# measure 7 and 5 bond segments for SO8 and SO7, and for an SO tail the
# terminal OH is a TEXT label with its C-O bond drawn, so segments = carbons
# + 1 -- 6 and 4 carbons, which is what the corpus says. The SC tails differ
# because their terminal "Me" IS a carbon, so there segments = carbons. Reading
# a chain by eye is error class 48 and it cost two blocks here.
NONCHAIN_TAILS = {
    "SN1": "CN(C)CCS", "SN2": "CCCCNCCS", "SN3": "CCN(CC)CCS",
    "SN4": "SCCN1CCCC1", "SN5": "SCCN1CCCCC1", "SN6": "SCCN1CCCCCC1",
    "SN7": "CCCN(CCC)CCS", "SN8": "CCCCN(CCCC)CCS",
    "SO1": "OC(=O)CCS", "SO2": "OC(=O)CCCCCS", "SO3": "OC(=O)CCCCCCCCCCS",
    "SO4": "OCCS", "SO5": "OCC(O)CS", "SO6": "OCCCS", "SO7": "OCCCCS",
    "SO8": "OCCCCCCS",
}
BRANCHED_SC = {"SC1": "CC(C)CS", "SC2": "CC(C)CCS"}

LINKER = {
    "name": "AEMA (G1)",
    "reagent": "2-(acryloyloxy)ethyl methacrylate",
    "smiles": "C=CC(=O)OCCOC(=O)C(C)=C",
    "arm_after_both_additions": "N-CH2CH2-C(=O)-O-CH2CH2-O-C(=O)-CH(CH3)-CH2-S-R",
    "source": "Scheme 1, and the Methods 'Synthesis of Modified Amine Core' / "
              "'Synthesis of the Ionizable Amino Lipid' paragraphs",
    "regiochemistry_note":
        "Read from the SYNTHESIS, not from a mass. The aza-Michael goes onto "
        "the ACRYLATE end (1.1 equiv per arm, BHT, 50 C) and the sulfa-Michael "
        "onto the METHACRYLATE end (DMPP, 50 C), so the thiol lands on the "
        "beta carbon and the alpha carbon keeps the methyl. Scheme 1 draws the "
        "product 3A4-SC5 with that connectivity. Molecular formula is identical "
        "for the swapped assignment, so HRMS or any mass check is blind to it "
        "-- this is the Lin regiochemistry class.",
}


def _runs(mask):
    out, s = [], None
    for i, v in enumerate(mask):
        if v and s is None:
            s = i
        elif not v and s is not None:
            out.append((s, i))
            s = None
    if s is not None:
        out.append((s, len(mask)))
    return out


def segment_tails(arr):
    """The 15 SC line drawings, as connected components of Figure 1b."""
    x0, x1, y0, y1 = TAIL_BOX
    ink = arr[y0:y1, x0:x1] < 140
    lab, n = ndimage.label(ink, structure=np.ones((3, 3)))
    out = []
    for i, sl in enumerate(ndimage.find_objects(lab), 1):
        comp = (lab[sl] == i)
        h = sl[0].stop - sl[0].start
        w = sl[1].stop - sl[1].start
        area = int(comp.sum())
        if area < 40 or w < 25 or h > 25:
            continue
        if area / float(h * w) > 0.30:        # text is dense; a line is not
            continue
        out.append({"x": x0 + sl[1].start, "y": y0 + sl[0].start,
                    "w": w, "h": h, "area": area, "mask": comp})
    # three columns, top to bottom within each
    def key(c):
        col = 0 if c["x"] < SC_COLUMN_SPLIT[0] else (
            1 if c["x"] < SC_COLUMN_SPLIT[1] else 2)
        return (col, c["y"])
    out.sort(key=key)
    return out


def trace(mask):
    """Ordered skeleton path and its branch degree. Refuses on a branch."""
    sk = skeletonize(mask)
    pts = np.argwhere(sk)
    nb = {}
    S = set(map(tuple, pts))
    for p in S:
        nb[p] = [q for q in
                 ((p[0] + dy, p[1] + dx)
                  for dy in (-1, 0, 1) for dx in (-1, 0, 1)
                  if (dy, dx) != (0, 0))
                 if q in S]
    ends = [p for p, v in nb.items() if len(v) == 1]
    branch = [p for p, v in nb.items() if len(v) > 2]
    if len(ends) != 2:
        return None, len(branch), f"{len(ends)} skeleton endpoints"
    # walk
    start = min(ends, key=lambda p: p[1])
    path, prev, cur = [start], None, start
    while True:
        nxt = [q for q in nb[cur] if q != prev]
        nxt = [q for q in nxt if q not in path[-3:]]
        if not nxt:
            break
        prev, cur = cur, nxt[0]
        path.append(cur)
    return [(p[1], p[0]) for p in path], len(branch), ""


def main():
    arr = np.asarray(Image.open(FIG).convert("L")).astype(float)
    comps = segment_tails(arr)
    if len(comps) != len(SC_ORDER):
        raise SystemExit(
            f"stage01: {len(comps)} SC line drawings segmented, expected "
            f"{len(SC_ORDER)}. REFUSING: if the segmentation is not exactly "
            "the 15 SC structures then the component-to-name mapping is "
            "guesswork, and a tail assigned to the wrong name changes 15 "
            "lipids at once.")

    qc = {"_what": "Chain lengths measured from Figure 1b by raster tracing.",
          "n_components": len(comps), "tails": []}
    measured = {}
    for name, c in zip(SC_ORDER, comps):
        path, nbranch, why = trace(c["mask"])
        token = int(name[2:])
        rec = {"tail": name, "bbox_w_px": c["w"], "bbox_h_px": c["h"],
               "area_px": c["area"], "skeleton_branch_nodes": nbranch}
        if nbranch or path is None:
            rec.update({"measured_carbons": None,
                        "why_not_measured":
                            "branched skeleton" if nbranch else why,
                        "structure_source": "drawn structure, Fig 1b"})
        else:
            v = chains.raster_chain_length(path, eps_px=3.0)
            n = v - 1                      # segments, not vertices
            got = chains.count_carbons(raster=n, name_token=token, label=name)
            rec.update({"rdp_vertices": v, "measured_carbons": got["n"],
                        "name_token_check":
                            chains.check_name_token(token, got["n"], name)})
            measured[name] = got["n"]
        qc["tails"].append(rec)

    # --- the blocks table ------------------------------------------------
    rows = []
    for name in SC_ORDER:
        if name in BRANCHED_SC:
            rows.append({"block": name, "kind": "tail",
                         "smiles": BRANCHED_SC[name],
                         "n_carbons": "", "measured": "no",
                         "read_by": "drawn structure, Fig 1b (branched: the "
                                    "skeleton has a degree-3 node, so a path "
                                    "trace does not describe it)",
                         "source": "Fig 1b"})
        else:
            n = measured[name]
            rows.append({"block": name, "kind": "tail",
                         "smiles": "C" * n + "S",
                         "n_carbons": n, "measured": "yes",
                         "read_by": f"raster trace of Fig 1b, {n} bond "
                                    "segments between the HS and Me labels",
                         "source": "Fig 1b"})
    for name, smi in NONCHAIN_TAILS.items():
        rows.append({"block": name, "kind": "tail", "smiles": smi,
                     "n_carbons": "", "measured": "no",
                     "read_by": "drawn structure, Fig 1b (not a chain: it "
                                "carries an amine, a ring, an acid or a diol)",
                     "source": "Fig 1b"})
    for name, (smi, plan, nh, src) in CORES.items():
        rows.append({"block": name, "kind": "core", "smiles": smi,
                     "n_carbons": "", "measured": "no",
                     "read_by": f"drawn structure; {sum(plan.values())} arms "
                                f"on {nh} N-H",
                     "source": src})
    rows.append({"block": LINKER["name"], "kind": "linker",
                 "smiles": LINKER["smiles"], "n_carbons": "", "measured": "no",
                 "read_by": "drawn reagent, Scheme 1", "source": "Scheme 1"})

    with open(OUT_BLOCKS, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    qc["linker"] = LINKER
    qc["name_token_mismatches"] = [
        t["tail"] for t in qc["tails"]
        if t.get("name_token_check") and not t["name_token_check"]["matches"]]
    import json
    with open(OUT_QC, "w") as fh:
        json.dump(qc, fh, indent=1)

    for t in qc["tails"]:
        print(f"  {t['tail']:<5} w={t['bbox_w_px']:>3} branch_nodes="
              f"{t['skeleton_branch_nodes']} carbons={t.get('measured_carbons')}"
              f" token_matches="
              f"{(t.get('name_token_check') or {}).get('matches')}")
    print(f"name-token mismatches: {qc['name_token_mismatches']}")
    print(f"wrote {len(rows)} blocks to {OUT_BLOCKS}")


if __name__ == "__main__":
    main()
