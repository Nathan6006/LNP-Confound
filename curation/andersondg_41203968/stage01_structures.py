#!/usr/bin/env python3
"""Reconstruct every AMG lipid the SI DRAWS, and gate each one on its own HRMS.

    python3 -I stage01_structures.py

Reads  source/41565_2025_2044_MOESM1_ESM.pdf  (pages 45-121, "Spectral
characterizations of the AMG ionizable lipids").
Writes 02_structures.csv  and  02_hrms.json.

The gate is the SI's printed HRMS m/z, not the printed formula string. The m/z
is a measurement; the formula beside it is a transcription, and at least one of
them (AMG101/AMG102) is the NEUTRAL formula labelled "[M+H]+", which would make
a formula-string comparison fail on a structure that is right. A compound whose
reconstruction does not reproduce its own m/z is recorded as FAIL and is not
carried forward.

Mass cannot see connectivity. It does see the hydrogen count, so it separates a
monocyclic headgroup from a spiro one, and it cannot tell 1,4,7-triazonane from
1-(2-aminoethyl)piperazine (both C6H15N3). Regiochemistry is read from the
drawing's own coordinates and from the synthesis in stage02.
"""
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import pymupdf                                            # noqa: E402
import chemdraw_vector as cdv                             # noqa: E402
from rdkit import Chem                                    # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt             # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula    # noqa: E402

SI = os.path.join(HERE, "source", "41565_2025_2044_MOESM1_ESM.pdf")
PAGE_FIRST, PAGE_LAST = 45, 121        # the characterization section
PROTON = 1.00727646
PPM_GATE = 25.0                        # printed m/z are 4 dp; 25 ppm is loose
TAG_MAX_PT = 70.0                      # a tag further than this is not a caption

LAB = re.compile(r"^(AMG\s?[0-9][0-9A-Za-z\-]*)\s*$")
# "calculated for" loses its "for" once (AMG49). A regex that required it
# dropped that compound's only external check -- and its structure matches the
# formula exactly, so the loss was silent.
HRMS = re.compile(r"HRMS[^\n]*?m/z\s*([0-9.]+)\s*\(calculated(?:\s+for)?\s*"
                  r"([A-Za-z0-9]+)\s*\[(M[^\]]*)\]")


def read_hrms(doc):
    """{label: {mz, formula_ion, ion, page}} from the SI's own text layer.

    The compound HEADING and the drawing's own tag carry the same string. They
    are told apart by font: the heading is Times/Cambria body text, the tag
    inside the drawing is Arial or Helvetica. Keying on the text alone would
    make every drawing tag reset the current compound.

    An HRMS line is paired with the nearest heading ON ITS OWN PAGE, preferring
    one above it, and only falls back to the last heading seen on an earlier
    page when its page carries none. "The last heading seen anywhere" put
    AMG1641's HRMS on AMG1541-OHD, because that compound's heading sits BELOW
    its own data -- and the mis-pairing then showed up as a failed gate on a
    structure that was right.
    """
    out, carried = {}, None
    for pno in range(PAGE_FIRST - 1, PAGE_LAST + 2):
        heads, hits = [], []
        seq = 0
        for blk in doc[pno - 1].get_text("dict")["blocks"]:
            if blk["type"] != 0:
                continue
            for ln in blk["lines"]:
                seq += 1
                txt = "".join(s["text"] for s in ln["spans"])
                if all(s["font"].startswith("Arial")
                       or s["font"].startswith("Helvetica")
                       for s in ln["spans"]):
                    continue
                m = LAB.match(txt.strip())
                if m:
                    heads.append((seq, m.group(1).replace(" ", "")))
                    continue
                h = HRMS.search(txt)
                if h:
                    hits.append((seq, h))
        for seq, h in hits:
            above = [lab for s2, lab in heads if s2 < seq]
            lab = above[-1] if above else carried
            if lab:
                out.setdefault(lab, {
                    "label": lab, "page": pno, "mz": float(h.group(1)),
                    "formula_ion": h.group(2), "ion": h.group(3)})
        if heads:
            carried = heads[-1][1]
    return out


def page_structures(page):
    """[(tag, smiles, formula, mass, n_atoms, dist_pt)] for one page.

    One tag takes one component: greedy nearest-pair matching, not
    nearest-tag-per-component. Nearest-tag-per-component let one tag claim four
    fragments of a reaction scheme while the structure beside it got none.
    """
    segs, _ = cdv.bond_segments(page)
    atoms, tags = cdv.glyph_boxes(page)
    if not segs or not tags:
        return []
    nodes, edges = cdv.graph(segs, atoms)
    comps = [c for c in cdv.components(nodes, edges) if len(c) >= 5]
    pairs = []
    for ci, c in enumerate(comps):
        xs = [nodes[i][1] for i in c]
        ys = [nodes[i][2] for i in c]
        bx = (min(xs), min(ys), max(xs), max(ys))
        for ti, (t, tb) in enumerate(tags):
            cx, cy = (tb[0] + tb[2]) / 2, (tb[1] + tb[3]) / 2
            dx = max(bx[0] - cx, 0.0, cx - bx[2])
            dy = max(bx[1] - cy, 0.0, cy - bx[3])
            pairs.append((math.hypot(dx, dy), ci, ti))
    pairs.sort()
    usedc, usedt, out = set(), set(), []
    for d, ci, ti in pairs:
        if d > TAG_MAX_PT or ci in usedc or ti in usedt:
            continue
        usedc.add(ci)
        usedt.add(ti)
        c = comps[ci]
        try:
            mol = cdv.to_mol(nodes, edges, c)
            out.append((tags[ti][0].replace(" ", ""), Chem.MolToSmiles(mol),
                        CalcMolFormula(mol), ExactMolWt(mol), len(c), d))
        except Exception as exc:                           # noqa: BLE001
            out.append((tags[ti][0].replace(" ", ""), "", f"FAIL:{exc}",
                        float("nan"), len(c), d))
    return out


def main():
    doc = pymupdf.open(SI)
    hr = read_hrms(doc)
    json.dump(hr, open(os.path.join(HERE, "02_hrms.json"), "w"), indent=1)

    rows = []
    for pno in range(PAGE_FIRST, PAGE_LAST + 1):
        for tag, smi, formula, mass, n, d in page_structures(doc[pno - 1]):
            rows.append({"label": tag, "page": pno, "n_heavy_atoms": n,
                         "tag_distance_pt": round(d, 1), "SMILES": smi,
                         "formula": formula, "exact_mass": mass})

    # --- the gate ---------------------------------------------------------
    by = {}
    for r in rows:
        by.setdefault(r["label"], []).append(r)
    npass = nfail = nohrms = 0
    for r in rows:
        h = hr.get(r["label"])
        r["hrms_mz"] = h["mz"] if h else ""
        r["hrms_ion"] = h["ion"] if h else ""
        r["hrms_formula_as_printed"] = h["formula_ion"] if h else ""
        if not h or not r["SMILES"]:
            r["ppm"] = ""
            r["status"] = "NO_HRMS" if not h else "UNPARSEABLE"
            nohrms += 1
            continue
        z = 2 if "2H" in h["ion"] else 1
        calc = (r["exact_mass"] + z * PROTON) / z
        ppm = (h["mz"] - calc) / calc * 1e6
        r["ppm"] = round(ppm, 1)
        r["status"] = "PASS" if abs(ppm) <= PPM_GATE else "FAIL"
        npass += r["status"] == "PASS"
        nfail += r["status"] == "FAIL"

    import csv
    cols = ["label", "page", "status", "ppm", "hrms_mz", "hrms_ion",
            "hrms_formula_as_printed", "formula", "exact_mass",
            "n_heavy_atoms", "tag_distance_pt", "SMILES"]
    with open(os.path.join(HERE, "02_structures.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in sorted(rows, key=lambda x: (x["page"], x["label"])):
            w.writerow({k: r.get(k, "") for k in cols})

    missing = sorted(set(hr) - {r["label"] for r in rows})
    print(f"HRMS entries in SI      : {len(hr)}")
    print(f"structures reconstructed: {len(rows)}")
    print(f"  PASS (<= {PPM_GATE:g} ppm)     : {npass}")
    print(f"  FAIL                  : {nfail}")
    print(f"  no HRMS / unparseable : {nohrms}")
    print(f"HRMS labels with no reconstructed structure: {len(missing)}")
    if missing:
        print("   ", missing)
    dup = [k for k, v in by.items() if len(v) > 1]
    if dup:
        print("labels reconstructed more than once:", dup)
    return 0


if __name__ == "__main__":
    sys.exit(main())
