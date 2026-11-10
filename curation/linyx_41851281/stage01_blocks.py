#!/usr/bin/env python3
"""Building blocks: 14 heads, 16 tails, 2 linkers -> 01_blocks.csv.

TWO INDEPENDENT METHODS, as curation/lib/chains requires.

  GEOMETRY  stage02_trace.py skeletonises Supplementary Fig. 1a and counts
            RDP vertices. The squiggle that marks the attachment bond is
            separated from the chain by BOND LENGTH, not by eye: the drawing's
            bonds measure 40-47 px and the squiggle's wiggles 5-10 px, with
            nothing between, so a >=25 px threshold is a gap, not a guess.
  HRMS      Supplementary Notes (MOESM1 pp.41-53) print `Mass calculated m/z`
            for 118 lipids. Subtracting the head's exact mass and dividing by
            the arm count gives the acrylate/acrylamide monomer's exact mass,
            hence its formula, hence the tail's carbon count -- with no
            picture involved. 22 of 22 (linker, tail) groups come out with a
            spread of 0.0000 Da across every lipid that uses them.

Mass cannot see connectivity, so the BRANCH POSITIONS of T10-T13 and the
spacer lengths of T14-T16 come from the trace, and the trace's totals are
checked against the formula. Where they disagree, chains.count_carbons raises.
"""
import csv
import json
import os
import re
import sys

import openpyxl
import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import chains                                              # noqa: E402

SI = os.path.join(HERE, "source", "41551_2026_1640_MOESM1_ESM.pdf")

EM = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048,
      "O": 15.9949146196, "S": 31.97207100}


def exact(formula):
    tot = 0.0
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if el:
            tot += EM[el] * (int(n) if n else 1)
    return tot


# ---------------------------------------------------------------- heads
# SMILES: connectivity read off Supplementary Fig. 1a (MOESM1 p.5, 2026-11-10).
# FORMULA: independently confirmed by HRMS -- see hrms_masses() and the
# back-solve in 07_qc.json. 13 of 14 heads reproduce to <0.0001 Da; H8 is
# confirmed by lipid 48 and its partner lipid 47's PRINTED mass is 0.0070 Da
# low (08_contradictions.csv C-03).
HEADS = {
    # name: (SMILES, formula, n_NH_sites, drawn_spacers, note)
    "H1":  ("OCCCC[NH2:1]", "C4H11NO", 2, [4],
            "4-aminobutan-1-ol; one primary amine, the OH does not undergo "
            "aza-Michael and is retained (mass confirms)."),
    "H2":  ("CN(C)CCC[NH2:1]", "C5H14N2", 2, [3],
            "N,N-dimethylpropane-1,3-diamine."),
    "H3":  ("OCC[NH:1]CC[NH2:1]", "C4H12N2O", 3, [2, 2],
            "2-[(2-aminoethyl)amino]ethan-1-ol."),
    "H4":  ("[NH2:1]CCN(C)CC[NH2:1]", "C5H15N3", 4, [2, 2],
            "N,N-bis(2-aminoethyl)methylamine."),
    "H5":  ("[NH2:1]CCN(CC)CC[NH2:1]", "C6H17N3", 4, [2, 2],
            "N,N-bis(2-aminoethyl)ethylamine."),
    "H6":  ("C[NH:1]CCN(C)CC[NH:1]C", "C7H19N3", 2, [2, 2],
            "N,N',N''-trimethyldiethylenetriamine; two secondary N-H."),
    "H7":  ("CN(C)CC[NH:1]CC[NH2:1]", "C6H17N3", 3, [2, 2],
            "N,N-dimethyldiethylenetriamine."),
    "H8":  ("CN(C)CCC[NH:1]CCC[NH2:1]", "C8H21N3", 3, [3, 3],
            "N,N-dimethyldipropylenetriamine."),
    "H9":  ("[NH2:1]CCCN(C)CCC[NH2:1]", "C7H19N3", 4, [3, 3],
            "N,N-bis(3-aminopropyl)methylamine."),
    "H10": ("[NH2:1]CCN1CC[NH:1]CC1", "C6H15N3", 3, [2],
            "1-(2-aminoethyl)piperazine."),
    "H11": ("[NH2:1]CCN(CCCO)CC[NH2:1]", "C7H19N3O", 4, [2, 2, 3],
            "3-{bis(2-aminoethyl)amino}propan-1-ol. The SI names it "
            "'2,2'-(hydroxyazanediyl)bis(ethan-1-amine)', an N-hydroxylamine, "
            "which contradicts both its own synthesis (trimethylene "
            "bromohydrin alkylates the central N) and the drawn structure; "
            "HRMS C7H19N3O settles it -- 08_contradictions.csv C-04."),
    "H12": ("CC(C)(C)OC(=O)N1CCN(CC[NH2:1])CC1", "C11H23N3O2", 2, [2],
            "tert-butyl 4-(2-aminoethyl)piperazine-1-carboxylate (Boc)."),
    "H13": ("[NH2:1]CCN(CC[NH2:1])CCCN1CCCCCC1", "C13H30N4", 4, [2, 2, 3],
            "N,N-bis(2-aminoethyl)-3-(azepan-1-yl)propan-1-amine."),
    "H14": ("[NH2:1]CCN(CC[NH2:1])CCCN1CCN(C)CC1", "C12H29N5", 4, [2, 2, 3],
            "N,N-bis(2-aminoethyl)-3-(4-methylpiperazin-1-yl)propan-1-amine."),
}

# ---------------------------------------------------------------- tails
# R-group SMILES with the attachment atom atom-mapped [*:3].
# `geom` is the carbon count the trace returns (see 02_trace.json); `branch`
# records how the trace partitions it at the branch point.
TAILS = {
    "T1":  ("[CH3:3]CCCCCCC/C=C\\C/C=C\\CCCCC", "C18H33", 18,
            "linoleyl, (9Z,12Z)-octadeca-9,12-dien-1-yl; the two horizontal "
            "bonds in the trace sit at C9-C10 and C12-C13 and both flanking "
            "vertices lie on the same side, so both alkenes are Z."),
    "T2":  ("[CH3:3]CCCCCCC/C=C\\CCCCCCCC", "C18H35", 18,
            "oleyl, (9Z)-octadec-9-en-1-yl; one horizontal bond at C9-C10, Z."),
    "T3":  ("[CH3:3]CCCCCCCCCCCCCCCCC", "C18H37", 18, "n-octadecyl (stearyl)."),
    "T4":  ("[CH3:3]CCCCCCCCCCCCCCC", "C16H33", 16, "n-hexadecyl (cetyl)."),
    "T5":  ("[CH3:3]CCCCCCCCCCCCC", "C14H29", 14, "n-tetradecyl."),
    "T6":  ("[CH3:3]CCCCCCCCCCC", "C12H25", 12, "n-dodecyl."),
    "T7":  ("[CH3:3]CCCCCCCCC", "C10H21", 10, "n-decyl."),
    "T8":  ("[CH3:3]CCCCCCC", "C8H17", 8, "n-octyl."),
    "T9":  ("[CH3:3]CCCCC", "C6H13", 6, "n-hexyl."),
    "T10": ("CC[CH2:3]CCCCCCCC", "C11H23", 11,
            "undecan-3-yl; the trace gives 2 carbons on one arm of the branch "
            "point and 8 on the other, +1 for the branch carbon. SECONDARY "
            "ester -- the 1H NMR O-CH of every OT10 lipid is at 4.8 ppm, not "
            "the 4.05 ppm triplet of the primary esters."),
    "T11": ("CCCCCC[CH2:3]CCCCCCCC", "C15H31", 15,
            "pentadecan-7-yl; trace 6 and 8 either side of the branch."),
    "T12": ("CCCCCCCC[CH2:3]CCCCCCCC", "C17H35", 17,
            "heptadecan-9-yl; trace 8 and 8, symmetric."),
    "T13": ("CC[CH2:3]CCCCCCCCCCCCC", "C16H33", 16,
            "hexadecan-3-yl; trace 2 and 13."),
    "T14": ("[CH3:3]CSSCCCCCCCCCCCCCC", "C16H33S2", 16,
            "2-(tetradecyldisulfanyl)ethyl. Confirmed in words by the SI's own "
            "synthesis: tetradecane-1-thiol + 2-mercaptoethanol + acryloyl "
            "chloride, 'Synthesis of 2-(tetradecyldisulfaneyl)ethyl acrylate'."),
    "T15": ("[CH3:3]CSCCCCCCCCCCCCCCC", "C17H35S", 17,
            "2-(pentadecylsulfanyl)ethyl. The alkyl is C15, ONE carbon longer "
            "than T14's and T16's C14 -- the trace says 15 long bonds and the "
            "HRMS monomer 342.2593 Da says C20H38O2S, independently."),
    "T16": ("[CH3:3]COCCCCCCCCCCCCCC", "C16H33O", 16,
            "2-(tetradecyloxy)ethyl."),
}
# NOTE on T16: the trace gives 2 carbons before the drawn O and 14 after.

LINKERS = {"O": ("ester", "[CH2:2]CC(=O)O{R}"),
           "N": ("amide", "[CH2:2]CC(=O)N{R}")}


def hrms_masses():
    d = pymupdf.open(SI)
    t = re.sub(r"\s+", " ", "".join(d[i].get_text() for i in range(40, 53)))
    out = {}
    pat = (r"Lipid (P?\d+|T\d+)\s*(?:\([^)]*\))?\s*: (.*?)"
           r"(?=Lipid (?:P?\d+|T\d+)\s*(?:\([^)]*\))?\s*:|$)")
    for m in re.finditer(pat, t):
        mm = re.search(r"Mass calculated m/z ([\d.]+)", m.group(2))
        if mm:
            out[m.group(1)] = float(mm.group(1))
    return out


# Two ML IDs where the correspondence table (Supplementary Data 4) and the
# PRINTED name table (Supplementary Fig. 1b) disagree, and HRMS settles it.
# Read off the rendered Fig. 1b on 2026-11-10 and confirmed by exact mass:
# 1064.9096 = H3 + 3 x OT1 and 1061.9575 = H3 + 3 x NT1, both to 0.0001 Da,
# whereas 2 arms is 320.27 / 319.29 Da short. 08_contradictions.csv C-02.
NAME_CORRECTIONS = {"5": ("H3-OT1(2)", "H3-OT1(3)"),
                    "6": ("H3-NT1(2)", "H3-NT1(3)")}

# One printed exact mass that is wrong. Its partner lipid 48 (H8-NT1(3),
# 1117.0361) reproduces to 0.0000 Da from the same head, so H8 = C8H21N3 is
# pinned and lipid 47's printed 1119.9811 is 0.0070 Da below the enumerated
# 1119.9881. 08_contradictions.csv C-03. Excluded from the monomer solve;
# NOT corrected in the paper's text and NOT used anywhere else.
MASS_EXCLUDE = {"47"}


def names():
    """ML-ID -> lipid name, with the two Fig.-1b corrections applied."""
    d = json.load(open(os.path.join(HERE, "00_names.json")))
    for k, (was, now) in NAME_CORRECTIONS.items():
        if d[k] != was:
            raise SystemExit("correction %s expected %r, found %r" % (k, was, d[k]))
        d[k] = now
    return d


def monomer_formula(linker, tail):
    """CH2=CH-C(=O)X-R as a formula string."""
    tf = TAILS[tail][1]
    c = int(re.search(r"C(\d+)", tf).group(1))
    h = int(re.search(r"H(\d+)", tf).group(1))
    s = re.search(r"S(\d*)", tf)
    o = re.search(r"O(\d*)", tf)
    ns = (int(s.group(1) or 1) if s else 0)
    no = (int(o.group(1) or 1) if o else 0)
    if linker == "O":
        f = "C%dH%dO%d" % (c + 3, h + 3, no + 2)
    else:
        f = "C%dH%dN1O%d" % (c + 3, h + 4, no + 1)
    if ns:
        f += "S%d" % ns
    return f


def main():
    masses = hrms_masses()
    nm = names()
    # independent HRMS solve of every monomer
    from collections import defaultdict
    solved = defaultdict(list)
    for mlid, name in nm.items():
        m = re.match(r"H(\d+)-([ON])T(\d+)\((\d)\)$", name)
        M = masses.get(mlid)
        if not m or M is None or mlid in MASS_EXCLUDE:
            continue
        h, lk, tl, n = "H" + m.group(1), m.group(2), "T" + m.group(3), int(m.group(4))
        solved[(lk, tl)].append(((M - exact(HEADS[h][1])) / n, mlid))

    rows = []
    for h, (smi, formula, sites, spacers, note) in HEADS.items():
        rows.append(dict(block=h, kind="head", smiles=smi, formula=formula,
                         exact_mass="%.5f" % exact(formula),
                         n_NH_sites=sites, geometry=";".join(map(str, spacers)),
                         hrms_check="", note=note))
    for t, (smi, tf, geom, note) in TAILS.items():
        for lk in ("O", "N"):
            key = (lk, t)
            if key not in solved:
                continue
            vals = [v for v, _ in solved[key]]
            f = monomer_formula(lk, t)
            calc = exact(f)
            spread = max(vals) - min(vals)
            if spread > 1e-3:
                raise SystemExit("monomer %s%s spread %.4f Da over %d lipids"
                                 % (lk, t, spread, len(vals)))
            dev = sum(vals) / len(vals) - calc
            if abs(dev) > 5e-4:
                raise SystemExit("monomer %s%s: HRMS %.4f vs formula %s %.4f"
                                 % (lk, t, sum(vals) / len(vals), f, calc))
            # the formula's carbon count for R, independent of the picture
            hrms_c = int(re.search(r"C(\d+)", f).group(1)) - 3
            # chains.count_carbons with one route returns `agreed: False` and
            # records which route was unavailable -- here the vector one, by
            # construction (the page has no vector paths). The SECOND method
            # is HRMS, which is not a geometry route and must not be passed as
            # `vector=`; a lie in that field would make the record claim a
            # measurement that was never taken. It is asserted explicitly.
            agree = chains.count_carbons(raster=geom, label="%s%s" % (lk, t))
            if agree["agreed"] or "vector" in agree["methods"]:
                raise SystemExit("chains API changed under this file")
            if hrms_c != geom:
                raise SystemExit(
                    "%s%s: the trace says R = %d C and the HRMS formula %s "
                    "says R = %d C. A chain length is an integer; one method "
                    "is wrong. Re-render, check the bbox, read the overlay."
                    % (lk, t, geom, f, hrms_c))
            rows.append(dict(
                block="%s%s" % (lk, t), kind="tail_monomer",
                smiles=LINKERS[lk][1].format(R=""), formula=f,
                exact_mass="%.5f" % calc, n_NH_sites="",
                geometry="R=%dC (trace)" % geom,
                hrms_check="R=%dC from %d lipids, spread %.4f Da, dev %+.5f Da"
                           % (hrms_c, len(vals), spread, dev),
                note=TAILS[t][3]))
            del agree
    for t, (smi, tf, geom, note) in TAILS.items():
        rows.append(dict(block=t, kind="tail_R", smiles=smi, formula=tf,
                         exact_mass="", n_NH_sites="",
                         geometry="%dC" % geom, hrms_check="", note=note))

    out = os.path.join(HERE, "01_blocks.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["block", "kind", "smiles", "formula",
                                           "exact_mass", "n_NH_sites",
                                           "geometry", "hrms_check", "note"])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print("wrote", out, len(rows), "rows")


if __name__ == "__main__":
    main()
