#!/usr/bin/env python3
"""The 48 building blocks of the PILL library, from SI page S5, VECTOR-first.

SOURCE. SI p.S5, section 2-1 "Chemical Reagents and Starting Materials": one
ChemDraw figure, three panels -- amine (A1-A26), boronic acid (B1-B8),
aldehyde (C1-C14). The page is FULLY VECTOR: 486 stroked paths carrying 483
line items, plus 219 filled paths that are glyph outlines (heteroatom labels
are drawn as outlines, not text, so the page's text layer names no block).

WHY A BOND-LINE COUNT AND NOT A CHAIN COUNT. curation/lib/chains.py exists
because counting a zig-zag by eye is this project's +/-1 failure mode (error
class 48). Its `vector_chain_length` counts vertices in one chain. This figure
lets us do something strictly stronger and uniform across all 48 blocks:

    drawn_lines(block) == heavy_atom_bonds(SMILES) + double_bonds(SMILES)

because ChemDraw draws one stroked line per heavy-atom bond and a second
parallel line for each double bond, and draws no bond to an implicit hydrogen
(the H in "NH2"/"OH" lives inside the label glyph). So every proposed structure
is tested against the drawing's own geometry, bond for bond -- branch points,
rings and double bonds included -- not just its longest chain.

THREE THINGS THIS CHECK IS BLIND TO, stated so nobody reads it as more than it
is:
  1. WHICH heteroatom a label is. N and O draw the same number of bonds. The
     element identities were read from a 600 dpi render (00_pages/).
  2. WHERE a branch sits when moving it preserves the bond count. 2-methyl vs
     3-methyl piperidine draw the same 11 lines. Also read from the render.
  3. Double-bond GEOMETRY (E/Z). Two parallel lines are two parallel lines.
     The blocks whose stereochemistry this leaves open are flagged in
     01_blocks.csv and queued in DECISIONS.md.

WHAT CLOSES THOSE GAPS for the nine blocks that actually enter a structure we
emit: the SI prints, for each of the nine purified lipids, the reagent amount
used at a stated molar amount. mass/mol (or volume x density/mol) pins the
reagent's MOLECULAR WEIGHT, independently of the drawing and of our reading of
it. See `AMOUNT_EVIDENCE` below and 07_qc.json.

Run:  python3 -I stage02_blocks.py        writes 01_blocks.csv
      python3 -I stage02_blocks.py --seed-defect <BLOCK>
                                          corrupts one block's SMILES by one
                                          carbon and confirms the check fires.
"""
import csv
import math
import os
import sys

import pymupdf
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
import chains            # noqa: E402  (vector_bonds, the PDF's own coordinates)
from errors import DetectionError, CurationLibError   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "ja6c13674_si_001.pdf")
PAGE = 4                      # 0-based; SI page S5
OUT = os.path.join(HERE, "01_blocks.csv")

# --- this paper's own geometry constants, kept here next to the figure ------
DIVIDER_MIN_LEN_PT = 100.0    # the two dashed group dividers are 261.4 pt;
                              # the longest BOND on the page is 9.0 pt.
CLUSTER_TOL_PT = 11.0         # endpoint linkage. A bond is 9.0 pt and the gap
                              # a heteroatom label leaves is ~7 pt, so 11 pt
                              # joins a molecule and 11 pt does not reach the
                              # next one (nearest inter-molecule gap: 13.0 pt,
                              # C2->C3). Asserted below, not assumed.
ANCHOR_TOL_PT = 4.0           # a cluster's (x0,y0) must land this close to its
                              # block's recorded anchor, and to exactly one.

# (label, x0, y0) of each cluster's top-left, measured from the page's own
# stroked-segment coordinates. These are POSITIONS, not measurements: they only
# say which drawing is which block. Every structural claim is checked.
ANCHORS = [
    ("A1", 186.0, 167.3), ("A2", 219.2, 168.1), ("A3", 289.6, 168.4),
    ("A4", 341.4, 168.2), ("A5", 386.7, 167.1),
    ("A6", 176.5, 218.8), ("A7", 254.4, 208.9), ("A8", 312.2, 213.6),
    ("A9", 361.0, 207.7), ("A10", 403.8, 209.6),
    ("A11", 172.2, 259.3), ("A12", 228.1, 255.4), ("A13", 276.6, 252.0),
    ("A14", 330.0, 259.6), ("A15", 404.3, 267.6),
    ("A16", 166.1, 315.5), ("A17", 199.4, 323.6), ("A18", 252.6, 324.1),
    ("A19", 321.7, 324.4), ("A20", 385.9, 319.0),
    ("A21", 170.5, 364.8), ("A22", 231.8, 356.1), ("A23", 292.5, 356.5),
    ("A24", 362.0, 357.8),
    ("A25", 241.7, 410.2), ("A26", 302.0, 404.2),
    ("B1", 103.6, 481.1), ("B2", 158.0, 480.8),
    ("B3", 106.0, 531.8), ("B4", 160.8, 531.5),
    ("B5", 106.0, 578.0), ("B6", 158.0, 580.5),
    ("B7", 106.0, 623.1), ("B8", 160.8, 620.6),
    ("C1", 214.7, 484.2), ("C2", 303.9, 482.5), ("C3", 399.4, 483.5),
    ("C4", 227.2, 514.7), ("C5", 302.6, 513.3), ("C6", 406.5, 514.5),
    ("C7", 247.8, 554.4), ("C8", 380.0, 550.4),
    ("C9", 223.4, 584.0), ("C10", 364.7, 586.1),
    ("C11", 219.5, 626.1), ("C12", 312.2, 624.0), ("C13", 227.9, 656.0),
    ("C14", 366.4, 629.3),
]

# The structures. Element identities and branch positions read from the 600 dpi
# renders in 00_pages/; every one is then tested bond-for-bond against the
# drawing by `drawn_lines == heavy_bonds + double_bonds`.
BLOCKS = {
 # --- amines, panel 1. The two dashed dividers on the page separate these
 #     into 1 reactive N-H (A1-A15), 2 (A16-A24) and 3 (A25-A26); that is the
 #     paper's own grouping and it matches the stoichiometry the SI uses
 #     (1 eq aldehyde/boronic acid for A2, 2 eq for A18).
 "A1":  ("morpholine", "C1COCCN1"),
 "A2":  ("N,N-dimethylpropane-1,3-diamine", "CN(C)CCCN"),
 "A3":  ("N,N-dimethylethane-1,2-diamine", "CN(C)CCN"),
 "A4":  ("N,N-diethylethane-1,2-diamine", "CCN(CC)CCN"),
 "A5":  ("N,N-diethylpropane-1,3-diamine", "CCN(CC)CCCN"),
 "A6":  ("N-(3-aminopropyl)-N'-(3-(dimethylamino)propyl)amine",
         "NCCCNCCCN(C)C"),
 "A7":  ("2,2'-((3-aminopropyl)azanediyl)diethanol", "NCCCN(CCO)CCO"),
 "A8":  ("2-(pyrrolidin-1-yl)ethan-1-amine", "NCCN1CCCC1"),
 "A9":  ("4-amino-2,2,6,6-tetramethylpiperidine", "CC1(C)CC(N)CC(C)(C)N1"),
 "A10": ("N,N-diisopropylethane-1,2-diamine", "CC(C)N(C(C)C)CCN"),
 "A11": ("3-(1H-imidazol-1-yl)propan-1-amine", "NCCCn1ccnc1"),
 "A12": ("2-morpholinoethan-1-amine", "NCCN1CCOCC1"),
 "A13": ("1-isopropylpiperidin-4-amine", "CC(C)N1CCC(N)CC1"),
 "A14": ("3-(2-methylpiperidin-1-yl)propan-1-amine", "CC1CCCCN1CCCN"),
 "A15": ("2-methoxyethan-1-amine", "NCCOC"),
 "A16": ("piperazine", "C1CNCCN1"),
 "A17": ("propane-1,3-diamine", "NCCCN"),
 "A18": ("pentane-1,5-diamine", "NCCCCCN"),
 "A19": ("N-(2-aminoethyl)ethane-1,2-diamine", "NCCNCCN"),
 "A20": ("N,N-bis(3-aminopropyl)methylamine", "CN(CCCN)CCCN"),
 "A21": ("N,N'-dimethylpropane-1,3-diamine", "CNCCCNC"),
 "A22": ("2,2'-(methylazanediyl)bis(ethan-1-amine)", "CN(CCN)CCN"),
 "A23": ("2,2'-(piperazine-1,4-diyl)bis(ethan-1-amine)",
         "NCCN1CCN(CCN)CC1"),
 "A24": ("N1-(2-(4-(2-aminoethyl)piperazin-1-yl)ethyl)ethane-1,2-diamine",
         "NCCN1CCN(CCNCCN)CC1"),
 "A25": ("tris(2-aminoethyl)amine", "NCCN(CCN)CCN"),
 "A26": ("tris(2-(methylamino)ethyl)amine", "CNCCN(CCNC)CCNC"),

 # --- boronic acids, panel 2 ---
 "B1": ("cyclopentylboronic acid", "OB(O)C1CCCC1"),
 "B2": ("cyclohexylboronic acid", "OB(O)C1CCCCC1"),
 "B3": ("butylboronic acid", "OB(O)CCCC"),
 "B4": ("isobutylboronic acid", "OB(O)CC(C)C"),
 "B5": ("propylboronic acid", "OB(O)CCC"),
 "B6": ("isopropylboronic acid", "OB(O)C(C)C"),
 "B7": ("hexylboronic acid", "OB(O)CCCCCC"),
 "B8": ("octylboronic acid", "OB(O)CCCCCCCC"),

 # --- aldehydes, panel 3 ---
 "C1":  ("nonanal", "CCCCCCCCC=O"),
 "C2":  ("undecanal", "CCCCCCCCCCC=O"),
 "C3":  ("dodecanal", "CCCCCCCCCCCC=O"),
 "C4":  ("(2E,6Z)-nona-2,6-dienal", r"CC/C=C\CC/C=C/C=O"),
 "C5":  ("(2E,4E)-deca-2,4-dienal", "CCCCC/C=C/C=C/C=O"),
 "C6":  ("undec-10-enal", "C=CCCCCCCCCC=O"),
 "C7":  ("(2E,6E)-3,7,11-trimethyldodeca-2,6,10-trienal (farnesal)",
         "CC(C)=CCC/C(C)=C/CC/C(C)=C/C=O"),
 "C8":  ("6-oxohexyl 4-methyloctanoate", "CCCCC(C)CCC(=O)OCCCCCC=O"),
 "C9":  ("6-oxohexyl nonanoate", "CCCCCCCCC(=O)OCCCCCC=O"),
 "C10": ("6-oxohexyl decanoate", "CCCCCCCCCC(=O)OCCCCCC=O"),
 "C11": ("decanal", "CCCCCCCCCC=O"),
 "C12": ("methyl 4-oxobutanoate", "COC(=O)CCC=O"),
 "C13": ("methyl 6-oxohexanoate", "COC(=O)CCCCC=O"),
 "C14": ("6-oxohexyl 2-hexyldecanoate",
         "CCCCCCCCC(CCCCCC)C(=O)OCCCCCC=O"),
}

# Independent, non-geometric evidence: the SI prints the amount of each reagent
# used at a stated molar amount, so amount/mol pins the MOLECULAR WEIGHT. For a
# liquid the printed volume needs a density, which is a literature constant and
# is therefore recorded as such -- it is the weakest link in these four rows and
# is labelled "volume" so a reviewer can see which rows lean on it.
# (block, SI page, printed amount, mmol, kind, density_g_per_mL_or_None)
AMOUNT_EVIDENCE = [
    ("A2",  "S12", 6.2,   0.05, "volume", 0.82),
    ("A18", "S14", 5.9,   0.05, "volume", 0.873),
    ("B1",  "S12", 5.7,   0.05, "mass",   None),
    ("B3",  "S12", 5.1,   0.05, "mass",   None),
    ("B8",  "S13", 7.9,   0.05, "mass",   None),
    ("C1",  "S12", 8.6,   0.05, "volume", 0.827),
    ("C8",  "S17", 25.6,  0.10, "mass",   None),
    ("C9",  "S18", 25.6,  0.10, "mass",   None),
    ("C14", "S19", 39.3,  0.10, "volume", 0.902),
]
AMOUNT_TOL_FRAC = 0.04        # 4 %: covers rounding of a printed microlitre
                              # and the density constant. Not a free parameter:
                              # the worst observed deviation is 1.0 %.

# Blocks whose drawn double-bond GEOMETRY the bond-line check cannot see.
STEREO_READ_FROM_RENDER = {"C4", "C5", "C7"}
# Blocks that enter a structure this folder emits.
USED = {"A2", "A18", "B1", "B3", "B8", "C1", "C8", "C9", "C14"}


def page_segments(page):
    """Every stroked bond line on the page, minus the two group dividers."""
    rect = page.rect
    segs = chains.vector_bonds(page, (rect.x0, rect.y0, rect.x1, rect.y1))
    short = [s for s in segs
             if math.hypot(s[2] - s[0], s[3] - s[1]) < DIVIDER_MIN_LEN_PT]
    dividers = len(segs) - len(short)
    if dividers != 2:
        raise DetectionError(
            "stage02: expected exactly 2 long segments (the dashed group "
            f"dividers at y=297.0 and y=396.7) but found {dividers}. The "
            "divider filter is a geometry constant of THIS figure; if the "
            "figure changed, every count below is measuring something else.")
    return short


def cluster(segs, tol):
    parent = list(range(len(segs)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def near(s, t):
        return min(math.hypot(ax - bx, ay - by)
                   for ax, ay in ((s[0], s[1]), (s[2], s[3]))
                   for bx, by in ((t[0], t[1]), (t[2], t[3])))

    for i in range(len(segs)):
        for j in range(i + 1, len(segs)):
            if near(segs[i], segs[j]) <= tol:
                a, b = find(i), find(j)
                if a != b:
                    parent[a] = b
    groups = {}
    for i, s in enumerate(segs):
        groups.setdefault(find(i), []).append(s)
    out = []
    for g in groups.values():
        xs = [v for s in g for v in (s[0], s[2])]
        ys = [v for s in g for v in (s[1], s[3])]
        out.append({"x0": min(xs), "y0": min(ys), "x1": max(xs),
                    "y1": max(ys), "n": len(g), "segs": g})
    return out


def cluster_separation(cl):
    """Smallest endpoint distance between two DIFFERENT clusters."""
    best = (1e9, None, None)
    for i in range(len(cl)):
        for j in range(i + 1, len(cl)):
            d = min(math.hypot(ax - bx, ay - by)
                    for s in cl[i]["segs"]
                    for ax, ay in ((s[0], s[1]), (s[2], s[3]))
                    for t in cl[j]["segs"]
                    for bx, by in ((t[0], t[1]), (t[2], t[3])))
            if d < best[0]:
                best = (d, i, j)
    return best


def expected_lines(smiles):
    """heavy-atom bonds + one extra line per double bond."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise CurationLibError(f"stage02: unparseable SMILES {smiles!r}")
    nb = m.GetNumBonds()
    # KEKULIZE FIRST. ChemDraw draws an aromatic ring as alternating single and
    # DOUBLE lines, so an aromatic bond order would undercount by one line per
    # formal double bond. Caught on this figure's first run: A11's imidazole
    # measured 11 lines against an aromatic-form prediction of 9.
    k = Chem.Mol(m)
    Chem.Kekulize(k, clearAromaticFlags=True)
    nd = sum(1 for b in k.GetBonds()
             if b.GetBondType() == Chem.BondType.DOUBLE)
    return nb + nd


def main():
    blocks = dict(BLOCKS)
    seeded = None
    if "--seed-defect" in sys.argv:
        seeded = sys.argv[sys.argv.index("--seed-defect") + 1]
        name, smi = blocks[seeded]
        m = Chem.MolFromSmiles(smi)
        # lengthen the molecule by one carbon: the exact +/-1 chain error this
        # check exists to catch (error class 48).
        blocks[seeded] = (name + " [SEEDED +1 CARBON]", "C" + smi)
        print("SEEDED: %s -> %s" % (seeded, blocks[seeded][1]))

    doc = pymupdf.open(SI)
    page = doc[PAGE]
    segs = page_segments(page)
    cl = cluster(segs, CLUSTER_TOL_PT)

    print("stroked bond lines on S5 (dividers removed): %d" % len(segs))
    print("clusters at tol %.1f pt: %d  (blocks declared: %d)"
          % (CLUSTER_TOL_PT, len(cl), len(blocks)))
    if len(cl) != len(blocks):
        raise DetectionError(
            f"stage02: {len(cl)} clusters for {len(blocks)} blocks. The "
            "linkage tolerance is either joining two molecules or splitting "
            "one. Refusing to guess which.")

    sep, _, _ = cluster_separation(cl)
    print("closest approach between two different clusters: %.1f pt "
          "(tolerance %.1f pt) -- the tolerance is below it, so no two "
          "molecules can have merged." % (sep, CLUSTER_TOL_PT))
    if sep <= CLUSTER_TOL_PT:
        raise DetectionError(
            f"stage02: two clusters come within {sep:.1f} pt, at or below the "
            f"{CLUSTER_TOL_PT} pt linkage tolerance. The clustering is not "
            "separating molecules; it got the right COUNT by luck.")

    # --- bijective anchor assignment -------------------------------------
    assign, used_cl = {}, set()
    for label, ax, ay in ANCHORS:
        hits = [i for i, c in enumerate(cl)
                if math.hypot(c["x0"] - ax, c["y0"] - ay) <= ANCHOR_TOL_PT]
        if len(hits) != 1:
            raise DetectionError(
                f"stage02: anchor {label} at ({ax},{ay}) matched {len(hits)} "
                "clusters, need exactly 1.")
        if hits[0] in used_cl:
            raise DetectionError(f"stage02: cluster reused by {label}.")
        used_cl.add(hits[0])
        assign[label] = cl[hits[0]]
    if len(used_cl) != len(cl):
        raise DetectionError("stage02: some cluster was never claimed.")

    # --- the check ---------------------------------------------------------
    rows, bad, total_expected = [], [], 0
    for label, (name, smi) in blocks.items():
        exp = expected_lines(smi)
        got = assign[label]["n"]
        total_expected += exp
        ok = (exp == got)
        if not ok:
            bad.append((label, name, exp, got))
        m = Chem.MolFromSmiles(smi)
        rows.append({
            "block": label,
            "panel": {"A": "amine", "B": "boronic_acid",
                      "C": "aldehyde"}[label[0]],
            "name": name,
            "SMILES": Chem.MolToSmiles(m),
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "MW": round(Descriptors.MolWt(m), 2),
            "drawn_lines_measured": got,
            "drawn_lines_expected": exp,
            "geometry_agrees": int(ok),
            "cluster_bbox": "%.1f,%.1f,%.1f,%.1f" % (
                assign[label]["x0"], assign[label]["y0"],
                assign[label]["x1"], assign[label]["y1"]),
            "used_in_an_emitted_structure": int(label in USED),
            "stereo_from_render_only": int(label in STEREO_READ_FROM_RENDER),
            "source": "SI p.S5 section 2-1 (vector)",
        })

    print("\nper-block drawn-line check: %d/%d agree"
          % (len(blocks) - len(bad), len(blocks)))
    for label, name, exp, got in bad:
        print("  MISMATCH %-4s %-45s expected %d drawn, measured %d"
              % (label, name[:45], exp, got))
    print("sum of expected lines over 48 blocks: %d   measured on page: %d"
          % (total_expected, len(segs)))

    # --- the reagent-amount cross-check -----------------------------------
    print("\nreagent amount -> molecular weight (independent of the drawing):")
    amt_bad = []
    for label, pg, amount, mmol, kind, dens in AMOUNT_EVIDENCE:
        mg = amount * (dens if kind == "volume" else 1.0)
        implied = mg / mmol
        mw = Descriptors.MolWt(Chem.MolFromSmiles(blocks[label][1]))
        dev = abs(implied - mw) / mw
        flag = "ok " if dev <= AMOUNT_TOL_FRAC else "BAD"
        if dev > AMOUNT_TOL_FRAC:
            amt_bad.append((label, implied, mw, dev))
        print("  %s %-4s %-6s %7.1f %-6s / %.3f mmol -> MW %7.2f   "
              "structure MW %7.2f   dev %5.2f %%"
              % (flag, label, pg, amount, kind, mmol, implied, mw, dev * 100))

    if seeded:
        if not bad and not amt_bad:
            raise SystemExit(
                "SEEDED DEFECT DID NOT FIRE on %s. A check that cannot fail "
                "is not a check (CLAUDE.md section 6)." % seeded)
        print("\nSEEDED DEFECT FIRED on %s -- the check has teeth." % seeded)
        return 0

    if bad or amt_bad:
        raise SystemExit("stage02: refusing to write 01_blocks.csv with a "
                         "failing check.")

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        for r in sorted(rows, key=lambda r: (r["block"][0],
                                             int(r["block"][1:]))):
            w.writerow(r)
    print("\nwrote %s (%d blocks)" % (OUT, len(rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
