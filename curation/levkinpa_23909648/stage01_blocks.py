#!/usr/bin/env python3
"""Stage 01 -- the library's BLOCKS, from Table S1 (SI pp.4-5), tier-4 text.

WHAT TABLE S1 IS. "ESI-MS of typical lipidoids chosen from the library": a
14-column text table indexed by amine number (1-20) and alkyl halide chain
length, giving calculated and observed m/z SEPARATELY FOR EACH NUMBER OF ALKYL
TAILS, 1 through 6.

WHAT IT PROVES, AND IT IS THE MOST IMPORTANT FACT IN THIS FOLDER.
One (amine, alkyl halide) pair carries up to FIVE different tail counts. So a
"lipidoid" named e.g. '1C14' is NOT one compound: it is the crude product of a
one-pot alkylation, a MIXTURE of mono-, di-, tri-, tetra- and penta-alkylated
species. The paper says so itself on p.7: "The 1C12-D product was isolated from
the crude lipidoid mixture by column chromatography". Trap 1 of the brief, in
its strongest form: 200 library members is 200 REACTION PRODUCTS, not 200
distinct structures.

WHAT IT DOES NOT CONTAIN. Any amine structure. The 'Amine' column holds the
bare integers 1-20. Only amine 1 is identified anywhere in this SI (p.3's
synthesis names 1-(2-aminoethyl)pyrrolidine, and the three purified hits'
IUPAC names on pp.7/10/13 all carry '(2-(pyrrolidin-1-yl)ethyl)').

THE UNPRINTED QUANTITY THIS STAGE CHECKS. The authors never printed the
difference between consecutive tail counts. Each added alkyl tail replaces one
N-H with one n-alkyl, so the calculated m/z must rise by exactly one CnH2n --
196.2191 Da for C14, 168.1878 for C12. That increment is external to anything
we compute from their numbers, it is pure chemistry, and it tests their calc
column 43 ways. It also independently CONFIRMS that the extra rows are
successive alkylations rather than, say, adducts or dimers.

Run:  python3 -I stage01_blocks.py
"""
import csv
import json
import os
import re
import sys

import pdfplumber
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "bc400158w_si_001.pdf")

# Table S1's 14-column rule, measured. Page index -> (rule y, body y-range).
#
# NOTE ON THE y=145.70 RULE. It is the SECOND header rule and it sits ABOVE
# the 'Amine | Alkyl | calc. | observed' sub-header line (y=147.5), not below
# it -- the sub-header's own second line ('halides') is at y=155.1. So the rule
# gives the column bands but is NOT the body boundary, and the body must start
# below the sub-header at y=160. The first version of this file took the rule
# as the body boundary and the sub-header line raised as an orphan, which is
# the refusal doing its job: a parser that had silently dropped it would have
# read the whole table with no complaint.
S1_PAGES = [(3, 145.70, (160.0, 760.0)), (4, 309.80, (55.0, 309.80))]
COLS = (["amine", "alkyl"]
        + [f"{k}_{n}" for n in range(1, 7) for k in ("calc", "obs")])
WATERMARK_X = 560.0
ROW_TOL = 3.0
PROTON = 1.007276            # the proton, for [M+H]+ / [M+2H]2+

# The one amine this SI identifies. Written from the chemical name printed on
# SI p.3, NOT taken from any product string: 1-(2-aminoethyl)pyrrolidine.
AMINE_1 = dict(
    number=1,
    name="1-(2-aminoethyl)pyrrolidine",
    source="SI p.3, 'Synthesis of lipidoid 1C12': 'The 1-(2-aminoethyl)"
           "pyrrolidine (1 eq., 0.25 mmol, 31.7 ul) and 1-iodododecane "
           "(2 eq., 0.5 mmol, 123.3 ul) were added to the vial.'",
    # Primary N tagged :1 for assemble.attach. Two alkylation sites (it is a
    # primary amine), which is why the crude product is a mono/di mixture.
    smiles_mapped="C1CCN(C1)CC[NH2:1]",
    smiles="NCCN1CCCC1",
    corroboration="the three purified hits' IUPAC names (SI pp.7/10/13) all "
                  "contain '(2-(pyrrolidin-1-yl)ethyl)', which is this amine "
                  "minus the alkylated primary N's hydrogens")


def column_bands(page, rule_y, tol=2.0, min_w=5.0):
    segs = sorted((round(r["x0"], 2), round(r["x1"], 2)) for r in page.rects
                  if (r["y1"] - r["y0"]) <= 1.0
                  and abs(r["top"] - rule_y) <= tol
                  and (r["x1"] - r["x0"]) >= min_w)
    if not segs:
        raise RuntimeError(f"no rule rects near y={rule_y}")
    return segs


def text_lines(page, bands, y0, y1):
    """Words -> lines of cells, placed by measured x. Straddles raise."""
    ws = [w for w in page.extract_words()
          if y0 < w["top"] <= y1 and w["x0"] < WATERMARK_X]
    lines = {}
    for w in ws:
        hit = [i for i, (bx0, bx1) in enumerate(bands)
               if w["x0"] >= bx0 - 1.0 and w["x1"] <= bx1 + 1.0]
        if len(hit) != 1:
            over = [i for i, (bx0, bx1) in enumerate(bands)
                    if w["x1"] > bx0 and w["x0"] < bx1]
            if len(over) > 1:
                raise RuntimeError(
                    f"Table S1: {w['text']!r} straddles columns {over}")
            continue
        key = next((k for k in lines if abs(k - w["top"]) <= ROW_TOL), None)
        if key is None:
            key = w["top"]
            lines[key] = [[] for _ in bands]
        lines[key][hit[0]].append((w["x0"], w["text"]))
    return [(round(k, 2), [" ".join(t for _, t in sorted(c)) for c in v])
            for k, v in sorted(lines.items())]


def main():
    doc = pdfplumber.open(PDF)
    lines = []
    for pi, rule_y, (y0, y1) in S1_PAGES:
        pg = doc.pages[pi]
        bands = column_bands(pg, rule_y)
        if len(bands) != 14:
            raise RuntimeError(
                f"Table S1 p.{pi + 1}: {len(bands)} column bands, expected 14 "
                "(amine, alkyl, and calc/observed for 1-6 tails).")
        for y, cells in text_lines(pg, bands, y0, y1):
            lines.append((pi + 1, y, cells))
    doc.close()

    # --- group lines into logical rows -------------------------------------
    # A logical row opens on a line carrying an amine number in column 0. Its
    # continuation lines hold the charge-state annotations ('(M+H+)') and, for
    # amine 11 and 16, a second [M+2H]2+ series. Grouping is structural; a
    # line that belongs to no row raises rather than being dropped.
    rows = []
    for pageno, y, c in lines:
        if re.fullmatch(r"\d+", c[0].strip()):
            rows.append(dict(page=pageno, y=y, amine=int(c[0]),
                             alkyl=c[1].strip(), main=c, extra=[]))
        else:
            if not rows:
                raise RuntimeError(f"Table S1: orphan line at p{pageno} y={y}")
            rows[-1]["extra"].append((pageno, y, c))

    # 23 is MEASURED, not expected: my first guess here was 25 and the
    # assertion caught it. The claims worth asserting are the library's SHAPE
    # (20 amines, 10 alkyl halides, 200 members), which is what the tracker
    # asserts and what the manuscript would cite -- not the number of rows the
    # authors chose to show ESI-MS for.
    if len(rows) != 23:
        raise RuntimeError(
            f"Table S1: {len(rows)} (amine, alkyl) rows, measured 23.")

    # --- the library's SHAPE ----------------------------------------------
    amines = sorted({r["amine"] for r in rows})
    alkyls = sorted({r["alkyl"] for r in rows},
                    key=lambda a: int(a.lstrip("C")))
    if amines != list(range(1, 21)):
        raise RuntimeError(
            f"Table S1: amine numbers {amines} are not exactly 1-20. The "
            "library's first dimension is what makes 20 x 10 = 200, so a gap "
            "here would mean the tracker's 200 rests on something else.")
    if len(alkyls) != 10:
        raise RuntimeError(
            f"Table S1: {len(alkyls)} distinct alkyl halides {alkyls}, "
            "expected 10.")
    shape = dict(
        distinct_amine_numbers=amines,
        n_amines=len(amines),
        distinct_alkyl_halides=alkyls,
        n_alkyl_halides=len(alkyls),
        product=len(amines) * len(alkyls),
        tracker_claim=200,
        agrees=len(amines) * len(alkyls) == 200,
        what_the_200_ARE="200 crude one-pot REACTION PRODUCTS, each a mixture "
                         "of mono- through penta-alkylated species. NOT 200 "
                         "distinct structures and NOT 200 individually "
                         "reported measurements.")

    # --- the UNPRINTED increment check ------------------------------------
    # m/z for tail count k -> k+1 must rise by exactly one CnH2n.
    def cnh2n(n):
        m = Chem.MolFromSmiles("C" * n)          # n-alkane, CnH(2n+2)
        return ExactMolWt(m) - 2 * 1.0078250319  # minus 2 H -> CnH2n

    inc = []
    for r in rows:
        n = int(r["alkyl"].lstrip("C"))
        want = cnh2n(n)
        vals = {}
        for k in range(1, 7):
            s = r["main"][2 + 2 * (k - 1)].strip()
            if s:
                vals[k] = float(s)
        ks = sorted(vals)
        for a, b in zip(ks, ks[1:]):
            if b != a + 1:
                continue                          # not consecutive
            got = vals[b] - vals[a]
            inc.append(dict(amine=r["amine"], alkyl=r["alkyl"],
                            tails=f"{a}->{b}", observed_step=round(got, 4),
                            expected_CnH2n=round(want, 4),
                            residual=round(got - want, 4)))
    worst = max(inc, key=lambda d: abs(d["residual"]))
    # Their calc column is printed to 1 decimal, so rounding alone permits
    # +/-0.1; allow 0.15 and report the worst rather than asserting a pass.
    inc_check = dict(
        n_consecutive_pairs=len(inc),
        tolerance_Da=0.15,
        rationale="their calc column is printed to one decimal, so two "
                  "roundings permit +/-0.1 before any chemistry is wrong",
        n_within_tolerance=sum(1 for d in inc if abs(d["residual"]) <= 0.15),
        worst=worst,
        all_pairs=inc)

    # --- amine 1: OUR exact masses vs THEIR calc column -------------------
    # Independent: we build the molecule and compute ExactMolWt; they printed a
    # number. Three-way split per the brief is in stage07_qc.py; here we only
    # record our side for amine 1's rows.
    a1 = Chem.MolFromSmiles(AMINE_1["smiles"])
    a1_mass = ExactMolWt(a1)
    # THE SI USES TWO DIFFERENT MASS CONVENTIONS AND LABELS ONLY ONE.
    # Table S1's '(M+H+)' annotation sits in the OBSERVED column only, and
    # measurement says that is literal: their 'calc.' column is the NEUTRAL
    # monoisotopic mass, not the protonated ion. Amine 1 + C12, 1 tail: their
    # calc 282.3 vs our neutral 282.3035 (d = +0.0035) and vs our [M+H]+
    # 283.3108 (d = -1.0073, i.e. exactly one proton). The HRMS section on
    # pp.9/12/15 uses the OTHER convention -- 'calculated m/z = 451.4986' for
    # 1C12-D is [M+H]+, matching our computed [M+H]+ to -0.00002 Da.
    #
    # This was caught because the first version of this check added a proton to
    # compare against Table S1's calc column and reported a 1.0421 Da residual.
    # Had the tolerance been 2 Da, every row would have "passed" one proton
    # out. Recorded in 08_contradictions.csv as C-2.
    ours = []
    for r in rows:
        if r["amine"] != 1:
            continue
        n = int(r["alkyl"].lstrip("C"))
        for k in range(1, 7):
            s = r["main"][2 + 2 * (k - 1)].strip()
            if not s:
                continue
            neutral = a1_mass + k * cnh2n(n)
            obs = float(r["main"][3 + 2 * (k - 1)].strip())
            ours.append(dict(
                amine=1, alkyl=r["alkyl"], tails=k,
                our_neutral_exact=round(neutral, 4),
                our_MH_plus=round(neutral + PROTON, 4),
                their_calc=float(s),
                their_observed=obs,
                # (a) ours vs their calc, on THEIR convention (neutral)
                residual_our_neutral_vs_their_calc=round(neutral - float(s), 4),
                # (b) ours vs their observed ion
                residual_our_MH_vs_their_observed=round(
                    neutral + PROTON - obs, 4),
                # (c) their own calc vs their own observed: must be ~1 proton
                their_obs_minus_their_calc=round(obs - float(s), 4)))

    # --- what the 'amine' integers can and cannot tell us -----------------
    # Solving each amine's neutral mass from the calc column is a real,
    # unprinted quantity -- but a MASS IS NOT A STRUCTURE, and this stage
    # refuses to turn one into the other.
    solved = {}
    for r in rows:
        n = int(r["alkyl"].lstrip("C"))
        for k in range(1, 7):
            s = r["main"][2 + 2 * (k - 1)].strip()
            if not s:
                continue
            m = float(s) - PROTON - k * cnh2n(n)
            solved.setdefault(r["amine"], []).append(round(m, 2))
    amine_masses = {}
    for a, ms in sorted(solved.items()):
        amine_masses[a] = dict(
            solved_neutral_masses=ms,
            spread=round(max(ms) - min(ms), 3),
            mean=round(sum(ms) / len(ms), 2),
            structure_known=(a == 1),
            note=("amine 1 = 1-(2-aminoethyl)pyrrolidine, named on SI p.3"
                  if a == 1 else
                  "STRUCTURE UNKNOWN. A solved mass constrains the formula "
                  "and nothing more; this curation does not guess a "
                  "structure from a mass. DECISIONS.md D-1."))

    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["block_type", "block_id", "IL_name_token", "structure_known",
                    "SMILES", "source", "note"])
        w.writerow(["amine_head", "amine_1", "1", "yes", AMINE_1["smiles"],
                    "SI p.3 synthesis of lipidoid 1C12",
                    "1-(2-aminoethyl)pyrrolidine; a PRIMARY amine, so the "
                    "one-pot alkylation gives a mono/di mixture"])
        for a in amines:
            if a == 1:
                continue
            w.writerow(["amine_head", f"amine_{a}", str(a), "no", "",
                        "SI Table S1 'Amine' column (an integer only)",
                        "NOT IN THIS SI. Solved neutral mass "
                        f"{amine_masses[a]['mean']} Da; a mass is not a "
                        "structure. DECISIONS.md D-1."])
        for alk in alkyls:
            n = int(alk.lstrip("C"))
            w.writerow(["alkyl_halide", alk, alk, "yes",
                        "C" * n + "I",
                        "SI Table S1 'Alkyl halides' column; SI p.3 names "
                        "1-iodododecane for C12",
                        f"1-iodo-n-alkane, {n} carbons, unbranched. Adds "
                        f"C{n}H{2 * n} per alkylation."])

    out = dict(
        _what="Stage 01. The library's blocks and shape, from Table S1's text "
              "layer. Tier 4; no pixels.",
        library_shape=shape,
        amine_1=AMINE_1,
        amine_1_our_mass_vs_their_calc=ours,
        unprinted_increment_check=inc_check,
        solved_amine_neutral_masses=amine_masses,
        multiple_alkylation_evidence=dict(
            rows_with_more_than_one_tail_count=sum(
                1 for r in rows
                if sum(1 for k in range(1, 7)
                       if r["main"][2 + 2 * (k - 1)].strip()) > 1),
            max_tail_count_observed=max(
                max((k for k in range(1, 7)
                     if r["main"][2 + 2 * (k - 1)].strip()), default=0)
                for r in rows),
            paper_says="SI p.7: 'The 1C12-D product was isolated from the "
                       "crude lipidoid mixture by column chromatography'"),
        rows=[dict(page=r["page"], y=r["y"], amine=r["amine"],
                   alkyl=r["alkyl"], cells=r["main"],
                   extra_lines=[c for _, _, c in r["extra"]]) for r in rows])
    with open(os.path.join(HERE, "01_blocks.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"OK  Table S1: {len(rows)} (amine, alkyl) rows")
    print(f"    library shape: {shape['n_amines']} amines x "
          f"{shape['n_alkyl_halides']} alkyl halides = {shape['product']} "
          f"(tracker says {shape['tracker_claim']}: "
          f"{'AGREES' if shape['agrees'] else 'DISAGREES'})")
    print(f"    amines {amines}")
    print(f"    alkyl halides {alkyls}")
    print(f"    multiple alkylation: max {out['multiple_alkylation_evidence']['max_tail_count_observed']} "
          f"tails; {out['multiple_alkylation_evidence']['rows_with_more_than_one_tail_count']}"
          f"/{len(rows)} rows show >1 tail count")
    print(f"    UNPRINTED increment check: "
          f"{inc_check['n_within_tolerance']}/{inc_check['n_consecutive_pairs']}"
          f" consecutive tail steps equal one CnH2n within 0.15 Da; "
          f"worst residual {worst['residual']} Da "
          f"(amine {worst['amine']} {worst['alkyl']} {worst['tails']})")
    print(f"    amine 1, 3-way mass split over {len(ours)} Table S1 points:")
    print(f"      (a) our neutral vs their calc.      max |d| "
          f"{max(abs(d['residual_our_neutral_vs_their_calc']) for d in ours):.4f} Da")
    print(f"      (b) our [M+H]+ vs their observed    max |d| "
          f"{max(abs(d['residual_our_MH_vs_their_observed']) for d in ours):.4f} Da")
    print(f"      (c) their observed - their calc.    "
          f"{[d['their_obs_minus_their_calc'] for d in ours]} "
          f"(one proton = {PROTON:.4f})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
