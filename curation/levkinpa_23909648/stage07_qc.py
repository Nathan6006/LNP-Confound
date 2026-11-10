#!/usr/bin/env python3
"""Stage 07 -- THE GATE. Every check, with numbers, and a non-zero exit.

This file is the only thing in the folder that is allowed to say "correct". It
re-derives rather than re-reads wherever it can: the three structures are
rebuilt from the printed IUPAC names by a route that does not import
stage03_enumerate, and the masses are recomputed from those builds.

EVERY GATE REPORTS AN AGREEMENT COUNT, and 10_seed_checks.py seeds a defect
into each one and requires BOTH a non-zero exit AND the count to move. A gate
whose count does not move under its own seed is not a check -- that is the
lesson of the ZT_2026 marginals that scored a literal 0 for twelve days while
the summary beside them said "ALL FOUR EXACT", and of the sibling whose +1-unit
seed was a no-op under a loose tolerance.

Run:  python3 -I stage07_qc.py   ->  exit 0 if and only if every gate passes
"""
import csv
import json
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem.Descriptors import ExactMolWt, MolWt
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
PROTON = 1.007276
H_ATOM = 1.0078250319

failures = []
gates = {}


def gate(name, ok, count, of, detail, tolerance=None):
    gates[name] = dict(pass_=bool(ok), agreement_count=count, of=of,
                       detail=detail, tolerance=tolerance)
    if not ok:
        failures.append(f"{name}: {detail}")
    flag = "PASS" if ok else "FAIL"
    print(f"  [{flag}] {name}: {count}/{of}  {detail}")


# ---------------------------------------------------------------------------
# Re-derivation, independent of stage03. Built by plain SMILES-string
# composition of fragments each written from the printed chemical name --
# deliberately a DIFFERENT code path from stage03's atom-mapped RDKit
# assembly, so a defect in one does not hide in the other.
# ---------------------------------------------------------------------------
AMINE_1 = "NCCN1CCCC1"            # 1-(2-aminoethyl)pyrrolidine, SI p.3


def rederive(chain, n_tails):
    """N-alkylated amine 1, composed as a string from the IUPAC name."""
    alk = "C" * chain
    if n_tails == 1:
        smi = f"{alk}NCCN1CCCC1"
    elif n_tails == 2:
        smi = f"{alk}N({alk})CCN1CCCC1"
    else:
        raise ValueError(n_tails)
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise RuntimeError(f"rederive({chain},{n_tails}) unparseable")
    return Chem.MolToSmiles(m), m


def cnh2n(n):
    return ExactMolWt(Chem.MolFromSmiles("C" * n)) - 2 * H_ATOM


PURIFIED = [
    # label, chain, tails, their HRMS calc, their HRMS obs, figure 1H sum,
    # text 1H sum, DRAWN formula, DRAWN exact mass, DRAWN molecular weight
    #
    # THE LAST THREE COME FROM THE STRUCTURE DRAWING, NOT THE TEXT LAYER.
    # Each of SI pp.7/10/13 carries a ChemDraw structure flattened into ~76
    # small image strips, with the compound label, Chemical Formula, Exact
    # Mass and Molecular Weight rendered INSIDE that raster -- pdfplumber
    # returns none of them (grep 'Formula' over the text layer: 0 hits). They
    # were read visually off a 400 dpi render, which is printed numerals, the
    # best case in tier 5, and they are STRUCTURE evidence, not values: no
    # emitted value depends on them. They are a FOURTH, independent mass
    # witness, and crucially a NEUTRAL one -- the HRMS section's number is
    # [M+H]+, so this leg can fail without that one.
    ("1C12-D", 12, 2, 451.4986, 451.4987, 63.20, 59,
     "C30H62N2", 450.49130, 450.82668),
    ("1C16-S", 16, 1, 339.3734, 339.3735, 44.45, 45,
     "C22H46N2", 338.36610, 338.61404),
    ("1C16-D", 16, 2, 563.62378, 563.62344, 75.76, 78,
     "C38H78N2", 562.61650, 563.03932),
]


def main():
    print("LL_2013 (Li & Levkin, Bioconjug Chem 2013, PMID 23909648) -- "
          "Stage 07 gate\n")

    s3 = {r["IL_name"]: r for r in csv.DictReader(
        open(os.path.join(HERE, "03_enumerated.csv")))}

    # =====================================================================
    # Q-1  TABLE GEOMETRY. Did we read the tables, and read them intact?
    # =====================================================================
    geo = json.load(open(os.path.join(HERE, "05_values_geometry.json")))
    vals = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))
    want_bands = {"S3": 3, "S4": 9, "S5a": 4, "S5b": 4, "S6": 3}
    want_rows = {"S3": 13, "S4": 13, "S5a": 8, "S5b": 22, "S6": 7}
    agree = sum(1 for t, n in want_bands.items()
                if len(geo["geometry"][t]["bands"]) == n
                and geo["geometry"][t]["n_rows"] == want_rows[t])
    gate("Q-1a_table_geometry", agree == 5 and len(vals) == 87, agree, 5,
         f"{len(vals)} values; column bands from the DRAWN rules match the "
         f"declared layout on 5/5 tables; S5 page-break band drift "
         f"{geo['geometry']['S5_page_break_band_drift_pt']} pt")

    # The dose blocks: 3 blocks of 10, each DNA label on row 6 of 10.
    blocks = geo["S5_dose_blocks"]
    ok = (len(blocks) == 3
          and all(len(b["rows"]) == 10 for b in blocks)
          and all(b["label_row_of_10"] == 6 for b in blocks)
          and [b["dna"] for b in blocks] == ["75 ng", "150 ng", "200 ng"])
    gate("Q-1b_S5_dose_blocks", ok,
         sum(1 for b in blocks if b["label_row_of_10"] == 6), 3,
         "Table S5's DNA amount is a MERGED CELL; its glyph sits on row 6 of "
         "each 10-row block. The first parser treated it as the block's LAST "
         "row and stamped 9 of 30 values with the wrong dose, silently")

    # =====================================================================
    # Q-2  LIBRARY SHAPE. 20 amines x 10 alkyl halides = 200.
    # =====================================================================
    b1 = json.load(open(os.path.join(HERE, "01_blocks.json")))
    sh = b1["library_shape"]
    gate("Q-2_library_shape",
         sh["n_amines"] == 20 and sh["n_alkyl_halides"] == 10
         and sh["product"] == 200,
         sh["product"], 200,
         f"{sh['n_amines']} amines x {sh['n_alkyl_halides']} alkyl halides "
         f"= {sh['product']}; the tracker says 200. These are 200 CRUDE "
         "REACTION PRODUCTS, not 200 structures and not 200 measurements")

    # =====================================================================
    # Q-3  THE LIBRARY IS MIXTURES. Table S1 reports several tail counts for
    #      a single (amine, alkyl) pair.
    # =====================================================================
    ma = b1["multiple_alkylation_evidence"]
    gate("Q-3_members_are_mixtures",
         ma["rows_with_more_than_one_tail_count"] == 23
         and ma["max_tail_count_observed"] == 6,
         ma["rows_with_more_than_one_tail_count"], 23,
         f"every one of the 23 (amine, alkyl) rows in Table S1 reports more "
         f"than one tail count, up to {ma['max_tail_count_observed']}; the "
         "screened material is a mono/di/tri/... mixture")

    # =====================================================================
    # Q-4  AN UNPRINTED QUANTITY: the per-tail increment in THEIR calc column
    #      must equal exactly one CnH2n. Their numbers only; no input of ours.
    # =====================================================================
    inc = b1["unprinted_increment_check"]
    gate("Q-4_unprinted_tail_increment",
         inc["n_within_tolerance"] == inc["n_consecutive_pairs"] == 38,
         inc["n_within_tolerance"], inc["n_consecutive_pairs"],
         f"consecutive tail steps in the authors' own calc column equal one "
         f"CnH2n; worst residual {inc['worst']['residual']} Da (amine "
         f"{inc['worst']['amine']} {inc['worst']['alkyl']} "
         f"{inc['worst']['tails']}); their column is printed to 1 dp so two "
         "roundings permit +/-0.1", tolerance="0.15 Da")

    # =====================================================================
    # Q-5  MASS, SPLIT FOUR WAYS. Each leg can fail without the others.
    # =====================================================================
    red = []
    for (lab, chain, tails, hcalc, hobs, hfig, htext,
         dform, dexact, dmw) in PURIFIED:
        smi, m = rederive(chain, tails)
        neutral = ExactMolWt(m)
        red.append(dict(label=lab, smiles=smi, formula=CalcMolFormula(m),
                        neutral=neutral, mh=neutral + PROTON,
                        hcalc=hcalc, hobs=hobs, hfig=hfig, htext=htext,
                        drawn_formula=dform, drawn_exact=dexact, drawn_mw=dmw,
                        our_mw=MolWt(m),
                        n_h=sum(a.GetTotalNumHs() for a in m.GetAtoms())))

    # (a) OURS vs THEIR CALCULATED m/z  -- our structure against their formula
    ppm = [(r["label"], (r["mh"] - r["hcalc"]) / r["hcalc"] * 1e6)
           for r in red]
    n_ok = sum(1 for _, p in ppm if abs(p) <= 5.0)
    gate("Q-5a_ours_vs_their_calculated_mz", n_ok == 3, n_ok, 3,
         "our recomputed [M+H]+ vs the authors' CALCULATED m/z: "
         + ", ".join(f"{l} {p:+.3f} ppm" for l, p in ppm),
         tolerance="5 ppm")

    # (b) THEIR CALCULATED vs THEIR OWN OBSERVED -- internal to the paper
    mda = [(r["label"], (r["hcalc"] - r["hobs"]) * 1000) for r in red]
    n_ok = sum(1 for _, d in mda if abs(d) <= 1.0)
    gate("Q-5b_their_calc_vs_their_observed", n_ok == 3, n_ok, 3,
         "the authors' own calculated vs observed HRMS: "
         + ", ".join(f"{l} {d:+.3f} mDa" for l, d in mda),
         tolerance="1 mDa")

    # (c) OURS vs THEIR OBSERVED -- the measurement, not the arithmetic
    ppm2 = [(r["label"], (r["mh"] - r["hobs"]) / r["hobs"] * 1e6) for r in red]
    n_ok = sum(1 for _, p in ppm2 if abs(p) <= 5.0)
    gate("Q-5c_ours_vs_their_observed_mz", n_ok == 3, n_ok, 3,
         "our recomputed [M+H]+ vs their OBSERVED m/z: "
         + ", ".join(f"{l} {p:+.3f} ppm" for l, p in ppm2),
         tolerance="5 ppm")

    # (d) AN UNPRINTED RELATION AMONG THEIR THREE PRINTED HRMS VALUES.
    #     Nothing of ours enters: two differences of their own numbers must
    #     equal exact alkyl increments. This is the check CLAUDE.md section 5
    #     prefers -- a quantity the authors never printed.
    d1 = 563.62378 - 339.3734          # 1C16-D minus 1C16-S = one C16H32
    d2 = 563.62378 - 451.4986          # 1C16-D minus 1C12-D = two C4H8
    rel = [("1C16-D - 1C16-S = 1 x C16H32", d1, cnh2n(16)),
           ("1C16-D - 1C12-D = 2 x C4H8", d2, 2 * cnh2n(4))]
    n_ok = sum(1 for _, got, want in rel if abs(got - want) <= 1e-4)
    gate("Q-5d_unprinted_HRMS_relations", n_ok == 2, n_ok, 2,
         "; ".join(f"{n}: {got:.5f} vs {want:.5f} "
                   f"({(got - want) * 1000:+.4f} mDa)"
                   for n, got, want in rel)
         + " -- these only close with tail multiplicities 2, 1, 2",
         tolerance="0.1 mDa")

    # (e) Table S1's amine-1 rows: OUR NEUTRAL vs THEIR CALC (their calc
    #     column is the NEUTRAL mass -- see 08_contradictions.csv C-2).
    a1 = ExactMolWt(Chem.MolFromSmiles(AMINE_1))
    s1 = b1["amine_1_our_mass_vs_their_calc"]
    res = [d["residual_our_neutral_vs_their_calc"] for d in s1]
    n_ok = sum(1 for r in res if abs(r) <= 0.06)
    gate("Q-5e_tableS1_amine1_neutral_masses", n_ok == len(s1) == 4, n_ok,
         len(s1),
         f"our neutral exact masses vs Table S1's calc column, residuals "
         f"{res} Da (amine 1 = {a1:.4f} Da)", tolerance="0.06 Da")

    # (f) THE FORMULA AND NEUTRAL EXACT MASS THE AUTHORS DREW. Independent of
    #     every leg above: a different number, in a different place, on a
    #     different convention (neutral, not [M+H]+).
    fm = [(r["label"], r["formula"], r["drawn_formula"]) for r in red]
    n_ok = sum(1 for _, a, b in fm if a == b)
    gate("Q-5f_drawn_molecular_formula", n_ok == 3, n_ok, 3,
         "our formula vs the formula printed under the structure drawing: "
         + ", ".join(f"{l} {a} vs {b}" for l, a, b in fm))

    ex = [(r["label"], (r["neutral"] - r["drawn_exact"]) * 1000) for r in red]
    n_ok = sum(1 for _, d in ex if abs(d) <= 0.2)
    gate("Q-5g_drawn_neutral_exact_mass", n_ok == 3, n_ok, 3,
         "our NEUTRAL exact mass vs the 'Exact Mass' printed under the "
         "drawing: " + ", ".join(f"{l} {d:+.4f} mDa" for l, d in ex)
         + " (the drawing's comma is a decimal separator: 450,49130)",
         tolerance="0.2 mDa")

    mw = [(r["label"], r["our_mw"] - r["drawn_mw"]) for r in red]
    n_ok = sum(1 for _, d in mw if abs(d) <= 0.05)
    gate("Q-5h_drawn_molecular_weight", n_ok == 3, n_ok, 3,
         "our average molecular weight vs the 'Molecular Weight' printed "
         "under the drawing: "
         + ", ".join(f"{l} {d:+.4f} Da" for l, d in mw), tolerance="0.05 Da")

    # =====================================================================
    # Q-6  1H INTEGRAL SUM vs MOLECULAR FORMULA. Mass cannot see
    #      connectivity; a methyl count can see tail NUMBER.
    # =====================================================================
    fig = [(r["label"], r["hfig"] - r["n_h"]) for r in red]
    n_ok = sum(1 for _, d in fig if abs(d) <= 3.0)
    gate("Q-6a_1H_figure_integral_sum", n_ok == 3, n_ok, 3,
         "spectrum integral sum minus formula H count: "
         + ", ".join(f"{l} {d:+.2f} H" for l, d in fig), tolerance="3 H")

    txt = [(r["label"], r["htext"] - r["n_h"]) for r in red]
    n_ok = sum(1 for _, d in txt if abs(d) <= 1.0)
    gate("Q-6b_1H_text_listing_sum", n_ok == 2, n_ok, 3,
         "the SI's TEXT listing minus formula H count: "
         + ", ".join(f"{l} {d:+d} H" for l, d in txt)
         + " -- 1C12-D is short by exactly one methyl (it prints 3H where the "
           "spectrum integrates 6.38 and C30H62N2 requires 6H). THIS GATE IS "
           "EXPECTED TO AGREE ON 2 OF 3: the third is a printing error in the "
           "paper, recorded as contradiction C-3, not a structural doubt",
         tolerance="1 H")

    # =====================================================================
    # Q-7  AN UNPRINTED QUANTITY: the DNA dose of Figures 2A, 3 and 8.
    #      Tables S3/S4/S6 do not print a dose. Table S5 does, per block, and
    #      its L2K control is 6.9+/-1.7 at 75 ng, 14.0+/-1.8 at 150 ng and
    #      22.4+/-4.0 at 200 ng -- three clearly distinct values. S3, S4 and
    #      S6 each print L2K = 14.0+/-1.8, matching the 150 ng block in BOTH
    #      the mean and the SD and matching neither other dose.
    # =====================================================================
    l2k = {}
    for v in vals:
        if v["lipidoid"] != "L2K":
            continue
        k = (v["source_table"], v["dna_ng"])
        l2k[k] = (v["value_as_printed"], v["sd_as_printed"])
    s5 = {k[1]: val for k, val in l2k.items() if k[0] == "Table S5"}
    undosed = {k[0]: val for k, val in l2k.items() if k[0] != "Table S5"}
    target = s5.get("150 ng")
    n_match = sum(1 for val in undosed.values() if val == target)
    distinct = len(set(s5.values())) == 3
    gate("Q-7_unprinted_dose_inference",
         n_match == 3 and distinct and target == ("14.0", "1.8"),
         n_match, 3,
         f"L2K at 75/150/200 ng = {s5.get('75 ng')}/{s5.get('150 ng')}/"
         f"{s5.get('200 ng')} (3 distinct); Tables S3, S4, S6 each print "
         f"{target} -- identical to the 150 ng block in mean AND SD. "
         "INFERENCE, not a printed fact: recorded as such, dose left as "
         "'not printed' on the rows and queued as DECISIONS.md D-4")

    # =====================================================================
    # Q-8  CROSS-PANEL AGREEMENT (trap 4). Tables S5 and S6 share two
    #      lipidoids. MEASURE the ratio; do not assume it is 1.
    # =====================================================================
    def get(tbl, lip, ratio=None, dose=None):
        for v in vals:
            if v["source_table"] == tbl and v["lipidoid"] == lip:
                if ratio and v["formulation_lipidoid_DOPE"] != ratio:
                    continue
                if dose and v["dna_ng"] != dose:
                    continue
                return float(v["value_as_printed"])
        return None

    # Robust to a missing cell on either panel: a None must FAIL this gate
    # cleanly, not crash it. Seed S21 removed one Table S6 row and the first
    # version of this block raised a TypeError on None/float -- which exits
    # non-zero, so the seed "passed", but the gate never reported and the
    # count never moved. A gate that crashes is not a gate that failed.
    cross = []
    for lip in ("1C16-S", "1C16-D"):
        s6v = get("Table S6", lip)
        s5c = [x for x in (get("Table S5", lip, r, "150 ng")
                           for r in ("2:1", "1:1", "1:2")) if x is not None]
        best = max(s5c) if s5c else None
        ratio = (s6v / best) if (s6v is not None and best) else None
        cross.append((lip, s6v, best, ratio))
    worst = max((c for c in cross if c[3] is not None),
                key=lambda t: t[3], default=(None, None, None, None))
    # This gate asserts that BOTH shared lipidoids were located on BOTH panels
    # and the ratio actually computed. It deliberately does NOT assert the
    # ratio is near 1 -- the measured gap is a finding, not a failure. But it
    # must not be allowed to pass by finding nothing: a report that cannot
    # fail is the ZT_2026 marginal that scored a literal 0 for twelve days.
    found = sum(1 for _, a, b, _r in cross if a is not None and b is not None)
    gate("Q-8_cross_panel_ratio", found == 2 and len(cross) == 2, found, 2,
         "; ".join(f"{l}: Table S6 {a}% vs the best Table S5 150 ng cell "
                   + (f"{b}% = {r:.2f}x" if r is not None
                      else f"{b}% = NOT COMPUTABLE (a panel cell is missing)")
                   for l, a, b, r in cross)
         + " -- REPORTED, NOT RECONCILED. Table S6 prints no DOPE ratio, so "
         + (f"{worst[0]}'s {worst[3]:.1f}x gap may be a formulation this "
            if worst[3] is not None else "the gap ")
         + "curation cannot identify. The 29 emitted rows keep each panel's "
           "own value and name its table on the row; nothing is averaged "
           "across panels. DECISIONS.md D-6")

    # =====================================================================
    # Q-12 DRAWN GEOMETRY. The only check here that is not arithmetic over a
    #      number the authors printed. See stage02_drawn.py.
    # =====================================================================
    dg = json.load(open(os.path.join(HERE, "02_drawn_geometry.json")))
    ds = dg["structures"]
    n_ok = sum(1 for r in ds if r["bond_count_agrees"]
               and r["closed_identity_holds"]
               and r["filled_segments_identity_holds"])
    gate("Q-12_drawn_bond_geometry", n_ok == 3, n_ok, 3,
         "drawn filled bond objects vs heavy bonds + double bonds "
         "(kekulised): "
         + ", ".join(f"{r['label']} {r['drawn_filled_bond_objects']}="
                     f"{r['expected_drawn_paths_value']}" for r in ds)
         + "; chains.vector_bonds(include_filled=True) returns exactly 2.0 "
           "segments per bond on all three pages (stroked mode returns 0 -- "
           "the bonds are filled outlines), and no glyph separation was "
           "needed because every glyph on these pages is a bitmap image "
           "rather than a path")

    # =====================================================================
    # Q-9  CORPUS OVERLAP (trap 2).
    # =====================================================================
    ov = json.load(open(os.path.join(HERE, "04_corpus_overlap.json")))
    nf = sum(s["n_full_key_hits"] for s in ov["our_structures"])
    nsk = sum(len(s["skeleton_only_hits"]) for s in ov["our_structures"])
    gate("Q-9_corpus_overlap", nf == 0 and nsk == 0, 3 - nf, 3,
         f"0 of 3 structures are in the corpus by FULL InChIKey and 0 by "
         f"skeleton block, over {ov['corpus_files']['unique_IL_SMILES_indexed']}"
         f" unique corpus IL_SMILES; 0 IL_name collisions. But the amine-1 "
         f"head substructure appears on "
         f"{ov['head_substructure_sweep']['n_hits']} corpus rows across 5 "
         "experiments (KZ_2016, YX_2024, YR_2024, XH_2025, YX_2026_SDL), on "
         "far more elaborate ester/amide lipids -- the head is not novel, "
         "these compounds are absent")

    # =====================================================================
    # Q-10 EMISSION RECONCILIATION.
    # =====================================================================
    em = json.load(open(os.path.join(HERE, "06_emit_summary.json")))
    ok = (em["n_emitted"] == 29 and em["n_skipped"] == 58
          and em["n_values_total"] == 87
          and em["n_emitted"] + em["n_skipped"] == em["n_values_total"])
    gate("Q-10_emission_reconciliation", ok, em["n_emitted"], 29,
         f"{em['n_emitted']} emitted + {em['n_skipped']} skipped = "
         f"{em['n_values_total']} values; confidence "
         f"{em['confidence_distribution']}; rows_with_a_caveat "
         f"{em['rows_with_a_caveat']}")

    # =====================================================================
    # Q-11 THE TWO BUILD ROUTES AGREE.
    # =====================================================================
    same = sum(1 for r in red
               if Chem.MolToSmiles(Chem.MolFromSmiles(s3[r["label"]]["IL_SMILES"]))
               == r["smiles"])
    gate("Q-11_two_build_routes_agree", same == 3, same, 3,
         "stage03's atom-mapped assemble.attach build and this file's "
         "independent string-composition re-derivation give identical "
         "canonical SMILES for all three compounds")

    # =====================================================================
    # Q-13 THE C10 ANCHOR. Does it match, and is it INDEPENDENT?
    #      Independence is tested on the GENERATOR'S COMPOSED TEMPLATE, not
    #      on the emitted anchor -- byte-identity of the emitted anchor is
    #      not the signal, because a correct re-derivation canonicalises to
    #      the same string as a correct enumeration. That is precisely why
    #      C10's own `anchor_may_be_a_copy` missed F-037.
    # =====================================================================
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "mkref", os.path.join(HERE, "make_reference.py"))
    mk = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mk)

    anchors = {r["lipid_name"]: r["reference_SMILES"] for r in csv.DictReader(
        open(os.path.join(HERE, "reference_structures.csv")))}
    n_match, n_indep, detail = 0, 0, []
    for a in mk.ANCHORS:
        lab = a["label"]
        template = mk.compose(a["chain"], a["tails"])
        enum_smiles = s3[lab]["IL_SMILES"]
        same_mol = (Chem.MolToSmiles(Chem.MolFromSmiles(anchors[lab]))
                    == Chem.MolToSmiles(Chem.MolFromSmiles(enum_smiles)))
        independent = template != enum_smiles
        n_match += same_mol
        n_indep += independent
        detail.append(f"{lab} match={'Y' if same_mol else 'N'} "
                      f"template!=enumeration={'Y' if independent else 'N'}")
    gate("Q-13_C10_anchor_match_and_independence",
         n_match == 3 and n_indep == 3, n_match + n_indep, 6,
         "; ".join(detail) + " -- the generator's HEAD-FIRST template (e.g. "
         f"'{mk.compose(12, 2)}') is byte-different from the enumeration's "
         f"canonical SMILES ('{s3['1C12-D']['IL_SMILES']}') for all three, "
         "and all three canonicalise to the same molecule")

    # ------------------------------------------------------------------
    out = dict(
        _what="Stage 07. THE GATE. Every check with numbers; exits non-zero "
              "on any failure. Each gate is seeded in 10_seed_checks.py and "
              "credited only when its agreement count moves.",
        paper="Li, Wang, Wu, Davidson, Levkin, Bioconjug Chem 2013, 24(9), "
              "1543-1551; PMID 23909648; DOI 10.1021/bc400158w",
        endpoint="average transfection efficiency, percentage of transfected "
                 "cells, plasmid DNA cargo, lipoplex with DOPE",
        value_tier="4_pdf_text_table -- zero values digitised",
        re_derived_structures=[dict(label=r["label"], smiles=r["smiles"],
                                    formula=r["formula"],
                                    neutral_exact=round(r["neutral"], 5),
                                    MH_plus=round(r["mh"], 5),
                                    formula_H=r["n_h"]) for r in red],
        gates=gates,
        n_gates=len(gates),
        n_failed=len(failures),
        failures=failures)
    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"\n{len(gates) - len(failures)}/{len(gates)} gates pass")
    if failures:
        print("FAILURES:")
        for f in failures:
            print("  -", f)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
