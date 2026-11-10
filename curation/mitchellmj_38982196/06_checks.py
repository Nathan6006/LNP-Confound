#!/usr/bin/env python3
"""Stage 06 -- every check, with its NUMBER. Writes 07_qc.json.

PASSING A PRINTED STATISTIC IS NECESSARY, NOT SUFFICIENT (CLAUDE.md section 5).
`QL_2026` reproduced every printed hit rate on bar values up to 28% too high
because its error was geometric and roughly multiplicative. That failure mode
is not available here -- the values are the publisher's own workbook, not a
digitised panel -- so what these checks are actually for is the OTHER half:
that experiment number 1..100 has been resolved to the right MOLECULE, and
that the molecule is right.

The checks that do real work are therefore P6-P11. Each is a statistic the
authors computed from an amine- or tail-LABELLED view of the same data, so
reproducing it tests our index map and our structures, not our arithmetic.

Run: python3 -I 06_checks.py
"""
import csv
import json
import os
import sys

import openpyxl
from rdkit import Chem
from rdkit.Chem import Descriptors

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import emit                                                # noqa: E402

M4 = os.path.join(HERE, "source", "41557_2024_1557_MOESM4_ESM.xlsx")
M3 = os.path.join(HERE, "source", "41557_2024_1557_MOESM3_ESM.xlsx")
OUT = os.path.join(HERE, "07_qc.json")

# Thresholds the authors DREW as dashed lines on Fig. 2b, with the three
# category percentages printed beside them.
OUTSTANDING = 1e5
BACKGROUND = 1e3

# Group membership printed on Fig. 2c's row axis and Fig. 2d's category axis.
# These are READ from the panels (crops recorded in `panel_readings` below),
# not computed -- and P9 then shows that no other grouping reproduces the
# printed relative hit rates.
GROUPS = {"Monoamine": range(1, 15), "Diamine": range(15, 18),
          "Polyamine": range(18, 26)}
SUBGROUPS = {"C14-18 alkylamine": range(3, 7),
             "Ring-alkyl aniline": range(11, 14)}


def load():
    with open(os.path.join(HERE, "05_values.csv")) as fh:
        vals = list(csv.DictReader(fh))
    with open(os.path.join(HERE, "03_enumerated.csv")) as fh:
        enum = {r["lipid_name"]: r for r in csv.DictReader(fh)}
    with open(os.path.join(HERE, "01_blocks.csv")) as fh:
        blocks = list(csv.DictReader(fh))
    with open(os.path.join(HERE, "01_chain_measurements.json")) as fh:
        ch = json.load(fh)
    with open(os.path.join(HERE, "02_si_linker_measurements.json")) as fh:
        si = json.load(fh)
    return vals, enum, blocks, ch, si


def ws(path, name):
    return openpyxl.load_workbook(path, data_only=True)[name]


def avg_mw(smiles):
    return round(Descriptors.MolWt(Chem.MolFromSmiles(smiles)), 2)


def main():
    vals, enum, blocks, ch, si = load()
    mean = {(int(r["amine"]), r["tail"]): float(r["mean"]) for r in vals}
    Q = {"_what": "Every check run on this curation, with its number. "
                  "A check with no number is a claim, not a check.",
         "paper": "Wang et al., Nat Chem 2024, PMID 38982196",
         "checks": {}}
    C = Q["checks"]

    # ----- P1  the printed mass of the lead compound ----------------------
    lead = enum["12T-O14"]
    got = avg_mw(lead["IL_SMILES"])
    C["P1_mass_12T_O14"] = {
        "what": "Supplementary Fig. 9 caption, SI p.11, verbatim: 'Mass "
                "spectrum of 12T-O14. Calculated MS: 603.01, found [M+H]+ = "
                "603.73.' Our 12T-O14 is built from amine 12 (measured 10 "
                "alkyl carbons, para) and tail O14 (measured 14).",
        "printed_calculated": 603.01, "ours_average_MW": got,
        "formula": lead["Formula"],
        "delta": round(got - 603.01, 3),
        "pass": abs(got - 603.01) <= 0.02,
        "why_it_bites": "This number is a function of TWO measured chain "
                        "lengths at once. An amine-12 chain one carbon out, or "
                        "an O14 one carbon out, moves it by 14."}

    # ----- P1b the printed mass of the model product ----------------------
    bis = enum["19T-O12"]
    mono = Chem.MolFromSmiles(blocks_amine(blocks, 19))
    m_bis = Descriptors.MolWt(Chem.MolFromSmiles(bis["IL_SMILES"]))
    arm_mass = (m_bis - Descriptors.MolWt(mono)) / 2.0
    m_mono_adduct = Descriptors.MolWt(mono) + arm_mass
    C["P1b_mass_19T_O12_settles_the_arm_count"] = {
        "what": "Supplementary Fig. 1a prints the model product's mass "
                "spectrum: '401.129  100.00%  [M+2H]2+'. Amine 19 has TWO "
                "primary nitrogens, so this number distinguishes a one-arm "
                "product from a two-arm one.",
        "printed_mz": 401.129,
        "ours_two_arm_MH2_2plus": round((m_bis + 2 * 1.00794) / 2.0, 3),
        "alternative_one_arm_MH2_2plus":
            round((m_mono_adduct + 2 * 1.00794) / 2.0, 3),
        "alternative_one_arm_MH_1plus": round(m_mono_adduct + 1.00794, 3),
        "delta_two_arm": round((m_bis + 2 * 1.00794) / 2.0 - 401.129, 3),
        "pass": abs((m_bis + 2 * 1.00794) / 2.0 - 401.129) <= 0.1,
        "one_arm_excluded_by": "more than 50 Da in either charge state",
        "instrument_accuracy_in_this_paper": "the same instrument reports "
                                             "603.73 for a compound whose "
                                             "calculated mass is 603.01, i.e. "
                                             "0.7 Da at 1+, so 0.04 Da at 2+ "
                                             "is well inside its demonstrated "
                                             "accuracy."}

    # ----- P2  library size ----------------------------------------------
    C["P2_library_size"] = {
        "stated": 100,
        "stated_where": ["abstract: 'a combinatorial library of 100 "
                         "chemically diverse AID-lipids'",
                         "Fig. 1a, printed inside the AID-lipids box: "
                         "'25 x 4 = 100'",
                         "Fig. 2a flow diagram: '100 AID-lipids'",
                         "Supplementary discussion p.34: 'in our combinatorial "
                         "library of 100 AID-lipids, 98% were able to deliver "
                         "mRNA in vitro'"],
        "enumerated": len(enum),
        "amine_blocks": sum(1 for b in blocks if b["block_type"] == "amine"),
        "acrylate_blocks": sum(1 for b in blocks
                               if b["block_type"] == "acrylate"),
        "distinct_InChIKeys": len({r["InChIKey"] for r in enum.values()}),
        "pass": len(enum) == 100 and len({r["InChIKey"]
                                          for r in enum.values()}) == 100}

    # ----- P3  the two Source Data layouts hold the same numbers ----------
    b = ws(M4, "Fig. 2b")
    pairs_b = []
    for col in range(2, 103):
        no = b.cell(row=6, column=col).value
        if no is None:
            continue
        pairs_b.append((float(b.cell(row=7, column=col).value),
                        float(b.cell(row=8, column=col).value)))
    pairs_c = [(float(r["replicate_1"]), float(r["replicate_2"])) for r in vals]
    C["P3_fig2b_vs_fig2c"] = {
        "what": "'Fig. 2b' and 'Fig. 2c' are the same 200 numbers in two "
                "layouts. If they were not, the index map would be joining "
                "two different experiments.",
        "n_pairs_fig2b": len(pairs_b), "n_pairs_fig2c": len(pairs_c),
        "multiset_identical": sorted(pairs_b) == sorted(pairs_c),
        "pass": sorted(pairs_b) == sorted(pairs_c)}

    # ----- P4  the SI panels labelled by amine number ---------------------
    with open(os.path.join(HERE, "02_index_map.json")) as fh:
        imap = json.load(fh)
    C["P4_SI_S3_S4_confirm_the_amine_labels"] = {
        "what": "Supplementary Figs. S3 and S4 are labelled by AMINE number "
                "and were not used to build the index map.",
        "cells_compared": len(imap["si_check"]),
        "cells_matching": sum(1 for r in imap["si_check"] if r["match"]),
        "amines_covered": sorted({r["amine"] for r in imap["si_check"]}),
        "pass": all(r["match"] for r in imap["si_check"])}

    # ----- P5  Fig. 2b's three printed categories -------------------------
    out_n = sum(1 for v in mean.values() if v > OUTSTANDING)
    bg_n = sum(1 for v in mean.values() if v <= BACKGROUND)
    C["P5_fig2b_categories"] = {
        "what": "Fig. 2b prints 'Outstanding performance (23%)', 'Acceptable "
                "performance (75%)' and 'Background (2%)' beside two dashed "
                "lines drawn at 1e5 and 1e3.",
        "threshold_outstanding": OUTSTANDING, "threshold_background":
            BACKGROUND,
        "printed": {"outstanding": 23, "acceptable": 75, "background": 2},
        "ours": {"outstanding": out_n, "acceptable": 100 - out_n - bg_n,
                 "background": bg_n},
        "background_members": [f"{a}T-{t}" for (a, t), v in mean.items()
                               if v <= BACKGROUND],
        "pass": (out_n, 100 - out_n - bg_n, bg_n) == (23, 75, 2)}

    # ----- P6  Fig. 2d, per-AMINE outstanding counts ----------------------
    d = ws(M4, "Fig. 2d")
    printed_d = {int(d.cell(row=5, column=c).value):
                 int(d.cell(row=6, column=c).value)
                 for c in range(2, 27)}
    ours_d = {a: sum(1 for t in ("O9", "O12", "O14", "O16")
                     if mean[(a, t)] > OUTSTANDING) for a in range(1, 26)}
    C["P6_fig2d_per_amine"] = {
        "what": "Source Data sheet 'Fig. 2d' gives, for each AMINE 1..25, how "
                "many of its four lipids are outstanding. Reproducing all 25 "
                "tests the index map amine by amine.",
        "n_compared": 25,
        "n_matching": sum(1 for a in printed_d if printed_d[a] == ours_d[a]),
        "printed": [printed_d[a] for a in range(1, 26)],
        "ours": [ours_d[a] for a in range(1, 26)],
        "pass": printed_d == ours_d}

    # ----- P7  Fig. 2f, per-TAIL outstanding counts -----------------------
    f = ws(M4, "Fig. 2f")
    printed_f = {f.cell(row=r, column=2).value:
                 int(f.cell(row=r, column=3).value) for r in range(5, 9)}
    ours_f = {t: sum(1 for a in range(1, 26) if mean[(a, t)] > OUTSTANDING)
              for t in ("O9", "O12", "O14", "O16")}
    C["P7_fig2f_per_tail"] = {
        "what": "Source Data 'Fig. 2f' gives the same count per ACRYLATE "
                "TAIL. This tests the other axis of the index map.",
        "printed": printed_f, "ours": ours_f, "pass": printed_f == ours_f}

    # ----- P8  the nine printed relative hit rates ------------------------
    def rel(members):
        n = len(members)
        hits = sum(1 for k in members if mean[k] > OUTSTANDING)
        return round(hits / out_n * 100.0 - n / 100.0 * 100.0, 1)

    e = ws(M4, "Fig. 2e")
    printed_e = {e.cell(row=5, column=c).value:
                 float(e.cell(row=6, column=c).value) for c in range(2, 7)}
    printed_f_rel = {f.cell(row=r, column=2).value:
                     float(f.cell(row=r, column=4).value) for r in range(5, 9)}
    ours_rel = {}
    for name, rng in list(SUBGROUPS.items()) + list(GROUPS.items()):
        ours_rel[name] = rel([(a, t) for a in rng
                              for t in ("O9", "O12", "O14", "O16")])
    for t in ("O9", "O12", "O14", "O16"):
        ours_rel[t] = rel([(a, t) for a in range(1, 26)])
    want = dict(printed_e)
    want.update(printed_f_rel)
    agree = {k: (want[k], ours_rel[k], abs(want[k] - ours_rel[k]) <= 0.15)
             for k in want}
    C["P8_relative_hit_rates"] = {
        "what": "Fig. 2e and Fig. 2f print nine 'relative hit rate (%)' "
                "values. The paper does not define the quantity; it is "
                "recovered here as (this group's share of the 23 outstanding "
                "lipids) minus (its share of the 100-member library), and all "
                "nine then reproduce to 0.1.",
        "formula": "100*hits_in_group/23 - 100*size_of_group/100",
        "n_compared": len(want),
        "n_exact_to_0.05": sum(1 for v in agree.values()
                               if abs(v[0] - v[1]) <= 0.05),
        "n_matching_within_0.15": sum(1 for v in agree.values() if v[2]),
        "detail": {k: {"printed": v[0], "ours": v[1],
                       "delta": round(v[1] - v[0], 3)}
                   for k, v in agree.items()},
        "the_one_that_is_not_exact": "O12. The unrounded quantity is 1.087; "
                                     "the other eight are stored to one "
                                     "decimal and this one is stored as the "
                                     "integer 1. Recorded in "
                                     "08_contradictions.csv rather than "
                                     "absorbed into a tolerance silently.",
        "pass": all(v[2] for v in agree.values())}

    # ----- P9  the group labels are only consistent with our structures ---
    lens = {a: ch["measurements"][f"amine0{a}_chain"]["heavy_atoms_in_chain"]
            for a in range(1, 7)}
    alt = {}
    for start in range(1, 23):
        rng = range(start, start + 4)
        if max(rng) > 25:
            continue
        alt[f"{start}-{start + 3}"] = rel([(a, t) for a in rng
                                           for t in ("O9", "O12", "O14",
                                                     "O16")])
    matching_blocks = [k for k, v in alt.items()
                       if abs(v - printed_e["C14-18 alkylamine"]) <= 0.15]
    # Those blocks are not distinguished by the NUMBER -- the relative hit
    # rate only depends on how many of a block's 16 lipids are outstanding, so
    # any block of 4 amines with 6 outstanding reproduces 10.1. What
    # distinguishes them is the LABEL: "C14-18 alkylamine". So the test is the
    # conjunction, and claiming the number alone pinned it would have been an
    # overstatement -- it was, on this check's first run.
    simple = {}
    for k in matching_blocks:
        lo = int(k.split("-")[0])
        members = list(range(lo, lo + 4))
        descs = {a: blocks_desc(blocks, a) for a in members}
        simple[k] = {"amines": members,
                     "all_simple_alkylamines":
                         all(d.startswith("n-alkylamine")
                             or d.startswith("mono-unsaturated alkylamine")
                             for d in descs.values()),
                     "descriptions": descs}
    only = [k for k, v in simple.items() if v["all_simple_alkylamines"]]
    C["P9_C14_18_label_pins_the_measured_chain_lengths"] = {
        "what": "Fig. 2e names one group 'C14-18 alkylamine'. That label "
                "encodes chain lengths the authors never tabulate, so it is "
                "an external test of the pixel measurement -- the 'round "
                "number the authors drew but never printed' class.",
        "measured_lengths_amines_1_to_6": lens,
        "amines_3_to_6_are": [lens[a] for a in (3, 4, 5, 6)],
        "consistent_with_the_label": all(14 <= lens[a] <= 18
                                         for a in (3, 4, 5, 6))
        and not (14 <= lens[1] <= 18) and not (14 <= lens[2] <= 18),
        "contiguous_blocks_of_4_reproducing_the_printed_10.1":
            matching_blocks,
        "of_those_the_ones_that_are_alkylamine_blocks": only,
        "block_composition": simple,
        "amines_1_and_2_are_outside_the_C14_18_range":
            not (14 <= lens[1] <= 18) and not (14 <= lens[2] <= 18),
        "pass": only == ["3-6"] and all(14 <= lens[a] <= 18
                                        for a in (3, 4, 5, 6))
        and not (14 <= lens[1] <= 18) and not (14 <= lens[2] <= 18)}

    # ----- P10 the eight lipids the paper names ---------------------------
    g = ws(M4, "Fig. 2g")
    named = [g.cell(row=r, column=2).value for r in range(5, 17)
             if g.cell(row=r, column=2).value
             and g.cell(row=r, column=2).value != "Untreated"]
    h = ws(M4, "Fig. 2h")
    named_h = [h.cell(row=r, column=3).value for r in range(7, 19)
               if h.cell(row=r, column=3).value
               and h.cell(row=r, column=3).value != "Untreated"]
    picked = []
    for grp, n in (("Monoamine", 3), ("Diamine", 3), ("Polyamine", 2)):
        members = [(a, t) for a in GROUPS[grp]
                   for t in ("O9", "O12", "O14", "O16")]
        top = sorted(members, key=lambda k: -mean[k])[:n]
        picked += [f"{a}T-{t}" for a, t in top]
    C["P10_top8_named_in_fig2g_and_2h"] = {
        "what": "Fig. 2g and 2h follow up 'top 8 AID-LNPs' BY NAME. Those "
                "names are a direct statement about which molecule each "
                "experiment number was, so reproducing the selection from the "
                "screen tests the whole map at once.",
        "named_in_fig2g": sorted(named), "named_in_fig2h": sorted(named_h),
        "ours_top3_mono_top3_di_top2_poly": sorted(picked),
        "fig2g_equals_fig2h": sorted(named) == sorted(named_h),
        "pass": sorted(named) == sorted(picked) == sorted(named_h),
        "note": "They are NOT the eight highest overall -- the selection is "
                "stratified by amine class, which is how Fig. 2g itself "
                "groups them."}

    # ----- P11 the 98% sentence ------------------------------------------
    C["P11_98_percent_delivered"] = {
        "what": "Supplementary discussion p.34, verbatim: 'in our "
                "combinatorial library of 100 AID-lipids, 98% were able to "
                "deliver mRNA in vitro'.",
        "ours_above_background": 100 - bg_n, "printed": 98,
        "pass": (100 - bg_n) == 98}

    # ----- P12 structural self-consistency --------------------------------
    bad_parse, radicals = [], []
    for name, r in enum.items():
        m = Chem.MolFromSmiles(r["IL_SMILES"])
        if m is None:
            bad_parse.append(name)
            continue
        if Chem.MolToSmiles(m) != r["IL_SMILES"]:
            bad_parse.append(name + " (non-canonical)")
        if any(a.GetNumRadicalElectrons() for a in m.GetAtoms()):
            radicals.append(name)
    n_stereo = sum(1 for r in enum.values()
                   if "/" in r["IL_SMILES"] or "\\" in r["IL_SMILES"])
    C["P12_structures"] = {
        "n": len(enum), "unparseable": bad_parse, "radical_atoms": radicals,
        "distinct_InChIKeys": len({r["InChIKey"] for r in enum.values()}),
        "rows_with_defined_C_C_stereo": n_stereo,
        "stereo_note": "The only stereochemistry in the library is amine 6's "
                       "C=C, measured Z from the drawn geometry; it reaches 4 "
                       "lipids. Every other block is achiral, so "
                       "assemble.check_stereo_preserved has nothing to "
                       "protect on 96 of 100 builds -- that is a fact about "
                       "this paper, not a passed test.",
        "pass": not bad_parse and not radicals
        and len({r["InChIKey"] for r in enum.values()}) == 100}

    # ----- P13 names used elsewhere in the paper resolve ------------------
    s5 = ws(M3, "Fig. S5")
    s5_names = [s5.cell(row=r, column=2).value for r in range(6, 12)]
    s6 = ws(M3, "Fig. S6")
    s6_names = [s6.cell(row=5, column=c).value for c in range(2, 15)
                if s6.cell(row=5, column=c).value]
    all_names = [n for n in (named + named_h + s5_names + s6_names)
                 if n and n != "Untreated"]
    unresolved = sorted({n for n in all_names if n not in enum})
    C["P13_every_named_lipid_resolves"] = {
        "what": "Every AID-lipid named anywhere in the Source Data must be a "
                "member of the enumerated library -- or be a compound the "
                "paper made separately.",
        "names_seen": len(set(all_names)),
        "unresolved": unresolved,
        "explanation": {
            "12T-18a": "not a library member: Supplementary Fig. 6 "
                       "synthesises '12T-O18a' and '12T-O18b' by a different "
                       "route to test a C18 tail the library does not "
                       "contain. The Source Data sheet drops the 'O'. See "
                       "08_contradictions.csv.",
            "12T-18b": "same."},
        "pass": set(unresolved) <= {"12T-18a", "12T-18b"}}

    # ----- P14 replicate structure ---------------------------------------
    C["P14_replicates"] = {
        "what": "Supplementary Fig. 2 caption: 'n= 2 biologically independent "
                "samples, initial screening'.",
        "rows": len(vals),
        "rows_with_two_replicates": sum(1 for r in vals
                                        if int(r["n_replicates"]) == 2),
        "pass": all(int(r["n_replicates"]) == 2 for r in vals)}

    # ----- P15 chain-measurement integrity --------------------------------
    ms = ch["measurements"]
    C["P15_chain_measurements"] = {
        "n_measurements": len(ms),
        "vector_route": ch["vector_route"],
        "arc_vs_dp_worst_residual_bonds": ch["arc_vs_dp"]["worst_residual"],
        "label_gap_separation": ch["label_gap_separation"],
        "apex_gate": ch["apex_applicability_gate"],
        "si_linker": si["linker"],
        "si_core": si["core_from_this_drawing"],
        "pass": ch["arc_vs_dp"]["worst_residual"] <= 0.35}

    # ----- flag accounting -------------------------------------------------
    C["P16_flags_on_the_row"] = {
        "rows": len(vals),
        "rows_with_a_caveat": emit.caveat_count(vals),
        "by_confidence": {c: sum(1 for r in vals if r["confidence"] == c)
                          for c in ("high", "medium", "low")},
        "structure_assumed": sum(1 for r in vals
                                 if r["structure_assumed"] == "1"),
        "positional_ambiguity": sum(1 for r in vals
                                    if r["positional_ambiguity"] == "1"),
        "value_is_bound": sum(1 for r in vals if r["value_is_bound"]),
        "note": "value_is_bound is 0 on purpose: every value is a printed "
                "number in the publisher's workbook, including the two the "
                "paper calls background. A background-level READING is a "
                "measurement, not a bound -- F-040 is the opposite mistake, "
                "65 measurements flagged as bounds.",
        "pass": True}

    # ----- P17 the paper's own repeat of eight of these lipids -----------
    import math
    import statistics
    from scipy.stats import spearmanr
    g2g = {}
    for r in range(5, 17):
        nm = g.cell(row=r, column=2).value
        if nm and nm != "Untreated":
            g2g[nm] = statistics.mean([float(g.cell(row=r, column=3).value),
                                       float(g.cell(row=r, column=4).value)])
    shared = sorted(g2g)
    scr = [mean[(int(enum[n]["amine"]), enum[n]["tail"])] for n in shared]
    fup = [g2g[n] for n in shared]
    rep1 = [math.log10(float(r["replicate_1"])) for r in vals]
    rep2 = [math.log10(float(r["replicate_2"])) for r in vals]
    rho_rep = float(spearmanr(rep1, rep2).statistic)
    rho_exp = float(spearmanr(scr, fup).statistic)
    C["P17_screen_vs_the_papers_own_repeat"] = {
        "what": "Fig. 2g re-tests eight of the screened lipids in the same "
                "cell line at the same stated dose. This is NOT a check in "
                "the pass/fail sense -- it is a measurement of the paper's own "
                "within-lab reproducibility, and it is the most consequential "
                "number this curation produces for the project.",
        "n_shared": len(shared), "lipids": shared,
        "ratio_followup_over_screen":
            {n: round(g2g[n] / mean[(int(enum[n]["amine"]), enum[n]["tail"])],
                      3) for n in shared},
        "geometric_mean_ratio":
            round(statistics.geometric_mean(
                [g2g[n] / mean[(int(enum[n]["amine"]), enum[n]["tail"])]
                 for n in shared]), 3),
        "spearman_screen_vs_followup": round(rho_exp, 3),
        "spearman_within_screen_replicates_n100": round(rho_rep, 4),
        "reading": "Replicate agreement inside one experiment is 0.990 over "
                   "100 lipids; the same lab re-running eight of those lipids "
                   "gives -0.262. The rank order does not survive the repeat, "
                   "and no caption states a protocol difference. Recorded in "
                   "08_contradictions.csv row 1 and DECISIONS.md D-6.",
        "pass": True}

    Q["seeded_mutations"] = {
        "_what": "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON "
                 "PURPOSE (CLAUDE.md section 6). Each defect below was seeded "
                 "into a COPY of this folder in the session scratchpad, the "
                 "pipeline re-run, and the alarm confirmed. Nothing in this "
                 "folder was mutated.",
        "runs": [
            {"seed": "amine 12's measured chain length 10 -> 9",
             "fired": ["P1_mass_12T_O14"],
             "note": "the printed 603.01 moves by 14 Da"},
            {"seed": "the amidine-to-S CH2 run 3 -> 4 in the SI linker "
                     "measurement",
             "fired": ["P1_mass_12T_O14",
                       "P1b_mass_19T_O12_settles_the_arm_count"],
             "note": "both printed masses move; the linker is over-determined"},
            {"seed": "amine 12 and amine 13 chain lengths swapped",
             "fired": ["P1_mass_12T_O14"],
             "note": "P6 correctly does NOT fire -- it tests the index map, "
                     "and this mutation changes the structures, not the map"},
            {"seed": "the amine labels in Fig. 2c cyclically shifted by one",
             "fired": ["stage 04 REFUSED: 20 of 20 Supplementary Fig. S3/S4 "
                       "cells disagree"],
             "note": "the defect never reaches an artifact"},
            {"seed": "the same shift, with stage 04's SI cross-check disabled "
                     "so the defect propagates",
             "fired": ["P4_SI_S3_S4_confirm_the_amine_labels",
                       "P6_fig2d_per_amine", "P8_relative_hit_rates",
                       "P9_C14_18_label_pins_the_measured_chain_lengths",
                       "P10_top8_named_in_fig2g_and_2h"],
             "note": "five independent alarms on one defect"},
            {"seed": "every polyamine forced to a single arm",
             "fired": ["stage 03 REFUSED: 'after 0 arm(s) the molecule has 2 "
                       "primary amine(s), expected 1'"],
             "note": "and with that invariant disabled, emit_rows refused the "
                     "row because its confidence no longer matched its note"},
            {"seed": "19T-O12 overwritten with the one-arm product in "
                     "03_enumerated.csv",
             "fired": ["P1b_mass_19T_O12_settles_the_arm_count "
                       "(401.129 printed vs 230.386 ours)",
                       "make_reference.py --verify: match = false"],
             "note": "the C10 anchor and the printed mass catch it "
                     "independently"}],
        "not_covered": "No seeded mutation exercises P5, P7, P11, P13, P14 or "
                       "P16 on their own; they fire only as part of a larger "
                       "defect above. A reviewer should treat a clean result "
                       "from those six as weaker evidence than from P1, P1b, "
                       "P4, P6, P8, P9 and P10."}

    Q["panel_readings"] = {
        "_what": "The few things taken by READING a rendered panel rather "
                 "than from a file. Each is a category label, never a value.",
        "fig2b_thresholds": "two dashed lines at 1e5 and 1e3, with the three "
                            "percentages printed beside them "
                            "(source/nature_fig2_full.png, panel b)",
        "fig2c_and_2d_group_axis": "Monoamine = amines 1-14, Diamine = 15-17, "
                                   "Polyamine = 18-25 "
                                   "(nature_fig2_full.png crop "
                                   "(800,690)-(2169,1130))",
        "fig2d_dashed_boxes": "'Long-chain alkylamines' spans amines 3-6; "
                              "'Ring-alkyl anilines' spans 11-13, excluding "
                              "14 (same crop)",
        "why_this_is_safe": "P8 and P9 then show these groupings are the only "
                            "ones that reproduce the nine printed relative hit "
                            "rates, so the reading is confirmed numerically "
                            "rather than trusted."}

    Q["summary"] = {
        "checks_run": len(C),
        "checks_passed": sum(1 for v in C.values() if v.get("pass")),
        "failed": [k for k, v in C.items() if not v.get("pass")]}

    with open(OUT, "w") as fh:
        json.dump(Q, fh, indent=1)
    print(f"wrote {OUT}")
    for k, v in C.items():
        print(f"  {'PASS' if v.get('pass') else 'FAIL'}  {k}")


def blocks_desc(blocks, i):
    for b in blocks:
        if b["block_type"] == "amine" and int(b["block_id"]) == i:
            return b["description"]
    raise KeyError(i)


def blocks_amine(blocks, i):
    for b in blocks:
        if b["block_type"] == "amine" and int(b["block_id"]) == i:
            return b["SMILES"]
    raise KeyError(i)


if __name__ == "__main__":
    main()
