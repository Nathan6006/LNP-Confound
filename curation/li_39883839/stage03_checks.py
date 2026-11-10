#!/usr/bin/env python3
"""Every structural check, with numbers. Writes 04_structure_qc.json.

Each check is a METRIC, not a boolean, so a seeded defect has to MOVE it --
CLAUDE.md section 6: "A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL
ON PURPOSE." stage03_seed.py seeds one defect per check and asserts the
corresponding metric changes.

Pass `--mutate <key>` to run with one seeded defect applied.
"""
import csv, json, os, re, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
from rdkit import Chem, RDLogger                              # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors          # noqa: E402
RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

# --- SI Section A, printed MS for the nine building blocks -----------------
# "MW calc'd for <neutral formula> (M + H+): <value>, found: <value>".
# The printed "calc'd" value is in fact the NEUTRAL MONOISOTOPIC mass, not
# [M+H]+; P5 tests that reading and reports the residual either way.
BLOCK_MS = {  # block: (printed formula, printed calc, printed found)
 "A1": ("C12H23N", 181.18, 181.68), "A2": ("C13H25N", 195.20, 195.35),
 "A3": ("C17H33N", 251.26, 251.67), "A4": ("C19H35N", 277.28, 277.47),
 "B1": ("C10H20O", 156.15, 156.22), "B2": ("C18H34O", 266.26, 266.42),
 "B3": ("C18H34O3", 298.25, 298.33), "B4": ("C22H42O3", 354.31, 354.58),
 "Bbu": ("C18H34O3", 299.47, 299.27),   # printed as the cation C18H35O3+
}

# --- SI Section A, printed MS for the six individually synthesised lipids ---
LIPID_MS = {  # lipid: (printed cation formula, printed calc, printed found)
 "H2A4B3":     ("C42H81N2O5+", 694.11, 694.80),
 "H2A4B4":     ("C46H88N2O5+", 750.22, 750.55),
 "H2A4-Ole":   ("C42H81N2O3+", 662.11, 662.34),
 "H2A4-Butyl": ("C42H81N2O5+", 694.11, 694.15),
 "A4B4-S3":    ("C47H90N2O5+", 763.68, 763.49),
 "A4B4-S4":    ("C48H92N2O5+", 777.70, 777.82),
}

# --- SI Section A, the (mass, amount) pair printed for every reagent --------
# This is the check that settles which aldehyde A4B4-S3/S4 actually used, and
# it is independent of the reagent NAME printed beside it.
STOICH = {  # lipid: [(printed reagent name, mg, umol, which block we claim)]
 "H2A4B3":     [("3-(dimethylamino)propanoic acid", 30, 256.09, "H2"),
                ("(Z)-1-isocyanooctadec-9-ene", 71.06, 256.09, "A4"),
                ("6-oxohexyl dodecanoate", 76.43, 256.09, "B3")],
 "H2A4B4":     [("3-(dimethylamine) propanoic acid", 30, 256.09, "H2"),
                ("(Z)-1-isocyanooctadec-9-ene", 71.06, 256.09, "A4"),
                ("6-oxohexyl 2-hexyl decanoate", 90.8, 256.09, "B4")],
 "H2A4-Ole":   [("3-(dimethylamine) propanoic acid", 30, 256.09, "H2"),
                ("(Z)-1-isocyanooctadec-9-ene", 71.06, 256.09, "A4"),
                ("olealdehyde", 68.24, 256.09, "B2")],
 "H2A4-Butyl": [("3-(dimethylamine) propanoic acid", 30, 256.09, "H2"),
                ("(Z)-1-isocyanooctadec-9-ene", 71.06, 256.09, "A4"),
                ("6-oxohexyl 2-butyl octanoate", 76.43, 256.09, "Bbu")],
 "A4B4-S3":    [("4-(dimethylamino)butanoic acid", 30, 228.7, "S3"),
                ("(Z)-1-isocyanooctadec-9-ene", 81.09, 228.7, "A4"),
                ("6-oxohexyl dodecanoate", 63.64, 228.7, "B4")],
 "A4B4-S4":    [("5-(dimethylamino)pentanoic acid", 30, 206.6, "S4"),
                ("(Z)-1-isocyanooctadec-9-ene", 57.33, 206.6, "A4"),
                ("6-oxohexyl dodecanoate", 73.26, 206.6, "B4")],
}

MUT = (sys.argv[sys.argv.index("--mutate") + 1]
       if "--mutate" in sys.argv else None)


def mw(smi):
    return Descriptors.MolWt(Chem.MolFromSmiles(smi))


def mono(smi):
    return Descriptors.ExactMolWt(Chem.MolFromSmiles(smi))


def reagent_mw(block, blocks):
    """The MW of the block AS THE PRINTED REAGENT (isocyanide / aldehyde /
    acid), reconstructed from the precursor SMILES in 01_blocks.csv."""
    smi = blocks[block]["smiles"]
    if block.startswith("A"):          # formamide -> isocyanide: lose H2O
        return mw(smi.replace(":1", "")) - 18.015
    if block.startswith("B"):          # alcohol -> aldehyde: lose H2
        return mw(smi.replace(":2", "")) - 2.016
    return mw(smi)                      # head acids are the reagent as written


def main():
    blocks = {r["block"]: r for r in
              csv.DictReader(open(os.path.join(HERE, "01_blocks.csv")))}
    enum = {r["IL_name"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
    tab = {r["Name"]: r for r in
           csv.DictReader(open(os.path.join(HERE, "02_tableS1.csv")))}

    if MUT == "P1":        # swap two head groups in the enumeration
        enum["H1A1B1"], enum["H2A1B1"] = (dict(enum["H2A1B1"]),
                                          dict(enum["H1A1B1"]))
    if MUT == "P2":        # nudge one MW
        enum["H5A2B3"] = dict(enum["H5A2B3"]);  enum["H5A2B3"]["MW_calc"] = "999.99"
    if MUT == "P4":        # B4 -> B3 for A4B4-S3, i.e. believe the printed name
        STOICH["A4B4-S3"][2] = ("6-oxohexyl dodecanoate", 63.64, 228.7, "B3")
    if MUT == "P5":
        BLOCK_MS["B4"] = ("C22H42O3", 340.0, 354.58)

    qc = {}

    # P1  every Table S1 formula reproduced by the enumeration -------------
    bad = [(n, tab[n]["Formula"], enum[n]["Formula"].rstrip("+"))
           for n in tab if enum[n]["Formula"].rstrip("+") != tab[n]["Formula"]]
    qc["P1_tableS1_formula"] = {
        "n_compared": len(tab), "n_exact": len(tab) - len(bad),
        "n_mismatch": len(bad), "mismatches": bad[:10],
        "what": "Each of the 144 enumerated library structures vs the molecular "
                "formula Table S1 prints for that name. External to the "
                "enumeration: the formula is the authors', the structure ours."}

    # P1b  Table S1's OWN internal consistency: does each printed formula
    #      weigh what the same row's printed MW says? Needs nothing from us.
    W = {"C": 12.011, "H": 1.008, "N": 14.007, "O": 15.999}

    def fw(f):
        return sum(W[e] * (int(n) if n else 1)
                   for e, n in re.findall(r"([CHNO])(\d*)", f))
    inc = [(n, tab[n]["Formula"], float(tab[n]["MW"]), round(fw(tab[n]["Formula"]), 2))
           for n in tab if abs(fw(tab[n]["Formula"]) - float(tab[n]["MW"])) > 0.06]
    qc["P1b_tableS1_self_consistency"] = {
        "n_rows": len(tab), "n_formula_disagrees_with_own_MW": len(inc),
        "rows": inc,
        "what": "A check on the PAPER alone -- it uses no structure of ours. "
                "It separates 'our structure is wrong' from 'their table is "
                "wrong' for every P1 mismatch."}

    # P1c  for the P1 mismatches, is the printed MW our formula's weight, and
    #      is the printed set a PERMUTATION of ours over the same names?
    mm = [n for n, _, _ in bad]
    ours_ms = sorted(round(float(enum[n]["MW_calc"]), 2) for n in mm)
    theirs_ms = sorted(round(float(tab[n]["MW"]), 2) for n in mm)
    agree_mw = [n for n in mm
                if abs(float(tab[n]["MW"]) - float(enum[n]["MW_calc"])) <= 0.05]
    qc["P1c_mismatch_diagnosis"] = {
        "names": mm,
        "n_where_printed_MW_equals_OUR_formula_weight": len(agree_mw),
        "which": agree_mw,
        "our_MW_multiset": ours_ms, "printed_MW_multiset": theirs_ms,
        "multisets_identical": ours_ms == theirs_ms,
        "what": "If the two multisets are identical the printed cells are the "
                "RIGHT values on the WRONG names -- a row-alignment slip in "
                "the table, not a missing or extra composition."}

    # P2  Table S1 molecular weights ---------------------------------------
    d = [(n, float(tab[n]["MW"]), float(enum[n]["MW_calc"]),
          round(float(enum[n]["MW_calc"]) - float(tab[n]["MW"]), 3)) for n in tab]
    worst = max(d, key=lambda t: abs(t[3]))
    qc["P2_tableS1_MW"] = {
        "n_compared": len(d), "max_abs_delta_Da": abs(worst[3]),
        "worst": {"name": worst[0], "printed": worst[1], "ours": worst[2]},
        "n_within_0.05_Da": sum(1 for t in d if abs(t[3]) <= 0.05),
        "mean_abs_delta_Da": round(sum(abs(t[3]) for t in d) / len(d), 4)}

    # P3  the library is 9 x 4 x 4 and the names say so ---------------------
    keys = {(r["head"], r["tailA"], r["tailB"]) for r in enum.values()
            if r["in_library"] == "1"}
    qc["P3_grid"] = {"n_library": len(keys),
                     "n_heads": len({k[0] for k in keys}),
                     "n_tailA": len({k[1] for k in keys}),
                     "n_tailB": len({k[2] for k in keys}),
                     "product": 9 * 4 * 4,
                     "complete": len(keys) == 144}

    # P4  reagent stoichiometry: mg / umol must give the block's own MW -----
    # Two metrics, because the single "n within tolerance" count is BLIND to
    # the defect this check exists for. If two masses inside one paragraph are
    # swapped, both of them fail individually whichever block we claim -- so a
    # per-reagent tolerance count does not move when the claim changes. The
    # multiset test does: it asks whether the paragraph's three masses are the
    # three masses our three claimed blocks should have, in any order.
    st, per_lipid = [], []
    for lip, rs in STOICH.items():
        implied, claimed = [], []
        for nm, mg, umol, blk in rs:
            iw = mg / (umol * 1e-6) / 1000.0
            bw = reagent_mw(blk, blocks)
            implied.append(iw); claimed.append(bw)
            st.append({"lipid": lip, "printed_name": nm, "block": blk,
                       "implied_MW": round(iw, 2), "block_MW": round(bw, 2),
                       "rel_err_pct": round(100 * (iw - bw) / bw, 2)})
        in_order = all(abs(i - c) / c <= 0.005 for i, c in zip(implied, claimed))
        si, sc = sorted(implied), sorted(claimed)
        as_set = all(abs(i - c) / c <= 0.005 for i, c in zip(si, sc))
        perm = None
        if as_set and not in_order:
            perm = [rs[implied.index(min(implied, key=lambda x: abs(x - c)))][3]
                    for c in claimed]
        per_lipid.append({"lipid": lip, "matches_as_printed": in_order,
                          "matches_as_a_multiset": as_set,
                          "implied_MWs": [round(x, 2) for x in implied],
                          "claimed_block_MWs": [round(x, 2) for x in claimed],
                          "mass_order_vs_name_order": perm})
    worst_s = max(st, key=lambda r: abs(r["rel_err_pct"]))
    qc["P4_reagent_stoichiometry"] = {
        "n_reagents": len(st),
        "n_within_0.5_pct": sum(1 for r in st if abs(r["rel_err_pct"]) <= 0.5),
        "n_lipids": len(per_lipid),
        "n_lipids_matching_as_printed":
            sum(1 for r in per_lipid if r["matches_as_printed"]),
        "n_lipids_matching_as_a_multiset":
            sum(1 for r in per_lipid if r["matches_as_a_multiset"]),
        "max_abs_rel_err_pct": abs(worst_s["rel_err_pct"]), "worst": worst_s,
        "per_lipid": per_lipid, "rows": st,
        "what": "mg / umol is the authors' own number pair and is independent "
                "of the reagent NAME printed beside it. This is what settles "
                "the A4B4-S3/S4 aldehyde against a wrong printed name."}

    # P5  building-block MS --------------------------------------------------
    bm = []
    for blk, (pf, calc, found) in BLOCK_MS.items():
        smi = blocks[blk]["smiles"]
        if blk.startswith("A"):
            rsmi = Chem.MolToSmiles(Chem.MolFromSmiles(smi.replace(":1", "")))
            neutral_mono = mono(rsmi) - 18.0106
        elif blk.startswith("B"):
            rsmi = Chem.MolToSmiles(Chem.MolFromSmiles(smi.replace(":2", "")))
            neutral_mono = mono(rsmi) - 2.0157
        else:
            neutral_mono = mono(smi)
        bm.append({"block": blk, "printed_formula": pf,
                   "printed_calc": calc,
                   "our_neutral_monoisotopic": round(neutral_mono, 3),
                   "delta_vs_neutral_mono": round(calc - neutral_mono, 3),
                   "delta_vs_MH_plus": round(calc - (neutral_mono + 1.00728), 3)})
    qc["P5_block_MS"] = {
        "n": len(bm),
        "n_matching_neutral_monoisotopic_within_0.02":
            sum(1 for r in bm if abs(r["delta_vs_neutral_mono"]) <= 0.02),
        "n_matching_MH_plus_within_0.02":
            sum(1 for r in bm if abs(r["delta_vs_MH_plus"]) <= 0.02),
        "rows": bm,
        "what": "The SI labels these '(M + H+)' but the numbers are the "
                "NEUTRAL monoisotopic masses. Testing both readings is the "
                "point: one of them fits 8/9 and the other fits none."}

    # P6  lipid MS, under both conventions -----------------------------------
    lm = []
    for lip, (pf, calc, found) in LIPID_MS.items():
        smi = enum[lip]["IL_SMILES"]
        lm.append({"lipid": lip, "printed_cation_formula": pf,
                   "printed_calc": calc, "printed_found": found,
                   "avg_MH": round(mw(smi) + 1.008, 2),
                   "mono_MH": round(mono(smi) + 1.00728, 2),
                   "d_avg": round(calc - (mw(smi) + 1.008), 2),
                   "d_mono": round(calc - (mono(smi) + 1.00728), 2)})
    qc["P6_lipid_MS"] = {
        "n": len(lm),
        "n_fit_average_MH_within_0.05":
            sum(1 for r in lm if abs(r["d_avg"]) <= 0.05),
        "n_fit_monoisotopic_MH_within_0.05":
            sum(1 for r in lm if abs(r["d_mono"]) <= 0.05),
        "n_fit_either": sum(1 for r in lm
                            if min(abs(r["d_avg"]), abs(r["d_mono"])) <= 0.05),
        "rows": lm,
        "what": "The SI uses TWO mass conventions without saying so: average "
                "[M+H]+ for the four H2A4* lipids, monoisotopic [M+H]+ for "
                "A4B4-S3/S4. Checked against only one, four or two of six "
                "would look wrong."}

    # P7  distinct structures vs names ---------------------------------------
    byk = {}
    for n, r in enum.items():
        byk.setdefault(r["InChIKey"], []).append(n)
    dups = {k: v for k, v in byk.items() if len(v) > 1}
    qc["P7_distinct_structures"] = {
        "n_names": len(enum), "n_distinct_full_InChIKeys": len(byk),
        "duplicate_groups": {k: sorted(v) for k, v in dups.items()},
        "what": "Trap 1 of the brief. The paper's own count of what it made is "
                "144 + 4; two of those names are ONE molecule."}

    out = os.path.join(HERE, "04_structure_qc.json" if not MUT
                       else "/dev/stdout")
    txt = json.dumps(qc, indent=1)
    if MUT:
        print(txt)
    else:
        open(out, "w").write(txt)
        for k, v in qc.items():
            print(k, {kk: vv for kk, vv in v.items()
                      if not isinstance(vv, (list, dict))})
    return 0


if __name__ == "__main__":
    sys.exit(main())
