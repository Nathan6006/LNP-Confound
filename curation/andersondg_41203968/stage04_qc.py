#!/usr/bin/env python3
"""Every check, with its number. Writes 07_qc.json and 08_contradictions.csv.

    python3 -I stage04_qc.py

Passing every printed statistic is NECESSARY, NOT SUFFICIENT (CLAUDE.md §5):
QL_2026 reproduced 48 printed hit rates on bars up to 28% too high. So the
checks below are of three different kinds and the QC file says which is which:

  EXTERNAL-MASS   the SI's own HRMS m/z against the reconstructed structure.
                  Blind to connectivity; sees the hydrogen count, which is
                  what separates a monocyclic headgroup from a spiro one.
  INTERNAL-GRAPH  the same block re-derived from up to 27 differently-drawn
                  molecules, and the blocks re-assembled into each lipid.
                  Blind to a systematic misreading of a drawing convention;
                  sees any one-off misread.
  PRINTED-TEXT    numbers the authors wrote in the article, recomputed from
                  the publisher's Source Data or from the structures.
"""
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from openpyxl import load_workbook                         # noqa: E402
from rdkit import Chem                                     # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt              # noqa: E402

SRC = os.path.join(HERE, "source")
PROTON = 1.00727646
QC = {"checks": [], "counts": {}}
CONTRA = []


def chk(cid, kind, what, result, **kw):
    row = {"id": cid, "kind": kind, "what": what, "result": result}
    row.update(kw)
    QC["checks"].append(row)
    return row


def load(name):
    return list(csv.DictReader(open(os.path.join(HERE, name))))


def xl(n, sh):
    wb = load_workbook(os.path.join(SRC, f"41565_2025_2044_MOESM{n}_ESM.xlsx"),
                       read_only=True, data_only=True)
    rows = [list(r) for r in wb[sh].iter_rows(values_only=True)]
    wb.close()
    return rows


def mean(v):
    v = [float(x) for x in v if isinstance(x, (int, float))]
    return sum(v)/len(v) if v else float("nan")


def main():
    structures = load("02_structures.csv")
    enum = load("03_enumerated.csv")
    blocks = load("01_blocks.csv")
    values = load("05_values.csv")
    cons = json.load(open(os.path.join(HERE, "02_blocks_consensus.json")))

    # ---- P1 library size ------------------------------------------------
    grid = [f"AMG{a}{g}" for a in range(1, 11) for g in range(1, 11)]
    built = {r["label"] for r in enum if r["SMILES"]}
    chk("P1", "PRINTED-TEXT",
        "Library size: the article states 'a combinatorial library of 100 "
        "unique ionizable lipids (AMG) by combining 10 amine headgroups (Am), "
        "9 glycidyl esters and 1 epoxy-ester tails (G)' and 'we then "
        "synthesized and purified (>90% purity) all 100'.",
        "PASS" if len([x for x in grid if x in built]) == 100 else "FAIL",
        stated=100, enumerated=len([x for x in grid if x in built]),
        drawn_in_SI=len([r for r in enum
                         if r["label"] in grid and r["drawn_in_SI"] == "True"]),
        enumerated_from_blocks_only=sorted(
            r["label"] for r in enum
            if r["label"] in grid and r["provenance"]
            == "enumerated_from_blocks"))

    # ---- P2 HRMS gate ---------------------------------------------------
    st = {}
    for r in structures:
        st.setdefault(r["status"], []).append(r["label"])
    chk("P2", "EXTERNAL-MASS",
        "Each reconstructed structure against the HRMS m/z printed beside it "
        "in the SI, as [M+H]+ or [M+2H]2+, gate 25 ppm.",
        "REPORTED", n_pass=len(st.get("PASS", [])),
        n_fail=len(st.get("FAIL", [])),
        n_without_printed_HRMS=len(st.get("NO_HRMS", [])),
        median_abs_ppm=round(sorted(abs(float(r["ppm"])) for r in structures
                                    if r["status"] == "PASS")
                             [len(st.get("PASS", []))//2], 2),
        failures=sorted(st.get("FAIL", [])))

    # ---- P3 block consensus --------------------------------------------
    reads = {r["block_id"]: int(r["n_independent_reads"]) for r in blocks}
    chk("P3", "INTERNAL-GRAPH",
        "Each building block is re-derived independently from every lipid "
        "that contains it; all reads of one block must give the same graph.",
        "PASS" if not cons["block_disagreements"] else "FAIL",
        n_blocks=len(blocks), n_disagreeing=len(cons["block_disagreements"]),
        max_independent_reads=max(reads.values()),
        median_independent_reads=sorted(reads.values())[len(reads)//2])

    # ---- P4 re-enumeration ----------------------------------------------
    m = [r for r in enum if r["matches_drawing"] == "True"]
    x = [r for r in enum if r["matches_drawing"] == "False"]
    chk("P4", "INTERNAL-GRAPH",
        "Blocks re-assembled with curation/lib/assemble must reproduce the "
        "drawn molecule exactly (canonical SMILES identity).",
        "PASS" if not x else "FAIL", n_match=len(m), n_differ=len(x),
        differ=[r["label"] for r in x],
        stereo_waivers=cons["n_stereo_waivers"],
        stereo_waiver_reason="assemble.check_stereo_preserved refuses a bond "
        "that CREATES a stereocentre; waived only where neither fragment nor "
        "product carries a single ASSIGNED centre, so nothing can be inverted")

    # ---- P5 structural statements in the article ------------------------
    bl = {r["block_id"]: r for r in blocks}
    hexacyclen = Chem.MolToSmiles(Chem.MolFromSmiles("C1CNCCNCCNCCNCCNCCN1"))
    am15 = Chem.MolToSmiles(Chem.MolFromSmiles(bl["Am15"]["free_amine_SMILES"]))
    chk("P5a", "PRINTED-TEXT",
        "'AMG1541 ... features a hexacyclen headgroup (Am15)' (Conclusions). "
        "Am15 is derived here only from AMG1541's and AMG15L-41's drawings.",
        "PASS" if am15 == hexacyclen else "FAIL",
        Am15_as_derived=am15, hexacyclen_written_from_its_name=hexacyclen)

    g41 = Chem.MolFromSmiles(bl["G41"]["SMILES_with_attachment_map"])
    acyl_c = sum(1 for a in g41.GetAtoms() if a.GetSymbol() == "C") - 3
    dbl = [b for b in g41.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE
           and b.GetBeginAtom().GetSymbol() == "C"
           and b.GetEndAtom().GetSymbol() == "C"]
    chk("P5b", "PRINTED-TEXT",
        "'... and monounsaturated 12 carbon tails' (Conclusions). G41's acyl "
        "chain, counted from the graph, not from the name token '41'.",
        "PASS" if acyl_c == 12 and len(dbl) == 1 else "FAIL",
        acyl_carbons=acyl_c, n_CC_double_bonds=len(dbl))

    tails = [r for r in blocks if r["kind"] != "amine_headgroup"]
    allcc, cis = 0, 0
    for r in tails:
        mol = Chem.MolFromSmiles(r["SMILES_with_attachment_map"])
        for b in mol.GetBonds():
            if (b.GetBondType() == Chem.BondType.DOUBLE
                    and b.GetBeginAtom().GetSymbol() == "C"
                    and b.GetEndAtom().GetSymbol() == "C"):
                allcc += 1
                cis += str(b.GetStereo()) in ("STEREOZ", "STEREOCIS")
    chk("P5c", "PRINTED-TEXT",
        "'Within the tails we included Z-alkenes ...' (Results). Geometry read "
        "from the drawing's own coordinates, not assumed.",
        "PASS" if allcc and cis == allcc else "FAIL",
        n_CC_double_bonds_in_tail_blocks=allcc, n_Z=cis, n_E=allcc-cis)

    # 9 glycidyl esters + 1 epoxy ester: a glycidyl ester opens at the CH2 of
    # the glycidyl group, so the arm starts CH2-CH(OH)-CH2-O-C(=O).
    glyc = Chem.MolFromSmarts("[CH2:2][CH](O)[CH2]OC=O")
    n_gly = n_epoxy = 0
    for i in range(1, 11):
        arm = Chem.MolFromSmiles(bl[f"G{i}"]["SMILES_with_attachment_map"])
        hit = arm.GetSubstructMatches(
            Chem.MolFromSmarts("[CH3,CH2][CH](O)[CH2]OC=O"))
        att = [a.GetIdx() for a in arm.GetAtoms() if a.GetAtomMapNum() == 2]
        if any(att[0] in h for h in hit):
            n_gly += 1
        else:
            n_epoxy += 1
    chk("P5d", "PRINTED-TEXT",
        "'... 9 glycidyl esters and 1 epoxy-ester tails (G)'. A glycidyl ester "
        "opens at the glycidyl CH2, so its arm reads CH2-CH(OH)-CH2-O-C(=O); "
        "the epoxy ester does not.",
        "PASS" if (n_gly, n_epoxy) == (9, 1) else "FAIL",
        n_glycidyl=n_gly, n_epoxy_ester=n_epoxy,
        the_epoxy_ester=[f"G{i}" for i in range(1, 11)
                         if not any(
                             [a.GetIdx() for a in Chem.MolFromSmiles(
                                 bl[f'G{i}']['SMILES_with_attachment_map']
                             ).GetAtoms() if a.GetAtomMapNum() == 2][0] in h
                             for h in Chem.MolFromSmiles(
                                 bl[f'G{i}']['SMILES_with_attachment_map']
                             ).GetSubstructMatches(Chem.MolFromSmarts(
                                 "[CH3,CH2][CH](O)[CH2]OC=O")))])

    # ---- P6 printed fold-changes, from the Source Data ------------------
    f2 = xl(4, "Figure 2")
    f5 = xl(7, "Figure 5")
    e5 = xl(12, "Extended Figure 5")

    liver = [mean(r[1:2]) for r in f2[36:39]], [mean(r[2:3]) for r in f2[36:39]]
    a_liv = mean([f2[i][1] for i in (36, 37, 38)])
    b_liv = mean([f2[i][2] for i in (36, 37, 38)])
    chk("P6a", "PRINTED-TEXT",
        "'AMG1541 LNPs exhibit 8.3-fold lower expression in the liver "
        "(Figure 2f)'. Figure 2f's two Source Data columns are UNLABELLED; "
        "the ratio identifies them -- column 1 / column 2 = the printed value, "
        "so column 1 is SM-102 and column 2 is AMG1541.",
        "PASS" if abs(a_liv/b_liv - 8.3) < 0.1 else "FAIL",
        recomputed=round(a_liv/b_liv, 3), printed=8.3,
        cells="MOESM4 'Figure 2' B38:C40")

    amg = [f2[i][4] for i in range(14, 21)]
    sm = [f2[i][1] for i in range(14, 17)]
    chk("P6b", "PRINTED-TEXT",
        "'... 3.5-fold higher expression at the injection site relative to "
        "SM-102 (Figure 2c)'.",
        "PASS" if abs(mean(amg)/mean(sm) - 3.5) < 0.1 else "FAIL",
        recomputed=round(mean(amg)/mean(sm), 3), printed=3.5,
        cells="MOESM4 'Figure 2' B15:B17 and E15:E21")

    f5c = {str(r[1]): r[2] for r in f5[23:37] if r[1] is not None}
    chk("P6c", "PRINTED-TEXT",
        "'Systematic optimization ... identified AMG1541 which gave 4.5-fold "
        "higher mRNA expression than our initial hit, AMG104'. Figure 5c's own "
        "'Fold Luciferase (AMG 104)' column.",
        "PASS" if abs(float(f5c["1541"]) - 4.5) < 0.05 else "FAIL",
        source_data_value=float(f5c["1541"]), printed=4.5,
        cells="MOESM7 'Figure 5' C31")

    ohd = mean(f5[48][1:4]), mean(f5[48][5:8])
    chk("P6d", "PRINTED-TEXT",
        "'LNPs formulated with AMG1541-OHD exhibited ~100-fold loss in IM "
        "FLuc expression (Figure 5f)'.",
        "PASS" if 50 <= ohd[0]/ohd[1] <= 200 else "FAIL",
        recomputed=round(ohd[0]/ohd[1], 1), printed="~100",
        cells="MOESM7 'Figure 5' B49:D49 and F49:H49")

    e5d = mean([e5[i][1] for i in (20, 21, 22)]), \
        mean([e5[i][2] for i in (20, 21, 22)])
    chk("P6e", "PRINTED-TEXT",
        "'AMG1559 ... showed a dramatic 34-fold drop in transfection "
        "(Extended Figure 5d)'.",
        "PASS" if abs(e5d[0]/e5d[1] - 34) < 1 else "FAIL",
        recomputed=round(e5d[0]/e5d[1], 2), printed=34,
        cells="MOESM12 'Extended Figure 5' B21:C23")

    # ---- P7 legend n against the Source Data ----------------------------
    legend_n = {
        "Figure 2a": {"AMG1042": 2, "AMG104": 3, "AMG1041": 3, "AMG1043": 3,
                      "AMG1044": 3, "AMG10C12": 3},
        "Figure 2b": {"AMG141": 2, "AMG1241": 2, "AMG1641": 2, "AMG1141": 2,
                      "AMG1441": 3, "AMG1841": 3, "AMG1941": 3, "AMG1341": 4,
                      "AMG1041": 5, "AMG1541": 5, "AMG1741": 5},
        "Extended Figure 2b": {"AMG1549": 3, "AMG1550": 3, "AMG1551": 3,
                               "AMG1553": 3, "AMG1555": 3, "AMG1541": 5,
                               "AMG1547": 5, "AMG1548": 5, "AMG1554": 5,
                               "AMG1556": 5}}
    pmap = {"2a": "Figure 2a", "2b": "Figure 2b",
            "Extended Figure 2b": "Extended Figure 2b"}
    agree = disagree = 0
    bad = []
    for v in values:
        panel = pmap.get(v["panel"])
        if not panel:
            continue
        want = legend_n[panel].get(v["lipid_label"])
        if want is None:
            continue
        if int(v["n_values"]) == want:
            agree += 1
        else:
            disagree += 1
            bad.append((panel, v["lipid_label"], int(v["n_values"]), want))
    chk("P7", "PRINTED-TEXT",
        "The per-lipid n printed in each figure legend against the number of "
        "values the publisher's Source Data actually contains for that lipid. "
        "27 statements, each naming a specific lipid.",
        "PASS" if disagree == 0 else "FAIL", n_agree=agree,
        n_disagree=disagree, disagreements=bad)

    # ---- P8 LC-MS purity ------------------------------------------------
    import pymupdf
    doc = pymupdf.open(os.path.join(SRC, "41565_2025_2044_MOESM1_ESM.pdf"))
    pur = []
    for p in range(44, 122):
        pur += [float(x) for x in
                re.findall(r"LC-MS purity:\s*([0-9.]+)\s*%", doc[p-1].get_text())]
    chk("P8", "PRINTED-TEXT",
        "'The purity of all AMG lipids was >=90% as judged by LC-MS' (SI p.5) "
        "against every 'LC-MS purity' value printed in the SI.",
        "PASS" if pur and min(pur) >= 90 else "FAIL",
        n_values=len(pur), min=min(pur) if pur else None,
        median=sorted(pur)[len(pur)//2] if pur else None)

    # ---- P9 the in vitro screen -----------------------------------------
    ev = [v for v in values if v["panel"] == "Extended Figure 1"]
    strong = {t: sum(1 for v in ev if float(v["mean"]) >= t)
              for t in (1e4, 2e4, 5e4, 1e5)}
    chk("P9", "PRINTED-TEXT",
        "'We found that >80 LNPs showed strong in vitro transfection' "
        "(Extended Figure 1). The article gives no threshold, so the count is "
        "reported against several and the claim is NOT scored.",
        "UNSCORED", n_lipids=len(ev), counts_at_RFU_threshold=strong)

    # ---- P10 orientation of the 10 x 10 matrix --------------------------
    # If rows and columns were read the wrong way round, all 100 values would
    # sit on the wrong lipids and NOTHING else in this curation would notice:
    # the structures are right, the labels are right, only the join is wrong.
    # The article states which factor dominates in vitro, so the data itself
    # can be asked.
    import math
    import statistics as stx
    emap = {r["label"]: (r["Am"], r["G"]) for r in enum}
    cell = {}
    for v in ev:
        if v["lipid_label"] in emap and float(v["mean"]) > 0:
            cell[emap[v["lipid_label"]]] = math.log10(float(v["mean"]))
    gm = stx.mean(cell.values())
    sst = sum((x-gm)**2 for x in cell.values())

    def explained(which):
        grp = {}
        for (am, g), x in cell.items():
            grp.setdefault(am if which == "Am" else g, []).append(x)
        return sum(len(xs)*(stx.mean(xs)-gm)**2 for xs in grp.values())/sst

    e_am, e_g = explained("Am"), explained("G")
    chk("P10", "PRINTED-TEXT",
        "'We chose to generate batches based on the tail structure, as this "
        "had the largest effect on mRNA delivery in vitro (Extended Figure "
        "1)'. One-way variance explained in log10(RFU) by the tail and by the "
        "headgroup. This is also the ORIENTATION check on the 10 x 10 matrix: "
        "transposing rows and columns would put all 100 values on the wrong "
        "lipids and would swap these two numbers, and nothing else here would "
        "notice.",
        "PASS" if e_g > e_am else "FAIL",
        variance_explained_by_tail_G=round(e_g, 3),
        variance_explained_by_headgroup_Am=round(e_am, 3),
        n_cells_used=len(cell),
        n_cells_excluded_zero_RFU=100-len(cell))

    QC["counts"] = {
        "AMG labels characterised in the SI": len(
            {r["label"] for r in structures}),
        "HRMS values printed in the SI": len(
            json.load(open(os.path.join(HERE, "02_hrms.json")))),
        "structures reconstructed from the drawings": len(structures),
        "lipids with a usable structure": len(built),
        "amine blocks": len([r for r in blocks
                             if r["kind"] == "amine_headgroup"]),
        "tail blocks": len([r for r in blocks
                            if r["kind"] != "amine_headgroup"]),
        "per-lipid values (tier 1, publisher Source Data)": len(values),
        "in vitro values": len([v for v in values if v["model"] == "in_vitro"]),
        "in vivo values": len([v for v in values if v["model"] == "in_vivo"]),
    }
    QC["value_tier"] = "1_publisher_source_data for every value; no pixel " \
                       "was read anywhere in this curation"
    json.dump(QC, open(os.path.join(HERE, "07_qc.json"), "w"), indent=1)

    for c in QC["checks"]:
        print(f'{c["id"]:<5} {c["kind"]:<15} {c["result"]:<9} '
              + "; ".join(f"{k}={v}" for k, v in c.items()
                          if k not in ("id", "kind", "what", "result")
                          and not isinstance(v, (list, dict)))[:110])
    return 0


if __name__ == "__main__":
    sys.exit(main())
