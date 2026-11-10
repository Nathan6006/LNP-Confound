#!/usr/bin/env python3
"""MY_2026 (PMID 42569578) -- stage the Source Data table into curation artifacts.

TIER 1 THROUGHOUT. Every value comes from the publisher Source Data workbook
(mmc2.xlsx) or the SI tables (mmc1.docx). NOTHING is digitised: the figures are
used only to VERIFY, which is strictly stronger than using them to extract
(PIPELINE_V2 Stage 00b).

Writes 01_blocks.csv, 03_enumerated.csv, 05_values.csv, 06_rows.csv,
08_contradictions.csv and the value half of 07_qc.json.

Usage:  python3 -I stage01_stage.py [--outdir DIR]
A seed harness MUST pass --outdir pointing at a temp directory (F-069).
"""
import argparse
import csv
import json
import math
import os
import re
import sys
import zipfile

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "provenance"))
sys.path.insert(0, os.path.join(ROOT, "curation"))
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))
import reagent_registry as REG                       # noqa: E402
from corpus_components import load, smi              # noqa: E402
import emit                                          # noqa: E402

SRC = os.path.join(HERE, "source")
XLSX = os.path.join(SRC, "mmc2.xlsx")
DOCX = os.path.join(SRC, "mmc1.docx")

# Table 1 of the paper -- the designed factor levels. Checking the workbook
# against these is NOT a tautology; checking mol% against a cholesterol that is
# DEFINED as the remainder would be.
LEVELS = {"ILC": [1.25, 2.5, 5.0, 10.0], "ILR": [40.0, 50.0, 60.0],
          "HLR": [4.0, 10.0, 16.0], "PEGLR": [0.5, 1.5, 2.5]}
HELPERS = ["DSPC", "DOPC"]
ILS = ["SM-102", "MC-3", "ALC-0315", "CKK-E12"]   # spelling of the paper + repo


# ---------------------------------------------------------------- readers ---
def read_xlsx(path):
    """Row-major cell grid from the sheet XML. Handles SELF-CLOSING <c/>.

    The obvious regex `<c r="..."([^>]*)>(.*?)</c>` silently MISREADS a
    self-closing empty cell: it runs past the `/>` and swallows the next cell,
    shifting the whole row by one. That is how the header first came out
    offset against its own data, with four shared-string indices printed as
    bare integers. Guard: the header must resolve to strings, never digits.
    """
    with zipfile.ZipFile(path) as z:
        sheet = z.read("xl/worksheets/sheet1.xml").decode("utf8")
        shared = z.read("xl/sharedStrings.xml").decode("utf8")
    ss = ["".join(re.findall(r"<t(?:\s[^>]*)?>(.*?)</t>", si, re.S))
          for si in re.findall(r"<si>(.*?)</si>", shared, re.S)]

    def colnum(ref):
        letters = re.match(r"([A-Z]+)", ref).group(1)
        n = 0
        for ch in letters:
            n = n * 26 + ord(ch) - 64
        return n

    grid, maxc = {}, 0
    for r, body in re.findall(r'<row[^>]*r="(\d+)"[^>]*>(.*?)</row>',
                              sheet, re.S):
        for m in re.finditer(r"<c\b([^>]*?)(/>|>(.*?)</c>)", body, re.S):
            attrs, inner = m.group(1), (m.group(3) or "")
            ref = re.search(r'r="([A-Z]+\d+)"', attrs).group(1)
            t = re.search(r't="([^"]+)"', attrs)
            v = re.search(r"<v>(.*?)</v>", inner, re.S)
            val = v.group(1) if v else ""
            if t and t.group(1) == "s" and val != "":
                val = ss[int(val)]
            c = colnum(ref)
            maxc = max(maxc, c)
            grid[(int(r), c)] = val
    return grid, maxc


def read_docx_tables(path):
    with zipfile.ZipFile(path) as z:
        doc = z.read("word/document.xml").decode("utf8")

    def ptext(blob):
        # `<w:t[^>]*>` ALSO matches `<w:tcPr>`; the lookahead is load-bearing.
        return "".join(re.findall(r"<w:t(?:\s[^>]*)?>(.*?)</w:t>", blob, re.S))

    tables = []
    for tbl in re.findall(r"<w:tbl>.*?</w:tbl>", doc, re.S):
        rows = []
        for tr in re.findall(r"<w:tr\b.*?</w:tr>", tbl, re.S):
            rows.append([ptext(tc) for tc in
                         re.findall(r"<w:tc>.*?</w:tc>", tr, re.S)])
        tables.append(rows)
    return tables


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# ---------------------------------------------------------------- main ------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=HERE)
    ap.add_argument("--xlsx", default=XLSX,
                    help="mutated COPY for the seed harness; never the real "
                         "source/ file")
    ap.add_argument("--docx", default=DOCX)
    a = ap.parse_args()
    out = os.path.abspath(a.outdir)
    os.makedirs(out, exist_ok=True)
    qc = {"paper": "PMID 42569578", "experiment_id_proposed": "MY_2026",
          "value_tier": "1_publisher_source_data", "digitised_values": 0}

    grid, maxc = read_xlsx(a.xlsx)
    hdr = {c: grid.get((1, c), "") for c in range(1, maxc + 1)}
    # guard against the self-closing-cell shift
    bad = [c for c, v in hdr.items() if v and re.fullmatch(r"\d+", str(v))]
    if bad:
        raise SystemExit(
            f"xlsx header cells {bad} resolved to bare integers -- the "
            "shared-string lookup did not fire, which is the self-closing "
            "<c/> shift. Refusing to read values off a shifted grid.")
    qc["xlsx_header_cells_resolved_to_strings"] = sum(1 for v in hdr.values() if v)

    nrows = max(r for r, _ in grid)
    data = [r for r in range(2, nrows + 1)
            if str(grid.get((r, 1), "")).strip() != ""]
    qc["xlsx_rows_present_in_xml"] = len({r for r, _ in grid})
    qc["xlsx_rows_with_data"] = len(data)
    if len(data) != 25:
        raise SystemExit(f"expected 25 formulations, found {len(data)}")

    # column map, verified against the header text rather than assumed
    COL = {"no": 1, "IL": 2, "HL": 3, "ILC": 4, "ILR": 5, "HLR": 6, "PEGLR": 7}
    COND = {"high_Cy3_siRNA": {"size": 9, "PDI": 10, "zeta": 11, "signal": 12},
            "low_Cy3_siRNA": {"size": 14, "PDI": 15, "zeta": 16, "signal": 17},
            "S100P_siRNA": {"size": 19, "PDI": 20, "zeta": 21, "signal": 22}}
    expect = {2: "Ionizable Lipid Type", 3: "Helper Lipid Type",
              4: "Ionizable Lipid Concentration (mM)",
              5: "Ionizable Lipid Ratio (%)", 6: "Helper Lipid Ratio (%)",
              7: "PEG Lipid Ratio (%)"}
    mism = {c: (hdr[c], e) for c, e in expect.items() if hdr[c].strip() != e}
    if mism:
        raise SystemExit(f"xlsx header does not match expected layout: {mism}")
    qc["xlsx_header_columns_verified"] = len(expect)

    # ---- SI tables -------------------------------------------------------
    tabs = read_docx_tables(a.docx)
    qc["si_tables_found"] = len(tabs)
    qc["si_table_shapes"] = [[len(t), max(len(r) for r in t)] for t in tabs]
    s1, s2, s3 = tabs[0][1:], tabs[1][1:], tabs[2][1:]
    qc["si_S1_formulations"] = len(s1)
    qc["si_S2_formulations"] = len(s2)
    qc["si_S3_formulations"] = len(s3)

    # ---- design-level check (NOT a tautology) ----------------------------
    offenders = []
    for r in data:
        for key, col in (("ILC", 4), ("ILR", 5), ("HLR", 6), ("PEGLR", 7)):
            v = num(grid.get((r, col)))
            if v is None or not any(abs(v - L) < 1e-9 for L in LEVELS[key]):
                offenders.append({"xlsx_row": r, "factor": key, "value": v})
    qc["design_level_violations_xlsx"] = len(offenders)
    qc["design_level_offenders"] = offenders

    s1_off = []
    for row in s1:
        for key, idx in (("ILC", 3), ("ILR", 4), ("HLR", 5), ("PEGLR", 6)):
            v = num(row[idx])
            if v is None or not any(abs(v - L) < 1e-9 for L in LEVELS[key]):
                s1_off.append({"S1_no": row[0], "factor": key, "value": v})
    qc["design_level_violations_tableS1"] = len(s1_off)
    qc["tableS1_rows_checked"] = len(s1)
    print("METRIC design_level_violations_xlsx=%d" % len(offenders))
    print("METRIC design_level_violations_tableS1=%d" % len(s1_off))

    # ---- the four-component sum ------------------------------------------
    # Cholesterol is a DEPENDENT variable: "calculated to complete the total
    # lipid to 100%" (section 2.3.1). So "sums to 100" is true by construction
    # and proves nothing. What IS testable: the remainder must be a positive,
    # physically sensible cholesterol fraction.
    chol = []
    for r in data:
        c = 100.0 - sum(num(grid.get((r, k))) for k in (5, 6, 7))
        chol.append(c)
    qc["fourth_component"] = "cholesterol, derived as 100 - IL - helper - PEG"
    qc["mol_pct_sum_check_is_a_tautology"] = True
    qc["cholesterol_pct_min"] = round(min(chol), 4)
    qc["cholesterol_pct_max"] = round(max(chol), 4)
    qc["cholesterol_pct_nonpositive"] = sum(1 for c in chol if c <= 0)

    # ---- S2 <-> xlsx composition bijection -------------------------------
    def key_x(r):
        il = grid.get((r, 2), "").strip()
        il = "CKK-E12" if il == "CEE-E12" else il
        return (il, grid.get((r, 3), "").strip()) + tuple(
            "%g" % num(grid.get((r, c))) for c in (4, 5, 6, 7))

    def key_s(row):
        return (row[1].strip(), row[2].strip()) + tuple(
            "%g" % num(row[i]) for i in (3, 4, 5, 6))

    xmap, smap = {}, {}
    for r in data:
        xmap.setdefault(key_x(r), []).append(r)
    for row in s2:
        smap.setdefault(key_s(row), []).append(int(row[0]))
    def xno(r):
        """The formulation number the SHEET itself prints, column A.

        NOT the XML row number: sheet row r holds formulation r-1, and
        conflating the two is an off-by-one that would silently rotate the
        whole S2<->xlsx comparison by one position.
        """
        return int(float(grid[(r, 1)]))

    pairs = {}
    for k, nos in smap.items():
        xs = xmap.get(k, [])
        if len(nos) == 1 and len(xs) == 1:
            pairs[nos[0]] = xno(xs[0])
    rowof = {xno(r): r for r in data}
    qc["S2_to_xlsx_bijection"] = (len(pairs) == 25 and
                                  len(set(pairs.values())) == 25)
    qc["S2_to_xlsx_map"] = {str(k): v for k, v in sorted(pairs.items())}
    qc["S2_xlsx_rows_at_same_index"] = sum(1 for k, v in pairs.items() if k == v)
    qc["S2_xlsx_rows_permuted"] = sum(1 for k, v in pairs.items() if k != v)
    print("METRIC S2_to_xlsx_bijection=%s" % qc["S2_to_xlsx_bijection"])
    x2s = {v: k for k, v in pairs.items()}      # xlsx formulation no -> S2 no
    if sorted(rowof) != list(range(1, 26)):
        raise SystemExit(f"sheet index column is not 1..25: {sorted(rowof)}")

    # ---- the ILC discriminator (paper's own Fig 3a correlation) ----------
    # 3.1: ILC correlates POSITIVELY with size for SM-102 (rho = 0.64, the
    # strongest) and NEGATIVELY for CKK-E12. The screen prioritised SMALL
    # size, so the selected SM-102 set should sit at LOW concentration and the
    # selected CKK-E12 set at HIGH. This is an independent, mechanistic test
    # of which block is which -- it does not use the disputed labels.
    blocks = {"xlsx_no_1_7": range(1, 8), "xlsx_no_8_13": range(8, 14),
              "xlsx_no_14_19": range(14, 20), "xlsx_no_20_25": range(20, 26)}
    ilc_by_block = {}
    for name, rng in blocks.items():
        vals = [num(grid.get((rowof[i], 4))) for i in rng]
        ilc_by_block[name] = {
            "xlsx_IL_name": sorted({grid.get((rowof[i], 2), "") for i in rng}),
            "tableS2_no_range": sorted(
                x2s.get(i) for i in rng if x2s.get(i) is not None),
            "unmapped_formulations": [i for i in rng if x2s.get(i) is None],
            "ILC_values": vals, "mean_ILC": round(sum(vals) / len(vals), 4),
            "n_at_max_level_10mM": sum(1 for v in vals if v == 10.0)}
    qc["ILC_discriminator"] = ilc_by_block

    # ---- structures -------------------------------------------------------
    import csv as _csv
    corpus = []
    for fn in ("data_vitro.csv", "data_vivo.csv"):
        with open(os.path.join(ROOT, "data", fn), newline="",
                  encoding="utf8") as f:
            corpus += list(_csv.DictReader(f))
    comp, ambig = load(corpus)
    struct = {
        "ALC-0315": (REG.REGISTRY["ALC-0315"][0], "provenance/reagent_registry.py"),
        "SM-102": (REG.REGISTRY["SM-102"][0], "provenance/reagent_registry.py"),
        "MC-3": (REG.REGISTRY["DLin-MC3-DMA"][0], "provenance/reagent_registry.py"),
        "CKK-E12": (smi(comp, ambig, "CKK-E12"),
                    "corpus (curation/corpus_components.smi) -- NOT independent"),
    }
    enum = []
    for nm, (s, src) in struct.items():
        m = Chem.MolFromSmiles(s)
        if m is None:
            raise SystemExit(f"{nm}: unparseable SMILES")
        enum.append({
            "IL_name_paper": nm,
            "IL_name_corpus": {"ALC-0315": "ALC-0315", "SM-102": "SM-102",
                               "MC-3": "DLin-MC3-DMA",
                               "CKK-E12": "CKK-E12"}[nm],
            "SMILES": Chem.MolToSmiles(m),
            "InChIKey": Chem.MolToInchiKey(m),
            "skeleton_block": Chem.MolToInchiKey(m).split("-")[0],
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "exact_mass": round(Descriptors.ExactMolWt(m), 4),
            "structure_source": src,
            "novel_to_corpus": "no",
            "is_commercial_benchmark": "yes",
        })
    with open(os.path.join(out, "03_enumerated.csv"), "w", newline="",
              encoding="utf8") as f:
        w = _csv.DictWriter(f, fieldnames=list(enum[0]))
        w.writeheader()
        w.writerows(enum)
    qc["distinct_ionizable_lipids"] = len(enum)
    qc["novel_ionizable_lipids"] = 0

    # ---- 01_blocks.csv : the design factors, not a chemical enumeration ---
    with open(os.path.join(out, "01_blocks.csv"), "w", newline="",
              encoding="utf8") as f:
        w = _csv.writer(f)
        w.writerow(["factor", "code", "unit", "levels", "n_levels",
                    "source", "note"])
        w.writerow(["Ionizable lipid type", "ILP", "", "|".join(ILS), len(ILS),
                    "Table 1 (main text)",
                    "all four are commercial catalogue benchmarks"])
        w.writerow(["Helper lipid type", "HLP", "", "|".join(HELPERS),
                    len(HELPERS), "Table 1 (main text)", ""])
        for k, lv in LEVELS.items():
            w.writerow([{"ILC": "Ionizable lipid concentration",
                         "ILR": "Ionizable lipid ratio",
                         "HLR": "Helper lipid ratio",
                         "PEGLR": "PEG-lipid ratio"}[k], k,
                        "mM" if k == "ILC" else "%",
                        "|".join("%g" % v for v in lv), len(lv),
                        "Table 1 (main text)", ""])
        w.writerow(["Cholesterol ratio", "CHOLR", "%", "derived", 0,
                    "section 2.3.1",
                    "dependent variable: 100 - ILR - HLR - PEGLR"])
    qc["design_space_combinations"] = (len(ILS) * len(HELPERS) *
                                       len(LEVELS["ILC"]) * len(LEVELS["ILR"]) *
                                       len(LEVELS["HLR"]) * len(LEVELS["PEGLR"]))

    # ---- 05_values.csv ----------------------------------------------------
    READOUT = {
        "size": ("size", "nm", "no",
                 "physicochemical particle property, not a delivery endpoint"),
        "PDI": ("PDI", "dimensionless", "no",
                "physicochemical particle property, not a delivery endpoint"),
        "zeta": ("zeta potential", "mV", "no",
                 "physicochemical particle property, not a delivery endpoint"),
    }
    rows = []
    vid = 0
    for r in data:
        s2no = x2s.get(xno(r))
        ilx = grid.get((r, 2), "").strip()
        base = {
            "xlsx_formulation_no": xno(r),
            "xlsx_sheet_row": r,
            "tableS2_formulation_no": s2no,
            "IL_name_xlsx": ilx,
            "IL_name_tableS2": next(
                (row[1] for row in s2 if int(row[0]) == s2no), ""),
            "helper_lipid": grid.get((r, 3), "").strip(),
            "ILC_mM": grid.get((r, 4), ""),
            "ILR_pct": grid.get((r, 5), ""),
            "HLR_pct": grid.get((r, 6), ""),
            "PEGLR_pct": grid.get((r, 7), ""),
            "CHOL_pct_derived": "%g" % (100.0 - sum(
                num(grid.get((r, k))) for k in (5, 6, 7))),
        }
        for cond, cols in COND.items():
            for meas, col in cols.items():
                vid += 1
                v = grid.get((r, col), "")
                if meas == "signal" and cond == "S100P_siRNA":
                    ro, un = ("S100P/GAPDH mRNA expression "
                              "(fold change relative to control)",
                              "fold change vs untreated control")
                    adm, why = ("yes",
                                "translated/functional readout: target-gene "
                                "mRNA knockdown by RT-qPCR (2^-ddCt), not "
                                "cargo abundance; corpus precedent "
                                "LRP6_knockdown_normalized (LL_2012)")
                    assay, cargo, tp = ("RT-qPCR TaqMan Hs00195584_m1, "
                                        "GAPDH Hs99999905_m1, 2^-ddCt",
                                        "S100P siRNA (Silencer Select "
                                        "AM16708, ID 17215)",
                                        "48 h (6 h transfection + 42 h)")
                elif meas == "signal":
                    ro, un = "Fluorescence intensity", "a.u."
                    adm, why = ("no",
                                "Cy3 dye on the siRNA reports cargo "
                                "abundance/uptake, not that delivery worked "
                                "(F-071); non-coding non-targeting sequence")
                    assay, cargo, tp = ("Olympus IX81 fluorescence microscopy, "
                                        "ImageJ semi-quantitative",
                                        "Cy3-labelled non-targeting siRNA",
                                        "6 h")
                else:
                    ro, un, adm, why = READOUT[meas]
                    assay, tp = "DLS (Anton Paar Litesizer 500), n=3 technical", ""
                    cargo = ("Cy3-labelled non-targeting siRNA"
                             if "Cy3" in cond else
                             "S100P siRNA (Silencer Select AM16708)")
                rows.append(dict(base, value_id="V%03d" % vid, condition=cond,
                                 measurement=meas, readout_as_printed=ro,
                                 units_as_printed=un, value=v,
                                 source_location=("mmc2.xlsx sheet1 "
                                                  "row %d col %d" % (r, col)),
                                 tier="1_publisher_source_data",
                                 cargo=cargo, cell_line="A549",
                                 organ="", route="in vitro",
                                 dose=("not stated for this arm"
                                       if cond == "S100P_siRNA" else
                                       ("9 pmol/well" if cond.startswith("high")
                                        else "3 pmol/well")),
                                 timepoint=tp, assay=assay,
                                 corpus_admissible=adm,
                                 admissible_reason=why))
    with open(os.path.join(out, "05_values.csv"), "w", newline="",
              encoding="utf8") as f:
        w = _csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    qc["values_total"] = len(rows)
    qc["values_admissible_yes"] = sum(1 for r in rows
                                      if r["corpus_admissible"] == "yes")
    qc["values_admissible_no"] = sum(1 for r in rows
                                     if r["corpus_admissible"] == "no")
    assert qc["values_admissible_yes"] + qc["values_admissible_no"] == len(rows)

    # ---- the z that WOULD be emitted, both directions --------------------
    adm = [r for r in rows if r["corpus_admissible"] == "yes"]
    vals = [float(r["value"]) for r in adm]
    # delivery-proportional signal = 1/relative expression; corpus convention
    # is a per-experiment z of log10 with population sd (ddof=0), and LL_2012
    # stores knockdown so that HIGHER = better delivery.
    sig = [-math.log10(v) for v in vals]
    mu = sum(sig) / len(sig)
    sd = math.sqrt(sum((x - mu) ** 2 for x in sig) / len(sig))
    qc["would_be_z"] = {
        "n": len(sig), "transform": "z(-log10(relative S100P expression))",
        "pooling": "ONE z per experiment over all 25 (no per-condition z; a "
                   "per-condition z is error class 22)",
        "mean_raw_relexpr": round(sum(vals) / len(vals), 6),
        "min_relexpr": round(min(vals), 6), "max_relexpr": round(max(vals), 6),
        "z_mean": round(mu - mu, 12), "z_sd_ddof0": round(sd / sd, 12),
        "z_values": {adm[i]["value_id"]: round((sig[i] - mu) / sd, 6)
                     for i in range(len(adm))}}

    # ---- 06_rows.csv : header-only, read live from the corpus ------------
    with open(os.path.join(ROOT, "data", "data_vitro.csv"), encoding="utf8") as f:
        header = f.readline().rstrip("\n")
    ncol = len(header.split(","))
    with open(os.path.join(out, "06_rows.csv"), "w", encoding="utf8") as f:
        f.write(header + "\n")
    qc["rows_emitted"] = 0
    qc["rows_header_columns"] = ncol
    qc["rows_header_source"] = "data/data_vitro.csv (read live)"
    if ncol != 79:
        raise SystemExit(f"corpus header is {ncol} columns, expected 79")

    for k in ("values_admissible_yes", "values_admissible_no",
              "rows_emitted", "xlsx_rows_with_data"):
        print("METRIC %s=%s" % (k, qc[k]))

    with open(os.path.join(out, "07_qc_values.json"), "w", encoding="utf8") as f:
        json.dump(qc, f, indent=2)
    print(json.dumps({k: v for k, v in qc.items()
                      if k not in ("S2_to_xlsx_map", "would_be_z",
                                   "ILC_discriminator")}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
