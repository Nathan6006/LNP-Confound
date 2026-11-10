#!/usr/bin/env python3
"""Build 05_values.csv -- EVERY value extracted from this paper, each with an
explicit corpus_admissible yes/no and the reason.

ADMISSIBILITY, stated once here and repeated per row.

data/'s Experiment_method vocabulary is DELIVERY-ONLY (luminescence /
protein-abundance / editing-efficiency / knockdown). Three classes of number in
this paper are therefore inadmissible regardless of quality:

  * formulation parameters (Table S1 molar ratios) -- an input, not a measurement;
  * physicochemical characterisation (Table S2 size / PDI / EE / RNA conc) --
    measured, but not delivery;
  * dispersion statistics (Tables S3, S4 coefficients of variation) -- a CV
    carries no mean and is not a measurement of anything delivered.

THE POOLED BARCODE SCREEN IS NOT INADMISSIBLE, AND AN EARLIER VERSION OF THIS
DOCSTRING SAID IT WAS. The paper's headline screen administers the b-LNP library
POOLED and reads mRNA barcodes by deep sequencing, so each value is a relative,
competitive share of its own pool. That is a real property -- but it is one the
corpus RECORDS rather than excludes: data/data_vivo.csv carries an
`Experiment_batching` column with vocabulary {barcoded, individual}, and
LR_2023 already contributes 407 barcoded in-vivo rows over 383 unique ionizable
lipids. The screen's values are simply NOT ACQUIRED -- they are in main-article
Figure 2, behind an ACS Cloudflare interstitial this run was not permitted to
pass. So 0 pooled values appear below for an ACQUISITION reason, not a design
one. See DECISIONS.md D-1 and D-3.

What IS admissible is the SEPARATE SINGLE-LNP VALIDATION: eight formulations
re-made individually, injected individually into their own mice, and read as
IVIS luciferase total flux. That is the same kind of measurement the corpus
already holds. SI Fig S2A carries it for six organs.
"""
import csv, json, os, re, sys
import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "nn5c19138_si_001.pdf")

# The 15 ILs in the order Tables S1/S2 lay them out, and the b-LNP numbers in
# each block. Block membership is 3 IL mol% x 2 PEG mol% = 6, in the printed
# order (35,2.5) (45,2.5) (55,2.5) (35,5) (45,5) (55,5). The GAPS are real:
# Table S1's caption says a greyed LNP "was not included in the pooled
# administration ... due to not meeting the selection criteria".
BLOCKS = [
    ("C10-306", [1, 2, 3, 4, 5, 6]),      ("C12-306", [7, 8, 9, 10, 11, 12]),
    ("C14-306", [13, 14, 15, 16, 17, 18]), ("C10-114", [20, 21, 22, 23, 24, 25]),
    ("C12-114", [26, 27, 28, 29, 30, 31]), ("C14-114", [32, 33, 34, 35, 36, 37]),
    ("C10-494", [38, 39, 40, 42, 43, 44]), ("C12-494", [45, 46, 47, 48, 49, 50]),
    ("C14-494", [51, 52, 53, 54, 58, 59]), ("C10-200", [60, 61, 62, 63, 64, 65]),
    ("C12-200", [66, 68, 69, 70, 71, 72]), ("C14-200", [73, 74, 75, 76, 77, 78]),
    ("C10-490", [79, 80, 81, 82, 83, 84]), ("C12-490", [85, 86, 87, 88, 89, 90]),
    ("C14-490", [91, 92, 94, 96, 97, 98]),
]
BENCH = {100: "DLin-MC3-DMA", 101: "ALC-0315", 102: "SM-102"}
VALID = ["13", "14", "70", "96", "97", "98", "45", "102"]   # Table S3 names these


def grab(pages):
    rows = []
    with pdfplumber.open(SI) as pdf:
        for p in pages:
            for tb in pdf.pages[p - 1].find_tables():
                for r in tb.extract():
                    rows.append(r)
    return rows


def parse(pages, cols):
    """Rows of (b_LNP:int, {col: str}), plus the sparse IL labels seen."""
    out, labels = {}, []
    for r in grab(pages):
        if not r or len(r) < 6:
            continue
        lab, num = (r[0] or "").strip(), (r[1] or "").strip()
        if lab and re.match(r"^C1[024]-\d+$", lab):
            labels.append((lab, num))
        if lab in ("MC3", "ALC-0315", "SM-102"):
            labels.append((lab, num))
        if not re.match(r"^\d+$", num):
            continue
        out[int(num)] = {c: (r[2 + i] or "").strip() for i, c in enumerate(cols)}
    return out, labels


def main():
    s1, lab1 = parse([22, 23, 24], ["IL", "DOPE", "Chol", "PEG"])
    s2, lab2 = parse([25, 26, 27], ["Size_nm", "PDI", "EE_pct", "RNA_ng_uL"])

    # b-LNP -> IL, from BLOCK STRUCTURE. Then CHECKED against the sparse labels
    # pdfplumber actually read, so the block map is tested, not trusted.
    il = {}
    for name, nums in BLOCKS:
        for n in nums:
            il[n] = name
    for n, nm in BENCH.items():
        il[n] = nm

    qc = {}
    seen = {}
    for lab, num in lab1 + lab2:
        if re.match(r"^\d+$", num or ""):
            seen.setdefault(lab, set()).add(int(num))
    mism = []
    for lab, nums in seen.items():
        key = "DLin-MC3-DMA" if lab == "MC3" else lab
        for n in nums:
            if il.get(n) != key:
                mism.append((lab, n, il.get(n)))
    qc["P4_block_map_vs_printed_labels"] = dict(
        n_labels_anchored=sum(len(v) for v in seen.values()),
        n_mismatched=len(mism), mismatches=mism,
        what=("The IL column is sparse -- one label per 6-row block, placed "
              "wherever the cell's text box sits. The IL for each b-LNP is "
              "therefore derived from BLOCK STRUCTURE (15 ILs x 6 "
              "formulations, with the greyed-out numbers skipped), and every "
              "label pdfplumber did read must land in the block it was "
              "assigned to."))

    # molar-ratio pattern implied by the design, checked against what is printed
    PAT = [("35", "2.5"), ("45", "2.5"), ("55", "2.5"),
           ("35", "5"), ("45", "5"), ("55", "5")]
    bad = []
    for name, nums in BLOCKS:
        for pos, n in enumerate(nums):
            if n not in s1:
                bad.append((n, "absent from Table S1")); continue
            want_il, want_peg = PAT[pos]
            got = (s1[n]["IL"], s1[n]["PEG"])
            if got != (want_il, want_peg):
                bad.append((n, f"printed {got}, design implies {(want_il, want_peg)}"))
    qc["P5_design_grid_matches_printed_molar_ratios"] = dict(
        n_checked=sum(len(b[1]) for b in BLOCKS), n_mismatched=len(bad),
        mismatches=bad,
        what=("Every block must run (35,2.5)(45,2.5)(55,2.5)(35,5)(45,5)(55,5). "
              "This is a 90-row test of the block map AND of the table read: a "
              "single mis-assigned b-LNP number breaks the pattern."))
    # DOPE and Chol are constant by design
    dope = sorted({s1[n]["DOPE"] for _, ns in BLOCKS for n in ns if n in s1})
    chol = sorted({s1[n]["Chol"] for _, ns in BLOCKS for n in ns if n in s1})
    qc["P6_helper_and_sterol_constant_across_library"] = dict(
        DOPE_values=dope, Chol_values=chol,
        constant=(len(dope) == 1 and len(chol) == 1))
    # the implied round number: the four mol% must sum to 100
    sums = {}
    for n, r in s1.items():
        try:
            sums[n] = round(sum(float(r[c]) for c in ("IL", "DOPE", "Chol", "PEG")), 6)
        except ValueError:
            pass
    qc["P7_molar_percents_sum_to_100"] = dict(
        n=len(sums), n_exactly_100=sum(1 for v in sums.values() if v == 100.0),
        distinct_sums=sorted(set(sums.values())),
        offenders={k: v for k, v in sums.items() if v != 100.0},
        what=("The authors never print this. A molar-ratio table that does not "
              "sum to 100 has a transcription error in it."))

    # ---------------------------------------------------------------- values
    dig = json.load(open(os.path.join(HERE, "04_figS1S2_measurements.json")))
    rows = []

    def add(**kw):
        rows.append(kw)

    SHARED = dict(cargo="luciferase mRNA", species="C57BL/6 mouse",
                  route="intravenous", dose="0.35 mg/kg mRNA",
                  timepoint="6 h", n_biological="4 mice",
                  pool_identity="NOT POOLED - formulated and injected individually")

    for b in dig["bars"]:
        if b["panel"] != "S2A":
            continue
        n = int(b["lnp"])
        add(value_id=f"figS2A/{b['lnp']}/{b['organ']}",
            source_locus="SI Figure S2A (p. 14)", tier="5_pixels",
            quantity_as_printed="Total flux (p/s)", units_as_printed="p/s",
            readout="IVIS bioluminescence from luciferase mRNA",
            assay="in vivo IVIS whole-organ bioluminescence imaging",
            b_LNP=b["lnp"], IL_name=il[n], organ=b["organ"],
            value=f"{b['value']:.6g}", **SHARED,
            corpus_admissible="yes",
            admissible_reason=("Single-LNP in-vivo luciferase luminescence: one "
                               "formulation per animal, not pooled, so the value "
                               "does not depend on other library members. Same "
                               "kind of measurement as the corpus's in-vivo "
                               "luminescence rows."))
    for b in dig["bars"]:
        if b["panel"] != "S1A":
            continue
        n = int(b["lnp"])
        add(value_id=f"figS1A/{b['lnp']}/{b['organ']}",
            source_locus="SI Figure S1A (p. 13)", tier="5_pixels",
            quantity_as_printed="% of total luminescence",
            units_as_printed="percent",
            readout="IVIS bioluminescence, organ signal / total organ signal x 100",
            assay="in vivo IVIS whole-organ bioluminescence imaging",
            b_LNP=b["lnp"], IL_name=il[n], organ=b["organ"],
            value=f"{b['value']:.6g}", **SHARED,
            corpus_admissible="no",
            admissible_reason=("SAME 48 measurements as Fig S2A, re-expressed as "
                               "a within-animal percentage. Emitting both would "
                               "enter one experiment twice. Used as the "
                               "independent cross-check P11/P12, not as a value."))

    # The three benchmark rows are printed under their OWN column headers in
    # Table S1 -- DSPC not DOPE, and DMG-PEG2k / ALC-0159 not C14-PEG2k. The
    # parse is positional, so the component NAME is restored here; the numbers
    # are untouched. Mislabelling a helper lipid is error class 30 territory.
    BENCH_HL = {100: "DSPC", 101: "DSPC", 102: "DSPC"}
    BENCH_PEG = {100: "DMG-PEG2000", 101: "ALC-0159", 102: "DMG-PEG2000"}
    for n, r in sorted(s1.items()):
        hl = BENCH_HL.get(n, "DOPE")
        pg = BENCH_PEG.get(n, "C14-PEG2000")
        for comp, col in (("ionizable lipid", "IL"), (hl, "DOPE"),
                          ("cholesterol", "Chol"), (pg, "PEG")):
            add(value_id=f"tableS1/{n}/{col}",
                source_locus="SI Table S1 (pp. 22-24)", tier="4_pdf_text_table",
                quantity_as_printed=f"Molar ratio, {comp}", units_as_printed="mol %",
                readout="formulation composition", assay="n/a - formulation input",
                b_LNP=str(n), IL_name=il.get(n, ""), organ="",
                value=r[col], cargo="", species="", route="", dose="",
                timepoint="", n_biological="", pool_identity="",
                corpus_admissible="no",
                admissible_reason=("A FORMULATION PARAMETER, not a measured "
                                   "endpoint. Carried into the emitted rows' "
                                   "formulation columns; never a delivery value."))
    for n, r in sorted(s2.items()):
        for col, unit in (("Size_nm", "nm"), ("PDI", "dimensionless"),
                          ("EE_pct", "%"), ("RNA_ng_uL", "ng/uL")):
            add(value_id=f"tableS2/{n}/{col}",
                source_locus="SI Table S2 (pp. 25-27)", tier="4_pdf_text_table",
                quantity_as_printed=col.replace("_", " "), units_as_printed=unit,
                readout="physicochemical characterisation",
                assay="DLS (size, PDI) / RiboGreen (EE, concentration)",
                b_LNP=str(n), IL_name=il.get(n, ""), organ="",
                value=r[col], cargo="", species="", route="", dose="",
                timepoint="", n_biological="", pool_identity="",
                corpus_admissible="no",
                admissible_reason=("NOT DELIVERY. Size, PDI, encapsulation "
                                   "efficiency and RNA concentration are "
                                   "physicochemical properties; data/'s "
                                   "Experiment_method vocabulary is delivery-only."))

    S3 = {"3C": [76.68, 34.85, 12.03, 54.99, 90.14, 36.57, 52.02, 58.67],
          "3D": [42.40, 67.00, 25.54, 94.73, 47.27, 31.02, 37.81, 21.98],
          "3E": [12.88, 59.18, 0.9918, 1.640, 1.071, 0.3346, 0.5283, 0.2108],
          "3F": [23.22, 31.16, 27.26, 45.37, 69.19, 67.33, 10.41, 18.34],
          "3G": [6.240, 18.39, 11.14, 23.34, 12.13, 12.80, 21.65, 16.82],
          "3H": [2.993, 6.702, 2.855, 2.530, 1.778, 0.3276, 0.8812, 6.160]}
    S4 = {"4A": [27.47, 6.842, 44.38, 41.35, 49.47, 29.11],
          "4B": [33.68, 13.88, 50.96, 35.53, 46.15, 28.69],
          "4C": [17.56, 26.13, 52.63, 26.93, 23.13, 36.74],
          "4D": [40.15, 32.36, 36.31, 41.50, 27.96, 48.54]}
    for panel, vals in S3.items():
        for lnp, cv in zip(VALID, vals):
            add(value_id=f"tableS3/{panel}/{lnp}",
                source_locus="SI Table S3 (p. 28)", tier="4_pdf_text_table",
                quantity_as_printed=f"Coefficient of variation, Figure {panel}",
                units_as_printed="%", readout="dispersion statistic",
                assay="n/a", b_LNP=lnp, IL_name=il[int(lnp)], organ="",
                value=str(cv), cargo="", species="", route="", dose="",
                timepoint="", n_biological="", pool_identity="",
                corpus_admissible="no",
                admissible_reason=("A CV IS NOT A MEASUREMENT. It carries no "
                                   "mean, so no delivery value can be recovered "
                                   "from it. Its use here is that it NAMES the "
                                   "single-LNP validation cohort."))
    for panel, vals in S4.items():
        for grp, cv in zip(["PBS", "13", "14", "97", "45", "102"], vals):
            add(value_id=f"tableS4/{panel}/{grp}",
                source_locus="SI Table S4 (p. 28)", tier="4_pdf_text_table",
                quantity_as_printed=f"Coefficient of variation, Figure {panel}",
                units_as_printed="%", readout="dispersion statistic", assay="n/a",
                b_LNP=grp, IL_name=(il[int(grp)] if grp.isdigit() else "PBS"),
                organ="brain", value=str(cv), cargo="", species="", route="",
                dose="", timepoint="", n_biological="", pool_identity="",
                corpus_admissible="no",
                admissible_reason="A CV is not a measurement (as Table S3).")

    cols = list(rows[0])
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})

    yes = sum(1 for r in rows if r["corpus_admissible"] == "yes")
    qc["value_counts"] = dict(
        total=len(rows), admissible=yes, inadmissible=len(rows) - yes,
        by_locus={})
    for r in rows:
        k = r["source_locus"]
        d = qc["value_counts"]["by_locus"].setdefault(
            k, {"total": 0, "admissible": 0})
        d["total"] += 1
        d["admissible"] += (r["corpus_admissible"] == "yes")
    assert yes + (len(rows) - yes) == len(rows)
    json.dump(qc, open(os.path.join(HERE, "05_values_qc.json"), "w"), indent=1)
    json.dump({"il_by_bLNP": il, "tableS1": s1, "tableS2": s2},
              open(os.path.join(HERE, "05_tables.json"), "w"), indent=1)
    print(json.dumps(qc, indent=1)[:2600])


if __name__ == "__main__":
    main()
