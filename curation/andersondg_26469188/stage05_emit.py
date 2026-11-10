#!/usr/bin/env python3
"""KK_2015 -- write 05_values.csv (every reported value) and 06_rows.csv.

ONE Z PER EXPERIMENT. The 36 EPO formulations are standardised together, as
one block: z of log10(serum EPO ng/mL) with POPULATION sd (ddof=0), the corpus
convention. There is one Model_type and one organ, so no per-cell or per-cohort
z arises and error class 22 cannot occur here.

WHAT IS NOT EMITTED, and why, is in DECISIONS.md: the 2 luciferase rows
(Table S3) and the 6 FVII siRNA rows (Table S4) are admissible readouts but a
different cargo and a different unit, and n=2 / n=2-formulations makes a z of
exactly +/-1 -- a number with no information in it.
"""
import csv
import json
import math
import os
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "provenance"))

import emit                                                       # noqa: E402
import make_reference                                             # noqa: E402
import reagent_registry                                           # noqa: E402
from corpus_components import load, smi                           # noqa: E402

OUTDIR = os.environ.get("KK2015_OUTDIR", HERE)
DATA = os.path.join(HERE, "..", "..", "data")
EXP = "KK_2015"
PMID = "26469188"
DOI = "10.1021/acs.nanolett.5b02497"
LINK = "https://doi.org/" + DOI

P = json.load(open(os.path.join(OUTDIR, "_parsed.json")))
S1 = {r["code"]: r for r in P["table_s1"]}
S3 = {r["code"]: r for r in P["table_s3"]}

# A-03's cholesterol mol% is printed as 31.5, which makes its four components
# sum to 110. See 08_contradictions.csv C-01: the value is OVER-DETERMINED at
# 21.5 -- 13 of 14 Library A rows have chol == 100 - (IL + phospholipid + PEG)
# exactly, and A-12 carries the IDENTICAL (60, 16, 2.5) levels and prints 21.5.
CHOL_CORRECTION = {"A-03": (31.5, 21.5)}


def corpus_rows():
    out = []
    for f in ("data_vitro.csv", "data_vivo.csv"):
        with open(os.path.join(DATA, f), newline="") as fh:
            out.extend(csv.DictReader(fh))
    return out


def header():
    """Read the 79-column header LIVE from the corpus, never restate it."""
    with open(os.path.join(DATA, "data_vivo.csv"), newline="") as fh:
        h = next(csv.reader(fh))
    with open(os.path.join(DATA, "data_vitro.csv"), newline="") as fh:
        h2 = next(csv.reader(fh))
    assert h == h2, "the two corpus headers differ; refusing to guess one"
    return h


# =============================================================== 05_values ===
VAL_COLS = ["value_id", "table", "page", "code", "readout_as_printed",
            "units_as_printed", "value", "sd", "n", "cargo", "cargo_detail",
            "cell_or_organ", "organ", "route", "dose", "timepoint", "assay",
            "il_name", "wt_ratio", "phospholipid", "il_molpct", "ph_molpct",
            "chol_molpct", "peg_molpct", "ee_pct", "size_nm", "pdi",
            "corpus_admissible", "admissible_reason", "emitted", "emit_reason"]

BASE = {"il_name": "C12-200", "route": "intravenous"}


def values():
    rows = []
    n = 0

    def add(**kw):
        nonlocal n
        n += 1
        r = {c: "" for c in VAL_COLS}
        r.update(BASE)
        r["value_id"] = "V-%02d" % n
        r.update(kw)
        rows.append(r)

    # ---- Table S1: serum EPO ------------------------------------------------
    for code in ["PBS", "Original"] + sorted(c for c in S1
                                             if c not in ("PBS", "Original")):
        r = S1[code]
        common = dict(table="S1", page=2, code=code,
                      readout_as_printed="Serum EPO",
                      units_as_printed="ng/mL", n=3, cargo="mRNA",
                      cargo_detail="erythropoietin (EPO) mRNA; species not stated",
                      cell_or_organ="mouse (strain not stated in the SI)",
                      organ="liver (per the abstract: 'delivery to the liver in "
                            "vivo'); the measurement itself is SERUM protein",
                      dose="15 ug total mRNA", timepoint="6 h",
                      assay="serum EPO ELISA (protein abundance)")
        if code == "PBS":
            add(**common, value="", sd="",
                corpus_admissible="no",
                admissible_reason="VEHICLE CONTROL. The readout (serum protein "
                "abundance by ELISA) IS an admissible delivery endpoint, but this "
                "row contains no ionizable lipid, so there is no IL structure to "
                "key the corpus on. Its printed value is 0.018 +/- 0.002 ng/mL.",
                emitted="no", emit_reason="no lipid")
            continue
        if r["unsynthesizable"]:
            add(**common, value="", sd="",
                wt_ratio=r["wt_ratio"], phospholipid=r["phospholipid"],
                il_molpct=r["il_molpct"], ph_molpct=r["ph_molpct"],
                chol_molpct=r["chol_molpct"], peg_molpct=r["peg_molpct"],
                corpus_admissible="no",
                admissible_reason="NO VALUE EXISTS. Table S1 prints '*' and its "
                "footnote states the LNP could not be synthesized due to "
                "insolubility of DSPE in ethanol at all concentrations tested.",
                emitted="no", emit_reason="no measurement")
            continue
        chol = r["chol_molpct"]
        extra = ""
        if code in CHOL_CORRECTION:
            printed, fixed = CHOL_CORRECTION[code]
            chol = fixed
            extra = (" NOTE: cholesterol is printed as %.1f, which makes the four "
                     "components sum to 110; the over-determined value is %.1f "
                     "(08_contradictions.csv C-01)." % (printed, fixed))
        add(**common, value=r["epo_mean"], sd=r["epo_sd"],
            wt_ratio=r["wt_ratio"], phospholipid=r["phospholipid"],
            il_molpct=r["il_molpct"], ph_molpct=r["ph_molpct"],
            chol_molpct=chol, peg_molpct=r["peg_molpct"],
            ee_pct=r["ee_pct"], size_nm=r["size_nm"], pdi=r["pdi"],
            corpus_admissible="yes",
            admissible_reason="Erythropoietin serum protein by ELISA is a "
            "PROTEIN-ABUNDANCE delivery readout." + extra,
            emitted="yes", emit_reason="in the KK_2015 z-block")

    # ---- Table S3: luciferase total flux -----------------------------------
    for code, r in sorted(S3.items()):
        add(table="S3", page=4, code=code + ", Luc",
            readout_as_printed="Total Flux", units_as_printed="(p/s) x 10^9",
            value=r["flux_mean_e9"], sd=r["flux_sd_e9"], n=3, cargo="mRNA",
            cargo_detail="luciferase mRNA",
            cell_or_organ="mouse (strain not stated in the SI)",
            organ="whole body (total flux; no organ stated)",
            dose="15 ug total mRNA", timepoint="6 h",
            assay="whole-animal bioluminescence imaging",
            wt_ratio=r["wt_ratio"], phospholipid=r["phospholipid"],
            il_molpct=r["il_molpct"], ph_molpct=r["ph_molpct"],
            chol_molpct=r["chol_molpct"], peg_molpct=r["peg_molpct"],
            ee_pct=r["ee_pct"], size_nm=r["size_nm"], pdi=r["pdi"],
            corpus_admissible="yes",
            admissible_reason="Luminescence is a delivery readout.",
            emitted="no",
            emit_reason="NOT EMITTED (D-1): only 2 formulations, a different cargo "
            "(Luc vs EPO mRNA) and a different unit ((p/s)x10^9 vs ng/mL). A 2-point "
            "z is exactly +/-1 by construction and carries no information, and "
            "pooling p/s with ng/mL in one z would standardise two incommensurable "
            "quantities.")

    # ---- Table S4: FVII siRNA knockdown ------------------------------------
    for r in P["table_s4"]:
        code = r["code"]
        if code == "PBS":
            add(table="S4", page=5, code="PBS",
                readout_as_printed="Serum FVII", units_as_printed="%",
                value=r["fvii_mean"], sd=r["fvii_sd"], n=3, cargo="siRNA",
                cargo_detail="Factor VII siRNA",
                cell_or_organ="mouse (strain not stated in the SI)",
                organ="liver (Factor VII is hepatocyte-expressed)",
                dose="n/a", timepoint="72 h",
                assay="serum Factor VII activity/protein (knockdown)",
                corpus_admissible="no",
                admissible_reason="VEHICLE CONTROL and the 100% normaliser; "
                "contains no ionizable lipid.",
                emitted="no", emit_reason="no lipid")
            continue
        add(table="S4", page=5, code="%s, siFVII @ %g mg/kg" % (code, r["dose_mg_kg"]),
            readout_as_printed="Serum FVII", units_as_printed="%", n=3,
            value=r["fvii_mean"], sd=r["fvii_sd"], cargo="siRNA",
            cargo_detail="Factor VII siRNA",
            cell_or_organ="mouse (strain not stated in the SI)",
            organ="liver (Factor VII is hepatocyte-expressed)",
            dose="%g mg/kg siRNA" % r["dose_mg_kg"], timepoint="72 h",
            assay="serum Factor VII activity/protein (knockdown)",
            wt_ratio=r["wt_ratio"], phospholipid=r["phospholipid"],
            il_molpct=r["il_molpct"], ph_molpct=r["ph_molpct"],
            chol_molpct=r["chol_molpct"], peg_molpct=r["peg_molpct"],
            ee_pct=r["ee_pct"], size_nm=r["size_nm"], pdi=r["pdi"],
            corpus_admissible="yes",
            admissible_reason="Target-protein knockdown is an admissible delivery "
            "readout.",
            emitted="no",
            emit_reason="NOT EMITTED (D-1): a different cargo (siRNA), a different "
            "readout (% of PBS, floored at 0 and inverted in direction), only 2 "
            "formulations, and the same formulation measured at 3 doses. 6 rows on "
            "2 formulations of ONE ionizable lipid.")
    return rows


# ================================================================ 06_rows ====
def emitted_rows():
    codes = sorted((c for c, r in S1.items()
                    if c not in ("PBS",) and not r["unsynthesizable"]),
                   key=lambda c: (c != "Original", c))
    assert len(codes) == 36, "expected 36 EPO LNPs, got %d" % len(codes)

    logs = {c: math.log10(S1[c]["epo_mean"]) for c in codes}
    mu = statistics.fmean(logs.values())
    sd = statistics.pstdev(list(logs.values()))          # POPULATION sd, ddof=0
    z = {c: (logs[c] - mu) / sd for c in codes}
    assert abs(statistics.fmean(z.values())) < 1e-12, "z block mean != 0"
    assert abs(statistics.pstdev(list(z.values())) - 1.0) < 1e-12, "z block sd != 1"

    comp, amb = load(corpus_rows())
    il = smi(comp, amb, "C12-200")
    helpers = {h: make_reference.build(h)[0] for h in ("DSPC", "DOPE", "DOPC")}
    chol = reagent_registry.REGISTRY["Cholesterol"][0]

    # ---- flags, through emit.emit_rows, with PER-ITEM CALLABLES ------------
    records = {c: dict(S1[c], _code=c,
                       _chol_corrected=c in CHOL_CORRECTION) for c in codes}

    def confidence(item, rec):
        # Derived from the record, never a literal (F-011). Every row shares the
        # same reason because every row shares the same ionizable lipid, and
        # that reason is a property of the record, not a constant in this file.
        return "medium" if rec["_code"] else ""

    def note(item, rec):
        n = ("C12-200's structure is NOT attested by this paper: no drawing, no "
             "formula, no HRMS, no MW, no synthesis anywhere in the SI, and the "
             "main text is closed-access. It is carried over from the corpus "
             "(AC_2025 498 rows + JW_2024 20 rows, InChIKey "
             "RVHYPUORVDKRTM-UHFFFAOYSA-N). The PEG-lipid is never named in the "
             "SI either and is emitted as 'Unknown' with no structure.")
        if rec["_chol_corrected"]:
            p, f = CHOL_CORRECTION[rec["_code"]]
            n += (" Cholesterol mol%% is printed as %.1f (four components sum to "
                  "110); emitted as the over-determined %.1f -- see "
                  "08_contradictions.csv C-01." % (p, f))
        return n

    def structure_assumed(item, rec):
        return "1"

    def positional_ambiguity(item, rec):
        return "0"

    def norm_unit(item, rec):
        return "experiment"

    flags = emit.emit_rows(
        codes, records,
        confidence=confidence, note=note,
        positional_ambiguity=positional_ambiguity,
        structure_assumed=structure_assumed,
        value_is_bound=None,          # no bound in this paper: every value is a
                                      # printed numeral with an SD
        value_normalization_unit=norm_unit,
        within_model=None,            # build_dataset.py computes it
        label="KK_2015 EPO")

    hdr = header()
    out = []
    for c, fl in zip(codes, flags):
        r = S1[c]
        chol_pct = CHOL_CORRECTION[c][1] if c in CHOL_CORRECTION else r["chol_molpct"]
        row = {k: "" for k in hdr}
        row.update({
            "Experiment_ID": EXP,
            "Formulation_ID": c,
            "IL_name": "C12-200", "IL_SMILES": il,
            "IL_molratio": "%g" % r["il_molpct"],
            "IL_to_nucleicacid_massratio": "%g" % r["wt_ratio"],
            "HL_name": r["phospholipid"], "HL_SMILES": helpers[r["phospholipid"]],
            "HL_molratio": "%g" % r["ph_molpct"],
            "CHL_name": "Cholesterol", "CHL_SMILES": chol,
            "CHL_molratio": "%g" % chol_pct,
            "PEG_name": "Unknown", "PEG_SMILES": "",
            "PEG_molratio": "%g" % r["peg_molpct"],
            "Model": "in_vivo", "Model_type": "Mouse_unspecified",
            "Model_target": "liver",
            "Route_of_administration": "intravenous",
            "Cargo": "mRNA", "Cargo_type": "EPO",
            "Dose_ug_nucleicacid": "15",
            "Experiment_method": "protein_adundance_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": "%.6f" % z[c],
            "Publication_link": LINK, "Publication_DOI": DOI,
            "Publication_PMID": PMID,
        })
        row.update(fl)
        out.append(row)
    return hdr, out, {"n": len(codes), "mu_log10": mu, "sd_log10": sd,
                      "z_min": min(z.values()), "z_max": max(z.values())}


def main():
    vals = values()
    with open(os.path.join(OUTDIR, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=VAL_COLS)
        w.writeheader()
        w.writerows(vals)

    hdr, rows, zinfo = emitted_rows()
    with open(os.path.join(OUTDIR, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr)
        w.writeheader()
        w.writerows(rows)

    adm = sum(1 for v in vals if v["corpus_admissible"] == "yes")
    nad = sum(1 for v in vals if v["corpus_admissible"] == "no")
    emd = sum(1 for v in vals if v["emitted"] == "yes")
    summary = {
        "value_rows_total": len(vals),
        "corpus_admissible_yes": adm, "corpus_admissible_no": nad,
        "counts_sum": adm + nad == len(vals),
        "emitted_yes": emd, "emitted_no": len(vals) - emd,
        "emitted_equals_06_rows": emd == len(rows),
        "by_table": {t: sum(1 for v in vals if v["table"] == t)
                     for t in ("S1", "S3", "S4")},
        "z_block": zinfo,
        "n_columns_emitted": len(hdr),
        "caveat_rows": emit.caveat_count([{k: r.get(k, "") for k in
                                           ("confidence", "positional_ambiguity",
                                            "structure_assumed", "value_is_bound")}
                                          for r in rows]),
    }
    with open(os.path.join(OUTDIR, "_emit.json"), "w") as fh:
        json.dump(summary, fh, indent=1, sort_keys=True, default=str)
    print(json.dumps(summary, indent=1, sort_keys=True, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
