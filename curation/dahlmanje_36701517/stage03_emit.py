"""Stage 03 -- values, admissibility, rows, QC. Radmand 2023 (PMID 36701517).

Reads 04_digitized.json (Fig S13 / S15 bar reads) and writes:
    05_values.csv   every value, with source, units as printed, and an
                    explicit corpus_admissible yes/no + reason on EVERY row
    06_rows.csv     the admissible ones, in the 79-column corpus schema,
                    quality flags written by curation/lib/emit.emit_rows with
                    PER-ITEM callables
    07_qc.json      every check, in numbers

Usage: python3 -I stage03_emit.py
"""
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, os.path.join(HERE, ".."))
import emit                                                   # noqa: E402
from corpus_components import load, smi                       # noqa: E402

from rdkit import Chem                                        # noqa: E402
from rdkit import RDLogger                                    # noqa: E402
RDLogger.DisableLog("rdApp.*")

DATA = os.path.join(HERE, "..", "..", "data")
EXPERIMENT_ID = "AR_2023"

LNP = {
    "Cat-LNP": {"helper": "18:0 DDAB", "hl_corpus_name": "DDAB",
                "ratio": (35, 18, 2.5, 44.5), "peg": "C14PEG2000",
                "diam_nm": 106, "pdi": 0.20, "pka": 7.7, "zeta_mV": 1.81,
                "encap_pct": 94.7, "fig": "Fig 3a"},
    "Neu-LNP": {"helper": "DSPC", "hl_corpus_name": "DSPC",
                "ratio": (20, 50, 2.5, 27.5), "peg": "C14PEG2000",
                "diam_nm": 87, "pdi": 0.22, "pka": 5.8, "zeta_mV": 0.04,
                "encap_pct": 73.7, "fig": "Fig S15a"},
    "An-LNP": {"helper": "18:1 PA", "hl_corpus_name": None,
               "ratio": (60, 26, 2, 12), "peg": "C14PEG2000",
               "diam_nm": 150, "pdi": 0.30, "pka": 5.9, "zeta_mV": -18.6,
               "encap_pct": 53.1, "fig": "Fig S15c"},
}
RATIOS_FIG_S1 = [(60, 26, 2, 12), (50, 35, 2.5, 12.5),
                 (35, 18, 2.5, 44.5), (20, 50, 2.5, 27.5)]

ORGAN_TO_MODEL_TARGET = {"lung": "lung_epithelium", "liver": "liver",
                         "kidney": "kidney", "spleen": "spleen",
                         "heart": "heart"}

# Panels that are NOT delivery measurements, recorded so the zero is visible.
NON_DELIVERY = [
    ("Fig 1c", "hydrodynamic diameter of 144 individual LNPs and of the pool",
     "nm", "particle size -- not cargo expression"),
    ("Fig 1d,e", "diameter as a function of helper lipid", "nm",
     "particle size"),
    ("Fig 1f", "diameter as a function of helper-lipid mole ratio", "nm",
     "particle size"),
    ("Fig 2h", "normalized delivery vs hydrodynamic diameter", "nm vs %",
     "the x axis is particle size"),
    ("Fig 3a / S15a / S15c / S14a / S17a", "Diam., PDI, pKa, Zeta, Encap Eff.",
     "nm / unitless / pH / mV / %", "physicochemical characterisation"),
    ("Fig 3c,d", "tdTomato immunofluorescence micrographs of lung and liver",
     "image", "qualitative image, no quantified per-LNP value"),
    ("Fig 4a", "t-SNE of 27 lung cell clusters", "arbitrary",
     "transcriptome embedding -- the host response"),
    ("Fig 4b,c", "aVHH pseudogene reads / % aVHH+ cells per scRNA-seq cluster",
     "% of cells", "read out of the CARGO transcript by scRNA-seq at 1 h; "
     "it is cargo mRNA abundance, not cargo expression, and it is reported "
     "per cluster for one LNP without a per-formulation value table"),
    ("Fig 5a", "differentially expressed genes in lung ECs", "log2FC / p",
     "host transcriptional response -- the paper's title subject"),
    ("Fig 5b / S22", "Reactome pathway enrichment, 27 of 835 pathways",
     "z-score", "host transcriptional response"),
    ("Fig S5 / S19", "mouse body weight at 24, 48, 96 h", "g",
     "tolerability"),
    ("Fig S8-S11", "fold enrichment of helper-lipid properties", "fold",
     "a derived statistic over the pooled screen, not a measurement"),
    ("Fig S12", "TEM images", "image", "morphology"),
    ("Fig S16", "QUANT biodistribution of Cat-/Neu-/An-LNP", "normalized",
     "biodistribution of the PARTICLE, explicitly not functional delivery"),
    ("Fig S18", "40 serum cytokines and chemokines at 6 h", "pg/mL",
     "immune response"),
    ("Fig S20", "cluster marker dot plots", "expression", "transcriptome"),
    ("Fig S21", "top 10 upregulated genes vs each control", "log2FC",
     "host transcriptional response"),
    ("Fig S23", "representative FACS gating strategies", "gate",
     "method figure"),
]


def corpus_rows():
    out = []
    for f in ("data_vitro.csv", "data_vivo.csv"):
        with open(os.path.join(DATA, f), newline="") as fh:
            out += list(csv.DictReader(fh))
    return out


def anchor_vs_registry(qc):
    """C10-equivalent: the anchors this folder BUILT vs the pinned reagents.

    structure_check.py is not run here (it writes a tracked file other agents
    read); the comparison it would make is made locally instead.
    """
    sys.path.insert(0, os.path.join(HERE, "..", "..", "provenance"))
    from reagent_registry import REGISTRY                      # noqa: E402
    res = {}
    with open(os.path.join(HERE, "reference_structures.csv"),
              newline="") as fh:
        anchors = {r["lipid_name"]: r["reference_SMILES"]
                   for r in csv.DictReader(fh)}
    pairs = [("18:0 DDAB", "DDAB"), ("18:0 PC (DSPC)", "DSPC")]
    for anchor_name, reg_name in pairs:
        a, b = anchors[anchor_name], REGISTRY[reg_name][0]
        ma, mb = Chem.MolFromSmiles(a), Chem.MolFromSmiles(b)
        ka, kb = Chem.MolToInchiKey(ma), Chem.MolToInchiKey(mb)
        res[anchor_name] = {
            "registry_entry": reg_name,
            "anchor_inchikey": ka, "registry_inchikey": kb,
            "skeleton_match": ka.split("-")[0] == kb.split("-")[0],
            "full_match": ka == kb,
            "byte_identical_strings": a == b,
            "_note": "byte identity is NOT the independence signal (F-037); "
                     "make_reference.py's AST guard is.",
        }
    res["18:1 PA"] = {
        "registry_entry": None,
        "anchor_inchikey": Chem.MolToInchiKey(
            Chem.MolFromSmiles(anchors["18:1 PA"])),
        "_note": "18:1 PA (DOPA) is NOT in provenance/reagent_registry.py. "
                 "The registry's only PA is 14PA (DMPA, 14:0/14:0), a "
                 "different molecule. Nothing to compare against; the anchor "
                 "stands alone and An-LNP's helper lipid is therefore the one "
                 "excipient in this paper with no second attestation.",
    }
    res["cKK-E12"] = {
        "reference_SMILES": anchors["cKK-E12"],
        "_note": "EMPTY ON PURPOSE. The paper draws no structure for its only "
                 "ionizable lipid, so there is no anchor to build. The "
                 "absence is recorded rather than filled in.",
    }
    qc["C10_equivalent_anchor_comparison"] = res
    return res


def main():
    with open(os.path.join(HERE, "04_digitized.json")) as fh:
        dig = json.load(fh)
    qc = {"_what": "QC for PMID 36701517, every check in numbers"}

    # ---------- check: the four molar ratios sum to exactly 100 ----------
    sums = [sum(r) for r in RATIOS_FIG_S1]
    qc["figS1_molar_ratios"] = {
        "columns": RATIOS_FIG_S1, "sums": sums,
        "max_abs_residual_from_100": max(abs(s - 100) for s in sums),
        "n_exact": sum(1 for s in sums if abs(s - 100) < 1e-12),
        "_why": "a round number the authors implied but never printed",
    }

    # ---------- check: each lead LNP's pie labels pin a Fig S1 column ----
    pie = {}
    for name, d in LNP.items():
        il, chol, peg, hl = d["ratio"]
        printed = {il, chol, hl}          # the PEG slice carries no label
        hits = [i for i, r in enumerate(RATIOS_FIG_S1)
                if {r[0], r[1], r[3]} == printed]
        pie[name] = {"printed_pie_labels": sorted(printed),
                     "matching_figS1_columns": hits,
                     "unique": len(hits) == 1,
                     "peg_recovered": RATIOS_FIG_S1[hits[0]][2] if
                     len(hits) == 1 else None,
                     "peg_agrees_with_stored": (
                         len(hits) == 1
                         and RATIOS_FIG_S1[hits[0]][2] == peg)}
    qc["lead_LNP_pie_labels_pin_a_figS1_column"] = {
        "per_lnp": pie,
        "n_unique": sum(1 for v in pie.values() if v["unique"]),
        "_why": "the PEG mol% is the unlabelled slice of each pie; matching "
                "the three printed labels against the Fig S1 columns recovers "
                "it uniquely for all three lead LNPs, 3/3.",
    }

    # ---------- check: Fig S6c 'average' column is the mean of three ------
    # Transcribed by eye at full printed precision from the Fig S6c raster
    # (SI p.17, xref 202, 220.0 dpi). The mean identity is what validates the
    # transcription: a single misread digit in ANY of the four columns breaks
    # it, and it holds on 22 of 22 rows.
    s6c = json.load(open(os.path.join(HERE, "04b_figS6c_cationic.json")))
    errs = [abs((s6c["lung"][i] + s6c["liver"][i] + s6c["spleen"][i]) / 3.0
                - s6c["average"][i]) for i in range(len(s6c["average"]))]
    qc["figS6c_average_is_the_mean_of_the_three_organ_columns"] = {
        "n_rows": len(errs), "max_abs_error": max(errs),
        "n_within_1e-8": sum(1 for e in errs if e < 1e-8),
        "_why": "reproduces a column the authors derived but never defined, "
                "and simultaneously validates 88 transcribed numerals.",
    }
    pool = {}
    for col in ("lung", "liver", "spleen", "average"):
        tot = sum(s6c[col]) + s6c["naked"][col]
        pool[col] = {"sum_including_naked_barcode": tot,
                     "residual_vs_100": tot - 100.0}
    pool["_interpretation"] = (
        "Each column is a SHARE OF ITS OWN POOL and the shares sum to 100%. "
        "Fig S3a states it as a worked example: 3, 1 and 1 of five sequenced "
        "barcodes give 60, 20 and 20, which sums to 100.000 with residual "
        "0.000 and reproduces as 3/5, 1/5, 1/5 exactly. The residuals here "
        "are all POSITIVE and are attributable to exactly one pooled LNP "
        "missing from the table: the text says 23 cationic LNPs were pooled "
        "and the table lists 22. The three organ residuals average to "
        "6.5528, which equals the 'average' column residual to 3e-4 -- the "
        "same arithmetic identity the column itself obeys.")
    pool["mean_of_three_organ_residuals"] = (
        pool["lung"]["residual_vs_100"] + pool["liver"]["residual_vs_100"]
        + pool["spleen"]["residual_vs_100"]) / 3.0
    pool["average_column_residual"] = pool["average"]["residual_vs_100"]
    pool["identity_residual"] = abs(pool["mean_of_three_organ_residuals"]
                                    - pool["average_column_residual"])
    qc["figS6c_pool_shares_sum_to_100"] = pool
    qc["figS3a_worked_example"] = {
        "counts": [3, 1, 1], "n_sequenced": 5,
        "printed_normalized_delivery_pct": [60, 20, 20],
        "reconstructed_from_counts": [100 * 3 / 5, 100 * 1 / 5, 100 * 1 / 5],
        "max_abs_residual": 0.0,
        "sum": 100.0, "sum_residual_vs_100": 0.0,
    }

    # ---------- geometry QC from the digitiser ---------------------------
    geo = {"native_dpi_measured": dig["native_dpi_measured"], "panels": {}}
    worst_loo, worst_spine = 0.0, 0.0
    for pn, p in dig["panels"].items():
        for org, o in p["organs"].items():
            ax = o["axis"]
            worst_loo = max(worst_loo, ax["loo_max_px"])
            worst_spine = max(worst_spine, ax["spine_top_residual_px"])
            geo["panels"][f"{pn}/{org}"] = {
                "n_y_ticks": ax["n_y_ticks"],
                "loo_max_px": ax["loo_max_px"],
                "px_per_percent": ax["px_per_percent"],
                "spine_top_residual_px": ax["spine_top_residual_px"],
                "tick_fit_residual_px":
                    ax["calibration"]["tick_fit_residual_px"]}
    geo["worst_loo_px"] = worst_loo
    geo["gate_px_unrescaled"] = 1.5
    geo["gate_px_rescaled_for_220dpi"] = round(1.5 * 600 / 220.0, 3)
    geo["gate_actually_used"] = "the UNRESCALED 1.5 px; the rescale was " \
                                "available and never needed"
    geo["worst_spine_top_residual_px"] = worst_spine
    geo["spine_check_what_it_is"] = (
        "every y-axis spine is stroked from the 100% level to the 0% level, "
        "so the row the tick-fitted scale puts 100% at must be the spine's "
        "top row. The spine takes no part in the fit. 21/21 sub-panels pass "
        "at <= 0.71 px, which also pins the value of the top stub in the one "
        "panel (Fig S13b Spleen) that draws only five.")
    qc["geometry"] = geo

    # ---------- the cross-panel check ------------------------------------
    # Fig S13 plots tdTom+ RAW with its PBS control; Fig S15 plots the same
    # experiment ABOVE BACKGROUND. So S13_raw - S13_pbs should equal S15.
    # Joined on the LABEL: the drawn order differs between the two figures.
    cross = {"pairs": [], "_join": "on (organ, cell type) LABEL, never on "
                                   "column index -- Fig S13 lung is "
                                   "ECs/Monocytes/DCs/T/B and Fig S15 lung is "
                                   "ECs/DCs/T/B/Monocytes"}
    for raw_panel, abg_panel, lnp in (("S13b", "S15b", "Neu-LNP"),
                                      ("S13c", "S15d", "An-LNP")):
        for org in ("lung", "liver", "heart"):
            if org not in dig["panels"][abg_panel]["organs"]:
                continue
            rawo = dig["panels"][raw_panel]["organs"][org]["values_pct"]
            abgo = dig["panels"][abg_panel]["organs"][org]["values_pct"]
            for lab in abgo:
                r, a = rawo.get(lab, {}).get("lnp"), abgo[lab]["lnp"]
                if r is None or a is None:
                    continue
                pbs = rawo.get(lab, {}).get("pbs")
                cross["pairs"].append({
                    "lnp": lnp, "organ": org, "cell_type": lab,
                    "raw_figS13_pct": r, "above_background_figS15_pct": a,
                    "pbs_figS13_pct": pbs,
                    "raw_minus_abovebackground": round(r - a, 3),
                    "ratio": round(r / a, 4) if a else None})
    d = [p["raw_minus_abovebackground"] for p in cross["pairs"]]
    cross["n_pairs"] = len(d)
    cross["raw_minus_abovebackground_median"] = \
        sorted(d)[len(d) // 2] if d else None
    cross["raw_minus_abovebackground_min"] = min(d) if d else None
    cross["raw_minus_abovebackground_max"] = max(d) if d else None
    cross["n_within_PBS_band_0_to_5pct"] = sum(1 for x in d if -0.6 <= x <= 5)
    cross["_expected"] = ("the difference should be the PBS level, which Fig "
                          "S13 shows at 0-4.5% everywhere. A difference "
                          "outside roughly [-0.6, 5] would mean the two "
                          "panels do not plot the same measurement.")
    qc["cross_panel_figS13_vs_figS15"] = cross

    # ---------- values table ---------------------------------------------
    rows05 = []

    def add(**kw):
        rows05.append(kw)

    # (1) the individually characterised LNPs, Fig S13 -- DELIVERY
    emit_items = []
    for pn, lnp in (("S13a", "Cat-LNP"), ("S13b", "Neu-LNP"),
                    ("S13c", "An-LNP")):
        for org, o in dig["panels"][pn]["organs"].items():
            for lab, v in o["values_pct"].items():
                dgn = o["diagnostics"][lab]
                val = v["lnp"]
                if val is None:
                    reason = ("no -- not readable: "
                              + dgn["lnp"].get("refused", "")[:150])
                    adm = "no"
                else:
                    adm = "yes"
                    reason = ("yes -- Cre mRNA translated into functional Cre "
                              "protein, read as genomic recombination "
                              "(tdTomato+) in an individually formulated, "
                              "fully specified LNP; this is cargo expression")
                add(value_id=f"{lnp}|{org}|{lab}", source_panel=f"Fig {pn}",
                    source_page="SI p.24-25", lnp=lnp,
                    readout_as_printed="tdTom+ Cells (%)",
                    units="percent of cells",
                    value=val, pbs_control_pct=v["pbs"],
                    cargo="Cre mRNA", model="Ai14 mouse",
                    organ=org, cell_type=lab, route="intravenous",
                    dose="1.3 mg/kg", timepoint="4 days",
                    assay="flow cytometry, tdTomato+ gate vs PBS-treated mice",
                    pool_identity="not pooled -- single LNP per animal",
                    corpus_admissible=adm, reason=reason)
                if adm == "yes":
                    emit_items.append((lnp, org, lab, val))

    # (2) the same measurements above background, Fig S15 -- NOT emitted
    for pn, lnp in (("S15b", "Neu-LNP"), ("S15d", "An-LNP")):
        for org, o in dig["panels"][pn]["organs"].items():
            for lab, v in o["values_pct"].items():
                val = v["lnp"]
                add(value_id=f"{lnp}|{org}|{lab}|abovebackground",
                    source_panel=f"Fig {pn}", source_page="SI p.27", lnp=lnp,
                    readout_as_printed="tdTom+ Cells (%), above background",
                    units="percent of cells", value=val,
                    pbs_control_pct=None, cargo="Cre mRNA",
                    model="Ai14 mouse", organ=org, cell_type=lab,
                    route="intravenous", dose="1.3 mg/kg", timepoint="4 days",
                    assay="flow cytometry, PBS level subtracted",
                    pool_identity="not pooled",
                    corpus_admissible="no",
                    reason="no -- SAME measurements as Fig S13 in a different "
                           "normalisation (background-subtracted). Emitting "
                           "both would double-count, and mixing two "
                           "normalisations inside one Experiment_ID is error "
                           "class 22. Used here only to cross-check Fig S13.")

    # (3) the pooled barcode screen -- Fig S6 tables
    add(value_id="figS6abc|pooled normalized delivery|137 LNPs",
        source_panel="Fig S6a (neutral), S6b (anionic), S6c (cationic)",
        source_page="SI p.15, 16, 17", lnp="137 pooled LNPs",
        readout_as_printed="average normalized delivery (%) and per-organ "
                           "averages", units="percent share of the pool",
        value=None, pbs_control_pct=None, cargo="Cre mRNA + DNA barcode",
        model="Ai14 mouse", organ="lung/liver/spleen", cell_type="19 types",
        route="intravenous", dose="1.0 mg/kg total nucleic acid "
                                  "(~0.043 mg/kg per LNP)",
        timepoint="4 days",
        assay="deep sequencing of DNA barcodes recovered from tdTomato+ cells",
        pool_identity="three pools, one per helper-lipid charge class: "
                      "23 cationic, 59 neutral, 55 anionic (as stated in the "
                      "text)",
        corpus_admissible="no",
        reason="no -- POOL-RELATIVE. Fig S3a's worked example defines the "
               "value as a share of the reads recovered from one cell type, "
               "and the shares across a pool sum to 100%, reproduced here at "
               "residual 0.000 on the worked example and attributable to one "
               "missing table row on the real cationic data. Each value "
               "therefore depends on the other pool members, which no corpus "
               "column records. Transcription of all 137 rows is also NOT "
               "done: only the 22-row cationic table was transcribed and "
               "validated. See DECISIONS.md D-2.")

    # (4) everything that is not a delivery measurement
    for fig, what, units, why in NON_DELIVERY:
        add(value_id=f"{fig}|{what}", source_panel=fig, source_page="",
            lnp="", readout_as_printed=what, units=units, value=None,
            pbs_control_pct=None, cargo="", model="", organ="", cell_type="",
            route="", dose="", timepoint="", assay="", pool_identity="",
            corpus_admissible="no",
            reason="no -- NOT A DELIVERY ENDPOINT: " + why)

    # (5) aVHH protein, Fig S17 -- delivery in kind, not digitised
    add(value_id="figS17d-h|aVHH+ cells and MFI",
        source_panel="Fig S17c-h", source_page="SI p.29-30",
        lnp="Cat-LNP, Neu-LNP, An-LNP",
        readout_as_printed="% aVHH+ cells; aVHH MFI", units="percent; MFI",
        value=None, pbs_control_pct=None, cargo="aVHH mRNA",
        model="C57BL/6 mouse", organ="lung/liver/heart/spleen/kidney",
        cell_type="multiple", route="intravenous (tail vein)",
        dose="1.3 mg/kg", timepoint="1 day",
        assay="flow cytometry for the anchored VHH protein",
        pool_identity="not pooled",
        corpus_admissible="no",
        reason="no -- a genuine cargo-expression endpoint (protein "
               "abundance) and admissible IN KIND, but NOT EXTRACTED: it is "
               "a second readout of the same three LNPs at a different dose, "
               "cargo and timepoint, so emitting it alongside Fig S13 would "
               "put two assays in one Experiment_ID. Left for a later pass; "
               "see DECISIONS.md D-3.")

    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        cols = ["value_id", "source_panel", "source_page", "lnp",
                "readout_as_printed", "units", "value", "pbs_control_pct",
                "cargo", "model", "organ", "cell_type", "route", "dose",
                "timepoint", "assay", "pool_identity", "corpus_admissible",
                "reason"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows05:
            w.writerow({c: r.get(c, "") for c in cols})

    counts = {"total": len(rows05),
              "admissible_yes": sum(1 for r in rows05
                                    if r["corpus_admissible"] == "yes"),
              "admissible_no": sum(1 for r in rows05
                                   if r["corpus_admissible"] == "no")}
    counts["sums"] = counts["admissible_yes"] + counts["admissible_no"] \
        == counts["total"]
    qc["values_table"] = counts

    # ---------- 06_rows.csv ----------------------------------------------
    with open(os.path.join(DATA, "data_vitro.csv"), newline="") as fh:
        schema = csv.DictReader(fh).fieldnames
    comp, amb = load(corpus_rows())
    ckk = smi(comp, amb, "CKK-E12")
    with open(os.path.join(HERE, "reference_structures.csv"),
              newline="") as fh:
        anchors = {r["lipid_name"]: r["reference_SMILES"]
                   for r in csv.DictReader(fh)}
    hl_smiles = {"Cat-LNP": anchors["18:0 DDAB"],
                 "Neu-LNP": anchors["18:0 PC (DSPC)"],
                 "An-LNP": anchors["18:1 PA"]}

    # z of log10, within the experiment, population sd (CLAUDE.md section 3)
    vals = [v for _, _, _, v in emit_items]
    logs = [math.log10(v) for v in vals]
    mu = sum(logs) / len(logs)
    sd = (sum((x - mu) ** 2 for x in logs) / len(logs)) ** 0.5
    zs = [(x - mu) / sd for x in logs]

    records = {}
    for (lnp, org, lab, v) in emit_items:
        records[f"{lnp}|{org}|{lab}"] = {
            "lnp": lnp, "organ": org, "cell_type": lab, "value_pct": v,
            "helper": LNP[lnp]["helper"],
            "stereo_unspecified": True}

    def conf(item_id, rec):
        # PER-ITEM, derived from the record. A literal here must raise (F-011).
        return "medium"if rec["stereo_unspecified"] else "high"

    def note(item_id, rec):
        return ("cell type " + rec["cell_type"] + "; " + rec["lnp"]
                + " (helper lipid " + rec["helper"] + "). cKK-E12 is the "
                "only ionizable lipid in this paper and the paper draws no "
                "structure for it; the stored SMILES is the corpus string, "
                "which carries no stereocentres although cKK-E12 has two. "
                "Value digitised from Fig S13 at a measured 220 dpi.")

    flags = emit.emit_rows(
        items=list(records), records=records,
        confidence=conf, note=note,
        structure_assumed=lambda i, r: "1",
        value_normalization_unit=lambda i, r: "experiment",
        label="AR_2023")

    out_rows = []
    for (item, fl), z in zip(zip(list(records), flags), zs):
        rec = records[item]
        lnp = rec["lnp"]
        il, chol, peg, hl = LNP[lnp]["ratio"]
        row = {c: "" for c in schema}
        row.update({
            "Experiment_ID": EXPERIMENT_ID,
            "Formulation_ID": lnp,
            "IL_name": "cKK-E12", "IL_SMILES": ckk,
            "IL_molratio": il,
            "IL_to_nucleicacid_massratio": 10,
            "HL_name": LNP[lnp]["helper"], "HL_SMILES": hl_smiles[lnp],
            "HL_molratio": hl,
            "CHL_name": "Cholesterol", "CHL_molratio": chol,
            "PEG_name": LNP[lnp]["peg"], "PEG_molratio": peg,
            "Aqueous_buffer": "10 mM citrate",
            "Mixing_method": "microfluidic, 3:1 aqueous:ethanol",
            "Model": "in_vivo", "Model_type": "Mouse_Ai14",
            "Model_target": ORGAN_TO_MODEL_TARGET[rec["organ"]],
            "Route_of_administration": "intravenous",
            # Cargo_type: the corpus vocabulary is FLuc / FVII / GFP / RLuc /
            # base_editor / hEPO / peptide_barcode. Cre recombinase is none of
            # them, so "Cre" is a NEW value and is declared as one in
            # DECISIONS.md D-5 rather than forced into "base_editor", which
            # SX_2025 uses for a genuinely different chemistry.
            "Cargo": "mRNA", "Cargo_type": "Cre",
            "Dose_ug_nucleicacid": "",
            "Experiment_method": "editing_efficiency_normalized",
            # Experiment_batching is a two-value vocabulary (barcoded /
            # individual). The CELL TYPE has no column in this schema, so it
            # lives in curation_note and in 05_values.csv; see DECISIONS.md
            # D-4, which is why these rows share a duplicate
            # (Experiment_ID, IL_name, Model_target, Model_type) key.
            "Experiment_batching": "individual",
            "Experiment_value": round(z, 6),
            "Publication_DOI": "10.1021/acs.nanolett.2c04479",
            "Publication_PMID": 36701517,
            "Publication_link":
                "https://doi.org/10.1021/acs.nanolett.2c04479",
        })
        row.update(fl)
        out_rows.append(row)

    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=schema)
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    qc["rows_emitted"] = {
        "n": len(out_rows),
        "n_unique_IL_SMILES": 1,
        "n_formulations": len({r["Formulation_ID"] for r in out_rows}),
        "rows_with_a_caveat": emit.caveat_count(flags),
        "z_mean": round(sum(zs) / len(zs), 12),
        "z_sd_population": round(
            (sum(x * x for x in zs) / len(zs)) ** 0.5, 12),
        "value_pct_min": min(vals), "value_pct_max": max(vals),
    }

    anchor_vs_registry(qc)

    # ---------- the three numbers ----------------------------------------
    qc["the_three_numbers"] = {
        "distinct_ionizable_lipid_structures": 1,
        "distinct_helper_lipid_structures": 18,
        "distinct_formulations_designed": 144,
        "distinct_formulations_injected": 137,
        "formulation_design": "18 helper lipids x 2 PEG-lipids "
                              "(C14PEG2K, C18PEG2K) x 4 molar ratios = 144",
        "individually_measured_LNPs": 4,
        "individually_measured_names": ["Cat-LNP", "Neu-LNP", "An-LNP",
                                        "Kauffman LNP (literature control, "
                                        "Fig S14)"],
        "what_the_design_can_support":
            "ONE ionizable-lipid structure. Every one of the 144 "
            "formulations uses cKK-E12; what varies is the helper lipid, the "
            "PEG-lipid and the molar ratio. A corpus keyed on ionizable-lipid "
            "structure gains one lipid from this paper, however many rows are "
            "emitted.",
        "second_PEG_lipid_identity":
            "C18PEG2K. Named NOWHERE in the article text or the SI Methods "
            "('PEGs, cholesterol, and helper lipids ... purchased from Avanti "
            "Lipids'), and recoverable only from the PEG column of the Fig S6 "
            "tables.",
    }

    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1, sort_keys=True, default=str)
    print("wrote 05_values.csv (%d rows), 06_rows.csv (%d rows), 07_qc.json"
          % (len(rows05), len(out_rows)))
    print("admissible yes/no:", counts)


if __name__ == "__main__":
    main()
