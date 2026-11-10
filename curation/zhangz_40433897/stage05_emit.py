#!/usr/bin/env python3
"""Stage 05 - 05_values.csv (every value), 06_rows.csv (the admissible ones),
01_blocks.csv (the rounds) and 03_enumerated.csv (the chemistry).

ADMISSIBILITY. data/'s Experiment_method vocabulary is delivery-only:
luminescence_normalized, luminescence_discretized_normalized and
LRP6_knockdown_normalized. A value is corpus_admissible only if it is a
nucleic-acid DELIVERY readout AND it belongs to a formulation whose ionizable
lipid has a structure. Particle size, encapsulation efficiency, haemolysis,
uptake, STING activation, body weight and tumour volume are not delivery
endpoints and are marked no, with the reason, rather than left out.

Quality flags are written ON THE ROW through curation/lib/emit.emit_rows with
PER-ITEM CALLABLES. A hard-coded literal raises there (F-011/F-025/F-046).
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import emit                                   # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

BARS = json.load(open(os.path.join(HERE, "02_bars.json")))["panels"]
STRUCT = json.load(open(os.path.join(HERE, "03_proposed_structures.json")))
PROP = {r["name"]: r for r in STRUCT["lipids"]}

# ---------------------------------------------------------------------------
# The rounds. "Iterative screening" means rounds, and rounds are different
# experiments: a different library, a different base formulation, a different
# cell line or a different cargo is a different standardisation block.
# ---------------------------------------------------------------------------
BLOCKS = [
    dict(block="R1_formulation_DoE", panels=["FigS4b"],
         what="60-prescription orthogonal formulation design on ONE ionizable "
              "lipid (VN2C2), siLuc knockdown in 4T1-Luc",
         n_distinct_structures=1, n_formulations=60, n_individually_measured=60,
         cargo="siLuc", cell="4T1-Luc",
         readout="Luciferase expression (%)",
         formulation_source="Supplementary Table 2 (PDF text layer, tier 2)"),
    dict(block="R2_structure_x_PEG_4T1", panels=["FigS5a"],
         what="6 ionizable-lipid structures x 3 DMG-PEG levels, siLuc "
              "knockdown in 4T1-Luc",
         n_distinct_structures=6, n_formulations=18, n_individually_measured=18,
         cargo="siLuc", cell="4T1-Luc",
         readout="Luciferase expression (%)",
         formulation_source="NOT RECOVERABLE - the three PEG levels are "
                            "labelled PEG-0.5 / PEG-1 / PEG-2 and no mol% for "
                            "any component is given for this figure"),
    dict(block="R3_structure_x_PEG_MC38", panels=["FigS5b"],
         what="the same 6 structures x 3 DMG-PEG levels, siLuc knockdown in "
              "MC38-Luc",
         n_distinct_structures=6, n_formulations=18, n_individually_measured=18,
         cargo="siLuc", cell="MC38-Luc",
         readout="Relative luciferase expression (%)",
         formulation_source="NOT RECOVERABLE, as R2"),
    dict(block="R4_lipid_to_pDNA_ratio", panels=["FigS7"],
         what="pGFP transfection vs lipid:pDNA mass ratio",
         n_distinct_structures=0, n_formulations=5, n_individually_measured=7,
         cargo="pGFP", cell=None,
         readout="Percentage of GFP+ cells (%)",
         formulation_source="NOT RECOVERABLE - the caption names no lipid and "
                            "no component mol%"),
    dict(block="R5_VE_TPGS_substitution", panels=["FigS8a", "FigS8b"],
         what="mFLuc delivery when VN2C1 is replaced by, or SM-102 is "
              "supplemented with, free vitamin E or TPGS",
         n_distinct_structures=2, n_formulations=8, n_individually_measured=10,
         cargo="mFLuc", cell=None,
         readout="Luciferase intensity (RLU)",
         formulation_source="NOT RECOVERABLE - no mol% stated"),
    dict(block="R6_storage", panels=["FigS13b"],
         what="mGFP delivery by the VN2C1 LNP before and after 7 days storage",
         n_distinct_structures=1, n_formulations=1, n_individually_measured=3,
         cargo="mGFP", cell="RAW264.7",
         readout="GFP+ cells (%)",
         formulation_source="NOT RECOVERABLE"),
]

CONTROLS = {"Lip3000", "Naked siRNA", "Naked pGFP", "PBS", "VE", "TPGS"}

NON_DELIVERY = [
    ("FigS3a", "Eff. Diam. (nm), 60 prescriptions", "particle size"),
    ("FigS3b", "gel retardation photographs", "no numeric scale; and gel "
     "retardation is a complexation readout, not delivery"),
    ("FigS4a", "four fitted GFP+% response surfaces", "a fitted surface over "
     "factor levels; no individual formulation is reported on it"),
    ("FigS6", "Eff. Diam. (nm) in DMEM/FBS", "colloidal stability"),
    ("FigS9", "p-STING western blot; type-I IFN transcript", "immune-pathway "
     "activation"),
    ("FigS10", "Encapsulation rate (%)", "encapsulation efficiency is a "
     "formulation property, not delivery"),
    ("FigS11", "haemolysis tubes at pH 7.4 / 5.5", "photograph, no numeric "
     "scale; membrane lysis is not delivery"),
    ("FigS12", "TEM micrographs", "morphology"),
    ("FigS13a", "particle size and PDI over 7 days", "particle size"),
    ("FigS14", "Cy5-mRNA uptake at 1/2/4 h", "UPTAKE is particle association, "
     "not cargo expression"),
    ("FigS15", "body weight", "toxicology"),
    ("FigS16", "tumour growth / survival", "therapeutic efficacy"),
    ("FigS17", "PD-L1 positive rate and MFI", "receptor expression"),
]


def lipid_of(block, series):
    if block == "R1_formulation_DoE":
        return "VN2C2" if series.startswith("LNP") else None
    if block in ("R2_structure_x_PEG_4T1", "R3_structure_x_PEG_MC38"):
        return series
    if block == "R5_VE_TPGS_substitution":
        if series == "VN2C1":
            return "VN2C1"
        if series in ("SM102", "SM-102"):
            return "SM-102"
        if series in ("VE+SM102", "TPGS+SM102"):
            return "SM-102"
        return None
    if block == "R6_storage":
        return "VN2C1" if series.startswith("VN2C1") else None
    return None


def admissible(block, series, lipid, value):
    if value is None:
        return "no", "the bar was REFUSED by curation/lib/bars.bar_top; see 07_qc.json"
    if series in CONTROLS:
        return "no", ("%s is a control arm, not an LNP of a library ionizable "
                      "lipid" % series)
    if lipid is None:
        return "no", "no ionizable lipid can be attached to this bar"
    rec = PROP.get(lipid)
    if lipid == "SM-102":
        return "no", ("SM-102 is a benchmark here and its formulation mol% is "
                      "not stated anywhere in the SI, so the row would carry a "
                      "structure and a value but no formulation")
    if rec is None or not rec.get("proposed_IL_SMILES"):
        return "no", ("%s has NO structure: %s" % (lipid, (rec or {}).get(
            "why", "not characterised in the SI")))
    return "yes", ("delivery endpoint, and %s has a single "
                   "connectivity-determinate core among the 34 drawn amines "
                   "(PROPOSED, not pinned - see DECISIONS.md D-1)" % lipid)


def main():
    # ---------------- 01_blocks.csv ----------------
    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(BLOCKS[0]))
        w.writeheader()
        for b in BLOCKS:
            r = dict(b)
            r["panels"] = ";".join(b["panels"])
            w.writerow(r)

    # ---------------- 03_enumerated.csv ----------------
    cols = ["IL_name", "IL_SMILES", "hrms_panel", "formula_ion_as_printed",
            "core_formula_derived", "core_CAS_in_FigS1a", "structure_pinned",
            "proposed_M_plus_H", "printed_calc_M_plus_H",
            "delta_vs_printed_mDa", "why"]
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in STRUCT["lipids"]:
            w.writerow({
                "IL_name": r["name"], "IL_SMILES": r.get("proposed_IL_SMILES", ""),
                "hrms_panel": r["panel"],
                "formula_ion_as_printed": r["formula_ion"],
                "core_formula_derived": r["core_formula"],
                "core_CAS_in_FigS1a": r["core_cas"] or "",
                "structure_pinned": "no",
                "proposed_M_plus_H": r.get("proposed_M_plus_H", ""),
                "printed_calc_M_plus_H": r.get("printed_calc_M_plus_H", ""),
                "delta_vs_printed_mDa": r.get("delta_vs_printed_mDa", ""),
                "why": r["why"]})

    # ---------------- 05_values.csv ----------------
    vcols = ["block", "panel", "source_location", "series", "IL_name",
             "readout_as_printed", "units", "value", "cargo", "cell_line",
             "organ", "route", "dose", "timepoint", "assay",
             "tier", "detector_status", "nrow_relaxed",
             "corpus_admissible", "admissibility_reason"]
    values = []
    for B in BLOCKS:
        for pid in B["panels"]:
            P = BARS[pid]
            for v in P["values"]:
                lip = lipid_of(B["block"], v["series"])
                adm, why = admissible(B["block"], v["series"], lip, v.get("value"))
                values.append({
                    "block": B["block"], "panel": pid,
                    "source_location": "%s, bar %d (x %d-%d)"
                                       % (P["values"][0].get("figure", pid),
                                          v["index"] + 1, v["x0"], v["x1"])
                    if False else "%s bar %d, embedded raster page %s"
                                  % (pid, v["index"] + 1,
                                     {"FigS4b": 6, "FigS5a": 7, "FigS5b": 7,
                                      "FigS7": 7, "FigS8a": 8, "FigS8b": 8,
                                      "FigS13b": 9}[pid]),
                    "series": v["series"], "IL_name": lip or "",
                    "readout_as_printed": B["readout"],
                    "units": "percent" if "%" in B["readout"] else "RLU",
                    "value": v.get("value", ""),
                    "cargo": B["cargo"], "cell_line": B["cell"] or "NOT STATED",
                    "organ": "in_vitro", "route": "in_vitro",
                    "dose": "NOT STATED IN THE SI",
                    "timepoint": "NOT STATED IN THE SI",
                    "assay": ("firefly luciferase, " +
                              ("knockdown" if B["cargo"] == "siLuc"
                               else "expression"))
                             if "uciferase" in B["readout"]
                             else "flow cytometry, GFP+ fraction",
                    "tier": "5_pixels",
                    "detector_status": v["status"],
                    "nrow_relaxed": v.get("nrow_relaxed") or "",
                    "corpus_admissible": adm,
                    "admissibility_reason": why})
    for pid, what, why in NON_DELIVERY:
        values.append({c: "" for c in vcols} | {
            "block": "not_a_delivery_endpoint", "panel": pid,
            "source_location": pid, "readout_as_printed": what,
            "tier": "not_extracted", "corpus_admissible": "no",
            "admissibility_reason": why})
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=vcols)
        w.writeheader()
        w.writerows(values)

    # ---------------- 06_rows.csv ----------------
    adm = [v for v in values if v["corpus_admissible"] == "yes"]
    items, records = [], {}
    for i, v in enumerate(adm):
        iid = "%s|%s|%s" % (v["block"], v["panel"], v["series"])
        items.append(iid)
        records[iid] = {
            "IL_name": v["IL_name"],
            "structure_record": PROP[v["IL_name"]],
            "nrow_relaxed": v["nrow_relaxed"],
            "block": v["block"],
            "panel": v["panel"],
            "dpi": {"FigS4b": "219.9", "FigS5a": "219.5", "FigS5b": "219.5",
                    "FigS7": "219.6", "FigS8a": "216.2", "FigS8b": "216.2",
                    "FigS13b": "190.0"}[v["panel"]],
        }

    def confidence(iid, rec):
        # derived from the item's OWN record, never a literal
        if not rec["structure_record"]["structure_pinned"]:
            return "low"
        return "medium" if rec["nrow_relaxed"] else "high"

    def note(iid, rec):
        """One sentence, built from THIS item's own enumeration record.

        The per-lipid reason comes out of 03_proposed_structures.json's `why`
        field, not from a shared sentence: VN2C2 and VN3C1 are the ONLY drawn
        amine with their formula, while VN3C2 has two and is selected by the
        CH2 ladder instead. A single shared sentence would have asserted
        uniqueness for VN3C2, which is false.
        """
        s = rec["structure_record"]
        why = s["why"].rstrip(".")
        base = ("Structure PROPOSED, not pinned: the paper draws no VNxCyOz "
                "structure anywhere, so this is inferred from the Fig S2 HRMS "
                "formula (%s, core %s) plus the CAS list drawn in Fig S1a - "
                "%s. Our build reproduces the printed [M+H]+ calc to %s mDa. "
                "Value digitised from the %s embedded raster at %s dpi "
                "(tier 5: no Source Data, no numeric value table, and zero "
                "vector paths on any figure page)."
                % (s["formula_ion"], s["core_formula"], why,
                   s.get("delta_vs_printed_mDa"), rec["panel"], rec["dpi"]))
        if rec["nrow_relaxed"]:
            base += (" Bar top accepted at a %s-row rule, that bar's own "
                     "measured rule thickness." % rec["nrow_relaxed"])
        return base

    def structure_assumed(iid, rec):
        return 1

    def value_norm(iid, rec):
        return "experiment"

    flags = emit.emit_rows(items, records,
                           confidence=confidence, note=note,
                           structure_assumed=structure_assumed,
                           value_normalization_unit=value_norm,
                           label="zhangz_40433897")

    rcols = (["item_id", "block", "panel", "IL_name", "IL_SMILES",
              "Experiment_method", "Cargo", "Cargo_type", "Model_type",
              "Experiment_value_raw", "units"] + emit.columns())
    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rcols)
        w.writeheader()
        for iid, v, f in zip(items, adm, flags):
            row = {"item_id": iid, "block": v["block"], "panel": v["panel"],
                   "IL_name": v["IL_name"],
                   "IL_SMILES": PROP[v["IL_name"]]["proposed_IL_SMILES"],
                   "Experiment_method": "siLuc_knockdown_normalized",
                   "Cargo": "siRNA", "Cargo_type": "FLuc",
                   "Model_type": v["cell_line"],
                   "Experiment_value_raw": v["value"], "units": v["units"]}
            row.update(f)
            w.writerow(row)

    counts = {
        "bars_digitised_total": sum(len(BARS[p]["values"]) for p in BARS),
        "bars_read": sum(1 for p in BARS.values()
                         for v in p["values"] if v["status"] == "ok"),
        "bars_refused": sum(p["n_refused"] for p in BARS.values()),
        "values_rows_in_05": len(values),
        "values_from_panels": len(values) - len(NON_DELIVERY),
        "non_delivery_panels_listed": len(NON_DELIVERY),
        "corpus_admissible_yes": len(adm),
        "corpus_admissible_no": len(values) - len(adm),
        "rows_in_06": len(items),
        "rows_with_a_caveat": emit.caveat_count(flags),
    }
    counts["sums"] = (counts["corpus_admissible_yes"] +
                      counts["corpus_admissible_no"] == counts["values_rows_in_05"]
                      and counts["bars_read"] + counts["bars_refused"] ==
                      counts["bars_digitised_total"])
    json.dump(counts, open(os.path.join(HERE, "05_counts.json"), "w"), indent=1)
    print(json.dumps(counts, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
