#!/usr/bin/env python3
"""Read the per-replicate delivery values out of the publisher Source Data.

    python3 -I stage05_values.py     ->  05_values.csv

TIER 1 THROUGHOUT. Every number here is a cell in
`source/41467_2026_69077_MOESM6_ESM.xlsx`. No panel is digitised; no pixel is
read. The exact sheet and cell range of each value is carried on its own row.

FOUR SEPARATE ASSAY BLOCKS, NOT ONE SCREEN. Fig 2c/4c are flow cytometry and
Fig 2d/4d are fluorescence microscopy on the same two reporter lines, and the
same compound disagrees across the pair by ~2x in both directions (CBB14:
58.21% flow vs 31.95% microscopy; CBB11: 28.61% flow vs 53.30% microscopy).
They get four distinct block ids and must not be pooled.

CARGO IS PROTEIN / RNP, NOT mRNA. There is no mRNA anywhere in this paper.
"""
import csv
import json
import os
import statistics
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "..", "lib"))
sys.path.insert(0, _HERE)
import emit                                                   # noqa: E402
import openpyxl                                               # noqa: E402

XLSX = os.path.join(_HERE, "source", "41467_2026_69077_MOESM6_ESM.xlsx")

CONTROLS = {"LF", "C", "NE", "A"}

# panel -> everything that is a property of the panel, not of a compound.
# Transcribed from the figure legends (JATS XML of PMC12992809) and the Methods
# sections "Delivery of Cre recombinase in vitro" / "Delivery of ABE in vitro".
PANELS = {
 "Fig2c": dict(
   sheet="Figure 2", header_row=1, first_data_row=2, n_rep=4,
   cells="C1:AA5",
   readout_as_printed="RFP conversion (%)",
   method="flow cytometry (Novocyte Quanteon, Agilent; FITC 530/30, PE 586/20)",
   cargo="Cre recombinase (protein)",
   cargo_detail="purified Cre recombinase enzyme; no nucleic acid cargo",
   cell_line="CS cells = HEK293-loxP-GFP-RFP colour-switching Cre reporter "
             "(GenTarget SC018-Bsd)",
   organism="human", route="in vitro, 48-well plate",
   ratio="[Cre]/[CBB(Z)] = 1:8",
   dose_as_printed="complexed at 10 uM Cre + 20-160 uM CBB in 12.5 uL OptiMEM "
                   "(80 uM at the emitted 8:1), then diluted to 250 uL and "
                   "applied to the cells",
   timepoint="3 h exposure, medium exchanged, read 24 h later",
   endpoint_direction="higher is better delivery"),
 "Fig2d": dict(
   sheet="Figure 2", header_row=6, first_data_row=7, n_rep=2,
   cells="C6:P8",
   readout_as_printed="RFP conversion (%)",
   method="fluorescence microscopy only (Keyence BZ-X810), NOT flow cytometry",
   cargo="Cre recombinase (protein)",
   cargo_detail="purified Cre recombinase enzyme; no nucleic acid cargo",
   cell_line="CS cells = HEK293-loxP-GFP-RFP colour-switching Cre reporter "
             "(GenTarget SC018-Bsd)",
   organism="human", route="in vitro, 48-well plate",
   ratio="[Cre]/[CBB(Z)] = 1:8",
   dose_as_printed="as Fig 2c",
   timepoint="3 h exposure, medium exchanged, read 24 h later",
   endpoint_direction="higher is better delivery"),
 "Fig4c": dict(
   sheet="Figure 4", header_row=1, first_data_row=2, n_rep=4,
   cells="C1:AA5",
   readout_as_printed="eGFP/mCherry (%)",
   method="flow cytometry (Novocyte Quanteon, Agilent)",
   cargo="ABE8e ribonucleoprotein (adenine base editor ABE8e protein + "
         "chemically modified synthetic sgRNA, 1:1.1)",
   cargo_detail="RNP, i.e. protein AND guide RNA; genome editing, not "
                "reporter translation",
   cell_line="rd12 reporter cells (NIH/3T3 background; mCherry-[Rpe65 rd12 "
             "204 bp, c.130C>T p.R44X]-eGFP colour switch)",
   organism="mouse", route="in vitro, 48-well plate",
   ratio="[ABE8e]/[CBB(Z)] = 1:8",
   dose_as_printed="10 uM final RNP in the complex; applied at 0.1 uM ABE in "
                   "OptiMEM",
   timepoint="3 h exposure, medium exchanged at 3 h and 24 h, analysed 48 h "
             "after applying ABE",
   endpoint_direction="higher is better editing"),
 "Fig4d": dict(
   sheet="Figure 4", header_row=8, first_data_row=9, n_rep=3,
   cells="C8:L11",
   readout_as_printed="eGFP/mCherry (%)",
   method="fluorescence microscopy only (Keyence BZ-X810), NOT flow cytometry",
   cargo="ABE8e ribonucleoprotein (ABE8e protein + synthetic sgRNA)",
   cargo_detail="RNP, i.e. protein AND guide RNA",
   cell_line="rd12 reporter cells (NIH/3T3 background)",
   organism="mouse", route="in vitro, 48-well plate",
   ratio="[ABE8e]/[CBB(Z)] = 1:8",
   dose_as_printed="as Fig 4c",
   timepoint="analysed 48 h after applying ABE",
   endpoint_direction="higher is better editing"),
}


def compound_id(header):
    """A Source Data column header -> our compound id, or None for a control."""
    h = str(header).strip()
    if h in CONTROLS:
        return None
    if h.lower() in ("ox11", "ox-11", "cbbox-11", "cbbox11"):
        return "CBBox11"
    if h.upper().startswith("Z"):
        return "CBBZ" + h[1:]
    # bare integers 1..22 are CBB1..CBB22
    try:
        n = int(float(h))
    except ValueError:
        raise SystemExit("unrecognised Source Data column header: " + repr(h))
    if not 1 <= n <= 22:
        raise SystemExit("compound index out of range: " + repr(h))
    return "CBB%d" % n


def read_panel(wb, p):
    ws = wb[p["sheet"]]
    grid = [list(r) for r in ws.iter_rows(values_only=True)]
    hdr = grid[p["header_row"] - 1]
    reps = grid[p["first_data_row"] - 1: p["first_data_row"] - 1 + p["n_rep"]]
    out, ctrl = [], {}
    for col in range(2, len(hdr)):
        h = hdr[col]
        if h is None or str(h).strip() == "":
            continue
        vals = []
        for r in reps:
            v = r[col] if col < len(r) else None
            if v is None or str(v).strip() == "":
                continue
            vals.append(float(v))
        if not vals:
            continue
        cid = compound_id(h)
        if cid is None:
            ctrl[str(h).strip()] = vals
        else:
            out.append((cid, vals, col))
    return out, ctrl


def main():
    enum = {r["compound"]: r
            for r in csv.DictReader(open(os.path.join(_HERE,
                                                      "03_enumerated.csv")))}
    # CBB18/CBB19 exist only as two sequence candidates; collapse to one record
    # whose SMILES is deliberately EMPTY (see DECISIONS.md D-1).
    for base in ("CBB18", "CBB19"):
        cands = [k for k in enum if k.startswith(base + "_")]
        r = dict(enum[cands[0]])
        r["compound"] = base
        r["canonical_SMILES"] = ""
        r["candidate_SMILES"] = " | ".join(
            enum[c]["compound"] + "=" + enum[c]["canonical_SMILES"]
            for c in sorted(cands))
        enum[base] = r

    wb = openpyxl.load_workbook(XLSX, read_only=True, data_only=True)
    rows, controls = [], {}
    for pid, p in PANELS.items():
        found, ctrl = read_panel(wb, p)
        controls[pid] = {k: dict(replicates=v, mean=round(statistics.fmean(v), 4))
                         for k, v in ctrl.items()}
        items = [c for c, _, _ in found]
        recs = {c: enum[c] for c in items}
        if set(items) - set(enum):
            raise SystemExit("no enumeration record for " +
                             str(set(items) - set(enum)))
        flags = emit.emit_rows(
            items, recs,
            # Every flag is read out of the compound's OWN enumeration record.
            # Not one literal: a hard-coded "high" is what F-011 was.
            confidence=lambda i, r: r["confidence"],
            note=lambda i, r: r["note"],
            structure_assumed=lambda i, r: r.get("structure_assumed") or "",
            positional_ambiguity=lambda i, r: r.get("positional_ambiguity") or "",
            # value_is_bound and value_normalization_unit are left unset ON
            # PURPOSE and the reason is in 07_qc.json: these are workbook
            # numbers, so there is no colormap floor or clipped axis to be
            # bounded by, and this folder does NOT z-score -- the corpus's
            # per-experiment standardisation happens at merge time.
            label=pid)
        for (cid, vals, col), fl in zip(found, flags):
            e = enum[cid]
            rows.append(dict(
                panel=pid, compound=cid,
                canonical_SMILES=e["canonical_SMILES"],
                candidate_SMILES=e.get("candidate_SMILES", ""),
                value_mean=round(statistics.fmean(vals), 4),
                value_sd=(round(statistics.stdev(vals), 4)
                          if len(vals) > 1 else ""),
                n_replicates=len(vals),
                replicates=";".join("%g" % v for v in vals),
                value_source_file="source/41467_2026_69077_MOESM6_ESM.xlsx",
                value_source_sheet=p["sheet"],
                value_source_cells=p["cells"],
                value_source_column=col + 1,
                value_source_tier="1_publisher_source_data",
                readout_as_printed=p["readout_as_printed"],
                method=p["method"], cargo=p["cargo"],
                cargo_detail=p["cargo_detail"], cell_line=p["cell_line"],
                organism=p["organism"], route=p["route"], ratio=p["ratio"],
                dose_as_printed=p["dose_as_printed"],
                timepoint=p["timepoint"],
                endpoint_direction=p["endpoint_direction"],
                scaffold=e["scaffold"], tail=e["tail"], **fl))

    cols = ["panel", "compound", "canonical_SMILES", "candidate_SMILES",
            "value_mean", "value_sd", "n_replicates", "replicates",
            "value_source_file", "value_source_sheet", "value_source_cells",
            "value_source_column", "value_source_tier", "readout_as_printed",
            "method", "cargo", "cargo_detail", "cell_line", "organism",
            "route", "ratio", "dose_as_printed", "timepoint",
            "endpoint_direction", "scaffold", "tail",
            "confidence", "positional_ambiguity", "structure_assumed",
            "value_is_bound", "value_normalization_unit",
            "Experiment_value_within_model", "curation_note"]
    with open(os.path.join(_HERE, "05_values.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)
    json.dump(controls, open(os.path.join(_HERE, "05_controls.json"), "w"),
              indent=1)

    n_struct = sum(1 for r in rows if r["canonical_SMILES"])
    print("rows: %d  (%d with a SMILES, %d deliberately without)"
          % (len(rows), n_struct, len(rows) - n_struct))
    for pid in PANELS:
        k = [r for r in rows if r["panel"] == pid]
        print("  %-7s %2d compounds, n=%s, controls %s"
              % (pid, len(k), k[0]["n_replicates"],
                 sorted(controls[pid])))
    print("rows_with_a_caveat =", emit.caveat_count(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
