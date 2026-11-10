#!/usr/bin/env python3
"""Stage 05 -- the values this paper actually puts a number on.

WHERE THE VALUES ARE. Not in the 384-lipid screen. Fig. 1c/h (circos) and
Fig. 1d/i (Sankey) are the only places per-member potencies appear, and neither
carries a numeric scale of any kind -- the heat-map legend is a ramp annotated
"Min" and "Max", the bar ring has no axis, and the Sankey's y axis is a bare
arrow. 00_value_sources.json records that refusal with its evidence.

What DOES have a scale is the follow-up: the named lipids the paper carried
into mice. Three panels, and they are THREE DIFFERENT EXPERIMENTS -- different
route, different dose, different tissue. They are emitted as three
Experiment_IDs and are never pooled.

  Fig. 1n   intravenous,   0.5  mg/kg, whole-body, 12 LNPs, n = 5
  Fig. 1o   intrathecal,   0.05 mg/kg, brain,       3 LNPs, n = 5
  Fig. 1p   intratracheal, 0.25 mg/kg, lung,        3 LNPs, n = 5

HOW. These panels are vector, not raster. Tick marks are stroked segments, the
group MEAN is a rule the figure itself draws, and every category label and tick
value is a PDF text span. Nothing here is estimated from a pixel. The mean is
taken from the drawn rule rather than averaged from the dots: the dots are
jittered and some coincide exactly, so an average of the markers we can SEE is
not the mean the authors plotted -- 55 of the 60 Fig. 1n markers are separable.

The scale comes from curation/lib/axis.calibrate, which fits scale AND origin
together from detected tick positions and has no API to supply an origin
(error class 33). Fig. 1p's 10^6 tick is ABSENT from the vector layer, so that
panel is fitted on its three contiguous ticks and the fourth is HELD OUT and
predicted -- which is a stronger test than including it would have been.
"""
import csv
import json
import math
import os
import sys

import pymupdf

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import axis                                                      # noqa: E402
import emit                                                      # noqa: E402
from errors import DetectionError                                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ART = os.path.join(HERE, "s41587-026-03298-8.pdf")
FIG1_PAGE = 2                                   # 0-based

# Panel bounds, in PDF points on that page. Each is the plot box of ONE
# quantification panel; the axis rule and the ticks are then found inside it.
PANELS = {
 "Fig_1n": {
   "experiment_id": "SD_2026",
   "axis_x": (149.0, 153.5), "y": (480.0, 530.0), "plot_x": (152.0, 292.0),
   "label_y": (526.0, 535.0),
   "route": "intravenous injection", "dose_mg_per_kg": 0.5,
   "organ": "whole body",
   "organ_note": "the panel quantifies the whole-body IVIS image; the Results "
                 "text states the expression is hepatic ('all formulations "
                 "produced hepatic expression'), and the per-organ breakdown "
                 "for spleen, heart, lung and kidney is Supplementary Fig. 53, "
                 "which is NOT extracted here.",
   "cargo": "ABE-NanoLuc mRNA", "cargo_kb": 5.7, "timepoint_h": 6,
   "n_mice": 5},
 "Fig_1o": {
   "experiment_id": "SD_2026_2",
   "axis_x": (503.0, 506.3), "y": (405.0, 456.0), "plot_x": (506.0, 548.5),
   "label_y": (452.0, 462.0),
   "route": "intrathecal injection", "dose_mg_per_kg": 0.05,
   "organ": "brain",
   "organ_note": "quantified from the head region of the whole-body image; "
                 "the Results text reports brain expression for this route.",
   "cargo": "ABE-NanoLuc mRNA", "cargo_kb": 5.7, "timepoint_h": 6,
   "n_mice": 5},
 "Fig_1p": {
   "experiment_id": "SD_2026_3",
   "axis_x": (503.0, 506.1), "y": (475.0, 526.0), "plot_x": (505.7, 548.5),
   "label_y": (522.0, 532.0),
   "route": "intratracheal instillation", "dose_mg_per_kg": 0.25,
   "organ": "lung",
   "organ_note": "the figure legend and the panel heading both say "
                 "intratracheal instillation; this panel's own y-axis label "
                 "text reads 'Radiance Intrathecal injection', which is a "
                 "copy of panel o's label. See 08_contradictions.csv C-2.",
   "cargo": "ABE-NanoLuc mRNA", "cargo_kb": 5.7, "timepoint_h": 6,
   "n_mice": 5},
}

MEAN_RULE_PT = (7.0, 9.0)      # the wide rule: the group mean
CAP_RULE_PT = (3.4, 4.6)       # the narrow rules: mean +/- s.d.


def page():
    return pymupdf.open(ART)[FIG1_PAGE]


def horizontals(p, xlo, xhi, ylo, yhi):
    out = []
    for q in p.get_drawings():
        if not q.get("width"):
            continue
        for it in q["items"]:
            if it[0] != "l":
                continue
            a, b = it[1], it[2]
            if abs(a.y - b.y) > 0.05:
                continue
            x0, x1 = min(a.x, b.x), max(a.x, b.x)
            if x0 >= xlo and x1 <= xhi and ylo <= a.y <= yhi:
                out.append((x0, x1, a.y, x1 - x0))
    return sorted(set(out))


def tick_values(p, panel):
    """The tick LABELS, read as PDF text: '10' + a superscript exponent."""
    xlo, xhi = panel["axis_x"][0] - 12.0, panel["axis_x"][0]
    ylo, yhi = panel["y"]
    spans = []
    for b in p.get_text("dict")["blocks"]:
        if b.get("type") != 0:
            continue
        for l in b.get("lines", []):
            if tuple(round(v, 2) for v in l.get("dir")) != (1.0, 0.0):
                continue
            for s in l["spans"]:
                x, y = s["bbox"][0], s["bbox"][1]
                if xlo <= x <= xhi and ylo - 6 <= y <= yhi + 6:
                    spans.append((round(y, 1), x, s["text"].strip()))
    # group by baseline: a '10' and its superscript share a line
    byline = {}
    for y, x, t in spans:
        key = min(byline, key=lambda k: abs(k - y)) if byline else None
        if key is None or abs(key - y) > 1.5:
            key = y
        byline.setdefault(key, []).append((x, t))
    out = []
    for y, parts in byline.items():
        parts.sort()
        txt = "".join(t for _, t in parts)
        if not txt.startswith("10"):
            continue
        exp = txt[2:]
        if not exp.lstrip("-").isdigit():
            raise DetectionError(
                f"stage05: tick label {txt!r} is not 10^n. Read it, do not "
                "guess it.")
        out.append((y, float(exp), txt))
    return sorted(out)


def calibrate_panel(p, name, panel):
    xlo, xhi = panel["axis_x"]
    ylo, yhi = panel["y"]
    ticks = [h for h in horizontals(p, xlo, xhi, ylo, yhi) if h[3] < 2.6]
    ys = sorted(t[2] for t in ticks)
    labels = tick_values(p, panel)
    if len(labels) < 3:
        raise DetectionError(f"stage05 {name}: {len(labels)} tick labels.")

    # Pair each DETECTED tick with the label nearest it, then keep only the
    # longest run whose VALUES step uniformly -- axis.calibrate requires that,
    # and a panel with a tick missing from the vector layer must not be forced.
    paired = []
    for y in ys:
        lab = min(labels, key=lambda L: abs(L[0] - y))
        paired.append((y, lab[1], lab[2]))
    paired.sort()
    best = []
    i = 0
    while i < len(paired):
        j = i + 1
        while (j < len(paired)
               and abs((paired[j][1] - paired[j - 1][1])
                       - (paired[i + 1][1] - paired[i][1])) < 1e-9):
            j += 1
        if j - i > len(best):
            best = paired[i:j]
        i += 1
    if len(best) < 3:
        raise DetectionError(
            f"stage05 {name}: only {len(best)} ticks step uniformly in value; "
            f"paired = {paired}.")
    heldout = [q for q in paired if q not in best]

    cal = axis.calibrate([q[0] for q in best],
                         [10.0 ** q[1] for q in best],
                         kind="log10", axis="y")
    qc = dict(cal.qc())
    qc["ticks_detected_px"] = [round(y, 3) for y in ys]
    qc["tick_labels_read_as_text"] = [[round(q[0], 1), q[2]] for q in paired]
    qc["ticks_used_in_fit"] = [[round(q[0], 3), q[2]] for q in best]
    qc["ticks_held_out"] = []
    for y, v, txt in heldout:
        pred = math.log10(cal.to_value(y))
        qc["ticks_held_out"].append({
            "label": txt, "value_log10": v, "pixel": round(y, 3),
            "predicted_log10": round(pred, 4),
            "error_decades": round(pred - v, 4),
            "error_px": round((pred - v) * cal.px_per_unit, 3)})
        if abs(pred - v) > 0.02:
            raise DetectionError(
                f"stage05 {name}: held-out tick {txt} predicts "
                f"10^{pred:.3f}, {pred - v:+.3f} decades out.")
    return cal, qc


def groups(p, name, panel):
    xlo, xhi = panel["plot_x"]
    ylo, yhi = panel["y"]
    H = horizontals(p, xlo, xhi, ylo, yhi)
    means = [h for h in H if MEAN_RULE_PT[0] <= h[3] <= MEAN_RULE_PT[1]]
    caps = [h for h in H if CAP_RULE_PT[0] <= h[3] <= CAP_RULE_PT[1]]
    out = []
    for x0, x1, y, w in sorted(means):
        xc = (x0 + x1) / 2.0
        mine = sorted(c[2] for c in caps if abs((c[0] + c[1]) / 2 - xc) < 0.6)
        if len(mine) != 2:
            raise DetectionError(
                f"stage05 {name}: group at x={xc:.2f} has {len(mine)} s.d. "
                "caps, expected 2 (mean + s.d. and mean - s.d.).")
        out.append({"x": xc, "mean_px": y,
                    "cap_hi_px": mine[0], "cap_lo_px": mine[1]})
    return out


def labels(p, panel, n):
    """Category names: the 45-degree rotated text spans under the axis."""
    ylo, yhi = panel["label_y"]
    xlo, xhi = panel["plot_x"]
    got = []
    for b in p.get_text("dict")["blocks"]:
        if b.get("type") != 0:
            continue
        for l in b.get("lines", []):
            d = tuple(round(v, 2) for v in l.get("dir"))
            if d == (1.0, 0.0) or d == (0.0, -1.0):
                continue
            for s in l["spans"]:
                x, y = s["bbox"][0], s["bbox"][1]
                if xlo - 12 <= x <= xhi and ylo <= y <= yhi:
                    got.append((x, s["text"].strip()))
    got.sort()
    if len(got) != n:
        raise DetectionError(
            f"stage05: found {len(got)} rotated category labels for {n} "
            f"groups: {got}")
    return [t for _, t in got]


def main():
    p = page()
    rows, qc = [], {"_what": "Stage 05: values from the three Fig. 1 panels "
                             "that carry a calibratable axis.",
                    "panels": {}}
    for name, panel in PANELS.items():
        cal, calqc = calibrate_panel(p, name, panel)
        gs = groups(p, name, panel)
        names = labels(p, panel, len(gs))
        vals = []
        for g, nm in zip(gs, names):
            v = cal.to_value(g["mean_px"])
            hi = cal.to_value(g["cap_hi_px"])
            lo = cal.to_value(g["cap_lo_px"])
            vals.append({"lipid": nm, "value": v, "sd_hi": hi, "sd_lo": lo,
                         "x": g["x"], "mean_px": g["mean_px"]})
        logs = [math.log10(v["value"]) for v in vals]
        mu = sum(logs) / len(logs)
        sd = math.sqrt(sum((x - mu) ** 2 for x in logs) / len(logs))
        for v, lg in zip(vals, logs):
            rows.append({
                "experiment_id": panel["experiment_id"],
                "source_panel": name,
                "source_locus": f"{os.path.basename(ART)} page "
                                f"{FIG1_PAGE + 1}, panel "
                                f"{name.split('_')[1]}, group centre "
                                f"x={v['x']:.2f} pt, mean rule "
                                f"y={v['mean_px']:.3f} pt",
                "lipid_name": v["lipid"],
                "is_library_member": "yes" if v["lipid"].startswith(
                    ("SC-", "LC-")) else "no (benchmark)",
                "cargo": panel["cargo"], "cargo_kb": panel["cargo_kb"],
                "route": panel["route"],
                "dose_mg_per_kg": panel["dose_mg_per_kg"],
                "organ": panel["organ"], "timepoint_h": panel["timepoint_h"],
                "readout": "whole-animal bioluminescence radiance",
                "units": "photons s^-1 cm^-2 sr^-1",
                "n_biological_replicates": panel["n_mice"],
                "value": f"{v['value']:.4g}",
                "value_log10": f"{lg:.4f}",
                "sd_upper": f"{v['sd_hi']:.4g}",
                "sd_lower": f"{v['sd_lo']:.4g}",
                "z_within_experiment": f"{(lg - mu) / sd:.4f}" if sd else "",
                "n_in_experiment": len(vals)})
        qc["panels"][name] = {
            "calibration": calqc,
            "n_groups": len(gs),
            "category_labels_read_as_text": names,
            "group_x_spacing_pt": [round(gs[i + 1]["x"] - gs[i]["x"], 4)
                                   for i in range(len(gs) - 1)],
            "log10_mean": round(mu, 4), "log10_sd_pop": round(sd, 4),
            "values": {v["lipid"]: round(math.log10(v["value"]), 4)
                       for v in vals}}

    # --- flags, through the ONE writer -------------------------------------
    recs = {f"{r['experiment_id']}|{r['lipid_name']}": r for r in rows}
    keys = list(recs)

    def conf(k, rec):
        # `confidence` is about the STRUCTURE on the row, not the value, and
        # that is the uncomfortable part here: for seven of these twelve names
        # this curation has NO structure at all. Derived from the record:
        #   benchmark  -> high   a catalogue reagent, unambiguous
        #   LC-1       -> medium pinned to H4+L4+Ta1+Tb1 by Supplementary
        #                        Table 5's CAS numbers, but its blocks are read
        #                        from a schematic figure
        #   the rest   -> low    the paper DRAWS them and the SI characterises
        #                        them, but which of the 384 enumerated members
        #                        each one is has not been determined (D-2), so
        #                        the row would carry a name and no structure.
        if rec["is_library_member"] == "no (benchmark)":
            return "high"
        return "medium" if rec["lipid_name"] == "LC-1" else "low"

    def note(k, rec):
        bits = [f"value digitised from the PDF's own vector geometry in "
                f"{rec['source_panel']} (mean rule as drawn, axis fitted from "
                f"detected tick marks); {rec['route']} at "
                f"{rec['dose_mg_per_kg']} mg/kg, {rec['organ']}, "
                f"{rec['cargo']}, {rec['timepoint_h']} h"]
        if rec["lipid_name"] == "LC-1":
            bits.append("structure pinned to the enumerated member "
                        "H4-L4-Ta1-Tb1 by Supplementary Table 5's CAS numbers")
        elif rec["is_library_member"] != "no (benchmark)":
            bits.append("NO STRUCTURE is assigned to this row: the paper draws "
                        "this lipid and the SI gives its formula, but which of "
                        "the 384 enumerated members it is has not been "
                        "determined and its printed formula matches between 2 "
                        "and 12 of them (DECISIONS.md D-2)")
        if int(rec["n_in_experiment"]) < 5:
            bits.append(f"this experiment has only "
                        f"{rec['n_in_experiment']} lipids, so its within-"
                        "experiment z-score is weakly determined")
        return "; ".join(bits) + "."

    flags = emit.emit_rows(
        keys, recs, confidence=conf, note=note,
        value_normalization_unit=lambda k, rec: "experiment",
        positional_ambiguity=lambda k, rec: "0",
        structure_assumed=lambda k, rec: "0",
        label="sd_2026 values")
    for r, f in zip(rows, flags):
        r.update(f)
    qc["rows_with_a_caveat"] = emit.caveat_count(flags)
    qc["n_rows"] = len(rows)

    cols = ([k for k in rows[0] if k not in emit.columns()] + emit.columns())
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})
    with open(os.path.join(HERE, "05_values_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(f"05_values.csv: {len(rows)} rows across "
          f"{len(set(r['experiment_id'] for r in rows))} experiments")
    for k, v in qc["panels"].items():
        print(" ", k, v["n_groups"], "groups;", v["category_labels_read_as_text"])
        print("     residual_px", v["calibration"].get("residual_px"),
              "held-out", v["calibration"]["ticks_held_out"])


if __name__ == "__main__":
    main()
