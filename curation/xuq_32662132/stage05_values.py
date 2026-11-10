"""Stage 05 - every delivery value in reach, with its exact source location.

Two origins, kept apart on purpose:
  * Figure S3  - PRINTED numerals on the raster. Transcribed, not estimated.
  * Figure S2  - bar geometry, measured by digitize_figS2.py through
                 curation/lib (axis.calibrate + bars.bar_top).

Figure S2b (CCK-8 viability) is carried here as a check, flagged
is_delivery_endpoint = 0, because the corpus endpoint set is delivery.
"""
import csv, json, os

HERE = os.path.dirname(os.path.abspath(__file__))

CELL = "primary human CD8+ T lymphocytes"

# Figure S3: read off source/figS3_page6_xref106.png at 6x. Glyph height ~22 px.
FIG_S3 = [
    ("93-O17S",   5,  2.68), ("93-O17S",  15, 4.16),
    ("93-O17S",  30,  3.19), ("93-O17S",  60, 7.05),
    ("9322-O17S", 5,  8.13), ("9322-O17S", 15, 6.83),
    ("9322-O17S", 30, 11.10), ("9322-O17S", 60, 0.74),
    ("Lipofectamine 2000", 5, 2.52), ("untreated", None, 0.44),
]

# Figure S2a category labels, in drawn left-to-right order. Confirmed BY the
# data, not assumed: 16:4:1:1 appears at positions 1, 6 and 10 and those three
# bars read to the same pixel row (07_qc.json repeated_formulation_identity).
S2A_LABELS = ["16:8:1:1", "16:4:1:1", "16:1:1:1", "16:0:1:1",
              "16:4:8:1", "16:4:4:1", "16:4:1:1", "16:4:0:1",
              "16:4:1:8", "16:4:1:4", "16:4:1:1", "16:4:1:0",
              "16:0:0:0", "UT"]
S2B_LABELS = ["6", "15", "30", "45", "60", "UT"]
S2C_LABELS = ["0.5", "1", "2", "5", "7", "18", "26", "UT"]
S2D_LABELS = ["6", "15", "30", "60", "120"]          # UT is drawn as no bar


def main():
    geo = json.load(open(os.path.join(HERE, "05_figS2_geometry.json")))
    rows = []

    for lip, dose, pct in FIG_S3:
        rows.append(dict(
            value_id="S3-%s-%s" % (lip, dose), source_figure="Figure S3",
            source_location="source/figS3_page6_xref106.png, histogram inset label",
            tier="5a_printed_numerals_on_a_raster",
            lipidoid=lip,
            readout="percent GFP-positive cells", units="percent", value=pct,
            cargo="EGFP mRNA", cargo_type="mRNA", model="in_vitro",
            model_target=CELL, route="cell culture",
            dose_lipidoid_ug_per_mL=dose, dose_nucleicacid="not stated for this panel",
            timepoint="not stated for this panel", assay="flow cytometry",
            replicates="not stated", is_delivery_endpoint=1,
            is_control=int(lip in ("Lipofectamine 2000", "untreated")),
            geometry_note="none - the number is printed"))

    def add_bars(pid, labels, readout, units, kind, extra):
        rec = geo[pid]
        if "REFUSED" in rec:
            raise SystemExit("panel %s was refused by the digitiser: %s" % (pid, rec["REFUSED"]))
        for b, lab in zip(rec["bars"], labels):
            if b is None:
                continue
            is_ctl = lab in ("UT", "16:0:0:0")
            d = dict(
                value_id="S2%s-%s" % (pid, lab), source_figure="Figure S2%s" % pid,
                source_location=("source/figS2_page5_xref101.png, bar window x=[%d,%d), "
                                 "top row %g" % (b["x0"], b["x1"], b["top_row"])),
                tier="5b_bar_geometry",
                lipidoid=("untreated" if lab == "UT" else "93-O17S"),
                readout=readout, units=units, value=b["value"],
                cargo="FLuc mRNA" if kind != "viability" else "FLuc mRNA",
                cargo_type="mRNA", model="in_vitro", model_target=CELL,
                route="cell culture", replicates="n=3 (S2a,b) / n=2 (S2c,d) per caption",
                is_delivery_endpoint=int(kind != "viability"),
                is_control=int(is_ctl),
                geometry_note=("flank delta %.1f px, coverage %.3f, width %d px"
                               % (b["detail"]["flank_delta_px"], b["detail"]["coverage"],
                                  b["detail"]["width_px"])))
            d.update(extra(lab))
            rows.append(d)

    add_bars("a", S2A_LABELS, "firefly luciferase luminescence", "RLU", "luminescence",
             lambda lab: dict(formulation_w_w=("not applicable" if lab == "UT" else lab),
                              dose_lipidoid_ug_per_mL="not stated for this panel",
                              dose_nucleicacid="not stated for this panel",
                              timepoint="not stated for this panel", assay="luciferase assay"))
    add_bars("b", S2B_LABELS, "cell viability (CCK-8), normalised to untreated", "percent",
             "viability",
             lambda lab: dict(formulation_w_w="16:4:1:1",
                              dose_lipidoid_ug_per_mL=(None if lab == "UT" else float(lab)),
                              dose_nucleicacid="not stated for this panel",
                              timepoint="6 h per the Fig S2 caption; the Experimental "
                                        "Section says 7 or 21 h (see 08_contradictions.csv C-05)",
                              assay="CCK-8"))
    add_bars("c", S2C_LABELS, "firefly luciferase luminescence", "RLU", "luminescence",
             lambda lab: dict(formulation_w_w="16:4:1:1",
                              dose_lipidoid_ug_per_mL="not stated for this panel",
                              dose_nucleicacid="not stated for this panel",
                              timepoint=("not applicable" if lab == "UT" else lab + " h"),
                              assay="luciferase assay"))
    add_bars("d", S2D_LABELS, "firefly luciferase luminescence", "RLU", "luminescence",
             lambda lab: dict(formulation_w_w="16:4:1:1",
                              dose_lipidoid_ug_per_mL=float(lab),
                              dose_nucleicacid="not stated for this panel",
                              timepoint="not stated for this panel", assay="luciferase assay"))

    # one printed outcome that is not in any panel we hold
    rows.append(dict(
        value_id="ABSTRACT-cre-recombination", source_figure="Abstract (Europe PMC core record)",
        source_location="abstractText, final sentence", tier="4_pdf_text_table",
        lipidoid="the lead lipidoid - NOT NAMED in the abstract",
        readout="percent gene recombination in mouse T lymphocytes", units="percent",
        value=8.2, cargo="Cre recombinase mRNA", cargo_type="mRNA", model="in_vivo",
        model_target="mouse T lymphocytes (Ai14 reporter, spleen)",
        route="intravenous", dose_lipidoid_ug_per_mL=None,
        dose_nucleicacid="15 ug Cre mRNA + 150 ug lipidoid per mouse, two IV injections "
                         "5 days apart (SI section 9)",
        timepoint="day 10 after the first injection", assay="flow cytometry of splenocytes",
        replicates="not stated", is_delivery_endpoint=1, is_control=0,
        formulation_w_w="16:4:1:1",
        geometry_note="printed statistic, recorded as evidence; not emitted as a row because "
                      "the abstract does not say WHICH lipidoid it belongs to"))

    cols = ["value_id", "source_figure", "source_location", "tier", "lipidoid", "readout",
            "units", "value", "cargo", "cargo_type", "model", "model_target", "route",
            "formulation_w_w", "dose_lipidoid_ug_per_mL", "dose_nucleicacid", "timepoint",
            "assay", "replicates", "is_delivery_endpoint", "is_control", "geometry_note"]
    for r in rows:
        for c in cols:
            r.setdefault(c, "")
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader(); w.writerows(rows)
    nd = sum(r["is_delivery_endpoint"] for r in rows)
    print("values:", len(rows), "| delivery-endpoint:", nd,
          "| controls:", sum(r["is_control"] for r in rows),
          "| distinct lipidoids:", len({r["lipidoid"] for r in rows}))

if __name__ == "__main__":
    main()
