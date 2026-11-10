#!/usr/bin/env python3
"""05_values.csv -- every value, with the exact panel or table cell it came from.

    python3 -I stage05_values.py

TWO CARGOES, AND THEY ARE NOT POOLED. The paper reports siRNA knockdown and
mRNA expression. They are different cargoes measured in different assays, and
Figure 1a draws them as two separate stripes with two separate legends and two
separate binning rules (0-25/25-70/70-100 % knockdown; log luminescence
0-1/1-2/>2). They are kept as separate records here and the siRNA half is not
emitted at all -- it is a re-binning of reference 11's data, which the corpus
already holds as KZ_2016. See 00_value_sources.json and DECISIONS.md D-2.

THREE RECORD TYPES, three different `source` strings:
  fig1a_mRNA    465 ordinal readings, IGROV-1, in vitro      <- emitted
  fig1a_siRNA   465 ordinal readings, prior data, re-binned  <- NOT emitted
  fig2b/fig3b   42 liver radiances, C57BL/6, in vivo         <- emitted
and the tier-4 SI table rides along on the in-vivo rows as context, never as a
value.
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
HEAT = os.path.join(HERE, "04_fig1_heatmap.json")
BARS = os.path.join(HERE, "04_invivo_bars.json")
SI = os.path.join(HERE, "02_si_table_S3.csv")
OUT = os.path.join(HERE, "05_values.csv")

FIELDS = ["lipid", "cargo", "record", "source_locator", "value_tier",
          "raw_value", "unit", "band_lo", "band_hi", "value_is_bound",
          "emitted", "si_diameter_nm", "si_PDI", "si_pKa", "si_EE_pct",
          "note"]


def main():
    heat = json.load(open(HEAT))
    bars = json.load(open(BARS))
    si = {r["lipid"]: r for r in csv.DictReader(open(SI))}

    rows = []
    for c in heat["cells"]:
        loc = (f"Fig 1a, row {c['row']} ({c['core']}), column {c['col']} "
               f"({c['tail']})")
        rows.append({
            "lipid": c["lipid"], "cargo": "mRNA", "record": "fig1a_mRNA",
            "source_locator": loc + ", LOWER stripe, cell centre "
                              f"({c['x_px']}, {c['mRNA_y_px']}) px",
            "value_tier": "5_pixels",
            "raw_value": c["mRNA_ordinal"],
            "unit": "ordinal 1-3 of log10 luminescence bands "
                    "(1 = log 0-1, 2 = log 1-2, 3 = > log 2)",
            "band_lo": c["mRNA_lo"], "band_hi": c["mRNA_hi"],
            "value_is_bound": c["mRNA_bound"], "emitted": "yes",
            "si_diameter_nm": "", "si_PDI": "", "si_pKa": "", "si_EE_pct": "",
            "note": f"nearest-swatch distance {c['mRNA_dist']} RGB",
        })
        rows.append({
            "lipid": c["lipid"], "cargo": "siRNA", "record": "fig1a_siRNA",
            "source_locator": loc + ", UPPER stripe, cell centre "
                              f"({c['x_px']}, {c['siRNA_y_px']}) px",
            "value_tier": "5_pixels",
            "raw_value": c["siRNA_ordinal"],
            "unit": "ordinal 1-3 of % knockdown bands "
                    "(1 = 0-25%, 2 = 25-70%, 3 = 70-100%)",
            "band_lo": c["siRNA_lo"], "band_hi": c["siRNA_hi"],
            "value_is_bound": c["siRNA_bound"], "emitted": "no",
            "si_diameter_nm": "", "si_PDI": "", "si_pKa": "", "si_EE_pct": "",
            "note": "NOT EMITTED: a re-binning of reference 11 (Zhou et al., "
                    "PNAS 2016, PMID 26729861), already in the corpus as "
                    "KZ_2016 at four levels. See DECISIONS.md D-2.",
        })

    for pname, panel in bars.items():
        fig = "Fig 2b" if pname == "fig2b" else "Fig 3b"
        for b in panel["bars"]:
            s = si.get(b["lipid"], {})
            rows.append({
                "lipid": b["lipid"], "cargo": "mRNA",
                "record": pname,
                "source_locator": f"{fig}, bar at x tick {b['x_tick_px']} px"
                                  + ("" if b["value"] is None else
                                     f", top row {b['bar_top_px']} px"),
                "value_tier": "5_pixels",
                "raw_value": "" if b["value"] is None else b["value"],
                "unit": "average liver radiance, p/sec/cm2/sr",
                "band_lo": "", "band_hi": "",
                "value_is_bound": "",
                "emitted": "no" if b["value"] is None else "yes",
                "si_diameter_nm": s.get("diameter_nm", ""),
                "si_PDI": s.get("PDI", ""), "si_pKa": s.get("pKa", ""),
                "si_EE_pct": s.get("encapsulation_efficiency_pct", ""),
                "note": b.get("note", ""),
            })

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    n_em = sum(1 for r in rows if r["emitted"] == "yes")
    print(f"wrote {len(rows)} value records to {OUT}; {n_em} emittable")
    for rec in ("fig1a_mRNA", "fig1a_siRNA", "fig2b", "fig3b"):
        sel = [r for r in rows if r["record"] == rec]
        print(f"  {rec:<12} {len(sel):>4} records, "
              f"{sum(1 for r in sel if r['emitted'] == 'yes'):>4} emittable")


if __name__ == "__main__":
    main()
