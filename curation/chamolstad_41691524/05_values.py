#!/usr/bin/env python3
"""Join the Fig S4A/S4C reads to Table S1's formulations and emit 05_values.csv.

Every quality flag is written on the ROW, through curation/lib/emit, from each
item's own record -- never a literal (error class 29).
"""
import sys, os, csv, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import emit

# ---- Table S1, SI pp. 27-28, from the PDF TEXT LAYER (tier 4, not pixels).
# fields: IL, IL mol%, phospholipid, PL mol%, cholesterol mol%, PEG lipid,
#         PEG mol%, N/P
T = """SM102|SM102|50|DSPC|10|38.5|DMG-PEG2k|1.5|6
1-F1|H9T6|25|DSPC|16|58.5|DMG-PEG2k|0.5|6
1-F2|H9T6|25|DOPE|16|58|DMG-PEG2k|1|7
1-F3|H9T6|25|DOPC|16|57.5|DMG-PEG2k|1.5|8
1-F4|H9T6|25|DOPG|16|57|DMG-PEG2k|2|9
1-F5|H9T6|30|DSPC|16|53|DMG-PEG2k|1|9
1-F6|H9T6|30|DOPE|16|53.5|DMG-PEG2k|0.5|8
1-F7|H9T6|30|DOPC|16|52|DMG-PEG2k|2|7
1-F8|H9T6|30|DOPG|16|52.5|DMG-PEG2k|1.5|6
1-F9|H9T6|35|DSPC|16|47.5|DMG-PEG2k|1.5|7
1-F10|H9T6|35|DOPE|16|47|DMG-PEG2k|2|6
1-F11|H9T6|35|DOPC|16|48.5|DMG-PEG2k|0.5|9
1-F12|H9T6|35|DOPG|16|48|DMG-PEG2k|1|8
1-F13|H9T6|40|DSPC|16|42|DMG-PEG2k|2|8
1-F14|H9T6|40|DOPE|16|42.5|DMG-PEG2k|1.5|9
1-F15|H9T6|40|DOPC|16|43|DMG-PEG2k|1|6
1-F16|H9T6|40|DOPG|16|43.5|DMG-PEG2k|0.5|7
2-F1|H9T6|35|DOPE|16|48.5|DMG-PEG2k|0.5|6
2-F2|H9T6|35|DOPE|16|48.5|DMG-PEG2k|0.5|8
2-F3|H9T6|35|DOPE|16|48.5|DMG-PEG2k|0.5|9
2-F4|H9T6|35|DOPE|16|48|DMG-PEG2k|1|6
2-F5|H9T6|35|DOPE|16|48|DMG-PEG2k|1|8
2-F6|H9T6|35|DOPE|16|48|DMG-PEG2k|1|9
2-F7|H9T6|35|DOPE|16|47|DMG-PEG2k|2|6
2-F8|H9T6|35|DOPE|16|47|DMG-PEG2k|2|8
2-F9|H9T6|35|DOPE|16|47|DMG-PEG2k|2|9
2-F10|H9T6|35|DOPC|16|48.5|DMG-PEG2k|0.5|6
2-F11|H9T6|35|DOPC|16|48.5|DMG-PEG2k|0.5|8
2-F12|H9T6|35|DOPC|16|48.5|DMG-PEG2k|0.5|9
2-F13|H9T6|35|DOPC|16|48|DMG-PEG2k|1|6
2-F14|H9T6|35|DOPC|16|48|DMG-PEG2k|1|8
2-F15|H9T6|35|DOPC|16|48|DMG-PEG2k|1|9
2-F16|H9T6|35|DOPC|16|47|DMG-PEG2k|2|6
2-F17|H9T6|35|DOPC|16|47|DMG-PEG2k|2|8
2-F18|H9T6|35|DOPC|16|47|DMG-PEG2k|2|9
2-F19|H9T6|35|DOPG|16|48.5|DMG-PEG2k|0.5|6
2-F20|H9T6|35|DOPG|16|48.5|DMG-PEG2k|0.5|8
2-F21|H9T6|35|DOPG|16|48.5|DMG-PEG2k|0.5|9
2-F22|H9T6|35|DOPG|16|48|DMG-PEG2k|1|6
2-F23|H9T6|35|DOPG|16|48|DMG-PEG2k|1|8
2-F24|H9T6|35|DOPG|16|48|DMG-PEG2k|1|9
2-F25|H9T6|35|DOPG|16|47|DMG-PEG2k|2|6
2-F26|H9T6|35|DOPG|16|47|DMG-PEG2k|2|8
2-F27|H9T6|35|DOPG|16|47|DMG-PEG2k|2|9"""
FORM = {}
for line in T.strip().split("\n"):
    f = line.split("|")
    FORM[f[0]] = dict(zip(["formulation", "IL", "IL_molpct", "phospholipid",
                           "PL_molpct", "chol_molpct", "PEG_lipid",
                           "PEG_molpct", "NP_ratio"], f))
for k, v in FORM.items():
    tot = sum(float(v[c]) for c in ("IL_molpct", "PL_molpct", "chol_molpct",
                                    "PEG_molpct"))
    assert abs(tot - 100.0) < 1e-9, (k, tot)   # P-check: mol% must sum to 100

geo = json.load(open(os.path.join(HERE, "04_figS4_geometry.json")))
enu = json.load(open(os.path.join(HERE, "02_enumerate.json")))
H9T6 = enu["series"]["6"]["smiles"]
# No SM-102 structure is written here on purpose. SM-102 appears only as the
# NORMALISER (its own bar is 1.000 by construction) and is not emitted as a
# row, so a SMILES for it would be an unused, unchecked claim sitting in a
# curation file -- and `reagents.smi` exists precisely so a name is never
# resolved to a structure by hand.

PANEL = {
    "S4A": {"cell": "DC2.4", "cell_type": "mouse dendritic cell line",
            "figure": "Figure S4A"},
    "S4C": {"cell": "C2C12", "cell_type": "mouse myoblast line",
            "figure": "Figure S4C"},
}
COMMON = dict(
    cargo="firefly luciferase (FLuc) mRNA",
    cargo_note="the SI figure caption says 'LNP-formulated luciferase mRNA' "
               "and does not state the UTR architecture or the nucleoside "
               "modification for this panel",
    assay="luciferase luminescence, microplate luminometer",
    readout_scale="relative to the SM102 benchmark LNP on the same plate "
                  "(the authors drew the 1.0 reference as a dashed line)",
    timepoint_h="24", route="in vitro transfection", organism="",
    dose="NOT STATED in the SI", LNP_preparation="pipette mixing",
    replicates="the caption says the data are representative of two "
               "independent experiments; n per bar is not stated",
    source_file="nn5c19036_si_002.pdf", source_page="11",
    value_tier="5_pixels")

rows, items, records = [], [], {}
for pid, P in PANEL.items():
    g = geo["panels"][pid]
    for name, b in g["bars"].items():
        if name == "SM102":
            continue            # the normaliser; 1.000 by construction
        it = f"{pid}:{name}"
        items.append(it)
        records[it] = {"panel": pid, "bar": b, "form": FORM[name],
                       "sm102_read": g["bars"].get("SM102", {}).get("value")}
        rows.append(dict(COMMON, item_id=it, panel=pid,
                         figure_panel=P["figure"], cell_line=P["cell"],
                         cell_type=P["cell_type"], **FORM[name],
                         IL_name="H9T6", IL_SMILES=H9T6,
                         value_rel_SM102=b["value"],
                         value_locus=(f"{P['figure']}, bar labelled {name} "
                                      f"(x {b['x0']}-{b['x1']} px, stroke top "
                                      f"row {b['row_stroke_top']}) in the "
                                      "embedded raster of SI p.11"),
                         bar_height_px=b["bar_height_px"],
                         crosscheck_delta_px=b["bar_top_crosscheck_delta_px"]))

# ---- quality flags, each derived from the item's own record
def _conf(it, rec):
    # Every row here is the same lipid, and that lipid is an unpurified
    # mixture: HPLC-ELSD puts the designed 6-tailed product's 4-tailed
    # congener at 55.26% of total area (SI p.5). The value is a tier-5 pixel
    # read as well. Neither is a "high".
    return "low" if rec["bar"]["stroke_width_px"] >= 5 else "medium"

def _note(it, rec):
    b, f = rec["bar"], rec["form"]
    base = ("Luciferase relative to SM102, read from the Fig S4 raster; the "
            "ionizable lipid H9T6 is a mixture (SI p.5 HPLC-ELSD: the "
            "4-tailed congener is 55.26% of total area) so the stored "
            "structure is the designed hexa-substituted product, not a pure "
            "compound.")
    if b["stroke_width_px"] >= 5:
        base += (f" The bar's top stroke measured {b['stroke_width_px']}px "
                 "against a 2-3px median, so the top rule may be merged with "
                 "an error-bar cap.")
    if b.get("top_within_6px_of_dashed_1.0"):
        base += (" The bar top is within 6px of the drawn 1.0 reference line, "
                 "which the detector treats as part of the silhouette.")
    return base

def _assumed(it, rec):
    return 1     # the 6-tail structure is the designed product, not isolated

flags = emit.emit_rows(items, records, confidence=_conf, note=_note,
                       structure_assumed=_assumed,
                       value_normalization_unit=lambda i, r: "",
                       label="chamolstad_41691524 Fig S4A/S4C")
for r, fl in zip(rows, flags):
    r.update(fl)

cols = list(rows[0].keys())
with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cols); w.writeheader(); w.writerows(rows)

print(f"wrote {len(rows)} rows, {len(cols)} columns")
print("rows_with_a_caveat =", emit.caveat_count(flags))
import collections
print("per panel:", collections.Counter(r["panel"] for r in rows))
print("confidence:", collections.Counter(r["confidence"] for r in rows))
