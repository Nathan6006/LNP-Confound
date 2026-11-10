#!/usr/bin/env python3
"""ZG_2020 -- the checks the authors left lying around, computed.

Run:  python3 -I stage03_figure_checks.py

NOTHING HERE DIGITISES A VALUE. The only article figures we can obtain are the
Europe PMC web JPEGs at 709 px on the long edge, with no page rect, so native
dpi is not derivable and axis.native_tol() would refuse them. Every check below
is on a COUNT, a SUM, a SIGN or a RANK -- quantities a low-resolution raster
can still settle. CLAUDE.md section 5's own rule: a low-resolution axis cannot
pin a value but can still test one.

The four checks:

  C-a  The Fig 2b molar-ratio table: four formulations, each must sum to 100.
       Drawn under the table as four columns of numerals and never summed by
       the authors -- a free check, and NOT a tautology of our own derivation
       because we transcribe four independent numerals per column.
  C-b  The factorial: 10 phospholipids x 3 sterols x 4 ratios must equal the
       120 the paper states twice (Fig 2c "120 total LNPs", Methods "formulate
       120 LNPs"). If a claimed factorial does not multiply out, it is not one.
  C-c  Fig 3h sign by COLOUR IDENTITY: the caption prints "five constrained
       phospholipids that perform well in the top 10%", and the panel draws
       enriched bars teal and depleted bars grey. Counting teal bars tests the
       printed five against the drawing.
  C-d  Fig 3h rank: the text names A-11 "most highly enriched". The tallest
       teal bar must therefore be A-11.
"""
import json, os, sys

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))

# Outputs go to OUTDIR, which is HERE for a normal run. The seed harnesses
# redirect it to a temp directory: a MUTATED copy must never overwrite the
# real curation outputs. It did, and the determinism check caught it --
# reference_structures.csv was left holding seed G3's REGIOISOMER anchors.
OUTDIR = HERE
FIG3 = os.path.join(HERE, "source", "epmc_suppl", "extracted",
                    "BTM2-5-e10161-g013.jpg")

# --- Fig 2b, transcribed by eye from the figure's own numerals -------------
# columns are the four formulation ratios; rows as printed in the panel.
FIG2B = {
    "Ionizable Lipid": [15, 15, 19, 26],
    "Cholesterol":     [28, 38.5, 44, 46.5],
    "PEG-Lipid":       [2, 1.5, 2, 2.5],
    "Phospholipid":    [55, 45, 35, 25],
}
N_PHOSPHOLIPID, N_STEROL, N_RATIO = 10, 3, 4
PAPER_SAYS_TOTAL, PAPER_SAYS_PASSED_QC = 120, 109

LABELS = ["A-9", "A-10", "A-11", "A-12", "A-13", "A-14", "A-15", "A-16",
          "A-17", "A-17-2Z"]
# Panel (h) bar-slot centres, in full-image x. Found by colour-run detection
# (not assumed): the ten slots sit on a 22.8 px pitch, residual < 1.1 px.
PANEL_H_BAND = (520, 730)
PANEL_H_CENTRES = [276, 300, 321, 344, 367, 390, 412, 435, 458, 481]


def main():
    qc = {}

    # ---- C-a: every formulation sums to 100 -------------------------------
    sums = [round(sum(FIG2B[k][i] for k in FIG2B), 6) for i in range(4)]
    qc["C_a_molar_ratio_sums"] = {
        "per_formulation_sum": sums,
        "n_summing_to_100": sum(1 for s in sums if abs(s - 100) < 1e-9),
        "n_formulations": 4,
        "offenders": [i + 1 for i, s in enumerate(sums) if abs(s - 100) > 1e-9],
        "_not_a_tautology": "Each column is four numerals transcribed "
                            "independently from the panel; nothing in our "
                            "pipeline derives one of them from the others.",
        "table": FIG2B,
    }

    # ---- C-b: the factorial actually multiplies out -----------------------
    prod = N_PHOSPHOLIPID * N_STEROL * N_RATIO
    qc["C_b_factorial"] = {
        "phospholipids": N_PHOSPHOLIPID, "sterols": N_STEROL,
        "molar_ratios": N_RATIO, "product": prod,
        "paper_states_total": PAPER_SAYS_TOTAL,
        "factorial_closes": prod == PAPER_SAYS_TOTAL,
        "passed_qc": PAPER_SAYS_PASSED_QC,
        "dropped_at_qc": PAPER_SAYS_TOTAL - PAPER_SAYS_PASSED_QC,
        "_which_were_dropped": "NOT IDENTIFIED ANYWHERE IN THE PAPER. The 11 "
                               "failures are shown only as a change in the "
                               "Fig 2c pie, so the 109 surviving members "
                               "cannot be enumerated, let alone addressed.",
    }

    # ---- C-c / C-d: Fig 3h sign by colour, and the rank of A-11 -----------
    a = np.asarray(Image.open(FIG3).convert("RGB")).astype(int)
    y0, y1 = PANEL_H_BAND
    sub = a[y0:y1, :]
    r, g, b = sub[:, :, 0], sub[:, :, 1], sub[:, :, 2]
    teal = (g > r + 18) & (g > b + 8) & (g < 235) & (g > 120)
    grey = (abs(r - g) < 14) & (abs(g - b) < 14) & (r > 140) & (r < 215)

    bars = []
    for lab, c in zip(LABELS, PANEL_H_CENTRES):
        col = slice(c - 6, c + 7)
        tn, qn = int(teal[:, col].sum()), int(grey[:, col].sum())
        if tn > qn and tn > 40:
            sign, h = "positive", int(
                np.ptp(np.nonzero(teal[:, col].any(1))[0]) + 1)
        elif qn > tn and qn > 40:
            sign, h = "negative", None       # see note below
        else:
            sign, h = "about_zero", 0
        bars.append({"label": lab, "teal_px": tn, "grey_px": qn,
                     "sign": sign, "teal_height_px": h})

    pos = [x for x in bars if x["sign"] == "positive"]
    qc["C_c_fig3h_sign_by_colour"] = {
        "_method": "Colour identity, not a bar-top detector. Enriched bars "
                   "are drawn teal and depleted bars grey, so the SIGN is "
                   "carried by hue and survives a 709 px raster. Separation "
                   "is categorical: every negative bar has EXACTLY 0 teal "
                   "pixels and every positive bar has <=44 grey pixels "
                   "against 260-572 teal.",
        "bars": bars,
        "n_positive": len(pos),
        "n_negative": sum(1 for x in bars if x["sign"] == "negative"),
        "n_about_zero": sum(1 for x in bars if x["sign"] == "about_zero"),
        "n_bars": len(bars),
        "caption_prints": 5,
        "matches_caption": len(pos) == 5,
        "positive_labels": [x["label"] for x in pos],
        "_NOT_MEASURED": "The MAGNITUDE of the negative bars. Our grey mask "
                         "over-captures below the baseline (all four returned "
                         "an identical 93 px extent, which the drawing "
                         "contradicts), so negative heights are reported as "
                         "None rather than as a number we do not trust. The "
                         "sign is unaffected: it comes from hue, not extent.",
    }
    tallest = max(pos, key=lambda x: x["teal_height_px"]) if pos else None
    qc["C_d_fig3h_rank"] = {
        "text_claims_most_enriched": "A-11",
        "tallest_teal_bar": tallest["label"] if tallest else None,
        "teal_heights_px": {x["label"]: x["teal_height_px"] for x in pos},
        "matches_text": bool(tallest and tallest["label"] == "A-11"),
        "_precision": "A RANK, not a value. 44 px vs the next 40 px is a 4 px "
                      "margin on a 709 px-wide figure; it orders A-11 first "
                      "but we do not convert any of these to enrichment "
                      "units and no emitted number depends on them.",
    }

    for k in ("C_a_molar_ratio_sums", "C_b_factorial",
              "C_c_fig3h_sign_by_colour", "C_d_fig3h_rank"):
        print(f"METRIC {k}=" + json.dumps(
            {kk: vv for kk, vv in qc[k].items() if not kk.startswith("_")
             and kk not in ("bars", "table", "teal_heights_px")}))

    fails = []
    if qc["C_a_molar_ratio_sums"]["offenders"]:
        fails.append("C-a molar ratios do not sum to 100")
    if not qc["C_b_factorial"]["factorial_closes"]:
        fails.append("C-b the claimed factorial does not multiply out")
    if not qc["C_c_fig3h_sign_by_colour"]["matches_caption"]:
        fails.append("C-c teal bar count != the caption's five")
    if not qc["C_d_fig3h_rank"]["matches_text"]:
        fails.append("C-d tallest enriched bar is not A-11")
    print("METRIC figure_check_failures=" + json.dumps(fails))

    with open(os.path.join(OUTDIR, "04_figure_checks.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    if fails:
        raise SystemExit("GATE FAILED: " + "; ".join(fails))
    return 0


if __name__ == "__main__":
    sys.exit(main())
