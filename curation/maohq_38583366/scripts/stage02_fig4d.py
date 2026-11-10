#!/usr/bin/env python3
"""Fig 4D -- flux per ORGAN for the six selected formulations, intraduodenal.

Tier 5. Same three-level method as stage01 (page / fill / stroke), same lib.

WHY THIS PANEL IS WORTH READING EVEN THOUGH NOTHING IS EMITTED FROM IT. Above
the bars the authors DREW six pie charts and printed one number on each: the
liver's share of the flux. They never computed that share from the bars. So
liver / sum(five organs) is a quantity external to this pipeline, six times
over, and one misread bar breaks it. This is the panel's own check on the
method stage01 uses.
"""
import json, os, sys
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "lib"))
import axis, bars                                        # noqa: E402

FIG = os.environ.get("STAGE02_FIG", os.path.join(
    HERE, "..", "source", "figs_elsevier", "gr4_lrg.jpg"))
OUT = os.environ.get("STAGE02_OUT", os.path.join(HERE, "..", "02_fig4d.json"))

AXIS_X = 1889
TICK_STRIP = (1880, 1887)
TICK_Y = (1000, 1620)
MAJOR_MIN_LEN = 15                 # measured: 19 (major) vs 12 (minor)
BASELINE_Y = 1566                  # rows 1564-1568
PLOT_X = (1890, 3290)
PAGE_ROWS = (1160, 1200)
MAJOR_VALUES = [1e7, 1e6, 1e5, 1e4]
MINOR_MULTIPLIERS = (8.2, 6.4, 4.6, 2.8)
CHROMA_MIN = 25
MIN_LEVEL_SEPARATION = 25

GATE_MINOR_RESID_PX = 1.5
GATE_LIVER_SHARE_PP = 1.0    # printed pie percentage vs the bars below it
GATE_N_BARS = 30

ORGANS = ["Heart", "Liver", "Spleen", "Lung", "Kidney"]
FORMS = ["AI-b3", "CII-b2", "CII-b3", "DI-a1", "DI-a3", "FIII-b3"]
# printed on the six pie charts directly above, left to right
PRINTED_LIVER_SHARE = [96.7, 91.4, 95.5, 95.2, 96.0, 87.5]

sys.path.insert(0, HERE)
from stage01_fig4b import stroke_centre, detect_bar_windows  # noqa: E402


def main():
    rgb = np.asarray(Image.open(FIG).convert("RGB")).astype(int)
    gray = rgb.mean(axis=2)
    mn = rgb.min(axis=2)

    ticks = sorted(t for t in axis.detect_ticks_1d(
        gray, "y", TICK_STRIP[0], TICK_STRIP[1], dark=140, min_len=3,
        max_len=40) if TICK_Y[0] < t < TICK_Y[1])

    def tick_len(t):
        band = gray[int(round(t)) - 3:int(round(t)) + 4, :AXIS_X].min(axis=0)
        x = AXIS_X - 1
        while x > AXIS_X - 200 and band[x] < 160:
            x -= 1
        return AXIS_X - 1 - x

    lens = {t: tick_len(t) for t in ticks}
    majors = sorted(t for t in ticks if lens[t] >= MAJOR_MIN_LEN)
    minors = sorted(t for t in ticks if lens[t] < MAJOR_MIN_LEN)
    if len(majors) != 4:
        raise SystemExit("stage02: expected 4 majors, got %r" % (majors,))
    cal = axis.calibrate(majors, MAJOR_VALUES, kind="log10", axis="y")

    predicted = {k * 10 ** (d - 1): cal.to_pixel(k * 10 ** (d - 1))
                 for d in (7, 6, 5, 4) for k in MINOR_MULTIPLIERS}
    minor_res = []
    for px in minors:
        v, p = min(predicted.items(), key=lambda kv: abs(kv[1] - px))
        minor_res.append({"measured_px": round(px, 3), "value": v,
                          "residual_px": round(px - p, 3)})
    mr = np.array([abs(m["residual_px"]) for m in minor_res])

    # bar windows, from a row just above the baseline
    row = rgb[BASELINE_Y - 8, PLOT_X[0]:PLOT_X[1]]
    on = (row.max(axis=1) - row.min(axis=1)) > CHROMA_MIN
    runs, s = [], None
    for i, v in enumerate(on):
        if v and s is None:
            s = i
        elif not v and s is not None:
            runs.append((PLOT_X[0] + s, PLOT_X[0] + i)); s = None
    if s is not None:
        runs.append((PLOT_X[0] + s, PLOT_X[1]))
    wins = [r for r in runs if r[1] - r[0] >= 10]
    if len(wins) != 30:
        raise SystemExit("stage02: expected 30 bars, found %d: %r"
                         % (len(wins), wins))

    page = float(np.median(mn[PAGE_ROWS[0]:PAGE_ROWS[1],
                              PLOT_X[0] + 20:PLOT_X[1] - 20]))
    out = []
    for i, (x0, x1) in enumerate(wins):
        interior = mn[BASELINE_Y - 25:BASELINE_Y - 5, x0 + 4:x1 - 4]
        fill = float(np.median(interior))
        strokev = float(np.percentile(mn[1150:BASELINE_Y - 2, x0:x1], 1))
        if page - fill < MIN_LEVEL_SEPARATION or \
           fill - strokev < MIN_LEVEL_SEPARATION:
            raise SystemExit("stage02: bar %d levels page=%.1f fill=%.1f "
                             "stroke=%.1f not separated" % (i, page, fill,
                                                            strokev))
        thr = (page + fill) / 2.0
        r = bars.bar_top(mn < thr, x0, x1, 1150, BASELINE_Y - 2,
                         baseline_row=BASELINE_Y)
        prof = (255.0 - mn[:, x0:x1]).mean(axis=1)
        centre, up, dn, peak = stroke_centre(prof, r["row"], 255.0 - page,
                                             255.0 - fill)
        out.append({"formulation": FORMS[i // 5], "organ": ORGANS[i % 5],
                    "x0": x0, "x1": x1, "bar_top_row": r["row"],
                    "flank_delta_px": r["flank_delta_px"],
                    "stroke_centre_row": round(centre, 3),
                    "value_photons_per_s": float(cal.to_value(centre))})

    checks = []
    for k, f in enumerate(FORMS):
        rows = out[5 * k:5 * k + 5]
        tot = sum(b["value_photons_per_s"] for b in rows)
        liver = [b for b in rows if b["organ"] == "Liver"][0]
        share = 100.0 * liver["value_photons_per_s"] / tot
        checks.append({"formulation": f, "printed_pct": PRINTED_LIVER_SHARE[k],
                       "computed_pct": round(share, 3),
                       "delta_pp": round(share - PRINTED_LIVER_SHARE[k], 3)})
    d = np.array([abs(c["delta_pp"]) for c in checks])

    failures = []
    if mr.max() > GATE_MINOR_RESID_PX:
        failures.append("held-out minor tick residual %.3f px > %.1f"
                        % (mr.max(), GATE_MINOR_RESID_PX))
    if d.max() > GATE_LIVER_SHARE_PP:
        failures.append("printed liver share vs bars: max |delta| %.3f pp > %.1f"
                        % (d.max(), GATE_LIVER_SHARE_PP))
    print("GATE minor_tick_residual_px_max %.4f" % mr.max())
    print("GATE liver_share_max_abs_delta_pp %.4f" % d.max())
    print("GATE liver_share_n_exact_to_1pp %d" % int((d <= 1.0).sum()))
    print("GATE n_bars %d" % len(out))

    ppd = cal.px_per_unit
    qc = {"panel": "Fig 4D",
          "calibration": cal.qc(),
          "major_ticks_px": [round(t, 3) for t in majors],
          "decade_gaps_px": [round(majors[i + 1] - majors[i], 3)
                             for i in range(3)],
          "minor_ticks_held_out_n": len(minor_res),
          "minor_tick_residual_px_max": round(float(mr.max()), 3),
          "minor_tick_residual_px_rms": round(float(np.sqrt((mr ** 2).mean())), 3),
          "precision": {"px_per_decade": round(ppd, 4),
                        "factor_per_px": round(10 ** (1.0 / ppd), 6)},
          "baseline_value_read": round(float(cal.to_value(BASELINE_Y)), 1),
          "liver_share_check": {
              "_what": "liver / sum(5 organ bars), against the percentage "
                       "PRINTED on the pie chart the authors drew above each "
                       "group. They never computed it from these bars, so it "
                       "is external to this read; one misread bar breaks it.",
              "rows": checks,
              "max_abs_delta_pp": round(float(d.max()), 3),
              "rms_delta_pp": round(float(np.sqrt((d ** 2).mean())), 3),
              "n_exact_to_1pp": int((d <= 1.0).sum()), "n": len(checks)},
          "text_claim_over_87_pct": {
              "claim": "'over 87% of the observed bioluminescent flux "
                       "localized in the liver' (Sec 3.3)",
              "min_printed": min(PRINTED_LIVER_SHARE),
              "min_computed": round(min(c["computed_pct"] for c in checks), 3)},
          "gate_failures": failures,
          "bars": out}
    with open(OUT, "w") as fh:
        json.dump(qc, fh, indent=1)
    print(json.dumps({k: v for k, v in qc.items() if k != "bars"}, indent=1))
    for b in out:
        print("  %-8s %-7s %.4e  flank=%.1f"
              % (b["formulation"], b["organ"], b["value_photons_per_s"],
                 b["flank_delta_px"]))
    if failures:
        raise SystemExit("stage02 GATE FAILED: " + "; ".join(failures))


if __name__ == "__main__":
    main()
