#!/usr/bin/env python3
"""Fig S3A (C2C12, first-round screen) measured as RATIOS OF BAR HEIGHTS.

Not emitted. Its only job is to test Fig S4C, which reports the SAME
formulations, the SAME cell line, the SAME readout and the SAME timepoint.

No y-axis calibration is needed and none is done: on a linear axis anchored at
zero, a ratio of bar heights IS a ratio of values, and the SM102 bar is 1.000
by construction. That makes this an independent read of a different panel with
a different scale, not a rerun of 04_digitize_figS4.py's calibration.
"""
import sys, os, json
import numpy as np
from PIL import Image
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import axis

# The raster comes out of the SI PDF in `source/`, not from a path handed in on
# the command line. A tracked output that depends on a session scratchpad is
# regenerable from nothing once the scratchpad is cleared -- error class 38,
# which has already cost this project `build_manifest.py` (F-023) and
# `lit_review/build_data.py`. The extracted PNG is written beside this script
# so a reviewer can look at exactly the pixels that were measured.
SI = os.path.join(HERE, "source", "nn5c19036_si_002.pdf")
PAGE, XREF = 9, 167
import pymupdf
_doc = pymupdf.open(SI)
_ids = [i[0] for i in _doc[PAGE - 1].get_images(full=True)]
if _ids != [XREF]:
    raise SystemExit(f"SI p.9 holds images {_ids}, expected [167]")
_px = pymupdf.Pixmap(_doc, XREF)
if _px.n > 3:
    _px = pymupdf.Pixmap(pymupdf.csRGB, _px)
IMG = os.path.join(HERE, "fig_p9_x167.png")
_px.save(IMG)
im = np.array(Image.open(IMG).convert("RGB")).astype(int)
g = im.mean(axis=2)
mx, mn = im.max(axis=2), im.min(axis=2)
sat = mx - mn
FILL = ((sat < 25) & (g > 70) & (g < 190)) | (sat >= 60)
DARK = (g < 90) & (sat < 60)
SOLID = FILL | DARK

BASE = 321                      # panel A's axis rule, x 58-551
LABELS = ["SM102"] + [f"1-F{i}" for i in range(1, 17)]

xt = [v for v in sorted(axis.detect_ticks_1d(g, "x", BASE + 3, BASE + 11,
                                             dark=140, min_len=2, max_len=14))
      if 100 < v < 560]
if len(xt) != len(LABELS):
    raise SystemExit(f"S3A: {len(xt)} category ticks, need {len(LABELS)}.")
loo = max(axis.loo_residuals(xt))
if loo > axis.LOO_TOL_PX:
    raise SystemExit(f"S3A: category ladder LOO {loo:.2f} px > tolerance.")

sil = SOLID[BASE - 2, :]
out = {"baseline_row": BASE, "n_ticks": len(xt),
       "x_tick_loo_max_px": round(loo, 3),
       "x_pitch_px": round(float(np.mean(np.diff(xt))), 3), "heights_px": {}}
for name, c in zip(LABELS, xt):
    ci = int(round(c))
    a = ci
    while a > 0 and sil[a - 1]: a -= 1
    b = ci
    while b < im.shape[1] - 1 and sil[b + 1]: b += 1
    cov = SOLID[:, a:b + 1].mean(axis=1)
    r = BASE - 2
    while r > 30 and cov[r - 1] >= 0.50:
        r -= 1
    out["heights_px"][name] = BASE - 2 - r
h0 = out["heights_px"]["SM102"]
if h0 < 20:
    raise SystemExit("S3A: the SM102 bar is too short to serve as the unit.")
out["ratio_to_SM102"] = {k: round(v / h0, 4) for k, v in out["heights_px"].items()}
with open(os.path.join(HERE, "04b_figS3A_ranks.json"), "w") as fh:
    json.dump(out, fh, indent=1)
print(json.dumps(out, indent=1))
