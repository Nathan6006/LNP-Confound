"""Run ONE seeded mutation, in its own process, and exit with a real code.

Why this file exists as well as stage04_selfcheck.py: that module decides a
seed's fate from a caught exception inside one interpreter. This one is the
harness the rule actually asks for -- one `python3 -I` process per seed, the
verdict taken from `returncode`, never from stdout and never through a pipe
(`python3 ... | tail` returns tail's exit code; that is F-012's shape).

Two refusals are built in, both because a seed that looks clean is worse than
no seed:

  * the LAST LINE of stdout must be `METRIC <n_bars> <v0> <v1> ...`. A seed
    that dies before producing its metric -- a TypeError, a bad index, a
    refusal upstream of the thing being tested -- has an absent metric line
    and is NOT CREDITED, however it exited.
  * exit 2 means "refused", exit 0 means "read", and the caller credits a
    seed only when the AGREEMENT COUNT MOVED against the control. A non-zero
    exit on its own earns nothing.

Usage: python3 -I seed_runner.py <seed-name> [arg]
       exit 0 = panel read, 2 = detector refused, 3 = metric never reached
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)
import axis                                                   # noqa: E402
import stage01_digitize as S                                  # noqa: E402
from errors import DetectionError, CalibrationError           # noqa: E402

PANEL = {"page": 23, "xref": 244, "spine_x": 597, "y0": 99, "y1": 338,
         "x_right": 1065, "paired": 0.22}
LABELS = S.S13_ORGANS[1][1]           # liver: 6 cell types


def read(img, derive_bounds=True):
    """Read the test panel.

    `derive_bounds` re-runs find_spines on the image being read, so the panel
    bounds come from the IMAGE and not from a constant. Without it the
    spine_erased seed is a NO-OP: the 100%-level check compares the fitted
    scale against PANEL["y0"], and a constant cannot notice that the spine it
    names has been rubbed out. That is the "removed the wrong object" shape.
    """
    y0, y1 = PANEL["y0"], PANEL["y1"]
    if derive_bounds:
        cands = [sp for sp in S.find_spines(img)
                 if abs(sp[0] - PANEL["spine_x"]) <= 3
                 and abs(sp[1] - PANEL["y0"]) <= 60]
        if not cands:
            raise DetectionError(
                "read: no y-axis spine found near x=%d, y=%d. The panel "
                "bounds are derived from the image, so an erased spine "
                "stops the read here." % (PANEL["spine_x"], PANEL["y0"]))
        _, y0, y1 = cands[0]
    vals, diag, cal = S.read_panel(
        img, PANEL["spine_x"], y0, y1, PANEL["x_right"],
        LABELS, S.CAL_TOL, S.LOO_TOL, PANEL["paired"])
    return vals, diag, cal


def mutate(img, name, arg):
    out = img.copy()
    base, _, cal = read(img)
    ticks = [int(round(t)) for t in cal["ticks_px"]]
    baseline = int(round(cal["baseline_row"]))
    if name == "control":
        return out
    if name == "tick_shift":
        dy = int(arg)
        x0, x1 = PANEL["spine_x"] - 5, PANEL["spine_x"] - 1
        strip = out[:, x0:x1].copy()
        r = ticks[0]
        out[r - 3:r + 4, x0:x1] = 254
        out[r + dy - 3:r + dy + 4, x0:x1] = strip[r - 3:r + 4]
        return out
    if name == "category_stub_deleted":
        ct = int(round(cal["category_ticks_px"][2]))
        out[baseline + 1:baseline + 6, ct - 3:ct + 4] = 254
        return out
    if name == "marker":
        w = int(arg)
        d = _diag(img)["Kupffer Cells"]["lnp"]
        cx = (d["x0"] + d["x1"]) // 2
        lo = cx - w // 2
        out[d["row"] - 30:d["row"] - 22, lo:lo + w] = [113, 90, 138]
        return out
    if name == "spine_erased":
        # removes the furniture the 100%-level check leans on
        out[PANEL["y0"] - 2:PANEL["y1"] + 2, PANEL["spine_x"] - 2:
            PANEL["spine_x"] + 3] = 254
        return out
    raise SystemExit("unknown seed " + name)


def _diag(img):
    return read(img)[1]


def main():
    name = sys.argv[1]
    arg = sys.argv[2] if len(sys.argv) > 2 else None
    img, _ = S.load(PANEL["page"], PANEL["xref"])
    bad = mutate(img, name, arg)
    try:
        vals, _, cal = read(bad)
    except (DetectionError, CalibrationError) as e:
        print("REFUSED", type(e).__name__, str(e).splitlines()[0][:160])
        print("METRIC 0")
        sys.exit(2)
    got = [vals[k]["lnp"] for k in sorted(vals)]
    n = sum(1 for v in got if v is not None)
    print("loo_max_px", cal["loo_max_px"],
          "spine_top_residual_px", cal["spine_top_residual_px"])
    print("METRIC", n, " ".join("None" if v is None else "%.3f" % v
                                for v in got))
    sys.exit(0)


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except BaseException as exc:                       # metric never reached
        print("CRASHED", type(exc).__name__, str(exc)[:160])
        sys.exit(3)
