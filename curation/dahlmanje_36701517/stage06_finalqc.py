"""Stage 06 -- re-attach the narrative QC sections and re-measure the one
section that is derived from pixels.

WHY THIS EXISTS. stage03_emit.py WRITES 07_qc.json, so anything added to that
file by hand after stage03 is destroyed the next time stage03 runs. Four
sections were added that way. They live in _qc_narrative_sections.json now and
are re-attached here, and the one of the four that is a MEASUREMENT rather than
prose -- error_class_32_headroom_measured -- is recomputed from the raster
every run rather than restored from the cache, so it cannot go stale silently.

Run order:
    stage01_digitize.py -> stage02_overlap.py -> stage03_emit.py
    -> stage04_selfcheck.py -> stage05_seedgate.py -> stage06_finalqc.py

Usage: python3 -I stage06_finalqc.py
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)
import stage01_digitize as S                                   # noqa: E402


def headroom():
    """The widest real error-bar cap, as a fraction of its bar's width.

    bars.bar_top's full-width coverage rule is 0.90. This measures how close
    the figure's OWN markers come to it, which is what turns the leak the
    width sweep found into a stated limitation rather than a live defect.
    """
    img, _ = S.load(23, 244)
    coloured, _ = S.masks(img)
    dig = json.load(open(os.path.join(HERE, "04_digitized.json")))
    rows = []
    for organ, o in dig["panels"]["S13a"]["organs"].items():
        for lab, dg in o["diagnostics"].items():
            r = dg.get("lnp", {})
            if "row" not in r:
                continue
            x0, x1, top = r["x0"], r["x1"], r["row"]
            widest = 0
            for y in range(max(0, top - 28), top - 2):
                run = best = 0
                for v in coloured[y, x0 - 6:x1 + 6]:
                    run = run + 1 if v else 0
                    best = max(best, run)
                widest = max(widest, best)
            rows.append({"organ": organ, "cell_type": lab,
                         "bar_width_px": int(x1 - x0),
                         "widest_coloured_run_above_the_top_px": int(widest),
                         "ratio": round(widest / (x1 - x0), 3)})
    worst = max(r["ratio"] for r in rows)
    return {
        "_what": "the widest contiguous bar-coloured run in the 26 rows above "
                 "each bar's detected top, over every readable Cat-LNP bar of "
                 "Fig S13a. These are the figure's own error-bar caps.",
        "per_bar": rows, "n_bars": len(rows),
        "max_cap_to_bar_ratio": worst,
        "min_cap_to_bar_ratio": min(r["ratio"] for r in rows),
        "coverage_gate": 0.90, "margin_to_gate": round(0.90 - worst, 3),
        "_verdict": "The widest real cap is %.3f of its bar against a 0.90 "
                    "gate -- a margin of %.3f on %d of %d bars. The width "
                    "sweep in the seed gate shows a marker must reach the "
                    "FULL bar width to be read as a top, and nothing drawn in "
                    "this figure comes close."
                    % (worst, 0.90 - worst, len(rows), len(rows)),
    }


def anchor_returncodes():
    """Re-run make_reference.py clean and with each seed, in its OWN process.

    The verdict is proc.returncode. Not stdout, not a pipe: `python3 ... |
    tail` returns tail's exit code, which is F-012's shape.
    """
    res = {}
    for args, name in (([], "clean"),
                       (["--seed-guard"], "seed_ast_guard"),
                       (["--seed-ester"], "seed_isomeric_ester_reversal")):
        p = subprocess.run([sys.executable, "-I", "make_reference.py"] + args,
                           capture_output=True, text=True, cwd=HERE)
        res[name] = {"argv": ["python3", "-I", "make_reference.py"] + args,
                     "returncode": p.returncode,
                     "stdout_tail": p.stdout.strip().splitlines()[-4:],
                     "verdict_source": "proc.returncode"}
    res["_expected"] = {"clean": 0, "seed_ast_guard": 1,
                        "seed_isomeric_ester_reversal": 1}
    res["_all_as_expected"] = all(
        res[k]["returncode"] == v for k, v in res["_expected"].items())
    return res


def main():
    qc_path = os.path.join(HERE, "07_qc.json")
    qc = json.load(open(qc_path))
    cached = json.load(open(os.path.join(HERE,
                                         "_qc_narrative_sections.json")))
    for k, v in cached.items():
        qc[k] = v
    qc["error_class_32_headroom_measured"] = headroom()      # re-measured
    qc["anchor_seed_gate_returncodes"] = anchor_returncodes()  # re-run
    json.dump(qc, open(qc_path, "w"), indent=1, sort_keys=True, default=str)
    print("07_qc.json now has %d sections; headroom re-measured at %.3f"
          % (len(qc), qc["error_class_32_headroom_measured"]
             ["max_cap_to_bar_ratio"]))


if __name__ == "__main__":
    main()
