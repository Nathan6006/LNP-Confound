#!/usr/bin/env python3
"""Stage 01 - Supplementary Table 3 (formulation details of the tested LNP library).

TIER 4 (pdf text table). The SI PDF carries this table as NATIVE TEXT, not as an
image: `page.get_images()` on SI pp. 11-13 returns only the 295x88 running
header logo, and every digit below is a text span. So no pixel is read here.

Words are taken WITH COORDINATES and clustered into rows by y, rather than
trusting the linear text stream, because a stream order is a rendering artefact
and a single mis-ordered span would silently shift a whole row.
"""
import csv, json, os, sys
import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "adma72049-sup-0001-SuppMat.pdf")
PAGES = [10, 11, 12]          # 0-based; SI pp. 11-13
HEADER = ["Index", "Ionizable lipid", "%mol", "Helper lipid", "%mol",
          "Cholesterol (%mol)", "PEG-lipid", "%mol"]


def rows_from_page(page):
    words = page.get_text("words")          # x0,y0,x1,y1,word,block,line,wordno
    words = [w for w in words if w[4].strip()]
    # cluster by y centre
    lines = {}
    for w in words:
        yc = round((w[1] + w[3]) / 2, 1)
        key = min(lines, key=lambda k: abs(k - yc)) if lines else None
        if key is None or abs(key - yc) > 3.0:
            key = yc
            lines.setdefault(key, [])
        lines[key].append(w)
    out = []
    for yc in sorted(lines):
        ws = sorted(lines[yc], key=lambda w: w[0])
        toks = [w[4] for w in ws]
        xs = [w[0] for w in ws]
        out.append((yc, toks, xs))
    return out


def main():
    doc = pymupdf.open(SI)
    recs, skipped = [], []
    for pno in PAGES:
        for yc, toks, xs in rows_from_page(doc[pno]):
            # a data row is exactly 8 tokens whose first is an integer index
            if len(toks) == 8 and toks[0].isdigit():
                recs.append({
                    "Index": int(toks[0]), "IL_name": toks[1],
                    "IL_molpct": float(toks[2]),
                    "HL_name": toks[3], "HL_molpct": float(toks[4]),
                    "CHL_molpct": float(toks[5]),
                    "PEG_name": toks[6], "PEG_molpct": float(toks[7]),
                    "si_page": pno + 1, "y": yc,
                })
            else:
                skipped.append({"page": pno + 1, "y": yc, "toks": toks})
    # --- gates -------------------------------------------------------
    idx = [r["Index"] for r in recs]
    assert idx == sorted(idx), "Index column is not ascending in file order"
    assert idx == list(range(1, len(idx) + 1)), (
        "Index is not contiguous 1..N: gaps at %s"
        % [a for a, b in zip(idx, range(1, len(idx) + 1)) if a != b][:5])
    names = [r["IL_name"] for r in recs]
    assert len(set(names)) == len(names), (
        "duplicate ionizable lipid name: %s"
        % [n for n in names if names.count(n) > 1][:5])
    # trap 8: mol% must sum to 100 on every formulation
    bad = []
    for r in recs:
        s = r["IL_molpct"] + r["HL_molpct"] + r["CHL_molpct"] + r["PEG_molpct"]
        r["molpct_sum"] = round(s, 6)
        if abs(s - 100.0) > 1e-9:
            bad.append((r["Index"], r["IL_name"], s))
    # every non-data line must be a header/caption/running head, never a digit row
    for s in skipped:
        assert not (s["toks"] and s["toks"][0].isdigit() and len(s["toks"]) > 4), \
            "a numeric-looking line was dropped: %r" % (s,)

    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "Index", "IL_name", "IL_molpct", "HL_name", "HL_molpct",
            "CHL_molpct", "PEG_name", "PEG_molpct", "molpct_sum", "si_page", "y"])
        w.writeheader()
        w.writerows(recs)

    rep = {
        "n_formulations": len(recs),
        "index_contiguous_1_to_N": True,
        "distinct_IL_names": len(set(names)),
        "molpct_sums_to_100": {"exact": len(recs) - len(bad),
                               "offenders": bad},
        "distinct_formulation_recipes": sorted(set(
            (r["IL_molpct"], r["HL_name"], r["HL_molpct"], r["CHL_molpct"],
             r["PEG_name"], r["PEG_molpct"]) for r in recs)),
        "lines_not_parsed_as_data": [s["toks"] for s in skipped],
    }
    print(json.dumps(rep, indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
