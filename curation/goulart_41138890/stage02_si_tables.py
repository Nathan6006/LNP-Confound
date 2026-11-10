#!/usr/bin/env python3
"""Tables S1 and S3 out of mmc1.docx, SYMBOLICALLY.

The .docx holds exactly ONE <w:tbl> (Table S4, pairwise statistics). Tables S1,
S2 and S3 are EMF images. A <w:tbl> count therefore reports one table for a
document holding four -- the brief's "read the primary object, not a proxy".

But an EMF is a Windows METAFILE, not a bitmap. image15.emf (Table S1) holds
694 EMR_EXTTEXTOUTW records, image17.emf (Table S3) 534: one Unicode glyph each
with an exact device-space reference point. So every digit below is a
TRANSCRIPTION of a text record, not a pixel measurement, and the only inference
is the grouping -- rows by equal y, cells by an x gap larger than any
within-cell advance.

The grouping threshold is MEASURED, not chosen: the within-line glyph advances
and the between-column gaps are printed, and the script asserts they are
cleanly separated before it groups anything.

MS-EMF RecordType numbering above 0x51 is off-by-one from the naive guess and
getting it wrong decodes EXTCREATEFONTINDIRECTW records as text (41 records of
binary garbage on the first attempt). Verified: STRETCHDIBITS=81,
EXTCREATEFONTINDIRECTW=82, EXTTEXTOUTA=83, EXTTEXTOUTW=84, POLYBEZIER16=85.
"""
import argparse
import csv
import hashlib
import json
import os
import struct
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DOCX = os.path.join(HERE, "source", "mmc1.docx")
DOCX_SHA256 = ("b4f0839295ca03107646ff2080a0f85d313b76238c11c6444afb9005"
               "bca977ab")
# image15 = the image immediately preceded by the "Table S1." caption,
# image17 = by "Table S3." -- established by walking word/document.xml's body
# in order and mapping each a:blip r:embed through word/_rels/document.xml.rels.
TABLES = {"S1": "word/media/image15.emf", "S3": "word/media/image17.emf"}
EXTTEXTOUTA, EXTTEXTOUTW = 83, 84
EMR_EOF = 14


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def emf_texts(data):
    """(x, y, string) for every EMR_EXTTEXTOUT record, in file order."""
    off, out = 0, []
    while off + 8 <= len(data):
        typ, size = struct.unpack_from("<II", data, off)
        if size < 8 or off + size > len(data):
            raise ValueError(f"EMF record stream desynchronised at {off}")
        if typ in (EXTTEXTOUTA, EXTTEXTOUTW):
            # Type(4) Size(4) Bounds(16) iGraphicsMode(4) exScale(4) eyScale(4)
            # then EmrText{ Reference(8) Chars(4) offString(4) Options(4) ... }
            rx, ry, chars, off_s, _ = struct.unpack_from("<iiIII", data,
                                                         off + 36)
            raw = data[off + off_s: off + off_s
                       + (chars * 2 if typ == EXTTEXTOUTW else chars)]
            s = (raw.decode("utf-16-le") if typ == EXTTEXTOUTW
                 else raw.decode("latin-1"))
            out.append((rx, ry, s))
        if typ == EMR_EOF:
            break
        off += size
    return out


def group(texts):
    """Rows by equal y; cells by an x gap wider than any within-cell advance.

    Returns (rows, diagnostics). Raises if the two gap populations are not
    cleanly separated, because then no threshold is defensible.
    """
    ys = sorted({t[1] for t in texts})
    gaps = []
    for y in ys:
        row = sorted([t for t in texts if t[1] == y])
        gaps += [b[0] - a[0] for a, b in zip(row, row[1:])]
    gaps = sorted(gaps)
    # The separation is the widest glyph advance INSIDE a cell vs the narrowest
    # between-column gap. It is the FIRST clean break scanning upward, not the
    # largest break overall: the largest break in this distribution is between
    # two column gaps (1000 -> 2214 on Table S1), and taking it collapses every
    # row into one cell. Criterion: the first adjacent pair whose ratio is
    # >= 1.5 AND whose absolute step is >= 100 device units.
    idx = None
    for i in range(len(gaps) - 1):
        if gaps[i + 1] >= 1.5 * gaps[i] and gaps[i + 1] - gaps[i] >= 100:
            idx = i
            break
    if idx is None:
        raise ValueError(
            "group: no adjacent pair in the glyph-gap distribution is both "
            ">=1.5x and >=100 device units apart. Within-cell advances and "
            "between-column gaps are not separable, so no threshold is "
            "defensible. REFUSAL.")
    lo_max, hi_min = gaps[idx], gaps[idx + 1]
    best_jump = hi_min - lo_max
    thresh = (lo_max + hi_min) / 2.0
    rows = []
    for y in ys:
        row = sorted([t for t in texts if t[1] == y])
        cells, cur = [], [row[0]]
        for a, b in zip(row, row[1:]):
            if b[0] - a[0] > thresh:
                cells.append(cur)
                cur = [b]
            else:
                cur.append(b)
        cells.append(cur)
        rows.append([(c[0][0], "".join(g[2] for g in c).strip())
                     for c in cells])
    diag = {"n_text_records": len(texts), "n_rows": len(ys),
            "gap_ratio_at_break": round(hi_min / lo_max, 3),
            "largest_gap_jump": best_jump,
            "widest_within_cell_advance": lo_max,
            "narrowest_column_gap": hi_min,
            "threshold_used": thresh}
    return rows, diag


def num(s):
    """The SI writes decimals with a COMMA (Brazilian locale): '47,5', '0,113'.
    A naive float() would raise; a naive strip of commas would turn 47,5 into
    475. Converted explicitly, and '-' stays None (not 0)."""
    s = s.strip()
    if s in ("-", ""):
        return None
    return float(s.replace(",", "."))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=HERE)
    ap.add_argument("--seed", default=None)
    a = ap.parse_args()
    if a.seed and os.path.abspath(a.outdir) == os.path.abspath(HERE):
        sys.exit("REFUSED: a seeded run must write to a temp --outdir (F-069).")
    os.makedirs(a.outdir, exist_ok=True)

    got = sha256(DOCX)
    if got != DOCX_SHA256:
        sys.exit(f"mmc1.docx sha256 {got} != pinned {DOCX_SHA256}")

    out, diags = {}, {}
    with zipfile.ZipFile(DOCX) as z:
        for tag, member in TABLES.items():
            data = z.read(member)
            if a.seed == "swap_digit" and tag == "S1":
                data = data.replace("4".encode("utf-16-le"),
                                    "9".encode("utf-16-le"), 1)
            rows, diag = group(emf_texts(data))
            out[tag] = rows
            diags[tag] = diag

    # ---------- Table S1: 18 b-DNA LNPs ----------
    s1 = []
    for r in out["S1"]:
        if not r or not r[0][1].startswith("LNP"):
            continue
        cells = [c[1] for c in r]
        # columns: name | IL | DOPE | DSPC | DOTAP | Chol | C14PEG | d | PDI | zeta
        assert len(cells) == 10, (len(cells), cells)
        helper = [("DOPE", cells[2]), ("DSPC", cells[3]), ("DOTAP", cells[4])]
        named = [(n, v) for n, v in helper if v.strip() != "-"]
        assert len(named) == 1, (cells,)
        s1.append({"formulation": cells[0],
                   "IL_molratio": num(cells[1]),
                   "helper_lipid": named[0][0],
                   "HL_molratio": num(named[0][1]),
                   "CHL_molratio": num(cells[5]),
                   "PEG_molratio": num(cells[6]),
                   "diameter_nm": num(cells[7]),
                   "PDI": num(cells[8]),
                   "zeta_mV": num(cells[9])})

    # ---------- Table S3: 11 pDNA-GFP LNPs ----------
    s3 = []
    for r in out["S3"]:
        if not r or not r[0][1].startswith("LNP"):
            continue
        cells = [c[1] for c in r]
        assert len(cells) == 11, (len(cells), cells)
        helper = [("DOPE", cells[2]), ("DSPC", cells[3]), ("DOTAP", cells[4])]
        named = [(n, v) for n, v in helper if v.strip() != "-"]
        assert len(named) == 1, (cells,)
        s3.append({"formulation": cells[0],
                   "IL_molratio": num(cells[1]),
                   "helper_lipid": named[0][0],
                   "HL_molratio": num(named[0][1]),
                   "CHL_molratio": num(cells[5]),
                   "PEG_molratio": num(cells[6]),
                   "diameter_nm": num(cells[7]),
                   "PDI": num(cells[8]),
                   "zeta_mV": num(cells[9]),
                   "pKa": num(cells[10])})

    # ---------- checks ----------
    checks = {}
    checks["S1_n_formulations"] = len(s1)
    checks["S3_n_formulations"] = len(s3)
    # the factorial the design claims: 3 helpers x 6 contents = 18 = row count
    helpers = sorted({r["helper_lipid"] for r in s1})
    contents = sorted({r["HL_molratio"] for r in s1})
    checks["S1_helpers"] = helpers
    checks["S1_helper_contents"] = contents
    checks["S1_factorial"] = (f"{len(helpers)} x {len(contents)} = "
                              f"{len(helpers) * len(contents)}")
    checks["S1_factorial_equals_row_count"] = (
        len(helpers) * len(contents) == len(s1))
    for h in helpers:
        got_c = sorted(r["HL_molratio"] for r in s1
                       if r["helper_lipid"] == h)
        checks[f"S1_contents_{h}"] = got_c
    # mol% sum. Cholesterol is almost certainly the authors' BALANCE, so this
    # is a weak check on its own (brief section 10) -- but all four percents
    # are PRINTED independently, so it is still a real transcription check.
    sums = [round(r["IL_molratio"] + r["HL_molratio"] + r["CHL_molratio"]
                  + r["PEG_molratio"], 6) for r in s1]
    checks["S1_molpct_sums"] = sorted(set(sums))
    checks["S1_molpct_sum_exact_100"] = sum(1 for v in sums if v == 100.0)
    checks["S1_molpct_sum_n"] = len(sums)
    checks["S1_molpct_sum_is_tautological"] = (
        "PARTLY: cholesterol is the residual by construction, so the identity "
        "can only fail on a transcription error in one of the four printed "
        "percents -- which is what it is here for. It is NOT evidence the "
        "design is correct.")
    sums3 = [round(r["IL_molratio"] + r["HL_molratio"] + r["CHL_molratio"]
                   + r["PEG_molratio"], 6) for r in s3]
    checks["S3_molpct_sum_exact_100"] = sum(1 for v in sums3 if v == 100.0)
    checks["S3_molpct_sum_n"] = len(sums3)
    checks["IL_molratio_constant"] = sorted(
        {r["IL_molratio"] for r in s1 + s3})
    checks["PEG_molratio_constant"] = sorted(
        {r["PEG_molratio"] for r in s1 + s3})
    # S3 must be a SUBSET of S1's formulations, by name AND by composition
    byname = {r["formulation"]: r for r in s1}
    agree, disagree = 0, []
    for r in s3:
        p = byname.get(r["formulation"])
        if p is None:
            disagree.append((r["formulation"], "absent from S1"))
            continue
        same = all(p[k] == r[k] for k in ("IL_molratio", "helper_lipid",
                                          "HL_molratio", "CHL_molratio",
                                          "PEG_molratio"))
        if same:
            agree += 1
        else:
            disagree.append((r["formulation"], "composition differs"))
    checks["S3_subset_of_S1_agree"] = agree
    checks["S3_subset_of_S1_disagree"] = disagree
    checks["glyph_grouping"] = diags

    with open(os.path.join(a.outdir, "02_si_tables.json"), "w") as fh:
        json.dump({"tableS1": s1, "tableS3": s3, "checks": checks}, fh,
                  indent=1)
    for tag, rowset in (("S1", s1), ("S3", s3)):
        path = os.path.join(a.outdir, f"02_table{tag}.csv")
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rowset[0].keys()))
            w.writeheader()
            w.writerows(rowset)
    print(json.dumps(checks, indent=1))
    print(f"wrote 02_si_tables.json, 02_tableS1.csv ({len(s1)}), "
          f"02_tableS3.csv ({len(s3)})")


if __name__ == "__main__":
    main()
