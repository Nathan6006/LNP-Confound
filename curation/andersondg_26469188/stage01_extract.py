#!/usr/bin/env python3
"""KK_2015 (Kauffman 2015, PMID 26469188) -- parse the SI's native text tables.

TIER 4 ONLY. Every value here is a text character inside the PDF's content
stream. No page is rendered; no raster is sampled. `curation/lib`'s axis/bars/
heatmap modules are deliberately NOT imported because there is nothing to
calibrate -- see 00_value_sources.json "pixels".

Writes (into OUTDIR, default = this folder):
    01_blocks.csv        the design factors and their levels
    05_values.csv        every reported value, with corpus_admissible + reason
    _parsed.json         the parsed tables, for stage02/stage03 to reuse
"""
import csv
import json
import os
import re
import sys

import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "nl5b02497_si_001.pdf")
# F-069: outputs go where the caller says, so a seed harness can never write
# into the curation folder. Default is the folder only for a normal run.
OUTDIR = os.environ.get("KK2015_OUTDIR", HERE)

PHOSPH = {"DSPC", "DOPE", "DOPC", "DSPE"}

# ---------------------------------------------------------------- parsing ----

def si_pages():
    with pdfplumber.open(SI) as pdf:
        return [(pg.extract_text() or "") for pg in pdf.pages]


_S1 = re.compile(
    r"^(?P<code>PBS|Original|[ABC]-\d{2})\s+"
    r"(?P<wt>n/a|[\d.]+)\s+"
    r"(?P<ph>n/a|DSPC|DOPE|DOPC|DSPE)\s+"
    r"(?P<il>n/a|[\d.]+)\s+(?P<phm>n/a|[\d.]+)\s+"
    r"(?P<chol>n/a|[\d.]+)\s+(?P<peg>n/a|[\d.]+)\s+"
    r"(?P<ee>n/a|\*|[\d.]+)\s+(?P<size>n/a|\*|[\d.]+)\s+"
    r"(?P<pdi>n/a|\*|[\d.]+)\s+"
    r"(?P<val>\*|[\d.]+\s*±\s*[\d.]+)\s*$"
)


def parse_table_s1(text):
    """SI Table S1, page 2. 41 data lines: PBS + Original + A-01..C-38."""
    rows = []
    for line in text.split("\n"):
        m = _S1.match(line.strip())
        if not m:
            continue
        d = m.groupdict()
        rec = {"code": d["code"]}
        if d["val"] == "*":
            rec["epo_mean"] = rec["epo_sd"] = None
            rec["unsynthesizable"] = True
        else:
            a, b = re.split(r"\s*±\s*", d["val"])
            rec["epo_mean"], rec["epo_sd"] = float(a), float(b)
            rec["unsynthesizable"] = False
        for k, src in (("wt_ratio", "wt"), ("il_molpct", "il"),
                       ("ph_molpct", "phm"), ("chol_molpct", "chol"),
                       ("peg_molpct", "peg"), ("ee_pct", "ee"),
                       ("size_nm", "size"), ("pdi", "pdi")):
            rec[k] = None if d[src] in ("n/a", "*") else float(d[src])
        rec["phospholipid"] = None if d["ph"] == "n/a" else d["ph"]
        rows.append(rec)
    return rows


_DSD = re.compile(
    r"^(?P<code>A-\d{2})\s+(?P<design>[+\-0—]{4,6})\s+"
    r"(?P<head>DO|DS)\s+(?P<tail>PE|PC)\s+"
    r"(?P<wt>[\d.]+)\s+(?P<il>[\d.]+)\s+(?P<phm>[\d.]+)\s+(?P<peg>[\d.]+)\s*$"
)
_FFD = re.compile(
    r"^(?P<code>B-\d{2})\s+(?P<design>[+\-0—]{3,5})\s+"
    r"(?P<ph>DSPC|DOPE)\s+"
    r"(?P<wt>[\d.]+)\s+(?P<il>[\d.]+)\s+(?P<phm>[\d.]+)\s+(?P<peg>[\d.]+)\s*$"
)


LEVEL_CHARS = set("+-0")


def _fix_emdash(s, want):
    """pdfplumber fuses a run of hyphens into U+2014. Expand it to reach `want`.

    This is a TEXT-EXTRACTION artifact, not a data question, and the expansion
    is only accepted when it is FORCED: exactly one em dash, and replacing it
    with (want - len(s) + 1) hyphens is the unique way to reach `want`.

    Its first version had a hole its own seed found: it returned early whenever
    `len(s) == want`, so an em-dash string that HAPPENED to be the right width
    was handed back with the em dash still in it -- and `check_coded_levels`
    would then index that character as if it were a design level, silently
    adding a fourth key to a three-key map and masking a violation. A guard
    that short-circuits on a length test has given up the ability to see a
    non-canonical character. So the character set is now validated on the way
    OUT, unconditionally, and length alone never buys an exit.
    """
    out = s
    if "\u2014" in out:
        if out.count("\u2014") != 1:
            raise ValueError("design string %r has %d em dashes; only a single "
                             "fused hyphen run can be expanded unambiguously"
                             % (s, out.count("\u2014")))
        need = want - (len(out) - 1)
        if need < 2:
            raise ValueError("em-dash expansion of %r to width %d would not "
                             "lengthen it, so the em dash is not a fused hyphen "
                             "run" % (s, want))
        out = out.replace("\u2014", "-" * need)
    if len(out) != want:
        raise ValueError("design string %r normalises to %r, width %d, expected %d"
                         % (s, out, len(out), want))
    bad = sorted(set(out) - LEVEL_CHARS)
    if bad:
        raise ValueError("design string %r normalises to %r which contains "
                         "non-level character(s) %r; a level is one of %s"
                         % (s, out, bad, sorted(LEVEL_CHARS)))
    return out, out != s


def parse_designs(p10, p11):
    out = {}
    fixed = []
    for line in p10.split("\n"):
        m = _DSD.match(line.strip())
        if m:
            d = m.groupdict()
            ds, was = _fix_emdash(d["design"], 6)
            if was:
                fixed.append((d["code"], d["design"], ds))
            out[d["code"]] = {
                "design": ds, "library": "A", "head": d["head"], "tail": d["tail"],
                "phospholipid": d["head"] + d["tail"],
                "wt_ratio": float(d["wt"]), "il_molpct": float(d["il"]),
                "ph_molpct": float(d["phm"]), "peg_molpct": float(d["peg"]),
            }
    for line in p11.split("\n"):
        m = _FFD.match(line.strip())
        if m:
            d = m.groupdict()
            ds, was = _fix_emdash(d["design"], 5)
            if was:
                fixed.append((d["code"], d["design"], ds))
            out[d["code"]] = {
                "design": ds, "library": "B", "head": d["ph"][:2], "tail": d["ph"][2:],
                "phospholipid": d["ph"],
                "wt_ratio": float(d["wt"]), "il_molpct": float(d["il"]),
                "ph_molpct": float(d["phm"]), "peg_molpct": float(d["peg"]),
            }
    return out, fixed


_S3 = re.compile(
    r"^(?P<code>Original|C-35),?\s+(?P<wt>[\d.]+)\s+(?P<ph>DSPC|DOPE)\s+"
    r"(?P<il>[\d.]+)\s+(?P<phm>[\d.]+)\s+(?P<chol>[\d.]+)\s+(?P<peg>[\d.]+)\s+"
    r"(?P<ee>[\d.]+)\s+(?P<size>[\d.]+)\s+(?P<pdi>[\d.]+)\s+"
    r"(?P<val>[\d.]+\s*±\s*[\d.]+)\s*$"
)


def parse_table_s3(text):
    """Luciferase mRNA. The code cell wraps ('Original,' / 'Luc'), so the
    numeric tail lands on the line that starts with the code."""
    rows = []
    for line in text.split("\n"):
        m = _S3.match(line.strip())
        if not m:
            continue
        d = m.groupdict()
        a, b = re.split(r"\s*±\s*", d["val"])
        rows.append({
            "code": d["code"], "wt_ratio": float(d["wt"]), "phospholipid": d["ph"],
            "il_molpct": float(d["il"]), "ph_molpct": float(d["phm"]),
            "chol_molpct": float(d["chol"]), "peg_molpct": float(d["peg"]),
            "ee_pct": float(d["ee"]), "size_nm": float(d["size"]),
            "pdi": float(d["pdi"]),
            "flux_mean_e9": float(a), "flux_sd_e9": float(b),
        })
    return rows


def parse_table_s4(pg):
    """FVII siRNA. Table S4's READING ORDER IS NOT ITS VISUAL ORDER: the dose
    cell is vertically centred, so each formulation's FIRST dose line
    ('0.1 mg/kg: ...') sorts ABOVE its own code label. A line-order parser
    silently drops those two lines (measured: 6 of 8 value rows recovered).

    So the three dose lines are assigned to a formulation by NATIVE TEXT
    y-coordinate -- the PDF's own character positions, not pixels. Each
    formulation is anchored on its phospholipid token, and every dose line goes
    to the nearest anchor. The assignment margin is asserted and reported, so a
    layout this does not fit raises instead of guessing.
    """
    words = pg.extract_words()
    text = pg.extract_text() or ""

    anchors = {}
    for w in words:
        if w["text"] in ("DSPC", "DOPE"):
            anchors[w["text"]] = w["top"]
    if set(anchors) != {"DSPC", "DOPE"}:
        raise ValueError("Table S4: expected one DSPC and one DOPE anchor, got %r"
                         % sorted(anchors))
    code_of = {"DSPC": "Original", "DOPE": "C-35"}

    # the formulation header cells, matched on their own line
    head_re = re.compile(
        r"^(?P<wt>[\d.]+)\s+(?P<ph>DSPC|DOPE)\s+(?P<il>[\d.]+)\s+(?P<phm>[\d.]+)\s+"
        r"(?P<chol>[\d.]+)\s+(?P<peg>[\d.]+)\s+(?P<ee>[\d.]+)\s+(?P<size>[\d.]+)\s+"
        r"(?P<pdi>[\d.]+)\s+[\d.]+\s*mg/kg:")
    formul = {}
    for line in text.split("\n"):
        m = head_re.match(line.strip())
        if m:
            d = m.groupdict()
            formul[code_of[d["ph"]]] = {
                "wt_ratio": float(d["wt"]), "phospholipid": d["ph"],
                "il_molpct": float(d["il"]), "ph_molpct": float(d["phm"]),
                "chol_molpct": float(d["chol"]), "peg_molpct": float(d["peg"]),
                "ee_pct": float(d["ee"]), "size_nm": float(d["size"]),
                "pdi": float(d["pdi"])}
    if set(formul) != {"Original", "C-35"}:
        raise ValueError("Table S4: formulation header cells not found: %r" % sorted(formul))

    # every 'X mg/kg: m +/- s' run, located by the y of its 'mg/kg:' token
    doses, margins = [], []
    for w in words:
        if w["text"] != "mg/kg:":
            continue
        y = w["top"]
        band = [x for x in words if abs(x["top"] - y) < 2.0]
        band.sort(key=lambda x: x["x0"])
        frag = " ".join(x["text"] for x in band)
        m = re.search(r"([\d.]+)\s*mg/kg:\s*([\d.]+)\s*\u00b1\s*([\d.]+)", frag)
        if not m:
            raise ValueError("Table S4: unparsed dose band at y=%.2f: %r" % (y, frag))
        d = sorted(anchors.items(), key=lambda kv: abs(kv[1] - y))
        near, far = d[0], d[1]
        margins.append((abs(far[1] - y), abs(near[1] - y)))
        doses.append({"code": code_of[near[0]], "dose_mg_kg": float(m.group(1)),
                      "fvii_mean": float(m.group(2)), "fvii_sd": float(m.group(3))})

    for code in ("Original", "C-35"):
        got = sorted(r["dose_mg_kg"] for r in doses if r["code"] == code)
        if got != [0.01, 0.03, 0.1]:
            raise ValueError("Table S4: %s got doses %r, expected [0.01, 0.03, 0.1]"
                             % (code, got))

    rows = []
    m = re.search(r"^PBS\s+(?:n/a\s+){9}([\d.]+)\s*\u00b1\s*([\d.]+)$",
                  text, re.M)
    if not m:
        raise ValueError("Table S4: PBS row not found")
    rows.append({"code": "PBS", "dose_mg_kg": None,
                 "fvii_mean": float(m.group(1)), "fvii_sd": float(m.group(2))})
    for r in doses:
        rec = dict(formul[r["code"]])
        rec.update(r)
        rows.append(rec)

    worst = min(margins, key=lambda t: t[0] / max(t[1], 1e-9))
    return rows, {"n_dose_bands": len(doses),
                  "anchor_separation_pt": round(abs(anchors["DOPE"] - anchors["DSPC"]), 2),
                  "worst_assignment_far_pt": round(worst[0], 2),
                  "worst_assignment_near_pt": round(worst[1], 2),
                  "worst_margin_ratio": round(worst[0] / max(worst[1], 1e-9), 2)}


def parse_table_s2(text):
    """Printed regression statistics used as external checks."""
    g = {}
    for pat, key in ((r"^R2\s+([\d.]+)$", "r2"),
                     (r"^Root Mean Square Error\s+([\d.]+)$", "rmse"),
                     (r"^Mean\s+([\d.]+)$", "mean"),
                     (r"^n\s+(\d+)$", "n")):
        for line in text.split("\n"):
            m = re.match(pat, line.strip())
            if m:
                g[key] = float(m.group(1))
                break
    m = re.search(r"Phospholipid\[DOPE\]\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", text)
    if m:
        g["dope_estimate"] = float(m.group(1))
        g["dope_stderr"] = float(m.group(2))
        g["dope_tratio"] = float(m.group(3))
    return g


# ----------------------------------------------------------------- output ----

def write_blocks(path):
    """01_blocks.csv -- the DESIGN FACTORS. This paper's 'blocks' are not
    chemical building blocks of a combinatorial library: there is one ionizable
    lipid and the library is a formulation design. Stating that here, in the
    file a reader opens first, is the point."""
    rows = [
        ("ionizable_lipid", "C12-200", 1,
         "FIXED across every formulation in the paper. This is the whole answer to "
         "'how many ionizable lipids': one.", "yes"),
        ("phospholipid", "DSPC | DOPE | DOPC | DSPE", 4,
         "The only CHEMICAL identity that varies. DSPE appears in 3 designs "
         "(A-07, A-12, A-14) that could not be synthesized, so only 3 of the 4 "
         "are ever measured.", "yes"),
        ("sterol", "Cholesterol", 1,
         "FIXED. Written only as 'Chol.' in Table S1; never drawn and never "
         "given a systematic name.", "no"),
        ("peg_lipid", "UNSTATED", 1,
         "FIXED but NEVER NAMED anywhere in the SI -- the column header is the "
         "bare word 'PEG'. The main text is not held (closed access), so the "
         "identity is unavailable. Emitted as PEG_name='Unknown'.", "no"),
        ("il_to_mrna_weight_ratio", "2.5 | 5 | 7.5 | 10 | 12.5 | 15 | 20 | 25", 8,
         "Quantitative factor. 3 levels per library; Library C extends it to 25.",
         "yes"),
        ("il_molpct", "30 | 35 | 40 | 50 | 60", 5, "Quantitative factor.", "yes"),
        ("phospholipid_molpct", "4 | 10 | 16 | 22 | 28", 5, "Quantitative factor.", "yes"),
        ("peg_molpct", "0.5 | 1.5 | 2.5 | 3 | 3.5", 5, "Quantitative factor.", "yes"),
        ("cholesterol_molpct", "balance to 100", 0,
         "NOT an independent factor: it is 100 minus the other three, so it is a "
         "dependent column. Checking that the four sum to 100 therefore tests OUR "
         "TRANSCRIPTION, not the authors' design -- see 07_qc.json.", "yes"),
    ]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["factor", "levels", "n_levels", "note", "varies_across_arms"])
        w.writerows(rows)
    return len(rows)


def main():
    pages = si_pages()
    s1 = parse_table_s1(pages[1])
    designs, emdash_fixed = parse_designs(pages[9], pages[10])
    s2 = parse_table_s2(pages[2])
    s3 = parse_table_s3(pages[3])
    with pdfplumber.open(SI) as pdf:
        s4, s4geom = parse_table_s4(pdf.pages[4])

    assert len(s1) == 40, "Table S1: expected 40 data lines (PBS + Original + 38 library codes), got %d" % len(s1)
    assert len(designs) == 32, "design tables: expected 32 rows, got %d" % len(designs)
    assert len(s3) == 2, "Table S3: expected 2 rows, got %d" % len(s3)
    assert len(s4) == 7, "Table S4: expected 7 value rows, got %d" % len(s4)

    os.makedirs(OUTDIR, exist_ok=True)
    nb = write_blocks(os.path.join(OUTDIR, "01_blocks.csv"))

    out = {"table_s1": s1, "designs": designs, "table_s2": s2,
           "table_s3": s3, "table_s4": s4, "emdash_fixed": emdash_fixed,
           "table_s4_geometry": s4geom}
    with open(os.path.join(OUTDIR, "_parsed.json"), "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)

    print("Table S1 rows          : %d (PBS + Original + 38 library codes)" % len(s1))
    codes = {r["code"] for r in s1}
    want = {"PBS", "Original"} | {"A-%02d" % i for i in range(1, 15)} | {"B-%02d" % i for i in range(15, 33)} | {"C-%02d" % i for i in range(33, 39)}
    assert codes == want, "Table S1 codes differ: missing %s extra %s" % (sorted(want - codes), sorted(codes - want))
    print("  with an EPO value    : %d" % sum(1 for r in s1 if r["epo_mean"] is not None))
    print("  unsynthesizable (*)  : %d  %s"
          % (sum(1 for r in s1 if r["unsynthesizable"]),
             [r["code"] for r in s1 if r["unsynthesizable"]]))
    print("design-table rows      : %d" % len(designs))
    print("em-dash normalisations : %d  %s" % (len(emdash_fixed), emdash_fixed))
    print("Table S2 statistics    : %s" % s2)
    print("Table S3 rows          : %d" % len(s3))
    print("Table S4 value rows    : %d  geometry %s" % (len(s4), s4geom))
    print("01_blocks.csv factors  : %d" % nb)
    return 0


if __name__ == "__main__":
    sys.exit(main())
