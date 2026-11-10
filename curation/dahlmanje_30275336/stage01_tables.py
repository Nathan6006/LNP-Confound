"""Extract the four LNP library formulation tables from the SI Appendix.

Tier 4 (pdfplumber native text cells). No pixels. Writes 01_formulations.csv.

The free external check this paper offers: the four mol% columns are a
composition and MUST sum to exactly 100. Nothing in the paper prints that sum.
"""
import sys, csv, json
import pdfplumber

PDF = sys.argv[1]
OUT = sys.argv[2]

# (page, table_index, library_label, has_helper_lipid_column, expected_last_lnp)
SPECS = [
    (4,  2, "invitro_54",  False, 54),
    (7,  0, "screen1_112", False, 112),
    (10, 0, "screen2",     False, None),
    (11, 0, "screen3",     True,  None),
]

def reflow(rows):
    """Page 11 packs several LNPs per visual row, newline-separated inside each
    cell. Split every cell on newline and emit one record per sub-line.
    A row whose cells disagree on sub-line count is REFUSED, not guessed."""
    out = []
    for r in rows:
        parts = [(c or "").split("\n") for c in r]
        n = max(len(p) for p in parts)
        if n == 1:
            out.append([(c or "").strip() for c in r])
            continue
        widths = {len(p) for p in parts if any(x.strip() for x in p)}
        if widths - {1, n}:
            raise ValueError("ragged multi-line row, refusing to guess: %r" % (r,))
        for i in range(n):
            out.append([(p[i] if len(p) == n else (p[0] if len(p) == 1 else "")).strip()
                        for p in parts])
    return out

records = []
meta = []
with pdfplumber.open(PDF) as pdf:
    for page, ti, label, has_hl, expect_last in SPECS:
        tbls = pdf.pages[page - 1].find_tables()
        rows = reflow(tbls[ti].extract())
        header = rows[0]
        # the last grid row on pages 7/10/11 is a trailing artifact whose every
        # cell is the literal "n"; it is not an LNP.
        body = [r for r in rows[1:] if r[0].isdigit() and r[1] != "n"]
        for r in body:
            d = dict(zip(header, r))
            rec = {
                "library": label,
                "si_page": page,
                "lnp_number": int(r[0]),
                "compound": d.get("Compound", ""),
                "peg_type": d.get("PEG Type", ""),
                "mass_ratio_compound_to_nucleicacid": d.get("Mass Ratio", ""),
                "helper_lipid_type": d.get("Helper Lipid Type", "DOPE" if not has_hl else ""),
                "compound_molpct": d.get("Compound Mole %", ""),
                "cholesterol_molpct": d.get("Cholesterol Mole %", ""),
                "peg_molpct": d.get("PEG Mole %", ""),
                # the 54-LNP in vitro library is a 3-component formulation: it has
                # no helper-lipid column at all, which is 0 mol%, not "unknown".
                "helper_molpct": d.get("Helper Lipid Mole %", d.get("DOPE Mole %", "0")),
            }
            records.append(rec)
        meta.append({"library": label, "si_page": page, "header": header,
                     "n_rows": len(body),
                     "lnp_numbers_min": min(int(r[0]) for r in body),
                     "lnp_numbers_max": max(int(r[0]) for r in body),
                     "contiguous_1_to_max": sorted(int(r[0]) for r in body) ==
                                            list(range(1, max(int(r[0]) for r in body) + 1)),
                     "expected_last_lnp_from_article_text": expect_last})

# ---- the free check: mol% must sum to exactly 100 ----
def molsum(rec):
    return sum(float(rec[k]) for k in
               ("compound_molpct", "cholesterol_molpct", "peg_molpct", "helper_molpct"))

sums = {}
bad = []
for rec in records:
    try:
        s = round(molsum(rec), 6)
    except ValueError:
        bad.append((rec["library"], rec["lnp_number"], "non-numeric"))
        continue
    rec["molpct_sum"] = s
    sums.setdefault(rec["library"], []).append(s)
    if s != 100.0:
        bad.append((rec["library"], rec["lnp_number"], s))

with open(OUT, "w", newline="") as fh:
    cols = list(records[0].keys())
    for extra in ("molpct_sum",):
        if extra not in cols:
            cols.append(extra)
    w = csv.DictWriter(fh, fieldnames=cols, restval="")
    w.writeheader()
    w.writerows(records)

report = {
    "n_formulations_total": len(records),
    "per_library": meta,
    "compounds_seen": sorted({r["compound"] for r in records}),
    "peg_types_seen": sorted({r["peg_type"] for r in records}),
    "helper_types_seen": sorted({r["helper_lipid_type"] for r in records}),
    "molpct_sums_to_100": {
        "n_checked": sum(len(v) for v in sums.values()),
        "n_exactly_100": sum(1 for v in sums.values() for s in v if s == 100.0),
        "n_off": len(bad),
        "offenders": bad[:20],
        "per_library_distinct_sums": {k: sorted(set(v)) for k, v in sums.items()},
    },
}
print(json.dumps(report, indent=1))
