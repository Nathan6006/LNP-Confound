#!/usr/bin/env python3
"""Tables S1 and S2 and the Figure S7 printed pKa labels -- tier 4, PDF text.

NEITHER TABLE CARRIES A DELIVERY VALUE. Both are physicochemical. They are
extracted because they supply the only in-paper cross-check available on
anything, and because of the trap below.

TWO pKa CONVENTIONS IN ONE SI -- the pKa analogue of the two-mass-conventions
trap. Table S1's "pKa" is the LNP SURFACE pKa, measured on assembled
nanoparticles by TNS fluorescence titration in aqueous buffer (SI p.3: "LNP
surface pKa was determined as previously described [20] ... the pH at which
50% amine protonation occurred"). Table S2's "Average pKa" is the LIPIDOID pKa,
measured on free lipidoid dissolved in ETHANOL and titrated with 0.1 M HCl
(SI p.2 and pp.16-17, Henderson-Hasselbalch from the titration curve). They are
different quantities on different materials in different solvents and they
DISAGREE by up to 0.76 pH units. Neither may be pooled with the other and
neither is a delivery value. Stated here so a later reader cannot discover the
gap by averaging them.
"""
import csv
import json
import os
import re

import pymupdf

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "smll201805097-sup-0001-S1.pdf")

ORDER = ["306O6", "306O8", "306O9", "306Oi10", "306O10", "306O11", "306O12",
         "306O13", "306O14", "306O16", "306O18"]


def page_lines(n):
    doc = pymupdf.open(PDF)
    t = doc[n - 1].get_text()
    doc.close()
    return [x.strip() for x in t.splitlines() if x.strip()]


def table_s1():
    """SI p.6. Eleven LNPs x four columns, emitted as pure text in row order."""
    ln = page_lines(6)
    i = ln.index("306O6")
    rows, j = [], i
    for name in ORDER:
        if ln[j] != name:
            raise SystemExit(f"Table S1: expected {name}, found {ln[j]!r}")
        vals = ln[j + 1:j + 5]
        rows.append({"LNP": name,
                     "z_avg_diameter_nm": float(vals[0]),
                     "PDI": float(vals[1]),
                     "mRNA_entrapment_pct": float(vals[2]),
                     "LNP_surface_pKa_TNS_aqueous": float(vals[3])})
        j += 5
    return rows


def table_s2():
    """SI p.15. Ten lipidoids; 306O18 is absent. 'x +/- y' per cell."""
    ln = page_lines(15)
    i = ln.index("306O6")
    rows, j = [], i
    for name in ORDER:
        if name == "306O18":
            continue
        if ln[j] != name:
            raise SystemExit(f"Table S2: expected {name}, found {ln[j]!r}")
        pka = re.match(r"([\d.]+)\s*±\s*([\d.]+)", ln[j + 1])
        buf = re.match(r"([\d.]+)\s*±\s*([\d.]+)", ln[j + 2])
        if not (pka and buf):
            raise SystemExit(f"Table S2 {name}: {ln[j+1]!r} {ln[j+2]!r}")
        rows.append({"Lipidoid": name,
                     "lipidoid_pKa_ethanol_titration": float(pka.group(1)),
                     "lipidoid_pKa_sd": float(pka.group(2)),
                     "buffering_capacity_uL_pH7.4_to_4": float(buf.group(1)),
                     "buffering_capacity_sd": float(buf.group(2))})
        j += 3
    return rows


def fig_s7_printed_pka():
    """SI p.13. Eight panels, each printing 'pKa= x' beside its lipidoid name.

    Paired GEOMETRICALLY, by nearest name to each pKa label, because the draw
    order is not consistent: in seven panels the name precedes the pKa line and
    in the 306O10 panel it follows it. A line-order rule silently left that one
    unpaired and reported 7 of 8 -- which is the shape of a check that scores
    less than it could and looks like a partial failure instead of a pairing
    bug (F-038: three marginal checks scored a literal 0 for twelve days).
    """
    doc = pymupdf.open(PDF)
    words = doc[13 - 1].get_text("words")
    doc.close()
    names, pkas = [], []
    for w in words:
        x, y, txt = (w[0] + w[2]) / 2.0, (w[1] + w[3]) / 2.0, w[4].strip()
        if txt in ORDER:
            names.append((x, y, txt))
        m = re.match(r"^([\d.]+)$", txt)
        if m and pkas and pkas[-1][2] is None:
            pkas[-1] = (pkas[-1][0], pkas[-1][1], float(m.group(1)))
        if txt.startswith("pKa="):
            tail = txt[4:].strip()
            pkas.append((x, y, float(tail) if tail else None))
    out, used = {}, set()
    for px, py, val in pkas:
        cands = sorted((((px - nx) ** 2 + (py - ny) ** 2) ** 0.5, nm)
                       for nx, ny, nm in names if nm not in used)
        if not cands:
            out[f"_unpaired_{len(out)}"] = val
            continue
        dist, nm = cands[0]
        used.add(nm)
        out[nm] = val
    return out


def main():
    s1, s2 = table_s1(), table_s2()
    s7 = fig_s7_printed_pka()

    with open(os.path.join(HERE, "04_table_s1.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(s1[0].keys()))
        w.writeheader(); w.writerows(s1)
    with open(os.path.join(HERE, "04_table_s2.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(s2[0].keys()))
        w.writeheader(); w.writerows(s2)

    # --- cross-check 1: Fig S7's printed pKa vs Table S1's pKa column -----
    cmp1 = []
    for k, v in s7.items():
        if k.startswith("_"):
            cmp1.append({"lipidoid": k, "fig_s7": v, "table_s1": None,
                         "match": False})
            continue
        t = next(r["LNP_surface_pKa_TNS_aqueous"] for r in s1 if r["LNP"] == k)
        cmp1.append({"lipidoid": k, "fig_s7": v, "table_s1": t,
                     "delta": round(v - t, 4), "match": v == t})
    # --- cross-check 2: the two pKa CONVENTIONS, deliberately not pooled ---
    cmp2 = []
    for r in s2:
        t = next(x["LNP_surface_pKa_TNS_aqueous"] for x in s1
                 if x["LNP"] == r["Lipidoid"])
        cmp2.append({"lipidoid": r["Lipidoid"],
                     "lipidoid_pKa_ethanol": r["lipidoid_pKa_ethanol_titration"],
                     "LNP_surface_pKa_TNS": t,
                     "delta": round(r["lipidoid_pKa_ethanol_titration"] - t, 4)})
    # --- cross-check 3: mole percents sum to exactly 100 -------------------
    molpct = {"lipidoid": 35.0, "second_listed": 16.0, "third_listed": 46.5,
              "C14-PEG2000": 2.5}
    qc = {
        "_what": "Tier-4 PDF-text tables and the checks they make possible.",
        "table_s1_rows": len(s1), "table_s2_rows": len(s2),
        "table_s2_missing_vs_s1": sorted(
            {r["LNP"] for r in s1} - {r["Lipidoid"] for r in s2}),
        "fig_s7_printed_pka_count": len(s7),
        "check_figS7_vs_tableS1_pKa": {
            "n": len([c for c in cmp1 if c["table_s1"] is not None]),
            "exact_matches": sum(1 for c in cmp1 if c["match"]),
            "max_abs_delta": max((abs(c.get("delta", 99)) for c in cmp1
                                  if c["table_s1"] is not None), default=None),
            "per_lipid": cmp1,
        },
        "check_two_pKa_conventions_DO_NOT_POOL": {
            "_what": ("Table S2 is the LIPIDOID pKa in ethanol; Table S1 is "
                      "the LNP SURFACE pKa in aqueous buffer by TNS. Different "
                      "material, different solvent, different method. This is "
                      "reported as a SPREAD, not as a disagreement to fix."),
            "n": len(cmp2),
            "max_abs_delta_pH": max(abs(c["delta"]) for c in cmp2),
            "mean_delta_pH": round(
                sum(c["delta"] for c in cmp2) / len(cmp2), 4),
            "n_with_delta_above_0.3_pH": sum(1 for c in cmp2
                                             if abs(c["delta"]) > 0.3),
            "per_lipid": cmp2,
        },
        "check_mole_percents_sum_to_100": {
            "_what": ("The SI prints 'a molar ratio of 35 : 16 : 46.5 : 2.5' "
                      "for four components. The sum is a round number the "
                      "authors never printed, and it confirms these are mol% "
                      "with no fifth component unaccounted for. WHICH "
                      "component takes 16 and which takes 46.5 is a separate "
                      "question and is DECISIONS.md D-1."),
            "values": molpct, "sum": round(sum(molpct.values()), 6),
            "equals_100_exactly": sum(molpct.values()) == 100.0,
        },
    }
    with open(os.path.join(HERE, "04_tables_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(json.dumps({k: v for k, v in qc.items() if not k.startswith("_")},
                     indent=1, default=str)[:3000])


if __name__ == "__main__":
    main()
