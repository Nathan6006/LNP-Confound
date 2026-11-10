#!/usr/bin/env python3
"""Which b-LNPs were EXCLUDED from the pooled administration?

Table S1/S2's caption: "Gray text indicates the LNP was not included in the
pooled administration for high throughput in vivo screening due to not meeting
the selection criteria (< 200 nm in size, < 0.3 PDI, > 50% encapsulation
efficiency, and > 15 ng/uL mRNA)."

Greying is a COLOUR property, invisible to a plain text extraction -- so the
pool size looked "not determinable" until the character colours were read.
pdfplumber exposes `non_stroking_color`, and the two values on these pages are
(0,0,0) and (0.7490196,)*3, i.e. 191/255 grey. So pool membership IS recoverable.

AND IT IS CHECKABLE AGAINST THE STATED CRITERIA, which is the point: the caption
names four numeric thresholds, Table S2 prints the four numbers, so "is this row
grey?" and "does this row fail a threshold?" are two independent statements that
must agree. A row grey for no stated reason, or failing a threshold while black,
is a finding.
"""
import json, os, re
import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "nn5c19138_si_001.pdf")
GREY = (0.7490196, 0.7490196, 0.7490196)


def row_colours(pages):
    """b-LNP number -> True if its row's digits are grey."""
    out = {}
    with pdfplumber.open(SI) as pdf:
        for pno in pages:
            page = pdf.pages[pno - 1]
            # group characters into visual rows by their top coordinate
            rows = {}
            for ch in page.chars:
                rows.setdefault(round(ch["top"] / 3.0), []).append(ch)
            for _, chs in rows.items():
                chs.sort(key=lambda c: c["x0"])
                txt = "".join(c["text"] for c in chs)
                m = re.match(r"^\s*(?:C1[024]-\d+\s*)?(\d{1,3})\b", txt)
                if not m:
                    continue
                n = int(m.group(1))
                if not (1 <= n <= 102):
                    continue
                digits = [c for c in chs if c["text"].isdigit()]
                if not digits:
                    continue
                grey = sum(1 for c in digits
                           if c.get("non_stroking_color") == GREY)
                out.setdefault(n, []).append(grey / len(digits))
    return {n: (max(v) > 0.5) for n, v in out.items()}


def main():
    tab = json.load(open(os.path.join(HERE, "05_tables.json")))
    s2 = tab["tableS2"]
    grey1 = row_colours([22, 23, 24])
    grey2 = row_colours([25, 26, 27])

    def fails(n):
        """The caption's four stated criteria, applied to Table S2's numbers."""
        r = s2.get(str(n))
        if r is None:
            return None, "absent"
        bad = []
        def num(x):
            x = (x or "").replace("%", "").strip()
            try:
                return float(x)
            except ValueError:
                return None
        size, pdi = num(r["Size_nm"]), num(r["PDI"])
        ee, conc = num(r["EE_pct"]), num(r["RNA_ng_uL"])
        if size is None: bad.append("size not measurable")
        elif size >= 200: bad.append(f"size {size} >= 200 nm")
        if pdi is None: bad.append("PDI not measurable")
        elif pdi >= 0.3: bad.append(f"PDI {pdi} >= 0.3")
        if ee is None: bad.append("EE not measurable")
        elif ee <= 50: bad.append(f"EE {ee}% <= 50%")
        if conc is None: bad.append("conc not measurable")
        elif conc <= 15: bad.append(f"RNA {conc} <= 15 ng/uL")
        return (len(bad) > 0), "; ".join(bad)

    rec, disagree = {}, []
    for n in sorted(set(grey1) | set(grey2) | {int(k) for k in s2}):
        g = grey1.get(n, grey2.get(n))
        f, why = fails(n)
        rec[n] = dict(grey_in_table=g, fails_stated_criteria=f, reason=why)
        if g is not None and f is not None and g != f:
            disagree.append((n, g, f, why))

    tabulated = sorted(int(k) for k in s2)
    greyed = sorted(n for n, v in rec.items() if v["grey_in_table"])
    out = dict(
        n_bLNP_numbers_assigned=98,
        n_bLNP_rows_tabulated=len(tabulated),
        numbers_assigned_but_never_tabulated=sorted(
            set(range(1, 99)) - set(tabulated)),
        n_greyed_excluded_from_pool=len(greyed),
        greyed_bLNPs=greyed,
        n_in_pooled_administration=len([n for n in tabulated
                                        if n <= 98 and not rec[n]["grey_in_table"]]),
        P14_grey_vs_stated_criteria=dict(
            n_compared=sum(1 for v in rec.values()
                           if v["grey_in_table"] is not None
                           and v["fails_stated_criteria"] is not None),
            n_disagreeing=len(disagree), disagreements=disagree,
            what=("The caption's four thresholds and the grey colouring are two "
                  "independent statements of the same fact. Agreement validates "
                  "both the colour read and the criteria transcription.")),
        per_bLNP=rec)
    json.dump(out, open(os.path.join(HERE, "09_pool_membership.json"), "w"), indent=1)
    print("tabulated b-LNP rows :", out["n_bLNP_rows_tabulated"])
    print("numbers never tabulated:", out["numbers_assigned_but_never_tabulated"])
    print("greyed (excluded)    :", len(greyed), greyed)
    print("IN THE POOL          :", out["n_in_pooled_administration"])
    print("P14 grey vs criteria : compared",
          out["P14_grey_vs_stated_criteria"]["n_compared"],
          "disagreeing", out["P14_grey_vs_stated_criteria"]["n_disagreeing"])
    for d in disagree:
        print("   disagree:", d)


if __name__ == "__main__":
    main()
