#!/usr/bin/env python3
"""Parse SI Tables S4 and S9 out of the PDF's text layer. Never transcribed.

Table S9 (pp 40-41) is the authoritative sequence list for peptoids 222-266 and
also carries the CALCULATED EXACT MASS for every one of them -- 45 independent
anchors against which an enumerated structure's composition can be checked.

Table S4 (p 14) covers the 34 lipid-block peptoids and carries the four DOE
parameters, including Total Lipid Carbons. It overlaps Table S9 on all 34
sequences, so the two tables check each other.

WHY A PARSER AND NOT A TRANSCRIPTION. A hand-typed sequence is a silent
single-character failure mode, and a peptoid IS its monomer order -- a formula
cannot see a transposition. Parsing makes the sequence a measured fact a
reviewer can re-derive with one command.

WHY WORD POSITIONS AND NOT extract_text(). Both tables wrap long sequence cells
across two or three physical lines, and extract_text() interleaves the wrapped
fragment with the NEXT row's numbers:

    Apd-Ehx-Hex-Ehx-Hex-Ehx-Hex-Ehx-
    246   25.2   6.455   93.19   1389.13   1390.14
    Hex

Read as lines, peptoid 246's sequence loses its final monomer and gains nothing
to say so. Clustering words by their y centre and then rejoining any fragment
that ends in a hyphen reconstructs the cell. The parser ASSERTS the rejoined
monomer count against the row's own printed monomer count where one exists, so
a mis-rejoin fails loudly instead of emitting a short peptoid.
"""
import os
import re

import pdfplumber

HERE = os.path.dirname(os.path.abspath(__file__))
SI = os.path.join(HERE, "source", "nn4c05513_si_001.pdf")

MONO = r"[A-Z][a-z]{2}"
SEQ_RE = re.compile(rf"^{MONO}(?:-{MONO})*-?$")
NUM_RE = re.compile(r"^\d+(?:\.\d+)?$|^n\.d$|^n\.d\.$|^Multimodal$")


class SeqParseError(Exception):
    """The table did not parse the way this file says it does."""


# The SI is stamped with a per-page download watermark running DOWN the right
# margin ("Downloaded / from / pubs.acs.org/... / by / PHILLIPS / EXETER /
# ACADEMY / user / on / 10 / November / 2026"), one word per line at x0 ~ 588.
# Those words land in the same y bands as table rows and broke the row
# classifier: a row carrying a wrapped sequence PLUS a watermark word is not
# "all sequence tokens", so it was silently dropped and peptoid 246 came out
# with the single monomer "Hex". Anything right of the table body is excluded.
X_TABLE_MAX = 560.0


def _rows(page, y_tol=2.5, x_max=None):
    """Words grouped into visual rows by y centre, watermark column excluded.

    x_max defaults to the MODULE GLOBAL, read at call time rather than bound as
    a default argument. That is not a style choice: with `x_max=X_TABLE_MAX` in
    the signature the value is frozen when the module is imported, so
    stage07_qc.py's seeded defect -- which sets seqtables.X_TABLE_MAX to 10000
    to re-admit the watermark -- changed nothing at all and the seed reported
    34/34 -> 34/34 while claiming to have broken the parser. A seed that
    mutates a shadowed name proves the check is fine when it has not tested it.
    """
    if x_max is None:
        x_max = X_TABLE_MAX
    words = [w for w in page.extract_words(use_text_flow=False)
             if w["x0"] < x_max]
    rows = []
    for w in sorted(words, key=lambda w: ((w["top"] + w["bottom"]) / 2, w["x0"])):
        yc = (w["top"] + w["bottom"]) / 2
        for r in rows:
            if abs(r["y"] - yc) <= y_tol:
                r["w"].append(w)
                break
        else:
            rows.append({"y": yc, "w": [w]})
    for r in rows:
        r["w"].sort(key=lambda w: w["x0"])
    return rows


def _parse_pages(pages, lo, hi):
    """Rows keyed by peptoid id, with sequence fragments rejoined."""
    recs = {}
    order = []
    pending = []          # sequence fragments seen before their id row
    current = None        # id whose sequence may still be continued
    with pdfplumber.open(SI) as pdf:
        for pno in pages:
            for r in _rows(pdf.pages[pno - 1]):
                toks = [w["text"] for w in r["w"]]
                if not toks:
                    continue
                # a pure sequence-fragment row
                if all(SEQ_RE.match(t) for t in toks):
                    frag = "".join(toks)
                    if current is not None and recs[current]["open"]:
                        recs[current]["seq"] += frag
                        recs[current]["open"] = frag.endswith("-")
                    else:
                        pending.append(frag)
                    continue
                # a record row: starts with an in-range integer
                if toks[0].isdigit() and lo <= int(toks[0]) <= hi:
                    pid = int(toks[0])
                    seqtoks = [t for t in toks[1:] if SEQ_RE.match(t)]
                    nums = [t for t in toks[1:] if NUM_RE.match(t)]
                    seq = "".join(seqtoks)
                    if pending:
                        seq = "".join(pending) + seq
                        pending = []
                    if pid in recs:
                        raise SeqParseError(f"peptoid {pid} seen twice")
                    recs[pid] = {"seq": seq, "nums": nums,
                                 "open": seq.endswith("-") or not seqtoks}
                    order.append(pid)
                    current = pid
                    continue
                current = None
    if pending:
        raise SeqParseError(f"unattached sequence fragments: {pending}")
    for pid, r in recs.items():
        r["seq"] = r["seq"].rstrip("-")
        r["monomers"] = r["seq"].split("-") if r["seq"] else []
        if r["open"] and not r["monomers"]:
            raise SeqParseError(f"peptoid {pid} parsed with no sequence")
        del r["open"]
    return recs, order


def table_s9():
    """{peptoid: {monomers, calc_exact_mass, found_exact_mass, ...}} for 222-266."""
    recs, order = _parse_pages([40, 41], 222, 266)
    out = {}
    for pid, r in recs.items():
        nums = [n for n in r["nums"] if NUM_RE.match(n) and "." in n or n.isdigit()]
        floats = [float(n) for n in r["nums"] if re.match(r"^\d+(\.\d+)?$", n)]
        # trailing two numbers are calculated then found exact mass
        if len(floats) < 2:
            raise SeqParseError(f"peptoid {pid}: fewer than 2 numeric cells")
        out[pid] = {"monomers": r["monomers"],
                    "calc_exact_mass": floats[-2],
                    "found_exact_mass": floats[-1],
                    "raw_nums": r["nums"]}
    return out, order


def table_s4():
    """{peptoid: {monomers, lipid_monomers, branched, lipid_carbons, bu_type}}."""
    recs, order = _parse_pages([14], 226, 266)
    out = {}
    for pid, r in recs.items():
        ints = [int(n) for n in r["nums"] if n.isdigit()]
        if len(ints) < 3:
            raise SeqParseError(f"peptoid {pid}: Table S4 row has {ints}")
        mons = r["monomers"]
        # THE GUARD THAT DID NOT EXIST. Table S4 wraps a long sequence across
        # three physical lines and the id row itself then carries NO sequence
        # token. When the watermark column defeated the row classifier, peptoid
        # 246 parsed as the single monomer ["Hex"] and nothing complained --
        # a one-monomer peptoid is not a peptoid. Every row of this table is a
        # headgroup plus its lipid block, so the parsed monomer count must be
        # exactly one more than the row's own printed lipid-monomer count.
        # That is the table checking its own rejoin, which is the only reason
        # to believe a rejoin at all.
        if len(mons) != ints[0] + 1:
            raise SeqParseError(
                f"peptoid {pid}: parsed {len(mons)} monomers {mons} but the "
                f"row prints {ints[0]} lipid monomers, so {ints[0] + 1} were "
                "expected (one headgroup + the lipid block). The wrapped "
                "sequence cell was not rejoined correctly.")
        out[pid] = {"monomers": mons,
                    "lipid_monomers": ints[0],
                    "branched_unsaturated": ints[1],
                    "lipid_carbons": ints[2]}
    return out, order
