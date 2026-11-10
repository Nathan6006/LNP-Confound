#!/usr/bin/env python3
"""Stage 04 -- put the quality flags ON THE ROW, through curation/lib/emit.

F-001: four of five curations computed their flags, wrote them into their own
QC JSON, and emitted rows without them. 1,424 corpus rows carry a convention
and say nothing about it. "A flag that lives only in a per-folder QC file is
not a flag, it is a memo."

This curation emits NO value rows (see 00_value_sources.json: the SI carries no
delivery endpoint at all), so there is nothing to merge today. The flags are
written onto 03_enumerated.csv anyway, because the facts that qualify these
structures were established today and the reader who eventually pairs them with
values will not re-derive them. The alternative -- "we will add the flags when
the values arrive" -- is exactly the shape F-001 took.

Every flag is derived from the lipid's OWN record by a callable; `emit_rows`
refuses a literal, which makes error class 29 unwritable rather than
discouraged.

    python3 -I stage04_flags.py        # rewrites 03_enumerated.csv in place
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import emit                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ENUM = os.path.join(HERE, "03_enumerated.csv")

# The SI contradicts itself on these, and the contradiction is recorded in
# 08_contradictions.csv with both loci. They are listed here as DATA, not as a
# judgement: each key says which record field is in conflict, and the callables
# below read that field rather than this comment.
P_SERIES_NOTE = (
    "The SI's printed formula and printed [M+H]+ agree with EACH OTHER and are "
    "one CH2 below this structure; the benzylic CH2 kept here is the one its "
    "own reagent name ('...benzyl 6-bromohexanoate'), its intermediate name "
    "('...oxy)methyl)phenoxy)propanoic acid', SI p.11), its 1H NMR 2H singlet "
    "at 5.04 and the Materials list ('4-hydroxybenzyl alcohol') all require. "
    "Open: DECISIONS.md D-3.")


def confidence(_id, rec):
    """high only when BOTH SI loci independently confirm the built structure."""
    if rec["core"] == "SP":
        return "low"
    if rec["linker"] == "P":
        return "low"
    if rec["formula_stated"] == "":
        return "low"
    if rec["formula_agrees"] == "1" and rec["mz_agrees"] == "1":
        return "high"
    if rec["si_self_consistent"] == "1":
        return "low"
    return "medium"


def note(_id, rec):
    if rec["core"] == "SP":
        return ("SP core: the SI names the tail alcohol twice and differently "
                "-- '2-heptylnonan-1-ol' (SI 3.1.47 body and the Materials "
                "list) vs '2-Heptyldecan-1-ol' (SI 3.2.2). The former is "
                "built here. No SP member has a formula or an m/z anywhere in "
                "the SI, so nothing in the paper discriminates them, and the "
                "two differ by one CH2 -- 2-heptyldecyl would make every SP "
                "lipid an ISOMER of its SS counterpart. Open: DECISIONS.md D-2."
                + (" " + P_SERIES_NOTE if rec["linker"] == "P" else ""))
    if rec["linker"] == "P":
        return P_SERIES_NOTE
    if rec["formula_stated"] == "":
        return "No formula or m/z for this lipid anywhere in the SI."
    if rec["formula_agrees"] == "1" and rec["mz_agrees"] == "1":
        return ""
    if rec["si_self_consistent"] == "1":
        return (f"The SI prints formula {rec['formula_stated']} and [M+H]+ "
                f"{rec['mz_calc_stated']}, which agree with each other and sit "
                f"{rec['mh_plus_delta_vs_SI']} Da below this structure "
                f"({rec['formula']}). The structure follows the SI's own "
                "reagent sentence and the rest of its series; the HRMS pair "
                "does not. Open: DECISIONS.md D-3.")
    return (f"The SI's own formula ({rec['formula_stated']}) and its own "
            f"printed [M+H]+ ({rec['mz_calc_stated']}) are mutually "
            f"inconsistent, so only one of them can be evidence. This "
            f"structure ({rec['formula']}, [M+H]+ {rec['mh_plus']}) is "
            "confirmed by whichever of the two is self-consistent with the "
            "series; see 07_qc.json check P3 and 08_contradictions.csv.")


def structure_assumed(_id, rec):
    return "1" if (rec["core"] == "SP" or rec["linker"] == "P"
                   or confidence(_id, rec) == "low") else "0"


def main():
    with open(ENUM) as fh:
        rows = list(csv.DictReader(fh))
    recs = {r["name"]: r for r in rows}
    flags = emit.emit_rows(
        items=[r["name"] for r in rows], records=recs,
        confidence=confidence, note=note,
        structure_assumed=structure_assumed,
        positional_ambiguity=lambda i, r: "0",
        # No value is emitted by this curation, so these two are deliberately
        # left as "not stated" rather than as a reassuring 0/"experiment".
        value_is_bound=None, value_normalization_unit=None, within_model=None,
        label="liu_41162324")
    qcols = emit.columns()
    out = []
    for r, f in zip(rows, flags):
        r = dict(r)
        r.update(f)
        out.append(r)
    cols = [c for c in rows[0] if c not in qcols] + qcols
    with open(ENUM, "w", newline="") as fh:
        w = csv.DictWriter(fh, cols)
        w.writeheader()
        w.writerows(out)
    n = emit.caveat_count(out)
    byconf = {}
    for r in out:
        byconf[r["confidence"]] = byconf.get(r["confidence"], 0) + 1
    print(f"flags written onto {len(out)} rows; rows_with_a_caveat = {n}")
    print("confidence:", json.dumps(byconf, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
