#!/usr/bin/env python3
"""Join digitised values to structures and emit 05_values.csv.

WHAT IS EMITTED AS A CORPUS-CANDIDATE ROW, AND WHAT IS NOT.

  S22C  4 rows  EMITTED.  In vitro luciferase (RLU), IGROV1, a delivery
                          magnitude of the kind the corpus already holds.
  S24   4 rows  HELD.     In vivo average radiance, but the ORGAN is never
                          stated anywhere in the SI, and Model_target cannot
                          be guessed. DECISIONS.md D-2.
  S5   11 rows  HELD.     Spleen/liver RATIO -- a selectivity contrast, not a
  S23A  4 rows  HELD.     delivery magnitude. DECISIONS.md D-1. S5 and S23A
                          also disagree with each other (08_contradictions).

Every row carries its quality columns through curation/lib/emit.emit_rows,
which is the one writer of those columns and which refuses a hard-coded
confidence. Confidence here is DERIVED from each lipid's own tail-count
anchor, because a literal is what turned 24 disclosed-low rows into the
project's highest confidence (F-011).
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__),
                                                "..", "lib")))
import emit                                                     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# Assay metadata, read off the SI. Each field names the page it came from.
ASSAY = {
    "S22C": {
        "pool": "in_vitro", "panel_page": 34,
        "Model_type": "IGROV1", "Model_target": "luminescence_normalized",
        "cargo": "Firefly luciferase (Fluc) mRNA",
        "cell_or_organ": "IGROV1 human ovarian cancer cells",
        "route": "in vitro transfection, 10% FBS medium",
        "dose": "25 ng mRNA per well (1.5e4 cells/well, 96-well)",
        "timepoint": "24 h",
        "assay": "Bio-Lumi II Firefly Luciferase Assay Kit, RLU",
        "n": 3, "unit": "RLU",
        "source_meta_page": 7,
    },
    "S24_left": {
        "pool": "in_vivo", "panel_page": 36,
        "Model_type": "C57BL/6 mouse", "Model_target": "UNSTATED",
        "cargo": "Fluc mRNA", "cell_or_organ": "UNSTATED (see D-2)",
        "route": "intravenous", "dose": "0.25 mg/kg", "timepoint": "6 h",
        "assay": "IVIS average radiance [p/sec/cm2/sr]", "n": 3,
        "unit": "p/sec/cm2/sr", "source_meta_page": 11,
    },
    "S24_right": {
        "pool": "in_vivo", "panel_page": 36,
        "Model_type": "C57BL/6 mouse", "Model_target": "UNSTATED",
        "cargo": "Fluc mRNA", "cell_or_organ": "UNSTATED (see D-2)",
        "route": "intravenous", "dose": "0.25 mg/kg", "timepoint": "6 h",
        "assay": "IVIS average radiance [p/sec/cm2/sr]", "n": 3,
        "unit": "p/sec/cm2/sr", "source_meta_page": 11,
    },
    "S5": {
        "pool": "in_vivo", "panel_page": 17,
        "Model_type": "C57BL/6 mouse", "Model_target": "spleen_over_liver_ratio",
        "cargo": "Fluc mRNA", "cell_or_organ": "spleen and liver, ex vivo",
        "route": "intravenous", "dose": "0.25 mg/kg", "timepoint": "6 h",
        "assay": "ratio of IVIS average radiance efficiency, spleen / liver",
        "n": "not stated (no error bars drawn)",
        "unit": "ratio", "source_meta_page": 11,
    },
    "S23A": {
        "pool": "in_vivo", "panel_page": 35,
        "Model_type": "C57BL/6 mouse", "Model_target": "spleen_over_liver_ratio",
        "cargo": "Fluc mRNA", "cell_or_organ": "spleen and liver, ex vivo",
        "route": "intravenous", "dose": "not stated on this panel",
        "timepoint": "6 h",
        "assay": "ratio of IVIS average radiance efficiency, spleen / liver",
        "n": 3, "unit": "ratio", "source_meta_page": 11,
    },
}

EMITTED_PANELS = {"S22C"}


def confidence_of(item_id, rec):
    """DERIVED, never a literal. The structure's own tail-count anchor decides."""
    return {"hrms": "high",
            "inferred_linear": "medium",
            "inferred_branched": "medium"}[rec["tail_count_anchor"]]


def structure_assumed_of(item_id, rec):
    return "0" if rec["tail_count_anchor"] == "hrms" else "1"


def note_of(item_id, rec):
    bits = []
    if rec["tail_count_anchor"] == "hrms":
        bits.append(
            f"Structure confirmed by the SI's own HRMS on p.{rec['si_page']}: "
            f"our independently assembled SMILES gives [M+H]+ "
            f"{float(rec['ours_MH_plus']):.4f} against the printed calculated "
            f"{rec['si_calc_MH_plus']} "
            f"({float(rec['ppm_ours_vs_si_calc']):+.2f} ppm).")
    else:
        bits.append(
            f"ASSUMED: the number of carbonate tails. {rec['IL_name']} is "
            f"named only in Figure S5 and has no HRMS, NMR or drawn structure "
            f"anywhere in the SI; {rec['n_tails']} tails is inferred from the "
            f"{'linear' if rec['tail_count_anchor'].endswith('linear') else 'branched'} "
            "members of its own series that do carry HRMS. The headgroup, the "
            "carbonate linkage and the alcohol identity are not assumed.")
    if rec.get("tie_note"):
        bits.append(rec["tie_note"])
    return " ".join(bits)


def main():
    recs = {r["IL_name"]: r for r in csv.DictReader(
        open(os.path.join(HERE, "03_enumerated.csv")))}
    dig = json.load(open(os.path.join(HERE, "01_digitised.json")))

    # Flag the S22C pixel tie on the rows it affects, from the measurement.
    s22 = dig["qc"]["S22C"]["bars"]
    rowmap = {}
    for b in s22:
        rowmap.setdefault(b["row"], []).append(b["IL_name"])
    for r, names in rowmap.items():
        if len(names) > 1:
            for nm in names:
                recs[nm]["tie_note"] = (
                    "VALUE TIE: in Figure S22C the top rules of "
                    + " and ".join(names) + " are detected on the same pixel "
                    f"row ({r}), so their stored values are identical by "
                    "construction; at 69.0 px per decade the panel resolves "
                    "no better than 3.4% and cannot order these two.")

    out, held = [], []
    for row in dig["rows"]:
        panel = row["panel"]
        a = ASSAY[panel]
        rec = recs[row["IL_name"]]
        base = {
            "IL_name": row["IL_name"], "IL_SMILES": rec["IL_SMILES"],
            "InChIKey": rec["InChIKey"], "formula": rec["formula"],
            "linkage": rec["linkage"], "n_tails": rec["n_tails"],
            "panel": panel, "si_page": row["page"],
            "source_location": f"Figure {panel}, SI p.{row['page']}, "
                               f"{'bar' if panel != 'S5' else 'point'} "
                               f"labelled {row['IL_name']}",
            "raw_value": round(row["value"], 6), "unit": row["unit"],
            "pool": a["pool"], "Model_type": a["Model_type"],
            "Model_target": a["Model_target"], "cargo": a["cargo"],
            "cell_or_organ": a["cell_or_organ"], "route": a["route"],
            "dose": a["dose"], "timepoint": a["timepoint"],
            "assay": a["assay"], "replicates_n": a["n"],
            "methods_page": a["source_meta_page"],
            "emitted_as_corpus_candidate": panel in EMITTED_PANELS,
        }
        (out if panel in EMITTED_PANELS else held).append((base, rec))

    # Quality columns, through the one writer, for every row we keep.
    allrows = out + held
    ids = [f"{b['panel']}::{b['IL_name']}" for b, _ in allrows]
    records = {i: r for i, (_, r) in zip(ids, allrows)}
    quality = emit.emit_rows(
        ids, records,
        confidence=confidence_of,
        note=note_of,
        structure_assumed=structure_assumed_of,
        positional_ambiguity=lambda i, rec: "0",
        value_is_bound=lambda i, rec: "",
        value_normalization_unit=lambda i, rec: "",
        label="lius_40500862")

    final = []
    for (base, _), q in zip(allrows, quality):
        base.update(q)
        final.append(base)

    cols = list(final[0].keys())
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in final:
            w.writerow(r)

    caveats = emit.caveat_count(final)
    print(f"{len(final)} value rows written to 05_values.csv")
    print(f"  emitted as corpus candidates : {len(out)}  (panel S22C)")
    print(f"  held, not emitted            : {len(held)}")
    print(f"  rows_with_a_caveat           : {caveats}")
    for r in final:
        if r["emitted_as_corpus_candidate"]:
            print(f"   EMIT {r['IL_name']:<12} {r['raw_value']:>12.4g} "
                  f"{r['unit']}  confidence={r['confidence']}")
    return caveats


if __name__ == "__main__":
    main()
