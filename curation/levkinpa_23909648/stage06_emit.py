#!/usr/bin/env python3
"""Stage 06 -- join values to structures and write the quality columns.

WHAT IS EMITTABLE. A row needs BOTH a determined single structure AND a value.
That intersection is Tables S5 and S6 restricted to the three purified
compounds: 27 rows from S5 (3 compounds x 3 DOPE ratios x 3 DNA doses) and 2
from S6 (the two unblended compounds) = 29.

WHAT IS NOT, AND WHY -- each a different reason, not one excuse:
  * 50 rows (Tables S3 and S4, the 12 crude lipidoids): the material is a
    one-pot MIXTURE, and for 9 of the 12 the amine head is also unknown.
  * 4 rows (Table S6's blends): a formulation of TWO ionizable lipids. There
    is no single IL_SMILES; emitting one of the two would assign the blend's
    value to a compound measured separately on the same panel.
  * 4 rows (the L2K controls): Lipofectamine 2000 is a proprietary mixture
    with no defined structure. curation/lib/reagents deliberately refuses to
    pin such a name, and PEG-lipids are absent from the registry for the same
    reason -- the composition is a distribution, so no single SMILES is right.

THE CONFIDENCE IS DERIVED, NOT WRITTEN. emit.emit_rows refuses a literal, and
the reason is F-011: a hard-coded "high" relabelled 24 rows whose own record
read "ASSUMED: the paper draws no structure" as the project's highest
confidence, inverting the labelling so the most uncertain rows carried the
strongest claim. So `_confidence()` below scores each structure against named,
checkable criteria and returns the note alongside, and the criteria are the
three independent witnesses this paper supplies: a printed IUPAC name, an HRMS
pair, and a 1H integral sum.

Run:  python3 -I stage06_emit.py
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import emit                                                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

PPM_TOL = 5.0          # routine HRMS acceptance
H_SUM_TOL = 3.0        # 1H integral sums against the formula H count

EMITTABLE_LIPIDS = {"1C12-D", "1C16-S", "1C16-D"}


def _confidence(item_id, rec):
    """high only if all three independent witnesses agree. Never a literal."""
    ok_name = bool(rec["iupac_name_as_printed"].strip())
    ppm = abs(float(rec["ppm_ours_vs_their_calc"]))
    ok_mass = ppm <= PPM_TOL
    # The FIGURE's integral sum is the primary datum; the text listing is a
    # transcription of it and 1C12-D's text is demonstrably wrong (it prints
    # 3H for a methyl the figure integrates at 6.38, on a structure whose
    # formula requires 6H). Scoring the figure, and recording the text
    # discrepancy as a contradiction, is the honest split.
    h_fig = float(rec["h_nmr_fig_integral_total"])
    h_formula = float(rec["formula_H_count"])
    ok_h = abs(h_fig - h_formula) <= H_SUM_TOL
    if ok_name and ok_mass and ok_h:
        return "high"
    return "medium"


def _note(item_id, rec):
    bits = []
    ppm = float(rec["ppm_ours_vs_their_calc"])
    h_text = int(rec["h_nmr_text_minus_formula"])
    bits.append(
        f"Structure from the full IUPAC name printed on SI p.{rec['si_page']}, "
        f"independently confirmed by HRMS ({ppm:+.3f} ppm, our computed "
        f"[M+H]+ vs the authors' calculated m/z)")
    if h_text != 0:
        bits.append(
            f"the SI's 1H NMR TEXT lists {h_text:+d} H against the formula "
            f"{rec['molecular_formula']}, which is a printing error in the "
            f"paper, not a structural doubt: the spectrum's own integral sum "
            f"is {rec['h_nmr_fig_integral_total']} H (see "
            f"08_contradictions.csv)")
    bits.append("the screened material for the unpurified library members of "
                "the same name is a one-pot mixture; THIS row is the "
                "column-purified single compound")
    return "; ".join(bits) + "."


def _bound(item_id, rec_and_value):
    return ""


def main():
    structs = {r["IL_name"]: r for r in csv.DictReader(
        open(os.path.join(HERE, "03_enumerated.csv")))}
    values = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))

    rows, skipped = [], []
    for v in values:
        lip = v["lipidoid"]
        if lip not in EMITTABLE_LIPIDS:
            if lip == "L2K":
                why = ("Lipofectamine 2000: a proprietary commercial mixture "
                       "with no defined structure. No corpus row can exist.")
            elif ":" in lip or "1C16-S:1C16-D" in v["value_id"]:
                why = "a blend of two ionizable lipids; no single IL_SMILES."
            else:
                why = ("crude one-pot mixture"
                       + ("" if lip.startswith("1C")
                          else " AND the amine head is unknown")
                       + "; see 03_not_determinable.csv and DECISIONS.md D-1.")
            skipped.append(dict(value_id=v["value_id"], lipidoid=lip,
                                table=v["source_table"], why=why))
            continue
        if v["source_table"] not in ("Table S5", "Table S6"):
            skipped.append(dict(
                value_id=v["value_id"], lipidoid=lip,
                table=v["source_table"],
                why="this panel's entry of that name is the CRUDE product; "
                    "only Tables S5/S6 report the purified compounds."))
            continue
        # Table S6's blends carry a ':' in the label and are caught above.
        rows.append(dict(value=v, struct=structs[lip]))

    # --- the quality columns, one call, per-item resolvers ----------------
    items = [r["value"]["value_id"] for r in rows]
    records = {r["value"]["value_id"]: r["struct"] for r in rows}

    # value_is_bound: a printed 0.0 with SD 0.0 is a non-detect, so the true
    # value is at or below the assay's detection limit -> an UPPER bound.
    # Derived from the value itself, per item, not asserted for the table.
    bound_map, note_extra = {}, {}
    for r in rows:
        vid = r["value"]["value_id"]
        val = r["value"]["value_as_printed"].strip()
        sd = r["value"]["sd_as_printed"].strip()
        if float(val) == 0.0 and float(sd or 0) == 0.0:
            bound_map[vid] = "upper"
            note_extra[vid] = (
                " This row's value is printed as 0 with an SD of 0, i.e. a "
                "non-detect: the true transfection efficiency is at or below "
                "the assay's detection limit, not exactly zero. It also "
                "cannot be log-transformed, which the corpus convention "
                "requires -- DECISIONS.md D-5.")
        else:
            bound_map[vid] = ""

    def note(item_id, rec):
        return _note(item_id, rec) + note_extra.get(item_id, "")

    qual = emit.emit_rows(
        items=items, records=records,
        confidence=_confidence,
        note=note,
        positional_ambiguity={i: "0" for i in items},
        structure_assumed={i: "0" for i in items},
        value_is_bound=bound_map,
        # value_normalization_unit and Experiment_value_within_model are left
        # UNASSESSED on purpose. No standardisation has been applied: the
        # values are raw printed percentages. Which unit the per-experiment
        # z-score should be taken over is a merge-time convention question
        # with a real answer at stake (these 29 rows span 3 DNA doses and 3
        # DOPE ratios), and one value is exactly 0 so log10 is undefined.
        # DECISIONS.md D-5. Writing a unit here would be asserting a
        # standardisation that was never performed.
        value_normalization_unit=None,
        within_model=None,
        label="LL_2013")

    out = []
    for r, q in zip(rows, qual):
        v, s = r["value"], r["struct"]
        rec = dict(
            Experiment_ID="LL_2013",
            IL_name=v["lipidoid"],
            IL_SMILES=s["IL_SMILES"],
            IL_molecular_formula=s["molecular_formula"],
            IL_InChIKey=s["inchikey_full"],
            Model="in_vitro",
            Model_target="NOT STATED IN THE SI",
            Model_type="NOT STATED IN THE SI",
            Cargo="plasmid DNA",
            Cargo_dose_ng=v["dna_ng"],
            Helper_lipid="DOPE",
            lipidoid_to_DOPE_molar_ratio=v["formulation_lipidoid_DOPE"],
            Experiment_value_raw_as_printed=v["value_as_printed"],
            Experiment_value_raw_sd=v["sd_as_printed"],
            Experiment_value_unit="percent of transfected cells",
            Experiment_value="",          # NOT standardised; see D-5
            Publication_PMID="23909648",
            Publication_DOI="10.1021/bc400158w",
            source_table=v["source_table"],
            source_page=v["source_page"],
            value_id=v["value_id"])
        rec.update(q)
        out.append(rec)

    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
        w.writeheader()
        for r in out:
            w.writerow(r)

    with open(os.path.join(HERE, "06_skipped.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["value_id", "lipidoid", "table",
                                           "why"])
        w.writeheader()
        for r in skipped:
            w.writerow(r)

    caveats = emit.caveat_count(qual)
    summary = dict(
        _what="Stage 06. Values joined to structures; quality columns written "
              "by curation/lib/emit, the single writer.",
        n_values_total=len(values),
        n_emitted=len(out),
        n_skipped=len(skipped),
        skipped_by_reason={},
        rows_with_a_caveat=caveats,
        confidence_distribution={},
        value_is_bound_distribution={},
        NOT_MERGED="This curation does NOT merge into data/. Two required "
                   "corpus fields are unresolved: Model_target/Model_type "
                   "(the cell line is not in the SI, D-2) and the "
                   "standardisation convention (D-5). build_dataset.py would "
                   "refuse the folder anyway while DECISIONS.md carries an "
                   "OPEN entry, which is the gate working.")
    for s in skipped:
        k = s["why"][:60]
        summary["skipped_by_reason"][k] = \
            summary["skipped_by_reason"].get(k, 0) + 1
    for q in qual:
        c = q["confidence"]
        summary["confidence_distribution"][c] = \
            summary["confidence_distribution"].get(c, 0) + 1
        b = q["value_is_bound"] or "(none)"
        summary["value_is_bound_distribution"][b] = \
            summary["value_is_bound_distribution"].get(b, 0) + 1
    with open(os.path.join(HERE, "06_emit_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)

    print(f"OK  {len(out)} rows emitted, {len(skipped)} skipped, of "
          f"{len(values)} values")
    print(f"    confidence: {summary['confidence_distribution']}")
    print(f"    value_is_bound: {summary['value_is_bound_distribution']}")
    print(f"    rows_with_a_caveat (F-024 definition): {caveats}")
    for k, n in summary["skipped_by_reason"].items():
        print(f"    skipped {n:>2}: {k}...")
    return 0


if __name__ == "__main__":
    sys.exit(main())
