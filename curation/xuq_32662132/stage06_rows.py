"""Stage 06 - corpus-schema rows, flags written ON THE ROW via emit.emit_rows.

WHAT IS EMITTED, AND WHY SO LITTLE. The only delivery values in reach belong to
one ionizable structure, 93-O17S. Of those, the set that the corpus schema can
represent without inventing a screen is Figure S2a: twelve distinct excipient
ratios plus the lipidoid-only formulation, each individually measured. That is a
formulation DoE on one lipid, which has direct precedent (ZT_2026, 486
formulations of one lipid). Everything else gets a per-row skip reason below and
the counts sum to the 44 rows of 05_values.csv.

Confidence is DERIVED per item from the enumeration record. A hard-coded literal
raises in emit.emit_rows (F-011).
"""
import csv, json, math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "lib")))
import emit                                              # noqa: E402

EXPERIMENT_ID = "XZ_2020"       # Xuewei Zhao, 2020. Free: not among the 54 live ids.

SKIP_REASONS = {
    "S3-93-O17S": "dose-response of ONE lipid (5/15/30/60 ug/mL). Four rows would "
                  "carry the same (Experiment_ID, IL_name) key with no chemical "
                  "difference between them. OPEN decision D-2.",
    "S3-9322-O17S": "dose-response of one lipid, same as D-2. NOTE: the reason changed on "
                    "2026-11-10. It used to read 'no structure' - head 9322 is now "
                    "recovered (corpus-attested as 1-(3-aminopropyl)-2-methylimidazole and "
                    "confirmed against Figure S1 at 0.20 px), so these four values DO have "
                    "a molecule and are blocked only by the dose-response shape. See "
                    "08_retracted_suspicions.csv R-06.",
    "S3-Lipofectamine 2000": "commercial transfection reagent control, not a library member.",
    "S3-untreated": "untreated control, not a library member.",
    "S2a-UT": "untreated control, not a library member.",
    "S2b": "cell viability (CCK-8). Not a delivery endpoint; the corpus excludes "
           "viability and hemolysis readouts.",
    "S2c": "time course of ONE lipid at ONE formulation. Same key collision as D-2, "
           "and the corpus has no time axis within an experiment.",
    "S2d": "dose response of ONE lipid at ONE formulation. Same as S2c.",
    "ABSTRACT": "the abstract does not say WHICH lipidoid achieved 8.2 % recombination, "
                "so the value has no molecule to attach to.",
}


def skip_reason(vid):
    for k in sorted(SKIP_REASONS, key=len, reverse=True):
        if vid.startswith(k):
            return SKIP_REASONS[k]
    return None


def _confidence(item_id, rec):
    """Derived from the record, never written as a literal.

    `high` would require this paper to have drawn or characterised the molecule.
    It drew nothing we can see and its one mass spectrum is unreadable, so the
    head identity rests entirely on records from OTHER papers. That is a real
    caveat and it is what keeps every row at `medium`.
    """
    head_drawn_here = rec.get("head_drawn_in_this_paper") == "1"
    tail_synth_here = rec.get("tail_synthesis_in_this_SI") == "1"
    mass_anchor_here = rec.get("mass_anchor_in_this_paper") == "1"
    if head_drawn_here and tail_synth_here and mass_anchor_here:
        return "high"
    if tail_synth_here:
        return "medium"
    return "low"


def _note(item_id, rec):
    return ("Head 93 is not drawn anywhere in this paper's Supporting Information and "
            "its one mass spectrum (Figure S1) is illegible at the published raster "
            "resolution; the head is 1-(3-aminopropyl)imidazole on the agreement of "
            "three external records (AA_2008 and KW_2014 corpus rows, and "
            "curation/xuq_32520527's read of a drawn figure). The tail is built from "
            "the named reagents of its own synthesis in SI section 1.2. The excipient "
            "ratio is the WEIGHT ratio %s (lipidoid:cholesterol:DOPE:DSPE-PEG); it is "
            "NOT converted to mole percent because DSPE-PEG's block length is a "
            "distribution, so HL/CHL/PEG_molratio are left blank."
            % rec["formulation_w_w"])


def main():
    vals = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))
    enum = {r["IL_name"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}

    emitted, skipped = [], []
    for v in vals:
        r = skip_reason(v["value_id"])
        if r is None and v["lipidoid"] == "93-O17S" and v["source_figure"] == "Figure S2a":
            emitted.append(v)
        else:
            skipped.append(dict(value_id=v["value_id"], lipidoid=v["lipidoid"],
                                reason=r or "not a Figure S2a 93-O17S formulation row"))
    if len(emitted) + len(skipped) != len(vals):
        raise SystemExit("counts do not sum")

    # --- Experiment_value: per-experiment z-score of log10, population sd ---
    y = [math.log10(float(v["value"])) for v in emitted]
    mu = sum(y) / len(y)
    sd = (sum((t - mu) ** 2 for t in y) / len(y)) ** 0.5
    z = [(t - mu) / sd for t in y]

    base = enum["93-O17S"]
    records, items = {}, []
    for v in emitted:
        iid = v["value_id"]
        items.append(iid)
        records[iid] = dict(base, formulation_w_w=v["formulation_w_w"],
                            head_drawn_in_this_paper="0",
                            tail_synthesis_in_this_SI="1",
                            mass_anchor_in_this_paper="0")

    quality = emit.emit_rows(
        items, records,
        confidence=_confidence, note=_note,
        positional_ambiguity=lambda i, r: 0,
        structure_assumed=lambda i, r: 0,
        value_is_bound=lambda i, r: "",
        value_normalization_unit=lambda i, r: "experiment",
        within_model=lambda i, r: "",         # filled below: one Model_type here
        label="XZ_2020 Figure S2a")

    COLS = open(os.path.join(HERE, "..", "..", "data", "data_vitro.csv")).readline().strip().split(",")
    rows = []
    for k, (v, zz, q) in enumerate(zip(emitted, z, quality)):
        row = {c: "" for c in COLS}
        row.update({
            "Index": k + 1, "LNP_ID": "%s_L%03d" % (EXPERIMENT_ID, k + 1),
            "Experiment_ID": EXPERIMENT_ID,
            "Formulation_ID": "%s_F%02d" % (EXPERIMENT_ID, k + 1),
            "IL_name": "93-O17S", "IL_SMILES": base["IL_SMILES"],
            "IL_head_name": "93", "IL_head_SMILES": "NCCCn1ccnc1",
            "IL_linker_name": "ester", "IL_tail1_name": "O17S", "IL_tail2_name": "O17S",
            "HL_name": "DOPE", "CHL_name": "Cholesterol", "PEG_name": "DSPE-PEG",
            "Aqueous_buffer": "25 mM sodium acetate pH 5.2",
            "Dialysis_buffer": "dialysed >=2 h, 3.5K MWCO (buffer not stated)",
            "Mixing_method": "handmixed (ethanol dropwise into 3 volumes of acetate buffer, vortex)",
            "Model": "in_vitro", "Model_type": "primary_human_CD8_T_lymphocyte",
            "Model_target": "in_vitro", "Route_of_administration": "in_vitro",
            "Cargo": "mRNA", "Cargo_type": "FLuc",
            "Experiment_method": "luminescence_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": round(zz, 6),
            "Experiment_value_within_model": round(zz, 6),
            "Publication_DOI": "10.1002/anie.202008082",
            "Publication_PMID": 32662132,
            "Publication_link": "https://doi.org/10.1002/anie.202008082",
        })
        for c, val in q.items():
            if c == "Experiment_value_within_model" and val == "":
                continue
            row[c] = val
        rows.append(row)

    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader(); w.writerows(rows)
    with open(os.path.join(HERE, "06_skipped.json"), "w") as f:
        json.dump({"n_values_total": len(vals), "n_emitted": len(emitted),
                   "n_skipped": len(skipped), "sums": len(emitted) + len(skipped) == len(vals),
                   "skipped": skipped}, f, indent=1)
    print("emitted", len(rows), "skipped", len(skipped), "of", len(vals),
          "| caveat rows", emit.caveat_count(rows),
          "| z mean %.2e sd %.6f" % (sum(z) / len(z), (sum((t - sum(z) / len(z)) ** 2
                                                           for t in z) / len(z)) ** 0.5))

if __name__ == "__main__":
    main()
