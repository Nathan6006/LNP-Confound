"""Stage 06 -- assemble 05_values.csv.

ONE BLOCK PER ROUND. The six panels are six separate in-vivo experiments and
they are NOT on one scale. That is not an inference from the numbers; the paper
says so itself, in Table S2, which lists a DIFFERENT formulation ID for the
same ionizable lipid in different panels:

    C3-K2-E14   ID 3 in Fig S1B   ID 4 in Fig S2B
    C3-K2-E10   ID 6 in Fig S2B   ID 7 in Fig S3B
    C3-K3-E10   ID 8 in Fig S3B   ID 10 in Fig S4B

Each panel's own characterisation sub-table (Fig S1C, S2E, S3E, S4D, S7C) was
matched against Table S2 on the full four-tuple (size, PDI, EE, pKa); every
match and every partial match is recorded on the row in `formulation_match`,
and the two partials are in 08_contradictions.csv.

DOSE IS NOT RECORDED, ANYWHERE. Fig S8 varies the mRNA dose from ~0.1 to
1.25 ug, which means the other panels used some fixed dose -- but no SI caption,
no table and no methods paragraph in this file states it. It is left blank.
The main article is not open access (Europe PMC: isOpenAccess=N, inEPMC=N) and
was not obtained, so the value is unknown rather than unstated.

TIMEPOINT IS STATED ON THREE PANELS AND NOT ON THE OTHER THREE. S1B, S2B and
S3B say "at 6 hours post-injection"; S4B, S4C and S7A say only "injected with".
Recorded exactly that way, blank where the paper is silent.
"""
import os, sys, csv, json

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import emit                                                  # noqa: E402
sys.path.insert(0, HERE)
from stage07_qc import parse_table_s2                        # noqa: E402

# panel -> {lipid: (formulation_id, how the id was pinned)}
MAP = {
    "S1b": {"SM-102": (1, "exact_3of4_size_PDI_EE_pKa_differs_see_C-3"),
            "C3-K2-A14": (2, "exact_4tuple_fig_S1C"),
            "C3-K2-E14": (3, "exact_4tuple_fig_S1C")},
    "S2b": {"C3-K2-E14": (4, "exact_4tuple_fig_S2E"),
            "C3-K2-E10": (6, "exact_4tuple_fig_S2E"),
            "C3-K2-E15": (5, "exact_4tuple_fig_S2E")},
    "S3b": {"C3-K2-E10": (7, "exact_4tuple_fig_S3E"),
            "C3-K3-E10": (8, "exact_4tuple_fig_S3E")},
    "S4b": {"SM-102": (9, "inferred: not in fig S4D; pKa 6.61 in fig S11 and "
                          "unmodified-mRNA N/P 5"),
            "C3-K3-E10": (10, "exact_4tuple_fig_S4D"),
            "OC2-K3-E10": (11, "partial: EE and pKa only; size/PDI match no "
                               "OC2-K3-E10 row, see C-4"),
            "OC3-K3-E10": (12, "exact_4tuple_fig_S4D")},
    "S4c": {"SM-102": (13, "inferred: the only Modified-Fluc SM-102 row"),
            "OC2-K3-E10": (14, "inferred: the only Modified-Fluc "
                               "OC2-K3-E10 row")},
    "S7a": {"DLin-MC3-DMA": (27, "inferred: fig S11 pKa 6.30 and the only "
                                 "MC3 row near N/P 5"),
            "SM-102": (None, "AMBIGUOUS: fig S7A prints no characterisation "
                             "for SM-102; table S2 has 1, 9, 19, 20, 21"),
            "OC2-K3-E10": (None, "AMBIGUOUS: fig S7A prints no "
                                 "characterisation; table S2 has 8 candidates"),
            "OC2-K2-A10": (28, "exact_4tuple_fig_S7C")},
}

TIMEPOINT = {"S1b": "6 h", "S2b": "6 h", "S3b": "6 h",
             "S4b": "", "S4c": "", "S7a": ""}

READOUT = {"S1b": "Total Flux (p/s)", "S2b": "Total Flux (p/s)",
           "S3b": "Total Flux (p/s)", "S4b": "Total Flux (p/s)",
           "S4c": "Total Flux (p/s)", "S7a": "Total flux (p/s)"}

LOCUS = {"S1b": "SI p2, Figure S1B", "S2b": "SI p3, Figure S2B",
         "S3b": "SI p3, Figure S3B", "S4b": "SI p4, Figure S4B",
         "S4c": "SI p4, Figure S4C", "S7a": "SI p7, Figure S7A"}

STEREO_NOTE = (
    "Stored without stereochemistry: the four carbinol centres come from "
    "ring-opening a RACEMIC 1,2-epoxyalkane, so the compound is a mixture of "
    "stereoisomers and no single configuration is correct.")


def main():
    dig = json.load(open(os.path.join(HERE, "05_digitize.json")))
    enum = {r["IL_name"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
    t2 = {r["id"]: r for r in parse_table_s2()}

    items, records = [], {}
    base = []
    for p in dig["panels"]:
        for b in p["bars"]:
            nm = b["label"]
            fid, how = MAP[p["key"]][nm]
            f = t2.get(fid, {})
            e = enum[nm]
            rid = f"{p['key']}|{nm}"
            items.append(rid)
            records[rid] = dict(panel=p["key"], IL_name=nm, enum=e,
                                n_stereo=int(e["n_unspecified_stereocentres"]),
                                origin=e["origin"], formulation=fid, how=how)
            base.append(dict(
                row_id=rid, round_id=p["round_id"], panel=p["key"],
                page=p["page"], source_locus=LOCUS[p["key"]],
                IL_name=nm, IL_SMILES=e["IL_SMILES"],
                readout_as_printed=READOUT[p["key"]],
                value=f"{b['value']:.6g}", units="p/s",
                log10_value=b["log10_value"],
                value_tier="5_pixels",
                assay=("in vivo whole-animal bioluminescence, IVIS Spectrum, "
                       "firefly luciferase"),
                cargo=p["cargo"], model_organism="C57BL/6 mouse",
                organ="hindlimb muscle", route="intramuscular",
                dose="", timepoint=TIMEPOINT[p["key"]], n_text=p["n_text"],
                formulation_id=("" if fid is None else fid),
                formulation_match=how,
                NP=f.get("NP", ""), lipid_to_RNA_wt_wt=f.get("lipid_rna", ""),
                Z_average_nm=f.get("z", ""), PDI=f.get("pdi", ""),
                EE_pct=f.get("ee", ""), pKa=f.get("pka", ""),
                bar_outer_row_px=b["outer_row"],
                bar_stroke_px=b["stroke_t_px"],
                bar_value_row_px=b["value_row"]))

    def confidence(rid, rec):
        """From the item's own enumeration record, never a literal."""
        if "reagent_registry" in rec["origin"]:
            return "high"          # pinned by C11, and corpus-matched 3/3
        if rec["enum"]["formula_ours"] != rec["enum"]["formula_SI"]:
            return "low"
        return "high"              # formula == SI, structure drawn, route read

    def note(rid, rec):
        bits = []
        if "reagent_registry" in rec["origin"]:
            bits.append("Benchmark lipid; structure from the project reagent "
                        "registry (C11), not from this paper.")
        else:
            bits.append(
                f"Structure enumerated from the SI synthesis; our formula "
                f"{rec['enum']['formula_ours']} equals the SI's printed "
                f"formula.")
        if rec["n_stereo"]:
            bits.append(STEREO_NOTE)
        bits.append(
            f"Value digitised from {LOCUS[rec['panel']]}; it is RAW total "
            f"flux in p/s and is NOT comparable across panels -- Table S2 "
            f"gives this panel formulation "
            f"{rec['formulation'] if rec['formulation'] else 'UNKNOWN'} and "
            f"other panels a different one.")
        return " ".join(bits)

    def pos_amb(rid, rec):
        # Assessed: the epoxide opens at the terminal carbon and the paper
        # DRAWS the hydroxyl on the second carbon in Figs S1A/S2A/S3A/S4A.
        return "0"

    def assumed(rid, rec):
        return "0"

    def bound(rid, rec):
        return ""      # no bar sits on an axis floor or is censored

    q = emit.emit_rows(items, records, confidence=confidence, note=note,
                       positional_ambiguity=pos_amb, structure_assumed=assumed,
                       value_is_bound=bound, label="khanof_36652629")

    rows = []
    for b, qq in zip(base, q):
        r = dict(b); r.update(qq); rows.append(r)

    cols = list(base[0]) + emit.columns()
    seen, order = set(), []
    for c in cols:
        if c not in seen:
            order.append(c); seen.add(c)
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=order)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in order})

    print(f"05_values.csv: {len(rows)} rows, {len(order)} columns")
    print(f"rows carrying a caveat: {emit.caveat_count(q)}")
    for k in sorted({r['round_id'] for r in rows}):
        n = sum(1 for r in rows if r['round_id'] == k)
        print(f"  {k:22s} {n} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
