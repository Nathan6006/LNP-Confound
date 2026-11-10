#!/usr/bin/env python3
"""06_rows.csv -- the corpus rows for KW_2019, in the 79-column merge schema.

WHAT IS EMITTED: the ELEVEN 306-family lipidoids of Figure S6A, at one dose,
one cell line, one readout. Nothing else in the paper resolves a delivery value
per lipidoid; 00_value_sources.json says so in the intake gate and
05_values.csv carries a per-row reason for every value NOT emitted.

FORMULATION, SI Methods p.2 verbatim: "mixing lipidoid, cholesterol (Sigma
Aldrich), 1,2-dioleoyl-sn-glycero-3-phosphoethanolamine (DOPE, Avanti Polar
Lipids), C14-PEG2000 (Avanti Polar Lipids), at a molar ratio of
35 : 16 : 46.5 : 2.5 ... the final weight ratio of lipidoid: mRNA was 10 : 1".
The four percentages sum to exactly 100 (04_tables_qc.json). WHICH component
takes 16 and which takes 46.5 is DECISIONS.md D-1: read positionally the
sentence gives cholesterol 16 / DOPE 46.5, which NO row in the 29,195-row
corpus does (20 experiments use this exact 35/x/y/2.5 ratio and every one of
them has helper lipid 16 / cholesterol 46.5; zero rows anywhere have
cholesterol 16 or a helper lipid 46.5). Emitted on the corpus convention and
flagged on every row, because a formulation field that differs from every other
experiment for a prose-ordering reason would manufacture exactly the kind of
source-separability signal this project exists to measure.

ENDPOINT, SI Methods p.3 verbatim: HeLa (ATCC), 15,000 cells/well in white
96-well plates, transfected 6-8 h after seeding with LNPs carrying UNMODIFIED
firefly-luciferase mRNA (TriLink) at 100 ng per well (0.81 nM); Bright-Glo
(Promega) read 24 h post-transfection; n = 4.

Run: python3 -I stage06_rows.py   (after stage01, stage02, stage03)
"""
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import assemble                                                 # noqa: E402
import emit                                                     # noqa: E402
import reagents                                                 # noqa: E402

EXPERIMENT_ID = "KW_2019"
DOI = "10.1002/smll.201805097"
PMID = "30637934"
LINK = "https://onlinelibrary.wiley.com/doi/10.1002/smll.201805097"
CORPUS = os.path.join(HERE, "..", "..", "data")

COLS = [
    "Index", "LNP_ID", "Experiment_ID", "Formulation_ID", "IL_name",
    "IL_SMILES", "IL_protonated_SMILES", "IL_head_name", "IL_head_SMILES",
    "IL_linker_name", "IL_linker_SMILES", "IL_tail1_name", "IL_tail1_SMILES",
    "IL_tail2_name", "IL_tail2_SMILES", "IL_tail3_name", "IL_tail3_SMILES",
    "IL_tail4_name", "IL_tail4_SMILES", "IL_molratio",
    "IL_to_nucleicacid_massratio", "IL_to_nucleicacid_chargeratio",
    "HL_name", "HL_SMILES", "HL_molratio", "CHL_name", "CHL_SMILES",
    "CHL_molratio", "PEG_name", "PEG_SMILES", "PEG_molratio",
    "fifthcomponent_name", "fifthcomponent_SMILES", "fifthcomponent_molratio",
    "Aqueous_buffer", "Dialysis_buffer", "Mixing_method", "Model",
    "Model_type", "Model_target", "Route_of_administration", "Cargo",
    "Cargo_type", "Dose_ug_nucleicacid", "Experiment_method",
    "Experiment_batching", "Experiment_value", "Publication_link",
    "Publication_DOI", "Publication_PMID",
    "Heavy.Atoms", "Rings", "Aromatic.Rings", "Rotatable.Bonds", "ExactMolWt",
    "Topological.Polar.Surface.Area", "Hydrogen.Bond.Donors",
    "Hydrogen.Bond.Acceptors", "LogP", "Molar.Refractivity",
    "Fraction.sp3.Carbons", "sp3.Carbons", "Nitrogen.Count",
    "Molecular.Weight", "has_ester", "has_carbonate", "has_disulfide",
] + emit.columns()


def corpus_rows():
    out = []
    for fn in ("data_vitro.csv", "data_vivo.csv"):
        with open(os.path.join(CORPUS, fn), newline="") as fh:
            out += list(csv.DictReader(fh))
    return out


# --- the per-item confidence and note rules ---------------------------------
# Both are CALLABLES over the item's own enumeration record. emit.emit_rows
# refuses a literal (error class 29: F-011, F-025, F-046 -- three times, and
# the first inverted the labelling so the most uncertain rows carried the
# strongest claim).

ISOMER_MIXTURE_TAILS = {"Oi10"}     # isodecyl acrylate is sold as a mixture

CRUDE_NOTE = ("In-vitro values were measured on CRUDE 1:4 reaction product "
              "(SI Methods p.2: 'In vitro experiments were conducted with "
              "crude lipidoids'); the stored structure is the nominal "
              "tetra-substituted species, so structure_assumed=1. The 306 "
              "core's four N-H sites are two symmetric pairs, so the "
              "tetra-product has no regiochemical choice and "
              "positional_ambiguity is 0 -- the uncertainty is substitution "
              "DEGREE, not position.")
RATIO_NOTE = ("Excipient mol% read as helper 16 / cholesterol 46.5 against "
              "the SI sentence's positional order (cholesterol 16 / DOPE "
              "46.5); see DECISIONS.md D-1.")


def confidence_of(item_id, rec):
    """From the item's own tail record, never from a literal."""
    if rec["tail_token"] in ISOMER_MIXTURE_TAILS:
        return "medium"
    return "high"


def note_of(item_id, rec):
    bits = []
    if rec["tail_token"] in ISOMER_MIXTURE_TAILS:
        bits.append(
            "Tail is ISODECYL acrylate, a commercial branched-isomer mixture; "
            "stored as 8-methylnonyl, matching the corpus's NC_2024 A306_Oi10 "
            "(full-InChIKey identical). The specific branch is not pinned by "
            "the paper.")
    bits.append(CRUDE_NOTE)
    bits.append(RATIO_NOTE)
    bits.append(f"Value: Fig S6A bar top, tier 4b PDF vector coordinates "
                f"({rec['value_RLU']:,} RLU, s.d. {rec['sd_RLU']:,}); "
                f"{rec['precision_note']}")
    return " ".join(bits)


def structure_assumed_of(item_id, rec):
    return "1"          # every row: the tetra-substitution is nominal (crude)


def positional_ambiguity_of(item_id, rec):
    return "0"


def value_is_bound_of(item_id, rec):
    # NOT a bound. Every one of the eleven bars is DRAWN with a measurable
    # height above the zero tick (smallest 0.233 pt = 4,751 RLU), so these are
    # low-precision measurements, not censored values. F-040 is the precedent:
    # 65 rows flagged value_is_bound=upper turned out to be measurements, and
    # over-flagging a measurement as a bound is as wrong as the reverse.
    return ""


def main():
    with open(os.path.join(HERE, "03_enumerated.csv"), newline="") as fh:
        enum = {r["tail_token"]: r for r in csv.DictReader(fh)}
    with open(os.path.join(HERE, "02_fig_s6_panelA.csv"), newline="") as fh:
        panelA = {r["category"]: r for r in csv.DictReader(fh)}

    with open(os.path.join(HERE, "01_blocks.csv"), newline="") as fh:
        alcohols = {b["block_name"]: b["alcohol_SMILES"]
                    for b in csv.DictReader(fh) if b["alcohol_SMILES"]}

    crows = corpus_rows()
    DOPE = reagents.smi("DOPE", corpus_rows=crows)
    CHOL = reagents.smi("Cholesterol", corpus_rows=crows)
    PEG = reagents.smi("C14-PEG2000", corpus_rows=crows)
    AMINE306 = "CN(CCCN)CCCN"

    # --- build the records the flag callables read ----------------------
    recs, skipped = {}, []
    px_per_5e5 = json.load(open(os.path.join(
        HERE, "02_fig_s6_geometry.json")))["panelA_px_per_5e5_RLU"]
    rlu_per_px = 5.0e5 / px_per_5e5
    frame_rlu = json.load(open(os.path.join(
        HERE, "02_fig_s6_geometry.json")))["panelA_frame_offset_worth_RLU"]
    for tok, e in enum.items():
        pa = panelA.get(tok)
        if pa is None or not pa["value_RLU"]:
            skipped.append({"IL_name": e["IL_name"],
                            "reason": "no bar drawn in Figure S6A"})
            continue
        v = float(pa["value_RLU"])
        rec = dict(e)
        rec["value_RLU"] = int(v)
        rec["sd_RLU"] = int(pa["sd_RLU"]) if pa["sd_RLU"] else 0
        rec["bar_height_px"] = round(v / rlu_per_px, 3)
        rec["precision_note"] = (
            f"the fill/stroke frame offset of 0.057 pt is worth {frame_rlu:,.0f} "
            f"RLU, i.e. {100*frame_rlu/v:.1f}% of this value")
        recs[tok] = rec

    order = ["O6", "O8", "O9", "Oi10", "O10", "O11", "O12", "O13", "O14",
             "O16", "O18"]
    items = [t for t in order if t in recs]
    if len(items) + len(skipped) != len(enum):
        raise SystemExit(f"counts do not sum: {len(items)} emitted + "
                         f"{len(skipped)} skipped != {len(enum)} enumerated")

    # --- per-experiment z-score of log10, population sd (ddof=0) --------
    logs = [math.log10(recs[t]["value_RLU"]) for t in items]
    mu = sum(logs) / len(logs)
    sd = math.sqrt(sum((x - mu) ** 2 for x in logs) / len(logs))
    z = {t: (math.log10(recs[t]["value_RLU"]) - mu) / sd for t in items}

    flags = emit.emit_rows(
        items, recs,
        confidence=confidence_of, note=note_of,
        structure_assumed=structure_assumed_of,
        positional_ambiguity=positional_ambiguity_of,
        value_is_bound=value_is_bound_of,
        value_normalization_unit=lambda i, r: "experiment",
        within_model=None, label="KW_2019")

    rows = []
    for tok, fl in zip(items, flags):
        e = recs[tok]
        # The tail acrylate (the reagent as bought) is built by ESTERIFYING
        # acrylic acid with the tail alcohol through curation/lib/assemble,
        # the same audited path stage01 uses for the Michael arm. The first
        # version of this line concatenated "C=CC(=O)" onto the alcohol
        # SMILES, whose O is terminal, and silently produced
        # `C=CC(=O)CCCCCCO` -- a ketone bearing a free alcohol, not an ester.
        # It parsed. That is the string-surgery hazard this file's own comment
        # warned about two lines later, committed anyway, and caught only by
        # printing the column.
        acrylate = assemble.esterify("C=CC(=O)O", alcohols[tok],
                                     label=f"{tok} acrylate reagent")
        r = {c: "" for c in COLS}
        r.update({
            "Experiment_ID": EXPERIMENT_ID,
            "Formulation_ID": f"{EXPERIMENT_ID}_F1",
            "IL_name": e["IL_name"],
            "IL_SMILES": e["IL_SMILES"],
            "IL_head_name": "A306",
            "IL_head_SMILES": AMINE306,
            "IL_tail1_name": tok,
            "IL_tail1_SMILES": acrylate,
            "IL_molratio": 35,
            "IL_to_nucleicacid_massratio": 10,
            "HL_name": "DOPE", "HL_SMILES": DOPE, "HL_molratio": 16,
            "CHL_name": "Cholesterol", "CHL_SMILES": CHOL,
            "CHL_molratio": 46.5,
            "PEG_name": "C14-PEG2000", "PEG_SMILES": PEG, "PEG_molratio": 2.5,
            "Aqueous_buffer": "citrate",
            "Dialysis_buffer": "",      # in vitro LNPs were NOT dialysed
            "Mixing_method": "handmixed",
            "Model": "in_vitro", "Model_type": "HeLa",
            "Model_target": "in_vitro", "Route_of_administration": "in_vitro",
            "Cargo": "mRNA", "Cargo_type": "FLuc",
            "Dose_ug_nucleicacid": 0.1,
            "Experiment_method": "luminescence_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": z[tok],
            "Publication_link": LINK, "Publication_DOI": DOI,
            "Publication_PMID": PMID,
        })
        r.update(fl)
        rows.append(r)

    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader(); w.writerows(rows)

    qc = {
        "experiment_id": EXPERIMENT_ID,
        "enumerated": len(enum), "emitted": len(rows),
        "skipped": skipped,
        "counts_sum": len(rows) + len(skipped) == len(enum),
        "zscore_mean": round(sum(z.values()) / len(z), 12),
        "zscore_sd_ddof0": round(math.sqrt(
            sum(v * v for v in z.values()) / len(z)), 12),
        "rows_with_a_caveat": emit.caveat_count(rows),
        "confidence_counts": {c: sum(1 for r in rows if r["confidence"] == c)
                              for c in ("high", "medium", "low")},
        "excipients_resolved_by_lib_reagents": {
            "DOPE": DOPE[:50] + "...", "Cholesterol": CHOL[:50] + "...",
            "C14-PEG2000": PEG[:50] + "..."},
        "dialysis_buffer_blank_because":
            ("SI Methods p.2 dialyses only the in-vivo LNPs ('For in vivo "
             "experiments, the LNPs were dialyzed against PBS'); the in-vitro "
             "LNPs were diluted into PBS and used. Blank is 'not stated to "
             "have happened', not 'clean'."),
    }
    with open(os.path.join(HERE, "06_rows_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=1)
    print(json.dumps(qc, indent=1))
    print()
    for r in rows:
        print(f"  {r['IL_name']:12s} z={r['Experiment_value']:+.4f} "
              f"conf={r['confidence']:7s} assumed={r['structure_assumed']} "
              f"tail1={r['IL_tail1_SMILES']}")


if __name__ == "__main__":
    main()
