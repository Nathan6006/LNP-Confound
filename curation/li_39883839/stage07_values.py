#!/usr/bin/env python3
"""Assemble 05_values.csv: one row per MEASUREMENT, with its exact locus.

One row per (lipid, panel), not per lipid: three panels measured overlapping
sets of lipids in three different cohorts, and collapsing them would hide the
disagreement V4 found (2.8x and 25.8x on the two shared lipids).

Benchmarks resolve through curation/lib/reagents, i.e. through
provenance/reagent_registry.py, which writes the structure FROM THE CHEMICAL
NAME. Error class 30 is an unguarded corpus name -> structure lookup; `smi`
raises on an ambiguous name rather than taking the first row.
"""
import csv, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import reagents                                              # noqa: E402
HERE = os.path.dirname(os.path.abspath(__file__))

COMMON = dict(
 Publication_PMID="39883839", Experiment_proposed_ID="YX_2025",
 cargo="firefly luciferase mRNA (mLuc), TriLink, unmodified",
 readout_as_printed="Total flux (p/s)",
 readout_full="whole-animal IVIS bioluminescence total flux, photons/s; "
              "D-luciferin 10 mg/mL 0.2 mL i.p., imaged 10 min after luciferin",
 instrument="IVIS (PerkinElmer), Living Image 4.8.0",
 biological_system="mouse, strain NOT STATED for the I.M. mLuc panels "
                   "(Methods lists C57BL/6, BALB/c and Ai9 for the study)",
 organ_or_site="injection site / muscle (whole-animal image)",
 route="intramuscular (I.M.)",
 dose="0.05 mg/kg mLuc per mouse",
 timepoint="6 h post injection",
 formulation="ionizable:DSPC:cholesterol:C14-PEG2000 = 50:10:38.5:1.5 mol%; "
             "0.1 ug/uL total mRNA; citrate pH 4.0 aqueous phase; T-junction "
             "mixing; dialysed 12 h vs 1x PBS, 10 kDa MWCO",
 NP_ratio="NOT STATED anywhere in the article or SI",
 purity="purified lipid (Methods: 'Except for the high-throughput screening "
        "and batch-testing stages, all other experiments utilized purified "
        "lipids') -- Fig 2H is a batch-testing stage panel, so its four LNPs "
        "are crude Passerini products; Fig 2K and Fig 3I are purified")

PANELS = {
 "Fig2H": dict(locus="Figure 2H, 'Batch 3 (Classified by tail B)', bars B1-B4",
               n="2", cohort="batch 3",
               caption_dose="0.05 mg kg-1 mLuc per mouse, n = 2 (Fig 2F caption)"),
 "Fig2K": dict(locus="Figure 2K, six bars", n="3", cohort="final comparison",
               caption_dose="0.05 mg kg-1 mLuc per mouse, n = 3"),
 "Fig3I": dict(locus="Figure 3I, three bars", n="3", cohort="spacer series",
               caption_dose="0.05 mg kg-1 mLuc per mouse, n = 3"),
}

BENCH = {"DLin-MC3-DMA": "DLin-MC3-DMA", "SM-102": "SM-102"}


def main():
    enum = {r["IL_name"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
    f2 = json.load(open(os.path.join(HERE, "05_fig2.json")))
    f3 = json.load(open(os.path.join(HERE, "05_fig3i.json")))
    src = [("Fig2H", f2["Fig2H"]["bars"]), ("Fig2K", f2["Fig2K"]["bars"]),
           ("Fig3I", f3["bars"])]

    rows = []
    for pan, bars in src:
        P = PANELS[pan]
        for b in bars:
            nm = b["lipid"]
            if nm in BENCH:
                smi = reagents.smi(BENCH[nm])
                key = reagents.full_key(smi)
                member = "benchmark, not a library member"
            else:
                smi = enum[nm]["IL_SMILES"]
                key = enum[nm]["InChIKey"]
                member = ("library member" if enum[nm]["in_library"] == "1"
                          else "made outside the 9x4x4 grid")
            r = dict(COMMON)
            r.update(panel=pan, locus=P["locus"], n_mice=P["n"],
                     cohort=P["cohort"], caption_dose=P["caption_dose"],
                     IL_name=nm, IL_SMILES=smi, InChIKey=key,
                     library_status=member,
                     value=b.get("value_p_per_s", "") or "",
                     value_units="photons/s (total flux)",
                     value_is_bound="", upper_bound="", lower_bound="",
                     bar_top_row_px=b.get("top_row_px", ""),
                     flank_delta_px=b.get("flank_delta_px", ""),
                     tier="5_pixels",
                     refused=b.get("REFUSED", ""))
            if b.get("upper_bound_p_per_s"):
                r["value_is_bound"] = "upper"
                r["upper_bound"] = b["upper_bound_p_per_s"]
            if b.get("lower_bound_p_per_s"):
                r["value_is_bound"] = "lower"
                r["lower_bound"] = b["lower_bound_p_per_s"]
            rows.append(r)

    cols = ["panel", "locus", "IL_name", "InChIKey", "library_status", "value",
            "value_units", "value_is_bound", "upper_bound", "lower_bound", "refused",
            "n_mice", "cohort", "caption_dose", "bar_top_row_px",
            "flank_delta_px", "tier", "cargo", "readout_as_printed",
            "readout_full", "instrument", "biological_system",
            "organ_or_site", "route", "dose", "timepoint", "formulation",
            "NP_ratio", "purity", "Publication_PMID", "Experiment_proposed_ID",
            "IL_SMILES"]
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    got = [r for r in rows if r["value"] != ""]
    bnd = [r for r in rows if r["value_is_bound"] == "upper"]
    ref = [r for r in rows if r["refused"] and r["value"] == ""
           and r["value_is_bound"] == ""]
    low = [r for r in rows if r["value_is_bound"] == "lower"]
    print("  lower bound only   : %d  %s" % (len(low),
          sorted({r["IL_name"] for r in low})))
    print("05_values.csv: %d measurement rows over %d distinct InChIKeys"
          % (len(rows), len({r["InChIKey"] for r in rows})))
    print("  with a point value : %d  %s" % (len(got),
          sorted({r["IL_name"] for r in got})))
    print("  upper bound only   : %d  %s" % (len(bnd),
          sorted({r["IL_name"] for r in bnd})))
    print("  refused, no number : %d  %s" % (len(ref),
          sorted({r["IL_name"] for r in ref})))
    return 0


if __name__ == "__main__":
    sys.exit(main())
