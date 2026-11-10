#!/usr/bin/env python3
"""KG_2025 -- blocks, enumeration, values and the emitted corpus rows.

Run:  python3 -I curation/dahlmanje_41183194/stage05_emit.py <repo_root>

Writes 01_blocks.csv, 03_enumerated.csv, 05_values.csv, 06_rows.csv.

THE HELPER LIPID'S SMILES HAS ONE SOURCE, AND THIS FILE SAYS SO.
make_reference.py builds the two glycolipids from IUPAC-verified fragments and
this emitter READS that file rather than constructing a second string. That is
deliberate. Writing a second derivation of the same molecule here would not be
independent -- it would be me transcribing the same decision twice, which is
precisely F-037: the file written to fix a copied anchor re-copied the other
half, and the project then spent a round discovering that byte-identity is not
the signal.

The independence that IS real here is in the INPUTS, not in the string:

  from PIXELS   stage04_chains.py measures the drawn N-acyl chain off Fig. 1A
                with three methods. For C16 all three agree at 15 bonds = 16
                carbons, and chains.check_name_token confirms the token "C16"
                means the acyl chain and not something else.
  from the NAME make_reference.py builds N-hexadecanoyl-1-O-(alpha-D-
                galactopyranosyl)-D-erythro-sphingosine from fragments whose
                stereodescriptors are checked against their IUPAC names, each
                against a negative control of identical formula and mass.

Those two routes share no input. What they confirm is the acyl length and the
alpha/4E features. What NEITHER pins is the sphingoid base's length: the raster
gives 14 bonds by mid-line crossing and 13 by apex counting, so the base comes
from the reagent's name alone and EVERY ROW CARRYING IT IS FLAGGED
structure_assumed = 1, with the reason in curation_note. That is the flag doing
the job a QC file cannot (F-001).
"""
import csv
import json
import math
import os
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ROOT = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else ".")
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))
sys.path.insert(0, os.path.join(ROOT, "curation"))
sys.path.insert(0, os.path.join(ROOT, "provenance"))
import emit as EM                            # noqa: E402
import reagent_registry as RR                # noqa: E402
from corpus_components import load, smi      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = "KG_2025"

# ---------------------------------------------------------------------------
# The library design, Fig. 1A + 1B + the Results text.
# ---------------------------------------------------------------------------
BLOCKS = [
    ("ionizable_lipid", "SM-102", "Fig. 1A box", "library 1"),
    ("ionizable_lipid", "cKK-E12", "Fig. 1A box", "library 2"),
    ("sterol", "Cholesterol", "Fig. 1A box", "both libraries"),
    ("PEG_lipid", "C14-PEG2000", "Fig. 1A box", "both libraries"),
    ("PEG_lipid", "DMG-PEG", "Fig. 1A box", "both libraries"),
    ("PEG_lipid", "C18-PEG2000", "Fig. 1A box", "both libraries"),
    ("PEG_lipid", "ALC-0159", "Fig. 1A box", "both libraries"),
    ("helper_glycolipid", "C8 Galactosyl(a) Ceramide", "Fig. 1A, drawn in full",
     "both libraries"),
    ("helper_glycolipid", "C16 Galactosyl(a) Ceramide", "Fig. 1A, drawn in full",
     "both libraries"),
    ("helper_lipid_control", "DSPC", "Fig. 2C", "Control LNP61/62 only"),
]

# Fig. S5 (lead formulations) and Fig. 2C (the DSPC control).
FORMULATIONS = {
    "GlycoLNP45": dict(il="cKK-E12", il_mol=45.0, chol=27.5, peg="DMG-PEG",
                       peg_mol=2.5, hl="C16 Galactosyl(a) Ceramide",
                       hl_mol=25.0, dh=113, pdi=0.29, ee=26, mass_ratio="10:1",
                       src="Fig. S5"),
    "GlycoLNP59": dict(il="cKK-E12", il_mol=40.0, chol=27.5, peg="ALC-0159",
                       peg_mol=2.5, hl="C16 Galactosyl(a) Ceramide",
                       hl_mol=30.0, dh=80, pdi=0.14, ee=28, mass_ratio="10:1",
                       src="Fig. S5"),
    "GlycoLNP60": dict(il="cKK-E12", il_mol=45.0, chol=22.5, peg="ALC-0159",
                       peg_mol=2.5, hl="C16 Galactosyl(a) Ceramide",
                       hl_mol=30.0, dh=132, pdi=0.22, ee=18, mass_ratio="10:1",
                       src="Fig. S5"),
    "GlycoLNP61": dict(il="SM-102", il_mol=45.0, chol=27.5, peg="ALC-0159",
                       peg_mol=2.5, hl="C16 Galactosyl(a) Ceramide",
                       hl_mol=25.0, dh=115, pdi=0.08, ee=58, mass_ratio="10:1",
                       src="Fig. S5 and Fig. 2C"),
    "GlycoLNP62": dict(il="SM-102", il_mol=45.0, chol=32.5, peg="ALC-0159",
                       peg_mol=2.5, hl="C16 Galactosyl(a) Ceramide",
                       hl_mol=20.0, dh=84, pdi=0.01, ee=54, mass_ratio="10:1",
                       src="Fig. S5"),
    "Control LNP61": dict(il="SM-102", il_mol=45.0, chol=27.5, peg="ALC-0159",
                          peg_mol=2.5, hl="DSPC", hl_mol=25.0, dh=60,
                          pdi=None, ee=76, mass_ratio="10:1", src="Fig. 2C"),
}

# ---------------------------------------------------------------------------
# Every value in the paper that could conceivably be a delivery measurement.
# corpus_admissible is on EVERY row, with its reason.
# ---------------------------------------------------------------------------
ADMIT = ("yes", "percent of cells positive for the Thy1.1 protein reporter -- "
                "a protein-abundance delivery readout, matching the corpus's "
                "protein_adundance_normalized vocabulary")

VALUES = [
    # formulation, organ, cell type, value, units as printed, locus, tier, admissible
    dict(f="GlycoLNP61", organ="spleen", cell="macrophages (F4/80+)",
         v=78.0, printed="78", locus="Fig. 3C numeral, teal (solid trace)",
         tier="4_printed_numerals", adm=ADMIT, emit=True),
    dict(f="GlycoLNP61", organ="liver", cell="endothelial cells (CD31+)",
         v=54.0, printed="54", locus="Fig. 3C numeral, teal (solid trace)",
         tier="4_printed_numerals", adm=ADMIT, emit=True),
    dict(f="Control LNP61", organ="spleen", cell="macrophages (F4/80+)",
         v=66.0, printed="66", locus="Fig. 3C numeral, green (dotted trace)",
         tier="4_printed_numerals", adm=ADMIT, emit=True),
    dict(f="Control LNP61", organ="liver", cell="endothelial cells (CD31+)",
         v=81.0, printed="81 (Fig. 3C) / 82 (Results text)",
         locus="Fig. 3C numeral, green (dotted trace)",
         tier="4_printed_numerals", adm=ADMIT, emit=True),
    dict(f="GlycoLNP45", organ="liver", cell="endothelial cells (CD31+)",
         v=69.41, printed="not printed; bar mean",
         locus="Fig. 2B Liver panel", tier="5_pixels", adm=ADMIT, emit=True),
    dict(f="GlycoLNP45", organ="spleen", cell="macrophages (F4/80+)",
         v=35.80, printed="not printed; bar mean",
         locus="Fig. 2B Spleen panel", tier="5_pixels", adm=ADMIT, emit=True),
    dict(f="GlycoLNP59", organ="liver", cell="endothelial cells (CD31+)",
         v=69.41, printed="not printed; bar mean",
         locus="Fig. 2B Liver panel", tier="5_pixels", adm=ADMIT, emit=True),
    dict(f="GlycoLNP59", organ="spleen", cell="macrophages (F4/80+)",
         v=51.14, printed="not printed; bar mean",
         locus="Fig. 2B Spleen panel", tier="5_pixels", adm=ADMIT, emit=True),
    dict(f="GlycoLNP62", organ="liver", cell="endothelial cells (CD31+)",
         v=73.37, printed="not printed; bar mean",
         locus="Fig. 2B Liver panel", tier="5_pixels", adm=ADMIT, emit=True),
    dict(f="GlycoLNP62", organ="spleen", cell="macrophages (F4/80+)",
         v=77.84, printed="not printed; bar mean",
         locus="Fig. 2B Spleen panel", tier="5_pixels", adm=ADMIT, emit=True),

    # ---- admissible in KIND but NOT emitted, each with its reason ----------
    dict(f="GlycoLNP61", organ="lung", cell="not stated", v=13.0, printed="13",
         locus="Fig. 3C numeral", tier="4_printed_numerals",
         adm=("no", "its PAIR value (Control LNP61, lung) is 0%, and a 0 "
                    "cannot enter the corpus's per-experiment z-score of "
                    "log10. Emitting one half of the pair would put a "
                    "one-sided organ into the standardisation. No cell type "
                    "is stated for this organ either."), emit=False),
    dict(f="Control LNP61", organ="lung", cell="not stated", v=0.0, printed="0",
         locus="Fig. 3C numeral", tier="4_printed_numerals",
         adm=("no", "0% cannot be log10-transformed; the paper's own text "
                    "says neither particle delivered strongly to this organ"),
         emit=False),
    dict(f="GlycoLNP61", organ="kidney", cell="not stated", v=1.0, printed="1",
         locus="Fig. 3C numeral", tier="4_printed_numerals",
         adm=("no", "printed to the nearest whole percent at the very bottom "
                    "of the range, where 1 vs 0.5 vs 1.4 changes log10 by "
                    "0.45 decades; no cell type stated"), emit=False),
    dict(f="Control LNP61", organ="kidney", cell="not stated", v=1.0,
         printed="1", locus="Fig. 3C numeral", tier="4_printed_numerals",
         adm=("no", "same as the GlycoLNP61 kidney value"), emit=False),
    dict(f="GlycoLNP61", organ="heart", cell="not stated", v=0.0, printed="0",
         locus="Fig. 3C numeral", tier="4_printed_numerals",
         adm=("no", "0% cannot be log10-transformed"), emit=False),
    dict(f="Control LNP61", organ="heart", cell="not stated", v=0.0,
         printed="0", locus="Fig. 3C numeral", tier="4_printed_numerals",
         adm=("no", "0% cannot be log10-transformed"), emit=False),
    dict(f="GlycoLNP60", organ="n/a", cell="n/a", v=None, printed="none",
         locus="Fig. S5 / Fig. 2A only",
         tier="not_a_source",
         adm=("no", "FORMULATED AND CHARACTERISED BUT NEVER DOSED -- the text "
                    "excluded it on its 18% encapsulation efficiency, so it "
                    "has a complete formulation and no delivery value at all"),
         emit=False),
    dict(f="GlycoLNP62 / Control LNP62 pair", organ="spleen:liver", cell="n/a",
         v=None, printed="ratio only", locus="Fig. S8",
         tier="not_a_source",
         adm=("no", "A RATIO. Fig. S8 publishes only the spleen-to-liver "
                    "delivery ratio for this pair, never the two absolute "
                    "per-organ values, so Control LNP62 has no emittable "
                    "value and GlycoLNP62's own absolutes come from Fig. 2B "
                    "instead. CLAUDE.md section 8 already rules a selectivity "
                    "contrast is not a poolable readout."), emit=False),
    dict(f="all 109 screened LNPs", organ="various", cell="25 cell types",
         v=None, printed="none", locus="Fig. 1G/1H, Fig. S11",
         tier="not_a_source",
         adm=("no", "THE POOLED SCREEN PUBLISHES NO PER-LNP VALUE. Fig. 1G/1H "
                    "show normalized barcode counts aggregated into two "
                    "groups (C8 vs C16 glycolipid); Fig. S11 plots 52 points "
                    "that are (a) spleen/liver RATIOS and (b) unlabelled. "
                    "Fig. S4, captioned as the ranking of GlycoLNPs by "
                    "spleen-to-liver ratio, is a cartoon with no number and "
                    "no LNP name on it."), emit=False),
    dict(f="GlycoLNP61 dose response", organ="spleen", cell="macrophages",
         v=None, printed="not printed", locus="Fig. 3D",
         tier="5_pixels",
         adm=("no", "NOT EXTRACTED. Three doses (0.8/0.4/0.2 mg/kg) of a "
                    "formulation already emitted at 0.8 mg/kg; the 0.8 bar "
                    "here and the Fig. 3A/3C value are the same condition "
                    "measured in a different cohort, so emitting both would "
                    "put one formulation into the same experiment twice. See "
                    "DECISIONS D-4."), emit=False),
]


def z_scores(vals):
    """Per-experiment z of log10, population sd (ddof=0) -- the CLAUDE.md
    section 3 convention. ONE z for the whole experiment, pooling organs:
    a per-organ z is error class 22 and deletes the between-organ contrast
    that is this paper's entire finding."""
    logs = [math.log10(v) for v in vals]
    m = sum(logs) / len(logs)
    sd = math.sqrt(sum((x - m) ** 2 for x in logs) / len(logs))
    return [(x - m) / sd for x in logs], m, sd


def main():
    # ---- component structures -------------------------------------------
    rows_corpus = []
    for f in ("data_vitro.csv", "data_vivo.csv"):
        with open(os.path.join(ROOT, "data", f), newline="",
                  encoding="utf-8") as fh:
            rows_corpus += list(csv.DictReader(fh))
    comp, ambiguous = load(rows_corpus)

    # The glycolipid comes from the anchor file -- see this module's docstring.
    ref = {}
    with open(os.path.join(HERE, "reference_structures.csv"),
              newline="") as fh:
        for r in csv.DictReader(fh):
            if r["reference_SMILES"]:
                ref[r["lipid_name"]] = r["reference_SMILES"]

    STRUCT = {
        # registry first: these are catalogue reagents whose structure is
        # written from the chemical name, not inherited from the corpus.
        "SM-102": (RR.REGISTRY["SM-102"][0], "provenance/reagent_registry.py"),
        "DSPC": (RR.REGISTRY["DSPC"][0], "provenance/reagent_registry.py"),
        "Cholesterol": (RR.REGISTRY["Cholesterol"][0],
                        "provenance/reagent_registry.py"),
        # not in the registry -> the guarded corpus loader, which RAISES on an
        # ambiguous name rather than letting row order decide (error class 30).
        "cKK-E12": (smi(comp, ambiguous, "CKK-E12"),
                    "corpus_components.smi('CKK-E12'), unambiguous"),
        "ALC-0159": (smi(comp, ambiguous, "ALC-0159"),
                     "corpus_components.smi('ALC-0159'), unambiguous"),
        "DMG-PEG": (smi(comp, ambiguous, "DMG-PEG2000"),
                    "corpus_components.smi('DMG-PEG2000'); the paper prints "
                    "'DMG-PEG' with NO PEG block length, so the 2000 Da form "
                    "is ASSUMED"),
        "C16 Galactosyl(a) Ceramide":
            (ref["C16 Galactosyl(a) Ceramide"],
             "make_reference.py, from the chemical name via IUPAC-verified "
             "fragments; acyl length independently measured off Fig. 1A"),
    }

    # ---- 01_blocks.csv ---------------------------------------------------
    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["block_role", "component_name", "source_locus",
                    "used_in", "SMILES", "structure_provenance"])
        for role, name, locus, used in BLOCKS:
            s, prov = STRUCT.get(name, ("", "NOT RESOLVED -- this component "
                                            "appears in the design but in no "
                                            "formulation this curation emits"))
            w.writerow([role, name, locus, used, s, prov])

    # ---- 03_enumerated.csv ----------------------------------------------
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["formulation_id", "IL_name", "IL_molratio", "CHL_molratio",
                    "PEG_name", "PEG_molratio", "HL_name", "HL_molratio",
                    "molpct_sum", "Dh_nm", "PDI", "EE_pct",
                    "lipid_to_mRNA_mass_ratio", "source", "dosed"])
        for fid, F in FORMULATIONS.items():
            tot = F["il_mol"] + F["chol"] + F["peg_mol"] + F["hl_mol"]
            w.writerow([fid, F["il"], F["il_mol"], F["chol"], F["peg"],
                        F["peg_mol"], F["hl"], F["hl_mol"], tot, F["dh"],
                        F["pdi"], F["ee"], F["mass_ratio"], F["src"],
                        "no" if fid == "GlycoLNP60" else "yes"])

    # ---- 05_values.csv ---------------------------------------------------
    with open(os.path.join(HERE, "05_values.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["formulation", "organ", "cell_type", "value_percent",
                    "readout_as_printed", "units_as_printed", "source_locus",
                    "tier", "cargo", "model", "route", "dose", "timepoint",
                    "assay", "pooled", "pool_size", "corpus_admissible",
                    "admissible_reason", "emitted"])
        for V in VALUES:
            pooled = V["f"] == "all 109 screened LNPs"
            w.writerow([
                V["f"], V["organ"], V["cell"],
                "" if V["v"] is None else V["v"], V["printed"],
                "% Thy1.1+ cells", V["locus"], V["tier"],
                "Thy1.1 mRNA", "C57BL/6J mouse", "intravenous",
                "1.5 / 1.6 mg/kg total pool" if pooled else "0.8 mg/kg",
                "16 h", "flow cytometry" if not pooled else
                "DNA barcode deep sequencing",
                "yes" if pooled else "no",
                "52 (SM-102 library) / 57 (cKK-E12 library)" if pooled else "1",
                V["adm"][0], V["adm"][1], "yes" if V["emit"] else "no"])

    # ---- 06_rows.csv -----------------------------------------------------
    em = [V for V in VALUES if V["emit"]]
    zs, mu, sd = z_scores([V["v"] for V in em])

    records = {}
    for i, V in enumerate(em):
        F = FORMULATIONS[V["f"]]
        records[i] = {"formulation": V["f"], "organ": V["organ"],
                      "tier": V["tier"], "hl": F["hl"], "peg": F["peg"],
                      "locus": V["locus"], "printed": V["printed"]}

    def confidence(i, rec):
        if rec["hl"] != "DSPC":
            return "medium"          # sphingoid base length not pinned
        if rec["organ"] == "liver":
            return "medium"          # 81 in Fig. 3C vs 82 in the text
        return "high"

    def note(i, rec):
        bits = []
        if rec["tier"] == "5_pixels":
            bits.append("value digitised from Fig. 2B (no numeral is printed); "
                        "the same panel's GlycoLNP61 bars reproduce Fig. 3C's "
                        "printed 54 and 78 to +0.67 and -0.16 points")
        else:
            bits.append("value is a numeral printed inside Fig. 3C")
        if rec["hl"] != "DSPC":
            bits.append("the helper lipid's sphingoid base length is taken "
                        "from the reagent name, not measured: the drawing "
                        "gives 14 bonds by mid-line crossing and 13 by apex "
                        "counting, though its N-acyl chain measures 16 "
                        "carbons by all three methods")
        if rec["peg"] == "DMG-PEG":
            bits.append("the paper prints 'DMG-PEG' with no PEG block length; "
                        "the 2000 Da form is assumed")
        if rec["hl"] == "DSPC" and rec["organ"] == "liver":
            bits.append("Fig. 3C prints 81 for this value and the Results "
                        "text says 82; the figure numeral is emitted and the "
                        "conflict is recorded in 08_contradictions.csv")
        return "; ".join(bits) + "."

    quality = EM.emit_rows(
        items=list(range(len(em))), records=records,
        confidence=confidence, note=note,
        positional_ambiguity=lambda i, r: "0",
        structure_assumed=lambda i, r: "0" if (r["hl"] == "DSPC" and
                                               r["peg"] != "DMG-PEG") else "1",
        value_is_bound=lambda i, r: "",
        value_normalization_unit=lambda i, r: "experiment",
        label="KG_2025")

    hdr = ["Experiment_ID", "Formulation_ID", "IL_name", "IL_SMILES",
           "IL_molratio", "IL_to_nucleicacid_massratio", "HL_name",
           "HL_SMILES", "HL_molratio", "CHL_name", "CHL_SMILES",
           "CHL_molratio", "PEG_name", "PEG_SMILES", "PEG_molratio",
           "Model", "Model_type", "Model_target", "Route_of_administration",
           "Cargo", "Cargo_type", "Dose_ug_nucleicacid", "Experiment_method",
           "Experiment_batching", "Experiment_value",
           "Experiment_value_within_model", "Publication_PMID",
           "Publication_DOI"] + EM.columns()

    out = []
    for i, V in enumerate(em):
        F = FORMULATIONS[V["f"]]
        q = dict(quality[i])
        # within_model: Model_type is Mouse_B6 on every row, so the
        # re-standardisation inside (experiment x Model_type) is the identity.
        q["Experiment_value_within_model"] = round(zs[i], 6)
        row = {
            "Experiment_ID": EXP, "Formulation_ID": V["f"],
            "IL_name": F["il"], "IL_SMILES": STRUCT[F["il"]][0],
            "IL_molratio": F["il_mol"],
            "IL_to_nucleicacid_massratio": "",
            "HL_name": F["hl"], "HL_SMILES": STRUCT[F["hl"]][0],
            "HL_molratio": F["hl_mol"],
            "CHL_name": "Cholesterol",
            "CHL_SMILES": STRUCT["Cholesterol"][0], "CHL_molratio": F["chol"],
            "PEG_name": F["peg"], "PEG_SMILES": STRUCT[F["peg"]][0],
            "PEG_molratio": F["peg_mol"],
            "Model": "in_vivo", "Model_type": "Mouse_B6",
            "Model_target": V["organ"],
            "Route_of_administration": "intravenous",
            "Cargo": "mRNA", "Cargo_type": "Thy1.1",
            "Dose_ug_nucleicacid": "",
            "Experiment_method": "protein_adundance_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": round(zs[i], 6),
            "Publication_PMID": "41183194",
            "Publication_DOI": "10.1073/pnas.2409569122",
        }
        row.update(q)
        out.append(row)

    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=hdr)
        w.writeheader()
        w.writerows(out)

    # ---- checks ----------------------------------------------------------
    molpct = {fid: F["il_mol"] + F["chol"] + F["peg_mol"] + F["hl_mol"]
              for fid, F in FORMULATIONS.items()}
    chk = {
        "rows_emitted": len(out),
        "formulations_emitted": len({r["Formulation_ID"] for r in out}),
        "values_catalogued": len(VALUES),
        "values_admissible_yes": sum(1 for V in VALUES if V["adm"][0] == "yes"),
        "values_admissible_no": sum(1 for V in VALUES if V["adm"][0] == "no"),
        "molpct_sums": molpct,
        "molpct_all_exactly_100": all(abs(v - 100.0) < 1e-9
                                      for v in molpct.values()),
        "z_mean": round(sum(zs) / len(zs), 12),
        "z_sd_pop": round(math.sqrt(sum(x * x for x in zs) / len(zs)), 12),
        "log10_mean": round(mu, 6), "log10_sd_pop": round(sd, 6),
        "smiles_parse": {n: Chem.MolFromSmiles(s) is not None
                         for n, (s, _) in STRUCT.items()},
        "inchikeys": {n: Chem.MolToInchiKey(Chem.MolFromSmiles(s))
                      for n, (s, _) in STRUCT.items()},
    }
    print(json.dumps(chk, indent=1))
    print(f"METRIC rows={len(out)} molpct_100={chk['molpct_all_exactly_100']} "
          f"z_mean={chk['z_mean']} z_sd={chk['z_sd_pop']}")
    with open(os.path.join(HERE, "05_emit_checks.json"), "w") as fh:
        json.dump(chk, fh, indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
