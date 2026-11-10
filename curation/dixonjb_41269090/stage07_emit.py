#!/usr/bin/env python3
"""Stage 07 -- emit the admissible rows in the 79-column corpus schema.

The quality columns are written by curation/lib/emit.emit_rows, the one audited
writer, with PER-ITEM callables. A literal confidence raises (error class 29,
F-011/F-025/F-046) and a caveat with no note raises (F-024). Both refusals are
exercised on purpose by --seed-defect, in their own processes, and the verdict
is read from the RETURN CODE, never from stdout.

ONE z-block, EM_2025, over all 28 rows.

The first draft of this file used two blocks, one per cohort. That is wrong for
the reason F-016 records as error class 22: standardising sub-blocks separately
DELETES the offset between them, and the offset here is measurable rather than
hypothetical. LNP7 is in both cohorts, and at the one site that drains the
injection in BOTH (the popliteal node) they agree to within 11.4% -- 37.34% vs
33.52%, ratio 1.114. So the two cohorts sit on the same scale, one z is correct,
and value_normalization_unit is 'experiment'.

What DOES differ between the cohorts is which nodes drain the injection: study A
injected BOTH paws, study B only one hindlimb. That is a biological fact about
the site, not a scale offset, so it is carried in Model_target
('..._nondraining') rather than deleted by a separate standardisation. Without
it, LNP7 would carry two axillary values 14.4x apart with nothing in any column
to tell them apart.

Experiment_value follows the CLAUDE.md section 3 convention: per-experiment
z-score of log10, population sd (ddof=0).

Run:  python3 -I stage07_emit.py
      python3 -I stage07_emit.py --seed-defect literal_confidence
      python3 -I stage07_emit.py --seed-defect caveat_without_note
"""
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))
import emit                                            # noqa: E402
from errors import EmitError                           # noqa: E402

from rdkit import Chem, RDLogger                       # noqa: E402
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors  # noqa: E402

RDLogger.DisableLog("rdApp.*")

PMID = "41269090"
DOI = "10.1021/acsnano.5c12080"
LINK = "https://doi.org/10.1021/acsnano.5c12080"

SITE_TARGET = {
    "ALN": "lymph_node_axillary",
    "BLN": "lymph_node_brachial",
    "PLN": "lymph_node_popliteal",
    "PLN Contralateral": "lymph_node_popliteal_contralateral",
    "PLN Injected": "lymph_node_popliteal",
    "PLV": "lymphatic_vessel_popliteal",
}
OFF_TARGET = {"ALN", "BLN", "PLN Contralateral"}   # study B only

def target(site, study):
    """Model_target, carrying whether the site drains the injection.

    Study A injected both paws, so ALN, BLN and PLN all drain an injected paw.
    Study B injected one hindlimb, so its ALN, BLN and contralateral PLN do
    not. Collapsing the two would give LNP7 two axillary values 14.4x apart
    under one key.
    """
    t = SITE_TARGET[site]
    if study == "study_B" and site in OFF_TARGET:
        return t + "_nondraining"
    return t

SCHEMA72 = [
    "Index", "LNP_ID", "Experiment_ID", "Formulation_ID", "IL_name",
    "IL_SMILES", "IL_protonated_SMILES", "IL_head_name", "IL_head_SMILES",
    "IL_linker_name", "IL_linker_SMILES", "IL_tail1_name", "IL_tail1_SMILES",
    "IL_tail2_name", "IL_tail2_SMILES", "IL_tail3_name", "IL_tail3_SMILES",
    "IL_tail4_name", "IL_tail4_SMILES", "IL_molratio",
    "IL_to_nucleicacid_massratio", "IL_to_nucleicacid_chargeratio", "HL_name",
    "HL_SMILES", "HL_molratio", "CHL_name", "CHL_SMILES", "CHL_molratio",
    "PEG_name", "PEG_SMILES", "PEG_molratio", "fifthcomponent_name",
    "fifthcomponent_SMILES", "fifthcomponent_molratio", "Aqueous_buffer",
    "Dialysis_buffer", "Mixing_method", "Model", "Model_type", "Model_target",
    "Route_of_administration", "Cargo", "Cargo_type", "Dose_ug_nucleicacid",
    "Experiment_method", "Experiment_batching", "Experiment_value",
    "Publication_link", "Publication_DOI", "Publication_PMID", "Heavy.Atoms",
    "Rings", "Aromatic.Rings", "Rotatable.Bonds", "ExactMolWt",
    "Topological.Polar.Surface.Area", "Hydrogen.Bond.Donors",
    "Hydrogen.Bond.Acceptors", "LogP", "Molar.Refractivity",
    "Fraction.sp3.Carbons", "sp3.Carbons", "Nitrogen.Count",
    "Molecular.Weight", "has_ester", "has_carbonate", "has_disulfide",
    "Num_carbon_in_tail", "num_protonatable_nitrogens",
    "num_permanent_cationic_N", "formal_net_charge", "num_unsaturated_cc_bonds",
]

NOTE_COMMON = (
    "Value is a RATIO OF MEANS derived from SI Table {t}: the average aVHH+ "
    "cell count divided by the average LEC count. The paper's own figure plots "
    "a MEAN OF RATIOS (one point per independent experiment); the two "
    "statistics differ by up to 2.5x on one cell of this paper and only the "
    "ratio of means is derivable from published material. "
)


def descriptors(smi):
    m = Chem.MolFromSmiles(smi)
    pat = {"has_ester": "[CX3](=O)[OX2H0][#6]",
           "has_carbonate": "[OX2][CX3](=O)[OX2]",
           "has_disulfide": "[#16X2][#16X2]"}
    d = {
        "Heavy.Atoms": m.GetNumHeavyAtoms(),
        "Rings": rdMolDescriptors.CalcNumRings(m),
        "Aromatic.Rings": rdMolDescriptors.CalcNumAromaticRings(m),
        "Rotatable.Bonds": rdMolDescriptors.CalcNumRotatableBonds(m),
        "ExactMolWt": round(Descriptors.ExactMolWt(m), 6),
        "Topological.Polar.Surface.Area": round(rdMolDescriptors.CalcTPSA(m), 4),
        "Hydrogen.Bond.Donors": rdMolDescriptors.CalcNumHBD(m),
        "Hydrogen.Bond.Acceptors": rdMolDescriptors.CalcNumHBA(m),
        "LogP": round(Crippen.MolLogP(m), 4),
        "Molar.Refractivity": round(Crippen.MolMR(m), 4),
        "Fraction.sp3.Carbons": round(rdMolDescriptors.CalcFractionCSP3(m), 6),
        "sp3.Carbons": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C"
                           and a.GetHybridization() == Chem.HybridizationType.SP3),
        "Nitrogen.Count": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N"),
        "Molecular.Weight": round(Descriptors.MolWt(m), 6),
        "num_unsaturated_cc_bonds": sum(
            1 for b in m.GetBonds()
            if b.GetBondType() == Chem.BondType.DOUBLE
            and b.GetBeginAtom().GetSymbol() == "C"
            and b.GetEndAtom().GetSymbol() == "C"),
        "formal_net_charge": Chem.GetFormalCharge(m),
    }
    for k, sm in pat.items():
        d[k] = bool(m.HasSubstructMatch(Chem.MolFromSmarts(sm)))
    d["num_permanent_cationic_N"] = sum(
        1 for a in m.GetAtoms() if a.GetSymbol() == "N"
        and a.GetFormalCharge() == 1 and a.GetTotalNumHs() == 0
        and a.GetDegree() == 4)
    d["num_protonatable_nitrogens"] = sum(
        1 for a in m.GetAtoms() if a.GetSymbol() == "N"
        and a.GetFormalCharge() == 0 and not a.GetIsAromatic()
        and not any(n.GetSymbol() == "C" and
                    any(b.GetBondType() == Chem.BondType.DOUBLE and
                        b.GetOtherAtom(n).GetSymbol() == "O"
                        for b in n.GetBonds())
                    for n in a.GetNeighbors()))
    d["Num_carbon_in_tail"] = ""
    return d


def zblock(vals):
    lg = [math.log10(v) for v in vals]
    mu = sum(lg) / len(lg)
    sd = math.sqrt(sum((x - mu) ** 2 for x in lg) / len(lg))
    return [(x - mu) / sd for x in lg], mu, sd


def main(seed=None):
    forms = {r["formulation_id"]: r
             for r in csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
    vals = [r for r in csv.DictReader(open(os.path.join(HERE, "05_values.csv")))
            if r["corpus_admissible"] == "yes"]

    # ONE block. See the module docstring.
    byblock = {"EM_2025": list(vals)}

    rows, quality_in, records, zstats = [], [], {}, {}
    for exp, rs in sorted(byblock.items()):
        rs.sort(key=lambda r: (r["study"], r["site"], r["administration"]))
        z, mu, sd = zblock([float(r["value_percent"]) for r in rs])
        zstats[exp] = {"n": len(rs), "mean_log10": mu, "sd_log10_ddof0": sd}
        for r, zi in zip(rs, z):
            f = forms[r["administration"]]
            item = f"{exp}|{r['study']}|{r['administration']}|{r['site']}"
            tbl = "S2" if r["study"] == "study_A" else "S3/S4"
            off = (r["study"] == "study_B" and r["site"] in OFF_TARGET)
            records[item] = {
                "item": item, "formulation": r["administration"],
                "IL_name": f["IL_name"],
                "IL_structure_source": f["IL_structure_source"],
                "site": r["site"], "study": r["study"], "table": tbl,
                "off_target_control": off,
                "value_percent": r["value_percent"],
                "numerator": r["numerator_aVHH_pos"],
                "denominator": r["denominator_cells"],
            }
            quality_in.append(item)
            row = {c: "" for c in SCHEMA72}
            row.update({
                "Experiment_ID": exp,
                "Formulation_ID": f"{exp}_{r['administration']}",
                "IL_name": f["IL_name"], "IL_SMILES": f["IL_SMILES"],
                "IL_molratio": f["IL_molratio"],
                "IL_to_nucleicacid_massratio": 10.0,
                "HL_name": f["HL_name"], "HL_SMILES": f["HL_SMILES"],
                "HL_molratio": f["HL_molratio"],
                "CHL_name": f["CHL_name"], "CHL_SMILES": f["CHL_SMILES"],
                "CHL_molratio": f["CHL_molratio"],
                "PEG_name": f["PEG_name_corpus"], "PEG_SMILES": f["PEG_SMILES"],
                "PEG_molratio": f["PEG_molratio"],
                "Aqueous_buffer": "citrate", "Dialysis_buffer": "PBS",
                "Mixing_method": "microfluidics",
                "Model": "in_vivo", "Model_type": "Mouse_B6",
                "Model_target": target(r["site"], r["study"]),
                "Route_of_administration": "intradermal",
                "Cargo": "mRNA", "Cargo_type": "aVHH",
                "Dose_ug_nucleicacid": "",
                "Experiment_method": "protein_adundance_normalized",
                "Experiment_batching": "individual",
                "Experiment_value": round(zi, 6),
                "Publication_link": LINK, "Publication_DOI": DOI,
                "Publication_PMID": PMID,
            })
            row.update(descriptors(f["IL_SMILES"]))
            rows.append(row)

    # ---- the quality columns, via the one audited writer ------------------
    def confidence(item_id, rec):
        if seed == "literal_confidence":
            raise AssertionError("unreachable")
        src = rec["IL_structure_source"]
        if src == "reagent_registry":
            return "high"
        if src.startswith("corpus_components"):
            return "medium"
        raise EmitError(f"no confidence rule for structure source {src!r}")

    def note(item_id, rec):
        s = NOTE_COMMON.format(t=rec["table"])
        if rec["IL_name"] == "cKK-E12":
            s += ("Structure is the corpus string already carried by AC_2025, "
                  "JW_2024 and BL_2023 for this same molecule; it declares no "
                  "stereochemistry, while cKK-E12 is lysine-derived and the "
                  "commercial reagent (BroadPharm BP-29590) is a single "
                  "diastereomer, so the two diketopiperazine centres are "
                  "undeclared. confidence=medium records that, not a doubt "
                  "about the constitution. ")
        else:
            s += "Benchmark structure taken from provenance/reagent_registry.py. "
        s += ("Cohort: %s (%s). " % (
            rec["study"],
            "both paws injected, Figure 2 / Table S2" if rec["study"] == "study_A"
            else "one hindlimb injected, Figure 4 / Tables S3-S4"))
        if rec["off_target_control"]:
            s += ("This node does NOT drain the injection site in that cohort, "
                  "which Model_target records as _nondraining; it is the "
                  "paper's own off-target control. ")
        if seed == "caveat_without_note":
            return ""
        return s.strip()

    conf_arg = "high" if seed == "literal_confidence" else confidence
    q = emit.emit_rows(quality_in, records, confidence=conf_arg, note=note,
                       positional_ambiguity=lambda i, r: 0,
                       structure_assumed=lambda i, r: 0,
                       value_is_bound=lambda i, r: "",
                       value_normalization_unit=lambda i, r: "experiment",
                       within_model=lambda i, r: "",
                       id_key="item", label="dixonjb_41269090")

    qcols = emit.columns()
    for row, qq in zip(rows, q):
        for c in qcols:
            v = qq.get(c, "")
            row[c] = "" if v is None else v
        # Model_type is a single value in each block, so the within-model
        # restandardisation is the identity. Stated, not assumed.
        row["Experiment_value_within_model"] = row["Experiment_value"]

    header = SCHEMA72 + qcols
    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header)
        w.writeheader()
        w.writerows(rows)

    print(f"06_rows.csv   {len(rows)} rows x {len(header)} columns")
    for e, s in sorted(zstats.items()):
        print(f"  {e}: n={s['n']}  mean(log10 %)={s['mean_log10']:.6f}  "
              f"sd(ddof=0)={s['sd_log10_ddof0']:.6f}")
    json.dump({"zblocks": zstats}, open(os.path.join(HERE, "06_rows_z.json"), "w"),
              indent=1)
    return 0


if __name__ == "__main__":
    s = None
    if "--seed-defect" in sys.argv:
        s = sys.argv[sys.argv.index("--seed-defect") + 1]
    try:
        sys.exit(main(s))
    except EmitError as e:
        print(f"emit REFUSED: {e}", file=sys.stderr)
        sys.exit(3)
