#!/usr/bin/env python3
"""Stage 03 -- the component inventory and the (one-member) IL enumeration.

RUN: python3 -I scripts/stage03_structures.py <vitro.csv> <01_blocks.csv> <03_enumerated.csv> <out.json>

THERE IS NO ENUMERATION TO DO. This paper contains exactly ONE ionizable-lipid
structure, DLin-MC3-DMA, on all 1,080 training rows and all 72 validation rows.
What varies is the FORMULATION. So 03_enumerated.csv has one row, and the real
structural work is the excipient inventory in 01_blocks.csv.

Every structure is resolved through curation/lib/reagents.smi(), whose
precedence rule 1 is provenance/reagent_registry.py -- a file that writes each
structure FROM ITS CHEMICAL NAME rather than inheriting whatever the corpus
happens to hold. No SMILES is typed into this script. That rule exists because
`18PG` became symmetric 18:0/18:0 when someone typed it, and this paper is
about to show that `18PG` is contested even between two papers from one lab.

The one structure NOT resolved by name is this paper's `18PG`. The registry
pins `18PG` to SOPG (18:0/18:1) on the authority of YZ_2022's Materials
sentence; this paper DRAWS 18:1/18:1 and prints Double_bonds = 2. So the
structure is taken from the registry entry whose CHEMICAL IDENTITY matches what
the paper draws -- `DOPG` -- and the name/structure disagreement is carried on
the row and raised as a decision, never silently resolved. That is the AC_2025
rule: when a name and a structure disagree, queue it.
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "lib"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "provenance"))

from rdkit import Chem                                     # noqa: E402
from rdkit import RDLogger                                 # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors       # noqa: E402

RDLogger.DisableLog("rdApp.*")

import reagents                                            # noqa: E402

# role, published name, the REGISTRY KEY used to resolve it, why that key
COMPONENTS = [
    ("IL", "DLin-MC3-DMA", "Dlin-MC3-DMA",
     "Named in Methods (MedKoo Biosciences) and drawn in SI Table 2; the "
     "registry alias maps the paper's casing to the registry key."),
    ("CHL", "Cholesterol", "Cholesterol",
     "Named in Methods (Sigma-Aldrich) and drawn in SI Table 2."),
    ("HL", "DOTAP", "DOTAP",
     "Named in Methods (Avanti) and drawn in SI Table 2 as 18:1/18:1 TAP "
     "chloride; the drawing agrees with the registry."),
    ("HL", "DDAB", "DDAB",
     "Named in Methods (Avanti) and drawn in SI Table 2 as "
     "dimethyldioctadecylammonium bromide; agrees with the registry."),
    ("HL", "DOPE", "DOPE",
     "Named in Methods (Avanti) and drawn in SI Table 2 as 18:1/18:1 PE; "
     "agrees with the registry."),
    ("HL", "DSPC", "DSPC",
     "Named in Methods (Avanti) and drawn in SI Table 2 as 18:0/18:0 PC; "
     "agrees with the registry."),
    ("HL", "14PA", "14PA",
     "Named in Methods (Avanti) and drawn in SI Table 2 as 14:0/14:0 PA, "
     "sodium salt; agrees with the registry."),
    ("HL", "18PG", "DOPG",
     "CONTESTED. Methods names 'Avanti 18PG'. SI Table 2 DRAWS a cis C=C in "
     "BOTH tails and SI Table 3 prints Double_bonds = 2, which is DOPG "
     "(18:1/18:1). The registry's `18PG` key is SOPG (18:0/18:1, one C=C) on "
     "the authority of YZ_2022's Materials sentence for THIS SAME 1,080-"
     "formulation library. Resolved to the structure the paper draws, via the "
     "registry's own DOPG entry WITH THE 1'-HEAD STEREOCENTRE CLEARED, "
     "because the drawing leaves the head glycerol undefined (rac) and the "
     "registry DOPG entry defines it; see drawn_18PG(). Flagged on every "
     "emitted row and raised as DECISIONS.md D-2."),
]

# DMG-PEG2000 is not given a SMILES *here*, in the block inventory, because the
# reagent registry excludes PEG lipids on purpose (CLAUDE.md section 3): the
# block length is a distribution and no single SMILES is the reagent. The
# emitted rows still carry a PEG_SMILES, because the corpus column does and
# because one structure per reagent name is what C17 gates -- and the structure
# the paper draws is byte-identical to the one YZ_2022/YZ_2024 already store.
PEG_NAME = "DMG-PEG2000"


def drawn_18PG(qc):
    """This paper's 18PG: dioleoyl backbone, and a 1'-RAC head glycerol.

    SI Table 2 draws a bold wedge to H on the glycerol BACKBONE C2 and a plain,
    unwedged OH on the HEAD glycerol C2' (renders/p2_x68_head.png, magnified
    4x from the native 300 dpi embedded image). So the head stereocentre is
    undefined in the drawing -- which is also what Avanti sells: their PG
    products are 1'-rac-glycerol.

    The registry's DOPG entry carries a DEFINED 1'-sn head, so emitting it
    verbatim would add a stereocentre the paper does not draw -- the same class
    of silent stereo error as F-013, in the opposite direction. So the head tag
    is CLEARED, and the clearing is checked rather than trusted: the skeleton
    InChIKey must be unchanged, the formula and exact mass must be unchanged,
    the count of specified stereocentres must fall by exactly one, and the
    survivor must be the backbone carbon. No SMILES is typed.
    """
    base = reagents.smi("DOPG")
    m = Chem.MolFromSmiles(base)
    before = Chem.FindMolChiralCenters(m, includeUnassigned=True,
                                       useLegacyImplementation=False)
    # the head stereocentre is the one whose neighbours include a free -OH and
    # a -CH2OH; identify it by substructure, never by atom index.
    patt = Chem.MolFromSmarts("[OX2]P(=O)([OX1,OX2])OC[CH1]([OX2H1])C[OX2H1]")
    hit = m.GetSubstructMatches(patt)
    if len(hit) != 1:
        raise SystemExit("drawn_18PG: expected one head-glycerol match, got %d"
                         % len(hit))
    head_c = [a for a in hit[0]
              if m.GetAtomWithIdx(a).GetSymbol() == "C"
              and m.GetAtomWithIdx(a).GetChiralTag()
              != Chem.ChiralType.CHI_UNSPECIFIED]
    if len(head_c) != 1:
        raise SystemExit("drawn_18PG: expected one tagged head carbon, got %d"
                         % len(head_c))
    m.GetAtomWithIdx(head_c[0]).SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    out = Chem.MolToSmiles(m)
    m2 = Chem.MolFromSmiles(out)
    after = Chem.FindMolChiralCenters(m2, includeUnassigned=True,
                                      useLegacyImplementation=False)
    nb = sum(1 for _, c in before if c != "?")
    na = sum(1 for _, c in after if c != "?")
    k1 = Chem.MolToInchiKey(Chem.MolFromSmiles(base))
    k2 = Chem.MolToInchiKey(m2)
    chk = {
        "base": base, "base_full_InChIKey": k1,
        "drawn": out, "drawn_full_InChIKey": k2,
        "skeleton_unchanged": k1.split("-")[0] == k2.split("-")[0],
        "specified_stereocentres_before": nb,
        "specified_stereocentres_after": na,
        "fell_by_exactly_one": nb - na == 1,
        "formula_unchanged": (rdMolDescriptors.CalcMolFormula(
            Chem.MolFromSmiles(base)) == rdMolDescriptors.CalcMolFormula(m2)),
        "exact_mass_unchanged": abs(
            Descriptors.ExactMolWt(Chem.MolFromSmiles(base))
            - Descriptors.ExactMolWt(m2)) < 1e-9,
        "surviving_centres": [c for _, c in after if c != "?"],
    }
    qc["drawn_18PG_derivation"] = chk
    for k in ("skeleton_unchanged", "fell_by_exactly_one",
              "formula_unchanged", "exact_mass_unchanged"):
        if not chk[k]:
            raise SystemExit("drawn_18PG: check %s failed: %s" % (k, chk))
    if na != 1:
        raise SystemExit("drawn_18PG: expected exactly one surviving "
                         "stereocentre (the backbone C2), got %d" % na)
    return out


def main(vitro, blocks_out, enum_out, qc_out):
    Q = {"resolved": [], "unresolved": []}
    rows = []
    for role, published, key, why in COMPONENTS:
        smi = drawn_18PG(Q) if published == "18PG" else reagents.smi(key)
        m = Chem.MolFromSmiles(smi)
        rows.append({
            "role": role,
            "published_name": published,
            "registry_key_used": key,
            "SMILES": Chem.MolToSmiles(m),
            "full_InChIKey": Chem.MolToInchiKey(m),
            "skeleton_InChIKey": Chem.MolToInchiKey(m).split("-")[0],
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "exact_mass": round(Descriptors.ExactMolWt(m), 4),
            "name_equals_registry_key": published == key,
            "resolution_note": why,
        })
        Q["resolved"].append(published)
    rows.append({
        "role": "PEG", "published_name": PEG_NAME, "registry_key_used": "",
        "SMILES": "", "full_InChIKey": "", "skeleton_InChIKey": "",
        "formula": "", "exact_mass": "",
        "name_equals_registry_key": "",
        "resolution_note": ("Deliberately left without a SMILES. CLAUDE.md "
                            "section 3: PEG lipids are absent from the "
                            "reagent registry on purpose because the block "
                            "length is a distribution, so no single SMILES is "
                            "right. SI Table 2 draws 45 ethylene oxides -- "
                            "the bracket subscript is 44 and one more is "
                            "drawn outside it -- which is ONE member of that "
                            "distribution and happens to be exactly the "
                            "structure the corpus already stores for "
                            "DMG-PEG2000 in YZ_2022/YZ_2024. 06_rows.csv "
                            "emits that structure so the reagent keeps one "
                            "structure across experiments (C17)."),
    })
    Q["unresolved"].append(PEG_NAME)

    with open(blocks_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---------------------------------------------------- the IL enumeration
    il = [r for r in rows if r["role"] == "IL"][0]
    # descriptor columns: taken from the corpus's OWN rows for this exact
    # structure, i.e. build_dataset.py's own output for this SMILES, rather
    # than re-implemented here. Re-implementing a descriptor definition is how
    # it drifts; and all rows carrying this structure must agree, which is
    # checked rather than assumed.
    DESC = ["Heavy.Atoms", "Rings", "Aromatic.Rings", "Rotatable.Bonds",
            "ExactMolWt", "Topological.Polar.Surface.Area",
            "Hydrogen.Bond.Donors", "Hydrogen.Bond.Acceptors", "LogP",
            "Molar.Refractivity", "Fraction.sp3.Carbons", "sp3.Carbons",
            "Nitrogen.Count", "Molecular.Weight", "has_ester",
            "has_carbonate", "has_disulfide", "Num_carbon_in_tail",
            "num_protonatable_nitrogens", "num_permanent_cationic_N",
            "formal_net_charge", "num_unsaturated_cc_bonds"]
    with open(vitro, newline="") as fh:
        corpus = list(csv.DictReader(fh))
    # Narrow by NAME first (cheap), then confirm every survivor by FULL
    # InChIKey. The name is only a filter -- CLAUDE.md section 3 is explicit
    # that 288 IL names carry different structures in different experiments,
    # so a name hit that fails the key check is dropped, not trusted.
    cand = [r for r in corpus if r["IL_name"].strip()
            in ("DLin-MC3-DMA", "Dlin-MC3-DMA") and r["IL_SMILES"].strip()]
    same, name_hit_wrong_structure = [], 0
    for r in cand:
        m = Chem.MolFromSmiles(r["IL_SMILES"])
        if m is not None and Chem.MolToInchiKey(m) == il["full_InChIKey"]:
            same.append(r)
        else:
            name_hit_wrong_structure += 1
    desc_sets = {c: sorted(set(r[c] for r in same)) for c in DESC}
    disagree = {c: v for c, v in desc_sets.items() if len(v) != 1}
    Q["descriptors_from_corpus"] = {
        "corpus_rows_with_this_exact_full_InChIKey": len(same),
        "name_hits_whose_structure_is_NOT_this_one": name_hit_wrong_structure,
        "experiments": sorted(set(r["Experiment_ID"] for r in same)),
        "columns_checked": len(DESC),
        "columns_with_more_than_one_value": disagree,
        "all_agree": len(disagree) == 0,
        "why": ("These are build_dataset.py's own descriptor values for this "
                "exact structure. Reusing them is not a shortcut: "
                "re-implementing a descriptor definition in a curation script "
                "is how it drifts from the builder's. Agreement across every "
                "corpus row carrying the structure is the check."),
    }
    if disagree:
        print(json.dumps({"FATAL": "corpus descriptors disagree",
                          "detail": disagree}, indent=1))
        return 1

    erow = {
        "IL_id": "MC3",
        "IL_name": "DLin-MC3-DMA",
        "IL_SMILES": il["SMILES"],
        "IL_full_InChIKey": il["full_InChIKey"],
        "IL_formula": il["formula"],
        "IL_exact_mass": il["exact_mass"],
        "n_formulations_training": 1080,
        "n_formulations_validation": 72,
        "structure_source": ("provenance/reagent_registry.py, written from "
                             "the chemical name "
                             "(6Z,9Z,28Z,31Z)-heptatriaconta-6,9,28,31-"
                             "tetraen-19-yl 4-(dimethylamino)butanoate"),
        "paper_draws_it": "SI Table 2, page 2, image xref 61 (renders/p2_x61.png)",
        "confidence_basis": ("A catalogue reagent (MedKoo), named in Methods, "
                             "drawn in full in the SI, and independently "
                             "re-derived as the C10 anchor by make_reference.py"),
    }
    for c in DESC:
        erow[c] = desc_sets[c][0]
    with open(enum_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(erow.keys()))
        w.writeheader()
        w.writerow(erow)

    Q["n_distinct_ionizable_lipid_structures"] = 1
    Q["n_distinct_helper_lipid_structures"] = 6
    Q["n_components_with_a_structure"] = len(rows) - 1
    Q["n_components_without_a_structure"] = 1
    with open(qc_out, "w") as fh:
        json.dump(Q, fh, indent=1)
    print(json.dumps({"blocks": blocks_out, "enumerated": enum_out,
                      "descriptors_all_agree": Q["descriptors_from_corpus"]["all_agree"],
                      "corpus_rows_matched": len(same)}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
