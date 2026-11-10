#!/usr/bin/env python3
"""Stage 02 -- what this paper shares with the live corpus. Writes 02_corpus_overlap.json.

    python3 -I 02_corpus_overlap.py

Read-only against data/data_vitro.csv and data/data_vivo.csv.

Three questions, kept separate because they answer differently:
  1. the TOKEN-PATTERN TEST -- does the paper's label '7C1' mean anything in the
     corpus's own naming? (It does not, and a substring join is a trap.)
  2. a SUBSTRUCTURE SWEEP, labelled honestly as a SUPERSET.
  3. the EXCIPIENTS, compared by InChIKey rather than by name.
"""
import json, os, sys
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))
import reagents                                            # noqa: E402

VITRO = os.path.join(ROOT, "data", "data_vitro.csv")
VIVO = os.path.join(ROOT, "data", "data_vivo.csv")

# The arm the paper DRAWS for 7C1 in Fig 3A / Fig 5A: a tertiary amine bearing
# CH2-CH(OH)-C13. Built here as a SMARTS, deliberately loose, and reported as a
# superset. The real 7C1 is a polymer; this pattern cannot identify it.
ARM_SMARTS = "[NX3;!$(N=*)][CH2][CH1]([OX2H1])[CH2]"
ARM_STRICT_SMARTS = ("[NX3;!$(N=*)][CH2][CH1]([OX2H1])"
                     "[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH3]")


def token_test(dv, dvv):
    """Would a substring join on the paper's own label '7C1' find anything?"""
    out = {}
    tot_rows = tot_struct = 0
    exps = set()
    for pool, d in (("vitro", dv), ("vivo", dvv)):
        m = d["IL_name"].astype(str).str.contains("7C1", case=False, na=False)
        sub = d[m]
        out[pool] = {
            "rows_matched": int(m.sum()),
            "distinct_structures": int(sub["IL_SMILES"].nunique()),
            "experiments": sorted(sub["Experiment_ID"].unique().tolist()),
            "examples": sub[["Experiment_ID", "IL_name", "IL_SMILES"]]
                        .drop_duplicates().head(4).values.tolist(),
        }
        tot_rows += int(m.sum())
        tot_struct += int(sub["IL_SMILES"].nunique())
        exps |= set(sub["Experiment_ID"].unique().tolist())
    out["total_rows"] = tot_rows
    out["total_distinct_structures"] = tot_struct
    out["total_experiments"] = sorted(exps)
    out["verdict"] = (
        "A substring join on the bare token '7C1' captures %d corpus rows spanning %d "
        "distinct structures across %d experiments, and NOT ONE of them is Dahlman's "
        "7C1 lipomer: they are combinatorial labels in which the characters happen to "
        "co-occur (KL_2010 '117C14' = head 117 + C14; LL_2012 'A7C10'; LM_2019 "
        "'A1Iso7C11'; ZH_2023 'A4I18R7C18-2'). The paper's own '7C1' is a compound "
        "index from Dahlman 2014, not a composition: the '7' is the amine and the "
        "'C1' is the epoxide index, while the drawn tail is C13 and the epoxide is "
        "C15. Three different numbers, none of them a chain length."
        % (tot_rows, tot_struct, len(exps)))
    return out


def substructure(dv, dvv):
    loose = Chem.MolFromSmarts(ARM_SMARTS)
    strict = Chem.MolFromSmarts(ARM_STRICT_SMARTS)
    assert loose is not None and strict is not None
    res = {}
    for pool, d in (("vitro", dv), ("vivo", dvv)):
        smis = sorted(set(d["IL_SMILES"].dropna().astype(str)))
        nl = ns = 0
        hits_strict = []
        for s in smis:
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                continue
            if mol.HasSubstructMatch(loose):
                nl += 1
                if mol.HasSubstructMatch(strict):
                    ns += 1
                    hits_strict.append(s)
        res[pool] = {"unique_structures_searched": len(smis),
                     "loose_hits_SUPERSET": nl,
                     "strict_hits_full_C13_arm": ns,
                     "strict_examples": hits_strict[:5]}
    res["_what"] = (
        "SUPERSET, not an identification. The loose pattern is 'any amine bearing "
        "CH2-CH(OH)-CH2', the amino-alcohol motif an epoxide ring-opening leaves; it "
        "matches a whole chemotype and says nothing about 7C1. The strict pattern "
        "additionally requires the full 13-carbon tail the paper draws. Neither can "
        "identify 7C1, because 7C1 is a polydisperse polymer with no single SMILES: a "
        "strict hit would be a small molecule carrying the same arm, not this lipid.")
    return res


def excipients(dv, dvv):
    """Compare by InChIKey, not by name (324 corpus structures carry >1 IL_name)."""
    drawn = {
        "Cholesterol": "C[C@H](CCCC(C)C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC"
                       "[C@]4(C)[C@H]3CC[C@]12C",
        "DSPC": "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)"
                "OC(=O)CCCCCCCCCCCCCCCCC",
    }
    out = {}
    for name, smi in drawn.items():
        ours_full = reagents.full_key(smi)
        ours_skel = reagents.skeleton(smi)
        found = {}
        for pool, d, col, scol in (("vitro", dv, "CHL_name", "CHL_SMILES"),
                                   ("vivo", dvv, "CHL_name", "CHL_SMILES"),
                                   ("vitro", dv, "HL_name", "HL_SMILES"),
                                   ("vivo", dvv, "HL_name", "HL_SMILES")):
            sub = d[d[col].astype(str) == name]
            if not len(sub):
                continue
            keys = {reagents.full_key(s) for s in set(sub[scol].dropna().astype(str))}
            found[f"{pool}:{col}"] = {
                "rows": int(len(sub)),
                "distinct_corpus_structures": len(keys),
                "full_inchikeys": sorted(keys),
                "matches_our_drawn_structure": ours_full in keys,
            }
        reg = reagents.registered(name)
        out[name] = {"our_smiles_from_the_drawn_figure": smi,
                     "our_full_inchikey": ours_full,
                     "our_skeleton": ours_skel,
                     "registry_entry": reg,
                     "matches_registry": (reg is not None
                                          and reagents.full_key(reg) == ours_full),
                     "corpus": found}
    return out


def peg_names(dv, dvv):
    names = sorted(set(dv["PEG_name"].dropna().astype(str))
                   | set(dvv["PEG_name"].dropna().astype(str)))
    return {"corpus_PEG_name_vocabulary": names,
            "this_paper_uses": ["C14PEG2000", "C16PEG2000", "C18PEG2000",
                                "C14PEG3000", "C18PEG3000"],
            "already_in_the_vocabulary": [],
            "note": "The paper names no catalogue product. Fig 3B draws a diacyl "
                    "glycerophosphoethanolamine-PEG, which for C14/C16/C18 would be "
                    "DMPE-/DPPE-/DSPE-PEG2000 -- two of those three names already "
                    "exist in the corpus. The mapping is an INFERENCE from a drawn "
                    "skeleton, so it is recorded and not applied; and PEG-lipids are "
                    "deliberately absent from provenance/reagent_registry.py anyway, "
                    "because the block length is a distribution."}


def shared_invivo(dvv):
    """The in-vivo pool's cross-experiment chemistry, for context."""
    g = dvv.groupby("IL_SMILES")["Experiment_ID"].nunique()
    shared = g[g > 1]
    rows = []
    for smi, n in shared.items():
        sub = dvv[dvv["IL_SMILES"] == smi]
        rows.append({"n_experiments": int(n),
                     "experiments": sorted(sub["Experiment_ID"].unique().tolist()),
                     "names": sorted(set(sub["IL_name"].astype(str))),
                     "smiles": smi})
    return {"n_structures_shared_across_experiments": int(len(shared)),
            "structures": rows,
            "does_this_paper_add_one": False,
            "why": "This paper's only ionizable lipid is 7C1, which has no SMILES and "
                   "is absent from the corpus under any name. It adds no shared "
                   "in-vivo chemistry. It also uses no benchmark ionizable lipid at "
                   "all -- no MC3, no SM-102, no ALC-0315 anywhere in the paper; its "
                   "only comparator is its own 80:20 7C1 formulation."}


def main():
    dv = pd.read_csv(VITRO, low_memory=False)
    dvv = pd.read_csv(VIVO, low_memory=False)
    out = {
        "_what": "What curation/dahlmanje_30394729 (Sago 2018, PMID 30394729) shares "
                 "with the live corpus. Produced by 02_corpus_overlap.py.",
        "corpus_shape": {"vitro": list(dv.shape), "vivo": list(dvv.shape),
                         "vitro_experiments": int(dv["Experiment_ID"].nunique()),
                         "vivo_experiments": int(dvv["Experiment_ID"].nunique())},
        "token_pattern_test": token_test(dv, dvv),
        "substructure_sweep_SUPERSET": substructure(dv, dvv),
        "excipients_by_inchikey": excipients(dv, dvv),
        "peg_lipid_names": peg_names(dv, dvv),
        "in_vivo_shared_chemistry": shared_invivo(dvv),
        "ionizable_lipids_this_paper_would_add": 0,
        "ionizable_lipids_this_paper_uses": 1,
    }
    with open(os.path.join(HERE, "02_corpus_overlap.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("excipients_by_inchikey",)}, indent=1)[:4000])
    print(json.dumps(out["excipients_by_inchikey"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
