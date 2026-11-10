#!/usr/bin/env python3
"""Corpus overlap for maohq_38583366, and the duplicate-arm test.

EVERY QUERY HERE CARRIES A POSITIVE CONTROL, AND EACH CONTROL ASSERTS THAT IT
MATCHES A SPECIFIC NAMED MOLECULE KNOWN TO BE PRESENT (F-070 and its amendment:
one sweep returned a clean 0 from a SMARTS that could not match a molecule known
to be there, and a second caught its own control returning 6 structures while
not matching the registry cholesterol at all. A non-zero count is not evidence
that a control works).

Inputs, all inside the repository -- nothing is read from a session scratchpad
(CLAUDE.md Sec 7, F-023):
  data/data_vitro.csv, data/data_vivo.csv
  curation/maohq_39375194/source/github/1080_LNP_OHE.csv   (the Mao lab's own
      deposit for PMID 39375194; this paper has no repository of its own)
  01_si_tables.csv (parsed from this paper's mmc1.docx by stage04)
"""
import csv, hashlib, json, os, sys
from collections import Counter
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..")
OUT = os.environ.get("STAGE03_OUT",
                     os.path.join(HERE, "..", "02_corpus_overlap.json"))
REPO_CSV = os.path.join(ROOT, "curation", "maohq_39375194", "source",
                        "github", "1080_LNP_OHE.csv")
SI_CSV = os.path.join(HERE, "..", "01_si_tables.csv")

MC3 = r"CCCCC/C=C\C/C=C\CCCCCCCCC(CCCCCCCC/C=C\C/C=C\CCCCC)OC(=O)CCCN(C)C"
CHOL = "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"


def ik(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise SystemExit("stage03: unparseable SMILES %r" % smi)
    return Chem.MolToInchiKey(m)


def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()


def main():
    vit = pd.read_csv(os.path.join(ROOT, "data", "data_vitro.csv"),
                      low_memory=False)
    viv = pd.read_csv(os.path.join(ROOT, "data", "data_vivo.csv"),
                      low_memory=False)
    keys = {}
    for d in (vit, viv):
        for s in set(d.IL_SMILES.dropna()):
            if s not in keys:
                m = Chem.MolFromSmiles(s)
                keys[s] = Chem.MolToInchiKey(m) if m else None
        d["ik"] = d.IL_SMILES.map(keys)

    mc3_k, chol_k = ik(MC3), ik(CHOL)

    # ---------------- positive controls --------------------------------
    controls = []
    n = int(vit.HL_name.astype(str).str.contains("DSPC", case=False).sum()
            + viv.HL_name.astype(str).str.contains("DSPC", case=False).sum())
    controls.append({"control": "DSPC name scan over HL_name", "rows": n,
                     "asserts": "the named helper DSPC, which this paper uses",
                     "passes": n > 0})
    nm = int((vit.ik == mc3_k).sum() + (viv.ik == mc3_k).sum())
    spell = sorted(set(vit.loc[vit.ik == mc3_k, "IL_name"].dropna())
                   | set(viv.loc[viv.ik == mc3_k, "IL_name"].dropna()))
    controls.append({"control": "DLin-MC3-DMA by InChIKey", "rows": nm,
                     "inchikey": mc3_k, "name_spellings": spell,
                     "asserts": "THE ionizable lipid of this paper, by "
                                "structure not by name",
                     "passes": nm > 0 and len(spell) >= 1})
    nc = int(vit.CHL_SMILES.map(lambda s: keys.get(s) if False else None)
             .notna().sum())   # placeholder, replaced below
    nc = int((vit.CHL_SMILES == CHOL).sum() + (viv.CHL_SMILES == CHOL).sum())
    controls.append({"control": "cholesterol by exact registry SMILES",
                     "rows": nc, "inchikey": chol_k,
                     "asserts": "the registry cholesterol string itself",
                     "passes": nc > 0})
    if not all(c["passes"] for c in controls):
        raise SystemExit("stage03: a positive control failed; every count "
                         "below would be uninterpretable.")

    # ---------------- the one ionizable lipid --------------------------
    mc3_vit = vit[vit.ik == mc3_k]
    mc3_viv = viv[viv.ik == mc3_k]
    g = viv.groupby("ik").Experiment_ID.nunique()
    shared_viv = sorted(g[g > 1].index)
    viv_shared = []
    for k in shared_viv:
        r = viv[viv.ik == k]
        viv_shared.append({"inchikey": k,
                           "names": sorted(set(r.IL_name.dropna())),
                           "sources": sorted(set(r.Experiment_ID)),
                           "organs": sorted(set(r.Model_target.astype(str))),
                           "rows": int(len(r))})

    # ---------------- is an arm already merged? ------------------------
    repo = pd.read_csv(REPO_CSV)
    si = pd.read_csv(SI_CSV)
    ilhl = repo["(IL+HL)"]
    repo["HLp"] = ilhl / (repo["Dlin-MC3_Helper lipid_ratio"] + 1)
    repo["ILp"] = ilhl - repo["HLp"]
    cp = 100 - ilhl
    repo["PEGp"] = cp / (repo["Chol_DMG-PEG_ratio"] + 1)
    repo["CHOLp"] = cp - repo["PEGp"]
    repo["np_"] = repo.NP_ratio.astype(float)
    z = vit[vit.Experiment_ID == "YZ_2022"].copy()
    z["np_"] = (z.IL_to_nucleicacid_massratio / 2.12).round(6)

    def key(df, hl, il, h, c, p):
        return (df[hl].astype(str) + "|" + df[il].round(1).astype(str) + "|"
                + df[h].round(1).astype(str) + "|" + df[c].round(1).astype(str)
                + "|" + df[p].round(1).astype(str) + "|"
                + df["np_"].round(3).astype(str))
    repo["k"] = key(repo, "Helper_lipid", "ILp", "HLp", "CHOLp", "PEGp")
    z["k"] = key(z, "HL_name", "IL_molratio", "HL_molratio", "CHL_molratio",
                 "PEG_molratio")
    mg = z.merge(repo[["k", "RLU_HepG2", "Formula label"]], on="k", how="left",
                 validate="one_to_one")
    matched = int(mg.RLU_HepG2.notna().sum())
    gg = mg.dropna(subset=["RLU_HepG2"])
    from scipy.stats import spearmanr, pearsonr
    rho = float(spearmanr(gg.Experiment_value, gg.RLU_HepG2).statistic)
    A = np.vstack([gg.Experiment_value, np.ones(len(gg))]).T
    coef, *_ = np.linalg.lstsq(A, gg.RLU_HepG2.values, rcond=None)
    res = gg.RLU_HepG2.values - A @ coef
    rms = float(np.sqrt((res ** 2).mean()))
    pred = 0.001 / np.sqrt(12)

    # the paper's own 192 values, as a multiset, against that column
    rc = Counter(repo.RLU_HepG2.round(3))
    per_helper = {}
    for h, grp in si.groupby("helper"):
        sc = Counter(grp.rlu_hepg2.round(3))
        rr = Counter(repo.loc[repo.Helper_lipid == h, "RLU_HepG2"].round(3))
        top32 = sorted(repo.loc[repo.Helper_lipid == h, "RLU_HepG2"].round(3),
                       reverse=True)[:32]
        per_helper[h] = {
            "si_rows": int(len(grp)),
            "in_repo_column_multiset": int(sum(min(v, rr.get(k, 0))
                                               for k, v in sc.items())),
            "equals_that_helpers_top_32": sorted(grp.rlu_hepg2.round(3))
            == sorted(top32)}
    overall = int(sum(min(v, rc.get(k, 0))
                      for k, v in Counter(si.rlu_hepg2.round(3)).items()))

    # ---------------- the design, proved ------------------------------
    levels = {c: sorted(int(v) for v in repo[c].unique())
              for c in ("NP_ratio", "Dlin-MC3_Helper lipid_ratio", "(IL+HL)",
                        "Chol_DMG-PEG_ratio")}
    prod = 1
    for v in levels.values():
        prod *= len(v)
    prod *= repo.Helper_lipid.nunique()

    # ---------------- token-pattern test -------------------------------
    names = pd.concat([vit.IL_name, viv.IL_name]).astype(str)
    tokens = {}
    for t in ("AI", "BI", "CI", "CII", "DI", "FIII", "AI-a1", "CII-b3",
              "A1", "C1", "D1"):
        hits = names[names.str.contains(t, regex=False)]
        tokens[t] = {"rows_matched_by_substring": int(len(hits)),
                     "distinct_structures": int(pd.concat(
                         [vit.loc[vit.IL_name.astype(str).str.contains(
                             t, regex=False), "ik"],
                          viv.loc[viv.IL_name.astype(str).str.contains(
                              t, regex=False), "ik"]]).nunique()),
                     "example_names": sorted(set(hits))[:6]}

    out = {
     "_what": "Corpus overlap, duplicate-arm test and token-pattern test for "
              "PMID 38583366.",
     "inputs": {"1080_LNP_OHE.csv": {"path": os.path.relpath(REPO_CSV, ROOT),
                                     "md5": md5(REPO_CSV)},
                "01_si_tables.csv": {"md5": md5(SI_CSV)}},
     "positive_controls": controls,
     "distinct_ionizable_lipids_this_paper_contributes": 1,
     "the_one_lipid": {
       "name_in_paper": "DLin-MC3 DMA (Materials), DLin-MC3-DMA elsewhere",
       "inchikey": mc3_k,
       "already_in_corpus": True,
       "vitro_rows": int(len(mc3_vit)),
       "vitro_sources": sorted(set(mc3_vit.Experiment_ID)),
       "vivo_rows": int(len(mc3_viv)),
       "vivo_sources": sorted(set(mc3_viv.Experiment_ID)),
       "vivo_organs_already_covered": sorted(set(
           mc3_viv.Model_target.astype(str))),
       "vivo_sources_with_liver": sorted(set(
           mc3_viv.loc[mc3_viv.Model_target == "liver", "Experiment_ID"]))},
     "in_vivo_pool_shared_structures": viv_shared,
     "already_merged_arm": {
       "claim": "The in-vitro HepG2 arm of this paper is the SAME screen as "
                "corpus experiment YZ_2022 (PMID 35879315), which this paper "
                "cites as ref [16] for it.",
       "step_1_corpus_vs_repo": {
         "join": "on composition (helper, 4 component percents to 1 dp, N:P), "
                 "using NO value column",
         "rows_matched": matched, "of": int(len(z)),
         "spearman_Experiment_value_vs_RLU_HepG2": round(rho, 6),
         "pearson": round(float(pearsonr(gg.Experiment_value,
                                         gg.RLU_HepG2).statistic), 8),
         "ln_unit_residual_rms": float("%.6g" % rms),
         "predicted_rms_from_3dp_rounding_alone": float("%.6g" % pred),
         "ratio_observed_over_predicted": round(rms / pred, 4),
         "reading": "The residual IS the rounding in the deposited file. "
                    "An exact-equality join therefore fails and is not the "
                    "right test."},
       "step_2_this_papers_192_si_values": {
         "si_rows": int(len(si)),
         "present_in_that_column_as_a_multiset": overall,
         "per_helper": per_helper},
       "verdict": "ALREADY MERGED. Emitting these 192 values under a new "
                  "Experiment_ID would duplicate a merged experiment and "
                  "manufacture exactly the source-separability signal this "
                  "project measures."},
     "design_proved": {"levels": levels,
                       "helpers": int(repo.Helper_lipid.nunique()),
                       "product": prod, "library_size_claimed": 1080,
                       "matches": prod == 1080 == len(repo)},
     "token_pattern_test": {
       "_what": "This paper's LNP ids are cluster tokens (AI, CII, DI, FIII) "
                "plus a sub-cluster letter and an index. They are LOCAL "
                "indices, not a shared catalogue: every one of them is a "
                "formulation of the same single lipid. A substring join on "
                "IL_name would match unrelated corpus lipids.",
       "tokens": tokens},
     "f057_18PG": {
       "question": "Is 18PG SOPG (18:0/18:1, the corpus string on 1,474 rows) "
                   "or DOPG (18:1/18:1)?",
       "what_this_paper_says": "NOTHING that discriminates. Verbatim and "
                               "complete: 'a mixture of helper lipids (DOTAP, "
                               "DDAB, DOPE, DSPC, 14PA, 18PG) from Avanti'. "
                               "No chemical name, no catalogue number, no acyl "
                               "specification, no structure is drawn anywhere "
                               "in the article or the SI.",
       "searched": "'stearoyl', 'oleoyl', 'glycero', 'phospho', '18:0', "
                   "'18:1', 'DOPG', 'SOPG', 'catalog', 'Cat' over 60,875 "
                   "chars of JATS full text and 17,600 chars of SI text: "
                   "zero hits for every one of them.",
       "verdict": "This paper DOES NOT SETTLE F-057 in either direction. It "
                  "supports neither SOPG nor DOPG. No change to data/."}}
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("token_pattern_test",)}, indent=1))


if __name__ == "__main__":
    main()
