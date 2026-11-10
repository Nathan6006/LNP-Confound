#!/usr/bin/env python3
"""Corpus overlap for the 15 Han et al. ionizable lipids, measured three ways
plus an honestly-labelled superset sweep.

Exact keys:  FULL InChIKey, SKELETON InChIKey (connectivity only), IL_name.
Superset:    a SMARTS for the paper's own arm motif, reported as a SUPERSET,
             because a sibling's first SMARTS returned 332 rows under a label
             claiming exactness when the truth was 128.

F-055 CROSS-LAB AMINE CATALOGUE. The amine numbers 200/306/114/490/494 are the
Anderson/Akinc lipidoid catalogue and are SHARED ACROSS LABS, so a head number
here can be pinned by another lab's rows. Matching is done on STRUCTURE, and the
name sweep below is reported separately and explicitly as a name sweep.

A LEADING TOKEN IS NOT AN IDENTITY. A sibling required a name to START with a
head number and so missed every row whose name began with a letter; and the same
digits occur in unrelated series indices. So the name sweep uses a DIGIT-BOUNDARY
search (not ^-anchored), and every hit is then adjudicated BY STRUCTURE, which is
how the false matches below were identified as false rather than guessed at.
"""
import csv, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
import pandas as pd, re
from rdkit import Chem
from rdkit.Chem import inchi

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
OURS = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))

ARM_SMARTS = "[NX3;!$(N=*)][CH2][CH1]([OX2H1])[CH2]"   # N-CH2-CH(OH)-CH2...
CORE_306   = "[CH3][NX3]([CH2][CH2][CH2][NX3])[CH2][CH2][CH2][NX3]"
CORE_PIPZ  = "[NX3]1[CH2][CH2][NX3][CH2][CH2]1"


def load(tag):
    d = pd.read_csv(os.path.join(ROOT, "data", f"data_{tag}.csv"), low_memory=False)
    d["_tag"] = tag
    return d


def keyify(s):
    m = Chem.MolFromSmiles(s) if isinstance(s, str) and s.strip() else None
    if m is None:
        return None, None
    k = inchi.MolToInchiKey(m)
    return k, k.split("-")[0]


def main():
    out = {"_what": __doc__.strip().splitlines()[0]}
    vit, viv = load("vitro"), load("vivo")
    allc = pd.concat([vit, viv], ignore_index=True)

    # index the corpus once, by structure
    uniq = allc[["IL_SMILES"]].dropna().drop_duplicates()
    full, skel = {}, {}
    bad = 0
    for s in uniq["IL_SMILES"]:
        k, sk = keyify(s)
        if k is None:
            bad += 1; continue
        full.setdefault(k, set()).add(s)
        skel.setdefault(sk, set()).add(s)
    out["corpus_index"] = dict(
        vitro_rows=int(len(vit)), vivo_rows=int(len(viv)),
        unique_IL_SMILES=int(len(uniq)), unparseable=bad,
        distinct_full_inchikeys=len(full), distinct_skeleton_inchikeys=len(skel))

    # ---------------------------------------------------- exact, by structure
    exact = []
    for r in OURS:
        hit_full = r["InChIKey"] in full
        hit_skel = r["InChIKey_skeleton"] in skel
        rec = dict(IL_name=r["IL_name"], InChIKey=r["InChIKey"],
                   full_inchikey_match=bool(hit_full),
                   skeleton_inchikey_match=bool(hit_skel),
                   corpus_rows_full=0, corpus_rows_skeleton=0,
                   corpus_experiments=[], corpus_names=[])
        if hit_skel:
            sm = skel[r["InChIKey_skeleton"]]
            m = allc[allc["IL_SMILES"].isin(sm)]
            rec["corpus_rows_skeleton"] = int(len(m))
            rec["corpus_experiments"] = sorted(set(m["Experiment_ID"]))
            rec["corpus_names"] = sorted(set(m["IL_name"].astype(str)))
            rec["corpus_pools"] = sorted(set(m["_tag"]))
        if hit_full:
            m = allc[allc["IL_SMILES"].isin(full[r["InChIKey"]])]
            rec["corpus_rows_full"] = int(len(m))
        exact.append(rec)
    out["exact_structure_match"] = exact
    out["exact_structure_match_summary"] = dict(
        n_ours=len(OURS),
        n_full_inchikey_match=sum(1 for e in exact if e["full_inchikey_match"]),
        n_skeleton_only_match=sum(1 for e in exact
                                  if e["skeleton_inchikey_match"] and not e["full_inchikey_match"]),
        n_no_match=sum(1 for e in exact if not e["skeleton_inchikey_match"]))

    # ------------------------------------------------------------ IL_name sweep
    name_hits = []
    for tok in ["200", "306", "114", "490", "494"]:
        pat = re.compile(r"(?<![0-9])" + tok + r"(?![0-9])")
        m = allc[allc["IL_name"].astype(str).str.contains(pat, regex=True, na=False)]
        for (nm, ex), g in m.groupby(["IL_name", "Experiment_ID"]):
            smis = sorted(set(g["IL_SMILES"].dropna()))
            # adjudicate BY STRUCTURE, not by the name
            carries_arm = carries_306 = carries_pipz = False
            for s in smis:
                mol = Chem.MolFromSmiles(s)
                if mol is None: continue
                if mol.HasSubstructMatch(Chem.MolFromSmarts(ARM_SMARTS)): carries_arm = True
                if mol.HasSubstructMatch(Chem.MolFromSmarts(CORE_306)):   carries_306 = True
                if mol.HasSubstructMatch(Chem.MolFromSmarts(CORE_PIPZ)):  carries_pipz = True
            name_hits.append(dict(
                amine_token=tok, corpus_IL_name=nm, Experiment_ID=ex,
                pool=sorted(set(g["_tag"])), rows=int(len(g)),
                n_distinct_SMILES=len(smis),
                carries_2hydroxyalkyl_arm=carries_arm,
                carries_amine306_core=carries_306,
                carries_piperazine=carries_pipz,
                SMILES_first=smis[0][:140] if smis else ""))
    out["IL_name_token_sweep"] = name_hits
    out["IL_name_token_sweep_note"] = (
        "A NAME SWEEP, NOT AN IDENTITY CLAIM. Each hit is adjudicated by the "
        "structure columns beside it. Search is digit-boundary, not ^-anchored, "
        "so names such as 'A200_C12' and '114-L1-C12' are both reached.")

    # ------------------------------------------- superset substructure sweep
    arm = Chem.MolFromSmarts(ARM_SMARTS)
    rows_arm = []
    for s in uniq["IL_SMILES"]:
        mol = Chem.MolFromSmiles(s)
        if mol is not None and mol.HasSubstructMatch(arm):
            rows_arm.append(s)
    marm = allc[allc["IL_SMILES"].isin(rows_arm)]
    out["superset_sweep_2hydroxyalkyl_amine"] = dict(
        SMARTS=ARM_SMARTS,
        LABEL="SUPERSET - NOT an identity match",
        why_superset=("This SMARTS matches ANY tertiary-amine 2-hydroxyalkyl motif, "
                      "including lipidoids from unrelated libraries and any lipid that "
                      "merely happens to carry an N-CH2-CH(OH)-CH2 fragment. It is "
                      "reported as an upper bound on how much of the corpus is "
                      "CHEMICALLY NEARBY, and must not be read as overlap."),
        unique_structures=len(rows_arm), rows=int(len(marm)),
        experiments=sorted(set(marm["Experiment_ID"])),
        rows_vitro=int((marm["_tag"] == "vitro").sum()),
        rows_vivo=int((marm["_tag"] == "vivo").sum()))


    # --------------------------------------- amine-core sweep (NEAR MISSES)
    # Does the corpus contain our AMINE CORE under a DIFFERENT tail chemistry?
    # That is the chemically interesting near miss: same head, different linker,
    # which is what makes a head number cross-lab-pinnable at all (F-055).
    CORE_SMARTS = {
        "200": "[NX3]([CH2][CH2][NX3]1[CH2][CH2][NX3]([CH2][CH2]1)[CH2][CH2][NX3][CH2][CH2][NX3])",
        "306": CORE_306,
        "114": "[NX3][CH2][CH2][NX3][CH2][CH2][NX3]",
        "494": "[NX3][CH2][CH2][OX2][CH2][CH2][NX3]1[CH2][CH2][NX3][CH2][CH2]1",
        "490": "[NX3][CH2][CH1]([OX2][CH2][CH3])[CH2][NX3]1[CH2][CH2][NX3][CH2][CH2]1",
    }
    core_sweep = {}
    ours_keys = {r["InChIKey"] for r in OURS}
    for amine, sma in CORE_SMARTS.items():
        q = Chem.MolFromSmarts(sma)
        if q is None:
            core_sweep[amine] = {"ERROR": "SMARTS did not parse"}; continue
        hits = [s for s in uniq["IL_SMILES"]
                if (Chem.MolFromSmiles(s) is not None
                    and Chem.MolFromSmiles(s).HasSubstructMatch(q))]
        mh = allc[allc["IL_SMILES"].isin(hits)]
        # exclude structures that ARE ours, so this reports only near misses
        not_ours = []
        for s in hits:
            mm = Chem.MolFromSmiles(s)
            if inchi.MolToInchiKey(mm) not in ours_keys:
                not_ours.append(s)
        mn = allc[allc["IL_SMILES"].isin(not_ours)]
        core_sweep[amine] = dict(
            SMARTS=sma, LABEL="SUPERSET of same-core corpus chemistry",
            unique_structures=len(hits), rows=int(len(mh)),
            experiments=sorted(set(mh["Experiment_ID"])),
            near_miss_unique_structures_excluding_ours=len(not_ours),
            near_miss_rows=int(len(mn)),
            near_miss_experiments=sorted(set(mn["Experiment_ID"])),
            near_miss_names_sample=sorted(set(mn["IL_name"].astype(str)))[:16])
    out["amine_core_near_miss_sweep"] = core_sweep

    # -------------------------- the KL_2010 / AC_2025 name-vs-structure split
    # MEASURED, not inferred from a grep: KL_2010 is Love et al. 2010, the SAME
    # epoxide-amine lipidoid chemistry, and its 114-series matches ours exactly.
    # Its 200-series does not. Both facts are computed here.
    split = []
    for nm in ["200C8", "200C10", "200C12", "200C14", "114C12", "114C14", "114C16"]:
        m = allc[allc["IL_name"].astype(str) == nm]
        if not len(m):
            continue
        for s in sorted(set(m["IL_SMILES"].dropna())):
            mm = Chem.MolFromSmiles(s)
            k = inchi.MolToInchiKey(mm)
            ourname = [r["IL_name"] for r in OURS if r["InChIKey"] == k]
            from rdkit.Chem.rdMolDescriptors import CalcMolFormula as CF
            from rdkit.Chem.Descriptors import MolWt as MW
            split.append(dict(
                corpus_IL_name=nm,
                Experiment_ID=sorted(set(m["Experiment_ID"]))[0],
                rows=int(len(m)), formula=CF(mm), MW=round(MW(mm), 2),
                n_N=sum(1 for a in mm.GetAtoms() if a.GetSymbol() == "N"),
                n_2hydroxyalkyl_arms=len(mm.GetSubstructMatches(
                    Chem.MolFromSmarts("[NX3][CH2][CH1]([OX2H1])[CH2]"))),
                matches_our_lipid=ourname or None))
    # and the structure we derive for the same name
    for r in OURS:
        if r["amine_core"] in ("200", "114") and r["tail"] in ("C10", "C12", "C14"):
            split.append(dict(
                corpus_IL_name=f"OURS:{r['IL_name']}", Experiment_ID="EH_2026(proposed)",
                rows=None, formula=r["formula"], MW=float(r["MW"]),
                n_N=int(r["formula"].split("N")[1].split("O")[0] or 1),
                n_2hydroxyalkyl_arms=int(r["n_arms_built"]),
                matches_our_lipid=[r["IL_name"]]))
    out["name_vs_structure_split_200_vs_114"] = split
    out["name_vs_structure_split_note"] = (
        "KL_2010 is Love et al. PNAS 2010, the ORIGINAL epoxide-amine lipidoid "
        "library that this paper's amine numbering comes from. Its 114C12 and "
        "114C14 match our C12-114 and C14-114 on the FULL InChIKey, which "
        "independently confirms the cross-lab amine catalogue (F-055). Its "
        "200C8/200C10/200C12/200C14 do NOT match our C*-200 and are not even the "
        "same size class. This is recorded as a CONTRADICTION IN data/, not as a "
        "doubt about our structure: AC_2025 (498 rows) and JW_2024 both carry our "
        "exact C12-200 InChIKey, and our MW is reproduced by the paper's own "
        "stoichiometry. We do not edit data/. See 08_contradictions.csv.")

    json.dump(out, open(os.path.join(HERE, "02_corpus_overlap.json"), "w"), indent=1)
    s = out["exact_structure_match_summary"]
    print("EXACT:", s)
    for e in exact:
        if e["skeleton_inchikey_match"]:
            print(f"  HIT {e['IL_name']}: full={e['full_inchikey_match']} "
                  f"rows={e['corpus_rows_skeleton']} exp={e['corpus_experiments']} "
                  f"names={e['corpus_names']} pools={e.get('corpus_pools')}")
    print("SUPERSET:", {k: v for k, v in out["superset_sweep_2hydroxyalkyl_amine"].items()
                        if k in ("unique_structures", "rows", "rows_vitro", "rows_vivo")})
    print("name sweep rows:", len(name_hits))


if __name__ == "__main__":
    main()
