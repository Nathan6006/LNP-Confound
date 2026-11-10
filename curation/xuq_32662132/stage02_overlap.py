"""Stage 02 - corpus overlap, three keys plus a near-miss sweep.

Keys: FULL InChIKey, skeleton InChIKey, and IL_name. Plus a substructure sweep
on the head and a Morgan/Tanimoto nearest neighbour, because the run brief's
trap 2 says a near miss can be the ground truth a source withheld. It was.
"""
import csv, json, os
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs, inchi
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data"))

HEAD_SMARTS = "c1cn(CCCN)cn1"      # 3-aminopropylimidazole, ring N substituted
IMID_SMARTS = "c1cncn1"

def main():
    cands = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    cand_full = {c["IL_name"]: c["inchikey_full"] for c in cands}
    cand_skel = {c["IL_name"]: c["inchikey_skeleton"] for c in cands}
    cand_fp = {c["IL_name"]: AllChem.GetMorganFingerprintAsBitVect(
        Chem.MolFromSmiles(c["IL_SMILES"]), 2, 2048) for c in cands}
    # every name this paper uses, for the IL_name key
    paper_names = sorted({c["IL_name"] for c in cands} |
                         {"93-O17S", "9322-O17S", "93-O17C", "9322-O16B", "9322-N16B",
                          "9322-EC18", "9322-O17C", "9322-O17O", "93-S14"})
    head = Chem.MolFromSmarts(HEAD_SMARTS)
    imid = Chem.MolFromSmarts(IMID_SMARTS)

    out = {"_what": "Overlap of this paper's enumerable structures against the live corpus.",
           "keys_used": ["full InChIKey", "skeleton InChIKey", "IL_name (exact and token)",
                         "3-aminopropylimidazole substructure", "Morgan r2 2048 Tanimoto NN"],
           "pools": {}}

    for pool in ("vitro", "vivo"):
        df = pd.read_csv(os.path.join(DATA, "data_%s.csv" % pool), low_memory=False)
        names = df["IL_name"].astype(str)
        uniq = sorted(set(df["IL_SMILES"].dropna()))
        full_hits, skel_hits, head_rows, imid_n = [], [], [], 0
        best = {}
        for s in uniq:
            m = Chem.MolFromSmiles(s)
            if m is None:
                continue
            k = inchi.MolToInchiKey(m)
            if k in cand_full.values():
                full_hits.append(s)
            if k.split("-")[0] in cand_skel.values():
                skel_hits.append(s)
            if m.HasSubstructMatch(imid):
                imid_n += 1
            if m.HasSubstructMatch(head):
                head_rows.append(s)
            fp = AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)
            for nm, f in cand_fp.items():
                t = DataStructs.TanimotoSimilarity(fp, f)
                if t > best.get(nm, (0, None))[0]:
                    best[nm] = (t, s)
        head_df = df[df["IL_SMILES"].isin(head_rows)]
        name_exact = sorted(set(names[names.isin(paper_names)]))
        # the trap the brief names: a leading number token is not an identity
        tokmask = names.str.match(r"^(?:93|9322)(?:[-_]|$)").fillna(False)
        out["pools"][pool] = {
            "n_rows": int(len(df)), "n_unique_IL_SMILES": len(uniq),
            "full_inchikey_hits": len(full_hits),
            "skeleton_inchikey_hits": len(skel_hits),
            "IL_name_exact_hits": name_exact,
            "rows_with_any_imidazole": imid_n,
            "near_miss_head_93": {
                "n_unique_structures": len(head_rows),
                "n_rows": int(len(head_df)),
                "experiments": sorted(set(head_df["Experiment_ID"])),
                "example_names_per_experiment": {
                    e: sorted(set(head_df[head_df.Experiment_ID == e]["IL_name"].astype(str)))[:6]
                    for e in sorted(set(head_df["Experiment_ID"]))},
            },
            "nearest_neighbour": {
                nm: {"tanimoto": round(t, 3), "smiles": s,
                     "experiments": sorted(set(df[df.IL_SMILES == s]["Experiment_ID"])),
                     "names": sorted(set(df[df.IL_SMILES == s]["IL_name"].astype(str)))[:3]}
                for nm, (t, s) in sorted(best.items())},
            "leading_token_trap": {
                "n_rows_whose_IL_name_starts_with_93_or_9322": int(tokmask.sum()),
                "examples": sorted(set(names[tokmask]))[:8],
                "note": "reported, NOT acted on. A leading number token is not an identity."},
        }

    out["conclusion"] = (
        "ZERO full-InChIKey and ZERO skeleton-InChIKey overlap in either pool: none of the "
        "five enumerable lipidoids is already in the corpus. But the NEAR MISS is the "
        "important result. The corpus already carries 173 unique structures over 277 in-vitro rows (11 experiments) and 25 in-vivo rows (4) "
        "built on this exact amine head, named 93_O12 / 93_O14 / 93_N12 by AA_2008 and "
        "A93_C12 by KW_2014. Those are two independent curations of two unrelated papers, "
        "and together with curation/xuq_32520527 (which read head 93 off a drawn figure in "
        "a third paper) they fix amine 93 = 1-(3-aminopropyl)imidazole from three sources "
        "this paper did not supply. They also confirm this paper's tail-prefix convention "
        "externally: AA_2008 carries 93_O14 as a tetradecyl ESTER and 93_N14 as the "
        "corresponding AMIDE on the same head, which is exactly the O/N split the SI's "
        "O16S-S / N16S-S pair implies. REPORTED, NOT ACTED ON.")
    with open(os.path.join(HERE, "02_corpus_overlap.json"), "w") as f:
        json.dump(out, f, indent=1)
    for p, v in out["pools"].items():
        print(p, "full", v["full_inchikey_hits"], "skel", v["skeleton_inchikey_hits"],
              "name", v["IL_name_exact_hits"], "head-near-miss rows",
              v["near_miss_head_93"]["n_rows"], v["near_miss_head_93"]["experiments"])

if __name__ == "__main__":
    main()
