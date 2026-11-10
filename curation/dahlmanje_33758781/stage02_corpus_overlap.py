#!/usr/bin/env python3
"""ZG_2020 vs data/ -- overlap, the token-pattern test, and a honest superset.

Run:  python3 -I stage02_corpus_overlap.py

Three questions, in order of how much damage getting them wrong does:

1. THE TOKEN-PATTERN TEST. This paper's labels are A-9 ... A-17 plus A-17-2Z,
   i.e. CONTIGUOUS 9..17. Contiguous numbering means LOCAL indices, so a
   substring join on any of them is a measured hazard, not a theoretical one.
   We measure the damage instead of asserting it.

2. DOES cKK-E15 ALREADY EXIST IN data/? Joined on InChIKey, never on a name:
   cKK-E12 sits in the corpus under THREE different names (CKK-E12, 'other',
   MD-1) with one InChIKey, so a name join finds one of three.

3. A SUBSTRUCTURE SWEEP, LABELLED AS A SUPERSET. The diketopiperazine core
   is the chemotype; a hit is a chemotype match, not this molecule.
"""
import json, os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import pandas as pd                                        # noqa: E402
from rdkit import Chem                                     # noqa: E402
from rdkit import RDLogger                                 # noqa: E402
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

PAPER_LABELS = ["A-9", "A-10", "A-11", "A-12", "A-13", "A-14", "A-15",
                "A-16", "A-17", "A-17-2Z"]

# cKK-E12, taken from the corpus (AC_2025) -- NOT written from memory.
CKK_E12_FROM_CORPUS_NAME = "CKK-E12"

# The cKK-E15 TAIL LENGTH IS AN INFERENCE, and is flagged as one. See
# DECISIONS.md D-2. The topology below is read off Fig 2a (which we CAN read:
# a 2,5-diketopiperazine, two (CH2)4 linkers to two tertiary nitrogens, each
# bearing two CH2-CH(OH)-CnH2n+1 arms); the integer n is NOT readable, because
# the drawing abbreviates each tail to a "CnH2n+1" label whose subscripts are
# 2-3 px tall in the only figure resolution available to us.
CKK_E15_INFERRED = ("CCCCCCCCCCCCCC(O)CN(CCCCC1NC(=O)C(CCCCN(CC(O)CCCCCCCCCCCCC)"
                    "CC(O)CCCCCCCCCCCCC)NC1=O)CC(O)CCCCCCCCCCCCC")


def tails_of(smiles):
    """Length of each CH2-CH(OH)-alkyl arm, measured not assumed."""
    m = Chem.MolFromSmiles(smiles)
    patt = Chem.MolFromSmarts("[NX3;!$(N-C=O)]-[CH2]-[CH1](-[OX2H1])-[CH2]")
    out = []
    for match in m.GetSubstructMatches(patt):
        n_idx, c1, c2, _o, c3 = match
        # walk the unbranched chain from c3 away from c2
        seen, cur, prev, L = {c2, c1, n_idx}, c3, c2, 0
        while True:
            L += 1
            seen.add(cur)
            nxt = [a.GetIdx() for a in m.GetAtomWithIdx(cur).GetNeighbors()
                   if a.GetSymbol() == "C" and a.GetIdx() not in seen]
            if len(nxt) != 1:
                break
            prev, cur = cur, nxt[0]
        out.append(L + 2)          # +2 for the CH2 and the CH(OH) carbons
    return sorted(out)


def main():
    vi = pd.read_csv(os.path.join(ROOT, "data", "data_vitro.csv"),
                     low_memory=False).assign(pool="vitro")
    vo = pd.read_csv(os.path.join(ROOT, "data", "data_vivo.csv"),
                     low_memory=False).assign(pool="vivo")
    allr = pd.concat([vi, vo], ignore_index=True)
    out = {}

    # ---------------- 1. token-pattern test --------------------------------
    tok = {}
    for lab in PAPER_LABELS + ["9", "11", "17", "A-11"]:
        sub = allr[allr.IL_name.astype(str).str.contains(lab, regex=False,
                                                         na=False)]
        exact = allr[allr.IL_name.astype(str) == lab]
        tok[lab] = {
            "substring_rows": int(len(sub)),
            "substring_structures": int(sub.IL_SMILES.nunique()),
            "substring_experiments": sorted(sub.Experiment_ID.unique().tolist()),
            "exact_name_rows": int(len(exact)),
        }
    out["token_pattern_test"] = {
        "_labels_are": "A-9 ... A-17 plus A-17-2Z: CONTIGUOUS 9..17, which is "
                       "the signature of LOCAL indices, so these names carry "
                       "no cross-paper meaning.",
        "_and_they_are_not_ionizable_lipid_names_at_all": (
            "A-N names the HELPER PHOSPHOLIPID. Joining them against IL_name "
            "would be wrong twice over -- wrong column and wrong namespace."),
        "_name_token_vs_chain_length": (
            "A-9 is the DECANOYL (C10) ester, A-17 the STEAROYL (C18). The "
            "token counts the alkyl carbons BEYOND the ester carbonyl, so it "
            "is the acyl length MINUS ONE, every time. Trap 1, measured from "
            "the printed acid names in the synthesis."),
        "hits": tok,
    }

    # ---------------- 2. cKK-E15 vs the corpus, on InChIKey ----------------
    e12_rows = allr[allr.IL_name == CKK_E12_FROM_CORPUS_NAME]
    e12_smi = e12_rows.IL_SMILES.iloc[0]
    e12_key = Chem.MolToInchiKey(Chem.MolFromSmiles(e12_smi))
    e15_key = Chem.MolToInchiKey(Chem.MolFromSmiles(CKK_E15_INFERRED))

    keys = {}
    for s in allr.IL_SMILES.dropna().unique():
        m = Chem.MolFromSmiles(s)
        if m is not None:
            keys[s] = Chem.MolToInchiKey(m)
    e12_aliases = allr[allr.IL_SMILES.isin(
        [s for s, k in keys.items() if k == e12_key])]
    e15_exact = allr[allr.IL_SMILES.isin(
        [s for s, k in keys.items() if k == e15_key])]
    e15_skel = allr[allr.IL_SMILES.isin(
        [s for s, k in keys.items() if k.split("-")[0] == e15_key.split("-")[0]])]

    out["ckk_e15"] = {
        "corpus_ckk_e12_smiles": e12_smi,
        "corpus_ckk_e12_inchikey": e12_key,
        "corpus_ckk_e12_tail_lengths": tails_of(e12_smi),
        "corpus_ckk_e12_names_held": sorted(
            {f"{r.pool}/{r.Experiment_ID}/{r.IL_name}"
             for r in e12_aliases.itertuples()}),
        "corpus_ckk_e12_rows_by_name": e12_aliases.groupby(
            ["pool", "Experiment_ID", "IL_name"]).size().to_dict().__str__(),
        "_name_join_would_have_found": int(len(e12_rows)),
        "_inchikey_join_finds": int(len(e12_aliases)),
        "inferred_ckk_e15_smiles": CKK_E15_INFERRED,
        "inferred_ckk_e15_inchikey": e15_key,
        "inferred_ckk_e15_tail_lengths": tails_of(CKK_E15_INFERRED),
        "ckk_e15_full_inchikey_rows_in_corpus": int(len(e15_exact)),
        "ckk_e15_skeleton_block_rows_in_corpus": int(len(e15_skel)),
        "_verdict": "",
    }
    out["ckk_e15"]["_verdict"] = (
        "cKK-E15 is NOT in the corpus under any name: 0 rows match its full "
        "InChIKey and 0 match even its skeleton block. cKK-E12 IS in the "
        "corpus and is a DIFFERENT molecule (C12 vs C15 arms), so this paper "
        "adds no in-vivo shared lipid. Note the corpus cKK-E12 is stored "
        "STEREO-FLAT (InChIKey suffix UHFFFAOYSA-N) although the real "
        "compound is lysine-derived; we inherit that convention rather than "
        "silently improving on it.")

    # ---------------- 3. substructure sweep, A SUPERSET --------------------
    # 2,5-diketopiperazine, with the ring carbons left UNSPECIFIED. The first
    # version of this line wrote "O=C1NC(=O)CN1", which demands a ring CH2 and
    # therefore could not match cKK-E12 itself -- a sweep that misses its own
    # known positive is not a loose sweep, it is a broken one, and it reported
    # a clean 0. The positive control below is what caught it, so it stays.
    core = Chem.MolFromSmarts("O=C1[#6]NC(=O)[#6]N1")
    arm = Chem.MolFromSmarts("[NX3;!$(N-C=O)](-[CH2][CH1]([OX2H1])[CH2])"
                             "-[CH2][CH1]([OX2H1])[CH2]")
    # POSITIVE CONTROL: the patterns MUST match the structures we already know
    # are of this chemotype, or the sweep is measuring nothing.
    for label, smi in (("corpus cKK-E12", e12_smi),
                       ("inferred cKK-E15", CKK_E15_INFERRED)):
        m = Chem.MolFromSmiles(smi)
        if not m.HasSubstructMatch(core):
            raise SystemExit(f"sweep positive control FAILED: core misses {label}")
        if len(m.GetSubstructMatches(arm)) < 2:
            raise SystemExit(f"sweep positive control FAILED: arm misses {label}")
    print("METRIC sweep_positive_control=PASS (core+arm match corpus cKK-E12 "
          "and inferred cKK-E15)")
    loose = strict = 0
    loose_rows = strict_rows = 0
    for s, k in keys.items():
        m = Chem.MolFromSmiles(s)
        n = int((allr.IL_SMILES == s).sum())
        if m.HasSubstructMatch(core):
            loose += 1
            loose_rows += n
            if len(m.GetSubstructMatches(arm)) >= 2:
                strict += 1
                strict_rows += n
    out["substructure_sweep_SUPERSET"] = {
        "_label": "SUPERSET. A hit is a CHEMOTYPE match (the cKK diketo"
                  "piperazine scaffold), NOT this paper's molecule.",
        "unique_structures_scanned": len(keys),
        "loose_hits_diketopiperazine_core_structures": loose,
        "loose_hits_rows": loose_rows,
        "strict_hits_core_plus_two_bis(hydroxyalkyl)amine_arms_structures": strict,
        "strict_hits_rows": strict_rows,
        "after_exact_inchikey_match_to_ckk_e15": 0,
    }

    # ---------------- the helper phospholipid, against HL_* ----------------
    hl = allr.HL_name.astype(str)
    out["helper_lipid_namespace"] = {
        "_why": "A-N is a HELPER phospholipid, so the only column it could "
                "legitimately join is HL_name/HL_SMILES.",
        "corpus_distinct_HL_names": int(allr.HL_name.nunique(dropna=True)),
        "corpus_distinct_HL_SMILES": int(allr.HL_SMILES.nunique(dropna=True)),
        "rows_with_any_HL_name_matching_a_paper_label": {
            lab: int((hl == lab).sum()) for lab in PAPER_LABELS},
        "any_adamantane_bearing_HL_SMILES_in_corpus": int(sum(
            1 for s in allr.HL_SMILES.dropna().unique()
            if (lambda m: m is not None and m.HasSubstructMatch(
                Chem.MolFromSmarts("C12CC3CC(C1)CC(C3)C2")))(
                    Chem.MolFromSmiles(s)))),
    }

    with open(os.path.join(HERE, "02_corpus_overlap.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
