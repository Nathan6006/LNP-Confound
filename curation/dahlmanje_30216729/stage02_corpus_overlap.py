#!/usr/bin/env python3
"""Stage 02 -- what, if anything, this paper's chemistry already has in the corpus.

    python3 -I stage02_corpus_overlap.py [OUTDIR]

This paper contributes ONE ionizable component, `7C1`, and 7C1 is a
polydisperse polymer conjugate (branched PEI600 opened with a C15 epoxide,
drawn as a generic R-group scheme on SI p14), so there is no single SMILES to
join on. The overlap question therefore splits in two:

  (a) the ionizable component -- answered structurally: it CANNOT be written,
      so no InChIKey join is possible, and the token joins below show what a
      name join would have done instead;
  (b) the excipients the paper DOES draw in full -- cholesterol and DSPC --
      which can be joined on InChIKey off provenance/reagent_registry.py.

EVERY QUERY CARRIES A POSITIVE CONTROL (F-070). A sweep that returns 0 because
its pattern cannot match anything is indistinguishable, in the report, from a
sweep that returns 0 because nothing is there.
"""
import json
import os
import sys

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
OUTDIR = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else HERE
sys.path.insert(0, os.path.join(REPO, "provenance"))
import reagent_registry as rr   # noqa: E402

NAME_COLS = ["IL_name", "HL_name", "CHL_name", "PEG_name",
             "fifthcomponent_name"]
SMI_COLS = ["IL_SMILES", "HL_SMILES", "CHL_SMILES", "PEG_SMILES",
            "fifthcomponent_SMILES"]


def ikey(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return Chem.MolToInchiKey(m)


def load_corpus():
    out = {}
    for pool, fn in (("vitro", "data_vitro.csv"), ("vivo", "data_vivo.csv")):
        out[pool] = pd.read_csv(os.path.join(REPO, "data", fn), low_memory=False)
    return out


# -------------------------------------------------- token-pattern test ------

def token_test(corpus):
    """MANDATORY. Measure what a naive substring join on this paper's name
    tokens would have returned. `7C1` packs three characters, none of which is
    a chain length; `C14`/`C18` here are PEG ALKYL TAILS, a different
    namespace from the ionizable lipid entirely.
    """
    tokens = ["7C1", "7", "C1", "C14", "C18", "CHOL", "Chol", "DSPC"]
    res = {}
    for tok in tokens:
        rows, structs, where = 0, set(), {}
        for pool, df in corpus.items():
            for nc, sc in zip(NAME_COLS, SMI_COLS):
                m = df[nc].astype(str).str.contains(tok, case=True, regex=False,
                                                    na=False)
                n = int(m.sum())
                if n:
                    rows += n
                    where[f"{pool}:{nc}"] = n
                    structs |= set(df.loc[m, sc].dropna().astype(str))
        res[tok] = {"substring_rows": rows, "distinct_structures": len(structs),
                    "columns_hit": where,
                    "example_names": sorted(
                        {str(x) for pool, df in corpus.items()
                         for nc in NAME_COLS
                         for x in df.loc[df[nc].astype(str).str.contains(
                             tok, case=True, regex=False, na=False), nc]
                         .astype(str).unique()[:6]})[:8]}
    # exact-match join, the only name join that would be defensible
    exact = {}
    for tok in ["7C1", "7c1", "7C1 ", "PEI600", "PEI", "C15 epoxide"]:
        rows = 0
        for pool, df in corpus.items():
            for nc in NAME_COLS:
                rows += int((df[nc].astype(str).str.strip() == tok).sum())
        exact[tok] = rows
    res["_exact_name_match"] = exact
    res["_positive_control"] = {
        "token": "DSPC",
        "expectation": "DSPC is a corpus helper lipid and MUST return many rows",
        "substring_rows": res["DSPC"]["substring_rows"],
        "passes": res["DSPC"]["substring_rows"] > 100,
        "why_this_control": ("if the name-column scan were reading nothing -- a "
                            "wrong column list, an all-NaN cast -- every token "
                            "would return 0 and the 7C1 zero would be vacuous"),
    }
    res["_interpretation"] = (
        "A substring join on this paper's own tokens is catastrophic and the "
        "numbers below quantify it. `7C1` matches nothing by any route, so the "
        "ionizable component is genuinely absent -- but `C14` and `C18`, which "
        "ARE tokens printed in this paper's formulation tables, match corpus "
        "rows in bulk while naming PEG ALKYL TAILS in this paper and something "
        "else there. That is the wrong-namespace failure: the token is real, "
        "the match is real, and the join is still meaningless.")
    return res


# ------------------------------------------------ structural attestation ----

def structural_join(corpus):
    """InChIKey join for the components whose structures are fixed and public.

    Cholesterol and DSPC are the two excipients SI p14 draws in full. 7C1 is
    absent from the registry on purpose and must stay absent: like a PEG-lipid,
    its block is a DISTRIBUTION, so no single SMILES is right.
    """
    res = {"queried": {}, "positive_control": {}}
    targets = {"Cholesterol": "drawn in full on SI p14",
               "DSPC": "drawn in full on SI p14",
               "DLin-MC3-DMA": "POSITIVE CONTROL -- not in this paper"}
    # Build the corpus InChIKey index once.
    index = {}
    for pool, df in corpus.items():
        for sc, nc in zip(SMI_COLS, NAME_COLS):
            sub = df[[sc, nc]].dropna(subset=[sc])
            for smi, nm in zip(sub[sc].astype(str), sub[nc].astype(str)):
                k = index.setdefault(smi, None)
                if k is None:
                    index[smi] = ikey(smi)
    key_to_rows, key_to_names = {}, {}
    for pool, df in corpus.items():
        for sc, nc in zip(SMI_COLS, NAME_COLS):
            sub = df[[sc, nc]].dropna(subset=[sc])
            for smi, nm in zip(sub[sc].astype(str), sub[nc].astype(str)):
                k = index.get(smi)
                if not k:
                    continue
                key_to_rows[k] = key_to_rows.get(k, 0) + 1
                key_to_names.setdefault(k, set()).add(nm)

    for name, why in targets.items():
        entry_raw = rr.REGISTRY.get(name)
        if entry_raw is None:
            res["queried"][name] = {"in_registry": False, "why": why}
            continue
        # REGISTRY values are (smiles, chemical_name); take the SMILES and
        # carry the registry's own chemical name through so the report states
        # WHICH molecule was pinned, not just that something was.
        smi, chem_name = entry_raw[0], entry_raw[1]
        k = ikey(smi)
        entry = {
            "in_registry": True, "why_queried": why,
            "registry_smiles": smi, "registry_chemical_name": chem_name,
            "inchikey": k,
            "corpus_rows_by_inchikey": key_to_rows.get(k, 0),
            "names_the_structure_sits_under": sorted(key_to_names.get(k, [])),
        }
        if "POSITIVE CONTROL" in why:
            entry["passes"] = entry["corpus_rows_by_inchikey"] > 1000
            entry["expectation"] = ("DLin-MC3-DMA is documented at ~3,087 corpus "
                                    "rows; a join returning 0 would mean the "
                                    "InChIKey index is broken and every other "
                                    "zero here is uninterpretable")
            res["positive_control"] = entry
        else:
            res["queried"][name] = entry

    res["ionizable_component_7C1"] = {
        "in_registry": "7C1" in rr.REGISTRY,
        "can_a_single_smiles_be_written": False,
        "why_not": ("SI p14 draws 7C1 as a GENERIC SCHEME: a C15 epoxide (drawn "
                    "as R-epoxide with R=C13) reacted with PEI600 in ethanol at "
                    "90 C for 72 h. PEI600 is branched, polydisperse "
                    "polyethylenimine of nominal MW 600, so the product is a "
                    "distribution over both backbone topology and degree of "
                    "substitution. No atom-resolved structure is drawn, stated, "
                    "or derivable, and the paper gives no formula, no HRMS and "
                    "no NMR for it."),
        "precedent": ("CLAUDE.md section 3: PEG-lipids are deliberately absent "
                      "from the reagent registry because the block length is a "
                      "distribution, so no single SMILES is right. 7C1 is the "
                      "same case on the ionizable component."),
        "must_not_be_guessed": ("Writing a representative monomer would put a "
                                "made-up structure in IL_SMILES. The brief's "
                                "rule is to flag, never guess -- and with zero "
                                "emittable rows there is nothing to flag it on."),
    }
    return res


# ----------------------------------------------------- substructure sweep ---

def substructure_sweep(corpus):
    """SUPERSET sweep: does the corpus already contain epoxide-derived
    amino-alcohol lipidoids of 7C1's chemotype? Reported as a SUPERSET, never
    as an overlap -- a chemotype match is not the same molecule.
    """
    queries = {
        "aminoalcohol_from_epoxide_ring_opening":
            ("[NX3;!$(N-C=O)]-[CH2]-[CH;!$(C=O)]-[OX2H1]",
             "the motif 7C1 would carry: a ring-opened 1,2-epoxide giving "
             "N-CH2-CH(OH)-alkyl"),
        "aminoalcohol_with_long_tail":
            ("[NX3;!$(N-C=O)]-[CH2]-[CH]([OX2H1])-[CH2][CH2][CH2][CH2][CH2]",
             "same motif, tightened to require >=5 further tail carbons"),
        "polyamine_backbone_3plus_N":
            ("[NX3]-[CH2]-[CH2]-[NX3]-[CH2]-[CH2]-[NX3]",
             "a PEI-like ethylenimine run of three nitrogens"),
    }
    # A POSITIVE CONTROL MUST ASSERT ITS OWN PREMISE.
    # First version of this function hand-wrote a steroid SMARTS, got 6 hits,
    # and reported `passes: True` -- while the pattern did NOT match the
    # registry's cholesterol at all (checked: HasSubstructMatch -> False). A
    # non-zero count made a mis-specified control look sound, which is the
    # F-070 failure one step further on than F-070 itself. So each control now
    # names the molecule it MUST match, the sweep verifies the match, and a
    # control whose own premise fails raises instead of reporting a number.
    controls = {
        "tertiary_aliphatic_amine": (
            "[NX3;H0;!$(N-C=O);!$(N-a)]([CX4])([CX4])[CX4]",
            "POSITIVE CONTROL -- ionizable lipids are full of these; must be "
            "in the thousands or the sweep is reading nothing",
            rr.REGISTRY["DLin-MC3-DMA"][0], "DLin-MC3-DMA"),
        "sterol_skeleton": (
            # Derived FROM the registry cholesterol by Murcko scaffold with
            # bond orders and element identities generified, rather than
            # hand-written: a hand-written ring system is exactly what failed.
            None,
            "POSITIVE CONTROL -- cholesterol is the corpus's near-universal "
            "sterol; must be non-zero AND must match the registry cholesterol",
            rr.REGISTRY["Cholesterol"][0], "Cholesterol"),
    }
    from rdkit.Chem.Scaffolds import MurckoScaffold
    _chol = Chem.MolFromSmiles(rr.REGISTRY["Cholesterol"][0])
    _scaf = MurckoScaffold.MakeScaffoldGeneric(
        MurckoScaffold.GetScaffoldForMol(_chol))
    # MakeScaffoldGeneric writes EXPLICIT SINGLE BONDS, so the emitted SMARTS
    # cannot match cholesterol's own delta-5 alkene -- the guard below caught
    # exactly that. Every atom in this scaffold is [#6], so there is no charge
    # or explicit-H hyphen to damage: replacing every '-' with '~' makes the
    # ring topology bond-order agnostic, which is what "sterol skeleton" means.
    controls["sterol_skeleton"] = (
        Chem.MolToSmarts(_scaf).replace("-", "~"),) + \
        controls["sterol_skeleton"][1:]
    # one pass over distinct structures, all component columns
    smiles = set()
    for pool, df in corpus.items():
        for sc in SMI_COLS:
            smiles |= set(df[sc].dropna().astype(str))
    mols = {}
    unparseable = 0
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        if m is None:
            unparseable += 1
        else:
            mols[s] = m

    def run(smarts):
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            raise SystemExit(f"bad SMARTS, refusing to report a 0: {smarts}")
        return [s for s, m in mols.items() if m.HasSubstructMatch(patt)]

    out = {"n_distinct_structures_scanned": len(mols),
           "n_unparseable_skipped": unparseable,
           "label": "SUPERSET -- a chemotype match is not the same molecule",
           "queries": {}, "positive_controls": {}}
    for name, (sm, why, must_match_smiles, must_match_name) in controls.items():
        patt = Chem.MolFromSmarts(sm)
        ref = Chem.MolFromSmiles(must_match_smiles)
        premise = bool(patt is not None and ref is not None
                       and ref.HasSubstructMatch(patt))
        if not premise:
            raise SystemExit(
                f"POSITIVE CONTROL {name!r} does not match the molecule it "
                f"names ({must_match_name}). Refusing to report any sweep "
                f"count: a control that cannot match its own reference "
                f"certifies nothing, and a non-zero corpus count would make "
                f"that invisible. SMARTS: {sm}")
        hits = run(sm)
        out["positive_controls"][name] = {
            "smarts": sm, "why": why,
            "premise_molecule": must_match_name,
            "premise_smiles": must_match_smiles,
            "premise_verified_smarts_matches_that_molecule": premise,
            "n_structures": len(hits),
            "reference_structure_found_in_corpus": must_match_smiles in mols,
            "passes": len(hits) > 0 and premise}
    for name, (sm, why) in queries.items():
        hits = run(sm)
        out["queries"][name] = {
            "smarts": sm, "why": why, "n_structures": len(hits),
            "example_structures": sorted(hits)[:3]}
    out["verdict"] = ("Even a non-zero chemotype count would not be an overlap: "
                      "7C1 has no atom-resolved structure, so it cannot be "
                      "shown identical to anything. The loose counts are "
                      "reported as a superset and claimed as nothing more.")
    return out


# --------------------------------------------------- in-vivo shared check ---

def invivo_shared(corpus):
    """The in-vivo pool has exactly two structures shared across experiments.
    Measure that, so the claim is this run's measurement and not a remembered
    number, and then say where this paper would have landed."""
    df = corpus["vivo"]
    sub = df[["Experiment_ID", "IL_SMILES", "IL_name"]].dropna(subset=["IL_SMILES"])
    keys = {}
    for s in sub["IL_SMILES"].astype(str).unique():
        keys[s] = ikey(s)
    sub = sub.assign(ikey=sub["IL_SMILES"].astype(str).map(keys))
    g = sub.dropna(subset=["ikey"]).groupby("ikey")
    shared = {k: {"n_sources": int(v.Experiment_ID.nunique()),
                  "sources": sorted(v.Experiment_ID.unique()),
                  "names": sorted(set(v.IL_name.astype(str))),
                  "rows": int(len(v))}
              for k, v in g if v.Experiment_ID.nunique() > 1}
    return {
        "n_structures_shared_across_invivo_experiments": len(shared),
        "shared": shared,
        "this_paper_would_add": 0,
        "why": ("7C1 has no writable structure, so it cannot join this set even "
                "in principle, and this paper emits no in-vivo row anyway. It "
                "neither adds a third shared in-vivo structure nor raises the "
                "source count on either existing one."),
    }


def main():
    corpus = load_corpus()
    out = {
        "_what": ("Corpus-overlap measurement for PMID 30216729 (Sago 2018, "
                  "Nano Letters). Every query carries a positive control."),
        "corpus_shape": {p: list(df.shape) for p, df in corpus.items()},
        "token_pattern_test": token_test(corpus),
        "structural_attestation": structural_join(corpus),
        "substructure_sweep_superset": substructure_sweep(corpus),
        "invivo_shared_structures": invivo_shared(corpus),
    }
    os.makedirs(OUTDIR, exist_ok=True)
    with open(os.path.join(OUTDIR, "02_corpus_overlap.json"), "w") as f:
        json.dump(out, f, indent=1)
    t = out["token_pattern_test"]
    print(json.dumps({
        "token_7C1_substring_rows": t["7C1"]["substring_rows"],
        "token_7_substring_rows": t["7"]["substring_rows"],
        "token_C14_substring_rows": t["C14"]["substring_rows"],
        "token_C18_substring_rows": t["C18"]["substring_rows"],
        "exact_name_match": t["_exact_name_match"],
        "token_positive_control_DSPC": t["_positive_control"],
        "inchikey_positive_control_MC3_rows":
            out["structural_attestation"]["positive_control"].get("corpus_rows_by_inchikey"),
        "inchikey_positive_control_passes":
            out["structural_attestation"]["positive_control"].get("passes"),
        "cholesterol_rows": out["structural_attestation"]["queried"]["Cholesterol"]["corpus_rows_by_inchikey"],
        "dspc_rows": out["structural_attestation"]["queried"]["DSPC"]["corpus_rows_by_inchikey"],
        "sweep_controls": {k: v["n_structures"] for k, v in
                           out["substructure_sweep_superset"]["positive_controls"].items()},
        "sweep_queries": {k: v["n_structures"] for k, v in
                          out["substructure_sweep_superset"]["queries"].items()},
        "invivo_shared": {k: v["names"] for k, v in
                          out["invivo_shared_structures"]["shared"].items()},
    }, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
