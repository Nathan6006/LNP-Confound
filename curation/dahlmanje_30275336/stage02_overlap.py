"""Corpus overlap for Sago 2018 (PMID 30275336).

Three levels, each labelled for what it actually proves:

  EXACT      full InChIKey identity   -- IMPOSSIBLE here: this paper yields no
                                         determinable SMILES for any of its 8
                                         lipid-amine compounds (see DECISIONS
                                         D-2), so there is nothing to hash.
  NAME       IL_name identity and containment against the live corpus.
  SUPERSET   SMARTS substructure sweep. Deliberately labelled a SUPERSET: it
             answers "does the corpus contain this CLASS of chemistry", not
             "does it contain these molecules". A sibling curation reported a
             loose SMARTS hit count under a label claiming exactness (332 vs a
             true 128); the label here is the fix.

Excipients ARE determinable and ARE compared, through the guarded lookup in
curation/corpus_components.py, which RAISES on an ambiguous name (error class
30) rather than resolving it by row order.
"""
import sys, os, json
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

ROOT = sys.argv[1]
sys.path.insert(0, os.path.join(ROOT, "curation"))
from corpus_components import load, smi            # noqa: E402

vi = pd.read_csv(os.path.join(ROOT, "data", "data_vitro.csv"), low_memory=False)
vv = pd.read_csv(os.path.join(ROOT, "data", "data_vivo.csv"), low_memory=False)

# ---------------------------------------------------------------- NAME level
COMPOUNDS = ["7C1", "73", "78", "80", "83", "88", "90", "92"]
LNP_NAMES = ["7C2", "7C3", "7C4"]
TOKENS = COMPOUNDS + LNP_NAMES

name_report = {}
for pool, d in (("vitro", vi), ("vivo", vv)):
    nm = d["IL_name"].astype(str).str.strip()
    for t in TOKENS:
        exact = int((nm == t).sum())
        # honest superset: any name CONTAINING the token anywhere
        contain = nm.str.contains(t, regex=False, na=False)
        key = f"{pool}:{t}"
        name_report[key] = {
            "exact_IL_name_matches": exact,
            "names_CONTAINING_token_SUPERSET": int(contain.sum()),
            "distinct_names_containing": sorted(set(nm[contain]))[:8],
            "experiments_containing": sorted(set(d.loc[contain, "Experiment_ID"])),
        }

# ------------------------------------------------- the name-collision warning
# The brief's F-055 note records that lipidoid AMINE catalogue numbers are
# shared across labs (93, 9322, 80). This paper's 73/78/80/83/88/90 look like
# the same catalogue. They are NOT: here they index COMPOUNDS (amine x tail
# chemistry), and the same index 80 denotes a PEI1200 + C12-acrylamide polymer.
# Measured, not eyeballed: how many corpus rows would a naive IL_name join on
# the bare token capture, and how many distinct structures would it span?
collision = {}
for pool, d in (("vitro", vi), ("vivo", vv)):
    nm = d["IL_name"].astype(str).str.strip()
    for t in ["73", "78", "80", "83", "88", "90"]:
        c = nm.str.contains(t, regex=False, na=False)
        collision[f"{pool}:{t}"] = {
            "rows_a_naive_substring_join_would_capture": int(c.sum()),
            "distinct_structures_those_rows_span": int(d.loc[c, "IL_SMILES"].nunique()),
            "experiments": sorted(set(d.loc[c, "Experiment_ID"])),
        }

# ------------------------------------------------------------ SUPERSET sweep
# Motifs this paper's chemistry would produce IF a structure were determinable.
SMARTS = {
  "epoxide_opened_2hydroxyalkyl_on_N": "[NX3;!$(N=*)][CH2][CH1]([OX2H1])[CX4]",
  "michael_adduct_3oxopropyl_ester_on_N": "[NX3;!$(N=*)][CH2][CH2]C(=O)[OX2][CX4]",
  "michael_adduct_3oxopropyl_amide_on_N": "[NX3;!$(N=*)][CH2][CH2]C(=O)[NX3]",
  "aminoethyl_aminoethyl_TETA_backbone": "[NX3][CH2][CH2][NX3][CH2][CH2][NX3]",
}
sweep = {}
for pool, d in (("vitro", vi), ("vivo", vv)):
    uniq = d[["IL_SMILES"]].drop_duplicates()["IL_SMILES"].tolist()
    mols = []
    unparsed = 0
    for s in uniq:
        m = Chem.MolFromSmiles(s)
        if m is None:
            unparsed += 1
        else:
            mols.append((s, m))
    for label, sma in SMARTS.items():
        q = Chem.MolFromSmarts(sma)
        assert q is not None, label
        hits = [s for s, m in mols if m.HasSubstructMatch(q)]
        rows = int(d["IL_SMILES"].isin(set(hits)).sum())
        sweep[f"{pool}:{label}"] = {
            "LABEL": "SUPERSET - class membership, NOT identity with this paper's compounds",
            "distinct_corpus_structures_matching": len(hits),
            "corpus_rows_matching": rows,
            "experiments": sorted(set(d.loc[d["IL_SMILES"].isin(set(hits)), "Experiment_ID"])),
        }
    sweep[f"{pool}:_unparsed_structures"] = unparsed

# ------------------------------------------------------------- EXCIPIENTS
# corpus_components.load() calls .strip() on every *_name cell; pandas reads a
# missing name as float NaN. Fill in OUR frame rather than touch the shared
# module (it is outside this curation folder and other agents read it).
_both = pd.concat([vi, vv], ignore_index=True)
_fillcols = [c for c in _both.columns if c.endswith("_name") or c.endswith("_SMILES")]
_both[_fillcols] = _both[_fillcols].fillna("")
rows_all = _both.to_dict("records")
comp, ambiguous = load(rows_all)
EXCIPIENTS = ["DOPE", "DOPC", "DSPC", "Cholesterol", "cholesterol",
              "18:1 Lyso PC", "C14-PEG2000", "DMG-PEG2000"]
exc = {}
for name in EXCIPIENTS:
    try:
        exc[name] = {"found": True, "corpus_SMILES": smi(comp, ambiguous, name)}
    except BaseException as e:
        exc[name] = {"found": False, "why": type(e).__name__ + ": " + str(e)[:200]}

out = {
  "_scope": "Sago et al. PNAS 2018, PMID 30275336. Compared against the LIVE corpus.",
  "corpus_shape": {"vitro": list(vi.shape), "vivo": list(vv.shape)},
  "EXACT_InChIKey_overlap": {
      "performed": False,
      "why": "This paper yields NO determinable SMILES for any of its 8 "
             "lipid-amine compounds (7C1, 73, 78, 80, 83, 88, 90, 92). Four "
             "are PEI600/PEI1200 conjugates (polydisperse branched polymers); "
             "three are triethylenetetramine conjugates and one a macrocyclic "
             "tetraamine conjugate, all drawn as a single-repeat-unit cartoon "
             "with the backbone truncated to two stubs and NO stated degree of "
             "substitution. There is nothing to canonicalise, so an InChIKey "
             "comparison would be a comparison of invented molecules.",
      "full_InChIKey_matches": None,
      "skeleton_InChIKey_matches": None},
  "NAME_overlap": name_report,
  "NAME_TOKEN_COLLISION_WARNING": collision,
  "SUPERSET_substructure_sweep": sweep,
  "EXCIPIENT_overlap_guarded_lookup": exc,
  "EXPERIMENT_ID_collision_check": {
      "candidate_ids_for_this_paper": ["CS_2018"],
      "vitro_ids_in_use": sorted(set(vi["Experiment_ID"])),
      "vivo_ids_in_use": sorted(set(vv["Experiment_ID"])),
      "CS_2018_in_use": "CS_2018" in set(vi["Experiment_ID"]) | set(vv["Experiment_ID"]),
      "note": "SS_2018 (PMID 29653760) already exists in the vivo pool and is a "
              "DIFFERENT paper; checked by PMID, not by the year token."},
  "POOLED_READOUT_PRECEDENT": {
      "finding": "The corpus schema already carries an Experiment_batching "
                 "column whose vocabulary is {barcoded, individual}.",
      "LR_2023": {"pmid": 37127709, "rows": int((vv["Experiment_ID"] == "LR_2023").sum()),
                  "Experiment_batching": sorted(set(vv.loc[vv["Experiment_ID"] == "LR_2023",
                                                           "Experiment_batching"].astype(str))),
                  "Cargo_type": sorted(set(vv.loc[vv["Experiment_ID"] == "LR_2023",
                                                  "Cargo_type"].astype(str)))},
      "consequence": "A POOLED multiplexed in vivo barcode screen is ALREADY in "
                     "the corpus, declared as such. Pooling is therefore a "
                     "labelled property, not grounds for refusal. See D-1."},
}
print(json.dumps(out, indent=1, default=str))
