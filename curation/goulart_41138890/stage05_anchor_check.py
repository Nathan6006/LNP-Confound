#!/usr/bin/env python3
"""C10 by hand, in this folder (provenance/structure_check.py is off limits
for this run), PLUS the two independence guards MADE TO FAIL ON PURPOSE.

What is being tested, and what is not:

  TESTED   whether the structure the corpus carries under the name `C12-200`
           is the molecule Figure 1A draws. The anchor is built by
           make_reference.py, which reads no corpus file.
  TESTED   whether the comparison can see a CONNECTIVITY error at identical
           composition -- seeded as a Markovnikov epoxide opening, which must
           come back Delta formula = none and Delta exact mass = 0.000000 Da
           with a different InChIKey.
  TESTED   whether the AST independence guard and the template-provenance
           guard can fail at all. Each is seeded in a CHILD PROCESS writing to
           a temp directory, and a clean control must come back clean.
  NOT      SM-102, ALC-0315, cholesterol, DOPE, DSPC, DOTAP. The paper draws
  TESTED   the excipients but they are pinned from provenance/reagent_registry,
           and comparing a registry pin against itself is not an anchor. SM-102
           and ALC-0315 are not drawn at all. Against an empty anchor the
           guards are VACUOUS, not passing.
"""
import csv
import json
import os
import subprocess
import sys
import tempfile

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
sys.path.insert(0, HERE)
import make_reference as mr                            # noqa: E402

VITRO = os.path.join(HERE, "..", "..", "data", "data_vitro.csv")


def corpus_c12_200():
    """Every distinct structure the corpus carries under the name C12-200."""
    out = {}
    with open(VITRO, newline="") as fh:
        for r in csv.DictReader(fh):
            if (r.get("IL_name") or "").strip() == "C12-200":
                s = r["IL_SMILES"].strip()
                k = Chem.MolToInchiKey(Chem.MolFromSmiles(s))
                e = out.setdefault(k, {"smiles": s, "rows": 0, "exps": set(),
                                       "stored_MW": set()})
                e["rows"] += 1
                e["exps"].add(r["Experiment_ID"])
                e["stored_MW"].add(r.get("Molecular.Weight", ""))
    for v in out.values():
        v["exps"] = sorted(v["exps"])
        v["stored_MW"] = sorted(v["stored_MW"])
    return out


def describe(smiles):
    m = Chem.MolFromSmiles(smiles)
    return {"inchikey": Chem.MolToInchiKey(m),
            "skeleton_block": Chem.MolToInchiKey(m).split("-")[0],
            "formula": rdMolDescriptors.CalcMolFormula(m),
            "exact_mass": round(Descriptors.ExactMolWt(m), 6)}


def child(seed=None, mutate=None):
    """Run make_reference.py in a temp dir, optionally on a MUTATED copy."""
    td = tempfile.mkdtemp(prefix="goulart_anchor_")
    script = os.path.join(HERE, "make_reference.py")
    if mutate is not None:
        src = open(script).read()
        old, new = mutate
        assert old in src, f"mutation anchor {old!r} not found"
        script = os.path.join(td, "make_reference_seeded.py")
        open(script, "w").write(src.replace(old, new, 1))
    cmd = [PY, "-I", script, "--outdir", td]
    if seed:
        cmd += ["--seed", seed]
    p = subprocess.run(cmd, capture_output=True, text=True, cwd=HERE)
    return p.returncode, p.stdout, p.stderr, td, script


def main():
    out = {}

    # ---- 1. the anchor vs the corpus ------------------------------------
    anchor, _ = mr.build()
    corpus = corpus_c12_200()
    a = describe(anchor)
    out["anchor"] = a
    out["corpus_C12_200"] = {k: {"rows": v["rows"], "experiments": v["exps"],
                                 "stored_Molecular_Weight": v["stored_MW"],
                                 **describe(v["smiles"])}
                             for k, v in corpus.items()}
    out["n_distinct_corpus_structures_under_this_name"] = len(corpus)
    match = [k for k in corpus if k == a["inchikey"]]
    out["C10_result"] = ("exact_match_full_inchikey" if match
                         else "MISMATCH")
    out["C10_matched_rows"] = (corpus[match[0]]["rows"] if match else 0)

    # ---- 2. can the comparison see connectivity at equal composition? ----
    seeded, _ = mr.build(seed="markovnikov")
    s = describe(seeded)
    out["connectivity_seed"] = {
        "what": ("Markovnikov opening of the C12 epoxide: "
                 "N-CH(C10H21)-CH2OH instead of N-CH2-CH(OH)-C10H21"),
        "seeded": s, "clean": a,
        "delta_exact_mass_Da": round(s["exact_mass"] - a["exact_mass"], 6),
        "formula_identical": s["formula"] == a["formula"],
        "skeleton_block_identical": s["skeleton_block"] == a["skeleton_block"],
        "inchikey_differs": s["inchikey"] != a["inchikey"]}
    out["connectivity_seed"]["credited"] = (
        out["connectivity_seed"]["formula_identical"]
        and out["connectivity_seed"]["delta_exact_mass_Da"] == 0.0
        and out["connectivity_seed"]["inchikey_differs"])

    # ---- 3. the two guards, made to fail on purpose ---------------------
    enumerated = [v["smiles"] for v in corpus.values()] + [anchor]
    guards = {}
    guards["clean_control"] = {
        "ast_independence_violations": mr.check_independence(),
        "lifted_template_literals": mr.check_no_lifted_template(enumerated)}
    guards["clean_control"]["clean"] = (
        not guards["clean_control"]["ast_independence_violations"]
        and not guards["clean_control"]["lifted_template_literals"])

    # 3a. seed an import of the corpus helper
    rc, so, se, td, script = child(
        mutate=('BANNED_IMPORTS = {"corpus_components"',
                'import corpus_components  # SEEDED\nBANNED_IMPORTS = '
                '{"corpus_components"'))
    v = mr.check_independence(script)
    guards["seed_banned_import"] = {
        "violations": v, "fires": bool(v),
        "note": ("the guard is run against the MUTATED FILE's AST; a guard "
                 "that could not see its own seeded import would report the "
                 "clean module and pass (the F-037 shape)")}

    # 3b. seed a long literal lifted out of the corpus structure
    lifted = sorted(corpus.values(), key=lambda x: -x["rows"])[0]["smiles"][:40]
    rc2, so2, se2, td2, script2 = child(
        mutate=('    core = "NCCNCCN1CCN(CCN)CC1"',
                f'    _lifted = "{lifted}"\n'
                '    core = "NCCNCCN1CCN(CCN)CC1"'))
    v2 = mr.check_no_lifted_template(enumerated, script2)
    guards["seed_lifted_template"] = {
        "literal_seeded_len": len(lifted), "violations": v2,
        "fires": bool(v2)}
    guards["both_guards_fire"] = (guards["seed_banned_import"]["fires"]
                                  and guards["seed_lifted_template"]["fires"]
                                  and guards["clean_control"]["clean"])
    out["independence_guards"] = guards

    # ---- 4. the emitted anchor file -------------------------------------
    ref = os.path.join(HERE, "reference_structures.csv")
    with open(ref, newline="") as fh:
        rows = list(csv.DictReader(fh))
    out["reference_structures_csv"] = {
        "header_exact": list(rows[0].keys()) == [
            "lipid_name", "reference_SMILES", "source_figure", "read_on",
            "read_by", "expected_mismatch", "note"],
        "n_rows": len(rows),
        "n_with_structure": sum(1 for r in rows if r["reference_SMILES"]),
        "n_empty_on_purpose": sum(1 for r in rows
                                  if not r["reference_SMILES"]),
        "expected_mismatches": [r["lipid_name"] for r in rows
                                if r["expected_mismatch"]]}

    with open(os.path.join(HERE, "05_anchor_check.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
    ok = (out["C10_result"] == "exact_match_full_inchikey"
          and out["connectivity_seed"]["credited"]
          and guards["both_guards_fire"]
          and out["reference_structures_csv"]["header_exact"])
    print("VERDICT:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
