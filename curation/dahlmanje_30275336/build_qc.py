"""Consolidate every check into 07_qc.json. Every entry carries NUMBERS, and
every check that was seeded records whether the AGREEMENT COUNT MOVED --
not merely whether the seed raised. Five no-op seeds were caught across this
project in one day by exactly that distinction, and two of them were in this
folder (see chain_check.seed_history).
"""
import json, subprocess, sys, csv, collections, os

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "pnas.1811276115.sapp.pdf")
ROOT = sys.argv[1]

def run(args):
    p = subprocess.run([sys.executable, "-I"] + args, capture_output=True,
                       text=True, cwd=HERE)
    return p.returncode, p.stdout, p.stderr

qc = {"paper": {"pmid": "30275336", "pmcid": "PMC6196543",
                "doi": "10.1073/pnas.1811276115",
                "title": "High-throughput in vivo screen of functional mRNA "
                         "delivery identifies nanoparticles for endothelial "
                         "cell gene editing.",
                "licence": "CC BY-NC-ND 4.0 (NOT CC0/CC-BY compatible)"}}

# ---- 1. the readout reconstruction ---------------------------------------
seeds04 = {}
for s in ["none", "wrong_raw_count", "drop_input_normalisation", "no_renormalisation"]:
    rc, out, err = run(["stage04_readout.py", s])
    d = json.loads(out)
    seeds04[s] = {"rc": rc,
        "step1_totals": d["step1_column_totals"]["AGREEMENT"],
        "step2_cells": d["step2_share_of_sample"]["AGREEMENT"],
        "step3_cells": d["step3_normalised_to_input_then_RENORMALISED"]["AGREEMENT"],
        "max_residual_from_100":
            d["THE_FREE_CHECK_the_authors_drew_but_did_not_compute"]["max_abs_residual_from_100"]}
base = seeds04["none"]
qc["readout_reconstruction"] = {
  "what": "Reconstruct the paper's own definition of 'normalized delivery' from "
          "its worked example (SI Fig S1b, p.2) and reproduce every printed cell.",
  "control": base,
  "FORMULA": json.loads(run(["stage04_readout.py"])[1])["FORMULA_ESTABLISHED"],
  "seeds": seeds04,
  "seed_verdict": {
     s: {"agreement_moved": (v["step1_totals"], v["step2_cells"], v["step3_cells"]) !=
                            (base["step1_totals"], base["step2_cells"], base["step3_cells"])
                            or v["max_residual_from_100"] != base["max_residual_from_100"]}
     for s, v in seeds04.items() if s != "none"},
  "THE_FREE_CHECK": "The authors DREW 'Total (%) = 100' under both normalised "
          "tables and never computed it. The unrounded reconstruction gives "
          "column sums of exactly 100.000000, residual 0.0. The readout is a "
          "COMPOSITIONAL share of its pool.",
}

# ---- 2. formulation tables + mol% sums -----------------------------------
rc, out, err = run(["stage01_tables.py", PDF, os.path.join(HERE, "01_formulations.csv")])
rep = json.loads(out)
qc["formulation_tables"] = {
  "tier": "4_pdfplumber_native_text_cells (no pixels)",
  "n_formulations": rep["n_formulations_total"],
  "per_library": {m["library"]: {"n": m["n_rows"], "lnp_numbers": "1..%d" % m["lnp_numbers_max"],
                                 "contiguous": m["contiguous_1_to_max"]}
                  for m in rep["per_library"]},
  "compounds_seen": rep["compounds_seen"],
  "peg_types_seen": rep["peg_types_seen"],
  "helper_types_seen": rep["helper_types_seen"],
  "MOLPCT_SUMS_TO_100": {
     "this_is_the_free_check": "A composition must sum to 100. The paper prints "
        "no such sum anywhere; the stacked-bar panels Fig 2C and SI Fig S1i draw "
        "it reaching exactly 100.",
     "n_checked": rep["molpct_sums_to_100"]["n_checked"],
     "n_exactly_100": rep["molpct_sums_to_100"]["n_exactly_100"],
     "n_off": rep["molpct_sums_to_100"]["n_off"],
     "offenders": rep["molpct_sums_to_100"]["offenders"],
     "verdict": "%d/%d exact. The 7 exceptions are the paper's, not the "
                "extraction's: pdfplumber's table grid and an independent "
                "extract_text() read of the same page return identical digits. "
                "See 08_contradictions.csv C-4 and C-5."
                % (rep["molpct_sums_to_100"]["n_exactly_100"],
                   rep["molpct_sums_to_100"]["n_checked"]),
  },
}

# ---- 3. chain check -------------------------------------------------------
chain = {}
for s in ["none", "drop_a_real_bond", "inject_a_bond_sized_glyph", "inject_a_small_glyph"]:
    rc, out, err = run(["stage03_chaincheck.py", PDF, s])
    d = json.loads(out)
    chain[s] = {"rc": rc, "bond_quads": d["n_full_bond_quads"],
                "CARBONS": d["MEASURED_acyl_carbons_including_carbonyl"],
                "contiguous": d["tiling_is_contiguous"],
                "token_matches": d["token_matches_measurement"]}
qc["chain_check"] = {
  "what": "The PEG column's token is C14PEG. 'A name token is not a chain length', "
          "so the drawn chain is measured.",
  "mode": "the shared curation/lib/chains.vector_bonds REFUSED in its default "
          "stroked-only mode and chains.vector_chain_length raised -- correctly. "
          "This page is a second measured instance of the class that module "
          "documents from andersondg_41203968: every bond is a FILLED polygon "
          "with width=None (603 of page 7's 915 paths). The refusal was "
          "respected, not widened silently; the measurement instead counts the "
          "PDF's own path bounding boxes under an explicit uniformity gate.",
  "control": chain["none"],
  "seeds": chain,
  "seed_history": "THE FIRST TWO SEEDS WRITTEN HERE WERE NO-OPS and both "
          "printed a clean 14: paths[:-1] dropped the short carbonyl-to-oxygen "
          "stub rather than a bond, and a small glyph stub is SUPPOSED to be "
          "excluded. They were replaced with seeds that move the count. "
          "inject_a_small_glyph is kept as the negative control.",
  "seed_verdict": {s: {"CARBONS_moved": chain[s]["CARBONS"] != chain["none"]["CARBONS"],
                       "expected_to_move": s in ("drop_a_real_bond", "inject_a_bond_sized_glyph")}
                   for s in chain if s != "none"},
  "RESULT": "14 carbons including the carbonyl = myristoyl; the C14 token is "
            "validated. Combined with the drawn phosphoethanolamine-carbamate "
            "head, the PEG-lipid is 14:0 PEG2000 PE (DMPE-PEG2000), which is "
            "NOT what the corpus stores under the name 'C14-PEG2000'. C-6.",
}

# ---- 4. anchors -----------------------------------------------------------
anch = {}
for s in ["none", "seed_epimer", "seed_banned_call"]:
    rc, out, err = run(["make_reference.py", s])
    d = json.loads(out)
    a = d["anchors"]["Cholesterol"]
    anch[s] = {"AST_CLEAN": d["AST_independence_gate"]["CLEAN"],
               "banned_calls": d["AST_independence_gate"]["banned_calls_found"],
               "CHOL_full_InChIKey_match": a["FULL_InChIKey_match"],
               "CHOL_skeleton_match": a["skeleton_match"],
               "CHOL_formula_match": a["formula_match"],
               "CHOL_exact_mass_match": a["exact_mass_match"],
               "DOPE_full_InChIKey_match": d["anchors"]["DOPE"]["FULL_InChIKey_match"]}
qc["c10_anchors"] = {
  "ionizable_anchors_with_a_SMILES": 0,
  "ionizable_anchors_EMPTY_ON_PURPOSE": 8,
  "why_empty": "The paper draws no complete ionizable-lipid structure - only a "
               "single repeat-unit cartoon with the backbone truncated. "
               "CLAUDE.md section 5: store an empty reference_SMILES so the "
               "ABSENCE of an anchor is visible.",
  "excipient_anchors_built": ["Cholesterol", "DOPE"],
  "independence_method": "AST walk of make_reference.py at run time. NOT a "
               "substring scan - a sibling's substring scan flagged its own "
               "docstring, and this module's docstring names every banned "
               "function. NOT byte-identity either: that is the F-037 trap, "
               "because a correct re-derivation canonicalises to the same "
               "bytes as a correct enumeration -- and it did here "
               "(canonical_SMILES_byte_identical = true for both anchors, "
               "which is the expected outcome and proves nothing on its own).",
  "control": anch["none"],
  "seeds": anch,
  "seed_verdict": {
    "seed_epimer": {"fires": anch["seed_epimer"]["CHOL_full_InChIKey_match"] is False,
       "and_crucially": "formula_match=%s, exact_mass_match=%s, skeleton_match=%s "
            "ALL REMAIN TRUE. The seeded defect has an IDENTICAL molecular "
            "formula (C27H46O) and an IDENTICAL exact mass (386.3549). This is "
            "F-013 exactly -- epicholesterol on 537 rows -- and only the FULL "
            "InChIKey sees it."
            % (anch["seed_epimer"]["CHOL_formula_match"],
               anch["seed_epimer"]["CHOL_exact_mass_match"],
               anch["seed_epimer"]["CHOL_skeleton_match"])},
    "seed_banned_call": {"fires": anch["seed_banned_call"]["AST_CLEAN"] is False}},
}

# ---- 5. emit guards -------------------------------------------------------
rc, out, err = run(["stage06_emit.py", ROOT])
qc["emit"] = json.loads(out)

# ---- 6. values ledger -----------------------------------------------------
vals = list(csv.DictReader(open(os.path.join(HERE, "05_values.csv"))))
ext = [v for v in vals if v["tier"] != "NOT_EXTRACTED"]
qc["values_ledger"] = {
  "total_rows_in_05_values": len(vals),
  "by_tier": dict(collections.Counter(v["tier"] for v in vals)),
  "numeric_values_extracted": len(ext),
  "of_which_duplicates_of_another_row": sum(1 for v in ext if v["duplicate_of"].startswith("V")),
  "distinct_numerals": len(ext) - sum(1 for v in ext if v["duplicate_of"].startswith("V")),
  "individually_measured": dict(collections.Counter(v["individually_measured"] for v in vals)),
  "corpus_admissible": dict(collections.Counter(v["corpus_admissible"] for v in vals)),
  "COUNTS_SUM": len(ext) + sum(1 for v in vals if v["tier"] == "NOT_EXTRACTED") == len(vals),
}

# ---- 7. THE THREE NUMBERS -------------------------------------------------
enum = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
forms = list(csv.DictReader(open(os.path.join(HERE, "01_formulations.csv"))))
qc["THE_THREE_NUMBERS"] = {
  "distinct_ionizable_STRUCTURES": {
     "named_in_the_paper": len(enum),
     "with_a_determinable_SMILES": sum(1 for e in enum if e["IL_SMILES"].strip()),
     "note": "8 compounds are NAMED (7C1, 73, 78, 80, 83, 88, 90, 92). ZERO "
             "have a determinable structure. Four are PEI600/PEI1200 "
             "conjugates (polydisperse branched polymers), three are "
             "triethylenetetramine conjugates and one a macrocyclic-tetraamine "
             "conjugate, every one drawn as a truncated repeat unit with no "
             "stated degree of substitution."},
  "distinct_FORMULATIONS": {
     "tabulated_in_the_SI": len(forms),
     "per_library": dict(collections.Counter(f["library"] for f in forms)),
     "note": "Fully recoverable, tier 4, mol% verified. Formulations are the "
             "one thing this paper documents completely."},
  "INDIVIDUALLY_MEASURED": {
     "LNPs": 2, "which": ["7C2", "7C3"],
     "ionizable_compound_they_share": "7C1 (both; they differ only in helper "
             "lipid and mol ratio, so the ionizable structural variance across "
             "the individually-measured set is ZERO)",
     "printed_values": len(ext),
     "note": "Out of 425 tabulated formulations and 255 pooled in vivo, exactly "
             "2 LNPs were measured individually, and they share one compound."},
  "WHAT_THE_BARCODE_DESIGN_CAN_SUPPORT": "A pooled competitive design can in "
     "principle support a per-LNP value for every pooled member -- 255 in vivo "
     "plus 54 in vitro. It does NOT here, because the paper publishes those "
     "values only as unlabelled points and an unlabelled clustered heat map. "
     "The design supports 309 values; the publication attributes 0 of them.",
}

with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
    json.dump(qc, fh, indent=1)
print(json.dumps({k: (v if k == "THE_THREE_NUMBERS" else "written")
                  for k, v in qc.items()}, indent=1)[:2600])
