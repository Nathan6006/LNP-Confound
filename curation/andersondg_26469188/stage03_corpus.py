#!/usr/bin/env python3
"""KK_2015 -- corpus overlap, the token-pattern test, and the component table.

READ-ONLY on data/. Nothing here writes to the corpus.

F-070: every corpus query carries a POSITIVE CONTROL -- something already known
to be present must come back. A clean zero from a query that cannot match is
indistinguishable from a true zero.
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from rdkit import Chem, RDLogger                                  # noqa: E402
from rdkit.Chem import rdMolDescriptors as rdm                    # noqa: E402

from corpus_components import load, smi                           # noqa: E402

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.environ.get("KK2015_OUTDIR", HERE)
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DATA = os.path.join(ROOT, "data")

# the ionizable lipid this paper uses, and the helper lipids it measures
IL_NAME = "C12-200"
HELPERS = ["DSPC", "DOPE", "DOPC"]          # DSPE is designed but never made
HELPERS_DESIGNED_ONLY = ["DSPE"]


def read(path, pool):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["_pool"] = pool
    return rows


def inchikey(s):
    m = Chem.MolFromSmiles(s)
    return (None, None) if m is None else (Chem.MolToInchiKey(m), m)


def main():
    vitro = read(os.path.join(DATA, "data_vitro.csv"), "vitro")
    vivo = read(os.path.join(DATA, "data_vivo.csv"), "vivo")
    corpus = vitro + vivo
    out = {"_what": "KK_2015 overlap against the live corpus. Read-only.",
           "corpus_shape": {"vitro_rows": len(vitro), "vivo_rows": len(vivo),
                            "vitro_experiments": len({r["Experiment_ID"] for r in vitro}),
                            "vivo_experiments": len({r["Experiment_ID"] for r in vivo}),
                            "total_experiment_IDs":
                                len({r["Experiment_ID"] for r in corpus})}}

    # ------------------------------------------------- experiment ID check --
    ids = sorted({r["Experiment_ID"] for r in corpus})
    pmids = {}
    for r in corpus:
        pmids.setdefault(r["Experiment_ID"], set()).add(r.get("Publication_PMID", ""))
    out["experiment_id"] = {
        "proposed": "KK_2015",
        "derivation": "FIRST author Kauffman, Kevin J. + 2015. The folder is named "
                      "andersondg_26469188 after the PI, which is not the corpus key.",
        "collides_with_existing_ID": "KK_2015" in ids,
        "n_existing_IDs": len(ids),
        "all_existing_IDs": ids,
        "checked_by_PMID_not_year_token": True,
        "pmid_26469188_already_in_corpus":
            any("26469188" in {str(x).split(".")[0] for x in v} for v in pmids.values()),
        "nearest_tokens_examined": {
            i: sorted({str(x).split(".")[0] for x in pmids[i] if x})
            for i in ids if i.startswith("K") or i.endswith("2015")},
    }

    # ------------------------------------------------------ InChIKey joins --
    uniq = sorted({r["IL_SMILES"] for r in corpus if r.get("IL_SMILES")})
    keys = {}
    for s in uniq:
        k, _ = inchikey(s)
        if k:
            keys[s] = k
    for r in corpus:
        r["_IK"] = keys.get(r.get("IL_SMILES"))
    out["inchikey_machinery"] = {"unique_IL_SMILES": len(uniq),
                                 "unparseable": len(uniq) - len(keys)}

    # POSITIVE CONTROL #1 -- cKK-E12, documented to sit under 3 names
    ck = [r for r in corpus if r["IL_name"] == "CKK-E12"]
    ckk = {r["_IK"] for r in ck}
    ckjoin = [r for r in corpus if r["_IK"] in ckk]
    out["positive_control_1_cKK_E12"] = {
        "why": "F-070. A structure can sit under several names; this one is "
               "documented as CKK-E12 / other / MD-1. If the InChIKey join cannot "
               "reproduce that, the join is broken and a zero means nothing.",
        "name_join_rows": len(ck),
        "inchikey_join_rows": len(ckjoin),
        "names_recovered": sorted({r["IL_name"] for r in ckjoin}),
        "by_experiment": {f'{r}': n for r, n in sorted(
            {(f'{x["_pool"]}/{x["Experiment_ID"]}/{x["IL_name"]}'): sum(
                1 for y in ckjoin
                if f'{y["_pool"]}/{y["Experiment_ID"]}/{y["IL_name"]}'
                == f'{x["_pool"]}/{x["Experiment_ID"]}/{x["IL_name"]}')
                for x in ckjoin}.items())},
        "expected_216_name_241_inchikey": (len(ck) == 216 and len(ckjoin) == 241),
        "PASSES": len(ckjoin) > len(ck),
    }

    # POSITIVE CONTROL #2 -- MC3, documented with two spellings
    mc = [r for r in corpus if "MC3" in r["IL_name"].upper()]
    out["positive_control_2_MC3"] = {
        "name_substring_rows": len(mc),
        "spellings": sorted({r["IL_name"] for r in mc
                             if r["IL_name"].upper().replace("-", "")
                             in ("DLINMC3DMA", "MC3")}),
        "PASSES": len(mc) > 0,
    }

    # --------------------------------------------- the paper's own lipid ---
    comp, amb = load(corpus)
    try:
        il_smiles = smi(comp, amb, IL_NAME)
        il_guard = "resolved (unambiguous in the corpus)"
    except Exception as e:                                   # pragma: no cover
        il_smiles, il_guard = None, "RAISED: %s" % e
    ik, mol = inchikey(il_smiles)
    byname = [r for r in corpus if r["IL_name"] == IL_NAME]
    byik = [r for r in corpus if r["_IK"] == ik]
    out["C12_200_in_corpus"] = {
        "guarded_lookup": il_guard,
        "SMILES": il_smiles,
        "formula": rdm.CalcMolFormula(mol),
        "exact_mass": round(rdm.CalcExactMolWt(mol), 4),
        "n_nitrogen": sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N"),
        "inchikey": ik,
        "name_join_rows": len(byname),
        "inchikey_join_rows": len(byik),
        "names_under_this_structure": sorted({r["IL_name"] for r in byik}),
        "by_pool_experiment": {k: v for k, v in sorted(
            {(f'{r["_pool"]}/{r["Experiment_ID"]}/{r["IL_name"]}'): sum(
                1 for y in byik
                if f'{y["_pool"]}/{y["Experiment_ID"]}/{y["IL_name"]}'
                == f'{r["_pool"]}/{r["Experiment_ID"]}/{r["IL_name"]}')
                for r in byik}.items())},
        "vitro_rows": sum(1 for r in byik if r["_pool"] == "vitro"),
        "vivo_rows": sum(1 for r in byik if r["_pool"] == "vivo"),
        "same_shape_as_positive_control": "yes -- 20 extra rows recovered under the "
                                          "name 'other' in JW_2024, exactly as "
                                          "cKK-E12 hides 19 there",
    }

    # ------------------------------------- the in-vivo shared-lipid question --
    vivo_ik = {}
    for r in vivo:
        if r["_IK"]:
            vivo_ik.setdefault(r["_IK"], set()).add(r["Experiment_ID"])
    shared = {k: sorted(v) for k, v in vivo_ik.items() if len(v) > 1}
    names = {}
    for r in vivo:
        if r["_IK"] in shared:
            names.setdefault(r["_IK"], set()).add(r["IL_name"])
    out["in_vivo_shared_structures"] = {
        "n_structures_shared_across_in_vivo_experiments": len(shared),
        "detail": {k: {"experiments": v, "names": sorted(names[k])}
                   for k, v in shared.items()},
        "C12_200_currently_in_how_many_in_vivo_experiments":
            len(vivo_ik.get(ik, set())),
        "would_KK_2015_add_an_in_vivo_SHARED_lipid": False,
        "why": "C12-200 appears in 0 in-vivo experiments today (518 rows, all "
               "in vitro, AC_2025 + JW_2024). Adding ONE in-vivo experiment takes "
               "it from 0 to 1, and a structure needs >=2 in-vivo experiments to be "
               "shared. So KK_2015 would create the FIRST in-vivo C12-200 rows and "
               "still leave the in-vivo shared set at exactly DLin-MC3-DMA and "
               "SM-102. It does however create a vitro<->vivo bridge for a "
               "518-row structure, which is a different and smaller thing.",
    }

    # --------------------------------------------- token-pattern test (F-055) --
    tok = {}
    for r in corpus:
        if "200" in r["IL_name"]:
            tok.setdefault(r["IL_name"], {"rows": 0, "smiles": set(),
                                          "exps": set()})
            tok[r["IL_name"]]["rows"] += 1
            tok[r["IL_name"]]["smiles"].add(r["IL_SMILES"])
            tok[r["IL_name"]]["exps"].add(r["Experiment_ID"])
    naive = [r for r in corpus if "200" in r["IL_name"]]
    core = Chem.MolFromSmarts("[#7](CCN1CCN(CC[#7])CC1)CC[#7]")
    pent = {}
    for n, d in sorted(tok.items()):
        ms = [Chem.MolFromSmiles(s) for s in d["smiles"]]
        pent[n] = {"rows": d["rows"], "experiments": sorted(d["exps"]),
                   "n_distinct_structures": len(d["smiles"]),
                   "formulas": sorted({rdm.CalcMolFormula(m) for m in ms if m}),
                   "n_nitrogen": sorted({sum(1 for a in m.GetAtoms()
                                             if a.GetSymbol() == "N")
                                         for m in ms if m}),
                   "has_amine200_pentaamine_core":
                       sorted({bool(m.HasSubstructMatch(core)) for m in ms if m})}
    out["token_pattern_test"] = {
        "the_name": "C12-200",
        "pattern": "C<tail length>-<amine catalogue index>. 'C12' is the tail (a "
                   "2-hydroxydodecyl arm from 1,2-epoxydodecane); '200' is the index "
                   "of the polyamine core in this lab's amine catalogue.",
        "which_token_is_which_PROVEN_BY": "structural attestation, not a substring "
            "guess: in KW_2014 the A200_* series holds '200' FIXED while the tail "
            "token runs C10..C14, and all five of those structures contain the same "
            "5-nitrogen pentaamine core; the carbon count is what changes. So '200' "
            "names the core and the C-number names the tail. KK_2015 writes the two "
            "tokens in the OPPOSITE order from KL_2010's '200C12'.",
        "is_the_numbering_contiguous": False,
        "contiguity_note": "F-055: the amine catalogue is sparse and non-contiguous, "
                           "so a contiguous 1..N name would mean LOCAL indices. '200' "
                           "is a catalogue index, not a local one.",
        "naive_substring_damage": {
            "query": "IL_name contains '200'",
            "rows_returned": len(naive),
            "distinct_names": len(tok),
            "distinct_structures": len({r["IL_SMILES"] for r in naive}),
            "verdict": "a naive substring join on '200' returns %d rows spanning %d "
                       "names and %d distinct structures with at least three "
                       "incompatible chemotypes. Never used; the join here is on "
                       "InChIKey." % (len(naive), len(tok),
                                      len({r["IL_SMILES"] for r in naive})),
        },
        "per_name": pent,
    }

    # ---------------------------------------------------- F-061, report only --
    kl = [r for r in corpus if r["IL_name"] in ("200C8", "200C10", "200C12", "200C14")]
    out["F061_known_conflict"] = {
        "status": "REPORTED, NOT REPAIRED. data/ is untouched.",
        "KL_2010_200C_series": {
            "rows": len(kl),
            "names": sorted({r["IL_name"] for r in kl}),
            "formulas": sorted({rdm.CalcMolFormula(Chem.MolFromSmiles(r["IL_SMILES"]))
                                for r in kl}),
            "n_nitrogen": sorted({sum(1 for a in Chem.MolFromSmiles(r["IL_SMILES"]).GetAtoms()
                                      if a.GetSymbol() == "N") for r in kl}),
            "has_pentaamine_core": sorted({bool(Chem.MolFromSmiles(r["IL_SMILES"])
                                                .HasSubstructMatch(core)) for r in kl}),
        },
        "AC_2025_C12_200": {"formula": rdm.CalcMolFormula(mol),
                            "n_nitrogen": 5, "rows": len(byname)},
        "which_form_does_THIS_paper_support": "NEITHER, by its own evidence. "
            "Kauffman 2015 cites C12-200 by name only: the SI contains no structure "
            "drawing, no molecular formula, no HRMS, no molecular weight and no "
            "synthesis. The main text is closed-access and not held. The weight "
            "ratio and the molar percentages are independent recipe parameters, so "
            "no molecular mass is derivable from them either. This paper therefore "
            "CANNOT adjudicate F-061 and is not evidence on either side.",
        "what_the_corpus_itself_says": "Two independent experiments agree the '200' "
            "core is a 5-nitrogen pentaamine bearing 5 arms: AC_2025's C12-200 "
            "(C70H145N5O5, 498 rows) and KW_2014's A200_C10..C14 acrylate series. "
            "KL_2010's 4 rows are the lone dissent at 3 N and 1 arm. That is "
            "corroboration of the majority form but it is NOT this paper's evidence, "
            "and the AC_2025 precedent is exactly why nothing is repaired here.",
    }

    # ------------------------------------- substructure sweep (A SUPERSET) --
    loose = Chem.MolFromSmarts("[#7]-[CH2]-[CH1](-[OX2H1])")   # N-CH2-CH(OH)-
    strict = core
    assert mol.HasSubstructMatch(loose), "POSITIVE CONTROL FAILED: C12-200 does " \
        "not match its own loose SMARTS"
    assert mol.HasSubstructMatch(strict), "POSITIVE CONTROL FAILED: C12-200 does " \
        "not match its own strict SMARTS"
    uniqmols = {s: Chem.MolFromSmiles(s) for s in keys}
    nl = {s for s, m in uniqmols.items() if m and m.HasSubstructMatch(loose)}
    ns = {s for s in nl if uniqmols[s].HasSubstructMatch(strict)}
    ne = {s for s in ns if keys[s] == ik}
    rowsof = lambda ss: sum(1 for r in corpus if r.get("IL_SMILES") in ss)
    out["substructure_sweep"] = {
        "LABEL": "SUPERSET. The loose and strict tiers are chemotype supersets, not "
                 "identity. Only the exact tier is identity.",
        "positive_control": "C12-200 itself matches at the loose, strict and exact "
                            "tiers. Asserted in code; the sweep raises if it does not.",
        "tier_loose": {"smarts": "[#7]-[CH2]-[CH1](-[OX2H1])",
                       "means": "any amine bearing a 2-hydroxyalkyl arm, i.e. the "
                                "epoxide-ring-opening chemotype",
                       "structures": len(nl), "rows": rowsof(nl),
                       "of_total_structures": len(keys)},
        "tier_strict": {"smarts": "[#7](CCN1CCN(CC[#7])CC1)CC[#7]",
                        "means": "the amine-200 pentaamine core AND a hydroxyalkyl arm",
                        "structures": len(ns), "rows": rowsof(ns)},
        "tier_exact": {"means": "full InChIKey identity with C12-200",
                       "structures": len(ne), "rows": rowsof(ne)},
        "experiments_at_strict_tier": sorted({r["Experiment_ID"] for r in corpus
                                              if r.get("IL_SMILES") in ns}),
    }

    # --------------------------------------------------- helper lipids, 3-way --
    hl = {}
    for r in corpus:
        if r.get("HL_name") and r.get("HL_SMILES"):
            hl.setdefault(r["HL_name"].strip(), set()).add(r["HL_SMILES"])
    helpers = {}
    for h in HELPERS + HELPERS_DESIGNED_ONLY:
        s = hl.get(h, set())
        rows = sum(1 for r in corpus if (r.get("HL_name") or "").strip() == h)
        rec = {"corpus_rows": rows, "n_distinct_corpus_SMILES": len(s)}
        if len(s) == 1:
            sm = next(iter(s))
            k, m = inchikey(sm)
            rec.update({"corpus_SMILES": sm, "corpus_inchikey": k,
                        "formula": rdm.CalcMolFormula(m)})
        helpers[h] = rec
    out["helper_lipids_in_corpus"] = {
        "note": "DSPE has 0 corpus rows. It is also never measured in this paper: "
                "all 3 DSPE designs (A-07, A-12, A-14) failed to formulate.",
        "per_helper": helpers,
    }

    # -------------------------------------------------- 03_enumerated.csv ----
    rows = [{
        "component": "ionizable_lipid", "name": IL_NAME, "SMILES": il_smiles,
        "formula": rdm.CalcMolFormula(mol), "inchikey": ik,
        "source_of_structure": "LIVE CORPUS via corpus_components.smi() (guarded; "
                               "unambiguous). NOT from this paper.",
        "drawn_in_this_paper": "no",
        "attested_in_this_paper_by": "name only",
        "confidence": "medium",
        "note": "This paper draws no structure, prints no formula, no HRMS and no "
                "MW for C12-200. The structure is carried over from the corpus, "
                "where AC_2025 (498 rows) and JW_2024 (20 rows) agree. See F-061.",
    }]
    for h in HELPERS:
        s = next(iter(hl[h]))
        k, m = inchikey(s)
        rows.append({
            "component": "helper_lipid", "name": h, "SMILES": s,
            "formula": rdm.CalcMolFormula(m), "inchikey": k,
            "source_of_structure": "provenance/reagent_registry.py (structure written "
                                   "from the chemical name), cross-checked against the "
                                   "corpus HL_SMILES and against this SI's own name "
                                   "decomposition -- see 07_qc.json C11_helper_3way",
            "drawn_in_this_paper": "no",
            "attested_in_this_paper_by": "systematic chemical name, decomposed in the "
                                         "Table S1/S3/S4 footnotes",
            "confidence": "high",
            "note": "catalogue phospholipid; the SI spells out DS/DO and PC/PE",
        })
    rows.append({
        "component": "helper_lipid", "name": "DSPE", "SMILES": "",
        "formula": "", "inchikey": "",
        "source_of_structure": "NOT RESOLVED and not needed",
        "drawn_in_this_paper": "no",
        "attested_in_this_paper_by": "abbreviation only",
        "confidence": "",
        "note": "Designed into A-07, A-12 and A-14 but 'LNP could not be synthesized "
                "due to insolubility of DSPE in ethanol at all concentrations tested'. "
                "No value exists, so no structure is emitted. Also note DSPE alone is "
                "an incomplete specification: a phosphoethanolamine with no stated "
                "head-group derivatisation.",
    })
    rows.append({
        "component": "sterol", "name": "Cholesterol",
        "SMILES": "", "formula": "", "inchikey": "",
        "source_of_structure": "DEFERRED to build_dataset.py / the registry at merge "
                               "time; deliberately NOT hand-written here",
        "drawn_in_this_paper": "no",
        "attested_in_this_paper_by": "the abbreviation 'Chol.' in a column header",
        "confidence": "",
        "note": "Trap 6: a hand-written sterol came out a diastereomer with 4 of 8 "
                "centres inverted, at identical formula, identical exact mass and an "
                "identical skeleton InChIKey block. This paper gives no systematic "
                "name and no drawing, so no sterol structure is asserted by this "
                "curation.",
    })
    rows.append({
        "component": "peg_lipid", "name": "Unknown", "SMILES": "",
        "formula": "", "inchikey": "",
        "source_of_structure": "NONE -- unavailable",
        "drawn_in_this_paper": "no",
        "attested_in_this_paper_by": "nothing: the column header is the bare word 'PEG'",
        "confidence": "",
        "note": "The PEG-lipid is NEVER NAMED in the SI. Emitted as PEG_name='Unknown' "
                "with a blank SMILES, which is an existing corpus convention (1,798 "
                "rows in KZ_2016 and RG_2023). Guessing was specifically refused: the "
                "briefing records C14-PEG2000 drawn as a phosphoethanolamine-carbamate-"
                "PEG but stored as the glycerol ether, so a plausible guess here is a "
                "wrong structure on 36 rows.",
    })

    with open(os.path.join(OUTDIR, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(OUTDIR, "02_corpus_overlap.json"), "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True, default=str)

    print("corpus:", out["corpus_shape"])
    print("KK_2015 collides:", out["experiment_id"]["collides_with_existing_ID"],
          "| PMID already present:", out["experiment_id"]["pmid_26469188_already_in_corpus"])
    print("POSCTRL cKK-E12: name", out["positive_control_1_cKK_E12"]["name_join_rows"],
          "-> IK", out["positive_control_1_cKK_E12"]["inchikey_join_rows"],
          out["positive_control_1_cKK_E12"]["names_recovered"],
          "PASS", out["positive_control_1_cKK_E12"]["PASSES"])
    c = out["C12_200_in_corpus"]
    print("C12-200:", c["formula"], c["inchikey"], "| name", c["name_join_rows"],
          "-> IK", c["inchikey_join_rows"], c["names_under_this_structure"],
          "| vitro", c["vitro_rows"], "vivo", c["vivo_rows"])
    print("in-vivo shared structures:",
          out["in_vivo_shared_structures"]["n_structures_shared_across_in_vivo_experiments"],
          {k: v["names"] for k, v in out["in_vivo_shared_structures"]["detail"].items()})
    print("naive '200' substring:", out["token_pattern_test"]["naive_substring_damage"])
    s = out["substructure_sweep"]
    print("sweep SUPERSET: loose %d structs/%d rows -> strict %d/%d -> exact %d/%d"
          % (s["tier_loose"]["structures"], s["tier_loose"]["rows"],
             s["tier_strict"]["structures"], s["tier_strict"]["rows"],
             s["tier_exact"]["structures"], s["tier_exact"]["rows"]))
    print("strict-tier experiments:", s["experiments_at_strict_tier"])
    for h, d in out["helper_lipids_in_corpus"]["per_helper"].items():
        print("  helper %-5s corpus_rows=%-6s distinct_SMILES=%s"
              % (h, d["corpus_rows"], d["n_distinct_corpus_SMILES"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
