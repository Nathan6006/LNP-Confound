#!/usr/bin/env python3
"""Stage 02 -- corpus overlap, name-token traps, and the one free external
check this paper does get.

WHAT CAN AND CANNOT BE COMPARED. Zero lipidoids are enumerable (stage 03), so
there is nothing to match on the key that matters -- a lipidoid's own full
InChIKey. What IS comparable:

  (a) the two tail REAGENTS, by full InChIKey and by skeleton block;
  (b) this paper's name tokens ('1'..'16', '16I', '16II') against every
      IL_name in the corpus;
  (c) a substructure sweep, reported as a SUPERSET, with a STRICT re-check --
      the loose version of probe (c) produced 17 + 2 false positives by
      matching a C18:1 chain read from the methyl end, which is the documented
      hazard (a sibling's first SMARTS returned 332 rows under a label
      claiming exactness when the exact number was 128);
  (d) THE FREE EXTERNAL CHECK: AA_2008's `*_N18` series is 44 bis-Michael
      adducts of N-OCTADECYLACRYLAMIDE -- independently curated by someone
      else, from a different paper -- which attests this paper's saturated
      tail reagent structure AND, by its amine numbering, DISPROVES the
      single most tempting way to invent 16 structures here.

Run:  python3 -I stage02_corpus_overlap.py
"""
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

import pandas as pd                                        # noqa: E402
from rdkit import Chem, RDLogger                           # noqa: E402

RDLogger.DisableLog("rdApp.*")

NAME_TOKENS = [str(i) for i in range(1, 17)] + ["16I", "16II"]

# What stage03 built, independently of anything in this file.
OLEYL_ACRYLAMIDE = "C=CC(=O)NCCCCCCCC/C=C\\CCCCCCCC"
OCTADECYL_ACRYLAMIDE = "C=CC(=O)NCCCCCCCCCCCCCCCCCC"


def canon(s):
    m = Chem.MolFromSmiles(s)
    return None if m is None else Chem.MolToSmiles(m)


def walk_linear_chain(mol, start_idx, from_idx):
    """Walk an UNBRANCHED carbon chain from `start_idx`, having come from
    `from_idx`. Returns (n_carbons, n_double_bonds, pos_of_first_double_bond,
    ok) where ok is False if the chain branches, hits a ring or a heteroatom,
    or does not end in a CH3.

    THIS IS THE STRICT VERSION. A SMARTS like `NCCCCCCCC/C=C\\CCCCCCCC` has
    no terminal-methyl and no degree constraint, so it matches a C18:1 chain
    read from the METHYL end into a substituted interior carbon. That gave 17
    YX_2026_SDL and 2 ZH_2023 false positives whose amide nitrogen carries
    tert-butyl, not oleyl.
    """
    n, prev, cur, ndb, dbpos = 0, from_idx, start_idx, 0, None
    while True:
        at = mol.GetAtomWithIdx(cur)
        if at.GetSymbol() != "C" or at.IsInRing():
            return (n, ndb, dbpos, False)
        n += 1
        heavy = [a.GetIdx() for a in at.GetNeighbors()
                 if a.GetIdx() != prev]
        if any(mol.GetAtomWithIdx(i).GetSymbol() != "C" for i in heavy):
            return (n, ndb, dbpos, False)
        for b in at.GetBonds():
            o = b.GetOtherAtomIdx(cur)
            if o != prev and b.GetBondType() == Chem.BondType.DOUBLE:
                ndb += 1
                if dbpos is None:
                    dbpos = n
        if len(heavy) == 0:                      # terminal CH3
            return (n, ndb, dbpos, at.GetTotalNumHs() == 3)
        if len(heavy) > 1:                       # branch
            return (n, ndb, dbpos, False)
        prev, cur = cur, heavy[0]


def strict_n_alkyl_amide_arms(smiles):
    """Every arm of the form  [N(tert, not amide)]-CH2-CH2-C(=O)-NH-[linear
    alkyl], i.e. a Michael adduct of an N-alkyl ACRYLAMIDE. Returns a list of
    dicts describing each arm's liberated acrylamide."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return []
    patt = Chem.MolFromSmarts(
        "[NX3;!$([NX3][CX3]=[OX1])]-[CH2]-[CH2]-[CX3](=[OX1])-[NX3H1]")
    out = []
    for match in m.GetSubstructMatches(patt, uniquify=True):
        n_basic, c1, c2, c_carbonyl, o_carbonyl, n_amide = match
        # the amide N's other substituent is the candidate alkyl chain
        others = [a.GetIdx() for a in m.GetAtomWithIdx(n_amide).GetNeighbors()
                  if a.GetIdx() != c_carbonyl]
        if len(others) != 1:
            continue
        nC, ndb, dbpos, ok = walk_linear_chain(m, others[0], n_amide)
        if not ok:
            continue
        # the liberated acrylamide, rebuilt from the MEASURED chain, not assumed
        if ndb == 0:
            chain = "C" * nC
            lib = f"C=CC(=O)N{chain}"
        elif ndb == 1:
            pre = "C" * (dbpos - 1)
            post = "C" * (nC - dbpos - 1)
            lib = f"C=CC(=O)N{pre}/C=C\\{post}"
        else:
            continue
        out.append({"chain_C": nC, "n_db": ndb, "db_pos_from_N": dbpos,
                    "liberated_acrylamide": canon(lib)})
    return out


# Seeds for the two external checks of section (d). Applied in their own
# process by stage07_qc.py.
SEEDS = {
    # the attested tail reagent, one carbon short. d1's agreement count must
    # collapse: if it does not, d1 was not comparing what it claims to.
    "seed_tail_structure": "octadecyl_17",
    # pretend this paper's labels were the 16 LOWEST catalogue numbers. The
    # catalogue-number DISPROOF must then stop holding -- a disproof that
    # holds whatever the tokens are is not a disproof.
    "seed_tokens_are_catalogue": "catalogue_tokens",
}


def main(seed=None):
    global OCTADECYL_ACRYLAMIDE
    enum = os.path.join(HERE, "03_enumerated.csv")
    if not os.path.exists(enum):
        print("run stage03_enumerate.py first", file=sys.stderr)
        return 1
    with open(enum) as fh:
        reagents = list(csv.DictReader(fh))

    if seed == "seed_tail_structure":
        OCTADECYL_ACRYLAMIDE = "C=CC(=O)NCCCCCCCCCCCCCCCCC"   # C17
    corpus = {p: pd.read_csv(os.path.join(ROOT, "data", f"data_{p}.csv"),
                             low_memory=False) for p in ("vitro", "vivo")}
    res = {"preamble": {
        "lipidoids_enumerated_by_this_curation": 0,
        "so_the_key_that_matters": "a lipidoid's own full InChIKey -- there is "
            "none to compare. Everything below concerns the tail reagents, "
            "names, or fragments.",
        "corpus_shape": {k: list(v.shape) for k, v in corpus.items()},
    }}
    fails = []

    # ---- (a) tail reagents ------------------------------------------------
    keys = {}
    for pool, df in corpus.items():
        d = {}
        for s in df["IL_SMILES"].dropna().astype(str).unique():
            m = Chem.MolFromSmiles(s)
            if m is not None:
                d.setdefault(Chem.MolToInchiKey(m), []).append(s)
        keys[pool] = d
    all_keys = {}
    for pool, d in keys.items():
        for k in d:
            all_keys.setdefault(k, set()).add(pool)
    reagent_hits = []
    for r in reagents:
        ik, skel = r["inchikey"], r["inchikey"].split("-")[0]
        reagent_hits.append({
            "reagent": r["role"], "inchikey": ik,
            "full_inchikey_match_in_corpus": sorted(all_keys.get(ik, [])) or None,
            "skeleton_block_matches_in_corpus":
                sorted({k for k in all_keys if k.split("-")[0] == skel}) or None,
        })
    res["a_tail_reagents"] = {
        "unique_parseable_corpus_IL_SMILES": {p: len(d) for p, d in keys.items()},
        "hits": reagent_hits,
        "interpretation": "0 matches, and 0 is the EXPECTED answer: a tail "
            "reagent is a monomer, not an ionizable lipid. A hit would have "
            "meant the corpus stored a reagent as a library member.",
    }

    # ---- (b) name tokens --------------------------------------------------
    name_rows = {}
    for pool, df in corpus.items():
        for nm in df["IL_name"].dropna().astype(str):
            name_rows.setdefault(nm.strip(), set()).add(pool)
    exact = {t: sorted(name_rows.get(t, [])) for t in NAME_TOKENS}
    leading = {t: sum(1 for n in name_rows if n.startswith(t))
               for t in NAME_TOKENS}
    n_leading = len({n for t in NAME_TOKENS for n in name_rows
                     if n.startswith(t)})
    res["b_name_tokens"] = {
        "distinct_IL_names_in_corpus": len(name_rows),
        "exact_matches": {k: v for k, v in exact.items() if v},
        "n_exact_matches": sum(1 for v in exact.values() if v),
        "distinct_corpus_names_merely_STARTING_with_one_of_our_tokens": n_leading,
        "per_token_leading_counts": leading,
        "interpretation": "A LEADING NUMBER TOKEN IS NOT AN IDENTITY. The "
            f"{n_leading} leading-prefix names are reported to SIZE the trap, "
            "never as candidates. Not acted on.",
    }

    # ---- (c) substructure sweep, loose then STRICT ------------------------
    probes = {
        "loose_oleyl_amide_chain": "C(=O)NCCCCCCCC/C=C\\CCCCCCCC",
        "loose_octadecyl_amide_chain": "C(=O)NCCCCCCCCCCCCCCCCCC",
        "loose_michael_adduct_of_oleyl_acrylamide":
            "[N;!$(N-C=O)]CCC(=O)NCCCCCCCC/C=C\\CCCCCCCC",
        "loose_michael_adduct_of_octadecyl_acrylamide":
            "[N;!$(N-C=O)]CCC(=O)NCCCCCCCCCCCCCCCCCC",
    }
    sweep = {}
    for label, smarts in probes.items():
        q = Chem.MolFromSmarts(smarts)
        if q is None:
            raise ValueError(f"bad SMARTS {label}")
        per = {}
        for pool, df in corpus.items():
            sub = df.dropna(subset=["IL_SMILES"])
            hit = set()
            for s in sub.IL_SMILES.astype(str).unique():
                m = Chem.MolFromSmiles(s)
                if m is not None and m.HasSubstructMatch(q):
                    hit.add(s)
            rows = sub[sub.IL_SMILES.astype(str).isin(hit)]
            per[pool] = {"unique_SMILES": len(hit), "rows": len(rows),
                         "experiments": sorted(rows.Experiment_ID.unique().tolist())}
        sweep[label] = {"smarts": smarts, "per_pool": per}

    # STRICT re-check of the two Michael probes
    strict = {}
    for tgt_label, tgt_smiles in (("oleyl", OLEYL_ACRYLAMIDE),
                                  ("octadecyl", OCTADECYL_ACRYLAMIDE)):
        tgt = canon(tgt_smiles)
        per = {}
        for pool, df in corpus.items():
            sub = df.dropna(subset=["IL_SMILES"])
            hit, byexp = set(), {}
            for s in sub.IL_SMILES.astype(str).unique():
                arms = strict_n_alkyl_amide_arms(s)
                if any(a["liberated_acrylamide"] == tgt for a in arms):
                    hit.add(s)
            rows = sub[sub.IL_SMILES.astype(str).isin(hit)]
            for e, g in rows.groupby("Experiment_ID"):
                byexp[e] = {"rows": len(g), "unique_SMILES": g.IL_SMILES.nunique()}
            per[pool] = {"unique_SMILES": len(hit), "rows": len(rows),
                         "by_experiment": byexp}
        strict[tgt_label] = per

    loose_ol = sum(sweep["loose_michael_adduct_of_oleyl_acrylamide"]
                   ["per_pool"][p]["unique_SMILES"] for p in corpus)
    strict_ol = sum(strict["oleyl"][p]["unique_SMILES"] for p in corpus)
    loose_oc = sum(sweep["loose_michael_adduct_of_octadecyl_acrylamide"]
                   ["per_pool"][p]["unique_SMILES"] for p in corpus)
    strict_oc = sum(strict["octadecyl"][p]["unique_SMILES"] for p in corpus)

    res["c_substructure_sweep"] = {
        "THIS_IS_A_SUPERSET": "the `loose_*` counts are counts of corpus "
            "lipids CONTAINING the queried fragment. They are NOT matches to "
            "any molecule in this paper and must not be quoted as such.",
        "loose_probes": sweep,
        "strict_recheck": {
            "what_strict_means": "the arm must be [tertiary non-amide N]-CH2-"
                "CH2-C(=O)-NH-[chain], and the chain must be walked atom by "
                "atom and found UNBRANCHED, ring-free, all-carbon and "
                "terminated by a CH3, with its length and double-bond "
                "position MEASURED. The liberated acrylamide is then rebuilt "
                "from the measured chain and compared as canonical SMILES.",
            "results": strict,
            "loose_vs_strict": {
                "oleyl_michael_adduct": {"loose": loose_ol, "strict": strict_ol,
                                         "false_positives": loose_ol - strict_ol},
                "octadecyl_michael_adduct": {"loose": loose_oc,
                                             "strict": strict_oc,
                                             "false_positives": loose_oc - strict_oc},
            },
            "why_the_oleyl_loose_probe_lied": "it matched a C18:1 chain read "
                "from the METHYL end into a SUBSTITUTED interior carbon. In "
                "YX_2026_SDL's A8+B11+C15+D3 the amide nitrogen carries "
                "tert-butyl, not oleyl. The loose SMARTS constrains neither "
                "the terminal methyl nor the interior degree.",
        },
    }

    # ---- (d) THE FREE EXTERNAL CHECK: AA_2008's N18 series ---------------
    vit = corpus["vitro"]
    aa = vit[vit.Experiment_ID == "AA_2008"]
    fam = {}
    for nm in aa.IL_name.dropna().astype(str).unique():
        mt = re.match(r"^(\d+)_([NO])(\d+)$", nm)
        if mt:
            fam.setdefault(mt.group(2) + mt.group(3), set()).add(int(mt.group(1)))
    n18_amines = sorted(fam.get("N18", []))
    all_aa_amines = sorted(set().union(*fam.values())) if fam else []

    # d1: do AA_2008's N18 lipids liberate EXACTLY our octadecyl acrylamide?
    tgt_oc = canon(OCTADECYL_ACRYLAMIDE)
    n18 = aa[aa.IL_name.astype(str).str.match(r"^\d+_N18$")]
    agree, disagree, armless = 0, [], 0
    for s in n18.IL_SMILES.astype(str).unique():
        arms = strict_n_alkyl_amide_arms(s)
        libs = {a["liberated_acrylamide"] for a in arms}
        if not libs:
            armless += 1
        elif tgt_oc in libs:
            agree += 1
        else:
            disagree.append({"smiles": s, "liberated": sorted(libs)})

    # d2: are this paper's tokens 1..16 a subset of the catalogue?
    ours = set(range(1, 17))
    if seed == "seed_tokens_are_catalogue":
        ours = set(sorted(all_aa_amines)[:16])
    in_cat = sorted(ours & set(all_aa_amines))
    not_in_cat = sorted(ours - set(all_aa_amines))

    res["d_free_external_check"] = {
        "what_it_is": "AA_2008 (Akinc/Anderson 2008 lipidoid library), curated "
            "into this corpus independently of this folder and from a "
            "different paper, contains tail families named `<amine>_<N|O><n>`. "
            "The `N` families are ACRYLAMIDES and the `O` families ACRYLATES; "
            "`n` is the alkyl carbon count. `N18` is therefore "
            "N-octadecylacrylamide -- this paper's own saturated comparator.",
        "d1_tail_structure_attestation": {
            "AA_2008_N18_unique_lipids_tested": int(n18.IL_SMILES.nunique()),
            "liberate_our_octadecyl_acrylamide_exactly": agree,
            "liberate_something_else": len(disagree),
            "no_strict_arm_found": armless,
            "disagreements": disagree[:5],
            "our_independently_built_SMILES": tgt_oc,
            "what_this_attests": "The structure stage03 built from the "
                "chemical name alone is byte-identical, as canonical SMILES, "
                "to the acrylamide liberated by retro-Michael from "
                f"{agree} AA_2008 N18 lipids that another curator entered from "
                "another paper. This UPGRADES 01_blocks.csv's "
                "`octadecyl_acrylamide` evidence: the abstract says only "
                "'n-octadecyl tails' and never prints 'octadecyl acrylamide', "
                "but the corpus shows the amide is this chemistry's standing "
                "convention and that the amide tail reagent is exactly this "
                "molecule. It does NOT prove this paper used the amide -- only "
                "that the structure we wrote for that name is the right one.",
            "what_it_does_NOT_attest": "the OLEYL reagent. The corpus contains "
                f"{strict_ol} strict Michael adducts of oleyl acrylamide -- "
                "i.e. none -- so the unsaturated tail, which is this paper's "
                "whole point, has no corpus attestation at all.",
        },
        "d2_catalogue_number_DISPROOF": {
            "AA_2008_amine_numbers_in_use": all_aa_amines,
            "n_catalogue_amines": len(all_aa_amines),
            "AA_2008_N18_amine_numbers": n18_amines,
            "this_papers_tokens": sorted(ours),
            "tokens_that_ARE_catalogue_numbers": in_cat,
            "tokens_that_are_NOT_catalogue_numbers": not_in_cat,
            "the_argument": "The Anderson/Akinc catalogue numbers actually in "
                "use are NOT CONTIGUOUS -- they run 1, 6, 7, 10, 11, 13, 15, "
                "17, 20, ... skipping 2, 3, 4, 5, 8, 9, 12, 14 and 16. This "
                "paper's Figure S3 plots a CONTIGUOUS run 1..16, and "
                f"{len(not_in_cat)} of those 16 ({not_in_cat}) are not "
                "catalogue numbers at all. Under the hypothesis 'these labels "
                "are catalogue amine numbers' all 16 would have to be; "
                f"{len(not_in_cat)} are not, so the hypothesis is REJECTED.",
            "why_this_matters": "Reading 1..16 as catalogue numbers is the "
                "single most tempting way to invent 16 structures for this "
                "paper, and it is available: AA_2008 stores a structure for "
                f"each of {len(in_cat)} of our tokens. Doing it would have "
                "assigned up to 16 measured values to wrong molecules -- the "
                "AC_2025 failure mode. This is a MEASURED disproof, not a "
                "stylistic preference, and it is the only external check this "
                "paper's name tokens admit.",
            "corroborating_siblings": "curation/xuq_24574196 states the same "
                "for its own Xu-lab paper: \"The paper's local label '1' is "
                "NOT catalogue amine 1.\" The catalogue numbers recorded in "
                "xuq_32520527 / xuq_32438474 / xuq_32662132 are all >= 75.",
        },
    }
    if armless or disagree:
        fails.append(f"d1: {armless} N18 lipids with no strict arm, "
                     f"{len(disagree)} liberating something else")
    if not not_in_cat:
        fails.append("d2: every token 1..16 is a catalogue number -- the "
                     "disproof does not hold and the labels must be "
                     "re-examined before anything is claimed about them")

    res["verdict"] = {
        "lipidoid_level_overlap": "NONE POSSIBLE -- nothing enumerated.",
        "tail_reagent_overlap": "0 full-InChIKey, 0 skeleton-block.",
        "free_external_check_obtained": True,
        "what_it_gave": [
            f"the octadecyl acrylamide structure attested by {agree} "
            "independently curated AA_2008 N18 lipids (retro-Michael, "
            "canonical-SMILES identity)",
            f"a measured DISPROOF that this paper's 1..16 are catalogue amine "
            f"numbers ({len(not_in_cat)} of 16 are not catalogue members)",
        ],
        "what_it_did_NOT_give": "any structure for this paper's 16 amine "
            "heads, and no attestation of the oleyl reagent.",
    }
    res["checks_failed"] = fails

    n_agree = agree + (1 if not_in_cat else 0)
    n_total = int(n18.IL_SMILES.nunique()) + 1
    res["seeded_check_counts"] = {"n_agree": n_agree, "n_total": n_total}
    if seed is None:
        with open(os.path.join(HERE, "02_corpus_overlap.json"), "w") as fh:
            json.dump(res, fh, indent=2)
    print(f"n_agree={n_agree} n_total={n_total}")

    print("corpus:", {k: list(v.shape) for k, v in corpus.items()})
    print(f"(a) tail reagents: full-InChIKey matches "
          f"{sum(1 for h in reagent_hits if h['full_inchikey_match_in_corpus'])}"
          f", skeleton matches "
          f"{sum(1 for h in reagent_hits if h['skeleton_block_matches_in_corpus'])}")
    print(f"(b) IL_name exact matches on our tokens: "
          f"{res['b_name_tokens']['n_exact_matches']}; corpus names merely "
          f"STARTING with one of them: {n_leading} (SUPERSET)")
    print(f"(c) oleyl Michael adduct   loose {loose_ol} -> STRICT {strict_ol} "
          f"({loose_ol - strict_ol} false positives)")
    print(f"    octadecyl Michael adduct loose {loose_oc} -> STRICT {strict_oc} "
          f"({loose_oc - strict_oc} false positives)")
    print(f"(d) AA_2008 N18: {agree}/{n18.IL_SMILES.nunique()} lipids liberate "
          f"EXACTLY our octadecyl acrylamide ({armless} armless, "
          f"{len(disagree)} other)")
    print(f"    catalogue amines in use: {len(all_aa_amines)}, non-contiguous")
    print(f"    our tokens {min(ours)}..{max(ours)}: {len(in_cat)} are "
          f"catalogue numbers, {len(not_in_cat)} are NOT -> catalogue "
          f"hypothesis "
          f"{'REJECTED' if not_in_cat else 'NOT REJECTED'} ({not_in_cat})")
    if fails:
        print("FAILED:", *fails, sep="\n  - ")
        return 1
    return 0


if __name__ == "__main__":
    _s = None
    if "--seed" in sys.argv:
        _s = sys.argv[sys.argv.index("--seed") + 1]
    sys.exit(main(seed=_s))
