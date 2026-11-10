#!/usr/bin/env python3
"""Stage 02 -- corpus overlap, at three levels, plus the shared-formulation
ratio that trap 4 demands.

RUN: python3 -I scripts/stage02_overlap.py <1080csv> <vitro.csv> <vivo.csv> <out.json>

Three independent keys, because they fail differently:

  1. FULL InChIKey   -- structure with stereo. The only key that sees the
                        epicholesterol class (F-013).
  2. skeleton block  -- connectivity without stereo.
  3. IL_name / HL_name -- the label. CLAUDE.md section 3: 288 IL names carry
                        different structures in different experiments, so a
                        name hit is a lead, never a conclusion.

And a fourth that is specific to this paper: the FORMULATION key. `YZ_2022`
already in the corpus is 1,080 rows, DLin-MC3-DMA only, the same six helper
lipids at 180 rows each, pDNA, HepG2. If those are the same formulations, the
corpus already holds this paper's HepG2 column and merging it would double-count
a published screen. Trap 4 says: if two sources share a formulation, MEASURE THE
RATIO. That is what this does, and it never averages or rescales the two.

READ-ONLY on data/. Nothing here is written outside this curation folder.
"""
import csv
import json
import math
import sys

from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

sys.path.insert(0, __file__.rsplit("/", 2)[0] + "/../lib")

OHE = ["18PG", "14PA", "DOPE", "DSPC", "DOTAP", "DDAB"]

# The ionizable lipid, taken VERBATIM from provenance/reagent_registry.py so
# the structure is the one written from the chemical name, not one typed into
# this script. reagents.smi() is the sanctioned path (lib/reagents.py
# precedence rule 1).
sys.path.insert(0, __file__.rsplit("/", 2)[0] + "/../../provenance")


_KCACHE = {}


def keys(smiles):
    """FULL InChIKey and skeleton block, memoised on the SMILES string.

    Memoised because the corpus holds ~31k rows x 5 component columns and the
    distinct structures number in the tens of thousands, not the hundreds of
    thousands; without the cache this recomputes the same key up to 30,000
    times and the stage does not finish. The cache is keyed on the exact
    string, so two spellings of one molecule are still computed separately --
    it is a speed-up, never a canonicalisation. (CLAUDE.md error class 37: a
    check that canonicalises its input has given up the ability to detect
    non-canonical input.)
    """
    if smiles in _KCACHE:
        return _KCACHE[smiles]
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        out = (None, None)
    else:
        k = Chem.MolToInchiKey(m)
        out = (k, k.split("-")[0])
    _KCACHE[smiles] = out
    return out


def spearman(a, b):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    da = math.sqrt(sum((x - ma) ** 2 for x in ra))
    db = math.sqrt(sum((y - mb) ** 2 for y in rb))
    return num / (da * db) if da and db else float("nan")


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    da = math.sqrt(sum((x - ma) ** 2 for x in a))
    db = math.sqrt(sum((y - mb) ** 2 for y in b))
    return num / (da * db) if da and db else float("nan")


def zscore(v):
    n = len(v)
    m = sum(v) / n
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / n)   # population sd, ddof=0
    return [(x - m) / sd for x in v]


def main(repo_csv, vitro, vivo, out_path):
    import reagents  # noqa: E402   lib/reagents.py

    with open(repo_csv, newline="", encoding="utf-8-sig") as fh:
        tr = list(csv.DictReader(fh))
    corpus = []
    for f, pool in ((vitro, "vitro"), (vivo, "vivo")):
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh):
                r["_pool"] = pool
                corpus.append(r)

    Q = {}

    # ---- our structures ------------------------------------------------
    il_smiles = reagents.smi("Dlin-MC3-DMA")
    ours = {"DLin-MC3-DMA (IL)": il_smiles}
    # helper lipids: five resolve from the registry; 18PG is the conflict and
    # is handled explicitly below, never silently.
    for n in ["DOTAP", "DDAB", "DOPE", "DSPC", "14PA"]:
        ours[n + " (HL)"] = reagents.smi(n)
    ours["18PG as DRAWN in this paper (HL)"] = reagents.smi("DOPG")
    ours["Cholesterol (CHL)"] = reagents.smi("Cholesterol")

    # ================================================== 1-3  structure/name
    hits = {}
    for label, smi in ours.items():
        fk, sk = keys(smi)
        fmatch, smatch = {}, {}
        for r in corpus:
            for pre in ("IL", "HL", "CHL", "PEG", "fifthcomponent"):
                s = (r.get(pre + "_SMILES") or "").strip()
                if not s:
                    continue
                k2, s2 = keys(s)
                if k2 == fk:
                    fmatch.setdefault((r["Experiment_ID"], pre), 0)
                    fmatch[(r["Experiment_ID"], pre)] += 1
                elif s2 == sk:
                    smatch.setdefault((r["Experiment_ID"], pre), 0)
                    smatch[(r["Experiment_ID"], pre)] += 1
        hits[label] = {
            "full_inchikey": fk,
            "skeleton": sk,
            "FULL_key_hits": {"%s/%s" % k: v for k, v in sorted(fmatch.items())},
            "n_rows_full": sum(fmatch.values()),
            "SKELETON_only_hits_stereo_differs":
                {"%s/%s" % k: v for k, v in sorted(smatch.items())},
            "n_rows_skeleton_only": sum(smatch.values()),
        }
    Q["structure_overlap"] = hits

    # name overlap
    names = {}
    for n in ["DLin-MC3-DMA", "Dlin-MC3-DMA"] + OHE + ["Cholesterol"]:
        e = {}
        for r in corpus:
            for pre in ("IL", "HL", "CHL", "PEG"):
                if (r.get(pre + "_name") or "").strip() == n:
                    k = "%s/%s" % (r["Experiment_ID"], pre)
                    e[k] = e.get(k, 0) + 1
        if e:
            names[n] = e
    Q["name_overlap"] = names

    # ---- the 18PG conflict, stated as a measurement ---------------------
    reg18 = reagents.smi("18PG")
    drawn18 = reagents.smi("DOPG")
    f1, s1 = keys(reg18)
    f2, s2 = keys(drawn18)
    rows18 = [r for r in corpus if (r.get("HL_name") or "").strip() == "18PG"]
    Q["reagent_18PG_conflict"] = {
        "corpus_registry_18PG": {"smiles": reg18, "full_inchikey": f1,
                                 "skeleton": s1,
                                 "identity": "SOPG, 18:0/18:1, ONE C=C"},
        "this_paper_drawn_18PG": {"smiles": drawn18, "full_inchikey": f2,
                                  "skeleton": s2,
                                  "identity": "DOPG, 18:1/18:1, TWO C=C"},
        "same_skeleton": s1 == s2,
        "same_full_key": f1 == f2,
        "corpus_rows_named_18PG": len(rows18),
        "corpus_experiments_named_18PG":
            sorted(set(r["Experiment_ID"] for r in rows18)),
        "paper_evidence": [
            "SI Table 2 draws 18PG with a cis C=C in BOTH tails "
            "(renders/p2_x68.png, extracted at its native 300 dpi).",
            "SI Table 3 prints Double_bonds = 2 for 18PG; SOPG has 1.",
            "The repo file carries Double_bonds = 2 on all 180 18PG rows.",
        ],
        "mass_cannot_see_this": ("Wrong -- here mass CAN see it: SOPG is "
                                 "C42H80O10P(-) and DOPG is C42H78O10P(-), two "
                                 "hydrogens apart. What mass could not see is "
                                 "WHICH tail carries the unsaturation. The "
                                 "discriminator used is the drawn structure "
                                 "plus the paper's own printed double-bond "
                                 "count, which agree."),
        "posture": ("REPORTED, NOT ACTED ON. Three independent corpus "
                    "experiments carry 18PG as SOPG on the authority of a "
                    "verbatim Materials sentence in YZ_2022. This paper draws "
                    "DOPG. Both could be right -- they are different purchases "
                    "-- or one curation is wrong. Deciding which half to "
                    "repair is exactly the AC_2025 question, whose first "
                    "confident repair put 419 measured values on the wrong "
                    "molecules. See DECISIONS.md D-2."),
    }

    # ================================================== 4  FORMULATION overlap
    # Derive the four molar percents, then join to YZ_2022 / YZ_2024.
    npmap = {}          # N/P -> IL_to_nucleicacid_massratio seen in the corpus
    yz = {}
    for exp in ("YZ_2022", "YZ_2024"):
        yz[exp] = [r for r in corpus if r["Experiment_ID"] == exp]

    def fkey(hl, il, hlp, chl, peg, mass):
        return (hl, round(il, 1), round(hlp, 1), round(chl, 1),
                round(peg, 1), mass)

    # the corpus mass ratios are 8.48 / 16.96 / 25.44 == 1:2:3, which is the
    # N/P 4:8:12 ladder; the map is DERIVED from that ratio, not assumed.
    masses = sorted(set(float(r["IL_to_nucleicacid_massratio"])
                        for r in yz["YZ_2022"]))
    nps = sorted(set(float(r["NP_ratio"]) for r in tr))
    npmap = dict(zip(nps, masses))
    ladder_ok = (len(masses) == len(nps)
                 and all(abs(masses[i] / masses[0] - nps[i] / nps[0]) < 1e-9
                         for i in range(len(nps))))

    ours_rows = []
    for r in tr:
        a, h, p = (float(r["(IL+HL)"]), float(r["HL_(IL+HL)"]),
                   float(r["PEG_(Chol+PEG)"]))
        hl = [c for c in OHE if r[c].strip() == "1"][0]
        ours_rows.append({
            "row": r,
            "hl": hl,
            "IL": a * (100 - h) / 100.0,
            "HL": a * h / 100.0,
            "CHL": (100 - a) * (100 - p) / 100.0,
            "PEG": (100 - a) * p / 100.0,
            "np": float(r["NP_ratio"]),
        })

    form = {}
    for exp in ("YZ_2022", "YZ_2024"):
        idx = {}
        dupes = 0
        for r in yz[exp]:
            k = fkey(r["HL_name"].strip(), float(r["IL_molratio"]),
                     float(r["HL_molratio"]), float(r["CHL_molratio"]),
                     float(r["PEG_molratio"]),
                     float(r["IL_to_nucleicacid_massratio"]))
            if k in idx:
                dupes += 1
            idx[k] = r
        matched, unmatched = [], 0
        for o in ours_rows:
            k = fkey(o["hl"], o["IL"], o["HL"], o["CHL"], o["PEG"],
                     npmap[o["np"]])
            if k in idx:
                matched.append((o, idx[k]))
            else:
                unmatched += 1
        form[exp] = {
            "corpus_rows": len(yz[exp]),
            "corpus_keys_distinct": len(idx),
            "corpus_key_collisions": dupes,
            "our_rows_matched": len(matched),
            "our_rows_unmatched": unmatched,
            "match_is_one_to_one":
                len(matched) == len(tr) and dupes == 0 and len(idx) == len(tr),
        }
        if matched:
            # ---- the ratio, per cell column
            cells = [("RLU_HepG2", "HepG2"), ("RLU_N2a", "N2a"),
                     ("RLU_ARPE19", "ARPE19"), ("RLU_HEK293", "HEK293"),
                     ("RLU_B16", "B16"), ("RLU_PC3", "PC3")]
            per = {}
            for col, nm in cells:
                pairs = [(float(o["row"][col]), float(c["Experiment_value"]))
                         for o, c in matched if o["row"][col].strip() != ""]
                if len(pairs) < 10:
                    continue
                a = [x for x, _ in pairs]
                b = [y for _, y in pairs]
                az = zscore(a)
                # the corpus value is a per-experiment z of log10 signal; a z
                # of ln is IDENTICALLY a z of log10, so this is a like-for-like
                # comparison with no rescaling of either side.
                resid = [abs(x - y) for x, y in zip(az, b)]
                per[nm] = {
                    "n": len(pairs),
                    "spearman_rank": round(spearman(a, b), 6),
                    "pearson_on_z": round(pearson(az, b), 6),
                    "max_abs_residual_after_z": round(max(resid), 6),
                    "mean_abs_residual_after_z":
                        round(sum(resid) / len(resid), 6),
                }
            form[exp]["value_agreement"] = per
    Q["formulation_overlap"] = {
        "np_to_massratio_map": {str(k): v for k, v in npmap.items()},
        "map_derived_from_the_1_2_3_ladder": ladder_ok,
        "per_experiment": form,
        "note": ("Trap 4: where two sources share a formulation, measure the "
                 "ratio. Neither side is averaged or rescaled onto the other; "
                 "the corpus side is a per-experiment z of log10 and our side "
                 "is z-scored ln, and a z of ln IS a z of log10, so the "
                 "comparison is exact rather than fitted."),
    }

    with open(out_path, "w") as fh:
        json.dump(Q, fh, indent=1)
    print(json.dumps({
        "formulation_overlap": {k: {kk: vv for kk, vv in v.items()
                                    if kk != "value_agreement"}
                                for k, v in form.items()},
        "value_agreement": {k: v.get("value_agreement") for k, v in form.items()},
        "18PG_same_full_key": Q["reagent_18PG_conflict"]["same_full_key"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:5]))
