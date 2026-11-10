#!/usr/bin/env python3
"""Stage 03 -- the structures this paper actually determines, and only those.

WHAT IS DETERMINABLE FROM THE MATERIAL HELD. Three compounds, and three only:

  1C12-D  N-dodecyl-N-(2-(pyrrolidin-1-yl)ethyl)dodecan-1-amine      SI p.7
  1C16-S  N-(2-(pyrrolidin-1-yl)ethyl)hexadecan-1-amine              SI p.10
  1C16-D  N-hexadecyl-N-(2-(pyrrolidin-1-yl)ethyl)hexadecan-1-amine  SI p.13

Each is a COLUMN-PURIFIED single compound with a full IUPAC name, a 1H NMR
integral set, a 13C peak list and an HRMS pair. The name is the connectivity
statement; the HRMS and the integral sum are the independent checks on it.

WHAT IS NOT DETERMINABLE. Everything else.
  * The 20 amine heads are integers in Table S1. Only amine 1 is identified.
  * The 12 lipidoids carrying Table S3/S4 values are CRUDE one-pot products,
    i.e. mixtures of mono- through penta-alkylated species (Table S1 shows up
    to 6 tails for a single (amine, alkyl) pair; 23/23 rows show >1). Even the
    three amine-1 crude entries -- 1C12, 1C14, 1C16 -- are mixtures whose
    composition the paper never reports, so they are NOT single structures.

So this stage builds 3 structures, not 200 and not 12. It does not enumerate
the library: enumerating 200 names against 1 known head would produce 200 rows
of which 190 carry no structure at all.

HOW THEY ARE BUILT. curation/lib/assemble.attach, atom-map to atom-map. No
string surgery on a canonical SMILES: a sibling did that and silently
mislocated every ester, with identical formula, exact mass identical to 8
decimal places, and both HRMS anchors still passing at 0.3-0.5 ppm. There is
no stereocentre anywhere in this chemistry (asserted below), which removes the
F-013 class but not the connectivity class -- mass cannot see connectivity, so
the 1H integral sum is carried as a separate witness.

Run:  python3 -I stage03_enumerate.py
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import assemble                                              # noqa: E402
from rdkit import Chem                                       # noqa: E402
from rdkit.Chem import inchi                                 # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt                # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PROTON = 1.007276

# The head, written from the chemical name printed on SI p.3
# ('1-(2-aminoethyl)pyrrolidine'), with the PRIMARY nitrogen atom-mapped :1.
# The ring N is deliberately NOT mapped: it is tertiary and cannot be
# alkylated by this chemistry without quaternising, which the paper does not
# claim and the HRMS would show as a permanent cation.
HEAD_MAPPED = "C1CCN(C1)CC[NH2:1]"

# A 1-iodo-n-alkane's alkyl group, attachment carbon mapped :2. Built as a
# string of 'C's with the first bracketed -- there is no canonical SMILES being
# edited here, the fragment is written from the reagent name.
def alkyl_mapped(n):
    if n < 2:
        raise ValueError("this library's shortest halide is C6")
    return "[CH2:2]" + "C" * (n - 1)


# The three purified compounds, keyed by the paper's own label, each carrying
# the IUPAC name the paper prints and the HRMS pair it prints.
PURIFIED = [
    dict(label="1C12-D", amine=1, chain=12, n_tails=2,
         iupac="N-dodecyl-N-(2-(pyrrolidin-1-yl)ethyl)dodecan-1-amine",
         si_page=7, hrms_page=9,
         hrms_calc_mz=451.4986, hrms_obs_mz=451.4987,
         yield_pct=26.7, chromatography="silica gel, chloroform/methanol 80/1",
         h_nmr_text_total=59, h_nmr_fig_integral_total=63.20),
    dict(label="1C16-S", amine=1, chain=16, n_tails=1,
         iupac="N-(2-(pyrrolidin-1-yl)ethyl)hexadecan-1-amine",
         si_page=10, hrms_page=12,
         hrms_calc_mz=339.3734, hrms_obs_mz=339.3735,
         yield_pct=15.6, chromatography="silica gel, chloroform/methanol 50/1",
         h_nmr_text_total=45, h_nmr_fig_integral_total=44.45),
    dict(label="1C16-D", amine=1, chain=16, n_tails=2,
         iupac="N-hexadecyl-N-(2-(pyrrolidin-1-yl)ethyl)hexadecan-1-amine",
         si_page=13, hrms_page=15,
         hrms_calc_mz=563.62378, hrms_obs_mz=563.62344,
         yield_pct=49.3, chromatography="silica gel, chloroform/methanol 50/1",
         h_nmr_text_total=78, h_nmr_fig_integral_total=75.76),
]

# The secondary, NON-RING amine created by the first alkylation. Used to
# re-tag for the second alkylation. Asserted UNIQUE at the call site: this is
# the only place a nitrogen is located by pattern rather than by being written
# into the fragment, and 'exactly one match or raise' is what makes it safe.
SECONDARY_AMINE = Chem.MolFromSmarts("[NX3;H1;!R]")


def retag(smiles, patt, mapno, label):
    """Set an atom map on the single atom matching `patt`. Raises otherwise."""
    m = Chem.MolFromSmiles(smiles)
    hits = m.GetSubstructMatches(patt)
    if len(hits) != 1:
        raise RuntimeError(
            f"retag({label}): {len(hits)} matches for the secondary amine "
            "pattern, need exactly 1. Refusing to pick one -- 'atom index 0 "
            "of the fragment' is forbidden for the same reason.")
    m.GetAtomWithIdx(hits[0][0]).SetAtomMapNum(mapno)
    return Chem.MolToSmiles(m)


def build(amine_head_mapped, chain, n_tails, label):
    """N-alkylate the primary amine `n_tails` times. Each step is one attach."""
    smi = assemble.attach(amine_head_mapped, alkyl_mapped(chain),
                          a_map=1, b_map=2, bond="single",
                          label=f"{label} tail 1")
    for t in range(2, n_tails + 1):
        smi = assemble.attach(retag(smi, SECONDARY_AMINE, 1,
                                    f"{label} before tail {t}"),
                              alkyl_mapped(chain), a_map=1, b_map=2,
                              bond="single", label=f"{label} tail {t}")
    return smi


def main():
    recs = []
    for p in PURIFIED:
        smi = build(HEAD_MAPPED, p["chain"], p["n_tails"], p["label"])
        m = Chem.MolFromSmiles(smi)

        # --- post-build assertions beyond what assemble already makes ------
        assemble.check_no_radicals(m, label=p["label"])
        sc = assemble.stereocentres(m)
        if sc:
            raise RuntimeError(
                f"{p['label']}: {len(sc)} stereocentre(s) {sc}. This chemistry "
                "is an N-alkylation of an achiral amine with n-alkyl iodides; "
                "a stereocentre here means the build is wrong.")
        # the pyrrolidine N must still be tertiary and in a ring
        ring_n = [a for a in m.GetAtoms()
                  if a.GetSymbol() == "N" and a.IsInRing()]
        if len(ring_n) != 1 or ring_n[0].GetTotalNumHs() != 0 \
                or ring_n[0].GetDegree() != 3 \
                or ring_n[0].GetFormalCharge() != 0:
            raise RuntimeError(
                f"{p['label']}: the pyrrolidine nitrogen is not a neutral "
                "tertiary ring N. Alkylating it would quaternise the head, "
                "which the paper does not claim and HRMS would show as a "
                "permanent cation.")
        # the alkylated N must carry exactly n_tails alkyl chains
        chain_n = [a for a in m.GetAtoms()
                   if a.GetSymbol() == "N" and not a.IsInRing()]
        if len(chain_n) != 1:
            raise RuntimeError(f"{p['label']}: {len(chain_n)} acyclic N.")
        want_h = {1: 1, 2: 0}[p["n_tails"]]
        if chain_n[0].GetTotalNumHs() != want_h:
            raise RuntimeError(
                f"{p['label']}: the alkylated N carries "
                f"{chain_n[0].GetTotalNumHs()} H, expected {want_h} for "
                f"{p['n_tails']} tail(s).")

        neutral = ExactMolWt(m)
        formula = CalcMolFormula(m)
        n_h = sum(a.GetTotalNumHs() for a in m.GetAtoms())
        recs.append(dict(
            IL_name=p["label"],
            IL_SMILES=smi,
            iupac_name_as_printed=p["iupac"],
            amine_head="amine 1 = 1-(2-aminoethyl)pyrrolidine",
            alkyl_tail=f"n-C{p['chain']}",
            n_tails=p["n_tails"],
            molecular_formula=formula,
            our_neutral_exact=round(neutral, 5),
            our_MH_plus=round(neutral + PROTON, 5),
            their_hrms_calc_mz=p["hrms_calc_mz"],
            their_hrms_obs_mz=p["hrms_obs_mz"],
            d_ours_vs_their_calc_mDa=round(
                (neutral + PROTON - p["hrms_calc_mz"]) * 1000, 4),
            d_ours_vs_their_obs_mDa=round(
                (neutral + PROTON - p["hrms_obs_mz"]) * 1000, 4),
            d_their_calc_vs_their_obs_mDa=round(
                (p["hrms_calc_mz"] - p["hrms_obs_mz"]) * 1000, 4),
            ppm_ours_vs_their_calc=round(
                (neutral + PROTON - p["hrms_calc_mz"])
                / p["hrms_calc_mz"] * 1e6, 3),
            formula_H_count=n_h,
            h_nmr_text_total=p["h_nmr_text_total"],
            h_nmr_text_minus_formula=p["h_nmr_text_total"] - n_h,
            h_nmr_fig_integral_total=p["h_nmr_fig_integral_total"],
            inchikey_full=inchi.MolToInchiKey(m),
            inchikey_skeleton=inchi.MolToInchiKey(m).split("-")[0],
            si_page=p["si_page"], hrms_page=p["hrms_page"],
            yield_pct=p["yield_pct"], chromatography=p["chromatography"],
            structure_source="full IUPAC name printed in the SI, corroborated "
                             "by HRMS and the 1H integral sum",
            confidence="high",
            is_a_mixture=0))

    # --- the crude library entries, recorded but NOT given a structure -----
    crude = []
    for lab, chain in [("1C12", 12), ("1C14", 14), ("1C16", 16)]:
        crude.append(dict(
            IL_name=lab, IL_SMILES="", amine_head="amine 1 (KNOWN)",
            alkyl_tail=f"n-C{chain}", n_tails="1-2+ (MIXTURE, not reported)",
            why_no_structure="Crude one-pot product. Table S1 reports BOTH a "
                             f"1-tail and a 2-tail species for amine 1 + C{chain}"
                             " and the paper purified 1C12-D and 1C16-S/D out "
                             "of exactly these mixtures, so the screened "
                             "material is a mixture of at least two compounds "
                             "in an unreported ratio. One SMILES would be a "
                             "fabrication.",
            is_a_mixture=1))
    for lab in ["3C12", "4C16", "4C18", "5C14", "10C14", "11C9", "12C16",
                "14C6", "16C14"]:
        crude.append(dict(
            IL_name=lab, IL_SMILES="",
            amine_head=f"amine {lab.split('C')[0]} (UNKNOWN -- not in this SI)",
            alkyl_tail="C" + lab.split("C")[1],
            n_tails="1-2+ (MIXTURE, not reported)",
            why_no_structure="The amine head is an integer in Table S1 and is "
                             "identified nowhere in this SI; the material is "
                             "also a crude mixture. Two independent reasons, "
                             "either one sufficient. DECISIONS.md D-1.",
            is_a_mixture=1))

    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        cols = list(recs[0].keys())
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in recs:
            w.writerow(r)

    with open(os.path.join(HERE, "03_not_determinable.csv"), "w",
              newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(crude[0].keys()))
        w.writeheader()
        for r in crude:
            w.writerow(r)

    out = dict(
        _what="Stage 03. The three structures this paper determines, built "
              "with curation/lib/assemble atom-map attachment.",
        n_structures_built=len(recs),
        n_labels_with_values_but_no_structure=len(crude),
        head_fragment=HEAD_MAPPED,
        build_route="assemble.attach(head[N:1], alkyl[CH2:2]) per tail; the "
                    "second tail re-tags the unique acyclic secondary amine "
                    "via SMARTS [NX3;H1;!R], asserted to match exactly once",
        no_stereocentres="asserted per structure; this chemistry is achiral, "
                         "so the F-013 chiral-tag class cannot arise -- but "
                         "the connectivity class still can, which is why the "
                         "1H integral sum is carried separately",
        structures=recs,
        not_determinable=crude)
    with open(os.path.join(HERE, "03_enumerated.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print(f"OK  {len(recs)} structures built; {len(crude)} valued labels "
          "left WITHOUT a structure")
    for r in recs:
        print(f"    {r['IL_name']:<7} {r['molecular_formula']:<10} "
              f"{r['IL_SMILES']}")
        print(f"            [M+H]+ ours {r['our_MH_plus']:.5f} vs their calc "
              f"{r['their_hrms_calc_mz']} "
              f"({r['d_ours_vs_their_calc_mDa']:+.3f} mDa, "
              f"{r['ppm_ours_vs_their_calc']:+.3f} ppm); "
              f"InChIKey {r['inchikey_full']}")
        print(f"            1H text total {r['h_nmr_text_total']}H vs formula "
              f"{r['formula_H_count']}H -> "
              f"{r['h_nmr_text_minus_formula']:+d}; figure integral sum "
              f"{r['h_nmr_fig_integral_total']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
