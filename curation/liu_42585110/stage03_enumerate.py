#!/usr/bin/env python3
"""Enumerate the 40 RAILs and check every one against the SI's own ESI-MS.

Run from this folder:   python3 -I stage03_enumerate.py   -> 03_enumerated.csv

THE ROUTE, in the order `assemble` is willing to run it:

  1. esterify(retinoid acid, linker-amine-alcohol)
         The linker stands in for the paper's alkyl bromide as the
         corresponding primary amine (stage02's docstring says why). It has
         exactly one free OH, so `esterify` knows which oxygen to use and does
         not have to be told -- and the retinoid goes on BEFORE the head arm,
         because once the head arm is present the molecule has two free
         hydroxyls and `esterify` refuses, correctly.
  2. attach(that amine, the ATX alkyl fragment)   -> the secondary amine
  3. attach(that secondary amine, the 4-hydroxybutyl arm) -> the tertiary amine

No RWMol, no RemoveAtom, no string concatenation. The only place this module
touches an atom directly is `_map_the_only_nitrogen`, which sets an atom map
and refuses if there is not exactly one nitrogen; it forms and breaks nothing.

WHAT THE MASS CHECK CAN AND CANNOT SEE. 39 of the 40 lipids carry a printed
ESI-MS [M+H]+ in the SI (AT3-Bex is the one that does not). Mass confirms the
molecular formula and therefore the arm COUNTS. It is blind to:
  * which nitrogen an arm sits on -- settled here by the synthesis being a
    double alkylation of ONE primary amine, drawn as a product in Schemes
    S5-S12;
  * the direction of an ester -- AT4 reverses it relative to AT1-AT3 and the
    formula is unchanged, so the InChIKey column below is what distinguishes
    them;
  * double-bond geometry -- ATRA, 13CRA and 9CRA are isobaric to 0.000000 Da.
    stage02 asserts their E/Z pattern AND the position of the Z bond, and the
    distinct-InChIKey assertion below is what makes the three stay distinct
    through 8 assembly steps each.
"""
import csv
import json
import os
import re
import sys

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
import assemble                                                  # noqa: E402

from rdkit import Chem                                           # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors             # noqa: E402
from rdkit import RDLogger                                       # noqa: E402
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from stage02_blocks import RETINOIDS, T_ALKYL, LINKER, HEAD_ARM  # noqa: E402

OUT = os.path.join(HERE, "03_enumerated.csv")
SI_TEXT = os.path.join(HERE, "00_pages", "si_text.txt")

PROTON = 1.00727646          # mass of H+ (H atom minus one electron)
H_ATOM = 1.00782503          # mass of a neutral H atom
NA_CATION = 22.98922         # Na+ : 22.98977 - 0.000549

# The panel order of Figure 1a/1c/1d, left to right. Read off the rotated
# axis labels of Figure 1a at native resolution (see 04_fig1a_heatmap.py).
RETINOID_ORDER = ["ATRA", "9CRA", "13CRA", "AVA", "Bex"]
LINKER_ORDER = ["c", ""]
SCAFFOLD_ORDER = ["AT1", "AT2", "AT3", "AT4"]


def panel_order():
    """The 40 library member names, in Figure 1's left-to-right order."""
    names = []
    for s in SCAFFOLD_ORDER:
        for lk in LINKER_ORDER:
            for r in RETINOID_ORDER:
                names.append(f"{s}-c-{r}" if lk == "c" else f"{s}-{r}")
    return names


def _map_the_only_nitrogen(smiles, mapno=1):
    """Return `smiles` with atom map `mapno` on its single nitrogen."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise SystemExit(f"stage03: will not parse: {smiles}")
    ns = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "N"]
    if len(ns) != 1:
        raise SystemExit(
            f"stage03: {len(ns)} nitrogens, need exactly 1, in {smiles}. "
            "Mapping 'the nitrogen' when there are two is how an arm lands on "
            "the wrong one, and mass cannot see it.")
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    m.GetAtomWithIdx(ns[0]).SetAtomMapNum(mapno)
    return Chem.MolToSmiles(m)


def build(scaffold, linker, retinoid):
    """One RAIL, by esterify -> attach -> attach."""
    ra_smiles = RETINOIDS[retinoid][0]
    label = f"{scaffold}-c-{retinoid}" if linker == "c" else \
            f"{scaffold}-{retinoid}"

    arm = assemble.esterify(ra_smiles, LINKER[linker][0],
                            label=f"{label}: retinoid onto linker")
    arm = _map_the_only_nitrogen(arm, 1)
    sec = assemble.attach(arm, T_ALKYL[scaffold][0], a_map=1, b_map=2,
                          label=f"{label}: first N-alkylation (ATX arm)")
    sec = _map_the_only_nitrogen(sec, 1)
    full = assemble.attach(sec, HEAD_ARM[0], a_map=1, b_map=2,
                           label=f"{label}: second N-alkylation (head arm)")
    assemble.check_no_radicals(full, label=label)
    return assemble.canonical(full)


def printed_hrms():
    """{compound: (adduct, calcd, found)} parsed out of the SI text."""
    if not os.path.exists(SI_TEXT):
        raise SystemExit(
            f"stage03: {SI_TEXT} is missing. Run stage01_si_text.py first; "
            "the HRMS table is parsed from the SI, never transcribed by hand.")
    t = open(SI_TEXT, encoding="utf-8").read()
    out = {}
    for m in re.finditer(r"Synthesis of ([A-Za-z0-9\-]+)\.(.*?)"
                         r"(?=Synthesis of |\Z)", t, re.S):
        hit = re.search(r"Calcd\.?\s*for\s*\[M([+][A-Za-z]+)\][+],?\s*"
                        r"([0-9]+\.[0-9]+),\s*Found,?\s*([0-9]+\.[0-9]+)",
                        m.group(2))
        if hit:
            out[m.group(1)] = (hit.group(1), float(hit.group(2)),
                               float(hit.group(3)))
    return out


def main():
    hrms = printed_hrms()
    rows, keys = [], {}
    # convention test: is the SI's "Calcd. for [M+H]+" a protonated mass or a
    # neutral-plus-H-atom mass? Determined from the data, with the margin
    # printed, rather than assumed. The two differ by one electron, 0.00055 Da.
    resid = {"proton": [], "h_atom": []}

    for name in panel_order():
        parts = name.split("-")
        scaffold = parts[0]
        linker = "c" if len(parts) == 3 else ""
        retinoid = parts[-1]
        smi = build(scaffold, linker, retinoid)
        m = Chem.MolFromSmiles(smi)
        formula = rdMolDescriptors.CalcMolFormula(m)
        neutral = Descriptors.ExactMolWt(m)
        key = Chem.MolToInchiKey(m)
        if key in keys:
            raise SystemExit(
                f"stage03: {name} and {keys[key]} share InChIKey {key}. The "
                "library would not be 40 distinct structures. Two geometric "
                "isomers collapsing here means a stereo marker was lost in "
                "assembly -- which mass cannot see.")
        keys[key] = name

        ad, calcd, found = hrms.get(name, ("", None, None))
        d_proton = d_hatom = None
        if calcd is not None:
            if ad != "+H":
                raise SystemExit(f"stage03: {name} adduct {ad}, expected +H.")
            d_proton = (neutral + PROTON) - calcd
            d_hatom = (neutral + H_ATOM) - calcd
            resid["proton"].append(abs(d_proton))
            resid["h_atom"].append(abs(d_hatom))

        rows.append({
            "IL_name": name,
            "IL_SMILES": smi,
            "inchikey": key,
            "formula": formula,
            "exact_mass_neutral": f"{neutral:.6f}",
            "atx_tail_motif": scaffold,
            "linker": linker or "(none)",
            "retinoid": retinoid,
            "si_adduct": ad,
            "si_calcd_mz": "" if calcd is None else f"{calcd:.3f}",
            "si_found_mz": "" if found is None else f"{found:.3f}",
            "our_mz_minus_si_calcd_mDa": "" if d_proton is None
                                         else f"{d_proton*1000:.2f}",
            "confidence": "high",
            "confidence_reason":
                ("structure built from the synthesis the SI describes in prose "
                 "AND draws as a product in Schemes S5-S12; molecular formula "
                 "confirmed by the SI's own printed ESI-MS")
                if calcd is not None else
                ("structure built from the synthesis the SI describes in prose "
                 "AND draws as a product in Scheme S8; the SI prints 1H NMR "
                 "for this compound but NO ESI-MS, so the formula has no mass "
                 "confirmation -- it is the only one of the 40 without one"),
        })

    # --- which mass convention? -------------------------------------------
    n = len(resid["proton"])
    worst = {k: max(v) for k, v in resid.items()}
    best = min(worst, key=worst.get)
    other = "h_atom" if best == "proton" else "proton"
    print(f"ESI-MS convention over {n} printed [M+H]+ values:")
    print(f"   as a PROTONATED ion (M + 1.00728): worst |delta| "
          f"{worst['proton']*1000:.2f} mDa")
    print(f"   as M + a neutral H atom (1.00783): worst |delta| "
          f"{worst['h_atom']*1000:.2f} mDa")
    print(f"   -> the SI reports {best}; margin over the next-best reading "
          f"{(worst[other]-worst[best])*1000:.2f} mDa")

    tol = 0.0015   # 1.5 mDa; the SI prints 3 decimals, so +-0.5 mDa is rounding
    bad = [r for r in rows if r["our_mz_minus_si_calcd_mDa"] not in ("",)
           and abs(float(r["our_mz_minus_si_calcd_mDa"])) > tol * 1000]
    print(f"\nHRMS agreement: {n - len(bad)}/{n} within {tol*1000:.1f} mDa")
    if bad:
        for r in bad:
            print(f"   MISMATCH {r['IL_name']} {r['formula']} "
                  f"{r['our_mz_minus_si_calcd_mDa']} mDa")
        raise SystemExit(
            "stage03: an enumerated structure disagrees with the SI's own "
            "printed exact mass. That is a formula error -- an arm count, a "
            "chain length or the head-arm carbon count. Do not widen the "
            "tolerance.")

    if len(keys) != 40:
        raise SystemExit(f"stage03: {len(keys)} distinct structures, need 40.")

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {len(rows)} lipids to {OUT}")
    print(f"distinct full InChIKeys: {len(keys)}/40")
    sk = {k.split('-')[0] for k in keys}
    # 4 ATX motifs x 2 linkers x 3 DISTINCT RETINOID CONNECTIVITIES. Three,
    # not five: ATRA, 13CRA and 9CRA are geometric isomers of one skeleton,
    # leaving {that skeleton, AVA, Bex}. My first version of this assertion
    # said 32 and was wrong -- it counted 4 connectivities. The check fired on
    # its own author, which is the point of writing it as an assertion.
    print(f"distinct InChIKey SKELETON blocks: {len(sk)} "
          f"(expected 24 = 4 motifs x 2 linkers x 3 distinct retinoid "
          f"connectivities, since ATRA/13CRA/9CRA share one skeleton)")
    if len(sk) != 24:
        raise SystemExit(
            f"stage03: {len(sk)} skeleton blocks, expected 24. The three "
            "retinoic acid geometric isomers must share a skeleton and AVA "
            "and Bex must each have their own.")
    import collections
    sizes = collections.Counter(collections.Counter(
        k.split('-')[0] for k in keys).values())
    print(f"skeleton multiplicities: {dict(sizes)} "
          "(8 skeletons of 3 = the RA isomer triples, 16 of 1 = AVA and Bex)")
    if sizes != collections.Counter({1: 16, 3: 8}):
        raise SystemExit(
            f"stage03: skeleton multiplicities {dict(sizes)}, expected "
            "{3: 8, 1: 16}. 8 x 3 + 16 x 1 = 40.")


if __name__ == "__main__":
    main()
