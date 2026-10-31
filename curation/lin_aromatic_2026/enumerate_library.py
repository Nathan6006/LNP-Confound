#!/usr/bin/env python3
"""Enumerate the 840 aromatic ionizable lipids of Lin et al., JACS 2026.

Chemistry, from the Results text and the SI Materials list:

  An amine precursor nA (8A-12A) is conjugated with an aromatic analog mB
  (1B-7B), then the REMAINING N-H sites are capped with a hydrophobic tail.
  "The notations nA1 and nA2 denoted reactions where one amine was conjugated
  with one and two equiv of aromatic analogs, respectively, while the remaining
  reactive [sites were capped with tails]."

  Both couplings are N-H additions, and which one applies is decided by the
  block's own functional group, not by a lookup table:
     epoxide   N-H + oxirane      -> N-CH2-CH(OH)-R      (ring-opening)
     acrylate  N-H + CH2=CR'-CO2R -> N-CH2-CHR'-CO2R     (Michael addition)
  1B-5B and EP6-EP18 are epoxides; 6B, 7B and Ac10/12/18/8-1/10-1 are
  acrylates. 7B is a METHACRYLATE, so its arm carries an alpha-methyl - handled
  by reading the double bond's substitution rather than assuming acrylate.

  5 amines x 7 aromatics x 12 tails x {1 equiv, 2 equiv} = 840

The epoxide is opened by RWMol surgery at the LESS SUBSTITUTED ring carbon
(the CH2), which is the regiochemistry of amine ring-opening, rather than by
string manipulation, so the aromatic epoxides with ether linkers (2B-5B) are
handled by the same code as the simple 1,2-epoxyalkanes.

CONFIDENCE. The amine, aromatic and linear-tail structures are named outright
in the SI Materials list and carry those names here. Ac8-1 and Ac10-1 are named
"iso-octyl acrylate" and "iso-decyl acrylate" - commercial isomer MIXTURES - so
a single SMILES is a drawing convention and those lipids are flagged medium.
Positional ambiguity is also flagged: where fewer arms are added than there are
N-H sites, the paper does not say which nitrogens react, so a convention is
applied (primary before secondary, RDKit index as tiebreak) exactly as for the
iPhos library.

Writes 03_enumerated.csv and 06_merged.csv (structures joined to Fig 2e values).
"""
import csv
import json

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors as D
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")

BLOCKS = {r["block_id"]: r for r in csv.DictReader(open("blocks.csv"))}

AMINES = ["8A", "9A", "10A", "11A", "12A"]
AROMATICS = [f"{m}B" for m in range(1, 8)]
TAILS = ["EP6", "EP8", "EP10", "EP12", "EP14", "EP16", "EP18",
         "Ac10", "Ac12", "Ac18", "Ac8-1", "Ac10-1"]

# named in the SI Materials list as commercial isomer mixtures
SUPPLIED = {
    "Ac8-1": ("C=CC(=O)OCCCCCC(C)C", "iso-octyl acrylate", "medium"),
    "Ac10-1": ("C=CC(=O)OCCCCCCCC(C)C", "iso-decyl acrylate", "medium"),
}


def block_smiles(bid):
    if bid in SUPPLIED:
        return SUPPLIED[bid][0]
    return BLOCKS[bid]["SMILES"]


def nh_sites(mol):
    """(atom index, free H count) for every N that can still take an arm."""
    out = []
    for a in mol.GetAtoms():
        if a.GetSymbol() != "N" or a.GetIsAromatic():
            continue
        h = a.GetTotalNumHs()
        if h:
            out.append((a.GetIdx(), h))
    return out


def add_epoxide(rw, n_idx, frag_smiles):
    """N-H + oxirane -> N-CH2-CH(OH)-R, opening at the less substituted carbon."""
    frag = Chem.RWMol(Chem.MolFromSmiles(frag_smiles))
    ring = None
    for atom in frag.GetAtoms():
        if atom.GetSymbol() != "O" or atom.GetDegree() != 2:
            continue
        nb = [x.GetIdx() for x in atom.GetNeighbors()]
        if len(nb) == 2 and frag.GetBondBetweenAtoms(nb[0], nb[1]) is not None:
            ring = (atom.GetIdx(), nb[0], nb[1])
            break
    if ring is None:
        return None
    o, c1, c2 = ring
    # the CH2 is the ring carbon with fewer heavy neighbours
    def heavy(i):
        return sum(1 for x in frag.GetAtomWithIdx(i).GetNeighbors()
                   if x.GetAtomicNum() > 1)
    ch2, ch = (c1, c2) if heavy(c1) <= heavy(c2) else (c2, c1)
    frag.RemoveBond(o, ch2)                      # open the ring at the CH2
    frag.GetAtomWithIdx(o).SetNoImplicit(False)  # O becomes the hydroxyl
    off = rw.GetNumAtoms()
    rw.InsertMol(frag.GetMol())
    rw.AddBond(n_idx, off + ch2, Chem.BondType.SINGLE)
    return rw


def add_acrylate(rw, n_idx, frag_smiles):
    """N-H + CH2=CR'-CO2R -> N-CH2-CHR'-CO2R (Michael addition)."""
    frag = Chem.RWMol(Chem.MolFromSmiles(frag_smiles))
    target = None
    for bond in frag.GetBonds():
        if bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetSymbol() != "C" or b.GetSymbol() != "C":
            continue
        # the terminal CH2 is the one with fewer heavy neighbours
        ha = sum(1 for x in a.GetNeighbors() if x.GetAtomicNum() > 1)
        hb = sum(1 for x in b.GetNeighbors() if x.GetAtomicNum() > 1)
        ch2, alpha = (a, b) if ha <= hb else (b, a)
        target = (bond.GetIdx(), ch2.GetIdx(), alpha.GetIdx())
        break
    if target is None:
        return None
    bidx, ch2, _alpha = target
    frag.GetBondWithIdx(bidx).SetBondType(Chem.BondType.SINGLE)
    off = rw.GetNumAtoms()
    rw.InsertMol(frag.GetMol())
    rw.AddBond(n_idx, off + ch2, Chem.BondType.SINGLE)
    return rw


def is_epoxide(smiles):
    m = Chem.MolFromSmiles(smiles)
    return m is not None and m.HasSubstructMatch(Chem.MolFromSmarts("C1CO1"))


def build(amine_id, arom_id, tail_id, n_equiv):
    base = Chem.MolFromSmiles(block_smiles(amine_id))
    sites = nh_sites(base)
    total = sum(h for _, h in sites)
    if n_equiv > total:
        return None, None, None

    # Aromatic arms go ONE PER NITROGEN, primary before secondary, only
    # doubling up on a nitrogen once every eligible one has an arm. This is not
    # a guess: the SI's "precise synthesis" of 10A2-3B-EP10 protects spermine's
    # primary amines as trifluoroacetamides, Boc-protects the secondaries, frees
    # the two PRIMARY amines, then adds 2.0 eq of 3B - one aromatic per primary
    # nitrogen - before deprotecting and capping the rest with EP10. An earlier
    # version put both aromatic arms on a single primary nitrogen; the molecular
    # formula is identical either way, so HRMS cannot distinguish them, but the
    # connectivity was wrong.
    order = sorted(sites, key=lambda s: (-s[1], s[0]))
    plan = []
    left = n_equiv
    for rnd in (1, 2):
        for idx, h in order:
            if left <= 0:
                break
            if h >= rnd and sum(1 for i, k in plan if i == idx) < rnd:
                plan.append((idx, "aromatic"))
                left -= 1
    used = {}
    for idx, _ in plan:
        used[idx] = used.get(idx, 0) + 1
    for idx, h in order:
        for _ in range(h - used.get(idx, 0)):
            plan.append((idx, "tail"))

    arom, tail = block_smiles(arom_id), block_smiles(tail_id)
    rw = Chem.RWMol(base)
    for idx, kind in plan:
        frag = arom if kind == "aromatic" else tail
        fn = add_epoxide if is_epoxide(frag) else add_acrylate
        if fn(rw, idx, frag) is None:
            return None, None, None
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None, None, None
    n_tails = sum(1 for _, k in plan if k == "tail")
    return Chem.MolToSmiles(out), n_tails, total


def main():
    rows, fails = [], []
    for amine in AMINES:
        for eq in (1, 2):
            for arom in AROMATICS:
                for tail in TAILS:
                    smi, n_tails, n_sites = build(amine, arom, tail, eq)
                    name = f"{amine}{eq}-{arom}-{tail}"
                    if smi is None:
                        fails.append(name)
                        continue
                    m = Chem.MolFromSmiles(smi)
                    conf = "high"
                    if tail in SUPPLIED:
                        conf = "medium"
                    rows.append({
                        "lipid_name": name,
                        "amine": amine, "equiv": eq, "aromatic": arom, "tail": tail,
                        "IL_SMILES": smi,
                        "n_NH_sites": n_sites,
                        "n_aromatic_arms": eq, "n_tail_arms": n_tails,
                        "positional_ambiguity": int(eq < n_sites),
                        "MolWt": round(D.MolWt(m), 2),
                        "ExactMolWt": round(D.ExactMolWt(m), 4),
                        "formula": CalcMolFormula(m),
                        "num_N": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N"),
                        "confidence": conf,
                        "amine_name": BLOCKS[amine]["name_in_paper"],
                        "aromatic_name": BLOCKS[arom]["name_in_paper"],
                        "tail_name": (SUPPLIED[tail][1] if tail in SUPPLIED
                                      else BLOCKS[tail]["name_in_paper"]),
                        "source": "Fig 1b-1d building blocks + SI Materials list",
                        "extraction_method": "enumerated_from_named_reagents",
                    })
    with open("03_enumerated.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    vals = {r["lipid_name"]: r for r in csv.DictReader(open("04_fig2e_840.csv"))}
    merged, missing = [], []
    for r in rows:
        v = vals.get(r["lipid_name"])
        if not v:
            missing.append(r["lipid_name"])
            continue
        out = dict(r)
        out.update({k: v[k] for k in ("bin_label", "bin_index", "rlu_midpoint",
                                      "cell_rgb", "match_distance", "assay",
                                      "readout_metric", "value_scale",
                                      "source_figure")})
        merged.append(out)
    with open("06_merged.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(merged[0]))
        w.writeheader()
        w.writerows(merged)

    uniq = len({r["IL_SMILES"] for r in rows})
    qc = {
        "paper": "Lin et al., JACS 2026, 148, 6960-6969",
        "enumerated": len(rows), "expected": 840, "failed": fails,
        "unique_SMILES": uniq,
        "merged_with_values": len(merged), "missing_values": missing,
        "positional_ambiguity_rows": sum(r["positional_ambiguity"] for r in rows),
        "medium_confidence_rows": sum(1 for r in rows if r["confidence"] == "medium"),
        "medium_confidence_reason": "Ac8-1 and Ac10-1 are named in the SI as "
                                    "iso-octyl / iso-decyl acrylate, commercial isomer "
                                    "MIXTURES; a single SMILES is a convention.",
        "MW_range": [min(r["MolWt"] for r in rows), max(r["MolWt"] for r in rows)],
    }
    json.dump(qc, open("09_enumeration_qc.json", "w"), indent=1)
    print(f"enumerated {len(rows)}/840  unique SMILES {uniq}  failed {len(fails)}")
    print(f"merged with Fig 2e values: {len(merged)}  missing {len(missing)}")
    print(f"MW range {qc['MW_range']}")
    print(f"positional ambiguity rows {qc['positional_ambiguity_rows']}, "
          f"medium confidence {qc['medium_confidence_rows']}")


if __name__ == "__main__":
    main()
