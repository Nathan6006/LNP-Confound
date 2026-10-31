#!/usr/bin/env python3
"""Enumerate the Siegwart esterified-cholesterol lipid library (Tian et al., Adv Mater 2026).

Design, read from Fig 1b and the Results text:
  step 1  amine head (C1-C6) + alkyl acrylate (O1/O2) -> Michael addition at EVERY N-H,
          giving a tertiary amino lipid core
  step 2  that core + a cholesterol-derived alkyl bromide (B1/B2/B3) -> quaternary
          ammonium cholesterol-conjugated lipid
  naming  CxOyBz;  lead compound is C6O2B2

Arithmetic check: (6 synthesized heads x 2 alkyls = 12 cores, plus 5 commercial cores
SM-102 / OF-2 / TT3 / MC3 / Lipid 14) x 3 cholesterol tails = 51, matching the stated
"library of 51 new lipids".

Stoichiometry independently confirmed by NMR rather than assumed: C6O2 shows
delta 4.04 (t, 10H) = five ester OCH2 groups and delta 0.87 (t, 15H) = five terminal CH3,
i.e. FIVE acrylate arms - exactly the five N-H sites in C6. C6O2B2 keeps both integrals
and adds the cholesterol signatures (5.30 d 1H vinyl, 4.54 m 1H C3-H, 0.61 s 3H C18-Me).

Confidence is recorded per building block; anything not independently confirmed is
flagged rather than presented as certain.
"""
import csv
import json

from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors as D
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

RDLogger.DisableLog("rdApp.*")

# --- amine heads, read from Fig 1b -------------------------------------------
HEADS = {
    "C1": ("CNCCN", "N-methylethylenediamine", "figure", 3),
    "C2": ("CN(C)CCN", "N,N-dimethylethylenediamine", "figure", 2),
    "C3": ("NCCN(CCN)CCN", "tris(2-aminoethyl)amine", "figure", 6),
    "C4": ("CNCCN(CCNC)CCNC", "tris(2-(methylamino)ethyl)amine", "figure", 3),
    "C5": ("NCCCN1CCN(CCCN)CC1", "1,4-bis(3-aminopropyl)piperazine", "figure", 4),
    "C6": ("NCCN1CCN(CCNCCN)CC1", "aminoethyl-piperazine-aminoethyl-aminoethyl",
           "figure + NMR (5 N-H confirmed by 10H OCH2 / 15H CH3 in C6O2)", 5),
}

# --- alkyl acrylates ----------------------------------------------------------
# O2 is pinned by the C6O2 integrals: the 1.27 envelope (92H over 5 chains, ~18H each)
# plus OCH2 and the beta-CH2 is consistent with a linear C12 (dodecyl) chain.
ALKYLS = {
    "O2": ("C=CC(=O)OCCCCCCCCCCCC", "dodecyl acrylate", "figure + NMR integrals", "high"),
    "O1": ("C=CC(=O)OCCCCCCCC(C)C", "isodecyl acrylate (branched terminus)",
           "figure only - branch position not independently confirmed", "medium"),
}

# --- cholesterol-derived alkyl bromides --------------------------------------
# Text: "tunable spacing between the cholesterol group and the amine core using
# 3-, 5-, or 7-carbon linkers". Modelled as the cholesteryl ester of a
# bromo-alkanoic acid. Linker length is stated; the ester connectivity is inferred
# from the C3-H at 4.54 ppm, so these carry a medium confidence flag.
CHOL = "CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
TAILS = {
    "B1": (3, "cholesteryl ester, 3-carbon linker"),
    "B2": (5, "cholesteryl ester, 5-carbon linker"),
    "B3": (7, "cholesteryl ester, 7-carbon linker"),
}


def michael(head_smiles, acrylate_smiles, n_arms):
    """Add acrylate across every N-H: N-H + CH2=CH-CO2R -> N-CH2CH2-CO2R."""
    arm = acrylate_smiles.replace("C=CC(=O)", "CCC(=O)", 1)  # saturated arm fragment
    mol = Chem.MolFromSmiles(head_smiles)
    rw = Chem.RWMol(mol)
    placed = 0
    for a in list(rw.GetAtoms()):
        if a.GetSymbol() != "N":
            continue
        nh = a.GetTotalNumHs()
        for _ in range(nh):
            if placed >= n_arms:
                break
            frag = Chem.MolFromSmiles(arm)
            off = rw.GetNumAtoms()
            rw.InsertMol(frag)
            rw.AddBond(a.GetIdx(), off, Chem.BondType.SINGLE)
            placed += 1
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return Chem.MolToSmiles(out), placed


ATTACH_MAP = 9          # atom map number tagging the fragment's free chain end


def _cholesteryl_ester(n_link):
    """N-(CH2)(n-1)-C(=O)-O-C3(cholesteryl), by esterifying the REAL C3-OH.

    Two bugs have lived in this function. Both are described here because each
    fix introduced the next, which is the actual lesson.

    1. The first version did `CHOL.replace("(O)", "", 1)` and concatenated the
       acyl group in front of the resulting string. Deleting the hydroxyl left
       C3 as a 2-coordinate carbon carrying a RADICAL, and because SMILES
       concatenation bonds to the FIRST atom of the appended fragment, the ester
       oxygen attached to `CC(C)CCC...` - a terminal methyl of the isooctyl side
       chain, the opposite end of the steroid from C3.

    2. The fix for that was RWMol fragment assembly: bond the acid's oxygen to
       C3, then `RemoveAtom` cholesterol's own hydroxyl oxygen. That silently
       INVERTED C3. `RWMol.RemoveAtom` drops the oxygen out of the middle of
       C3's neighbour list without correcting C3's chiral tag, so the tag then
       means the opposite of what it meant. Every one of the 537 corpus rows
       carried epicholesterol (3alpha, InChIKey HVYWMOMLDIMFJA-VEIPTCAHSA-N)
       instead of cholesterol (3beta, -DPAQBDIFSA-N). Formula, exact mass and
       all 17 descriptor columns are bit-identical between the two epimers, so
       C1-C9 were blind, and C10 was blind too because its reference structure
       was a byte-for-byte copy of this function's own output.

    This version does the condensation the way the chemistry does it: the ACID
    loses its OH and the ALCOHOL's oxygen becomes the ester oxygen. No atom
    adjacent to a stereocentre is deleted, so no chiral tag is left stale. The
    only atom removed is the acid's terminal oxygen, which has no stereochemistry
    anywhere near it.

    The paper is unambiguous that this is 3beta: Fig 1b (p.3) draws B1/B2/B3 with
    the C3-O bond as a bold wedge cofacial with the C18/C19 angular methyls,
    Fig 2a draws the same relation in a different orientation, the 1H C3-H
    appears at 4.71-4.58 (the broad axial tt of a 3beta-acyloxy sterol; a 3alpha
    ester would be near 5.0 and far narrower), and the Fig S7 13C spectrum has
    an isolated peak at 73.9 ppm where a 3alpha ester would sit near 70.9.
    Nothing in the paper or SI mentions an epi, allo or 3alpha sterol at all.

    The free chain end is tagged with ATTACH_MAP rather than left to be found by
    position: `quaternize` used to bond the amine to atom index 0 of this
    fragment's canonical SMILES, which happens to be the chain terminus for all
    three linkers but is an unstated invariant of exactly the kind that produced
    bug 1 above.
    """
    frag = Chem.RWMol(Chem.MolFromSmiles("C" * (n_link - 1) + "C(=O)O"))
    acid_oh = frag.GetNumAtoms() - 1          # the ACID's terminal -OH oxygen
    acyl_c = frag.GetAtomWithIdx(acid_oh).GetNeighbors()[0].GetIdx()
    frag.GetAtomWithIdx(0).SetAtomMapNum(ATTACH_MAP)    # the free chain end
    frag.RemoveAtom(acid_oh)                  # condensation: the ACID loses OH
    chol = Chem.MolFromSmiles(CHOL)
    oh = next((a for a in chol.GetAtoms()
               if a.GetSymbol() == "O" and a.GetDegree() == 1), None)
    if oh is None:
        return None, None
    off = frag.GetNumAtoms()
    frag.InsertMol(chol)
    # bond to the ALCOHOL's oxygen, which is retained -- C3's neighbour list is
    # never touched, so its chiral tag stays valid
    frag.AddBond(acyl_c, off + oh.GetIdx(), Chem.BondType.SINGLE)
    out = frag.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None, None
    return Chem.MolToSmiles(out), 0


def _quaternization_site(mol):
    """Pick the nitrogen that is actually alkylated.

    Fig 2a shows the N+ of C6O2B2 is a PIPERAZINE RING nitrogen, and the
    cholesteryl-ester linker hangs off it. The previous code took
    `max(atom index)`, an artifact of SMILES writing order, which for C6 landed
    on the terminal nitrogen bearing two beta-amino-ester arms - the LEAST
    nucleophilic nitrogen in the molecule, and not the one drawn.

    Rule, in order: tertiary nitrogens only; prefer a ring (piperazine) N;
    among those prefer the fewest electron-withdrawing beta arms (the most
    nucleophilic); then prefer the least sterically hindered; break remaining
    ties by atom index. For C5/C6 this selects a piperazine N, matching Fig 2a.
    For C1-C4, which have no ring, the paper draws no structure and the site is
    an ASSUMPTION - flagged per row.

    BOTH beta-ester AND beta-acyloxy arms are counted, which the first repair of
    this function did not do, and the omission mattered. The old
    `beta_ester_arms` walked only N-CH2-CH2-C(=O), so it could not see an arm of
    the form N-CH2CH2-O-C(=O)R. On `Lipid 14` both candidate nitrogens therefore
    scored 0 and the rule fell through to the atom-index tiebreak -- making the
    answer a function of how the SMILES happened to be written, which is the
    exact defect this function was rewritten to eliminate. Rewriting the same
    molecule tail-first made the two rules swap answers. Counting beta-acyloxy
    as well makes the choice order-independent: the terminal dimethylamino wins
    from either spelling, and C5/C6 still resolve to the piperazine N of Fig 2a.

    The hindrance tiebreak is secondary to the beta count on purpose -- an
    electron-poor nitrogen is a worse nucleophile than a merely crowded one -
    and it only ever decides a tie the beta count leaves open.
    """
    cands = [a for a in mol.GetAtoms()
             if a.GetSymbol() == "N" and a.GetTotalNumHs() == 0
             and sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1) == 3]
    if not cands:
        return None, None

    def beta_ewg_arms(a):
        """Electron-withdrawing arms two bonds out from this nitrogen.

        Counts N-CH2-CH2-C(=O)-  (beta-ester, the Michael adduct arm) and
        N-CH2-CH2-O-C(=O)-      (beta-acyloxy, e.g. Lipid 14's ester arm).
        Both pull density off the nitrogen and both make it a worse nucleophile,
        so for the purpose of choosing an alkylation site they are one thing.
        """
        n = 0
        for nb in a.GetNeighbors():
            if nb.GetSymbol() != "C":
                continue
            for nb2 in nb.GetNeighbors():
                if nb2.GetIdx() == a.GetIdx():
                    continue
                if nb2.GetSymbol() == "C":          # beta-ester: ...C-C(=O)
                    for nb3 in nb2.GetNeighbors():
                        if nb3.GetSymbol() == "C" and any(
                                b.GetBondType() == Chem.BondType.DOUBLE
                                and b.GetOtherAtom(nb3).GetSymbol() == "O"
                                for b in nb3.GetBonds()):
                            n += 1
                elif nb2.GetSymbol() == "O":        # beta-acyloxy: ...C-O-C(=O)
                    for nb3 in nb2.GetNeighbors():
                        if nb3.GetIdx() == nb.GetIdx() or nb3.GetSymbol() != "C":
                            continue
                        if any(b.GetBondType() == Chem.BondType.DOUBLE
                               and b.GetOtherAtom(nb3).GetSymbol() == "O"
                               for b in nb3.GetBonds()):
                            n += 1
        return n

    def hindrance(a):
        """Heavy atoms on the carbons alpha to this nitrogen: a methyl on N
        contributes 0, a long alkyl chain contributes 1 each. Distinguishes a
        terminal N(CH3)2 from a trialkylamine flanked by two C8 chains."""
        return sum(1 for nb in a.GetNeighbors()
                   for nb2 in nb.GetNeighbors()
                   if nb2.GetIdx() != a.GetIdx() and nb2.GetAtomicNum() > 1)

    ranked = sorted(cands, key=lambda a: (not a.IsInRing(), beta_ewg_arms(a),
                                          hindrance(a), a.GetIdx()))
    chosen = ranked[0]
    return chosen, ("drawn in Fig 2a (piperazine N)" if chosen.IsInRing()
                    else "ASSUMED: no ring N and the paper draws no structure")


def quaternize(core_smiles, n_link):
    """Alkylate the selected tertiary N with cholesteryl-O2C-(CH2)n-."""
    tail_frag, _ = _cholesteryl_ester(n_link)
    if tail_frag is None:
        return None, None
    mol = Chem.MolFromSmiles(core_smiles)
    target, site_note = _quaternization_site(mol)
    if target is None:
        return None, None
    rw = Chem.RWMol(mol)
    frag = Chem.MolFromSmiles(tail_frag)
    if frag is None:
        return None, None
    tidx = target.GetIdx()
    off = rw.GetNumAtoms()
    rw.InsertMol(frag)
    # Bond to the atom TAGGED as the free chain end, not to index 0 of the
    # fragment's canonical SMILES. Index 0 is the chain terminus for all three
    # linkers today, but that is an unstated invariant and this function has
    # already shipped one bug that was exactly an unstated invariant about
    # which atom a concatenation attaches to.
    attach = next((a.GetIdx() for a in frag.GetAtoms()
                   if a.GetAtomMapNum() == ATTACH_MAP), None)
    if attach is None:
        return None, None
    # An atom-mapped atom round-trips through SMILES in bracket form, `[CH3:9]`,
    # which pins its hydrogen count; adding a bond to it without clearing that
    # makes it hypervalent and SanitizeMol rejects the whole molecule.
    at = rw.GetAtomWithIdx(off + attach)
    at.SetAtomMapNum(0)
    at.SetNoImplicit(False)
    at.SetNumExplicitHs(0)
    rw.AddBond(tidx, off + attach, Chem.BondType.SINGLE)
    rw.GetAtomWithIdx(tidx).SetFormalCharge(1)
    rw.GetAtomWithIdx(tidx).SetNoImplicit(True)
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None, None
    if any(a.GetNumRadicalElectrons() for a in out.GetAtoms()):
        return None, None          # hard guard: this is how the old bug shipped
    return Chem.MolToSmiles(out), site_note


def main():
    rows = []
    for hk, (hs, hname, hsrc, narm) in HEADS.items():
        for ak, (as_, aname, asrc, aconf) in ALKYLS.items():
            core, placed = michael(hs, as_, narm)
            for bk, (nlink, bname) in TAILS.items():
                smi, site_note = quaternize(core, nlink)
                if smi is None:
                    continue
                mol = Chem.MolFromSmiles(smi)
                rows.append({
                    "lipid_name": f"{hk}{ak}{bk}",
                    "amine_head": hk, "head_desc": hname,
                    "alkyl": ak, "alkyl_desc": aname,
                    "chol_tail": bk, "chol_desc": bname,
                    "n_acrylate_arms": placed,
                    "IL_SMILES": smi,
                    "MolWt": round(D.MolWt(mol), 2),
                    "formula": CalcMolFormula(mol),
                    "num_N": sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "N"),
                    "quaternization_site": site_note,
                    "quaternization_site_assumed": int(site_note.startswith("ASSUMED")),
                    "confidence": ("high" if (hk in ("C5", "C6") and ak == "O2")
                                   else ("low" if site_note.startswith("ASSUMED")
                                         else aconf)),
                    "head_evidence": hsrc,
                    "source_pdf": "Adv Mater 2026 Tian (main)",
                    "figure": "Fig 1b", "extraction_method": "structure_diagram+NMR",
                })
    with open("03_enumerated.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"enumerated {len(rows)} synthesized-core lipids "
          f"(6 heads x 2 alkyls x 3 tails = 36 expected)")
    print(f"unique SMILES: {len({r['IL_SMILES'] for r in rows})}")
    lead = [r for r in rows if r["lipid_name"] == "C6O2B2"]
    if lead:
        r = lead[0]
        print(f"\nlead compound C6O2B2: {r['formula']}  MW {r['MolWt']}  "
              f"arms={r['n_acrylate_arms']} (NMR says 5)")


if __name__ == "__main__":
    main()
