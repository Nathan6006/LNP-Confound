#!/usr/bin/env python3
"""Combine blocks into the full library -> 03_enumerated.csv

Run from this folder:   python3 -I stage03_enumerate.py   (after stage02)

WHAT IS ENUMERATED. 8 headgroups x 13 tails = 104 members, which is the number
the abstract claims ("Synthesized from combinations of 8 different polyamines
as headgroups and 13 acrylates as tails, 104 different lipidoids"). SI Table S3
prints a score grid of exactly 8 x 13 and SI Table S1 lists 13 members for the
2N headgroup alone, so the claim is a plain product and is reconciled, not
assumed -- see 07_qc.json check P1.

TWO STRUCTURES PER MEMBER, and the distinction matters:

  parent      the Michael adduct. Every N-H of the polyamine is consumed by one
              acrylate, giving 4/5/6 arms depending on the headgroup. The
              article states those counts independently and they match the N-H
              count on all 8 headgroups (stage02).

  methylated  what was actually TESTED. "lipidoid dispersion in aqueous media
              was facilitated by introducing permanent positive charges into
              the headgroup using methyl iodide." Figure 1A draws the product
              of that step in full for the 3NB (TREN) case: ALL FOUR nitrogens
              carry a methyl and a '+', the central tertiary nitrogen included.
              So methylation is exhaustive, and the tested species is an
              N-fold quaternary ammonium cation with N iodide counterions.

HOW THEY ARE BUILT. Arms go on through curation/lib/assemble.attach, by atom
map, one N-H at a time; the loop terminates when no nitrogen has a hydrogen
left, and the number of iterations is asserted against the headgroup's stated
tail count. Because every arm on a given member is identical, the order the
nitrogens are filled in cannot change the product.

THE ONE STEP THAT IS NOT assemble's. Quaternisation cannot go through `attach`:
adding a fourth substituent to a neutral nitrogen makes RDKit's sanitiser
reject the valence before the charge can be set, so `attach` would raise on
every member. `_quaternise` below does it with AddAtom/AddBond and
SetFormalCharge -- and with NO RemoveAtom anywhere, so the F-013 mechanism is
structurally absent rather than merely avoided. It then re-runs the library's
own assertions (check_no_radicals, check_stereo_preserved) plus an explicit
atom-count and formula ledger, because a local helper that nobody audited is
exactly the thing error class 49 is about.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import assemble                                        # noqa: E402
from errors import StructureError                      # noqa: E402
from rdkit import Chem                                 # noqa: E402
from rdkit.Chem import Descriptors, rdMolDescriptors   # noqa: E402

BLOCKS = os.path.join(HERE, "01_blocks.csv")
OUT = os.path.join(HERE, "03_enumerated.csv")

# Items whose structure rests on a convention rather than on a measurement, and
# the sentence that says WHICH feature. These drive the quality columns at emit
# time; they are NOT literals applied to every row (error class 29).
TAIL_CAVEATS = {
    "13B": ("structure_assumed",
            "Tail 13B is isodecyl acrylate, which is sold as a mixture of "
            "branched C10 isomers; Figure 1A draws one representative "
            "(8-methylnonyl) and the paper gives no NMR, HRMS or supplier "
            "isomer statement that would pin the branch. CLAUDE.md section 6 "
            "already lists 'the exact branch of iso-decyl acrylate' as "
            "unverifiable from the paper alone."),
    "9C": ("structure_assumed",
           "Tail 9C is the tricyclo[5.2.1.0^2,6]decanyl (hydrogenated "
           "dicyclopentadienyl) ester. Figure 1A draws the cage and puts the "
           "attachment on the two-carbon bridge next to a bridgehead, which "
           "is the commercial 8-yl; the 8-yl/9-yl and endo/exo isomers are "
           "not distinguishable from the drawing and the commercial monomer "
           "is a mixture of them."),
    "8CB": ("structure_assumed",
            "Tail 8CB is isobornyl acrylate. The constitution is fixed by the "
            "drawn cage and by its 208.3 Da entry in Table S4, but the paper "
            "states no stereochemistry and commercial isobornyl acrylate "
            "contains bornyl isomer; emitted without stereocentres assigned."),
}


def _quaternise(parent_smiles, label=""):
    """Methylate EVERY nitrogen and give each a +1 charge. No atom is removed.

    Returns (cation SMILES, n_methyls). Raises unless the product's ledger is
    exactly parent + one CH3 per nitrogen and every nitrogen ends up
    four-coordinate, charge +1, with no hydrogens left.
    """
    parent = Chem.MolFromSmiles(parent_smiles)
    if parent is None:
        raise StructureError(f"_quaternise({label!r}): unparseable parent.")
    n_idx = [a.GetIdx() for a in parent.GetAtoms() if a.GetSymbol() == "N"]
    if not n_idx:
        raise StructureError(f"_quaternise({label!r}): no nitrogen to methylate.")
    for i in n_idx:
        a = parent.GetAtomWithIdx(i)
        if a.GetTotalNumHs() != 0 or a.GetDegree() != 3:
            raise StructureError(
                f"_quaternise({label!r}): N{i} is {a.GetDegree()}-coordinate "
                f"with {a.GetTotalNumHs()} H. Every nitrogen must already be "
                "a fully substituted tertiary amine before methylation; one "
                "that is not means an arm was not attached, and quaternising "
                "it would hide that.")

    rw = Chem.RWMol(parent)
    for i in n_idx:
        c = rw.AddAtom(Chem.Atom(6))
        rw.AddBond(i, c, Chem.BondType.SINGLE)
        rw.GetAtomWithIdx(i).SetFormalCharge(1)
        rw.GetAtomWithIdx(i).SetNoImplicit(False)
    prod = rw.GetMol()
    Chem.SanitizeMol(prod)

    # --- the ledger: parent + one CH3 per N, nothing else ---------------
    if prod.GetNumAtoms() != parent.GetNumAtoms() + len(n_idx):
        raise StructureError(
            f"_quaternise({label!r}): heavy-atom count "
            f"{prod.GetNumAtoms()} != parent {parent.GetNumAtoms()} + "
            f"{len(n_idx)} methyl carbons.")
    if Chem.GetFormalCharge(prod) != len(n_idx):
        raise StructureError(
            f"_quaternise({label!r}): net charge "
            f"{Chem.GetFormalCharge(prod)} != {len(n_idx)} nitrogens.")
    for a in prod.GetAtoms():
        if a.GetSymbol() == "N" and not (a.GetDegree() == 4
                                         and a.GetFormalCharge() == 1
                                         and a.GetTotalNumHs() == 0):
            raise StructureError(
                f"_quaternise({label!r}): N{a.GetIdx()} ended up "
                f"degree {a.GetDegree()}, charge {a.GetFormalCharge()}, "
                f"{a.GetTotalNumHs()} H.")
    assemble.check_no_radicals(prod, label=label or "quaternise")
    assemble.check_stereo_preserved(prod, [parent],
                                    label=label or "quaternise")
    return Chem.MolToSmiles(prod), len(n_idx)


def build_parent(head_smiles, arm_smiles, expect_arms, label=""):
    """Attach one arm per N-H until no nitrogen has a hydrogen left."""
    cur = Chem.MolFromSmiles(head_smiles)
    n = 0
    while True:
        target = None
        for a in cur.GetAtoms():
            if a.GetSymbol() == "N" and a.GetTotalNumHs() > 0:
                target = a.GetIdx()
                break
        if target is None:
            break
        if n >= 12:
            raise StructureError(f"build_parent({label!r}): runaway loop.")
        for a in cur.GetAtoms():
            a.SetAtomMapNum(0)
        cur.GetAtomWithIdx(target).SetAtomMapNum(1)
        cur = Chem.MolFromSmiles(
            assemble.attach(cur, arm_smiles, a_map=1, b_map=2,
                            label=f"{label} arm {n + 1}"))
        n += 1
    if n != expect_arms:
        raise StructureError(
            f"build_parent({label!r}): attached {n} arms, the headgroup's "
            f"N-H count and the article's stated tail count say {expect_arms}.")
    return Chem.MolToSmiles(cur), n


def main():
    heads, tails = [], []
    with open(BLOCKS) as fh:
        for r in csv.DictReader(fh):
            (heads if r["kind"] == "headgroup" else tails).append(r)
    if len(heads) != 8 or len(tails) != 13:
        raise SystemExit(f"blocks: {len(heads)} heads x {len(tails)} tails, "
                         "expected 8 x 13")

    rows = []
    for h in heads:
        for t in tails:
            name = f"{h['name']}{t['name']}"
            # the arm = propanoate ester, atom-mapped where it meets the N
            alcohol = t["note"].split("alcohol R-OH as built: ")[1]
            arm = assemble.esterify("[CH3:2]CC(=O)O", alcohol,
                                    label=f"arm {t['name']}")
            parent, n_arms = build_parent(
                h["smiles"], arm, int(h["n_tails_computed"]), label=name)
            cation, n_me = _quaternise(parent, label=name)

            pm = Chem.MolFromSmiles(parent)
            cm = Chem.MolFromSmiles(cation)
            flag, why = TAIL_CAVEATS.get(t["name"], ("", ""))
            rows.append({
                "IL_name": name,
                "headgroup": h["name"], "tail": t["name"],
                "n_tails": n_arms,
                "n_nitrogens": int(h["n_nitrogens"]),
                "n_methyl_added": n_me,
                "net_charge_methylated": Chem.GetFormalCharge(cm),
                "counterion": f"{n_me} x iodide (CH3I quaternisation)",
                # `IL_SMILES` is the column name C10 (provenance/
                # structure_check.py) looks for, and it is set to the
                # METHYLATED CATION on purpose: that is the species the assay
                # was run on and the species Figure 1A draws in full, so it is
                # what a C10 anchor can be read against. The Michael-adduct
                # parent is kept beside it because DECISIONS.md D-2 -- which of
                # the two is "the lipid" if this paper is ever merged -- is
                # OPEN, and because it is the parent, not the cation, that the
                # corpus already contains for 11 members.
                "IL_SMILES": cation,
                "parent_SMILES": parent,
                "parent_formula": rdMolDescriptors.CalcMolFormula(pm),
                "parent_MW": round(Descriptors.MolWt(pm), 3),
                "methylated_cation_SMILES": cation,
                "methylated_cation_formula": rdMolDescriptors.CalcMolFormula(cm),
                "methylated_cation_MW": round(Descriptors.MolWt(cm), 3),
                "parent_InChIKey": Chem.MolToInchiKey(pm),
                "methylated_InChIKey": Chem.MolToInchiKey(cm),
                "structure_caveat_kind": flag,
                "structure_caveat_why": why,
                "source_figure": "Figure 1A (blocks + the drawn methylated "
                                 "hexa-substituted product)",
            })

    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    uniq_p = len({r["parent_InChIKey"] for r in rows})
    uniq_m = len({r["methylated_InChIKey"] for r in rows})
    print(f"wrote {len(rows)} members -> {OUT}")
    print(f"  unique parent InChIKeys     {uniq_p}/{len(rows)}")
    print(f"  unique methylated InChIKeys {uniq_m}/{len(rows)}")
    print(f"  arms per member             "
          f"{sorted({r['n_tails'] for r in rows})}")
    print(f"  charge on methylated cation "
          f"{sorted({r['net_charge_methylated'] for r in rows})}")
    print(f"  members carrying a structure caveat "
          f"{sum(1 for r in rows if r['structure_caveat_kind'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
