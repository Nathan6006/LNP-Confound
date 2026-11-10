#!/usr/bin/env python3
"""Enumerate peptoids 222-266 from SI Table S9's sequences. Writes 03_enumerated.csv.

Run from this folder:   python3 -I stage03_enumerate.py

HOW A PEPTOID IS BUILT HERE. N-terminus to C-terminus, one residue at a time,
each join made by curation/lib assemble.attach at ATOM-MAPPED attachment atoms.
No RWMol.RemoveAtom anywhere: nothing is deleted, so the F-013 mechanism (an
atom dropped out of a stereocentre's neighbour list without correcting its
chiral tag, 537 rows of epicholesterol) cannot occur. attach also asserts
no-radicals and stereo preservation on every single build.

    N-terminal residue   N(R)C[C:2]=O          nitrogen keeps one implicit H
    internal residue     [N:1](R)C[C:2]=O      nitrogen becomes tertiary
    C-terminal residue   [N:1](R)CC(N)=O       Rink amide primary carboxamide

and for Dea / Eae, whose submonomer amine is SECONDARY, the N-terminal residue
is N(A1)(A2)C[C:2]=O with no hydrogen at all. See blocks.py.

THE CHECK THAT MAKES THIS WORTH ANYTHING. Every one of the 45 peptoids has a
CALCULATED EXACT MASS printed in SI Table S9. The enumerated structure's
monoisotopic mass is compared against it. This is an external anchor: the
authors computed it from their intended structure, we compute it from ours, and
they were produced by different people from different objects.

WHAT THE MASS CHECK CANNOT SEE, stated so nobody reads a clean run as more than
it is: the ORDER of the monomers. Peptoids 236, 237, 259, 264 and 266 ALL have
calculated mass 1275.09 and all five are genuinely different molecules -- 266
is 236's deliberate sequence isomer, and the paper made it precisely to show
that composition rather than order drives performance. So mass confirms
composition and nothing else. Order is confirmed separately, two ways: SI Table
S4 independently prints the same sequences (34/34 agree), and SI Figure S7
DRAWS 236 and 266 in full, which is the C10 anchor in make_reference.py.

PEPTOID 231 IS NOT ENUMERATED. Its Dme monomer is not in Figure 1B's drawn pool
and C4H10N has several connectivities; see blocks.py. Excluded rather than
guessed.
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)

import assemble                                    # noqa: E402
from rdkit import Chem, RDLogger                   # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt      # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula  # noqa: E402

import blocks                                      # noqa: E402
import seqtables                                   # noqa: E402

RDLogger.DisableLog("rdApp.*")

# The paper prints the calculated exact mass to 2 decimals, so agreement is
# asserted at half a unit in the last printed place. A one-hydrogen error is
# 1.008 Da and a one-CH2 error is 14.016 Da; both are ~20x and ~2800x this
# tolerance, so the check has no room to pass a real structural error.
MASS_TOL = 0.005

SKIP = {231: "Dme monomer connectivity is not determined by any source we hold"}


def residue(code, position, n_residues):
    """SMILES for one residue, with the maps attach() needs."""
    m = blocks.MONOMERS[code]
    first, last = position == 0, position == n_residues - 1
    if first:
        if m["amine_order"] == "secondary":
            a1, a2 = m["arms"]
            head = f"N({a1})({a2})"
        else:
            head = f"N({m['sidechain']})"
    else:
        if m["amine_order"] == "secondary":
            raise ValueError(
                f"{code} is a secondary amine and can only be the N-terminal "
                f"residue; it was asked for at position {position}. A "
                "secondary-amine residue's backbone nitrogen already carries "
                "two arms, so it cannot be bromoacetylated again.")
        head = f"[N:1]({m['sidechain']})"
    # The carbonyl is written [CH:2]=O, not [C:2]=O. A bracket atom carries
    # noImplicit=True, so an UNBONDED [C:2] with only CH2 and =O sanitises to a
    # RADICAL -- and this fragment's map-2 carbon is deliberately left unbonded
    # for one iteration, so attach's own clearing (which touches only the two
    # atoms it joins) never reaches it. assemble.check_no_radicals caught this
    # on the first build, which is exactly the Siegwart radical-carbon class:
    # it would have parsed, canonicalised and given self-consistent descriptors.
    tail = "CC(N)=O" if last else "C[CH:2]=O"
    return head + tail


def build(monomers, label=""):
    """Canonical SMILES for a peptoid given its monomer sequence, N -> C."""
    n = len(monomers)
    if n < 2:
        raise ValueError(f"{label}: {n} residues is not a peptoid chain.")
    chain = residue(monomers[0], 0, n)
    for i in range(1, n):
        nxt = residue(monomers[i], i, n)
        # attach(a, b): a's map-1 atom bonds to b's map-2 atom. Here a is the
        # incoming residue's nitrogen and b is the growing chain's carbonyl
        # carbon, so the chain extends toward the C-terminus.
        chain = assemble.attach(nxt, chain, a_map=1, b_map=2,
                                label=f"{label} residue {i + 1}/{n}")
    if ":1]" in chain or ":2]" in chain:
        raise ValueError(f"{label}: unconsumed atom map in {chain}")
    return assemble.canonical(chain)


def main():
    s9, order = seqtables.table_s9()
    s4, _ = seqtables.table_s4()

    rows, failures = [], []
    for pid in order:
        rec = s9[pid]
        mons = rec["monomers"]
        if pid in SKIP:
            rows.append(dict(
                peptoid=pid, sequence="-".join(mons), n_residues=len(mons),
                smiles="", formula="", exact_mass="",
                paper_calc_exact_mass=rec["calc_exact_mass"],
                mass_delta="", mass_ok="", enumerated="no",
                reason=SKIP[pid]))
            continue
        smi = build(mons, label=f"peptoid {pid}")
        mol = Chem.MolFromSmiles(smi)
        mw = ExactMolWt(mol)
        delta = mw - rec["calc_exact_mass"]
        ok = abs(delta) <= MASS_TOL
        if not ok:
            failures.append((pid, "-".join(mons), round(mw, 4),
                             rec["calc_exact_mass"], round(delta, 4)))
        rows.append(dict(
            peptoid=pid, sequence="-".join(mons), n_residues=len(mons),
            smiles=smi, formula=CalcMolFormula(mol), exact_mass=round(mw, 4),
            paper_calc_exact_mass=rec["calc_exact_mass"],
            mass_delta=round(delta, 4), mass_ok="yes" if ok else "NO",
            enumerated="yes", reason=""))

    out = os.path.join(HERE, "03_enumerated.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    built = [r for r in rows if r["enumerated"] == "yes"]
    print(f"parsed from SI Table S9 : {len(s9)} peptoids (222-266)")
    print(f"cross-checked vs Table S4: {len(s4)} sequences, "
          f"{sum(1 for p in s4 if s4[p]['monomers'] == s9[p]['monomers'])} agree")
    print(f"enumerated              : {len(built)}")
    print(f"skipped                 : {len(rows) - len(built)} "
          f"({', '.join(str(p) for p in SKIP)})")
    print(f"exact-mass agreement    : "
          f"{sum(1 for r in built if r['mass_ok'] == 'yes')}/{len(built)} "
          f"within {MASS_TOL} Da")
    print(f"unique canonical SMILES : {len(set(r['smiles'] for r in built))}")
    if failures:
        print("\nMASS MISMATCHES:")
        for f in failures:
            print("   ", f)
        return 1
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
