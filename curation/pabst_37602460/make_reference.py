#!/usr/bin/env python3
"""The C10 anchor GENERATOR for pabst_37602460 -- not a stored string.

Run from this folder:   python3 -I make_reference.py

WHAT THE PAPER DRAWS. Figure 1A is a three-panel scheme. Its top row draws the
reaction in full: tris(2-aminoethyl)amine, plus six acrylates, giving the
hexa-substituted adduct, then CH3I giving a product in which ALL FOUR
nitrogens -- the three that reacted and the central tertiary one that did not
-- carry a methyl and a '+'. Its lower panels draw each of the eight
polyamines and each of the thirteen R groups as a full structure. So the
figure does state a complete connectivity for any (headgroup, tail) pair; it
simply does so in two pieces.

WHY THESE TWO ANCHORS.
  3NB8   is the member whose product the scheme itself draws. Its headgroup is
         the TREN in the top row and its tail is the n-butyl of 'R = 8', so
         every bond below was read off a drawn bond.
  4NC6B  exercises what 3NB8 cannot: a RING headgroup (cyclam, where an arm
         and a methyl go onto a nitrogen that already has two ring bonds) and a
         BRANCHED tail (tert-butyl). A single acyclic anchor would leave the
         cyclic headgroup untested, and the cyclic headgroup is the one whose
         identity SI Table S5 gets wrong (08_contradictions.csv C-2).

WHY THIS IS A DIFFERENT CODE PATH, STATED AS SOMETHING EXECUTED.
CLAUDE.md section 6 and PIPELINE_V2 section 3 both say the independence has to
be executed, not asserted, because it was asserted twice and false twice:

  F-013  the ZT_2026 anchor was a BYTE COPY of the enumeration -- same sha256,
         223 bytes -- while its note claimed independent fragment assembly
         "through the C3 oxygen (3-beta)" and the string was 3-alpha.
  F-037  the generator written to cure F-013 hard-coded an AMINE_HALF template
         that was the enumeration's canonical SMILES with the sterol swapped
         for '%' and one ring digit renumbered: byte-identical over 141
         characters.

So this file contains NO SMILES STRING AT ALL -- not a template, not a
fragment, not a placeholder. There is nothing here for a copy to hide in. The
structures are written as explicit atom lists and 0-based bond tables, in the
order a reader walks the drawing: the headgroup first, then one arm at a time
away from each nitrogen, then the methyls. The enumerator
(stage03_enumerate.py) shares none of this: it starts from SMILES for the
amine and the alcohol, makes the ester with assemble.esterify, attaches arms
with assemble.attach by atom map, and quaternises last. A bond table also has
no ring-closure digits, which is the mechanism that once fused piperazine to
cholesterol.

WHAT AN AGREEMENT HERE DOES AND DOES NOT PROVE. It proves the two routes agree
on connectivity, atom count and charge. It does NOT prove the reading of
Figure 1A was right -- that is what the printed cross-checks in 07_qc.json are
for (tail MW against Table S4, headgroup N/MW against Table S5, and the eleven
members whose parent already sits in data/data_vitro.csv under KW_2014,
AA_2008 and NC_2024 with a matching FULL InChIKey).
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

import assemble                                        # noqa: E402
from rdkit import Chem                                 # noqa: E402
from rdkit import RDLogger                             # noqa: E402
RDLogger.DisableLog("rdApp.*")

OUT = os.path.join(HERE, "reference_structures.csv")

SINGLE = Chem.BondType.SINGLE
DOUBLE = Chem.BondType.DOUBLE


class Draw:
    """A drawing pad: add an atom, bond it to one already there. No SMILES."""

    def __init__(self):
        self.rw = Chem.RWMol()
        self.charges = []

    def atom(self, element, bond_to=None, order=SINGLE):
        i = self.rw.AddAtom(Chem.Atom(element))
        if bond_to is not None:
            self.rw.AddBond(bond_to, i, order)
        return i

    def bond(self, i, j, order=SINGLE):
        self.rw.AddBond(i, j, order)

    def quaternary(self, i):
        """Mark a nitrogen as the figure draws it: a methyl and a '+'."""
        self.atom(6, bond_to=i)
        self.charges.append(i)

    def finish(self, label):
        for i in self.charges:
            self.rw.GetAtomWithIdx(i).SetFormalCharge(1)
            self.rw.GetAtomWithIdx(i).SetNoImplicit(False)
        m = self.rw.GetMol()
        Chem.SanitizeMol(m)
        assemble.check_no_radicals(m, label=label)
        return m


def arm(d, n_idx, r_shape):
    """One drawn arm: N-CH2-CH2-C(=O)-O-R, with R given as a shape callback.

    Read off the figure's product panel, left to right along the pink
    highlight: the nitrogen, two methylenes, the carbonyl carbon with its
    doubly bonded oxygen above it, the ester oxygen, then the R group.
    """
    c_alpha = d.atom(6, bond_to=n_idx)
    c_beta = d.atom(6, bond_to=c_alpha)
    c_carbonyl = d.atom(6, bond_to=c_beta)
    d.atom(8, bond_to=c_carbonyl, order=DOUBLE)
    o_ester = d.atom(8, bond_to=c_carbonyl)
    r_shape(d, o_ester)


def r_n_butyl(d, o):
    """'R = 8' in Figure 1A: a four-vertex zigzag, no branch."""
    a = d.atom(6, bond_to=o)
    b = d.atom(6, bond_to=a)
    c = d.atom(6, bond_to=b)
    d.atom(6, bond_to=c)


def r_n_pentyl(d, o):
    """'R = 9' in Figure 1A, and the arms drawn in SI Figure S1 panel (E)."""
    prev = o
    for _ in range(5):
        prev = d.atom(6, bond_to=prev)


def r_tert_butyl(d, o):
    """'6B' in Figure 1A: one bond to a carbon with three methyls on it."""
    q = d.atom(6, bond_to=o)
    d.atom(6, bond_to=q)
    d.atom(6, bond_to=q)
    d.atom(6, bond_to=q)


def draw_3NB8():
    """TREN hexa-substituted and fully methylated, exactly as the scheme draws it."""
    d = Draw()
    n_centre = d.atom(7)
    outer = []
    for _ in range(3):                       # three 2-carbon arms off the centre
        c1 = d.atom(6, bond_to=n_centre)
        c2 = d.atom(6, bond_to=c1)
        outer.append(d.atom(7, bond_to=c2))
    for n in outer:                          # two tails on each outer nitrogen
        arm(d, n, r_n_butyl)
        arm(d, n, r_n_butyl)
    for n in [n_centre] + outer:             # every nitrogen: a methyl and a +
        d.quaternary(n)
    return d.finish("3NB8 anchor")


def draw_4NC6B():
    """Cyclam, one tail on each ring nitrogen, fully methylated.

    The ring is walked the way Figure 1A draws it: starting at the top-left NH
    and going clockwise, the spacers between nitrogens are 3, 2, 3, 2 carbons.
    """
    d = Draw()
    spacers = [3, 2, 3, 2]
    first_n = d.atom(7)
    nitrogens = [first_n]
    prev = first_n
    for k, n_carbons in enumerate(spacers):
        for _ in range(n_carbons):
            prev = d.atom(6, bond_to=prev)
        if k < len(spacers) - 1:
            prev = d.atom(7, bond_to=prev)
            nitrogens.append(prev)
        else:
            d.bond(prev, first_n)            # close the 14-membered ring
    for n in nitrogens:
        arm(d, n, r_tert_butyl)
    for n in nitrogens:
        d.quaternary(n)
    return d.finish("4NC6B anchor")


def draw_4NC9():
    """Cyclam with four n-pentyl propanoate arms, from SI Figure S1 panel (E).

    Panel (E) draws the macrocycle with its headgroup protons labelled, and the
    labels themselves give the spacers away: between the two upper nitrogens
    the sequence is d, b', d -- three carbons, the middle one distinct -- while
    down the left side it is d, d, two carbons. Alternating 3 and 2 is cyclam;
    cyclen would be 2, 2, 2, 2. The ring is walked here from that panel, one
    labelled carbon at a time, and the methylation is taken from Figure 1A.
    """
    d = Draw()
    first_n = d.atom(7)
    nitrogens = [first_n]
    prev = first_n
    for k, n_carbons in enumerate([3, 2, 3, 2]):
        for _ in range(n_carbons):
            prev = d.atom(6, bond_to=prev)
        if k < 3:
            prev = d.atom(7, bond_to=prev)
            nitrogens.append(prev)
        else:
            d.bond(prev, first_n)
    for n in nitrogens:
        arm(d, n, r_n_pentyl)
    for n in nitrogens:
        d.quaternary(n)
    return d.finish("4NC9 anchor")


ANCHORS = [
    ("3NB8", draw_3NB8,
     "Figure 1A: the reaction scheme's own product panel (tris(2-aminoethyl)"
     "amine, six acrylates, then CH3I) combined with the 'R = 8' n-butyl "
     "drawing in the Linear tail box",
     "Walked off the drawing as an explicit atom/bond table: centre N, three "
     "2-carbon arms to three outer N, two N-CH2-CH2-C(=O)-O-butyl tails on "
     "each outer N, then a methyl and a '+' on all four nitrogens. No SMILES "
     "string exists anywhere in make_reference.py, so there is nothing for a "
     "transcription to hide in (F-013, F-037)."),
    ("4NC6B", draw_4NC6B,
     "Figure 1A: the 'Cyclic' headgroup box (cyclam, ring spacers 3,2,3,2 "
     "carbons) combined with the '6B' tert-butyl drawing in the "
     "Branched/cyclic tail box",
     "Same atom/bond-table route, chosen to exercise a RING nitrogen bearing "
     "two ring bonds plus an arm plus a methyl, and a branched tail -- "
     "neither of which the 3NB8 anchor touches."),
    ("4NC9", draw_4NC9,
     "SI Figure S1 panel (E), page S5: the full structure of 4NC9 drawn with "
     "NMR assignments, whose d/b'/d and d/d proton labels state the cyclam "
     "ring spacers explicitly; methylation from Figure 1A's product panel",
     "Read from a DIFFERENT figure in a DIFFERENT document than the other two "
     "anchors, and walked the way that panel labels its carbons. SI pages "
     "S3-S8 carry no text layer, so this structure exists only as a raster "
     "and could not have been copied from anything machine-readable."),
]


def build():
    enum = {r["IL_name"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
    rows = []
    for name, fn, figure, how in ANCHORS:
        mol = fn()
        smi = Chem.MolToSmiles(mol)
        e = enum[name]
        same = (Chem.MolToInchiKey(mol)
                == Chem.MolToInchiKey(Chem.MolFromSmiles(e["IL_SMILES"])))
        rows.append([name, smi, figure, "2026-11-10",
                     "curation agent, from the PMC figure raster at "
                     "4183 x 1847 px", "",
                     how + f" | anchor InChIKey {'==' if same else '!='} the "
                     f"enumeration's | anchor formula "
                     f"{Chem.rdMolDescriptors.CalcMolFormula(mol)}"])
        print(f"{name}: {'MATCH' if same else 'MISMATCH'}  "
              f"{Chem.rdMolDescriptors.CalcMolFormula(mol)}  "
              f"{mol.GetNumAtoms()} heavy atoms")
        if not same:
            print(f"   anchor      {smi}")
            print(f"   enumeration {e['IL_SMILES']}")
    return rows


if __name__ == "__main__":
    rows = build()
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "expected_mismatch", "note"])
        w.writerows(rows)
    print(f"wrote {len(rows)} anchor(s) to {OUT}")
