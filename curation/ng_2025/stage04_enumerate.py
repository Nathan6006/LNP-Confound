"""Enumerate the 310 Mannich lipids of the six libraries.

Building blocks come from 02_blocks_hand.csv (Fig 1b/1c, machine-verified
against the page's vector line art by stage02_blocks.py) and from the
Supplementary rasters, where chain lengths were counted by stage03_zigzag.py
rather than by eye.

Degree of substitution is NOT one rule. It differs by chemistry, and the two
dialkyl-ketone libraries are the reason this file exists rather than a one-line
product template:

  phenol / thiol / pyrrole / acetophenone
      every PRIMARY NH2 becomes a tertiary N bearing two tail-derived arms.
      Secondary amine N-H, tertiary N and the a16 amide N-H do not react
      (drawn intact in Supplementary Figs 2-3 and 5-6; C-a16's 1H NMR
      integrates the N-CH2-Ar multiplet for 8H = 4 arms; Methods run C-a16 at
      4.5 phenol : 4.2 HCHO : 1.0 amine).

  dialkyl ketone 1:2:1  (Supplementary Fig 7a)
      a DOUBLE Mannich that closes a ring: the two alpha carbons of the
      dialkyl 3-oxoglutarate each take one CH2 and the amine nitrogen bridges
      them, giving 1-R-4-oxopiperidine-3,5-dicarboxylate. One amine per lipid.

  dialkyl ketone 1:4:2  (Supplementary Fig 10a)
      the same reaction run twice over, giving the bicyclic bispidinone
      3,7-diaza-bicyclo[3.3.1]nonan-9-one-1,5-dicarboxylate. TWO amine
      residues per lipid, both the same amine.

The ketone products are therefore not R-N(CH2-tail)2 at all, and the
stoichiometries printed in the schemes (1:2:1 and 1:4:2) are exactly the ring
and double-ring closures. See 04_notes.md for the drawing inconsistency that
makes these two libraries confidence=low.
"""
import csv, os, itertools, json
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
RDLogger.DisableLog('rdApp.*')

D = os.path.dirname(os.path.abspath(__file__))

# --- tails, as the ARM that ends up bonded to nitrogen (atom 0 = that bond) ---
ARM = {
    'A': 'Cc1ccccc1O',                    # 2-hydroxybenzyl        (phenol)
    'B': 'Cc1cc(CCCCCC)ccc1O',            # 4-hexylphenol,  ortho
    'C': 'Cc1cc(CCCCCCCC)ccc1O',          # 4-octylphenol,  ortho
    'D': 'Cc1cc(CCCCCCCCC)ccc1O',         # 4-nonylphenol,  ortho
    'E': 'CSCCCCCC',                      # hexanethiol            (thiol)
    'F': 'CSCCCCCCCC',                    # octanethiol
    'G': 'CSCCCCCCCCCC',                  # decanethiol
    'H': 'CSCCCCCCCCCCCC',                # dodecanethiol
    'M': 'Cc1cccn1C',                     # 1-methylpyrrole, C2    (pyrrole)
    'N': 'Cc1cccn1CC',                    # 1-ethylpyrrole,  C2
    'O': 'Cc1cccn1CCCCCCCC',              # 1-octylpyrrole,  C2
    'P': 'CCC(=O)c1ccccc1',               # acetophenone           (aceto)
    'Q': 'CCC(=O)c1ccc(CCCCCC)cc1',       # 4'-hexylacetophenone
    'R': 'CCC(=O)c1ccc(CCCCCCCC)cc1',     # 4'-octylacetophenone
}
# dialkyl 3-oxoglutarate esters: the O-alkyl of both esters
ESTER = {'I': 'CC', 'J': 'CCCC', 'K': 'CCCCCC', 'L': 'CCCCCCCC'}

LIB = {
    'phenol':       dict(tails=list('ABCD'), amines='a', fig='Fig. 1e'),
    'thiol':        dict(tails=list('EFGH'), amines='a', fig='Fig. 1f'),
    'ketone_1_2_1': dict(tails=list('IJKL'), amines='b', fig='Fig. 1g'),
    'ketone_1_4_2': dict(tails=list('IJKL'), amines='b', fig='Fig. 1h'),
    'pyrrole':      dict(tails=list('MNO'),  amines='b', fig='Fig. 1i'),
    'acetophenone': dict(tails=list('PQR'),  amines='b', fig='Fig. 1j'),
}


def dialkylate(amine, arm):
    """Two tail arms onto every primary NH2; everything else untouched."""
    m = Chem.RWMol(Chem.MolFromSmiles(amine))
    prim = [a.GetIdx() for a in m.GetAtoms()
            if a.GetSymbol() == 'N' and a.GetTotalNumHs() == 2 and a.GetDegree() == 1]
    for n in prim:
        for _ in range(2):
            off = m.GetNumAtoms()
            m.InsertMol(Chem.MolFromSmiles(arm))
            m.AddBond(n, off, Chem.BondType.SINGLE)
    mol = m.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol), len(prim)


def piperidinone(amine, ester):
    rxn = AllChem.ReactionFromSmarts(
        f'[NX3;H2:1]>>[N:1]1CC(C(=O)O{ester})C(=O)C(C(=O)O{ester})C1')
    p = rxn.RunReactants((Chem.MolFromSmiles(amine),))
    mol = p[0][0]
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def bispidinone(amine, ester):
    rxn = AllChem.ReactionFromSmarts(
        f'[NX3;H2:1].[NX3;H2:2]>>'
        f'[N:1]1CC2(C(=O)O{ester})C(=O)C(C(=O)O{ester})(C1)C[N:2]C2')
    a = Chem.MolFromSmiles(amine)
    p = rxn.RunReactants((a, a))
    mol = p[0][0]
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def main():
    blocks = {r['block_id']: r['SMILES'] for r in csv.DictReader(open(os.path.join(D, '02_blocks_hand.csv')))}
    amines = {'a': [f'a{i}' for i in range(1, 17)],
              'b': [f'b{i}' for i in range(1, 14)]}
    rows = []
    for lib, cfg in LIB.items():
        for tail in cfg['tails']:
            for am in amines[cfg['amines']]:
                asmi = blocks['H_' + am]
                if lib == 'ketone_1_2_1':
                    name = f'{tail}-a-{am}'
                    alias = f'{tail}1-{am}'
                    smi = piperidinone(asmi, ESTER[tail])
                    narm = 1
                elif lib == 'ketone_1_4_2':
                    name = f'{tail}-b-{am}'
                    alias = f'{tail}2-{am}'
                    smi = bispidinone(asmi, ESTER[tail])
                    narm = 2
                else:
                    name = f'{tail}-{am}'
                    alias = name
                    smi, nprim = dialkylate(asmi, ARM[tail])
                    narm = 2 * nprim
                rows.append(dict(lipid_name=name, heatmap_label=alias, library=lib,
                                 tail=tail, amine=am, n_tail_arms=narm,
                                 IL_SMILES=smi, source_figure=cfg['fig']))
    with open(os.path.join(D, '03_enumerated.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    smis = [r['IL_SMILES'] for r in rows]
    print(f'enumerated {len(rows)} lipids   unique SMILES {len(set(smis))}')
    for lib in LIB:
        n = sum(r['library'] == lib for r in rows)
        print(f'  {lib:14} {n:3d}')
    dup = {s for s in smis if smis.count(s) > 1}
    print('duplicate structures:', len(dup))
    for s in sorted(dup):
        print('   ', [r['lipid_name'] for r in rows if r['IL_SMILES'] == s], s[:70])


if __name__ == '__main__':
    main()
