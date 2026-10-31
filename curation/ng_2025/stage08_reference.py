"""C10 anchors: the structures this paper DRAWS, built by a different code path.

Six anchors, one per chemistry plus both ketone ring systems. Each is assembled
from named fragments with RDKit's RWMol/CombineMols here, whereas the
enumeration in stage04 builds its products by attaching arms to the amine or by
running a reaction template - so a shared bug would have to occur twice, in two
different ways, to pass unnoticed.

B-a1 is the strongest of the six: the paper draws that exact product in main
Fig. 1d, which is VECTOR line art, so its reference string is DERIVED HERE from
the page's own bond coordinates (01_fig1d_graph.json, written by
stage01_vector.py) rather than typed by a human who has already seen the
enumerated answer. It used to be a string literal with a note claiming the same
thing, which is the "anchor that can only pass" failure ERROR_CLASSES warns
about; vector_smiles() below is the fix, so the claim is now executed.

That panel needs the glyph merge radius at 3.0 pt, not the 3.6 used for
Fig. 1b/1c: the product is drawn compactly and at 3.6 two separate bond
endpoints fall inside one "HO" label box, giving that oxygen three bonds and an
AtomValenceException. The radius is an argument to stage01_vector.py for
exactly this reason.
"""
import csv, json, os
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

D = os.path.dirname(os.path.abspath(__file__))


def vector_smiles(graph_path, min_atoms):
    """Largest structure in a stage01_vector.py graph, as SMILES.

    Same element rule as stage02_blocks.py: a label's path signature gives the
    atom (curves only, or the 12-line+20-curve "HO", is oxygen; anything else
    with straight strokes is nitrogen).
    """
    g = json.load(open(graph_path))
    nc = g['n_carbon_vertices']

    def element(c):
        return 'O' if (c[5] and not c[4]) or (c[4] == 12 and c[5] == 20) else 'N'

    comp = max((c for c in g['components'] if len(c) >= min_atoms), key=len)
    idx, m = {}, Chem.RWMol()
    for v in comp:
        idx[v] = m.AddAtom(Chem.Atom(element(g['glyphs'][v - nc]) if v >= nc else 'C'))
    for u, w, o in g['bonds']:
        if u in idx and w in idx:
            m.AddBond(idx[u], idx[w],
                      Chem.BondType.DOUBLE if o == 2 else Chem.BondType.SINGLE)
    mol = m.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def join(core, attachments):
    """core + [(core_atom_idx, fragment_smiles), ...] -> SMILES.

    Each fragment carries a single dummy atom '*' marking where it bonds, so
    the attachment point is written in the fragment itself and never as a bare
    index that could silently drift when a fragment is edited.
    """
    m = Chem.RWMol(Chem.MolFromSmiles(core))
    for ci, frag in attachments:
        f = Chem.MolFromSmiles(frag)
        star = [a.GetIdx() for a in f.GetAtoms() if a.GetAtomicNum() == 0]
        assert len(star) == 1, frag
        anchor = f.GetAtomWithIdx(star[0]).GetNeighbors()[0].GetIdx()
        off = m.GetNumAtoms()
        m.InsertMol(f)
        m.AddBond(ci, off + anchor, Chem.BondType.SINGLE)
        m.RemoveAtom(off + star[0])
    mol = m.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


def amine_n(smiles):
    """index of the single primary-amine nitrogen."""
    m = Chem.MolFromSmiles(smiles)
    n = [a.GetIdx() for a in m.GetAtoms()
         if a.GetSymbol() == 'N' and a.GetTotalNumHs() == 2 and a.GetDegree() == 1]
    assert len(n) == 1, smiles
    return n[0]


# --- fragments, named as the paper names them -------------------------------
HYDROXY_OCTYLBENZYL = 'Oc1ccc(CCCCCCCC)cc1C*'    # 2-hydroxy-5-octylbenzyl
PYRROLYLMETHYL_M = 'Cn1cccc1C*'                  # (1-methylpyrrol-2-yl)methyl
PHENACYLETHYL_P = 'O=C(c1ccccc1)CC*'             # 3-oxo-3-phenylpropyl

ANCHORS = []

# C-a16: cyclo(Lys-Lys) core, four 2-hydroxy-5-octylbenzyl arms on the two
# side-chain nitrogens. Drawn in Supplementary Figs 3a and 19a (the inset of
# the mass spectrum), whose own caption reads "Molecular Weight 1129.71".
cyclo_lys = 'NCCCCC1NC(=O)C(CCCCN)NC1=O'
m = Chem.MolFromSmiles(cyclo_lys)
prim = [a.GetIdx() for a in m.GetAtoms()
        if a.GetSymbol() == 'N' and a.GetTotalNumHs() == 2 and a.GetDegree() == 1]
ANCHORS.append(('C-a16',
                join(cyclo_lys, [(n, HYDROXY_OCTYLBENZYL) for n in prim for _ in (0, 1)]),
                'Supplementary Fig 3a; Supplementary Fig 19a (MS inset)',
                'four 2-hydroxy-5-octylbenzyl arms on the two lysine side-chain NH2 of '
                'cyclo(Lys-Lys); the diketopiperazine N-H stay intact'))

# M-b1 and P-b1: methylamine plus two arms, the generic Mannich product
ANCHORS.append(('M-b1', join('CN', [(amine_n('CN'), PYRROLYLMETHYL_M)] * 2),
                'Supplementary Fig 14a',
                'N-methyl, two (1-methylpyrrol-2-yl)methyl arms'))
ANCHORS.append(('P-b1', join('CN', [(amine_n('CN'), PHENACYLETHYL_P)] * 2),
                'Supplementary Fig 16a',
                'N-methyl, two 3-oxo-3-phenylpropyl arms (the acetophenone methyl '
                'becomes the CH2 beta to nitrogen)'))

# the two ketone ring systems, both with methylamine and the diethyl ester
ANCHORS.append(('I-a-b1', Chem.CanonSmiles('O=C1C(C(=O)OCC)CN(C)CC1C(=O)OCC'),
                'Supplementary Fig 8a',
                'diethyl 1-methyl-4-oxopiperidine-3,5-dicarboxylate: the 1:2:1 double '
                'Mannich closes a ring, it is not an acyclic bis-adduct'))
ANCHORS.append(('I-b-b1', Chem.CanonSmiles('O=C1C2(C(=O)OCC)CN(C)CC1(C(=O)OCC)CN(C)C2'),
                'Supplementary Fig 11a',
                'diethyl 3,7-dimethyl-9-oxo-3,7-diazabicyclo[3.3.1]nonane-1,5-'
                'dicarboxylate (bispidinone): the 1:4:2 ratio closes a second ring'))

# B-a1, derived at run time from the vector line art of main Fig. 1d
ANCHORS.append(('B-a1', vector_smiles(os.path.join(D, '01_fig1d_graph.json'), 20),
                'Fig. 1d',
                'DERIVED at run time from the page vector bond coordinates in '
                '01_fig1d_graph.json (stage01_vector.py, glyph merge radius 3.0): glyph '
                'signatures give the atoms, geometry gives the chain lengths. No human '
                'typed this string, so it cannot have been copied from the enumeration'))


def main():
    enum = {r['lipid_name']: r['IL_SMILES'] for r in csv.DictReader(open(os.path.join(D, '03_enumerated.csv')))}
    out = []
    for name, smi, fig, note in ANCHORS:
        out.append(dict(lipid_name=name, reference_SMILES=Chem.CanonSmiles(smi),
                        source_figure=fig, read_on='2026-10-29', read_by='claude-code',
                        expected_mismatch='', note=note))
        same = Chem.CanonSmiles(enum[name]) == Chem.CanonSmiles(smi)
        print(f"{'OK  ' if same else 'DIFF'} {name:8} {fig}")
        if not same:
            print('   ref :', Chem.CanonSmiles(smi))
            print('   enum:', Chem.CanonSmiles(enum[name]))
    with open(os.path.join(D, 'reference_structures.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)


if __name__ == '__main__':
    main()
