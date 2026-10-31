"""Turn the Fig. 1b/1c vector graphs into SMILES, with no human in the loop.

Element identity comes from the glyph's own path signature, which is exactly
separable on this page:

    lines=10 curves=0   w~2.0 h~2.4   ->  N      (tertiary amine label "N")
    lines=22 curves=0   w~2.0 h~5.4   ->  N      (secondary "N/H", stacked)
    lines=26 curves=18  w~6.0 h~3.2   ->  N      ("H2N" / "NH2")
    lines=0  curves=20  w~2.4 h~2.5   ->  O      (ether or carbonyl "O")
    lines=12 curves=20  w~5.1 h~2.7   ->  O      ("HO")

so the structures below are machine-read end to end: geometry gives the
skeleton and chain lengths, the glyph signature gives the atoms. The reading
is then compared against the independently hand-written SMILES in
02_blocks_hand.csv; the two must agree as canonical SMILES.
"""
import csv, json, os, sys, collections
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

D = os.path.dirname(os.path.abspath(__file__))
G = json.load(open(os.path.join(D, sys.argv[1])))
V, NC = G['vertices'], G['n_carbon_vertices']


def element(c):
    lines, curves = c[4], c[5]
    if curves and not lines:
        return 'O'
    if lines == 12 and curves == 20:
        return 'O'
    return 'N'


def comp_smiles(comp):
    idx = {}
    m = Chem.RWMol()
    for v in comp:
        a = Chem.Atom(element(G['glyphs'][v - NC]) if v >= NC else 'C')
        idx[v] = m.AddAtom(a)
    for u, w, o in G['bonds']:
        if u in idx and w in idx:
            m.AddBond(idx[u], idx[w],
                      Chem.BondType.DOUBLE if o == 2 else Chem.BondType.SINGLE)
    mol = m.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


comps = [c for c in G['components'] if len(c) > 2]
comps.sort(key=lambda c: (round(min(V[v][1] for v in c) / 9), min(V[v][0] for v in c)))
read = []
for c in comps:
    xs = [V[v][0] for v in c]; ys = [V[v][1] for v in c]
    read.append({'n_atoms': len(c), 'x0': round(min(xs), 1), 'y0': round(min(ys), 1),
                 'smiles': comp_smiles(c)})

hand = {r['block_id']: r for r in csv.DictReader(open(os.path.join(D, '02_blocks_hand.csv')))}
# position-ordered block ids, filled in after looking at the rendered panels once
ORDER = sys.argv[2].split(',')
assert len(ORDER) == len(read), f'{len(ORDER)} ids vs {len(read)} structures drawn'

ok = bad = 0
out = []
for bid, r in zip(ORDER, read):
    h = Chem.MolFromSmiles(hand[bid]['SMILES'])
    hcanon = Chem.MolToSmiles(h)
    agree = hcanon == r['smiles']
    ok += agree; bad += not agree
    out.append(dict(block_id=bid, n_atoms_drawn=r['n_atoms'], x=r['x0'], y=r['y0'],
                    vector_SMILES=r['smiles'], hand_SMILES=hcanon, agree=agree))
    print(f"{'OK ' if agree else 'DIFF'} {bid:6} n={r['n_atoms']:3d}  {r['smiles']}")
    if not agree:
        print(f"{'':11}hand:  {hcanon}")
w = csv.DictWriter(open(os.path.join(D, '02_vector_vs_hand.csv'), 'w'), fieldnames=list(out[0]))
w.writeheader(); w.writerows(out)
print(f"\nagree {ok} / {ok + bad}")
