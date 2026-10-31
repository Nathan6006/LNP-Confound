"""Per-row confidence flags. The corpus schema has no confidence column, so the
flags live beside the rows and nothing is merged until they are resolved.

Two issues, both structural, both traced to the paper contradicting itself:

 1. b3 and b4. The four amine-core panels (Supplementary Figs 7c, 10c, 13c,
    15c) all draw one CH2 between the primary nitrogen and the tertiary one,
    i.e. H2N-CH2-NMe2 and H2N-CH2-(4-methylpiperazin-1-yl), which are aminals.
    The pyrrole and acetophenone PRODUCT panels (Supplementary Figs 14 and 16)
    draw the same two amines with a direct N-N bond, i.e. 1,1-dimethylhydrazine
    and 1-amino-4-methylpiperazine - both ordinary catalogue reagents, whereas
    the aminals are not isolable. The ketone product panels (Figs 8, 9, 11, 12)
    side with the amine-core panels. Emitted AS DRAWN in the amine-core panels,
    because those panels are the paper's statement of which reagents were used
    and all four agree; the alternative SMILES is carried here so the choice is
    one line to reverse.

 2. the two dialkyl-ketone libraries. Enumerated from the reaction schemes
    (Supplementary Figs 7a and 10a), which draw the 1:2:1 product as the cyclic
    4-oxopiperidine-3,5-dicarboxylate and the 1:4:2 product as the bicyclic
    bispidinone. The per-lipid drawings do not all agree with their own scheme:
    of tail I's 13 lipids at 1:2:1, five are drawn as the open-chain
    mono-adduct, while tail J's are drawn cyclic for 12 of 13 - the same amine
    drawn both ways for two different tails, which is a drawing inconsistency
    and not chemistry. The stoichiometries themselves only make sense as the
    ring (2 CH2 + 1 amine) and the double ring (4 CH2 + 2 amines).
"""
import csv, os
from rdkit import Chem

D = os.path.dirname(os.path.abspath(__file__))
ALT = {'b3': ('NN(C)C', '1,1-dimethylhydrazine'),
       'b4': ('NN1CCN(C)CC1', '1-amino-4-methylpiperazine')}


def main():
    enum = list(csv.DictReader(open(os.path.join(D, '03_enumerated.csv'))))
    out = []
    for r in enum:
        if r['amine'] in ALT:
            alt, altname = ALT[r['amine']]
            out.append(dict(
                lipid_name=r['lipid_name'], confidence='low', issue='reagent_identity',
                detail=f"amine {r['amine']} is drawn as an aminal in the amine-core panels "
                       f"and as a hydrazine ({altname}) in Supplementary Figs 14 and 16",
                alternative_amine_SMILES=alt))
        if r['lipid_name'] == 'O-b10':
            out.append(dict(
                lipid_name='O-b10', confidence='low', issue='value_disagrees_between_figures',
                detail='Fig. 1i Source Data give log10 = 4.6097 (4.07e4), the highest value '
                       'in the whole pyrrole panel; Supplementary Fig 17b, which re-measures '
                       'six library members individually, plots this lipid at ~1.7e3 - about '
                       '24x lower. The other five lipids in Supp Fig 17b agree with Fig. 1 to '
                       'within 1.9x, so this is not generic run-to-run spread. Which figure '
                       'is wrong cannot be decided from the paper.',
                alternative_amine_SMILES=''))
        if r['library'].startswith('ketone'):
            out.append(dict(
                lipid_name=r['lipid_name'], confidence='low', issue='product_ring_system',
                detail='enumerated from the reaction scheme (cyclic for 1:2:1, bicyclic for '
                       '1:4:2); a minority of the per-lipid drawings show the open-chain '
                       'mono-adduct instead',
                alternative_amine_SMILES=''))
    with open(os.path.join(D, '07_flags.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader(); w.writerows(out)
    names = {r['lipid_name'] for r in out}
    print(f'flagged rows: {len(out)} covering {len(names)} of {len(enum)} lipids')
    for issue in ('reagent_identity', 'product_ring_system'):
        n = {r['lipid_name'] for r in out if r['issue'] == issue}
        print(f'  {issue:22} {len(n):3d} lipids')
    both = {r['lipid_name'] for r in out if r['issue'] == 'reagent_identity'} & \
           {r['lipid_name'] for r in out if r['issue'] == 'product_ring_system'}
    print(f'  both issues            {len(both):3d} lipids')
    print(f'  unflagged (confidence high) {len(enum) - len(names):3d} lipids')


if __name__ == '__main__':
    main()
