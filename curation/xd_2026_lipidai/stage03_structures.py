"""E1-E9 + benchmarks read from Fig 5A by geometric vertex tracing (stage02).

CORRECTED 2026-11-09 (F-045). This docstring used to say "Every chain length is
counted from the stage02 polylines, not by eye." That overstated what this file
does. The trace WAS run -- `02_*.json` exists and its vertex counts are right --
but the twelve SMILES below are **typed literals** in `DRAWN`, transcribed from
those polylines by hand. No code in this folder reads `02_*.json`, and nothing
asserts that a typed chain length equals the traced count. So the transcription,
which is the step the old wording claimed was machine-made, was unchecked.

That is error class 24 (documented as machine-derived, actually a typed literal)
one slot over from where it was first found -- there on a C10 anchor, here on the
enumerator. It cost nothing in the end: a round-2 reviewer re-traced all twelve
from the native raster and matched every one by FULL InChIKey, stereo included.
But the claim is precisely what would have stopped a reader re-checking.

So, stated accurately: chain lengths are TRANSCRIBED from the stage02 polylines.
If this is ever regenerated, the cheap improvement is an assert that loads each
`02_*.json` and compares its per-path vertex count against the typed SMILES.

Benchmarks are the C2 positive control: corpus SMILES is the ground truth.
`ALC-0315` is deliberately written AS DRAWN (3-carbon hydroxyalkyl) so the
mismatch against the corpus is recorded rather than silently corrected -- and
that mismatch is now independently confirmed twice, with the drawn head measured
at three vertices and full-length bonds, against four correct heads elsewhere in
the same panel.
"""
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

C12 = 'C'*12          # dodecyl
C14 = 'C'*14          # tetradecyl
OHT = 'CC(O)' + 'C'*12   # 2-hydroxytetradecyl, written N->CH2->CH(OH)->C12H25
HD  = 'CC(CCCCCC)' + 'C'*8  # 2-hexyldecyl, written N->CH2->CH(hexyl)->C8H17

def amine(a, b, rest):           # N(a)(b)rest
    return f'N({a})({b}){rest}'

DRAWN = {
 'ALC-0315': 'CCCCCCCCC(CCCCCC)C(=O)OCCCCCCN(CCCO)CCCCCCOC(=O)C(CCCCCC)CCCCCCCC',
 'SM-102'  : 'CCCCCCCCCCCOC(=O)CCCCCN(CCO)CCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCC',
 'MC3'     : r'CCCCC/C=C\C/C=C\CCCCCCCCC(CCCCCCCC/C=C\C/C=C\CCCCC)OC(=O)CCCN(C)C',
 'E1': f'OCCCCN({C12[:6]}NC(=O)N{HD})' + f'{C12[:6]}NC(=O)N{HD}',
 'E2': amine(OHT, OHT, 'CCNC(=O)CCN1CCN(CCC(=O)NCC' + f'N({OHT})({OHT})' + ')CC1'),
 'E3': amine(OHT, OHT, 'CCCN1CCN(CCC' + f'N({OHT})({OHT})' + ')CC1'),
 'E4': amine(C14, C14, 'CCCN1CCN(CCC' + f'N({C14})({C14})' + ')CC1'),
 'E5': amine(C14, C14, 'CCNC(=O)CCN(C)CCN(C)CCC(=O)NCC' + f'N({C14})({C14})'),
 'E6': amine(OHT, OHT, 'CCN(C)CCC(=O)N1CCN(C(=O)CCN(C)CC' + f'N({OHT})({OHT})' + ')CC1'),
 'E7': amine(C14, C14, 'CCN(C)CCC(=O)N1CCN(C(=O)CCN(C)CC' + f'N({C14})({C14})' + ')CC1'),
 'E8': amine(OHT, OHT, 'CCN(C)CCC(=O)NCCNC(=O)CCN(C)CC' + f'N({OHT})({OHT})'),
 'E9': f'N(CC(=O)N{C12})(CC(=O)N{C12})CCCN1CCN(CCC' + f'N(CC(=O)N{C12})(CC(=O)N{C12})' + ')CC1',
}

rows=[]
for k,s in DRAWN.items():
    m=Chem.MolFromSmiles(s)
    assert m is not None, (k,s)
    rows.append(dict(name=k, smiles_drawn=Chem.MolToSmiles(m),
                     formula=rdMolDescriptors.CalcMolFormula(m),
                     MolWt=round(Descriptors.MolWt(m),2),
                     ExactMolWt=round(Descriptors.ExactMolWt(m),4),
                     nN=sum(a.GetSymbol()=='N' for a in m.GetAtoms()),
                     nO=sum(a.GetSymbol()=='O' for a in m.GetAtoms()),
                     nC=sum(a.GetSymbol()=='C' for a in m.GetAtoms())))
df=pd.DataFrame(rows)
df.to_csv('/Users/nathanliu/LNP/confound/curation/xd_2026_lipidai/03_structures.csv',index=False)
pd.set_option('display.width',260); pd.set_option('display.max_colwidth',130)
print(df.to_string(index=False))
