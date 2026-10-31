import pandas as pd, os
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
D='/Users/nathanliu/LNP/confound/curation/xd_2026_lipidai'
P=os.path.join(D,'repo/dataset/LipidDB-Ionizable lipid Library.xlsx')
xl=pd.ExcelFile(P)
def canon(s):
    m=Chem.MolFromSmiles(str(s))
    return Chem.MolToSmiles(m) if m else None
store={}
for s in xl.sheet_names:
    d=xl.parse(s)
    cs=d.iloc[:,1].map(canon)
    store[s]=(d,cs)
    print(s,'n=',len(d),'unparseable=',cs.isna().sum(),'unique=',cs.nunique(),'dups=',cs.duplicated().sum())
# cross-sheet relation
a=set(store['Sheet1'][1].dropna()); b=set(store['sheet1+tail'][1].dropna()); c=set(store['sheet1-tail'][1].dropna())
print('Sheet1 n uniq',len(a),'+tail',len(b),'-tail',len(c))
print('+tail inside Sheet1:',len(b&a),'  -tail inside Sheet1:',len(c&a))
# corpus overlap
vit=pd.read_csv('/Users/nathanliu/LNP/confound/data/data_vitro.csv',low_memory=False)
viv=pd.read_csv('/Users/nathanliu/LNP/confound/data/data_vivo.csv',low_memory=False)
print([x for x in vit.columns if 'SMILES' in x.upper()])
sc=[x for x in vit.columns if x.upper()=='IL_SMILES'] or [x for x in vit.columns if 'SMILES' in x.upper()]
col=sc[0]; print('using',col)
cv=set(vit[col].dropna().map(canon).dropna()); cz=set(viv[col].dropna().map(canon).dropna())
print('corpus vitro uniq lipids',len(cv),'vivo',len(cz))
print('Sheet1 overlap vitro:',len(a&cv),'/',len(a))
print('Sheet1 overlap vivo :',len(a&cz),'/',len(a))
print('Sheet1 overlap union:',len(a&(cv|cz)),'/',len(a))
store['Sheet1'][0].assign(canonical_SMILES=store['Sheet1'][1]).to_csv(os.path.join(D,'01_sheet1_canon.csv'),index=False)
