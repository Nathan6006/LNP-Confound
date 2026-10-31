"""Cross-check the xlsx LipidDB against Table 1, the 26 rows the PAPER prints.

Table 1 is independent of the repo file, so it is the only external ground
truth available for LipidDB.

Parsing note: SMILES wrap over several lines, and the luminescence value is set
at the VERTICAL CENTRE of its row, so neither text order nor "nearest value"
assigns the wrap lines correctly. Rows are therefore cut 7 pt above each value
anchor, which brackets every block exactly; a DP split on "must parse as a
molecule" was tried first and rejected, because SMILES prefixes parse too and
it silently chose a shifted split. Every assembled SMILES is checked to parse.
'=' renders as U+FFFD in the extracted text.
"""
import re, collections
import pandas as pd, pdfplumber
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
D='/Users/nathanliu/LNP/confound/curation/xd_2026_lipidai/'
PDF=D+'1-s2.0-S2211383525007221-main.pdf'
SMI=re.compile('^[CNOSPHBrIFcnosp0-9()\\[\\]@+\\-=#/\\\\.�]+$')
# page index -> (smiles col x0, value col x0, table top, table bottom)
COL={3:(195,425,100,999), 4:(185,416,95,317)}
PAD=7.0                      # pt above a value anchor where its row starts

rows=[]
with pdfplumber.open(PDF) as pdf:
    for pno,(cx0,cxv,t0,t1) in COL.items():
        w=pdf.pages[pno].extract_words()
        frag=[x for x in w if cx0<=x['x0']<cxv and t0<=x['top']<t1
              and SMI.match(x['text']) and len(x['text'])>=2
              and not x['text'].isdigit()]
        val=[x for x in w if cxv<=x['x0']<cxv+80 and t0<=x['top']<t1
             and re.fullmatch(r'\d{4,}',x['text'])]
        val.sort(key=lambda x:x['top'])
        edges=[v['top']-PAD for v in val]+[10**6]
        for i,v in enumerate(val):
            f=[x for x in frag if edges[i]<=x['top']<edges[i+1]]
            f.sort(key=lambda x:(round(x['top'],1),x['x0']))
            cls=[x['text'] for x in w if edges[i]<=x['top']<edges[i+1]
                 and cxv+45<=x['x0']<cxv+75 and x['text'] in 'ABCD']
            ref=[x['text'] for x in w if edges[i]<=x['top']<edges[i+1]
                 and cxv+75<=x['x0']<cxv+120 and x['text'].isdigit()]
            rows.append(dict(page=pno+1, lum=int(v['text']),
                             table_class=cls[0] if cls else '',
                             table_ref=int(ref[0]) if ref else None,
                             n_lines=len(f),
                             smiles=''.join(x['text'] for x in f).replace('�','=')))

for i,r in enumerate(rows,1):
    r['table_row']=i
    m=Chem.MolFromSmiles(r['smiles'])
    r['parsed']=m is not None
    r['canonical']=Chem.MolToSmiles(m) if m else None
t1=pd.DataFrame(rows)
print(f"Table 1 rows: {len(t1)}   SMILES parsed: {t1.parsed.sum()}")
if (~t1.parsed).any():
    print(t1[~t1.parsed][['table_row','lum','smiles']].to_string(index=False))

db=pd.read_csv(D+'06_lipiddb_original.csv')
g=collections.defaultdict(list)
gc=collections.defaultdict(list)
gr=collections.defaultdict(list)
for sm,v,b,rf in zip(db.canonical_SMILES, db.total_expr, db.bin3, db.reference):
    g[sm].append(v); gc[sm].append(b); gr[sm].append(rf)
t1['in_xlsx']=t1.canonical.map(lambda s: bool(g.get(s)))
t1['xlsx_values']=t1.canonical.map(lambda s: g.get(s,[]))
t1['value_match']=[any(pd.notna(v) and abs(v-l)/max(l,1)<0.01 for v in vs)
                   for l,vs in zip(t1.lum,t1.xlsx_values)]
# the CLASS is the label the model actually trains on, so a class disagreement is
# a worse defect than a value disagreement and is reported separately
t1['xlsx_classes']=t1.canonical.map(lambda s: gc.get(s,[]))
t1['class_match']=[(c in cs) if cs else None for c,cs in zip(t1.table_class,t1.xlsx_classes)]
t1['xlsx_reference']=t1.canonical.map(lambda s: (gr.get(s) or [''])[0][:70])
print(f"\nstructures found in xlsx : {t1.in_xlsx.sum()}/{len(t1)}")
print(f"luminescence agreement   : {t1.value_match.sum()}/{t1.in_xlsx.sum()}")
print(f"class agreement          : {int((t1.class_match==True).sum())}/{t1.in_xlsx.sum()}")
print(t1[['table_row','lum','table_class','table_ref','parsed','in_xlsx',
          'value_match','class_match','xlsx_values','xlsx_classes']].to_string(index=False))
t1.to_csv(D+'05_table1_crosscheck.csv',index=False)
