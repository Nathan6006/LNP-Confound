"""Emit the LipidDB extraction as a standalone artifact (NOT merged into data/).

Decision 2026-10-29: LipidDB is a literature compilation - 485 rows whose
luminescence values were lifted from the 47 published studies the paper names
(§2.1.1), each on its own dose / route / timepoint / instrument scale. It is not
one experiment, so it cannot carry one Experiment_ID without fabricating a lab.
Three further facts make the pooling non-recoverable rather than merely untidy:
`descriptors_1107.csv` - the file the model actually trains on - carries an
`inject_way` feature with FOUR levels, so at least four administration routes are
pooled into the single "Total expression" column; 7 of the 24 duplicate
structures carry CONFLICTING class labels, with value disagreements up to five
orders of magnitude; and the released xlsx is not the complete released dataset,
since `inject_way` appears nowhere in it. It is kept here as a
second, independent literature-pooled corpus: useful as external evidence for
Results 1/2/5 and as a cross-compiler concordance test on the lipids it shares
with our corpus.

Also emitted: the MTA (Methyl Tail Augmentation) sheets, flagged `synthetic` -
those 720 rows are generated structures carrying an assumed label, not
measurements, and must never be treated as data.
"""
import re, json, collections
import pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
D='/Users/nathanliu/LNP/confound/curation/xd_2026_lipidai/'
XL=D+'repo/dataset/LipidDB-Ionizable lipid Library.xlsx'
canon=lambda s:(lambda m: Chem.MolToSmiles(m) if m else None)(Chem.MolFromSmiles(str(s)))
flat =lambda s:(lambda m:(Chem.RemoveStereochemistry(m), Chem.MolToSmiles(m))[1])(Chem.MolFromSmiles(s))

xl=pd.ExcelFile(XL)
s1=xl.parse('Sheet1')
s1.columns=['serial','smiles','ecfp1024','ecfp2048','total_expr','bin3','bin4','bin2',
            'refs','aux_a','aux_b','aux_c']
s1=s1[['serial','smiles','total_expr','bin3','bin4','bin2','refs']].copy()
for c in ('bin3','bin4','bin2'): s1[c]=s1[c].astype(str).str.strip()
s1['canonical_SMILES']=s1.smiles.map(canon)
# The authors appended an incrementing counter to each row so every reference
# string is unique, but NOT always to a trailing `doi:` token - the counter also
# rides on a PMCID, on a page range, or on a doi.org URL. Cutting only at `doi:`
# leaves 299 apparently-distinct strings, 6x the truth. Cut at the first
# identifier of ANY kind, which collapses them to ~50, against the 47 studies
# §2.1.1 states. This is the paper's own cross-check and it is the authoritative
# number; the string count is reported beside it, never instead of it.
IDENT=(r'\bdoi\s*:', r'https?://doi\.org/', r'\bPMID\s*:', r'\bPMCID\s*:')
def ref_key(x):
    x=re.sub(r'\s+',' ',str(x).replace('\xa0',' ')).strip()
    for pat in IDENT:
        m=re.search(pat,x,flags=re.I)
        if m: x=x[:m.start()]
    return x.strip(' .;')
s1['reference']=s1.refs.map(ref_key)
s1['reference_letters_only']=s1.refs.map(
    lambda x: re.sub(r'[^A-Za-z ]','',str(x))[:120].strip())
s1['source']='LipidDB_literature_compilation'
s1['value_kind']='in_vivo_total_luminescence_as_reported_by_each_source_publication'

aug=[]
for sh,tag in (('sheet1+tail','MTA_add_methyl'),('sheet1-tail','MTA_remove_methyl')):
    d=xl.parse(sh); d.columns=['serial','smiles','e1','e2','bin3','bin4','bin2','x']
    d=d[['serial','smiles','bin3','bin4','bin2']].copy()
    for c in ('bin3','bin4','bin2'): d[c]=d[c].astype(str).str.strip()
    d['canonical_SMILES']=d.smiles.map(canon); d['source']=tag
    d['value_kind']='SYNTHETIC - label inherited from the parent lipid, not measured'
    aug.append(d)
aug=pd.concat(aug,ignore_index=True)

# corpus overlap
vit=pd.read_csv(D+'../../data/data_vitro.csv',low_memory=False)
viv=pd.read_csv(D+'../../data/data_vivo.csv',low_memory=False)
cv={canon(x) for x in vit.IL_SMILES.dropna()}; cz={canon(x) for x in viv.IL_SMILES.dropna()}
s1['in_corpus_vitro']=s1.canonical_SMILES.isin(cv)
s1['in_corpus_vivo'] =s1.canonical_SMILES.isin(cz)

# duplicate structures whose class labels disagree between their source papers
conflicts=[]
for smi,g in s1.groupby('canonical_SMILES'):
    if len(g)>1 and g.bin3.nunique()>1:
        conflicts.append(dict(classes=sorted(set(g.bin3)),
                              values=sorted(v for v in g.total_expr.dropna()),
                              n_rows=int(len(g))))
n_conflict=len(conflicts)

# the feature file the model actually trains on, which the xlsx does not cover
desc=pd.read_csv(D+'repo/dataset/descriptors_1107.csv',low_memory=False,encoding='latin-1')
inject_levels=({int(k):int(v) for k,v in desc.inject_way.value_counts().items()}
               if 'inject_way' in desc.columns else None)

s1.to_csv(D+'06_lipiddb_original.csv',index=False)
aug.to_csv(D+'06_lipiddb_mta_synthetic.csv',index=False)

t1=pd.read_csv(D+'05_table1_crosscheck.csv')
s1f={flat(x) for x in s1.canonical_SMILES.dropna()}
t1['flat_in_xlsx']=t1.canonical.map(lambda s: flat(s) in s1f)

qc=dict(
 sheet1_rows=int(len(s1)),
 sheet1_unique_canonical=int(s1.canonical_SMILES.nunique()),
 sheet1_duplicate_rows=int(len(s1)-s1.canonical_SMILES.nunique()),
 sheet1_unparseable=int(s1.canonical_SMILES.isna().sum()),
 rows_with_numeric_value=int(s1.total_expr.notna().sum()),
 rows_without_value=int(s1.total_expr.isna().sum()),
 paper_states_n_lipids=475, xlsx_rows=485,
 distinct_source_publications_stated_by_paper=47,   # §2.1.1, verbatim
 distinct_reference_strings_after_normalisation=int(s1.reference.nunique()),
 distinct_reference_strings_letters_only=int(s1.reference_letters_only.nunique()),
 raw_reference_strings=int(s1.refs.nunique()),
 largest_single_source_rows=int(s1.reference.value_counts().iloc[0]),
 top10_sources_rows=int(s1.reference.value_counts().head(10).sum()),
 mta_rows=int(len(aug)), mta_unique=int(aug.canonical_SMILES.nunique()),
 mta_rows_identical_to_an_original_lipid=int(aug.canonical_SMILES.isin(
     set(s1.canonical_SMILES.dropna())).sum()),
 duplicate_structures_with_conflicting_class=int(n_conflict),
 duplicate_class_conflicts=conflicts,
 model_feature_file_has_inject_way_levels=inject_levels,
 overlap_with_corpus_vitro=int(s1.in_corpus_vitro.sum()),
 overlap_with_corpus_vivo=int(s1.in_corpus_vivo.sum()),
 overlap_with_corpus_union=int((s1.in_corpus_vitro|s1.in_corpus_vivo).sum()),
 table1_rows=int(len(t1)), table1_parsed=int(t1.parsed.sum()),
 table1_structure_in_xlsx=int(t1.in_xlsx.sum()),
 table1_structure_in_xlsx_ignoring_stereo=int(t1.flat_in_xlsx.sum()),
 table1_value_exact=int(t1.value_match.sum()),
 table1_class_agreement=int((t1.class_match==True).sum()),
 # the CLASS is the model's actual training label, so a class disagreement is a
 # worse defect than a value disagreement; both are listed per row
 table1_disagreements=[dict(table_row=int(r.table_row),
                            table_value=int(r.lum), table_class=r.table_class,
                            table_ref=None if pd.isna(r.table_ref) else int(r.table_ref),
                            xlsx_values=r.xlsx_values, xlsx_classes=r.xlsx_classes,
                            class_also_disagrees=bool(r.class_match is False))
                       for _,r in t1[t1.in_xlsx & (~t1.value_match | (t1.class_match==False))].iterrows()],
 table1_not_found=[int(r.table_row) for _,r in t1[~t1.flat_in_xlsx].iterrows()],
 merged_into_corpus=False,
 reason=("literature compilation across the 47 studies the paper names; a single "
         "Experiment_ID would fabricate a source and corrupt the source key "
         "Results 1-4 rest on. Independently fatal to pooling: the model's own "
         "feature file encodes FOUR administration routes (inject_way) that the "
         "released value column does not distinguish, and 7 duplicate structures "
         "carry conflicting class labels across their source papers."),
)
json.dump(qc,open(D+'06_lipiddb_qc.json','w'),indent=1)
print(json.dumps(qc,indent=1))
