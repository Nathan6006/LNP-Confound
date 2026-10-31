"""Emit XD_2026: the 12-lipid in-vivo screen the LipidAI paper ran itself.

Source: Duan X et al., Acta Pharm Sin B 2026;16(16):5053-5069, PMID 42592089,
doi 10.1016/j.apsb.2025.10.040. Fig 5A (structures), Fig 5C (luciferase
expression), Fig 5D (the A/B/C class of each lipid).

What is and is not in here
--------------------------
IN  : 12 LNPs - E1-E9 plus ALC-0315, SM-102 and MC3 - one formulation, BALB/c
      male mice, tail-vein IV, whole-body IVIS at 6 h, p/sec/cm^2/sr.
OUT : the 485-row LipidDB and its 720 MTA rows. LipidDB is a compilation across
      299 publications, not an experiment (see stage06_lipiddb.py); the MTA rows
      are generated structures with inherited labels, not measurements.

Values
------
Fig 5C is a log-axis bar chart, digitised in stage04_fig5c.py against the
10^6-10^9 major ticks (tick residuals < 0.002 decades). Every digitised value
falls in the A/B/C band Fig 5D independently assigns to that lipid, 12/12 - two
figures agreeing is the check that the digitisation is right.
Experiment_value is then z(log10(value)) over this experiment, population sd,
which is the corpus convention.

Structures
----------
E1-E9 were read off Fig 5A geometrically (stage02_trace.py: skeletonise, trace,
RDP-simplify, count vertices), never by eye. The same procedure reproduces
SM-102 and MC3 EXACTLY against the corpus - but ALC-0315 is DRAWN WITH A
3-CARBON HYDROXYALKYL HEAD where the real compound has 4 (hydroxybutyl). That
is a drawing error in the paper, confirmed by the other two controls being
exact. Consequence, recorded honestly: the three benchmarks are emitted with
their CORPUS SMILES, and E1-E9 carry confidence=medium, because they come from
a figure with a demonstrated +-1 carbon error rate of 1 in 3 on its controls.

Formulation (Methods 2.3.3): IL:DSPC:Chol:DMG-PEG2000 = 50:10:38.5:1.5 mol%,
microfluidic mixing, aqueous phase 0.01 M citrate pH 3, dialysed into 0.01 M
citrate pH 6. mRNA dose and N:P are NOT stated in the paper, so they stay blank.
"""
import csv, json, math, os
import pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

D   = os.path.dirname(os.path.abspath(__file__))
ROOT= os.path.join(D, '..', '..')
EID = 'XD_2026'
FID = 'XD_2026_Fvivo'
LINK= 'https://www.sciencedirect.com/science/article/pii/S2211383525007221'
DOI = '10.1016/j.apsb.2025.10.040'
PMID= '42592089'

corpus = pd.read_csv(os.path.join(ROOT,'data','data_vivo.csv'), low_memory=False)
def comp(col_name, col_smiles, name):
    s = corpus.loc[corpus[col_name]==name, col_smiles].dropna()
    assert len(s), name
    return s.iloc[0]
DSPC = comp('HL_name','HL_SMILES','DSPC')
CHOL = comp('CHL_name','CHL_SMILES','Cholesterol')
PEG  = comp('PEG_name','PEG_SMILES','DMG-PEG2000')
def il(name):
    s = corpus.loc[corpus.IL_name==name,'IL_SMILES'].dropna()
    if not len(s):
        s = pd.read_csv(os.path.join(ROOT,'data','data_vitro.csv'),
                        low_memory=False).pipe(lambda d: d.loc[d.IL_name==name,'IL_SMILES']).dropna()
    assert len(s), name
    return Chem.MolToSmiles(Chem.MolFromSmiles(s.iloc[0]))

drawn = pd.read_csv(os.path.join(D,'03_structures.csv')).set_index('name')
fig5c = {b['name']: b['value'] for b in json.load(open(os.path.join(D,'04_fig5c.json')))['bars']}
fig5c['DLin-MC3-DMA'] = fig5c.pop('MC3')
# Fig 5D "Reality" column, read off the rendered panel; A >5e7, B 5e6-5e7, C <5e6
CLASS = dict(zip(['ALC-0315','SM-102','DLin-MC3-DMA','E1','E2','E3','E4','E5',
                  'E6','E7','E8','E9'], list('AAAAABBBCCCC')))
# Fig 5A labels the benchmarks ALC-0315 / SM-102 / MC3; the corpus name for the
# third is DLin-MC3-DMA and that is what is EMITTED, so the corpus keeps one
# alias per molecule and the C10 anchor can be matched by name.
BENCH = {'ALC-0315':'ALC-0315', 'SM-102':'SM-102', 'DLin-MC3-DMA':'DLin-MC3-DMA'}
FIG_LABEL = {'DLin-MC3-DMA':'MC3'}

def band(v):
    return 'A' if v > 5e7 else ('B' if v > 5e6 else 'C')

names = list(CLASS)
vals  = [fig5c[n] for n in names]
logs  = [math.log10(v) for v in vals]
mu    = sum(logs)/len(logs)
sd    = (sum((x-mu)**2 for x in logs)/len(logs))**0.5         # population sd, ddof=0
z     = [(x-mu)/sd for x in logs]

COLS = open(os.path.join(ROOT,'data','data_vivo.csv')).readline().strip().split(',')
for c in (x['name'] for x in json.load(open(os.path.join(
        ROOT,'provenance','schema_quality.json')))['columns']):
    if c not in COLS:
        COLS.append(c)
rows, meta, mism = [], [], []
for n, v, zz in zip(names, vals, z):
    is_bench = n in BENCH
    key = FIG_LABEL.get(n, n)                      # the label Fig 5A uses
    smi = il(BENCH[n]) if is_bench else drawn.loc[key,'smiles_drawn']
    if is_bench and il(BENCH[n]) != drawn.loc[key,'smiles_drawn']:
        mism.append(dict(lipid=n, drawn=drawn.loc[key,'smiles_drawn'], corpus=il(BENCH[n])))
    assert band(v) == CLASS[n], (n, v, band(v), CLASS[n])
    r = {c:'' for c in COLS}
    r.update(Index='', LNP_ID='', Experiment_ID=EID, Formulation_ID=FID,
             IL_name=n, IL_SMILES=smi, IL_molratio='50',
             HL_name='DSPC', HL_SMILES=DSPC, HL_molratio='10',
             CHL_name='Cholesterol', CHL_SMILES=CHOL, CHL_molratio='38.5',
             PEG_name='DMG-PEG2000', PEG_SMILES=PEG, PEG_molratio='1.5',
             Aqueous_buffer='citrate', Dialysis_buffer='citrate',
             Mixing_method='microfluidics',
             Model='in_vivo', Model_type='Mouse_BALBc', Model_target='whole_body',
             Route_of_administration='intravenous', Cargo='mRNA', Cargo_type='FLuc',
             Experiment_method='luminescence_normalized',
             Experiment_batching='individual', Experiment_value=repr(zz),
             Publication_link=LINK, Publication_DOI=DOI, Publication_PMID=PMID)
    # Confidence onto the ROW, 2026-10-30 (review finding F-001). It was always
    # computed -- see `meta` below, and 07_qc.json's E_lipids_confidence -- but
    # only ever written to this folder's own files, so data_vivo.csv carried
    # nine vendor-supplied structures with no external anchor of any kind and
    # said nothing about it. E1-E9 have no name, no formula, no HRMS and no SI
    # structure: the ONLY source is a geometric trace of the Fig 5A drawing.
    r.update(confidence='high' if is_bench else 'medium')
    if not is_bench:
        r.update(structure_assumed='1', curation_note=(
            'Structure recovered only by geometric trace of the Fig 5A drawing '
            '(stage02_trace.py). The paper gives this lipid no name, no '
            'molecular formula, no HRMS and no SI structure, so there is no '
            'external anchor against which to check it.'))
    rows.append(r)
    meta.append(dict(IL_name=n, source='benchmark_corpus_SMILES' if is_bench
                                        else 'Fig5A_geometric_trace',
                     confidence='high' if is_bench else 'medium',
                     fig5c_value_p_per_s_per_cm2_per_sr=round(v),
                     fig5c_log10=round(math.log10(v),4),
                     fig5d_reality_class=CLASS[n], class_from_value=band(v),
                     Experiment_value=zz))

assert abs(sum(z)/len(z))<1e-12 and abs((sum(x*x for x in z)/len(z))**0.5-1)<1e-12
with open(os.path.join(D,'07_XD_2026_rows.csv'),'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=COLS); w.writeheader(); w.writerows(rows)
pd.DataFrame(meta).to_csv(os.path.join(D,'07_XD_2026_values.csv'),index=False)

# C10 inputs: what we emit, and the externally-known benchmark structures
with open(os.path.join(D,'03_enumerated.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['lipid_name','IL_SMILES'])
    for r in rows: w.writerow([r['IL_name'], r['IL_SMILES']])

qc=dict(experiment_id=EID, model='in_vivo', rows=len(rows),
        unique_canonical_SMILES=len({r['IL_SMILES'] for r in rows}),
        molratio_sum=50+10+38.5+1.5,
        zscore_mean=sum(z)/len(z), zscore_sd=(sum(x*x for x in z)/len(z))**0.5,
        value_range=[min(vals),max(vals)],
        fig5c_vs_fig5d_class_agreement=f'{len(names)}/{len(names)}',
        # C2 in audit.py is head-series mass monotonicity; comparing an emitted
        # structure against an external reference is a C10-class check.
        benchmark_C10_external=dict(checked=3, exact=3-len(mism), mismatches=mism),
        benchmarks_emitted_with='corpus SMILES',
        E_lipids_confidence='medium',
        # fields that are READINGS, not facts printed in the paper
        assumptions=dict(
          Model_target=("'whole_body'. Methods 2.3.4 describes BOTH whole-animal "
                        "imaging at 6 h AND ex-vivo heart/liver/spleen/lung/kidney "
                        "on the same instrument; neither the Fig 5C caption nor the "
                        "text says which ROI the histogram plots. whole_body is the "
                        "defensible reading for one number per lipid, but it is a "
                        "reading."),
          Route_of_administration=("'intravenous'. Methods 2.3.4 says only "
                        "'post-vaccination'; IV / tail-vein comes from the Fig 5 "
                        "caption and the body text, which agree with each other."),
          Dialysis_buffer=("'citrate' - faithful to Methods 2.3.3 (0.01 M citric "
                        "acid pH 6, 2 h) but a category the corpus has never "
                        "carried; every other source dialyses into PBS or saline."),
          Dose_ug_nucleicacid="absent from the paper; left blank, not imputed.",
          N_to_P_ratio="absent from the paper; the corpus has no column for it."),
        lipiddb_merged=False)
json.dump(qc, open(os.path.join(D,'07_qc.json'),'w'), indent=1)
print(json.dumps(qc, indent=1))
print()
print(pd.DataFrame(meta).to_string(index=False))
