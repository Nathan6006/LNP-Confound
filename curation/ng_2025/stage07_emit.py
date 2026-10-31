"""Emit NG_2025 rows in the 72-column corpus schema. NOT merged.

Three blocks, each z-scored on its own (population sd, ddof=0) exactly as the
corpus convention requires, and blocked the way BL_2023 and ZH_2023 already
are - by (Experiment_ID, endpoint), not by Experiment_ID alone:

  NG_2025     310 in-vitro lipids, MC38, Fig. 1e-j. The workbook's numbers are
              ALREADY log10(luminescence), so they are z-scored as they stand;
              taking log10 again would be wrong.
  NG_2025_2   the in-vivo follow-up of the 11 hits, two endpoints:
              liver after IV (Fig. 1m) and muscle after IM (Fig. 1p).
              Radiance is linear, so it is log10'd first, then z-scored.

Benchmark rows use the corpus's own canonical SMILES for DLin-MC3-DMA and the
corpus excipient strings, so C11 and the shared-lipid analysis see the same
molecule the rest of the corpus does. The benchmark LNP is the one formulation
here that uses DSPC rather than DOPE (Methods), and that is carried per row.

Held back, with reasons, in 07_held_back.csv.
"""
import csv, json, math, os
import numpy as np

D = os.path.dirname(os.path.abspath(__file__))
COLS = [c.strip() for c in open(os.path.join(D, '..', '..', 'data', 'data_vitro.csv')).readline().split(',')]

DOPE = 'CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[NH3+])OC(=O)CCCCCCC/C=C\\CCCCCCCC'
DSPC = 'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCCCCCC'
CHOL = 'CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C'
MC3 = 'CCCCC/C=C\\C/C=C\\CCCCCCCCC(CCCCCCCC/C=C\\C/C=C\\CCCCC)OC(=O)CCCN(C)C'
DOI = '10.1038/s41551-025-01422-8'
PMID = 40681859
LINK = 'https://www.nature.com/articles/s41551-025-01422-8'


def peg_smiles():
    import pandas as pd
    d = pd.read_csv(os.path.join(D, '..', '..', 'data', 'data_vivo.csv'), low_memory=False)
    s = d.loc[d.PEG_name == 'DMG-PEG2000', 'PEG_SMILES'].dropna()
    return s.value_counts().index[0]


def blank_row():
    return {c: '' for c in COLS}


def z(values):
    a = np.asarray(values, float)
    return (a - a.mean()) / a.std(ddof=0)


def main():
    PEG = peg_smiles()
    enum = {r['lipid_name']: r for r in csv.DictReader(open(os.path.join(D, '03_enumerated.csv')))}
    V = json.load(open(os.path.join(D, '04_values.json')))
    rows, held = [], []

    # ---------------- in vitro: 310 lipids, one block ----------------
    vit = V['vitro']
    missing = [v['lipid_name'] for v in vit if v['lipid_name'] not in enum]
    assert not missing, missing[:5]
    zs = z([v['log10_luminescence'] for v in vit])
    for v, zv in zip(vit, zs):
        e = enum[v['lipid_name']]
        r = blank_row()
        r.update({
            'Experiment_ID': 'NG_2025', 'Formulation_ID': 'NG_2025_Fvitro',
            'IL_name': v['lipid_name'], 'IL_SMILES': e['IL_SMILES'],
            'IL_molratio': 50, 'HL_name': 'DOPE', 'HL_SMILES': DOPE, 'HL_molratio': 10,
            'CHL_name': 'Cholesterol', 'CHL_SMILES': CHOL, 'CHL_molratio': 38.5,
            'PEG_name': 'DMG-PEG2000', 'PEG_SMILES': PEG, 'PEG_molratio': 1.5,
            'Aqueous_buffer': 'citrate', 'Mixing_method': 'handmixed',
            'Model': 'in_vitro', 'Model_type': 'MC38', 'Model_target': 'in_vitro',
            'Route_of_administration': 'in_vitro', 'Cargo': 'mRNA', 'Cargo_type': 'FLuc',
            'Dose_ug_nucleicacid': 0.01,
            'Experiment_method': 'luminescence_normalized',
            'Experiment_batching': 'individual',
            'Experiment_value': round(float(zv), 6),
            'Publication_link': LINK, 'Publication_DOI': DOI, 'Publication_PMID': PMID,
        })
        rows.append(r)

    # ---------------- in vivo: two endpoint blocks ----------------
    VIVO = {
        'Fig.1m': dict(target='liver', route='intravenous', panel='Fig. 1m'),
        'Fig.1p': dict(target='muscle', route='intramuscular', panel='Fig. 1p'),
    }
    for key, cfg in VIVO.items():
        blk = V['vivo'][key]
        names, keep = blk['names'], []
        for nm in names:
            if nm.lower().startswith('dlin') and key == 'Fig.1m':
                # the benchmark is the one formulation here built on DSPC, so it is
                # a different formulation and gets its own Formulation_ID
                keep.append((nm, 'DLin-MC3-DMA', MC3, DSPC, 'DSPC', 'NG_2025_Fvivo_MC3'))
            elif nm in enum:
                keep.append((nm, nm, enum[nm]['IL_SMILES'], DOPE, 'DOPE', 'NG_2025_Fvivo'))
            else:
                held.append(dict(panel=cfg['panel'], label=nm,
                                 reason='benchmark identity contradicts itself: the Fig. 1p '
                                        'x-axis labels this bar DLin-MC3-DMA, the Source Data '
                                        'column for the same three numbers is headed SM-102. '
                                        'The two are different molecules, so emitting either '
                                        'would put measured values on a possibly wrong structure.'))
        vals = [blk['mean'][nm] for nm, *_ in keep]
        zs = z([math.log10(v) for v in vals])
        for (label, name, smi, hl, hlname, fid), zv, raw in zip(keep, zs, vals):
            r = blank_row()
            r.update({
                'Experiment_ID': 'NG_2025_2', 'Formulation_ID': fid,
                'IL_name': name, 'IL_SMILES': smi,
                'IL_molratio': 50, 'HL_name': hlname, 'HL_SMILES': hl, 'HL_molratio': 10,
                'CHL_name': 'Cholesterol', 'CHL_SMILES': CHOL, 'CHL_molratio': 38.5,
                'PEG_name': 'DMG-PEG2000', 'PEG_SMILES': PEG, 'PEG_molratio': 1.5,
                'Aqueous_buffer': 'citrate', 'Dialysis_buffer': 'PBS',
                'Mixing_method': 'microfluidics',
                'Model': 'in_vivo', 'Model_type': 'Mouse_B6',
                'Model_target': cfg['target'],
                'Route_of_administration': cfg['route'],
                'Cargo': 'mRNA', 'Cargo_type': 'FLuc',
                'Dose_ug_nucleicacid': 2.5,
                'Experiment_method': 'luminescence_normalized',
                'Experiment_batching': 'individual',
                'Experiment_value': round(float(zv), 6),
                'Publication_link': LINK, 'Publication_DOI': DOI, 'Publication_PMID': PMID,
            })
            rows.append(r)

    with open(os.path.join(D, '07_NG_2025_rows.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        w.writeheader(); w.writerows(rows)
    with open(os.path.join(D, '07_held_back.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['panel', 'label', 'reason'])
        w.writeheader(); w.writerows(held)

    import collections
    c = collections.Counter((r['Experiment_ID'], r['Model_target']) for r in rows)
    print('rows:', len(rows))
    for k, n in sorted(c.items()):
        sub = [float(r['Experiment_value']) for r in rows
               if (r['Experiment_ID'], r['Model_target']) == k]
        print(f'  {k[0]:10} {k[1]:10} n={n:3d}  mean={np.mean(sub):+.2e}  sd={np.std(sub):.6f}')
    print('held back:', len(held))
    for h in held:
        print('  ', h['panel'], h['label'])


if __name__ == '__main__':
    main()
