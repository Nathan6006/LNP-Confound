"""All NG_2025 checks in one place -> 09_qc.json, and the flag file.

C1 is the one that would catch a miscounted chain, and this paper supplies it
for its lead compound only: the mass spectrum of C-a16 (Supplementary Fig 19a)
prints its own "Molecular Weight 1129.71" and an [M+H]+ at m/z 1130.4, and its
1H NMR integrates to a hydrogen count that is checked here against the
enumerated formula.
"""
import csv, json, math, os, collections
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors
RDLogger.DisableLog('rdApp.*')

D = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(D, '..', '..')

# the paper's own printed facts, transcribed once
PAPER = dict(
    library_size=310, in_vivo_hits=11,
    hits_per_phenol_tail_as_printed={'A': 0, 'B': 3, 'C': 6, 'D': 2},
    c_a16_mw_printed=1129.71, c_a16_mz_observed=1130.4,
    # 1H NMR (400 MHz, CDCl3): the integrals as printed
    c_a16_nmr_integrals=[8, 4, 2, 8, 12, 8, 12, 40, 12],
    c_a16_nmr_exchangeable_not_integrated=['4 phenol OH', '2 amide NH'],
)


def main():
    enum = list(csv.DictReader(open(os.path.join(D, '03_enumerated.csv'))))
    vals = json.load(open(os.path.join(D, '04_values.json')))
    rows = list(csv.DictReader(open(os.path.join(D, '07_NG_2025_rows.csv'))))
    qc = {}

    # C5 parse + nitrogen present
    mols = {r['lipid_name']: Chem.MolFromSmiles(r['IL_SMILES']) for r in enum}
    bad = [k for k, m in mols.items() if m is None]
    noN = [k for k, m in mols.items() if m is not None and
           not any(a.GetSymbol() == 'N' for a in m.GetAtoms())]
    rad = [k for k, m in mols.items() if m is not None and
           any(a.GetNumRadicalElectrons() for a in m.GetAtoms())]
    qc['C5_parse'] = dict(n=len(enum), unparseable=bad, no_nitrogen=noN, radicals=rad,
                          result='pass' if not (bad or noN or rad) else 'FAIL')

    # C3 library size, C4 duplicates
    canon = {k: Chem.MolToSmiles(m) for k, m in mols.items()}
    dup = [s for s, n in collections.Counter(canon.values()).items() if n > 1]
    qc['C3_library_size'] = dict(enumerated=len(enum), stated=PAPER['library_size'],
                                 result='pass' if len(enum) == PAPER['library_size'] else 'FAIL')
    qc['C4_duplicates_within_source'] = dict(n_duplicate_structures=len(dup),
                                             result='pass' if not dup else 'FAIL')
    qc['C4_noncanonical'] = dict(
        n=sum(1 for k, r in zip(canon, enum) if canon[k] != r['IL_SMILES']),
        result='pass' if all(canon[r['lipid_name']] == r['IL_SMILES'] for r in enum) else 'FAIL')

    # C1 mass of the one compound the paper characterises
    m = mols['C-a16']
    formula = rdMolDescriptors.CalcMolFormula(m)
    mw = Descriptors.MolWt(m)
    exact = Descriptors.ExactMolWt(m)
    nH = sum(a.GetTotalNumHs() for a in m.GetAtoms())
    nmr = sum(PAPER['c_a16_nmr_integrals'])
    qc['C1_mass_C_a16'] = dict(
        formula=formula, average_MW=round(mw, 2), printed_MW=PAPER['c_a16_mw_printed'],
        delta_MW=round(mw - PAPER['c_a16_mw_printed'], 3),
        exact_mass=round(exact, 4), M_plus_H_exact=round(exact + 1.00728, 4),
        M_plus_H_average=round(mw + 1.008, 2),
        observed_mz=PAPER['c_a16_mz_observed'],
        doubly_charged_expected=round((exact + 2 * 1.00728) / 2, 1),
        doubly_charged_observed=565.8,
        total_H_in_structure=nH, nmr_integral_total=nmr,
        exchangeable_H=nH - nmr,
        result='pass' if abs(mw - PAPER['c_a16_mw_printed']) < 0.05 and nH - nmr == 6 else 'FAIL',
        note='nH - NMR total must equal 6: the 4 phenol OH and 2 amide NH are not '
             'integrated in the printed list, which is itself a structural statement '
             '(four phenols, and the diketopiperazine NH unsubstituted)')

    # C6 molar ratios
    ratios = {(r['IL_molratio'], r['HL_molratio'], r['CHL_molratio'], r['PEG_molratio'])
              for r in rows}
    sums = {round(sum(float(x) for x in t), 3) for t in ratios}
    qc['C6_molar_ratios'] = dict(distinct=sorted(ratios), sums=sorted(sums),
                                 result='pass' if sums == {100.0} else 'FAIL')

    # C9 values: the paper's own hit rule, reconstructed
    ph = {v['lipid_name']: v['log10_luminescence'] for v in vals['vitro'] if v['block'] == 'Fig.1e'}
    panel = ['B-a11', 'B-a16', 'C-a9', 'C-a10', 'C-a11', 'C-a12', 'C-a15', 'C-a16',
             'D-a9', 'D-a11', 'D-a16']
    best = None
    for t in np.arange(1.5, 7.0, 0.001):
        hit = sorted(k for k, v in ph.items() if v > t)
        agree = len(set(hit) & set(panel)) - abs(len(hit) - len(panel))
        if best is None or agree > best[0]:
            best = (agree, float(t), hit)
    _, thr, hit = best
    per_tail = {L: sum(1 for k in hit if k.startswith(L + '-')) for L in 'ABCD'}
    qc['C2_hit_rule_reconstruction'] = dict(
        rule="the paper defines a hit as luminescence > 5x DLin-MC3-DMA",
        scope=('PHENOL LIBRARY ONLY (Fig. 1e). The paper states the rule inside its '
               'phenol-library paragraph and reports hits per phenol tail, so the '
               'reconstruction is scoped the same way. Applied to all 310 lipids the '
               'same threshold also catches G-a12 (thiol, log10 = 6.074), which clears '
               'it yet was not taken into mice - so no single threshold over the whole '
               'library reproduces exactly the 11.'),
        thiol_lipid_also_above_threshold=dict(name='G-a12', log10=6.074034),
        phenol_threshold_window='(5.937736, 6.036851]',
        threshold_log10=round(thr, 3),
        implied_MC3_luminescence=round(10 ** thr / 5, 1),
        hits=hit, in_vivo_panel=panel,
        reproduced=hit == sorted(panel),
        hits_per_tail_from_data=per_tail,
        hits_per_tail_as_printed=PAPER['hits_per_phenol_tail_as_printed'],
        result='pass' if hit == sorted(panel) else 'FAIL',
        discrepancy=('the Results text prints 0, 3, 6 and 2 hits for phenols A-D; the '
                     'Source Data give 0, 2, 6 and 3, and the paper\'s OWN in-vivo panel '
                     '(Fig. 1k-p) carries 2 B-lipids and 3 D-lipids, so the printed B and '
                     'D counts are transposed. Which 11 lipids are hits is not in doubt: '
                     'one threshold reproduces all 11 exactly.'))

    # C9 distribution
    v = [float(r['Experiment_value']) for r in rows]
    blocks = collections.defaultdict(list)
    for r in rows:
        blocks[(r['Experiment_ID'], r['Model_target'])].append(float(r['Experiment_value']))
    qc['C9_value_distribution'] = dict(
        rows=len(rows),
        blocks={f'{k[0]}|{k[1]}': dict(n=len(x), mean=round(float(np.mean(x)), 12),
                                       sd=round(float(np.std(x)), 9),
                                       min=round(min(x), 3), max=round(max(x), 3))
                for k, x in blocks.items()},
        result='pass' if all(abs(np.mean(x)) < 1e-6 and abs(np.std(x) - 1) < 1e-6
                             for x in blocks.values()) else 'FAIL')

    # C8 overlap with the corpus as it stands
    import pandas as pd
    have = set()
    for f in ('data_vitro.csv', 'data_vivo.csv'):
        d = pd.read_csv(os.path.join(ROOT, 'data', f), low_memory=False)
        have |= set(d.IL_SMILES.dropna())
    ov = {k for k, s in canon.items() if s in have}
    bench = {r['IL_name'] for r in rows if r['IL_SMILES'] in have} - set(canon)
    qc['C8_overlap_with_corpus'] = dict(
        novel_lipids=len(canon) - len(ov), already_in_corpus=sorted(ov),
        benchmark_rows_matching_corpus=sorted(bench),
        note='zero overlap among the 310 novel lipids is the expected result for a new '
             'chemistry; the benchmark rows are expected to match')

    # heat-map cross-check, carried through
    qc['C7_heatmap_vs_source_data'] = json.load(open(os.path.join(D, '06_heatmap_check.json')))
    qc['C10_drawn_structures'] = dict(
        anchors=[r['lipid_name'] for r in csv.DictReader(open(os.path.join(D, 'reference_structures.csv')))],
        note='run provenance/structure_check.py; all six were exact matches on 2026-10-29')

    qc['held_back'] = list(csv.DictReader(open(os.path.join(D, '07_held_back.csv'))))
    json.dump(qc, open(os.path.join(D, '09_qc.json'), 'w'), indent=1)

    for k, val in qc.items():
        if isinstance(val, dict) and 'result' in val:
            print(f"  {val['result']:5} {k}")
    print()
    print('C1  C-a16', qc['C1_mass_C_a16']['formula'],
          'MW', qc['C1_mass_C_a16']['average_MW'],
          'vs printed', qc['C1_mass_C_a16']['printed_MW'],
          '| [M+H]+ calc', qc['C1_mass_C_a16']['M_plus_H_average'],
          'obs', qc['C1_mass_C_a16']['observed_mz'],
          '| [M+2H]2+ calc', qc['C1_mass_C_a16']['doubly_charged_expected'], 'obs 565.8')
    print('C2  hit threshold log10 =', qc['C2_hit_rule_reconstruction']['threshold_log10'],
          '-> reproduces', len(hit), 'of 11 in-vivo lipids exactly:',
          qc['C2_hit_rule_reconstruction']['reproduced'])
    print('C8  novel lipids', qc['C8_overlap_with_corpus']['novel_lipids'], '/', len(canon),
          '| benchmarks matching corpus:', qc['C8_overlap_with_corpus']['benchmark_rows_matching_corpus'])


if __name__ == '__main__':
    main()
