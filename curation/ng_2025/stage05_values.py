"""Read the 310 in-vitro values and the two in-vivo panels from Source Data.

The publisher's Source Data workbook (MOESM3, sheet 'Fig1') holds the numbers
BEHIND the Fig. 1 heat maps, so nothing here is digitised from pixels. Each
block is a labelled grid: rows are the amine codes the sheet itself prints,
columns are the tails in the order the figure's own axis prints them, which
stage05_heatmap.py confirms against the drawn cell colours.

Fig.1m = isolated-liver radiance 6 h after IV injection (n = 3 mice).
Fig.1p = whole-mouse radiance after intramuscular injection (n = 3 mice).
"""
import json, os, sys
import openpyxl

D = os.path.dirname(os.path.abspath(__file__))
XL = os.path.join(D, '41551_2025_1422_MOESM3_ESM.xlsx')

# block header -> (tail labels left-to-right, lipid-name template)
BLOCK = {
    'Fig.1e': (['A', 'B', 'C', 'D'], '{t}-{a}'),
    'Fig.1f': (['E', 'F', 'G', 'H'], '{t}-{a}'),
    'Fig.1g': (['I', 'J', 'K', 'L'], '{t}-a-{a}'),
    'Fig.1h': (['I', 'J', 'K', 'L'], '{t}-b-{a}'),
    'Fig.1i': (['M', 'N', 'O'], '{t}-{a}'),
    'Fig.1j': (['P', 'Q', 'R'], '{t}-{a}'),
}


def cells(ws):
    return [[c for c in row] for row in ws.iter_rows(values_only=True)]


def main():
    ws = openpyxl.load_workbook(XL, data_only=True)['Fig1']
    g = cells(ws)
    # locate each block header and read the grid under it
    vitro, meta = [], {}
    for r, row in enumerate(g):
        for c, v in enumerate(row):
            if v in BLOCK:
                tails, tmpl = BLOCK[v]
                n = 0
                for rr in range(r + 1, len(g)):
                    lab = g[rr][c]
                    if not isinstance(lab, str) or not lab[:1] in 'ab':
                        break
                    vals = g[rr][c + 1: c + 1 + len(tails)]
                    assert all(isinstance(x, (int, float)) for x in vals), (v, lab, vals)
                    for t, val in zip(tails, vals):
                        vitro.append(dict(block=v, lipid_name=tmpl.format(t=t, a=lab),
                                          tail=t, amine=lab, log10_luminescence=float(val)))
                    n += 1
                meta[v] = dict(rows=n, cols=len(tails), n=n * len(tails))

    vivo = {}
    for key in ('Fig.1m', 'Fig.1p'):
        for r, row in enumerate(g):
            for c, v in enumerate(row):
                if v == key:
                    names = [x for x in g[r + 1][c:] if isinstance(x, str)]
                    reps = []
                    for rr in range(r + 2, len(g)):
                        vals = g[rr][c: c + len(names)]
                        if not all(isinstance(x, (int, float)) for x in vals):
                            break
                        reps.append([float(x) for x in vals])
                    vivo[key] = dict(names=names,
                                     replicates=reps,
                                     mean={nm: sum(rep[i] for rep in reps) / len(reps)
                                           for i, nm in enumerate(names)})

    out = dict(blocks=meta, n_vitro=len(vitro), vitro=vitro, vivo=vivo)
    json.dump(out, open(os.path.join(D, '04_values.json'), 'w'), indent=1)
    print('in-vitro values:', len(vitro))
    for k, m in meta.items():
        print(f'  {k}: {m["rows"]} x {m["cols"]} = {m["n"]}')
    print('sum:', sum(m['n'] for m in meta.values()))
    for k, v in vivo.items():
        print(f'  {k}: {len(v["names"])} groups x {len(v["replicates"])} replicates')
        print('   ', {nm: round(mv, 1) for nm, mv in v['mean'].items()})


if __name__ == '__main__':
    main()
