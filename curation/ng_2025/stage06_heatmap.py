"""Check the Source Data blocks against the drawn heat maps, cell by cell.

Fig. 1e-j are vector: every cell is a filled rectangle whose colour the PDF
states exactly. If the workbook block were transposed, shifted by a row, or
assigned to the wrong panel, the colour order would not track the numbers.
This maps each rectangle to its (row, column) from the axis label positions and
reports the Spearman correlation between the drawn colour and the number the
workbook gives for that same cell - per panel, so one good panel cannot hide a
bad one.
"""
import json, os, sys, collections
import pymupdf
from scipy.stats import spearmanr

D = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(D, 's41551-025-01422-8.pdf')
VALS = json.load(open(os.path.join(D, '04_values.json')))

# panel -> (block key, x of each tail column, y of each amine row) read from the
# figure's own axis labels (pymupdf text), not from guesswork
PANELS = {
    'Fig.1e': ('e', (305, 395), (168, 292)),
    'Fig.1f': ('f', (414, 505), (168, 292)),
    'Fig.1g': ('g', (85, 190), (314, 436)),
    'Fig.1h': ('h', (196, 300), (314, 436)),
    'Fig.1i': ('i', (301, 405), (314, 436)),
    'Fig.1j': ('j', (405, 515), (314, 436)),
}


def main():
    doc = pymupdf.open(PDF)
    page = doc[2]
    rects = []
    for p in page.get_drawings():
        for it in p['items']:
            if it[0] != 're':
                continue
            r = it[1]
            if p.get('fill') is None:
                continue
            rects.append((r.x0, r.y0, r.x1, r.y1, p['fill']))
    report = {}
    for key, (panel, (x0, x1), (y0, y1)) in PANELS.items():
        sel, seen = [], set()
        for r in rects:
            if not (x0 <= r[0] and r[2] <= x1 and y0 <= r[1] and r[3] <= y1):
                continue
            if not (15 < r[2] - r[0] < 30 and 5 < r[3] - r[1] < 12):
                continue
            k = (round(r[0], 0), round(r[1], 0))
            if k in seen:
                continue
            seen.add(k); sel.append(r)
        xs = sorted({round(r[0], 1) for r in sel})
        ys = sorted({round(r[1], 1) for r in sel})
        # cluster columns / rows
        def group(v, tol=3.0):
            out = [[v[0]]]
            for a in v[1:]:
                if a - out[-1][-1] < tol:
                    out[-1].append(a)
                else:
                    out.append([a])
            return [sum(g) / len(g) for g in out]
        cols, rows = group(xs), group(ys)
        grid = {}
        for r in sel:
            ci = min(range(len(cols)), key=lambda i: abs(cols[i] - r[0]))
            ri = min(range(len(rows)), key=lambda i: abs(rows[i] - r[1]))
            # the scale is a single hue (R=G varying, B pinned at 1), and a
            # DARKER cell is a HIGHER value, so negate to make the check monotone up
            grid[(ri, ci)] = -sum(r[4][:2]) / 2.0
        vals = [v for v in VALS['vitro'] if v['block'] == key]
        tails = sorted({v['tail'] for v in vals}, key=lambda t: [x['tail'] for x in vals].index(t))
        amines = []
        for v in vals:
            if v['amine'] not in amines:
                amines.append(v['amine'])
        got, exp, missing = [], [], 0
        for v in vals:
            ri, ci = amines.index(v['amine']), tails.index(v['tail'])
            if (ri, ci) not in grid:
                missing += 1
                continue
            got.append(grid[(ri, ci)])
            exp.append(v['log10_luminescence'])
        rho = spearmanr(got, exp).statistic if len(got) > 3 else float('nan')
        report[key] = dict(cells_drawn=len(grid), cells_in_workbook=len(vals),
                           matched=len(got), missing=missing,
                           n_cols=len(cols), n_rows=len(rows),
                           spearman_colour_vs_value=round(float(rho), 4))
        print(f'{key}: grid {len(rows)}x{len(cols)}  matched {len(got)}/{len(vals)}  '
              f'rho(colour, workbook value) = {rho:+.4f}')
    json.dump(report, open(os.path.join(D, '06_heatmap_check.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
