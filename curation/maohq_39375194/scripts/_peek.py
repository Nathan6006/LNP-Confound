import sys, csv, collections
p = sys.argv[1]
with open(p, newline='', encoding='utf-8-sig') as fh:
    rows = list(csv.DictReader(fh))
print("nrows", len(rows))
cols = list(rows[0].keys())
print("ncols", len(cols))
for c in cols:
    vals = [r[c] for r in rows]
    u = sorted(set(vals))
    print(f"{c:28s} nuniq={len(u):5d} ", (u[:8] if len(u)<=8 else u[:6]+['...']))
