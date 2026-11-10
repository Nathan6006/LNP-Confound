import sys
from openpyxl import load_workbook
for p in sys.argv[1:]:
    wb = load_workbook(p, read_only=True, data_only=True)
    print("="*70); print(p.split("/")[-1])
    for ws in wb.worksheets:
        vals = 0
        for row in ws.iter_rows(values_only=True):
            vals += sum(1 for c in row if isinstance(c, (int, float)))
        print(f"  sheet {ws.title!r}: {ws.max_row} rows x {ws.max_column} cols, "
              f"{vals} numeric cells")
    wb.close()
