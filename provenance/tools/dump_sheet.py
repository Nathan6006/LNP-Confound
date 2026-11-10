import sys
from openpyxl import load_workbook
wb = load_workbook(sys.argv[1], read_only=True, data_only=True)
for ws in wb.worksheets:
    if len(sys.argv) > 2 and ws.title != sys.argv[2]:
        continue
    print("### sheet", repr(ws.title))
    for i, row in enumerate(ws.iter_rows(values_only=True)):
        if i > 12:
            print("   ...")
            break
        print("  ", [("" if c is None else (round(c,3) if isinstance(c,float) else c)) for c in row][:12])
wb.close()
