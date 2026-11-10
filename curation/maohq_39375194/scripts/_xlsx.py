import sys
from openpyxl import load_workbook
wb = load_workbook(sys.argv[1], data_only=True)
for ws in wb.worksheets:
    print("=== SHEET", repr(ws.title), ws.max_row, "x", ws.max_column)
    for i, row in enumerate(ws.iter_rows(values_only=True), 1):
        if i > 8 and i < ws.max_row - 2:
            continue
        print(" ", i, [("" if v is None else v) for v in row])
