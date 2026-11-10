import sys, pdfplumber
with pdfplumber.open(sys.argv[1]) as pdf:
    for i, pg in enumerate(pdf.pages, 1):
        t = pg.extract_text() or ""
        print(f"\n########## PAGE {i}  ({pg.width:.0f}x{pg.height:.0f}) ##########")
        print(t)
