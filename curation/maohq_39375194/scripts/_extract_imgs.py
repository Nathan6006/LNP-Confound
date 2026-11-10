import sys, pymupdf
doc = pymupdf.open(sys.argv[1]); pg = doc[int(sys.argv[2]) - 1]
items = []
for im in pg.get_images(full=True):
    xref = im[0]
    for r in pg.get_image_rects(xref):
        items.append((round(r.y0, 1), round(r.x0, 1), xref, r))
items.sort()
# text labels in the left column, to pair with rows
words = pg.get_text("words")
for y0, x0, xref, r in items:
    labs = [w[4] for w in words if abs((w[1]+w[3])/2 - (r.y0+r.y1)/2) < r.height/2 + 6 and w[2] < 200]
    info = doc.extract_image(xref)
    out = f"{sys.argv[3]}/p{sys.argv[2]}_x{xref}.png"
    open(out, "wb").write(info["image"])
    print(f"y0={y0:7.1f} xref={xref} {info['width']}x{info['height']} -> {out}   LABEL={' '.join(labs)}")
