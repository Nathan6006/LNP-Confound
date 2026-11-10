import sys, fitz
doc = fitz.open(sys.argv[1])
pg = doc[int(sys.argv[2]) - 1]
# report embedded image native dpi before rendering anything
for im in pg.get_images(full=True):
    xref = im[0]
    info = doc.extract_image(xref)
    rects = pg.get_image_rects(xref)
    for r in rects:
        dpi_x = info["width"] / (r.width / 72.0)
        dpi_y = info["height"] / (r.height / 72.0)
        print(f"EMBEDDED xref={xref} {info['width']}x{info['height']}px "
              f"placed {r.width:.1f}x{r.height:.1f}pt -> native {dpi_x:.1f}x{dpi_y:.1f} dpi")
print("n_drawings", len(pg.get_drawings()))
if len(sys.argv) > 4:
    dpi = int(sys.argv[4])
    pg.get_pixmap(dpi=dpi).save(sys.argv[3])
    print("rendered", sys.argv[3], "at", dpi, "dpi")
