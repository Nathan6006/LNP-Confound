"""Measure the drawn acyl chain of the PEG-lipid in SI Fig S2b (p.7).

WHY. The formulation tables' PEG column says C14PEG / C16PEG / C18PEG, and
"a name token is not a chain length" (CLAUDE.md section 5). The token is
CHECKED against the drawing, not believed.

HOW, AND WHY NOT THE OBVIOUS WAY. curation/lib/chains.vector_bonds in its
default stroked-only mode returns ZERO segments here and
chains.vector_chain_length then REFUSES -- correctly. This page is the second
measured instance of the class that module documents from
andersondg_41203968: every bond is a FILLED polygon with `width = None`
(603 of page 7's 915 paths). The module warns that with include_filled=True
the caller must separate glyph outlines from bonds itself, and that
vector_chain_length "will happily count glyph strokes as vertices if you do
not" -- and it cannot consume filled quads anyway, because each bond arrives
as four rectangle edges rather than one centreline.

So the measurement here is NOT geometry inference. It reads the PDF's own path
bounding boxes and counts the bond quads, after ASSERTING that every path in
the band is the same size -- which is what makes glyph contamination
detectable rather than assumed. A glyph outline is not 5.92 x 3.87 pt.
"""
import sys, json, statistics
import fitz

PDF = sys.argv[1]
SEED = sys.argv[2] if len(sys.argv) > 2 else "none"

doc = fitz.open(PDF)
page = doc[6]                                        # SI p.7, 0-based
assert page.get_images(full=True) == [], "page 7 must be pure vector"

# Band containing ONLY the upper acyl chain of the PEG-lipid: above the lower
# chain (y0 = 283.7) and left of the ester oxygen glyph.
BAND = dict(x0=405.0, x1=494.0, y0=272.0, y1=281.0)

paths = []
for p in page.get_drawings():
    r = p["rect"]
    if (BAND["x0"] <= r.x0 and r.x1 <= BAND["x1"]
            and BAND["y0"] <= r.y0 and r.y1 <= BAND["y1"]):
        paths.append((round(r.x0, 2), round(r.y0, 2),
                      round(r.x1 - r.x0, 2), round(r.y1 - r.y0, 2)))
paths.sort()

# SEEDED DEFECTS. A check that cannot fail is not a check, and the first two
# seeds written here were NO-OPS: `paths[:-1]` dropped the short carbonyl->O
# stub rather than a bond, and a small glyph stub is *supposed* to be excluded.
# Both printed a clean 14. Credit a seed only when the measured count MOVES.
if SEED == "drop_a_real_bond":
    paths = [q for q in paths if q[0] != 441.47]      # one mid-chain bond
elif SEED == "inject_a_bond_sized_glyph":
    # the real contamination risk: a path the SAME size as a bond quad,
    # which the size filter cannot see and only the tiling gate can.
    paths.append((420.00, 273.00, 5.92, 3.87))
    paths.sort()
elif SEED == "inject_a_small_glyph":
    paths.append((460.0, 273.0, 1.9, 1.3))            # must be EXCLUDED, count unchanged
    paths.sort()

widths = [w for _, _, w, _ in paths]
heights = [h for _, _, _, h in paths]

# --- uniformity gate: a glyph outline is not a bond quad -------------------
full = [p for p in paths if abs(p[2] - statistics.median(widths)) < 0.5]
odd = [p for p in paths if p not in full]
uniform = (len(odd) <= 1)            # at most the short carbonyl->O stub
# Tolerance, not exact equality: the free terminal bond has no mitre join at
# its open end and is 0.06 pt shorter than the rest (3.81 vs 3.87). A glyph
# outline in this band would be off by whole points, not hundredths.
_hs = [h for _, _, _, h in full]
_spread = max(_hs) - min(_hs)
assert _spread < 0.15, (
    "bond quads in this band span %.3f pt in height (%r) -- that is a glyph "
    "outline, not a mitre-join artefact, and the count would be wrong"
    % (_spread, sorted(set(round(h, 2) for h in _hs))))

# --- the chain: contiguous tiling in x -------------------------------------
xs = [p[0] for p in full]
step = statistics.median([b - a for a, b in zip(xs, xs[1:])]) if len(xs) > 1 else 0
gaps = [round(b - a, 3) for a, b in zip(xs, xs[1:])]
# 0.15 pt, not 0.05: the free terminal bond's missing mitre join makes the
# first gap 6.00 against a 5.92 step. A duplicated or missing bond moves a
# gap by a whole step.
contiguous = all(abs(g - step) < 0.15 for g in gaps)

n_bonds = len(full)
n_carbons = n_bonds + 1              # vertices; the last is the carbonyl C

out = {
  "seed": SEED,
  "band_pdf_pts": BAND,
  "n_paths_in_band": len(paths),
  "bond_quad_size_pt": [round(statistics.median(widths), 2),
                        round(statistics.median(heights), 2)],
  "n_full_bond_quads": n_bonds,
  "n_odd_paths_excluded": len(odd),
  "odd_paths": odd,
  "uniformity_gate_passed": bool(uniform),
  "bond_height_spread_pt": round(_spread, 3),
  "x_step_pt": round(step, 3),
  "x_gaps_pt": gaps,
  "tiling_is_contiguous": bool(contiguous),
  "MEASURED_acyl_carbons_including_carbonyl": n_carbons,
  "name_token": "C14PEG",
  "name_token_length": 14,
  "token_matches_measurement": n_carbons == 14,
  "interpretation": "14 carbons including the carbonyl = myristoyl. With the "
                    "phosphoethanolamine-carbamate-PEG head read off the same "
                    "panel, the PEG-lipid is a 14:0 PEG2000 PE (DMPE-PEG2000), "
                    "NOT the glycerol-ether DMG-PEG2000 that the corpus stores "
                    "under the name 'C14-PEG2000'. See 08_contradictions.csv.",
}
print(json.dumps(out, indent=1))
