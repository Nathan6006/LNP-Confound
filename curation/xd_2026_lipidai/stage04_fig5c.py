"""Digitise Fig 5C (log-scale bar chart, 12 LNPs) from the 600-dpi page render.

Calibration uses the y-axis MAJOR tick marks (10^6..10^9), located as the long
horizontal ticks left of the axis spine. Replicate dots are read separately so
the bar value can be cross-checked against the spread of its own points.

BAR TOPS: FLANK COLUMNS, not a width-fraction rule. Corrected 2026-10-30,
review finding F-020.

The old rule declared the top at the first row where more than 60% of the inner
bar width is fill. The three replicate discs straddle the top edge and occlude
the MIDDLE of the width, so that rule fires several pixels INSIDE the bar. At
the true top row the inner fill fraction is 0.21-0.57 on the nine affected bars
and 0.68-0.97 on the three that are not -- the 0.60 threshold is exactly the
discriminator, so the rule IS the error rather than merely correlating with it.
Re-applying it reproduces the stored tops on 12 of 12 at exactly 0 px.

The bias is one-sided and not constant, so it does not cancel in the z-score:
every one of the 12 `Experiment_value`s moves (rms |dz| 0.034, max 0.058 on E4),
including the three whose pixel top does not change, because the block mean and
sd shift under them. Mean bias +0.051 decades, max +0.117 (E4).

The fix reads the FLANK columns, which a centred disc cannot reach, and requires
the left and right flanks to agree. Note the trap: a global `min` over all
columns is NOT safe -- on E1 the anti-aliased halo of a disc floating above the
bar makes columns 23-25 report a top 8 px too high.

UNITS. The calibration is 179.96 px/decade at this 600-dpi render and 89.98
px/decade in the figure's native 300 dpi (the page embeds a single 1541x1864
PNG that the render upsamples 2x). Pixel offsets quoted below are 600-dpi.
The calibration itself is sound -- 179.960 px/decade, major-tick residuals
<=4e-4 decades, validated against 31 unlabelled minor ticks at rms 0.53 px --
so this is NOT error class 20, and there is only one error here, not two.
"""
from collections import Counter

import numpy as np, json
from PIL import Image

P='/Users/nathanliu/LNP/confound/curation/xd_2026_lipidai/00_pages/png/p13_600.png'
im=np.array(Image.open(P).convert('RGB')).astype(int)
# Fig 5C plot region (found by inspection of the 600-dpi render)
Y0,Y1,X0,X1 = 4850, 5620, 900, 2000
sub=im[Y0:Y1, X0:X1]

r,g,b = sub[:,:,0], sub[:,:,1], sub[:,:,2]
fill = (abs(r-205)<28)&(abs(g-205)<28)&(b>r+12)&(b<255)      # lavender bar fill
dark = (r<110)&(g<110)&(b<110)                                # ink

# --- axis spine: the longest vertical dark run
colcount = dark.sum(0)
spine = int(np.argmax(colcount))
print('spine x (crop) =', spine, 'height', colcount[spine])

# --- major ticks: dark pixels immediately left of the spine, long runs
tick_rows=[]
for y in range(dark.shape[0]):
    run = dark[y, max(0,spine-15):spine].sum()
    if run >= 10: tick_rows.append(y)
# group
groups=[];cur=[tick_rows[0]]
for y in tick_rows[1:]:
    if y-cur[-1]<=2: cur.append(y)
    else: groups.append(sum(cur)/len(cur)); cur=[y]
groups.append(sum(cur)/len(cur))
print('major tick rows:', groups)

assert len(groups)==4, groups
decades=[9,8,7,6]
ys=np.array(groups); ex=np.array(decades,float)
A=np.polyfit(ys,ex,1)
print('calibration: log10(value) = %.8f*y + %.6f'%(A[0],A[1]))
resid=ex-np.polyval(A,ys); print('tick residuals (decades):', np.round(resid,4))

# --- bars
colfill = fill.sum(0)
on = colfill>40
segs=[];s=None
for i,v in enumerate(on):
    if v and s is None: s=i
    if not v and s is not None:
        if i-s>18: segs.append([s,i])
        s=None
merged=[segs[0]]
for a,b2 in segs[1:]:
    if a-merged[-1][1] < 12: merged[-1][1]=b2
    else: merged.append([a,b2])
segs=[tuple(x) for x in merged]
print('n bars:', len(segs))
names=['ALC-0315','SM-102','MC3','E1','E2','E3','E4','E5','E6','E7','E8','E9']
assert len(segs)==12
def flank_top(x0, x1, n_edge=5, inset=1, run=6):
    """Top row from the bar's left and right flank columns, which the replicate
    discs never cover. Requires the two flanks to agree within 1 px -- if they
    do not, the bar edge itself is not what is being measured."""
    def side(cols):
        tops = []
        for x in cols:
            c = fill[:, x]
            # first row starting a run of `run` consecutive fill pixels, so a
            # stray anti-aliased pixel or a disc halo cannot set the top
            for y in range(len(c) - run):
                if c[y:y + run].all():
                    tops.append(y)
                    break
        return Counter(tops).most_common(1)[0][0] if tops else None
    L = side(range(x0 + inset, x0 + inset + n_edge))
    R = side(range(x1 - inset - n_edge, x1 - inset))
    assert L is not None and R is not None, f'no flank fill for bar x[{x0},{x1}]'
    assert abs(L - R) <= 1, (f'bar x[{x0},{x1}]: left flank {L} and right flank '
                             f'{R} disagree by {abs(L-R)} px, so the top edge is '
                             f'not what is being measured')
    return min(L, R)


out=[]
for nm,(x0,x1) in zip(names,segs):
    col=fill[:, x0+6:x1-6]
    old_top=int(np.nonzero(col.sum(1) > 0.6*(x1-x0-12))[0].min())
    top=flank_top(x0, x1)
    occluded = col.sum(1)[top] / (x1-x0-12)
    val=10**np.polyval(A,top)
    out.append(dict(name=nm, x0=int(x0+X0), x1=int(x1+X0), top_y=int(top+Y0),
                    value=float(val), log10=float(np.polyval(A,top)),
                    top_y_old_60pct_rule=int(old_top+Y0),
                    px_above_old_rule=int(old_top-top),
                    inner_fill_fraction_at_true_top=round(float(occluded),3),
                    decades_above_old_rule=round(float(
                        np.polyval(A,top)-np.polyval(A,old_top)),4)))
    print(f'{nm:10s} x[{x0:4d},{x1:4d}] top_y={top:4d} '
          f'(60%-rule said {old_top:4d}, {old_top-top:+3d} px, '
          f'fill at true top {occluded:.2f}) value={val:.3e}')
json.dump(dict(calibration=dict(slope=A[0],intercept=A[1],tick_rows=groups,
                                decades=decades,residual_decades=list(np.round(resid,5))),
               bars=out),
          open('/Users/nathanliu/LNP/confound/curation/xd_2026_lipidai/04_fig5c.json','w'),indent=1)
