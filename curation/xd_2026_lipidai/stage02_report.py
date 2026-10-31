"""Vertex + glyph report per structure crop.

Chains are reported as merged polylines (bond length ~54 px in these crops);
glyph components (letters) are reported separately so heteroatom labels can be
placed. Endpoints that sit <0.6 bond from a glyph are bond stubs, not atoms.
"""
import json, sys, subprocess, math
D='/Users/nathanliu/LNP/confound/curation/xd_2026_lipidai'
P=D+'/00_pages/png'

def merge(pl, tol):
    out=[tuple(pl[0])]
    for p in pl[1:]:
        p=tuple(p)
        if math.dist(p,out[-1])<tol: out[-1]=((p[0]+out[-1][0])//2,(p[1]+out[-1][1])//2)
        else: out.append(p)
    return out

def run(name, tol=22):
    subprocess.run(['python3','-I',D+'/stage02_trace.py',f'{P}/str_{name}.png',
                    f'{D}/02_{name}.json',f'{P}/ov_{name}.png'],capture_output=True)
    comps=json.load(open(f'{D}/02_{name}.json'))
    print(f'===== {name} =====')
    glyphs=[]; chains=[]
    for c in comps:
        r0,c0,r1,c1=c['bbox']; w,h=c1-c0,r1-r0
        segs=[]
        for pl in c['polylines']:
            if len(pl)<2: continue
            m=merge(pl,tol)
            span=max(x for x,y in m)-min(x for x,y in m)
            vspan=max(y for x,y in m)-min(y for x,y in m)
            if max(span,vspan)<45: continue
            segs.append(m)
        if w<48 and h<48:
            glyphs.append(((c0+c1)//2,(r0+r1)//2,w,h,c['area'])); continue
        if not segs:
            glyphs.append(((c0+c1)//2,(r0+r1)//2,w,h,c['area'])); continue
        for m in segs:
            d=[round(math.dist(m[i],m[i+1])) for i in range(len(m)-1)]
            chains.append((m,d))
    print(' GLYPHS (cx,cy,w,h,area):')
    for g in sorted(glyphs): print('   ',g)
    print(' CHAINS:')
    for m,d in sorted(chains,key=lambda t:-len(t[0])):
        print(f'   n={len(m)} bondlen={d}')
        print(f'      {m}')
    print()

for n in sys.argv[1:]: run(n)
