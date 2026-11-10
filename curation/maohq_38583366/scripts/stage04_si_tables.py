"""Parse Supplementary Tables 1-6 out of mmc1.docx -> tidy CSV."""
import sys, re, csv
from xml.etree import ElementTree as ET
W='{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'
docxml, out = sys.argv[1], sys.argv[2]
root=ET.parse(docxml).getroot(); body=root.find(W+'body')
def ptext(p): return ''.join(t.text or '' for t in p.iter(W+'t'))
caps=[]; tbls=[]
pending=None
for el in body:
    if el.tag==W+'p':
        s=ptext(el).strip()
        if re.match(r'^Supplementary Table \d+\.', s): pending=s
    elif el.tag==W+'tbl':
        rows=[]
        for tr in el.findall(W+'tr'):
            rows.append([' '.join(ptext(p).strip() for p in tc.findall(W+'p')).strip()
                         for tc in tr.findall(W+'tc')])
        tbls.append((pending, rows)); pending=None
assert len(tbls)==6, f"expected 6 tables, got {len(tbls)}"
recs=[]
for cap, rows in tbls:
    m=re.search(r'top 32 (\S+)-based', cap); assert m, cap
    helper=m.group(1)
    hdr0, hdr1 = rows[0], rows[1]
    assert hdr0[0]=='LNP ID' and hdr0[2]=='RLU in HepG2', hdr0
    assert hdr1[1]==helper, (hdr1, helper)
    assert hdr1[1:6]==[helper,'DLin-MC3','Chol.','DMG-PEG','pDNA'], hdr1
    body_rows=rows[2:]
    assert len(body_rows)==32, (cap, len(body_rows))
    for r in body_rows:
        assert len(r)==9, r
        recs.append(dict(si_table=cap.split('.')[0], helper=helper, lnp_id=r[0],
                         wt_helper=r[1], wt_mc3=r[2], wt_chol=r[3], wt_peg=r[4],
                         wt_pdna=r[5], rlu_hepg2=r[6], z_average=r[7], pdi=r[8]))
with open(out,'w',newline='') as fh:
    w=csv.DictWriter(fh, fieldnames=list(recs[0].keys())); w.writeheader(); w.writerows(recs)
print(f"wrote {len(recs)} rows to {out}")
