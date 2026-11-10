import sys, re, xml.etree.ElementTree as ET
t = ET.parse(sys.argv[1]); r = t.getroot()
def txt(e):
    return re.sub(r'\s+', ' ', ''.join(e.itertext())).strip()
# permissions / license
for p in r.iter('permissions'):
    print("[PERMISSIONS]", txt(p)[:800]); print()
for sec in r.iter('sec'):
    ti = sec.find('title')
    name = txt(ti) if ti is not None else "(untitled)"
    print("\n===== SEC:", name, "=====")
    print(txt(sec)[:6000])
