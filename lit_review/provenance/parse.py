import xml.etree.ElementTree as ET, glob, re, json, os
recs=[]
NUM=r'(\d{1,3}(?:,\d{3})+|\d{2,6})'
# patterns: number followed by lipid-ish noun, or "library of N"
PATS=[
 re.compile(NUM+r'[\s\-]+(?:structurally[\s\-]+\w+[\s]+)?(?:novel[\s]+|distinct[\s]+|unique[\s]+|different[\s]+|new[\s]+|candidate[\s]+|chemically[\s\-]+\w+[\s]+)?(?:ionizable[\s\-]+|cationic[\s\-]+|amino[\s\-]+|degradable[\s\-]+|biodegradable[\s\-]+)?(lipid|lipids|lipidoid|lipidoids|LNP|LNPs|lipid nanoparticles?|formulations?|nanoparticles?|compounds?|members?|analogues?|analogs?|candidates?)\b',re.I),
 re.compile(r'librar(?:y|ies)\s+(?:of|containing|comprising|consisting of|with)\s+(?:over\s+|more than\s+|>\s*|~\s*|approximately\s+)?'+NUM,re.I),
 re.compile(r'(?:screen(?:ed|ing)?|synthesi[sz]ed|generated|constructed|evaluated|tested|assessed|prepared)\s+(?:a\s+librar\w+\s+of\s+)?(?:over\s+|more than\s+|>\s*|~\s*|approximately\s+|up to\s+)?'+NUM,re.I),
 re.compile(NUM+r'[\s\-]*member',re.I),
]
CTX=re.compile(r'lipid|LNP|nanoparticle|lipidoid',re.I)
for fn in sorted(glob.glob(os.path.join(SP:=os.path.dirname(os.path.abspath(__file__)),'xml','*.xml'))):
    try: t=ET.parse(fn)
    except Exception as e: print('parse fail',fn,e); continue
    for a in t.iter('PubmedArticle'):
        pmid=a.findtext('./MedlineCitation/PMID')
        ti=''.join((a.find('.//ArticleTitle').itertext()) if a.find('.//ArticleTitle') is not None else '')
        j=a.findtext('.//Journal/ISOAbbreviation') or a.findtext('.//Journal/Title') or ''
        yr=a.findtext('.//JournalIssue/PubDate/Year') or (a.findtext('.//JournalIssue/PubDate/MedlineDate') or '')[:4]
        aus=[f"{x.findtext('LastName')} {(x.findtext('Initials') or '')}".strip() for x in a.iter('Author') if x.findtext('LastName')]
        ab=' '.join(''.join(x.itertext()) for x in a.iter('AbstractText'))
        pt=[x.text for x in a.iter('PublicationType')]
        doi=''
        pd=a.find('./PubmedData/ArticleIdList')
        if pd is not None:
            for i in pd.findall('ArticleId'):
                if i.get('IdType')=='doi': doi=i.text or ''
        hits=[]
        for p in PATS:
            for m in p.finditer(ab):
                g=[x for x in m.groups() if x and re.fullmatch(r'[\d,]+',x)]
                if not g: continue
                v=int(g[0].replace(',',''))
                s=max(0,m.start()-90); e=min(len(ab),m.end()+90)
                win=ab[s:e]
                if v>=20 and v<=200000 and CTX.search(win):
                    hits.append((v,win.strip()))
        mx=max([h[0] for h in hits],default=0)
        recs.append(dict(pmid=pmid,title=ti,journal=j,year=yr,first=aus[0] if aus else '',last=aus[-1] if aus else '',
                         doi=doi,abstract=ab,maxlib=mx,hits=sorted(set(hits),reverse=True)[:6],ptypes=pt))
json.dump(recs,open(os.path.join(SP,'recs.json'),'w'))
print('records:',len(recs))
rev=[r for r in recs if 'Review' in r['ptypes']]
print('reviews:',len(rev))
cand=[r for r in recs if r['maxlib']>=50 and 'Review' not in r['ptypes']]
print('candidates maxlib>=50, non-review:',len(cand))
cand.sort(key=lambda r:-r['maxlib'])
for r in cand[:120]:
    print(f"{r['maxlib']:>7}  {r['pmid']}  {r['year']}  {r['journal'][:24]:24s}  {r['last'][:16]:16s}  {r['title'][:78]}")
