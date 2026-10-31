import xml.etree.ElementTree as ET, glob, re, json, os
D=os.path.dirname(os.path.abspath(__file__))
NUM=r'(\d{1,3}(?:,\d{3})+|\d{2,6})'
# phrases that genuinely denote a library/screen size
PATS=[
 (re.compile(r'librar(?:y|ies)\s+(?:of|containing|comprising|consisting of|with|totall?ing)\s+(?:over\s+|more than\s+|>\s*|~\s*|nearly\s+|approximately\s+|up to\s+|some\s+)?'+NUM,re.I),'library of N'),
 (re.compile(r'(?:screen(?:ed|ing|s)?|synthesi[sz]ed|generated|constructed|evaluated|tested|assessed|prepared|produced|created|curated|compiled|assembled|designed)\s+(?:a\s+(?:combinatorial\s+|chemically\s+diverse\s+|diverse\s+)?librar\w+\s+of\s+|data\s+from\s+|a\s+total\s+of\s+)?(?:over\s+|more than\s+|>\s*|~\s*|nearly\s+|approximately\s+|up to\s+)?'+NUM+r'(?=[\s,]*(?:distinct|unique|novel|different|new|structurally|chemically|candidate|ionizable|cationic|amino|degradable|biodegradable|lipid|lipids|lipidoid|lipidoids|LNP|LNPs|nanoparticle|nanoparticles|formulation|formulations|compound|compounds|member|members))',re.I),'verb N'),
 (re.compile(NUM+r'[\s\-]+(?:(?:structurally|chemically)[\s\-]+\w+[\s]+)?(?:novel\s+|distinct\s+|unique\s+|different\s+|new\s+|candidate\s+|diverse\s+)?(?:ionizable\s+|cationic\s+|amino\s+|degradable\s+|biodegradable\s+|aromatic\s+|branched[\s\-]+tail\s+)?(?:lipids|lipidoids|LNPs|lipid nanoparticles|formulations|nanoparticles)\b',re.I),'N lipids'),
 (re.compile(NUM+r'[\s\-]*member',re.I),'N-member'),
 (re.compile(r'(?:comprising|consisting of|spanning|covering|across)\s+(?:over\s+|more than\s+|>\s*|~\s*)?'+NUM+r'(?=[\s,]*(?:ionizable|lipid|lipids|lipidoid|LNP|LNPs|formulation|nanoparticle))',re.I),'spanning N'),
]
# numbers that are parts of known lipid / reagent names -> reject
BLOCK=re.compile(r'(?:ALC[\s\-]?0?|SM[\s\-]?|C12[\s\-]|C14[\s\-]|5A2[\s\-]SC|MC|KC|DLin[\s\-]?|LP[\s\-]?0?|cKK[\s\-]E|7C|306Oi|200|NTLA[\s\-]|BNT|mRNA[\s\-]1273|304O|A18B8C|O[\s\-])$',re.I)
CTX=re.compile(r'lipid|LNP|nanoparticle|lipidoid|formulation',re.I)
BAD=re.compile(r'sphingolipid|metabolomic|proteins?\b.*structur|phytochem|detergent screening|fermented|asthma|transporter|GPR|bacteri|fatty acid|untargeted metabolom',re.I)
recs=[]
for fn in sorted(glob.glob(os.path.join(D,'xml','*.xml'))):
    for a in ET.parse(fn).iter('PubmedArticle'):
        pmid=a.findtext('./MedlineCitation/PMID')
        el=a.find('.//ArticleTitle'); ti=''.join(el.itertext()) if el is not None else ''
        j=a.findtext('.//Journal/ISOAbbreviation') or a.findtext('.//Journal/Title') or ''
        yr=a.findtext('.//JournalIssue/PubDate/Year') or (a.findtext('.//JournalIssue/PubDate/MedlineDate') or '')[:4]
        aus=[f"{x.findtext('LastName')} {(x.findtext('Initials') or '')}".strip() for x in a.iter('Author') if x.findtext('LastName')]
        ab=' '.join(''.join(x.itertext()) for x in a.iter('AbstractText'))
        pt=[x.text for x in a.iter('PublicationType')]
        doi=''; pd=a.find('./PubmedData/ArticleIdList')
        if pd is not None:
            for i in pd.findall('ArticleId'):
                if i.get('IdType')=='doi': doi=i.text or ''
        hits=[]
        for p,lab in PATS:
            for m in p.finditer(ab):
                g=[x for x in m.groups() if x and re.fullmatch(r'[\d,]+',x)]
                if not g: continue
                raw=g[0]; v=int(raw.replace(',',''))
                if not (20<=v<=200000): continue
                pre=ab[max(0,m.start(1)-10):m.start(1)]
                if BLOCK.search(pre): continue
                s=max(0,m.start()-110); e=min(len(ab),m.end()+110); win=ab[s:e]
                if CTX.search(win): hits.append((v,lab,win.strip()))
        mx=max([h[0] for h in hits],default=0)
        recs.append(dict(pmid=pmid,title=ti,journal=j,year=yr,authors=aus,first=aus[0] if aus else '',last=aus[-1] if aus else '',
                         doi=doi,abstract=ab,maxlib=mx,hits=sorted(set(hits),reverse=True)[:5],
                         review='Review' in pt, offtopic=bool(BAD.search(ti+' '+ab))))
json.dump(recs,open(os.path.join(D,'recs2.json'),'w'))
cand=[r for r in recs if r['maxlib']>=40 and not r['review'] and not r['offtopic']]
cand.sort(key=lambda r:-r['maxlib'])
print('records',len(recs),'| clean candidates >=40:',len(cand))
for r in cand:
    print(f"{r['maxlib']:>7}  {r['pmid']}  {r['year']}  {r['journal'][:23]:23s}  {r['last'][:15]:15s}  {r['title'][:72]}")
