#!/usr/bin/env python3
"""Probe, without any login, whether each candidate's data can actually be reached.

Two independent routes, because they fail independently:

  SI    Many publishers paywall the article but serve supplementary files openly
        (Wiley, ACS, RSC and MDPI often do; Nature's static-content host frequently
        does). The earlier grading only checked Europe PMC, so it under-reports.
        Here: Europe PMC supplementaryFiles, then the publisher landing page parsed
        for SI links, then those links HEAD-ed to see if bytes come back.

  PRE   bioRxiv / ChemRxiv / Research Square preprints carry the same SI, free. Found
        via Crossref posted-content and Europe PMC's PPR source, matched back to the
        paper by title similarity so a same-topic different-paper hit is not counted.

Writes JSON so the result can be merged into data.json. Nothing here needs credentials;
anything that 403s is reported as such rather than retried, since the fix for those is
the user's own institutional session, not a different user-agent.

Usage: python3 -I probe_access.py <in.json> <out.json> [limit]
"""
import difflib
import json
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import zipfile

UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/120.0 Safari/537.36")
HDRS = {"User-Agent": UA, "Accept": "*/*"}

EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
CROSSREF = "https://api.crossref.org/works"

# file extensions that could hold a table of per-lipid values
DATA_EXT = re.compile(r"\.(xlsx?|csv|tsv|zip|txt|json)(\?|$)", re.I)
DOC_EXT = re.compile(r"\.(pdf|docx?|pptx?)(\?|$)", re.I)
# Publisher SI URL shapes. The first pass missed two big ones: Elsevier names its files
# mmc1/mmc2 with no "supplement" anywhere in the path, and Wiley's
# action/downloadSupplement?doi=...&file=... carries no file extension at all.
SI_HINT = re.compile(r"(supplement|supporting|supinfo|suppl|si[_-]|esm|moesm|mmc\d|"
                     r"additional[_-]file|appendix|data[_-]s\d|table[_-]s\d|/bin/)", re.I)


def fetch(url, method="GET", timeout=45, maxbytes=400_000):
    """Return (status, final_url, headers, body_prefix). Never raises."""
    try:
        req = urllib.request.Request(url, headers=HDRS, method=method)
        with urllib.request.urlopen(req, timeout=timeout) as r:
            body = b"" if method == "HEAD" else r.read(maxbytes)
            return r.status, r.geturl(), dict(r.headers), body
    except urllib.error.HTTPError as e:
        return e.code, url, dict(getattr(e, "headers", {}) or {}), b""
    except Exception as e:
        return 0, url, {"_err": type(e).__name__}, b""


def norm(s):
    return " ".join(re.sub(r"[^a-z0-9 ]", " ", (s or "").lower()).split())


def epmc_pmcid(pmid):
    """data.json carries no PMCID, so look it up - it is the key for the SI endpoint."""
    q = urllib.parse.urlencode({"query": f"EXT_ID:{pmid}", "format": "json",
                                "pageSize": "1", "resultType": "lite"})
    st, _u, _h, body = fetch(f"{EPMC}/search?{q}", maxbytes=100_000)
    if st != 200 or not body:
        return ""
    try:
        res = json.loads(body)["resultList"]["result"]
        return (res[0].get("pmcid") or "") if res else ""
    except Exception:
        return ""


def epmc_suppl(pmcid, outdir):
    """Download Europe PMC's supplementary archive and list what is inside it.

    Listing beats inspecting the article HTML: it names every file and its size, so a
    spreadsheet of per-lipid values can be told from a methods PDF without guessing.
    Each archive goes into its own directory and is only LISTED, never extracted - these
    are untrusted downloads.

    Returns (bytes, note, [{name, size, is_data}]).
    """
    st, _u, h, body = fetch(f"{EPMC}/{pmcid}/supplementaryFiles", maxbytes=80_000_000)
    if st != 200:
        return 0, f"HTTP {st}", []
    n = len(body)
    # a ~300-byte zip is a VALID BUT EMPTY archive - the trap that made the first pass
    # report 24 papers as having SI when most had nothing
    if n < 1200:
        return n, "empty archive", []
    d = os.path.join(outdir, pmcid)
    os.makedirs(d, exist_ok=True)
    p = os.path.join(d, "suppl.zip")
    with open(p, "wb") as f:
        f.write(body)
    try:
        with zipfile.ZipFile(p) as z:
            items = [{"name": i.filename, "size": i.file_size,
                      "is_data": bool(DATA_EXT.search(i.filename))}
                     for i in z.infolist() if not i.is_dir()]
    except zipfile.BadZipFile:
        return n, "not a zip", []
    return n, f"{len(items)} file(s)", sorted(items, key=lambda x: -x["size"])


def pmc_si(pmcid):
    """SI files declared in the PMC full text, which is free even for non-OA deposits.

    More reliable than the supplementaryFiles zip: the XML names each file, so the exact
    extension is visible and a spreadsheet can be told from a PDF before downloading.
    """
    st, _u, _h, body = fetch(f"{EPMC}/{pmcid}/fullTextXML", maxbytes=3_000_000)
    body = body.lstrip()
    if st != 200 or not body.startswith(b"<"):
        return []
    xml = body.decode("utf-8", "replace")
    files = set()
    for m in re.finditer(r"<supplementary-material[^>]*>(.*?)</supplementary-material>",
                         xml, re.S):
        for h in re.finditer(r'xlink:href\s*=\s*["\']([^"\']+)["\']', m.group(1)):
            files.add(h.group(1))
    for m in re.finditer(r'<supplementary-material[^>]*xlink:href\s*=\s*["\']([^"\']+)["\']',
                         xml):
        files.add(m.group(1))
    return [f"https://europepmc.org/articles/{pmcid}/bin/{f}" for f in sorted(files)]


def landing_si(doi):
    """Resolve the DOI, fetch the landing page, return candidate SI links."""
    st, final, _h, body = fetch(f"https://doi.org/{urllib.parse.quote(doi)}")
    if st != 200 or not body:
        return st, final, []
    html = body.decode("utf-8", "replace")
    links = set()
    for m in re.finditer(r'(?:href|content)\s*=\s*["\']([^"\']+)["\']', html):
        u = m.group(1)
        if not SI_HINT.search(u):
            continue
        # no extension requirement: Wiley's downloadSupplement endpoint has none
        links.add(urllib.parse.urljoin(final, u))
    return st, final, sorted(links)[:12]


def mdpi_si(doi):
    """MDPI serves every paper's SI at a fixed path."""
    return [f"https://www.mdpi.com/article/{doi}/s1"]


def check_links(urls):
    out = []
    for u in urls:
        st, _f, h, _b = fetch(u, method="HEAD", timeout=30)
        if st in (403, 405, 0):  # some hosts refuse HEAD; try a ranged GET
            st, _f, h, _b = fetch(u, timeout=30, maxbytes=2048)
        n = h.get("Content-Length")
        out.append({
            "url": u, "status": st,
            "bytes": int(n) if (n or "").isdigit() else None,
            "type": (h.get("Content-Type") or "").split(";")[0],
            "data_file": bool(DATA_EXT.search(u)),
        })
        time.sleep(0.3)
    return out


def preprints(title, doi, last_author=""):
    """Crossref posted-content + Europe PMC PPR.

    Matched on title BUT at a deliberately loose threshold, and also by author: a
    preprint is frequently retitled before publication, so title-only matching at a
    strict threshold misses real hits. Anything at or above LOOSE is reported with its
    similarity score for a human to confirm rather than silently accepted.
    """
    STRICT, LOOSE = 0.80, 0.55
    hits = []
    t = norm(title)
    q = urllib.parse.urlencode({"query.bibliographic": title[:300],
                                "filter": "type:posted-content", "rows": "20",
                                "select": "DOI,title,created,institution,container-title"})
    st, _u, _h, body = fetch(f"{CROSSREF}?{q}", maxbytes=200_000)
    if st == 200 and body:
        try:
            for it in json.loads(body)["message"]["items"]:
                cand = norm(" ".join(it.get("title") or []))
                sim = difflib.SequenceMatcher(None, t, cand).ratio()
                if sim >= LOOSE and (it.get("DOI") or "").lower() != (doi or "").lower():
                    inst = it.get("institution") or []
                    hits.append({"src": "crossref", "doi": it["DOI"], "sim": round(sim, 3),
                                 "title": cand[:110],
                                 "confident": sim >= STRICT,
                                 "server": (inst[0].get("name") if inst else
                                            (it.get("container-title") or [""])[0]),
                                 "date": "-".join(str(x) for x in
                                                  (it.get("created", {}).get("date-parts", [[]])[0]))})
        except Exception:
            pass
    time.sleep(0.3)
    # Europe PMC: search the preprint source by AUTHOR, not title, so a retitled
    # preprint is still found; rank what comes back by title similarity.
    queries = [f'TITLE:"{title[:120]}" AND SRC:PPR']
    if last_author:
        surname = last_author.split()[0]
        queries.append(f'AUTH:"{surname}" AND SRC:PPR AND (lipid OR mRNA OR nanoparticle)')
    for qq in queries:
        q2 = urllib.parse.urlencode({"query": qq, "format": "json",
                                     "pageSize": "25", "resultType": "lite"})
        st, _u, _h, body = fetch(f"{EPMC}/search?{q2}", maxbytes=400_000)
        if st != 200 or not body:
            continue
        try:
            for r in json.loads(body)["resultList"]["result"]:
                sim = difflib.SequenceMatcher(None, t, norm(r.get("title"))).ratio()
                if sim >= LOOSE:
                    hits.append({"src": "epmc", "doi": r.get("doi"), "sim": round(sim, 3),
                                 "title": (r.get("title") or "")[:110],
                                 "confident": sim >= STRICT,
                                 "server": r.get("journalTitle") or "preprint",
                                 "id": r.get("id"), "date": r.get("firstPublicationDate")})
        except Exception:
            pass
        time.sleep(0.3)
    # de-duplicate on DOI
    seen, out = set(), []
    for h in sorted(hits, key=lambda x: -x["sim"]):
        k = (h.get("doi") or h.get("id") or "").lower()
        if k and k in seen:
            continue
        seen.add(k)
        out.append(h)
    return out[:4]


def main():
    global DL
    DL = os.path.join(os.path.dirname(os.path.abspath(sys.argv[2])), "si_archives")
    os.makedirs(DL, exist_ok=True)
    rows = json.load(open(sys.argv[1]))["rows"]
    cands = [r for r in rows if r["status"] == "new"]
    if len(sys.argv) > 3:
        cands = cands[:int(sys.argv[3])]
    out = {}
    for i, r in enumerate(cands, 1):
        pmid, doi = r["pmid"], (r.get("doi") or "")
        rec = {"pmid": pmid, "doi": doi, "lib": r["lib"], "prio": r["prio"],
               "grade": r["grade"], "lab": r["lab"], "title": r["title"]}
        print(f"[{i}/{len(cands)}] {pmid} {r['lab'][:22]:22} n={r['lib']}", flush=True)

        # --- route 1: supplementary files, no login
        si = []
        pmcid = r.get("pmcid") or epmc_pmcid(pmid)
        rec["pmcid"] = pmcid
        pmc_links = []
        if pmcid:
            n, note, items = epmc_suppl(pmcid, DL)
            rec["epmc_suppl"] = {"bytes": n, "note": note, "files": items}
            if items:
                nd = sum(1 for x in items if x["is_data"])
                print(f"      SI ARCHIVE {n/1e6:.1f} MB, {len(items)} files, "
                      f"{nd} machine-readable", flush=True)
                for x in items[:6]:
                    flag = "  <-- DATA" if x["is_data"] else ""
                    print(f"         {x['size']/1024:8.0f} KB  {x['name'][:72]}{flag}",
                          flush=True)
        if doi:
            if doi.startswith("10.3390"):
                urls = mdpi_si(doi)
                st, final = 200, f"https://www.mdpi.com/{doi}"
            else:
                st, final, urls = landing_si(doi)
            rec["landing"] = {"status": st, "url": final, "n_si_links": len(urls)}
            si = check_links(sorted(set(urls) | set(pmc_links)))
        elif pmc_links:
            si = check_links(pmc_links)
        rec["si"] = si
        got = [s for s in si if s["status"] == 200]
        rec["si_reachable"] = len(got)
        rec["si_data_files"] = len([s for s in got if s["data_file"]])
        print(f"      landing HTTP {rec.get('landing', {}).get('status')} | "
              f"{len(si)} SI links | {rec['si_reachable']} reachable | "
              f"{rec['si_data_files']} data files", flush=True)

        # --- route 2: preprint with the same SI
        pp = preprints(r["title"], doi, r.get("last") or "")
        rec["preprints"] = pp
        if pp:
            for p in pp:
                print(f"      PREPRINT {p['src']}: {p.get('doi')} "
                      f"({p.get('server')}, {p.get('date')}) sim={p['sim']}", flush=True)

        out[pmid] = rec
        json.dump(out, open(sys.argv[2], "w"), indent=1)
        time.sleep(0.5)

    arch = [v for v in out.values() if (v.get("epmc_suppl") or {}).get("files")]
    adata = [v for v in arch
             if any(x["is_data"] for x in v["epmc_suppl"]["files"])]
    ns = sum(1 for v in out.values() if v["si_reachable"])
    npp = [v for v in out.values() if v["preprints"]]
    nppc = [v for v in npp if any(p.get("confident") for p in v["preprints"])]
    print(f"\n==== {len(out)} probed ====")
    print(f"  SI archive with contents : {len(arch)}")
    print(f"  ...of which hold a machine-readable file : {len(adata)}")
    print(f"  SI links reachable off the landing page  : {ns}")
    print(f"  preprint candidates : {len(npp)} ({len(nppc)} high-confidence)")


if __name__ == "__main__":
    main()
