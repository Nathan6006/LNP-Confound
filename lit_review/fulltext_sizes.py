#!/usr/bin/env python3
"""Resolve library sizes from Europe PMC full text where the abstract does not state one.

For each PMID: ask Europe PMC whether a full text is available, pull it, and print every
sentence that both contains a number in a plausible range and mentions a counted thing
(lipids, LNPs, formulations, library members). Prints sentences, not a verdict - deciding
what a number counts needs a human, because papers routinely state a synthesis count, a
screening count and a formulation count that all differ.

Usage: python3 -I fulltext_sizes.py <pmid> [<pmid> ...]
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request

UA = {"User-Agent": "lit-review-size-check/1.0 (academic use)"}
SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
# The endpoint takes the PMCID alone - a /{source}/{id}/ segment 404s.
FULL = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pid}/fullTextXML"

COUNTED = re.compile(
    r"\b(lipid|lipidoid|lipids|lipidoids|ionizable|LNP|LNPs|nanoparticle|nanoparticles|"
    r"formulation|formulations|library|libraries|member|members|compound|compounds|"
    r"amphiphile|amphiphiles|dendrimer|dendrimers|IAJD|IAJDs|analog|analogue)s?\b", re.I)
# a number, as digits or as a written-out word
NUMTOK = re.compile(r"\b(\d{1,3}(?:,\d{3})+|\d{2,5})\b")
CUE = re.compile(
    r"\b(librar|synthesi[sz]|screen|generat|prepar|construct|comprising|consisting|"
    r"total of|in total|a set of|panel of|array of)", re.I)


def get(url):
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=60) as r:
            return r.read()
    except Exception as e:
        return b"" if "HTTP Error 404" in str(e) else b""


def locate(pmid):
    q = urllib.parse.urlencode({"query": f"EXT_ID:{pmid}", "format": "json",
                                "resultType": "core", "pageSize": "1"})
    raw = get(f"{SEARCH}?{q}")
    if not raw:
        return None
    try:
        res = json.loads(raw)["resultList"]["result"]
    except Exception:
        return None
    if not res:
        return None
    r = res[0]
    return {
        "pmcid": r.get("pmcid"), "oa": r.get("isOpenAccess"),
        "inEPMC": r.get("inEPMC"), "src": r.get("source"), "id": r.get("id"),
        "title": (r.get("title") or "")[:80],
    }


def sentences(xml):
    txt = re.sub(r"<[^>]+>", " ", xml.decode("utf-8", "replace"))
    txt = re.sub(r"&[a-z]+;|&#\d+;", " ", txt)
    txt = re.sub(r"\s+", " ", txt)
    return re.split(r"(?<=[.!?])\s+(?=[A-Z0-9])", txt)


def main():
    for pmid in sys.argv[1:]:
        info = locate(pmid)
        print("=" * 104)
        if not info:
            print(f"{pmid}  NOT FOUND in Europe PMC")
            time.sleep(0.4)
            continue
        print(f"{pmid}  oa={info['oa']} inEPMC={info['inEPMC']} pmcid={info['pmcid']}  {info['title']}")
        xml = b""
        if info["pmcid"]:
            # the body begins with a newline before <!DOCTYPE, so strip before testing
            xml = get(FULL.format(pid=info["pmcid"])).lstrip()
        if not xml.startswith(b"<"):
            print("   no full text served")
            time.sleep(0.4)
            continue
        shown = 0
        seen = set()
        for s in sentences(xml):
            if len(s) > 420 or not CUE.search(s):
                continue
            nums = [int(n.replace(",", "")) for n in NUMTOK.findall(s)]
            nums = [n for n in nums if 10 <= n <= 60000]
            if not nums or not COUNTED.search(s):
                continue
            key = (tuple(sorted(set(nums))), s[:40])
            if key in seen:
                continue
            seen.add(key)
            print(f"   {sorted(set(nums))}  {s.strip()[:300]}")
            shown += 1
            if shown >= 10:
                break
        if not shown:
            print("   full text served, but no size-like sentence matched")
        time.sleep(0.4)


if __name__ == "__main__":
    main()
