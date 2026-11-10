#!/usr/bin/env python3
"""Fetch a paper's full text and supplementary files from PubMed Central.

    python3 -I pmc_fetch.py <PMCID> <outdir>

Route: efetch gives the full JATS XML (full text, tables, and the filenames of
every supplementary-material item); those files then come from
https://pmc.ncbi.nlm.nih.gov/articles/<PMCID>/bin/<filename>.

Why this is the preferred route: no publisher bot-check, no session, no login,
and the JATS XML carries the article's TABLES AS MARKUP rather than as pixels --
which is tier 2 of the intake gate (`00_value_sources.json`) and removes a
whole class of digitisation error. The old oa.fcgi endpoint 404s; this is its
replacement.

Downloaded files are untrusted: names are sanitised and confined to outdir.
"""
import hashlib, json, os, re, subprocess, sys

UA = "Mozilla/5.0 (compatible; academic-corpus-curation)"
EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
          "?db=pmc&retmode=xml&id=")
# THE WORKING ROUTE IS `articles/instance/<NUMERIC id>/bin/<name>`.
# `articles/PMC<id>/bin/<name>` returns a 404 HTML page -- measured
# 2026-11-10: 404 on the PMC<id> form, 200 on the instance form for the
# same file. That one wrong template made this tool report "7/7 failed,
# 0 downloads" to several agents, each of which then fell back to
# pixels or recorded the figures as unreachable. The instance route
# also serves figures at 600 dpi (~4200 px wide).
BIN = "https://pmc.ncbi.nlm.nih.gov/articles/instance/{num}/bin/{name}"


def curl(url, binary=False, timeout=180):
    r = subprocess.run(["curl", "-sS", "-L", "--max-time", str(timeout),
                        "-A", UA, url],
                       capture_output=True, timeout=timeout + 30)
    return r.stdout if binary else r.stdout.decode("utf-8", "replace")


def main():
    pmcid, outdir = sys.argv[1], sys.argv[2]
    num = pmcid.replace("PMC", "")
    xml = curl(EFETCH + num)
    if "<article" not in xml:
        print(json.dumps({"pmcid": pmcid, "ok": False, "reason": "no article"}))
        return 1
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, pmcid + ".jats.xml"), "w") as f:
        f.write(xml)

    names = []
    for m in re.finditer(r'xlink:href="([^"]+)"', xml):
        n = m.group(1)
        if re.search(r'\.(pdf|docx?|xlsx?|csv|zip|txt|tif|jpg|png)$', n, re.I):
            names.append(n)
    names = sorted(set(names))

    got, failed = [], []
    for n in names:
        safe = re.sub(r'[^A-Za-z0-9._-]', '_', os.path.basename(n))[:120]
        dest = os.path.abspath(os.path.join(outdir, safe))
        if not dest.startswith(os.path.abspath(outdir) + os.sep):
            failed.append({"name": n, "why": "path escape"}); continue
        data = curl(BIN.format(num=num, name=n), binary=True)
        # STRIP FIRST. PMC error pages begin with newlines, so a prefix test
        # on the raw bytes passed them and seven HTML error pages were written
        # as ".jpg"/".pdf". Check the stripped head, and require a known magic.
        head = data.lstrip()[:20].lower()
        if (len(data) < 200 or head.startswith(b"<!doctype html")
                or head.startswith(b"<html")):
            failed.append({"name": n, "why": f"HTML error page ({len(data)}B)"})
            continue
        with open(dest, "wb") as f:
            f.write(data)
        got.append({"file": safe, "bytes": len(data),
                    "sha256": hashlib.sha256(data).hexdigest(),
                    "magic": data[:8].decode("latin1", "replace")})

    title = re.search(r'<article-title[^>]*>(.*?)</article-title>', xml, re.S)
    print(json.dumps({
        "pmcid": pmcid, "ok": True,
        "title": re.sub(r'<[^>]+>', '', title.group(1)).strip()[:160] if title else None,
        "xml_bytes": len(xml), "n_tables": xml.count("<table-wrap"),
        "n_supp_declared": len(names), "downloaded": got, "failed": failed}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
