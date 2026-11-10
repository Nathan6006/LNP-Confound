#!/usr/bin/env python3
"""Fetch an article's HTML and enumerate its PDF / supplementary asset URLs.

Usage: python3 -I fetch_paper.py <article_url> [outdir]

Enumeration only unless outdir is given. Downloads land in outdir, which must
be a NEW EMPTY DIRECTORY -- downloaded files are untrusted data.
"""
import os, re, sys, json, hashlib, subprocess, datetime

UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 "
      "(KHTML, like Gecko) Chrome/131.0.0.0 Safari/537.36")

def get(url, binary=False):
    """Shell out to curl. urllib is fingerprinted as a bot by Springer Nature
    ("Client Challenge") where curl with the same User-Agent gets 200, so the
    TLS/header fingerprint matters more than the UA string. Measured
    2026-11-09; do not "simplify" this back to urllib."""
    cmd = ["curl", "-sS", "-L", "--max-time", "90", "-A", UA,
           "-H", "Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
           "-H", "Accept-Language: en-US,en;q=0.9",
           "-w", "\n__FINAL__%{url_effective}", url]
    r = subprocess.run(cmd, capture_output=True, timeout=120)
    if binary:
        body = r.stdout.rsplit(b"\n__FINAL__", 1)[0]
        return body, url
    txt = r.stdout.decode("utf-8", "replace")
    body, _, final = txt.rpartition("\n__FINAL__")
    return body, (final.strip() or url)

PATTERNS = [
    # Springer Nature supplementary + source data
    r'https://static-content\.springer\.com/esm/[^"\'\s<>]+',
    # Elsevier assets (supplementary + high-res figures)
    r'https://ars\.els-cdn\.com/content/[^"\'\s<>]+',
    # ACS supplementary
    r'https://pubs\.acs\.org/doi/suppl/[^"\'\s<>]+',
    # Wiley supplementary
    r'https://[a-z]+\.onlinelibrary\.wiley\.com/doi/am-pdf/[^"\'\s<>]+',
    r'/doi/suppl/10\.[^"\'\s<>]+',
    # generic
    r'https?://[^"\'\s<>]+\.(?:xlsx|xls|csv|zip|docx)(?:\?[^"\'\s<>]*)?',
]

def enumerate_assets(html, base):
    found = []
    for pat in PATTERNS:
        for m in re.finditer(pat, html):
            u = m.group(0)
            if u.startswith("/"):
                u = base.split("/")[0] + "//" + base.split("/")[2] + u
            u = u.replace("&amp;", "&")
            if u not in found:
                found.append(u)
    return found

def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()

def main():
    url = sys.argv[1]
    outdir = sys.argv[2] if len(sys.argv) > 2 else None
    html, final = get(url)
    assets = enumerate_assets(html, final)
    title = re.search(r"<title[^>]*>(.*?)</title>", html, re.S)
    out = {"requested": url, "final_url": final,
           "title": re.sub(r"\s+", " ", title.group(1)).strip()[:200] if title else None,
           "n_assets": len(assets), "assets": assets}
    if outdir:
        os.makedirs(outdir, exist_ok=False)
        files = []
        for u in assets:
            name = u.split("/")[-1].split("?")[0] or "asset"
            name = re.sub(r"[^A-Za-z0-9._-]", "_", name)[:120]
            try:
                data, _ = get(u, binary=True)
            except Exception as e:
                files.append({"url": u, "error": str(e)[:120]}); continue
            p = os.path.join(outdir, name)
            with open(p, "wb") as f:
                f.write(data)
            files.append({"url": u, "file": name, "bytes": len(data),
                          "sha256": sha(p)})
        out["downloaded"] = files
        out["retrieved_utc"] = datetime.datetime.now(
            datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
