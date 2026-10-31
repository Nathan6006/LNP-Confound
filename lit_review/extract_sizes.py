#!/usr/bin/env python3
"""Pull every candidate-size claim out of the abstracts, with context, for manual review.

The `lib` field in data.json came from a regex over abstracts and was hand-checked only
for the priority 1-2 rows. This re-derives it for all 61 candidates and prints each hit
with its surrounding words, so a human can see what the number actually counts. It
deliberately does NOT pick a winner - several abstracts state a lipid count and a
formulation count that differ, and which one matters depends on the row.

Usage: python3 -I extract_sizes.py <pubmed.xml>
"""
import re
import sys
import xml.etree.ElementTree as ET

WORDS = {
    # what is being counted -> the bucket it belongs in
    "lipid": "lipid", "lipids": "lipid", "lipidoid": "lipid", "lipidoids": "lipid",
    "ionizable": "lipid", "amphiphile": "lipid", "amphiphiles": "lipid",
    "dendrimer": "lipid", "dendrimers": "lipid", "iajd": "lipid", "iajds": "lipid",
    "lnp": "lnp", "lnps": "lnp", "nanoparticle": "lnp", "nanoparticles": "lnp",
    "formulation": "lnp", "formulations": "lnp", "particle": "lnp", "particles": "lnp",
    "candidate": "lnp", "candidates": "lnp",
    "library": "library", "libraries": "library",
    "member": "lipid", "members": "lipid", "compound": "lipid", "compounds": "lipid",
    "structure": "lipid", "structures": "lipid", "analog": "lipid", "analogs": "lipid",
    "analogue": "lipid", "analogues": "lipid",
}

# number words that appear as library sizes in older abstracts
NUMWORD = {
    "one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7,
    "eight": 8, "nine": 9, "ten": 10, "eleven": 11, "twelve": 12, "dozen": 12,
    "hundred": 100, "hundreds": None, "thousand": 1000, "thousands": None,
}

# reagent names whose digits are not library sizes
BLOCK = re.compile(
    r"(?:alc|sm|mc|dlin|c12|c14|c16|c18|ckk|lipid|peg|dmg|dspc|dope|dopc|dotap|"
    r"raw|hek|hela|a549|hepg2|dc2|b16|4t1|cho|il|tnf|cd|balb|c57bl|pbs|rpmi|dmem|"
    r"covid|sars|cov|mrna-|arcturus|ph|pka|nm|kda|kb|rna|dna|cas|sg)[\s-]*$",
    re.I,
)
NUM = re.compile(r"\b(\d{1,3}(?:,\d{3})+|\d{2,5})\b")


def abstracts(path):
    root = ET.parse(path).getroot()
    for a in root.findall(".//PubmedArticle"):
        pmid = a.findtext("./MedlineCitation/PMID") or ""
        art = a.find(".//Article")
        txt = " ".join(
            " ".join(e.itertext()) for e in art.findall("./Abstract/AbstractText")
        )
        title = " ".join(art.find("./ArticleTitle").itertext())
        yield pmid, title, re.sub(r"\s+", " ", txt).strip()


def hits(text):
    """Every number in `text` followed within 4 words by a counted-thing word."""
    out = []
    for m in NUM.finditer(text):
        pre = text[max(0, m.start() - 12):m.start()]
        if BLOCK.search(pre):
            continue
        n = int(m.group(1).replace(",", ""))
        if not (10 <= n <= 60000):
            continue
        tail = text[m.end():m.end() + 70]
        nxt = re.findall(r"[A-Za-z]+", tail)[:4]
        for i, w in enumerate(nxt):
            b = WORDS.get(w.lower())
            if b:
                ctx = text[max(0, m.start() - 60):m.end() + 60]
                out.append((n, b, i, " ".join(w for w in nxt[:i + 1]), ctx))
                break
    return out


def main():
    for pmid, title, ab in abstracts(sys.argv[1]):
        hs = hits(ab)
        print("=" * 108)
        print(f"{pmid}  {title[:95]}")
        if not ab:
            print("   (no abstract)")
            continue
        if not hs:
            # fall back: show any bare numbers so a human can judge
            bare = sorted({int(m.group(1).replace(",", "")) for m in NUM.finditer(ab)
                           if 10 <= int(m.group(1).replace(",", "")) <= 60000})
            print(f"   NO SIZE PHRASE. bare numbers present: {bare[:14]}")
        seen = set()
        for n, b, dist, phrase, ctx in sorted(hs, key=lambda x: -x[0]):
            if (n, b) in seen:
                continue
            seen.add((n, b))
            print(f"   {n:>6}  [{b:7}] {phrase}")
            print(f"           …{ctx.strip()}…")


if __name__ == "__main__":
    main()
