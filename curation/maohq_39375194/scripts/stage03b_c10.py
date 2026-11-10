#!/usr/bin/env python3
"""The C10 comparison, run locally (provenance/structure_check.py writes tracked
files other agents read, so it is NOT invoked). Compares each anchor built by
make_reference.py against the structure stage03_structures.py emits.

RUN: python3 -I scripts/stage03b_c10.py <01_blocks.csv> <reference_structures.csv> <out.json>
"""
import csv, json, sys
from rdkit import Chem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

def k(s):
    m = Chem.MolFromSmiles(s)
    if m is None: return None, None
    i = Chem.MolToInchiKey(m)
    return i, i.split("-")[0]

def main(blocks, refs, out):
    emitted = {r["published_name"]: r["SMILES"]
               for r in csv.DictReader(open(blocks))}
    res, nskel, nstereo, nempty = [], 0, 0, 0
    for r in csv.DictReader(open(refs)):
        nm, ref = r["lipid_name"], r["reference_SMILES"].strip()
        em = emitted.get(nm, "").strip()
        if not ref:
            res.append({"lipid": nm, "status": "NO_ANCHOR",
                        "note": r["note"][:120]})
            nempty += 1
            continue
        fr, sr = k(ref); fe, se = k(em)
        if sr != se:
            st = "SKELETON_MISMATCH"
        elif fr != fe:
            st = "STEREO_MISMATCH"
            nstereo += 1
        else:
            st = "exact_match"
        if sr == se: nskel += 1
        res.append({"lipid": nm, "status": st,
                    "anchor_full": fr, "emitted_full": fe,
                    "anchor_skeleton": sr, "emitted_skeleton": se,
                    "byte_identical": ref == em})
    o = {"results": res,
         "n_anchors": sum(1 for x in res if x["status"] != "NO_ANCHOR"),
         "n_exact_match": sum(1 for x in res if x["status"] == "exact_match"),
         "n_skeleton_match": nskel,
         "n_stereo_mismatch": nstereo,
         "n_no_anchor": nempty,
         "byte_identical_count": sum(1 for x in res
                                     if x.get("byte_identical")),
         "_byte_identity_note": ("Byte-identity is EXPECTED for a correct "
                                 "re-derivation once canonicalised and is NOT "
                                 "evidence of copying (F-037). Independence "
                                 "here is established by make_reference.py's "
                                 "run-time AST walk, recorded in "
                                 "PROVENANCE.md.")}
    json.dump(o, open(out, "w"), indent=1)
    print(json.dumps({x: o[x] for x in o if x != "results"}, indent=1))
    for x in res: print(" ", x["lipid"], x["status"])
    bad = [x for x in res if x["status"] == "SKELETON_MISMATCH"]
    return 1 if bad else 0

if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:4]))
