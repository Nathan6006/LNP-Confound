#!/usr/bin/env python3
"""Stage 01 -- extract every native SI table and write the formulation enumeration.

    python3 -I 01_si_tables.py

Tier 4 (pdfplumber native tables). No pixels anywhere in this file.

Writes:
    04_si_tables.json   every extracted table, verbatim, plus the arithmetic checks
    01_blocks.csv       the components the two libraries are built from
    03_enumerated.csv   all 160 formulated LNPs (120 library 1 + 40 library 2)

Reads ONLY source/ja8b08976_si_001.pdf.
"""
import csv, json, os, sys
import pdfplumber

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                    "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))
import reagents                                                # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PDF = os.path.join(HERE, "source", "ja8b08976_si_001.pdf")

# SI page numbers (1-based in the PDF). The SI's own page labels are S1..S15,
# offset by one because the cover page is unlabelled -- cite by CONTENT, which
# is what the assertions below pin, not by the contents page (see
# 08_contradictions.csv, the contents page is off against its own captions).
P_NORM = 5      # Figure S1E, the worked normalisation example
P_BARCODE = 6   # Figure S1F, the 156 barcodes
P_LIB1 = 9      # Figure S2A, Library 1 Composition
P_LIB2 = 12     # Figure S2F, Library 2 Composition


def tables(pdf, page_no):
    pg = pdf.pages[page_no - 1]
    return [t.extract() for t in pg.find_tables()]


def num(x):
    if x is None or x == "":
        raise ValueError("empty numeric cell -- refusing to guess")
    return float(x)


def parse_lib1(rows):
    """61x13: two side-by-side 6-column blocks, LNP 1-60 and 61-120."""
    hdr = rows[0]
    assert hdr[:6] == ["LNP #", "7C1 Mole %", "Chol %", "PEG %", "DSPC %", "PEG Type"], hdr
    assert hdr[7:] == ["LNP #", "7C1 Mole %", "Chol %", "PEG %", "DSPC %", "PEG Type"], hdr
    out = {}
    for r in rows[1:]:
        for off in (0, 7):
            blk = r[off:off + 6]
            n = int(blk[0])
            out[n] = {"lnp": n, "il_molpct": num(blk[1]), "chol_molpct": num(blk[2]),
                      "peg_molpct": num(blk[3]), "dspc_molpct": num(blk[4]),
                      "peg_type": blk[5]}
    assert sorted(out) == list(range(1, 121)), "library 1 is not 1..120"
    return out


def parse_lib2(rows):
    """41x6. The PEG-type column header is BLANK in the PDF; pinned by content."""
    hdr = rows[0]
    assert hdr[:5] == ["LNP #", "7C1 Mole %", "Cholesterol Mole %",
                       "PEG Mole %", "DSPC Mole %"], hdr
    assert hdr[5] in ("", None), "column 6 was expected to be the unlabelled PEG-type column"
    out = {}
    for r in rows[1:]:
        n = int(r[0])
        out[n] = {"lnp": n, "il_molpct": num(r[1]), "chol_molpct": num(r[2]),
                  "peg_molpct": num(r[3]), "dspc_molpct": num(r[4]),
                  "peg_type": r[5]}
    assert sorted(out) == list(range(1, 41)), "library 2 is not 1..40"
    return out


def ratio_key(f):
    return (f["il_molpct"], f["chol_molpct"], f["peg_molpct"], f["dspc_molpct"])


def factorial_check(lib, label, expect_pegs, expect_ratios):
    """PROVE the factorial: every PEG variant must carry the SAME ratio set, once."""
    by_peg = {}
    for f in lib.values():
        by_peg.setdefault(f["peg_type"], []).append(ratio_key(f))
    sets = {p: sorted(v) for p, v in by_peg.items()}
    first = sorted(sets[sorted(sets)[0]])
    identical = all(sorted(v) == first for v in sets.values())
    dupes = {p: len(v) - len(set(v)) for p, v in by_peg.items()}
    return {
        "label": label,
        "n_peg_variants": len(by_peg),
        "n_ratios_per_variant": {p: len(v) for p, v in by_peg.items()},
        "distinct_ratios_per_variant": {p: len(set(v)) for p, v in by_peg.items()},
        "every_variant_carries_the_same_ratio_set": identical,
        "duplicate_ratios_within_a_variant": dupes,
        "product": len(by_peg) * len(first),
        "n_rows": len(lib),
        "product_equals_rows": len(by_peg) * len(first) == len(lib),
        "claimed": f"{expect_pegs} PEG types x {expect_ratios} formulation ratios = {expect_pegs * expect_ratios}",
        "claim_holds": (len(by_peg) == expect_pegs and len(first) == expect_ratios
                        and len(by_peg) * len(first) == len(lib)),
    }


def molpct_check(lib, label):
    bad = []
    for n, f in sorted(lib.items()):
        s = f["il_molpct"] + f["chol_molpct"] + f["peg_molpct"] + f["dspc_molpct"]
        if abs(s - 100.0) > 1e-9:
            bad.append({"lnp": n, "sum": round(s, 4), "row": f})
    return {"label": label, "n": len(lib), "n_sum_to_100": len(lib) - len(bad),
            "n_violations": len(bad), "violations": bad}


def main():
    with pdfplumber.open(PDF) as pdf:
        t_norm = tables(pdf, P_NORM)
        t_bc = tables(pdf, P_BARCODE)
        t1 = tables(pdf, P_LIB1)
        t2 = tables(pdf, P_LIB2)

    lib1 = parse_lib1(t1[0])
    lib2 = parse_lib2(t2[0])

    # --- the 156 barcodes -------------------------------------------------
    bc = {}
    hdr = t_bc[0][0]
    assert hdr == ["BC #", "Seq", "BC #", "Seq", "BC #", "Seq", "BC #", "Seq"], hdr
    for r in t_bc[0][1:]:
        for off in (0, 2, 4, 6):
            n, s = r[off], r[off + 1]
            if n in (None, "") and s in (None, ""):
                continue
            bc[int(n)] = s
    barcode = {
        "n_cells_in_grid": 4 * (len(t_bc[0]) - 1),
        "n_barcodes": len(bc),
        "numbering_is_1_to_n": sorted(bc) == list(range(1, len(bc) + 1)),
        "all_8_nt": all(len(v) == 8 for v in bc.values()),
        "n_distinct_sequences": len(set(bc.values())),
        "claimed_in_text": 156,
        "matches_claim": len(bc) == 156 and len(set(bc.values())) == 156,
    }

    out = {
        "_what": "Every native (tier-4) table in the SI, verbatim, plus the arithmetic "
                 "checks run on them. Produced by 01_si_tables.py.",
        "source_pdf_sha256": "c09c8d95ba38590beb15d3d9b0a68e7d164a81fea94b31d4684cae6bbd03296a",
        "library_1": lib1,
        "library_2": lib2,
        "barcodes": barcode,
        "normalisation_tables_fig_S1E": t_norm,
        "checks": {
            "factorial_library_1": factorial_check(lib1, "library 1 (Fig S2A)", 5, 24),
            "factorial_library_2": factorial_check(lib2, "library 2 (Fig S2F)", 2, 20),
            "molpct_library_1": molpct_check(lib1, "library 1"),
            "molpct_library_2": molpct_check(lib2, "library 2"),
            "total_formulated": {
                "library_1": len(lib1), "library_2": len(lib2),
                "total": len(lib1) + len(lib2),
                "drawn_in_Fig_1D": 160,
                "matches": len(lib1) + len(lib2) == 160,
                "note": "Fig 1D draws '2x; 160 LNPs tested in vivo' and never computes it; "
                        "120 + 40 reproduces it exactly.",
            },
        },
    }
    with open(os.path.join(HERE, "04_si_tables.json"), "w") as f:
        json.dump(out, f, indent=1)

    # --- 01_blocks.csv ----------------------------------------------------
    blocks = [
        ["component", "role", "name_in_paper", "where_drawn", "structure_known",
         "smiles_or_why_not", "note"],
        ["ionizable", "ionizable lipomer", "7C1", "Fig 3A, Fig 5A (one arm only); "
         "Fig S1B (the synthesis)", "NO",
         "",
         "Branched PEI Mn~600 (1 molar eq) + a 1,2-epoxide with R = C13 (14 molar eq), "
         "90C/72h/EtOH. A polydisperse, statistically substituted polymer: no single "
         "SMILES is correct. Fig 3A/5A draw ONE representative arm with the backbone "
         "nitrogen replaced by a diethylamino group."],
        ["helper", "phospholipid", "DSPC", "Fig 5A", "YES",
         reagents.registered("DSPC"),
         "Library 2 only, at 0 or 10 mol%. Drawn in full in Fig 5A. SMILES from "
         "provenance/reagent_registry.py, never from memory; the independently "
         "assembled anchor in reference_structures.csv agrees by full InChIKey."],
        ["sterol", "cholesterol", "Cholesterol", "Fig 3A, Fig 5A", "YES",
         reagents.registered("Cholesterol"),
         "3-beta OH, drawn with wedges in both panels. SMILES from "
         "provenance/reagent_registry.py; the hand re-traversal written for the "
         "anchor gave a DIASTEREOMER (see reference_structures.csv), which is "
         "why no excipient SMILES in this folder is typed from memory."],
        ["peg", "PEG-lipid", "C14PEG2000", "Fig 3B (generic skeleton)", "ASSUMED",
         "", "Drawn as a diacyl-glycero-phosphoethanolamine-PEG with (OCH2CH2)45OCH3. "
         "C14 => DMPE-PEG2000 by that skeleton, but the paper never names a catalogue "
         "product and never states the PEG dispersity; not emitted."],
        ["peg", "PEG-lipid", "C16PEG2000", "Fig 3B, Fig 5A", "ASSUMED", "",
         "C16 => DPPE-PEG2000 by the drawn skeleton. Not named by the paper."],
        ["peg", "PEG-lipid", "C18PEG2000", "Fig 3B, Fig 5A", "ASSUMED", "",
         "C18 => DSPE-PEG2000 by the drawn skeleton. Not named by the paper."],
        ["peg", "PEG-lipid", "C14PEG3000", "Fig 3B", "ASSUMED", "", "Library 1 only."],
        ["peg", "PEG-lipid", "C18PEG3000", "Fig 3B", "ASSUMED", "", "Library 1 only."],
        ["cargo", "siRNA", "siICAM-2 / siLuciferase", "Fig S1C", "sequence given", "",
         "2'-O-methyl modified, dTdT overhangs, phosphorothioate."],
        ["cargo", "DNA barcode", "156 x 91-nt barcodes", "Fig S1F", "sequence given", "",
         "siRNA : DNA barcode mass ratio 10:1."],
    ]
    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as f:
        csv.writer(f).writerows(blocks)

    # --- 03_enumerated.csv ------------------------------------------------
    rows = [["library", "lnp_number", "il_name", "il_smiles", "il_molpct",
             "chol_molpct", "peg_molpct", "dspc_molpct", "peg_type",
             "molpct_sum", "pooled_and_measured", "individually_measured",
             "individual_label", "structure_resolved", "note"]]
    # The only formulations the paper measures ONE AT A TIME, from Fig S3A / Fig 6A.
    named = {
        ("1", 25): "Screen 1 Winner (Fig S3A)",
        ("2", 6): "Screen 2 LNP1 (Fig S3A)",
        ("2", 18): "BM1 -- Screen 2 winner (Fig S3A, Fig 6A)",
        ("2", 19): "Screen 2 LNP3 (Fig S3A)",
    }
    for lab, lib in (("1", lib1), ("2", lib2)):
        for n, f in sorted(lib.items()):
            s = f["il_molpct"] + f["chol_molpct"] + f["peg_molpct"] + f["dspc_molpct"]
            ind = named.get((lab, n), "")
            rows.append([lab, n, "7C1", "", f["il_molpct"], f["chol_molpct"],
                         f["peg_molpct"], f["dspc_molpct"], f["peg_type"],
                         round(s, 4), "unknown -- see note", "yes" if ind else "no",
                         ind, "NO",
                         "IL structure unresolved (polydisperse PEI600/C15-epoxide "
                         "polymer); pooled membership not published."])
    # The 'original' 80:20 7C1 formulation is in neither library table.
    rows.append(["0 (reference)", "", "7C1", "", 80.0, 0.0, 20.0, 0.0, "C14PEG2000",
                 100.0, "no", "yes", "Original 80:20 7C1 (Fig 6A, Fig S3A)", "NO",
                 "The comparator formulation. Stated in the Results text, Fig 6A and "
                 "Fig S3A; absent from both library composition tables."])
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as f:
        csv.writer(f).writerows(rows)

    print(json.dumps(out["checks"], indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
