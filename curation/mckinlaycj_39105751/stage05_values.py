#!/usr/bin/env python3
"""The five extractable values. Writes 05_values.csv. Flags go ON THE ROW.

Run from this folder:   python3 -I stage05_values.py

WHERE THESE VALUES COME FROM. Main-article Table 1, parsed from the PMC JATS
record as structured <table-wrap> markup. That is the publisher's own cell
geometry, so there is no row/column inference and nothing is digitised. It is
the ONLY per-member delivery value source in everything we hold; see
00_value_sources.json for the exhaustive inventory that establishes that.

THE SUPERSCRIPT. Every cell reads "2.46 +/- 0.55 x 10<sup>7</sup>". Flattening
the element tree with itertext() turns that into "2.46 +/- 0.55 x 107", which
parses cleanly as the number 107 and is wrong by fourteen orders of magnitude
while looking entirely plausible. <sup> is therefore handled explicitly and the
parser REFUSES a cell whose exponent it cannot locate, rather than falling
through to a plain-number reading. The same trap flattened Tesseract's
"1.2x10^5" into "109" in an earlier curation (CLAUDE.md section 5).

WHAT THE VALUE IS. Total photon flux, p/s, from IVIS at 6 h after a 0.125 mg/kg
IV dose of firefly-luciferase mRNA in female Balb/c mice, n = 3, mean +/- SD.
Formulation and dosing are quoted verbatim from the main-article Experimental
Methods in the same XML file.

WHAT THE VALUE IS NOT PINNED ON. Table 1's column header says only "luciferase
expression (p/s)" -- it does not name an organ. The surrounding text cites
"Table 1 and Figure 5D" together and Figure 5D is captioned "Total Fluc
bioluminescence in balb/c liver", so these are almost certainly liver values;
but "total flux" elsewhere in this paper means whole-animal. The two readings
differ by whatever fraction of signal is extrahepatic, which for these
liver-selective peptoids the paper puts above 97% at higher dose -- small, but
not zero and not measured here. Recorded on every row and in
08_contradictions.csv rather than silently resolved.
"""
import csv
import math
import os
import re
import sys
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, HERE)

import emit                        # noqa: E402

import blocks                      # noqa: E402

XML = os.path.join(HERE, "source", "PMC11342374.jats.xml")

# Verbatim from the main-article Experimental Methods, same XML file.
METADATA = {
    "model": "in_vivo",
    "species_strain": "Balb/c, female, 6-8 weeks",
    "route": "IV (tail vein), 100 uL injection volume",
    "dose": "0.125 mg/kg Fluc mRNA",
    "cargo": "firefly luciferase mRNA, N1-methylpseudouridine, CleanCap",
    "timepoint": "6 h post-injection (imaged 10 min after IP D-luciferin)",
    "readout": "total photon flux in the region of interest, IVIS Spectrum",
    "units_as_printed": "p/s",
    "n": 3,
    "statistic": "mean +/- SD",
    "formulation": "peptoid:DSPC:cholesterol:DMG-PEG2000 = 20:1.79:7.16:1.84 "
                   "by WEIGHT; microfluidic mixing against 10 mM sodium "
                   "citrate pH 5.0 at 1:3 v/v; dialysed into tris-sucrose "
                   "pH 7.4",
    "np_ratio": "NOT STATED. The paper gives a weight ratio of the four "
                "components and no N:P. Left blank rather than derived -- a "
                "peptoid's protonatable-nitrogen count at pH 5 is itself a "
                "prediction here (SI: Apd peptoids have predicted pKa ~5.5), "
                "so any N:P would be our arithmetic, not the paper's.",
}

ORGAN_AS_PRINTED = "not stated in Table 1"
ORGAN_INFERRED = ("liver -- the text cites 'Table 1 and Figure 5D' together "
                  "and Figure 5D is captioned 'Total Fluc bioluminescence in "
                  "balb/c liver'")

VALUE_RE = re.compile(r"([0-9.]+)\s*±\s*([0-9.]+)\s*×\s*10$")


def table1():
    """[(peptoid, label, mean, sd, printed)] from the JATS <table-wrap>."""
    root = ET.parse(XML).getroot()
    tw = None
    for t in root.iter("table-wrap"):
        lab = t.find("label")
        if lab is not None and "".join(lab.itertext()).strip() == "Table 1":
            tw = t
            break
    if tw is None:
        raise SystemExit("Table 1 not found in the JATS record.")

    out = []
    for tr in tw.iter("tr"):
        cells = [c for c in tr if c.tag in ("td", "th")]
        if len(cells) < 6 or cells[0].tag == "th":
            continue
        label = "".join(cells[0].itertext()).strip()
        m = re.match(r"^(\d+)", label)
        vcell = cells[5]
        # Explicit <sup> handling. The mantissa is the text BEFORE the <sup>;
        # the exponent is the <sup>'s own text. Anything else is refused.
        sup = vcell.find("sup")
        head = (vcell.text or "")
        if sup is None or not (sup.text or "").strip().lstrip("-").isdigit():
            if m and "".join(vcell.itertext()).strip():
                raise SystemExit(
                    f"row {label!r}: value cell "
                    f"{''.join(vcell.itertext())!r} has no parsable <sup> "
                    "exponent. Refusing to read it as a plain number -- "
                    "'2.46 +/- 0.55 x 10<sup>7</sup>' flattens to '...x 107'.")
            continue        # the 'model predicted' row has an empty value cell
        if not m:
            continue
        mm = VALUE_RE.search(head.strip())
        if not mm:
            raise SystemExit(f"row {label!r}: cannot parse mantissa {head!r}")
        exp = int(sup.text.strip())
        mean = float(mm.group(1)) * 10 ** exp
        sd = float(mm.group(2)) * 10 ** exp
        printed = f"{mm.group(1)} ± {mm.group(2)} × 10^{exp}"
        out.append((int(m.group(1)), label, mean, sd, printed,
                    [("".join(c.itertext()).strip()) for c in cells[1:5]]))
    return out


def load_enumerated():
    recs = {}
    with open(os.path.join(HERE, "03_enumerated.csv")) as fh:
        for r in csv.DictReader(fh):
            recs[int(r["peptoid"])] = r
    return recs


# --- the flag rules. Each is a function of the ITEM'S OWN record, because
#     emit.emit_rows refuses a literal: error class 29, found three times
#     (F-011, F-025, F-046), and in F-011 a hard-coded "high" relabelled 24
#     rows whose record read "ASSUMED: the paper draws no structure".

DRAWN_IN_FULL = {236, 266}      # SI Figure S7A draws these two complete

# A monomer whose structure rests on a NAME rather than on a bond-line drawing.
# Ole is the word "Oleyl" in Figure 1B with no structure beside it; Eae is a
# code that never appears in Figure 1B's pool at all and is resolved from the
# SI Materials list.
NAME_ONLY_MONOMERS = {"Ole", "Eae"}


def _confidence(pid, rec):
    mons = rec["sequence"].split("-")
    if any(blocks.MONOMERS[m]["sidechain"] is None
           and not blocks.MONOMERS[m].get("arms") for m in mons):
        return "low"
    if any(m in NAME_ONLY_MONOMERS for m in mons):
        return "medium"
    return "high"


def _structure_assumed(pid, rec):
    mons = rec["sequence"].split("-")
    return "1" if any(m in NAME_ONLY_MONOMERS for m in mons) else "0"


def _note(pid, rec):
    mons = rec["sequence"].split("-")
    bits = []
    name_only = sorted({m for m in mons if m in NAME_ONLY_MONOMERS})
    if name_only:
        bits.append(
            f"{'/'.join(name_only)} is defined in Figure 1B by a NAME, not a "
            "bond-line drawing, so its connectivity is the standard meaning "
            "of the word rather than a reading")
    if pid in DRAWN_IN_FULL:
        bits.append("structure drawn in full in SI Figure S7A and matched "
                    "monomer-for-monomer")
    if "Ehx" in mons:
        bits.append("Ehx (2-ethylhexyl) is a racemic centre the paper draws "
                    "without a wedge, stored flat")
    bits.append(f"value is Table 1's luciferase expression (p/s); the organ is "
                f"{ORGAN_AS_PRINTED} ({ORGAN_INFERRED})")
    return "; ".join(bits) + "."


def main():
    t1 = table1()
    enum = load_enumerated()

    ids = [p for p, *_ in t1]
    missing = [p for p in ids if p not in enum or not enum[p]["smiles"]]
    if missing:
        raise SystemExit(f"no enumerated structure for {missing}")
    records = {p: enum[p] for p in ids}

    # per-experiment z-score of log10, the corpus convention (CLAUDE.md s3),
    # population sd (ddof=0).
    logs = [math.log10(mean) for _, _, mean, _, _, _ in t1]
    mu = sum(logs) / len(logs)
    sd = math.sqrt(sum((x - mu) ** 2 for x in logs) / len(logs))

    flags = emit.emit_rows(
        ids, records,
        confidence=_confidence,
        note=_note,
        positional_ambiguity=lambda p, r: "0",
        structure_assumed=_structure_assumed,
        value_is_bound=lambda p, r: "",
        value_normalization_unit=lambda p, r: "experiment",
        within_model=None,
        label="mckinlaycj_39105751")

    qcols = emit.columns()
    fields = ["peptoid", "IL_name", "sequence", "n_residues", "IL_SMILES",
              "formula", "exact_mass", "source_location",
              "readout_as_printed", "units_as_printed", "value_printed",
              "value_mean_p_per_s", "value_sd_p_per_s", "log10_mean",
              "Experiment_value_z", "n_replicates", "statistic",
              "organ_as_printed", "organ_inferred", "model", "species_strain",
              "route", "dose", "cargo", "timepoint", "formulation",
              "np_ratio", "table1_lipid_monomers", "table1_lipid_carbons",
              "table1_branched_unsat", "table1_bu_identity"] + qcols

    out = os.path.join(HERE, "05_values.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for (pid, label, mean, sdv, printed, t1cols), fl in zip(t1, flags):
            rec = enum[pid]
            row = {
                "peptoid": pid,
                "IL_name": f"peptoid_{pid}",
                "sequence": rec["sequence"],
                "n_residues": rec["n_residues"],
                "IL_SMILES": rec["smiles"],
                "formula": rec["formula"],
                "exact_mass": rec["exact_mass"],
                "source_location": "main-article Table 1, column 'luciferase "
                                   "expression (p/s)', via PMC11342374 JATS "
                                   "<table-wrap id=tbl1>",
                "readout_as_printed": "luciferase expression (p/s)",
                "units_as_printed": METADATA["units_as_printed"],
                "value_printed": printed,
                "value_mean_p_per_s": f"{mean:.6g}",
                "value_sd_p_per_s": f"{sdv:.6g}",
                "log10_mean": round(math.log10(mean), 6),
                "Experiment_value_z": round((math.log10(mean) - mu) / sd, 6),
                "n_replicates": METADATA["n"],
                "statistic": METADATA["statistic"],
                "organ_as_printed": ORGAN_AS_PRINTED,
                "organ_inferred": ORGAN_INFERRED,
                "model": METADATA["model"],
                "species_strain": METADATA["species_strain"],
                "route": METADATA["route"],
                "dose": METADATA["dose"],
                "cargo": METADATA["cargo"],
                "timepoint": METADATA["timepoint"],
                "formulation": METADATA["formulation"],
                "np_ratio": METADATA["np_ratio"],
                "table1_lipid_monomers": t1cols[0],
                "table1_lipid_carbons": t1cols[1],
                "table1_branched_unsat": t1cols[2],
                "table1_bu_identity": t1cols[3],
            }
            row.update(fl)
            w.writerow(row)

    print(f"Table 1 rows parsed with an explicit <sup> exponent : {len(t1)}")
    for pid, label, mean, sdv, printed, _ in t1:
        print(f"  {pid:>4}  {printed:<22}  = {mean:.3g} p/s")
    print(f"\nz-score basis: n={len(t1)}, mean(log10)={mu:.4f}, "
          f"population sd={sd:.4f}")
    print(f"rows_with_a_caveat (F-024 definition): "
          f"{emit.caveat_count(flags)}/{len(flags)}")
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
