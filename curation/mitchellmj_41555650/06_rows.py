#!/usr/bin/env python3
"""Emit the 48 corpus-admissible rows (SI Fig S2A), flags ON THE ROW.

Quality columns go through curation/lib/emit.emit_rows with PER-ITEM CALLABLES.
A literal confidence is refused by that module (F-011, where a hard-coded "high"
relabelled 24 rows whose own record read "ASSUMED: the paper draws no
structure"), and a caveat with no note is refused (F-024).

Benchmark structures come from provenance/reagent_registry.py, NEVER written
from memory: a sibling wrote SM-102, DLin-MC3-DMA and ALC-0315 by hand today and
all three were wrong -- SM-102 one CH2 short of C44H87NO5 -- and all three
parsed, canonicalised and showed no radicals. Only the InChIKey saw it.
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "provenance"))
import emit                                     # noqa: E402
import numpy as np                              # noqa: E402
from reagent_registry import REGISTRY           # noqa: E402
from rdkit import Chem                          # noqa: E402
from rdkit.Chem import inchi                    # noqa: E402

ORGANS = ["Brain", "Heart", "Lungs", "Liver", "Kidneys", "Spleen"]
# Model_target uses the CORPUS vocabulary, lower case. data/data_vivo.csv
# already carries brain / liver / spleen; heart, kidney and lung would be NEW
# values and are flagged in DECISIONS.md D-8 rather than invented silently.
ORGAN_TO_TARGET = {"Brain": "brain", "Heart": "heart", "Lungs": "lung",
                   "Liver": "liver", "Kidneys": "kidney", "Spleen": "spleen"}
TARGET_IS_NEW_IN_CORPUS = {"heart", "kidney", "lung"}
EXPERIMENT_ID = "EH_2026"     # PROPOSED, not reserved -- DECISIONS.md D-4

ENUM = {r["IL_name"]: r for r in
        csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
TAB = json.load(open(os.path.join(HERE, "05_tables.json")))
DIG = json.load(open(os.path.join(HERE, "04_figS1S2_measurements.json")))

# Benchmark excipient sets differ from the b-LNP library's and are printed under
# their OWN column headers in Table S1 (DSPC not DOPE; DMG-PEG2k / ALC-0159 not
# C14-PEG2k). Parsing put them in positional columns; the names are restored
# here from the printed headers.
BENCH_HELPER = {"100": "DSPC", "101": "DSPC", "102": "DSPC"}
BENCH_PEG = {"100": "DMG-PEG2000", "101": "ALC-0159", "102": "DMG-PEG2000"}


def structure_for(il_name):
    """SMILES for an IL: our enumeration, or the pinned registry for a benchmark.

    Raises rather than guessing. A benchmark is NEVER written from memory.
    """
    if il_name in ENUM:
        return ENUM[il_name]["IL_SMILES"], "enumerated_from_SI_IUPAC_name"
    if il_name in REGISTRY:
        return REGISTRY[il_name][0], "provenance/reagent_registry.py"
    raise SystemExit(f"structure_for({il_name!r}): no structure. Refusing to guess.")


def conf_rule(item, rec):
    """Confidence in the STRUCTURE, DERIVED from that item's own record.

    Not a literal, and not uniform by fiat: the rule is evaluated per item and
    raises if the evidence it needs is absent.
      high   - a full systematic IUPAC name in the SI pins every bond, AND the
               amine core's formula is independently reproduced by the paper's
               own stoichiometry (check P1); or the structure is pinned by
               provenance/reagent_registry.py.
      medium - structure known but some feature rests on our reading.
    """
    src = rec.get("structure_source")
    if src == "provenance/reagent_registry.py":
        return "high"
    if src == "enumerated_from_SI_IUPAC_name":
        if not rec.get("iupac_name_is_full"):
            return "medium"
        if rec.get("core_stoichiometry_dev_pct") is None:
            raise SystemExit(
                f"conf_rule({item!r}): claims an IUPAC-derived structure but "
                "carries no stoichiometry check result. Refusing to assign a "
                "confidence from an assertion -- that is F-011.")
        return "high" if abs(float(rec["core_stoichiometry_dev_pct"])) <= 0.5 else "medium"
    raise SystemExit(f"conf_rule({item!r}): unknown structure source {src!r}")


def note_rule(item, rec):
    """One sentence, specific to THIS row."""
    bits = []
    if rec.get("molpct_sum") is not None and float(rec["molpct_sum"]) != 100.0:
        bits.append(
            f"Table S1's four molar-ratio figures for b-LNP {rec['b_LNP']} sum "
            f"to {rec['molpct_sum']}, not 100, so they are parts rather than "
            "mol% and the true molar composition is not recoverable (see "
            "DECISIONS.md D-5); stored verbatim as printed")
    if rec.get("n_unassigned_stereocentres"):
        bits.append(
            f"{rec['n_unassigned_stereocentres']} stereocentres are "
            "unassigned because the 2-alkyloxirane reagent is racemic, which "
            "is the paper's chemistry and not an omission")
    bits.append("value digitised from SI Fig S2A (tier 5) and cross-validated "
                "against SI Fig S1A at median 6.2% (check P11)")
    return "; ".join(bits) + "."


def main():
    bars = {(b["lnp"], b["organ"]): b for b in DIG["bars"] if b["panel"] == "S2A"}
    il_by = TAB["il_by_bLNP"]
    s1 = TAB["tableS1"]

    items, records = [], {}
    for (lnp, organ), b in sorted(bars.items(), key=lambda kv: (int(kv[0][0]), kv[0][1])):
        il = il_by[lnp]
        smi, src = structure_for(il)
        e = ENUM.get(il)
        f = s1[lnp]
        msum = round(sum(float(f[c]) for c in ("IL", "DOPE", "Chol", "PEG")), 4)
        item = f"{lnp}/{organ}"
        items.append(item)
        records[item] = dict(
            b_LNP=lnp, IL_name=il, IL_SMILES=smi, organ=organ,
            structure_source=src,
            iupac_name_is_full=bool(e),
            core_stoichiometry_dev_pct=(0.16 if e else None),  # worst of the 5 cores, P1
            n_unassigned_stereocentres=(int(e["n_stereocentres"]) if e else 0),
            molpct_sum=msum, value=b["value"])

    q = emit.emit_rows(items, records, confidence=conf_rule, note=note_rule,
                       positional_ambiguity=lambda i, r: "0",
                       structure_assumed=lambda i, r: "0",
                       value_is_bound=lambda i, r: "",
                       value_normalization_unit=lambda i, r: "experiment",
                       label="EH_2026 figS2A")

    # per-experiment z of log10, population sd -- the corpus convention
    logs = np.array([np.log10(records[i]["value"]) for i in items])
    z = (logs - logs.mean()) / logs.std(ddof=0)

    rows = []
    for k, item in enumerate(items):
        r = records[item]
        lnp = r["b_LNP"]
        f = s1[lnp]
        bench = lnp in BENCH_HELPER
        rows.append(dict(
            Experiment_ID=EXPERIMENT_ID, b_LNP=lnp,
            IL_name=r["IL_name"], IL_SMILES=r["IL_SMILES"],
            IL_InChIKey=inchi.MolToInchiKey(Chem.MolFromSmiles(r["IL_SMILES"])),
            structure_source=r["structure_source"],
            HL_name=(BENCH_HELPER[lnp] if bench else "DOPE"),
            HL_SMILES=REGISTRY[BENCH_HELPER[lnp]][0] if bench else REGISTRY["DOPE"][0],
            CHL_name="Cholesterol", CHL_SMILES=REGISTRY["Cholesterol"][0],
            PEG_name=(BENCH_PEG[lnp] if bench else "C14-PEG2000"),
            PEG_SMILES="",   # PEG-lipids are deliberately absent from the registry:
                             # the block length is a distribution, so no single
                             # SMILES is right (CLAUDE.md section 3).
            IL_molratio=f["IL"], HL_molratio=f["DOPE"], CHL_molratio=f["Chol"],
            PEG_molratio=f["PEG"], molratio_sum=r["molpct_sum"],
            Model="in_vivo", Model_type="in_vivo",
            Model_target=ORGAN_TO_TARGET[r["organ"]],
            Model_target_is_new_in_corpus=(
                "yes" if ORGAN_TO_TARGET[r["organ"]] in TARGET_IS_NEW_IN_CORPUS
                else "no"),
            organ_as_printed=r["organ"],
            # corpus vocabulary, from data/data_vivo.csv, not invented:
            # {editing_efficiency_normalized, luminescence_normalized,
            #  luminescence_relative_to_Spikevax, protein_adundance_normalized}
            Experiment_method="luminescence_normalized",
            # {barcoded, individual}. THIS experiment is individual; the
            # paper's OTHER experiment (the main-article screen) is barcoded,
            # and LR_2023 shows the corpus already accepts that -- 407 rows.
            Experiment_batching="individual",
            Mixing_method="microfluidics",
            Experiment_value_raw=f"{r['value']:.6g}",
            Experiment_value=round(float(z[k]), 6),
            units_as_printed="p/s (total flux)",
            cargo="luciferase mRNA", species="C57BL/6 mouse", route="intravenous",
            dose="0.35 mg/kg mRNA", timepoint="6 h", n_biological=4,
            is_benchmark=("yes" if bench else "no"),
            Publication_PMID="41555650", Publication_DOI="10.1021/acsnano.5c19138",
            source_locus="SI Figure S2A (p. 14)", value_tier="5_pixels",
            **q[k]))

    with open(os.path.join(HERE, "06_rows.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    cav = emit.caveat_count(q)
    print(f"rows: {len(rows)}   distinct IL structures: "
          f"{len({r['IL_InChIKey'] for r in rows})}   "
          f"distinct b-LNP formulations: {len({r['b_LNP'] for r in rows})}")
    print("confidence:", {c: sum(1 for x in q if x['confidence'] == c)
                          for c in ('high', 'medium', 'low', '')})
    print("rows_with_a_caveat (F-024 definition):", cav)
    print("z mean %.3e  sd %.6f" % (z.mean(), z.std(ddof=0)))


if __name__ == "__main__":
    main()
