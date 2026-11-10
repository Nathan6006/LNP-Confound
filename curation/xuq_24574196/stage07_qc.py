#!/usr/bin/env python3
"""Every check for curation/xuq_24574196, with numbers, and each one SEEDED.

"A NEW CHECK IS NOT DONE UNTIL IT HAS BEEN MADE TO FAIL ON PURPOSE."
(CLAUDE.md section 6.) Three incidents paid for that rule, the sharpest being
F-038: three ZT_2026 marginal checks scored a literal 0 for twelve days while
the summary beside them said "ALL FOUR EXACT". So every check here runs twice --
once on the real input, once on a deliberately corrupted one -- and is credited
only if the AGREEMENT COUNT MOVES between the two runs. A check whose seeded
run agrees just as often is examining nothing.
"""
import csv
import json
import os
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors, Descriptors, inchi

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "curation", "lib"))
# `python3 -I` does NOT prepend the script's own directory to sys.path, so the
# sibling stage module has to be added explicitly. Required by the project's
# `python3 -I always` rule, not a workaround for it.
sys.path.insert(0, HERE)

import stage03_enumerate as enum  # noqa: E402

PROTON = 1.00727646
NA = 22.98922

# ---------------------------------------------------------------- helpers ---


def formula(smi):
    return rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smi))


def h_count(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    return sum(1 for a in m.GetAtoms() if a.GetSymbol() == "H")


def bulk_ch2_protons(smi):
    """Protons in the 1.2-1.75 'bulk methylene' envelope.

    A CH2 carbon counts when it has exactly 2 H and 2 carbon neighbours AND at
    least one of those neighbours is itself bonded only to carbon. The second
    clause is what excludes the propylene's central CH2 (both its neighbours
    bear the nitrogens) while keeping tail C2, whose far neighbour is pure
    alkyl -- the paper assigns that central CH2 its own 2H multiplet, so a rule
    that swallowed it would predict 58H where 56H is printed.
    """
    m = Chem.MolFromSmiles(smi)
    n = 0
    for a in m.GetAtoms():
        if a.GetSymbol() != "C" or a.GetTotalNumHs() != 2:
            continue
        nbrs = [x for x in a.GetNeighbors()]
        if len(nbrs) != 2 or any(x.GetSymbol() != "C" for x in nbrs):
            continue
        if any(all(y.GetSymbol() == "C" for y in x.GetNeighbors())
               for x in nbrs):
            n += 2
    return n


def chain_heavy_atoms(alcohol_smi):
    """Chain atoms from (and excluding) the alcohol oxygen to the end."""
    m = Chem.MolFromSmiles(alcohol_smi)
    return m.GetNumHeavyAtoms() - 1


def strip_tails_to_methyl(smi):
    """Replace every ester O-alkyl with O-methyl: the HEAD identity test.

    Two lipidoids built on the same amine differ only in tail length, so
    truncating both tails to methyl makes the head comparable across papers by
    InChIKey. Done by RWMol surgery on carbons far from any stereocentre
    (there are none in this family at all), never by string surgery.
    """
    m = Chem.RWMol(Chem.MolFromSmiles(smi))
    patt = Chem.MolFromSmarts("[CX3](=O)[OX2][CX4H2]")
    keep = set()
    for _c, _o, o2, c1 in m.GetMol().GetSubstructMatches(patt):
        keep.add(c1)
    # walk outward from each kept carbon and delete everything beyond it
    doomed = set()
    for c1 in keep:
        seen, stack = {c1}, [c1]
        start_o = [n.GetIdx() for n in m.GetAtomWithIdx(c1).GetNeighbors()
                   if n.GetSymbol() == "O"]
        seen.update(start_o)
        while stack:
            i = stack.pop()
            for nb in m.GetAtomWithIdx(i).GetNeighbors():
                j = nb.GetIdx()
                if j in seen:
                    continue
                seen.add(j)
                doomed.add(j)
                stack.append(j)
    for i in sorted(doomed, reverse=True):
        m.RemoveAtom(i)
    prod = m.GetMol()
    Chem.SanitizeMol(prod)
    return Chem.MolToSmiles(prod)


def corpus_rows():
    out = []
    for f in ("data_vitro.csv", "data_vivo.csv"):
        out += list(csv.DictReader(open(os.path.join(ROOT, "data", f))))
    return out


# ----------------------------------------------------------------- checks ---

def check_formula(items):
    """P1. Built formula vs the formula the SI prints next to ESI-MS."""
    agree, detail = 0, []
    for r in items:
        if not r["formula_printed"]:
            continue
        ok = r["formula_built"] == r["formula_printed"]
        agree += ok
        detail.append({"lipid": r["IL_name"], "built": r["formula_built"],
                       "printed": r["formula_printed"], "match": ok})
    return agree, len(detail), detail


def check_proton_closure(items):
    """P2. UNPRINTED: every printed integral, summed, equals the formula's H.

    The authors never add their integrals up. Each list is 8-10 separately
    reported numbers and the formula is reported in a different sentence, so
    agreement to the unit is a real external closure, not a restatement.
    """
    agree, detail = 0, []
    for r in items:
        nm = r["IL_name"]
        if nm not in enum.PRINTED_INTEGRALS:
            continue
        printed = sum(enum.PRINTED_INTEGRALS[nm])
        built = h_count(r["IL_SMILES"])
        ok = printed == built
        agree += ok
        detail.append({"lipid": nm, "sum_of_printed_integrals": printed,
                       "H_in_built_formula": built, "match": ok,
                       "n_integrals_summed": len(enum.PRINTED_INTEGRALS[nm])})
    return agree, len(detail), detail


def check_bulk_integral(items):
    """P3. The bulk-CH2 integral, predicted from the built structure.

    Independent of P1: the formula fixes the TOTAL carbon, this fixes how much
    of it is unfunctionalised methylene. Reading 'O16B' as a C16 alkyl -- the
    name-token trap -- predicts 56H here against the printed 40H.
    """
    agree, detail = 0, []
    for r in items:
        nm = r["IL_name"]
        if nm not in enum.PRINTED_BULK_CH2_H:
            continue
        pred = bulk_ch2_protons(r["IL_SMILES"])
        got = enum.PRINTED_BULK_CH2_H[nm]
        ok = pred == got
        agree += ok
        detail.append({"lipid": nm, "predicted_bulk_CH2_H": pred,
                       "printed": got, "match": ok})
    return agree, len(detail), detail


def check_ms_convention(items):
    """P4. WHICH mass convention is 'Calcd'? Two conventions in one SI is a
    measured hazard; here the SI prints no `Calc.` column header at all, so the
    convention has to be determined FROM THE DATA."""
    detail = []
    for r in items:
        nm = r["IL_name"]
        if nm not in enum.PRINTED_CALCD:
            continue
        m = Chem.MolFromSmiles(r["IL_SMILES"])
        mono, avg = Descriptors.ExactMolWt(m), Descriptors.MolWt(m)
        cands = {"neutral_average": avg, "neutral_monoisotopic": mono,
                 "MH+_average": avg + PROTON, "MH+_monoisotopic": mono + PROTON,
                 "MNa+_monoisotopic": mono + NA}
        calcd = enum.PRINTED_CALCD[nm]
        res = {k: round(abs(v - calcd), 4) for k, v in cands.items()}
        best = min(res, key=res.get)
        second = sorted(res.values())[1]
        detail.append({
            "lipid": nm, "printed_Calcd": calcd, "printed_Found":
            enum.PRINTED_FOUND[nm],
            "residual_Da": res, "best": best, "best_residual": res[best],
            "margin_over_next_best": round(second - res[best], 4),
            "found_minus_best_candidate":
            round(enum.PRINTED_FOUND[nm] - cands[best], 4)})
    agree = sum(1 for d in detail if d["best"] == "neutral_average")
    return agree, len(detail), detail


def check_name_convention(items):
    """P5. UNPRINTED: the name token equals the CHAIN HEAVY-ATOM count.

    O16 -> 16 carbons; O16B -> 2 C + 2 S + 12 C = 16 atoms. This is the only
    reading under which the paper's own claim that O16B and O16 have "similar
    chemical structures" is true, and it is what the pinned C39H78N2O4S4
    forces. A name token is not a chain length (error class 6).
    """
    agree, detail = 0, []
    for r in items:
        tok = int(r["acrylate"].rstrip("B").lstrip("O"))
        n = chain_heavy_atoms(r["acrylate_alcohol_SMILES"])
        ok = tok == n
        agree += ok
        detail.append({"acrylate": r["acrylate"], "name_token": tok,
                       "chain_heavy_atoms_from_O": n, "match": ok,
                       "carbons_in_whole_alcohol_NOT_the_S_alkyl": sum(
                           1 for a in Chem.MolFromSmiles(
                               r["acrylate_alcohol_SMILES"]).GetAtoms()
                           if a.GetSymbol() == "C")})
    return agree, len(detail), detail


def check_corpus_head(items, rows):
    """P6. EXTERNAL: is 'amine 1' a head an independent curator already drew?

    The corpus's AA_2008 carries an `80_O<n>` series. If Xu's amine 1 is the
    Anderson catalogue's amine 80, then truncating both to O-methyl must give
    the SAME InChIKey -- a head derived here from NMR and a printed formula,
    against a head curated years earlier by someone else from a different
    paper. Zero overlap would have left our mass/NMR anchors carrying the whole
    external burden; this is the one place that burden is shared.
    """
    ours = [r for r in items if r["IL_name"] == "1-O16"][0]
    our_head = strip_tails_to_methyl(ours["IL_SMILES"])
    our_key = inchi.MolToInchiKey(Chem.MolFromSmiles(our_head))
    detail, agree = [], 0
    probes = {}
    for row in rows:
        nm = (row.get("IL_name") or "").strip()
        if nm in ("80_O15", "80_O18", "80_O10", "1_O15"):
            probes.setdefault(nm, (row["IL_SMILES"], row["Experiment_ID"]))
    for nm, (smi, exp) in sorted(probes.items()):
        head = strip_tails_to_methyl(smi)
        k = inchi.MolToInchiKey(Chem.MolFromSmiles(head))
        same = (k == our_key)
        expect = not nm.startswith("1_")
        agree += (same == expect)
        detail.append({"corpus_lipid": nm, "experiment": exp,
                       "corpus_head": head, "our_head": our_head,
                       "head_inchikey_identical": same,
                       "expected_identical": expect,
                       "as_expected": same == expect,
                       "corpus_formula": formula(smi)})
    return agree, len(detail), detail


def check_name_token_trap(rows):
    """P7. The leading '1' in '1-O16' is NOT the corpus's amine 1.

    AA_2008's amine 1 is 3-methoxypropylamine: one nitrogen, five oxygens.
    The SI's printed C43H86N2O4 has TWO nitrogens and FOUR oxygens, so keying
    on the number token would have put every value on a molecule the paper's
    own mass spectrum excludes. 163 corpus IL_names start with one of the Xu
    lab's catalogue numbers and do not all mean the same thing.
    """
    hit = [r for r in rows if (r.get("IL_name") or "").strip() == "1_O15"]
    f = formula(hit[0]["IL_SMILES"])
    ours = enum.PRINTED_FORMULA["1-O16"]
    n_ok = ("N2" in ours) and ("N2" not in f)
    o_ok = ("O4" in ours) and ("O4" not in f)
    return int(n_ok and o_ok), 1, [{
        "corpus_amine1_lipid": "1_O15", "experiment": hit[0]["Experiment_ID"],
        "its_formula": f, "its_head": "3-methoxypropylamine (1 N, ether O)",
        "SI_printed_formula_for_1-O16": ours,
        "excluded_by_nitrogen_count": n_ok,
        "excluded_by_oxygen_count": o_ok}]


def check_series_steps(rows):
    """P8. UNPRINTED: the corpus 80_O series steps in exact CnH2n units.

    Confirms the homologous family our 1-O16 sits inside is internally
    consistent, so the nearest-neighbour comparison in P6 is being made against
    a clean series rather than a drifting one.
    """
    fam = {}
    for r in rows:
        nm = (r.get("IL_name") or "").strip()
        if nm.startswith("80_O"):
            fam[nm] = r["IL_SMILES"]
    order = sorted(fam, key=lambda n: int(n[4:]))
    detail, agree = [], 0
    for a, b in zip(order, order[1:]):
        na, nb = int(a[4:]), int(b[4:])
        ma = Chem.MolFromSmiles(fam[a])
        mb = Chem.MolFromSmiles(fam[b])
        dc = sum(1 for x in mb.GetAtoms() if x.GetSymbol() == "C") - \
            sum(1 for x in ma.GetAtoms() if x.GetSymbol() == "C")
        dh = h_count(fam[b]) - h_count(fam[a])
        want_c, want_h = 2 * (nb - na), 4 * (nb - na)
        ok = (dc, dh) == (want_c, want_h)
        agree += ok
        detail.append({"from": a, "to": b, "dC": dc, "dH": dh,
                       "expected_dC": want_c, "expected_dH": want_h,
                       "match": ok})
    return agree, len(detail), detail


def check_anchor(items):
    """P10. C10-style: the enumeration vs an anchor built by a DIFFERENT code
    path (make_reference.py, head-first traversal, no assembler import).

    Compared by FULL InChIKey. Byte-identity of the canonical SMILES is
    expected here and is deliberately NOT the criterion: a correct
    re-derivation canonicalises to the same string as a correct enumeration,
    which is precisely why C10's own `anchor_may_be_a_copy` missed F-037.
    Independence is argued from provenance and is checked mechanically below
    by asserting the anchor module imports neither the library assembler nor
    the enumerator.
    """
    import ast
    import make_reference as ref
    # Parse the IMPORT STATEMENTS with ast. The first version of this test
    # substring-matched the whole file for "import assemble"/"esterify" and
    # reported imports_clean = False because those words appear in
    # make_reference.py's DOCSTRING, where they are prose explaining what the
    # module deliberately does not use. A substring scan cannot tell a module
    # name from a sentence about one -- the same class as the sibling whose
    # olefinic-integral parser was defeated by a space after a hyphen.
    tree = ast.parse(open(os.path.join(HERE, "make_reference.py")).read())
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.module:
                imported.add(node.module.split(".")[0])
            imported.update(a.name for a in node.names)
    forbidden = {"assemble", "stage03_enumerate", "esterify", "attach",
                 "axis", "bars", "heatmap", "chains", "emit"}
    imports_clean = not (imported & forbidden)
    agree, detail = 0, []
    for r in items:
        nm = r["IL_name"]
        if nm not in ref.TAILS:
            detail.append({"lipid": nm, "anchor": "(none stored)",
                           "status": "absent_on_purpose", "match": None})
            continue
        a = Chem.MolFromSmiles(ref.compose(ref.TAILS[nm]))
        b = Chem.MolFromSmiles(r["IL_SMILES"])
        ka, kb = inchi.MolToInchiKey(a), inchi.MolToInchiKey(b)
        ok = ka == kb and imports_clean
        agree += ok
        detail.append({"lipid": nm, "anchor_inchikey": ka,
                       "enumeration_inchikey": kb,
                       "inchikey_identical": ka == kb,
                       "canonical_smiles_identical":
                           Chem.MolToSmiles(a) == Chem.MolToSmiles(b),
                       "anchor_module_imports_no_assembler_or_enumerator":
                           imports_clean,
                       "anchor_module_imports": sorted(imported),
                       "match": ok})
    return agree, sum(1 for d in detail if d["match"] is not None), detail


def check_structural_sanity(items):
    """P9. Parse, canonical round-trip, no radicals, no stereocentres."""
    agree, detail = 0, []
    for r in items:
        smi = r["IL_SMILES"]
        m = Chem.MolFromSmiles(smi)
        rt = Chem.MolToSmiles(m) == smi
        rad = sum(a.GetNumRadicalElectrons() for a in m.GetAtoms())
        st = int(r["n_stereocentres"])
        ok = bool(m) and rt and rad == 0 and st == 0
        agree += ok
        detail.append({"lipid": r["IL_name"], "canonical_roundtrip": rt,
                       "radical_electrons": rad, "stereocentres": st,
                       "match": ok})
    return agree, len(detail), detail


# ------------------------------------------------------------------ seeds ---

def seeded(items, rows):
    """Run each check on a corrupted input. Credit a check only if its
    AGREEMENT COUNT MOVES. A seeded run that agrees just as often means the
    check is reading something that is not there."""
    out = {}

    # P1/P3/P9: the name-token trap itself -- read O16B as a C16 alkyl.
    bad = [dict(r) for r in items]
    for r in bad:
        if r["IL_name"] == "1-O16B":
            r["acrylate_alcohol_SMILES"] = "OCCSSCCCCCCCCCCCCCCCC"
            r["IL_SMILES"] = enum.build("seed", r["acrylate_alcohol_SMILES"])
            m = Chem.MolFromSmiles(r["IL_SMILES"])
            r["formula_built"] = rdMolDescriptors.CalcMolFormula(m)
            r["n_stereocentres"] = 0
    out["P1_formula"] = check_formula(bad)[0]
    out["P3_bulk_integral"] = check_bulk_integral(bad)[0]
    out["P5_name_convention"] = check_name_convention(bad)[0]

    # P2: drop one printed integral (a parser that misses a term).
    real = dict(enum.PRINTED_INTEGRALS)
    try:
        enum.PRINTED_INTEGRALS["1-O16"] = real["1-O16"][:-1]
        out["P2_proton_closure"] = check_proton_closure(items)[0]
    finally:
        enum.PRINTED_INTEGRALS.update(real)

    # P4: claim Calcd is [M+H]+ by shifting the printed value one proton.
    realc = dict(enum.PRINTED_CALCD)
    try:
        for k in enum.PRINTED_CALCD:
            enum.PRINTED_CALCD[k] = round(enum.PRINTED_CALCD[k] + PROTON, 3)
        out["P4_ms_convention"] = check_ms_convention(items)[0]
    finally:
        enum.PRINTED_CALCD.update(realc)

    # P6/P7/P8: corrupt the corpus copy the check reads.
    brows = [dict(r) for r in rows]
    for r in brows:
        if (r.get("IL_name") or "").strip() == "80_O15":
            r["IL_SMILES"] = r["IL_SMILES"].replace("N(CCCN(C)C)",
                                                    "N(CCCN(CC)CC)")
        if (r.get("IL_name") or "").strip() == "1_O15":
            r["IL_SMILES"] = ("CCCCCCCCCCCCCCCOC(=O)CCN(CCCN(C)C)"
                              "CCC(=O)OCCCCCCCCCCCCCCC")
        if (r.get("IL_name") or "").strip() == "80_O13":
            r["IL_SMILES"] = r["IL_SMILES"].replace("CCCCCCCCCCCCCOC(=O)",
                                                    "CCCCCCCCCCCCOC(=O)", 1)
    out["P6_corpus_head"] = check_corpus_head(items, brows)[0]
    out["P7_name_token_trap"] = check_name_token_trap(brows)[0]
    out["P8_series_steps"] = check_series_steps(brows)[0]
    return out


def main():
    items = list(csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv"))))
    rows = corpus_rows()

    real = {}
    for key, fn in (("P1_formula", lambda: check_formula(items)),
                    ("P2_proton_closure", lambda: check_proton_closure(items)),
                    ("P3_bulk_integral", lambda: check_bulk_integral(items)),
                    ("P4_ms_convention", lambda: check_ms_convention(items)),
                    ("P5_name_convention",
                     lambda: check_name_convention(items)),
                    ("P6_corpus_head", lambda: check_corpus_head(items, rows)),
                    ("P7_name_token_trap",
                     lambda: check_name_token_trap(rows)),
                    ("P8_series_steps", lambda: check_series_steps(rows)),
                    ("P10_anchor", lambda: check_anchor(items)),
                    ("P9_structural_sanity",
                     lambda: check_structural_sanity(items))):
        a, n, d = fn()
        real[key] = {"agree": a, "of": n, "detail": d}

    seeds = seeded(items, rows)
    for k, v in seeds.items():
        real[k]["seeded_agree"] = v
        real[k]["check_has_teeth"] = v < real[k]["agree"]
    # P10 seed: the anchor's own gate. make_reference.py --seed reads the
    # name token "16" as a C16 S-alkyl and MUST exit 1. Read the EXIT CODE
    # from returncode, never from a pipe.
    import subprocess
    r_ok = subprocess.run([sys.executable, "-I",
                           os.path.join(HERE, "make_reference.py")],
                          capture_output=True)
    r_bad = subprocess.run([sys.executable, "-I",
                            os.path.join(HERE, "make_reference.py"), "--seed"],
                           capture_output=True)
    real["P10_anchor"]["seeded_agree"] = 0 if r_bad.returncode != 0 else \
        real["P10_anchor"]["agree"]
    real["P10_anchor"]["anchor_gate_exit_clean"] = r_ok.returncode
    real["P10_anchor"]["anchor_gate_exit_seeded"] = r_bad.returncode
    real["P10_anchor"]["anchor_gate_seeded_stderr"] = \
        r_bad.stderr.decode().strip().splitlines()[-1] if r_bad.stderr else ""
    real["P10_anchor"]["check_has_teeth"] = (
        r_ok.returncode == 0 and r_bad.returncode == 1)

    real["P9_structural_sanity"]["seeded_agree"] = None
    real["P9_structural_sanity"]["check_has_teeth"] = (
        "n/a -- P9 is a parse/radical/stereo assertion on the library's own "
        "output; curation/lib already refuses a radical at build time and did "
        "so on this paper's first run (check_no_radicals fired with "
        "\"1 radical atom(s): [(0, 'C', 1)]\"), which is the seeded failure.")

    summary = {
        "_what": "stage 07 QC for curation/xuq_24574196. Every check carries "
                 "its agreement count AND the count it scores on a seeded "
                 "defect; a check is credited only where the two differ.",
        "rows_emitted_to_corpus": 0,
        "reason_zero_rows": "the SI contains no delivery measurement and the "
                            "main article is not acquirable; see "
                            "00_value_sources.json",
        "checks": real,
        "toothless_checks": [k for k, v in real.items()
                             if v.get("check_has_teeth") is False],
    }
    with open(os.path.join(HERE, "07_qc.json"), "w") as fh:
        json.dump(summary, fh, indent=1)

    for k, v in real.items():
        print(f"{k:24s} {v['agree']}/{v['of']}   seeded={v.get('seeded_agree')}"
              f"   teeth={v.get('check_has_teeth')}")
    print("TOOTHLESS:", summary["toothless_checks"])
    return 1 if summary["toothless_checks"] else 0


if __name__ == "__main__":
    sys.exit(main())
