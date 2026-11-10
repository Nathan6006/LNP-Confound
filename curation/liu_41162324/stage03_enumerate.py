#!/usr/bin/env python3
"""Stage 03 -- enumerate the 74 FILs and check every one against the SI's HRMS.

THE LIBRARY'S LOGIC, read out of the SI's own synthesis text (section 3.2), not
inferred from the names. Every FIL is one tertiary amine:

    N( hydroxyalkyl )( lipid-tail arm )( fluorinated-ester arm )

  * the hydroxyalkyl and lipid-tail arms arrive together, as a SECONDARY AMINE
    intermediate that section 3.1 synthesises and section 3.2 names per
    compound: `Heptadecan-9-yl 8-((2-hydroxyethyl)amino)octanoate` for the
    SS series (3.2.1.1), `6-((4-Hydroxybutyl)amino)hexyl 2-hexyldecanoate` for
    AL and `... 2-heptylnonanoate` for AS (3.2.3.3-3.2.3.10);
  * the fluorinated arm is a 6-BROMOHEXANOATE, so the nitrogen alkylates at the
    carbon that bore the bromine:
        N-(CH2)5-C(=O)O-[linker]-C(=O)O-[fluoroalkyl]

NAME TOKENS ARE NOT CHAIN LENGTHS, and here they are nearly -- but not quite.
`Cn` is the number of main-chain atoms between the two ester carbonyls, `Bn` is
that same C4 backbone carrying an n-carbon branch at the carbon bearing the
N-side ester oxygen, `O` swaps a methylene for an ether oxygen and `P` for a
para-phenylene reached through a benzylic CH2. `chains.check_name_token`
compares the token to the structure we built and REPORTS; it never corrects.

WHY THE ARMS ARE WRITTEN OUT PER LINKER rather than concatenated from a
`linker` string. String concatenation of SMILES fused piperazine to cholesterol
through a ring-closure digit collision (CLAUDE.md section 6). The `P` linker is
the only fragment here carrying a ring, so each arm is one explicit template
with a single `{TAIL}` slot, the tails carry no ring digit at all, and
`_check_no_ring_digit_collision` asserts that before anything is parsed.

THE EXTERNAL CHECK. The SI prints, for 42 of the 74, a molecular formula and a
calculated [M+H]+. Both are compared against the structure we built:
  - formula equality is an atom-count test the drawing cannot fake;
  - [M+H]+ recomputed from OUR structure against the SI's own printed number is
    a second, independent comparison, and it catches the case where the SI's
    formula and the SI's m/z disagree with EACH OTHER -- which they do, for
    every member of the P series and for five others.
Mass still cannot see connectivity: a formula match says nothing about which
nitrogen or which oxygen. That is C10's job.

    python3 -I stage03_enumerate.py
"""
import csv
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "lib"))
import assemble                                              # noqa: E402
import chains                                                # noqa: E402
from errors import StructureError                            # noqa: E402

from rdkit import Chem                                       # noqa: E402
from rdkit.Chem.Descriptors import ExactMolWt                # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
PROTON = 1.00727646                     # mass of H+, not of an H atom

# Monoisotopic masses, for recomputing [M+H]+ from the SI's OWN PRINTED
# FORMULA. That comparison is entirely internal to the SI and does not involve
# our structure at all, which is exactly why it is worth making: when the
# printed formula and the printed m/z disagree with each other, the HRMS line
# is not evidence about anything and must not be used to adjudicate a
# structure. When they agree with each other and disagree with us, the paper is
# making one coherent claim and the disagreement is ours to explain.
MONO = {"C": 12.0, "H": 1.00782503, "F": 18.99840322,
        "N": 14.0030740, "O": 15.9949146}

# ---------------------------------------------------------------- blocks ----
# Each core is the SECONDARY amine the SI names, with its nitrogen atom-mapped
# for assemble.attach. The map is the only edit; the rest is the SI's own name
# transcribed to SMILES.
CORES = {
 "SS": {"smiles": "CCCCCCCCC(CCCCCCCC)OC(=O)CCCCCCC[NH:1]CCO",
        "si_name": "heptadecan-9-yl 8-((2-hydroxyethyl)amino)octanoate",
        "si_locus": "SI 3.2.1.1 (SI PDF p.24), reagent sentence",
        "head": "2-hydroxyethyl", "tail": "heptadecan-9-yl octanoate"},
 "SP": {"smiles": "CCCCCCCC(CCCCCCC)COC(=O)CCCCCCC[NH:1]CCO",
        "si_name": "2-heptylnonyl 8-((2-hydroxyethyl)amino)octanoate",
        "si_locus": "SI 3.1.47 heading (SI PDF p.22)",
        "head": "2-hydroxyethyl", "tail": "2-heptylnonyl octanoate"},
 "AL": {"smiles": "CCCCCCCCC(CCCCCC)C(=O)OCCCCCC[NH:1]CCCCO",
        "si_name": "6-((4-hydroxybutyl)amino)hexyl 2-hexyldecanoate",
        "si_locus": "SI 3.2.3.3 (SI PDF p.36), reagent sentence",
        "head": "4-hydroxybutyl", "tail": "hexyl 2-hexyldecanoate"},
 "AS": {"smiles": "CCCCCCCC(CCCCCCC)C(=O)OCCCCCC[NH:1]CCCCO",
        "si_name": "6-((4-hydroxybutyl)amino)hexyl 2-heptylnonanoate",
        "si_locus": "SI 3.2.3.4 (SI PDF p.36), reagent sentence",
        "head": "4-hydroxybutyl", "tail": "hexyl 2-heptylnonanoate"},
}

# The four fluorinated alcohols, all four named in the SI Materials list
# (SI PDF p.3). These are the point of the paper: F-count and F-position are
# invisible to a chain-length count, so they are pinned to the purchased
# reagent's own IUPAC name, not to the `Fn` token.
TAILS = {
 "F1": {"ester": "CCCCF",                   "n_F": 1,
        "alcohol": "4-fluoro-1-butanol", "alcohol_smiles": "OCCCCF"},
 "F3": {"ester": "CCCC(F)(F)F",             "n_F": 3,
        "alcohol": "4,4,4-trifluoro-1-butanol", "alcohol_smiles": "OCCCC(F)(F)F"},
 "F5": {"ester": "CCC(F)(F)C(F)(F)F",       "n_F": 5,
        "alcohol": "3,3,4,4,4-pentafluorobutan-1-ol",
        "alcohol_smiles": "OCCC(F)(F)C(F)(F)F"},
 "F7": {"ester": "CC(F)(F)C(F)(F)C(F)(F)F", "n_F": 7,
        "alcohol": "2,2,3,3,4,4,4-heptafluoro-1-butanol",
        "alcohol_smiles": "OCC(F)(F)C(F)(F)C(F)(F)F"},
}

# One explicit arm template per linker. `{TAIL}` is the only slot.
# `[CH2:2]` is the carbon that bore the bromine.
ARMS = {
 "C3": "[CH2:2]CCCCC(=O)OCCCC(=O)O{TAIL}",
 "C4": "[CH2:2]CCCCC(=O)OCCCCC(=O)O{TAIL}",
 "C5": "[CH2:2]CCCCC(=O)OCCCCCC(=O)O{TAIL}",
 "C6": "[CH2:2]CCCCC(=O)OCCCCCCC(=O)O{TAIL}",
 "C7": "[CH2:2]CCCCC(=O)OCCCCCCCC(=O)O{TAIL}",
 "B1": "[CH2:2]CCCCC(=O)OC(C)CCCC(=O)O{TAIL}",
 "B2": "[CH2:2]CCCCC(=O)OC(CC)CCCC(=O)O{TAIL}",
 "B4": "[CH2:2]CCCCC(=O)OC(CCCC)CCCC(=O)O{TAIL}",
 "O":  "[CH2:2]CCCCC(=O)OCCCCOCCC(=O)O{TAIL}",
 "P":  "[CH2:2]CCCCC(=O)OCc1ccc(OCCC(=O)O{TAIL})cc1",
}
# The SI's own name for the bromide reagent each arm comes from, so a reader
# can check the template against a string the authors wrote.
ARM_SI_NAME = {
 "C3": "4-(fluoroalkoxy)-4-oxobutyl 6-bromohexanoate",
 "C4": "5-(fluoroalkoxy)-5-oxopentyl 6-bromohexanoate",
 "C5": "6-oxo-6-(fluoroalkoxy)hexyl 6-bromohexanoate",
 "C6": "fluoroalkyl 7-((6-bromohexanoyl)oxy)heptanoate",
 "C7": "fluoroalkyl 8-((6-bromohexanoyl)oxy)octanoate",
 "B1": "6-oxo-6-(fluoroalkoxy)hexan-2-yl 6-bromohexanoate",
 "B2": "fluoroalkyl 5-((6-bromohexanoyl)oxy)heptanoate",
 "B4": "fluoroalkyl 5-((6-bromohexanoyl)oxy)nonanoate",
 "O":  "4-(3-(fluoroalkoxy)-3-oxopropoxy)butyl 6-bromohexanoate",
 "P":  "4-(3-(fluoroalkoxy)-3-oxopropoxy)benzyl 6-bromohexanoate",
}
# CARBON atoms strictly between the two arm ester carbonyl carbons, along the
# shortest path -- hand-counted off the template ONCE, so that
# chains.check_name_token has something to compare the `Cn` token to, and so
# that `measured_main_chain` has an expectation to contradict. Heteroatoms on
# the path (the bridging ester O, the `O` linker's ether, the `P` linker's
# phenol O) are deliberately NOT counted: the token is a carbon count.
# Verified below against the built molecule; neither side is trusted alone.
LINKER_MAIN_CHAIN = {"C3": 3, "C4": 4, "C5": 5, "C6": 6, "C7": 7,
                     "B1": 4, "B2": 4, "B4": 4, "O": 6, "P": 7}

# The roster, read off Figure S51's x-axis labels at 900 dpi (75 ticks:
# SM102 plus 74 FILs) and cross-checked against the SI text layer's names.
ROSTER = (
    [f"SSC{n}F{f}" for f in (1, 3, 5, 7) for n in (3, 4, 5, 6)]
    + [f"SSC{t}F{f}" for t in ("O", "P", "B2", "B4") for f in (1, 3, 5, 7)]
    + [f"SPC{n}F{f}" for f in (1, 3, 5, 7) for n in (3, 4, 5, 6)]
    + [f"SPC{t}F{f}" for t in ("O", "P", "B2", "B4") for f in (1, 3, 5, 7)]
    + ["SSC7F5", "SSCB1F5", "ALC4F3", "ASC4F3", "ALC6F5", "ASC6F5",
       "ALCB1F5", "ALCB2F5", "ASCB1F5", "ASCB2F5"]
)


def formula_mass(formula):
    """Monoisotopic neutral mass of a CxHyFzNaOb string, or None if unparsed."""
    import re as _re
    toks = _re.findall(r"([CHFNO])(\d*)", formula)
    if not toks:
        return None
    seen = "".join(e + n for e, n in toks)
    if seen != formula:
        return None
    m = 0.0
    for el, n in toks:
        m += MONO[el] * (int(n) if n else 1)
    return m


def parse_name(name):
    """(core, linker, fluoro) from a FIL name, or raise. Never guesses."""
    if "F" not in name[2:]:
        raise StructureError(f"parse_name: no F token in {name!r}")
    core = name[:2]
    if core not in CORES:
        raise StructureError(f"parse_name: unknown core {core!r} in {name!r}")
    rest = name[2:]
    if not rest.startswith("C"):
        raise StructureError(f"parse_name: {name!r} has no 'C' linker marker")
    body, _, ftok = rest.rpartition("F")
    # The name writes the linker as the literal "C" followed by the token, so
    # "SSC3F1" carries linker C3 while "SSCOF1" carries linker O and "SSCB2F1"
    # carries B2. Reconstructing "C"+digits rather than stripping the leading
    # "C" unconditionally keeps C3 and a hypothetical bare 3 distinguishable.
    link = body[1:]
    if link.isdigit():
        link = "C" + link
    if link not in ARMS:
        raise StructureError(f"parse_name: unknown linker {link!r} in {name!r}")
    fl = "F" + ftok
    if fl not in TAILS:
        raise StructureError(f"parse_name: unknown tail {fl!r} in {name!r}")
    return core, link, fl


def _check_no_ring_digit_collision(template, tail):
    """Refuse a composition that could fuse two fragments through a ring digit.

    This is the piperazine/cholesterol class made unwritable: the only ring in
    this library is the `P` linker's benzene, and the fluoroalkyl tails are
    acyclic, so a digit appearing in BOTH halves would mean the templates have
    drifted.
    """
    dt = {c for c in template if c.isdigit()}
    dl = {c for c in tail if c.isdigit()}
    if dt & dl:
        raise StructureError(
            "stage03: ring-closure digit(s) "
            f"{sorted(dt & dl)} appear in both the arm template and the tail. "
            "Concatenating them would fuse the fragments through a ring "
            "closure -- the defect that fused piperazine to cholesterol "
            "(CLAUDE.md section 6). Write the arm out in full instead.")


def build(name):
    core, link, fl = parse_name(name)
    tmpl, tail = ARMS[link], TAILS[fl]["ester"]
    _check_no_ring_digit_collision(tmpl, tail)
    arm = tmpl.replace("{TAIL}", tail)
    smi = assemble.attach(CORES[core]["smiles"], arm,
                          a_map=1, b_map=2, label=name)
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise StructureError(f"build({name}): product does not parse")
    assemble.check_no_radicals(mol, label=name)
    return {"name": name, "core": core, "linker": link, "fluoro": fl,
            "smiles": assemble.canonical(smi),
            "formula": CalcMolFormula(mol),
            "exact_mass": ExactMolWt(mol),
            "mh_plus": ExactMolWt(mol) + PROTON,
            "n_F": sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "F")}


def measured_main_chain(name):
    """Count the linker's main-chain atoms ON THE BUILT MOLECULE.

    Walks the shortest path between the two ester carbonyl carbons and
    subtracts them, so the number compared to the `Cn` token is measured from
    the structure rather than copied from LINKER_MAIN_CHAIN.
    """
    rec = build(name)
    m = Chem.MolFromSmiles(rec["smiles"])
    patt = Chem.MolFromSmarts("[CX3](=O)[OX2][#6]")
    cs = sorted({h[0] for h in m.GetSubstructMatches(patt)})
    if len(cs) != 3:
        raise StructureError(
            f"measured_main_chain({name}): found {len(cs)} ester carbonyls, "
            "expected 3 (lipid tail, N-side, fluoroalkyl).")
    # The two that are NOT the lipid-tail ester are the arm's pair: they are
    # the pair whose shortest path contains no nitrogen.
    best = None
    for i in range(3):
        for j in range(i + 1, 3):
            path = Chem.GetShortestPath(m, cs[i], cs[j])
            if any(m.GetAtomWithIdx(k).GetSymbol() == "N" for k in path):
                continue
            if best is None or len(path) < len(best):
                best = path
    if best is None:
        raise StructureError(
            f"measured_main_chain({name}): every carbonyl pair is separated "
            "by the nitrogen; the arm template has drifted.")
    inner = best[1:-1]            # drop the two carbonyl carbons themselves
    return sum(1 for k in inner if m.GetAtomWithIdx(k).GetSymbol() == "C")


def main():
    si = json.load(open(os.path.join(HERE, "01_si_records.json")))
    # name -> SI record. 3.2.1.7's heading reads SSCO4F1 and 3.2.1.9's reads
    # SSC4F3 (a duplicate of 3.2.1.10); both are recorded in
    # 08_contradictions.csv and keyed here to the roster name they belong to,
    # which the formula and the series pattern both pin.
    HEADING_FIX = {"SSCO4F1": "SSCOF1", "3.2.1.9": "SSC3F3"}
    byname = {}
    for r in si["records"]:
        nm = HEADING_FIX.get(r["section"]) or HEADING_FIX.get(
            r["name_in_heading"], r["name_in_heading"])
        byname.setdefault(nm, []).append(r)

    rows, blocks = [], []
    for k, v in CORES.items():
        blocks.append({"block_type": "amine_core", "token": k,
                       "smiles": v["smiles"].replace("[NH:1]", "N"),
                       "si_name": v["si_name"], "si_locus": v["si_locus"],
                       "detail": f"head={v['head']}; tail={v['tail']}"})
    for k, v in TAILS.items():
        blocks.append({"block_type": "fluoroalcohol", "token": k,
                       "smiles": v["alcohol_smiles"], "si_name": v["alcohol"],
                       "si_locus": "SI section 1 Materials (SI PDF p.3)",
                       "detail": f"n_F={v['n_F']}"})
    for k, v in ARMS.items():
        blocks.append({"block_type": "bromide_arm", "token": k,
                       "smiles": v, "si_name": ARM_SI_NAME[k],
                       "si_locus": "SI section 3.2 reagent sentences",
                       "detail": f"main_chain_between_carbonyls="
                                 f"{LINKER_MAIN_CHAIN[k]}"})

    n_ok_formula = n_ok_mz = n_si = 0
    mismatches = []
    for nm in ROSTER:
        rec = build(nm)
        mc = measured_main_chain(nm)
        if mc != LINKER_MAIN_CHAIN[rec["linker"]]:
            raise StructureError(
                f"{nm}: the arm template's main chain measures {mc} atoms but "
                f"LINKER_MAIN_CHAIN says {LINKER_MAIN_CHAIN[rec['linker']]}. "
                "One of the two is wrong; neither is allowed to win silently.")
        tok = rec["linker"]
        rec["name_token_check"] = (
            chains.check_name_token(int(tok[1:]), mc, label=nm)
            if tok[0] == "C" else
            {"name_token_length": None, "measured_length": mc,
             "matches": None,
             "note": f"token {tok!r} is not a chain length: "
                     + ("a branch size at the N-side ester carbon"
                        if tok[0] == "B" else
                        "an ether-containing linker" if tok == "O" else
                        "a benzylic/para-phenylene linker")})
        srec = byname.get(nm, [None])[0]
        rec["si_section"] = srec["section"] if srec else ""
        rec["si_pdf_page"] = srec["si_pdf_page"] if srec else ""
        rec["formula_stated"] = srec["formula_stated"] if srec else ""
        rec["mz_calc_stated"] = srec["mz_calc_stated"] if srec else ""
        rec["mz_found_stated"] = srec["mz_found_stated"] if srec else ""
        if srec and srec["formula_stated"]:
            n_si += 1
            ok_f = (srec["formula_stated"] == rec["formula"])
            n_ok_formula += ok_f
            cands = [float(x) for x in
                     srec["mz_calc_stated"].replace(" ", "").split(",")
                     if x.replace(".", "").isdigit() and x.count(".") <= 1]
            dmin = min((abs(c - rec["mh_plus"]) for c in cands),
                       default=float("nan"))
            ok_m = dmin <= 0.05
            n_ok_mz += bool(ok_m)
            rec["formula_agrees"] = int(ok_f)
            rec["mh_plus_delta_vs_SI"] = round(dmin, 4)
            rec["mz_agrees"] = int(ok_m)
            fm = formula_mass(srec["formula_stated"])
            if fm is None:
                rec["si_formula_mh"] = ""
                rec["si_self_consistent"] = ""
            else:
                simh = fm + PROTON
                sidelta = min((abs(c - simh) for c in cands),
                              default=float("nan"))
                rec["si_formula_mh"] = round(simh, 4)
                rec["si_self_consistent"] = int(sidelta <= 0.05)
            if not (ok_f and ok_m):
                mismatches.append(
                    {"name": nm, "section": srec["section"],
                     "built_formula": rec["formula"],
                     "si_formula": srec["formula_stated"],
                     "built_mh": round(rec["mh_plus"], 4),
                     "si_mz_calc": srec["mz_calc_stated"],
                     "delta": round(dmin, 4),
                     "si_formula_mh": rec.get("si_formula_mh"),
                     "si_self_consistent": rec.get("si_self_consistent")})
        else:
            rec["formula_agrees"] = ""
            rec["mh_plus_delta_vs_SI"] = ""
            rec["mz_agrees"] = ""
            rec["si_formula_mh"] = ""
            rec["si_self_consistent"] = ""
        rows.append(rec)

    # ---- write ----
    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, ["block_type", "token", "smiles", "si_name",
                                "si_locus", "detail"])
        w.writeheader()
        w.writerows(blocks)
    cols = ["name", "core", "linker", "fluoro", "smiles", "formula",
            "exact_mass", "mh_plus", "n_F", "si_section", "si_pdf_page",
            "formula_stated", "formula_agrees", "mz_calc_stated",
            "mz_found_stated", "mh_plus_delta_vs_SI", "mz_agrees",
            "si_formula_mh", "si_self_consistent"]
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            r = dict(r)
            r["exact_mass"] = round(r["exact_mass"], 5)
            r["mh_plus"] = round(r["mh_plus"], 5)
            w.writerow(r)
    with open(os.path.join(HERE, "03_name_token_checks.json"), "w") as fh:
        json.dump({r["name"]: r["name_token_check"] for r in rows}, fh, indent=1)

    smis = {r["smiles"] for r in rows}
    print(f"enumerated {len(rows)} FILs, {len(smis)} unique canonical SMILES")
    print(f"SI records matched: {n_si}; formula agrees {n_ok_formula}/{n_si}; "
          f"[M+H]+ within 0.05 Da {n_ok_mz}/{n_si}")
    for m in mismatches:
        print("  MISMATCH", m)
    return 0


if __name__ == "__main__":
    sys.exit(main())
