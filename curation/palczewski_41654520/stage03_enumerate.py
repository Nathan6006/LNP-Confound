#!/usr/bin/env python3
"""Build every CBB / CBBZ structure and check it against the SI's own HRMS.

    python3 -I stage03_enumerate.py

Writes 03_enumerated.csv. (01_blocks.csv belongs to stage01_blocks.py.)

WHAT IS CHECKED, AND WHAT CANNOT BE. For each compound the SI prints both a
molecular formula and a calculated m/z. Those are two statements, and they can
disagree with each other -- so this script tests THREE things separately:

  A  our built SMILES' neutral formula, plus z protons, == the SI formula string
  B  the m/z recomputed from the SI FORMULA STRING  == the SI's printed m/z
  C  the m/z recomputed from OUR OWN built formula  == the SI's printed m/z

A failing with C passing means the SI's formula STRING is the typo and its m/z
is right. A and C failing together means we are wrong. Collapsing these into
one "HRMS passes" boolean is how a formula typo gets silently adopted.

NONE OF THIS SEES CONNECTIVITY. Mass and formula are blind to which nitrogen an
arm sits on; that came from reading the IUPAC names, and the C10 anchor in
make_reference.py is the only external test of topology here.
"""
import csv
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
# `python3 -I` implies -P, so the script's own directory is NOT on sys.path.
sys.path.insert(0, os.path.join(_HERE, "..", "lib"))
sys.path.insert(0, _HERE)
import assemble                                              # noqa: E402
import reagents                                              # noqa: E402
import paper_library as P                                    # noqa: E402

from rdkit import Chem, RDLogger                              # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula        # noqa: E402

RDLogger.DisableLog("rdApp.*")
HERE = os.path.dirname(os.path.abspath(__file__))

M_E = 0.00054857990907
MONO = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048,
        "O": 15.9949146196, "S": 31.97207100, "F": 18.99840322}


def parse_formula(s):
    import re
    out = {}
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", s):
        if el:
            out[el] = out.get(el, 0) + (int(n) if n else 1)
    return out


def mz(counts, z):
    """m/z of a cation whose ATOM counts are `counts` and whose charge is z."""
    m = sum(MONO[e] * n for e, n in counts.items())
    return (m - z * M_E) / z


def neutral_counts(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise SystemExit("unparseable: " + smiles[:120])
    mh = Chem.AddHs(m)
    c = {}
    for a in mh.GetAtoms():
        c[a.GetSymbol()] = c.get(a.GetSymbol(), 0) + 1
    return c, Chem.GetFormalCharge(m)


def amidate(core, acyl, a_map, label):
    """R-CH=O + core's mapped primary amine -> amide. No atom is removed."""
    return assemble.attach(core, acyl, a_map=a_map, b_map=1, label=label)


def build():
    """compound -> (smiles, record). `record` is the row's own provenance."""
    out = {}

    def rec(**kw):
        return kw

    # ---- the two bis-arm cores, both terminal amines mapped (:2 and :3)
    core_k1 = P.dye(P.HEAD_ETHOXY, P.K1_A, P.K1_B)
    core_k2 = P.dye(P.HEAD_ETHOXY, P.K2_A, P.K2_B)

    def bis_acyl(core, acyl, label):
        one = amidate(core, acyl, 2, label + ":arm1")
        return amidate(one, acyl, 3, label + ":arm2")

    # ---------------------------------------------------------- CBB21, CBB22
    out["CBB21"] = (assemble.canonical(P.dye(P.HEAD_ETHOXY,
                                             "NCCCN(C)CCCN", "NCCCN(C)CCCN")),
        rec(scaffold="CBB-K1 bis-sulfonamide", tail="none (free primary amines)",
            confidence="high",
            note="CBB21: Coomassie Brilliant Blue G-250 bis-sulfonamide of "
                 "3,3'-diamino-N-methyl-dipropylamine; built from its own SI "
                 "IUPAC name, HRMS C61H84N9O5S2+ exact."))
    out["CBB22"] = (assemble.canonical(P.dye(P.HEAD_ETHOXY, P.SO3H, P.SO3H)),
        rec(scaffold="CBB G-250 free acid", tail="none",
            confidence="high",
            note="CBB22 is commercial Coomassie Brilliant Blue G-250, emitted "
                 "as the neutral free acid C47H49N3O7S2 (832.00); the sodium "
                 "salt is 854.01 and the SI weighs 'brilliant G250 (854 mg, "
                 "1.0 mmol)', which corroborates the core independently."))

    # ------------------------------------------------- the bis-amide series
    simple = [("CBB3",  "dha"),            ("CBB6",  "palmitoyl"),
              ("CBB8",  "15-hexadecynoyl"),("CBB10", "triacontanoyl"),
              ("CBB11", "linoleoyl"),      ("CBB12", "lauroyl"),
              ("CBB13", "butyryl"),        ("CBB17", "trifluoroacetyl"),
              ("CBB20", "retinoyl")]
    for name, acyl in simple:
        smi, why = P.ACYL[acyl]
        out[name] = (bis_acyl(core_k1, smi, name),
            rec(scaffold="CBB-K1 bis-sulfonamide", tail=acyl, tail_source=why,
                confidence="high",
                note=f"{name}: CBB21 bis-acylated on both terminal primary "
                     f"amines with {acyl} ({why}); connectivity from the SI "
                     f"IUPAC name, arm count from its HRMS."))

    # --------------------------------- the quaternised (K2) pairs 2/4/7/14
    for name, parent, acyl in [("CBB2", "CBB1", None),
                               ("CBB4", "CBB3", "dha"),
                               ("CBB7", "CBB6", "palmitoyl"),
                               ("CBB14", "CBB13", "butyryl")]:
        if acyl is None:
            continue                                   # CBB2 handled with CBB1
        smi, why = P.ACYL[acyl]
        out[name] = (bis_acyl(core_k2, smi, name),
            rec(scaffold="CBB-K2 bis-sulfonamide (both linker tertiary amines "
                         "N-methylated to quaternary dimethylammonium)",
                tail=acyl, tail_source=why, confidence="high",
                note=f"{name}: {parent} exhaustively methylated with methyl "
                     f"iodide; the SI name says 'dimethylammonio' and "
                     f"'N,N-dimethylpropan-1-aminium', so both linker "
                     f"nitrogens are quaternary. Emitted as the dication."))

    # ------------------------------------- CBB1 / CBB2, cholesteryl carbamate
    # formic acid onto each terminal amine -> carbamic acid, then esterify with
    # the REGISTRY cholesterol. Cholesterol keeps its own oxygen; its C3 is
    # never a neighbour of a deleted atom, which is the F-013 rule.
    chol = reagents.smi("Cholesterol")
    formyl = P.ACYL["formyl_for_carbamate"][0]
    for name, core in (("CBB1", core_k1), ("CBB2", core_k2)):
        m = amidate(core, formyl, 2, name + ":carbamic1")
        m = assemble.esterify(m, chol, label=name + ":chol1")
        # re-map: the second arm's map :3 survived both operations
        m2 = assemble.attach(m, formyl, a_map=3, b_map=1,
                             label=name + ":carbamic2")
        m2 = assemble.esterify(m2, chol, label=name + ":chol2")
        out[name] = (m2, rec(
            scaffold=("CBB-K1 bis-sulfonamide" if name == "CBB1"
                      else "CBB-K2 bis-sulfonamide (quaternised)"),
            tail="cholesteryl carbamate", tail_source="cholesteryl chloroformate",
            confidence="high",
            note=f"{name}: both terminal amines carbamoylated with cholesteryl "
                 "chloroformate. Cholesterol is taken verbatim from "
                 "provenance/reagent_registry.py and esterified by "
                 "assemble.esterify, so the acid's OH leaves and C3-beta is "
                 "untouched; the SI name's (3S,...)-cyclopenta[a]phenanthren-3-yl "
                 "is cholesteryl, 3-beta."))

    # ------------------------------------------------- CBB5, the mono-adduct
    dha = P.ACYL["dha"][0]
    core5 = P.dye(P.HEAD_ETHOXY, P.K1_A, P.SO3H)
    out["CBB5"] = (amidate(core5, dha, 2, "CBB5"), rec(
        scaffold="CBB G-250 mono-sulfonamide / mono-sulfonic acid",
        tail="dha (one arm only)", confidence="high",
        note="CBB5: the silica-column byproduct of the CBB3 synthesis -- only "
             "ONE sulfonic acid was amidated, the other stays as free "
             "benzenesulfonic acid. The SI name and its HRMS N6 (not N9) both "
             "say so."))

    # ------------------------------------------------ CBB15 / CBB16, TEG esters
    out["CBB15"] = (assemble.canonical(
        P.dye(P.HEAD_ETHOXY, P.TEG_PALMITATE, P.TEG_PALMITATE)), rec(
        scaffold="CBB G-250 bis(triethylene-glycol sulfonic ester)",
        tail="palmitoyl on both TEG termini", confidence="high",
        note="CBB15: made by sulfonic ESTERIFICATION of CBB G-250 with "
             "triethylene glycol, then palmitoylation of the free OH. No "
             "nitrogen is added, and its HRMS N3 confirms that."))
    out["CBB16"] = (assemble.canonical(
        P.dye(P.HEAD_ETHOXY, P.TEG_PALMITATE, P.SO3H)), rec(
        scaffold="CBB G-250 mono(TEG sulfonic ester) / mono-sulfonic acid",
        tail="palmitoyl on the one TEG terminus", confidence="high",
        note="CBB16: the mono-esterified byproduct isolated alongside CBB15; "
             "HRMS C69H92N3O11S2+."))

    # ------------------------------------------------------ CBB9, macrocycle
    out["CBB9"] = (assemble.canonical(P.CBB9_MACROCYCLE), rec(
        scaffold="CBB-K1 bis-sulfonamide, macrocyclised",
        tail="oxalyl bridge (intramolecular)", confidence="high",
        note="CBB9: the byproduct of the CBB10 synthesis, in which excess "
             "oxalyl chloride bridged both terminal amines into a macrocyclic "
             "oxalamide ('17,18-dione' on adjacent skeletal atoms in the SI "
             "name). Closing a macrocycle needs two bonds, so the whole "
             "molecule is transcribed from the name rather than assembled, "
             "then canonicalised and radical-checked."))

    # --------------------------------------------- CBB18 / CBB19, CPP5 peptide
    for name in ("CBB18", "CBB19"):
        for seq, (pep, why) in P.CPP5.items():
            if name == "CBB19":
                m = bis_acyl(core_k1, pep, name + ":" + seq)
            else:
                m = amidate(core_k1, pep, 2, name + ":" + seq + ":pep")
                m = assemble.attach(m, P.ACYL["trifluoroacetyl"][0],
                                    a_map=3, b_map=1,
                                    label=name + ":" + seq + ":tfa")
            out[name + "_" + seq] = (m, rec(
                scaffold="CBB-K1 bis-sulfonamide",
                tail=("CPP5 x2" if name == "CBB19"
                      else "CPP5 x1 + trifluoroacetyl x1"),
                tail_source=why, cpp5_sequence=seq,
                confidence="low", structure_assumed=1,
                note=f"{name}: the CPP5 residue COMPOSITION is fixed by HRMS, "
                     f"but the SEQUENCE is not -- the SI heading writes KLVPM "
                     f"and the main text writes KLPVM, and the two permutations "
                     f"are isomers that no mass can separate. This candidate is "
                     f"{seq}. See DECISIONS.md D-1; nothing is emitted as THE "
                     f"structure until the user settles it."))

    # ----------------------------------------------------- CBBZ0 and CBBZ1-10
    out["CBBZ0"] = (assemble.canonical(
        P.dye(P.HEAD_PROPYNYL, P.SO3H, P.SO3H)), rec(
        scaffold="CBBZ propargyl-ether dye", tail="none (free alkyne)",
        confidence="medium",
        note="CBBZ0: the click handle, made by the Takahagi route with "
             "4-(prop-2-yn-1-yloxy)aniline in place of 4-ethoxyaniline. "
             "confidence=medium because it is the ONE library member the SI "
             "gives no HRMS for -- 1H/13C NMR only -- so its atom count has no "
             "independent check, and its SI name is written in the "
             "aminophenyl/ethanaminium tautomer while every other name in the "
             "paper uses the quinoid-imine form we emit."))
    phenol_core = P.dye("[OH:1]c4ccc(cc4)", P.SO3H, P.SO3H)
    for zn, (tail, why) in P.AZIDE_TAIL.items():
        name = "CBBZ" + zn[1:]
        out[name] = (assemble.attach(phenol_core, P.triazole_frag(tail),
                                     a_map=1, b_map=2, label=name), rec(
            scaffold="CBBZ propargyl-ether dye", tail="triazole-" + zn,
            tail_source=why, confidence="high",
            note=f"{name}: CuAAC of CBBZ0 with {why}. Built by etherifying the "
                 f"dye's phenolic oxygen with the 1,4-disubstituted "
                 f"triazol-4-ylmethyl the cycloaddition forms, which "
                 f"reconstitutes Ar-O-CH2-triazole without needing a "
                 f"ring-forming step."))

    # ---------------------------------------------------------- CBBox11
    out["CBBox11"] = ("", rec(
        scaffold="CBB-K1 bis-sulfonamide", tail="linoleoyl, mono-hydroxylated",
        confidence="low", structure_assumed=1, positional_ambiguity=1,
        note="CBBox11: the DMSO storage oxidation product of CBB11. The SI "
             "reports two chromatographically separated isomers with similar "
             "MS1 and says only that MS2 'suggested hydroxylation events on one "
             "of the polyunsaturated fatty acid moieties' -- which arm, and "
             "which carbon, are not determined. NO SMILES IS EMITTED: writing "
             "one would be a guess at an allylic position the paper does not "
             "assign. Its Fig 2d value is kept. See DECISIONS.md D-2."))
    return out


def main():
    P.triazole_frag = P._triazole
    built = build()

    rows, probs = [], []
    for name in sorted(built, key=lambda s: (len(s), s)):
        smi, r = built[name]
        if not smi:
            rows.append(dict(compound=name, canonical_SMILES="",
                             formula="", neutral_formula="", formal_charge="",
                             hrms_formula="", hrms_charge="", hrms_mz_printed="",
                             check_A_formula="NO STRUCTURE EMITTED",
                             check_B_mz_from_si_formula="",
                             check_C_mz_from_our_formula="", **r))
            continue
        counts, chg = neutral_counts(smi)
        mol = Chem.MolFromSmiles(smi)
        assemble.check_no_radicals(mol, label=name)
        # CBB18/CBB19 carry a sequence suffix; the SI's HRMS is for the
        # compound, and BOTH sequence candidates have the same formula -- which
        # is the point: running the check on both is what demonstrates that
        # mass cannot separate them.
        key = name.split("_")[0]
        A = B = C = "no HRMS in SI"
        hf = hz = hm = ""
        if key in P.HRMS:
            hf, hz, hm = P.HRMS[key]
            si = parse_formula(hf)
            # our cation at the SI's charge: neutral + (z - formal charge) H
            ours = dict(counts)
            add = hz - chg
            ours["H"] = ours.get("H", 0) + add
            A = "exact" if ours == si else ("MISMATCH " + json.dumps(
                {e: (ours.get(e, 0), si.get(e, 0)) for e in set(ours) | set(si)
                 if ours.get(e, 0) != si.get(e, 0)}))
            B = round(mz(si, hz) - hm, 4)
            C = round(mz(ours, hz) - hm, 4)
            if A != "exact" or abs(B) > 0.002 or abs(C) > 0.002:
                probs.append((name, A, B, C))
        rows.append(dict(
            compound=name, canonical_SMILES=smi, formula=CalcMolFormula(mol),
            neutral_formula="".join(f"{e}{counts[e]}" for e in
                                    sorted(counts, key=lambda x: "CHNOFS".index(x[0]))),
            formal_charge=chg, hrms_formula=hf, hrms_charge=hz,
            hrms_mz_printed=hm, check_A_formula=A,
            check_B_mz_from_si_formula=B, check_C_mz_from_our_formula=C, **r))

    cols = ["compound", "canonical_SMILES", "formula", "neutral_formula",
            "formal_charge", "scaffold", "tail", "tail_source",
            "cpp5_sequence", "hrms_formula", "hrms_charge", "hrms_mz_printed",
            "check_A_formula", "check_B_mz_from_si_formula",
            "check_C_mz_from_our_formula", "confidence", "structure_assumed",
            "positional_ambiguity", "note"]
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    # 01_blocks.csv is NOT written here. It was, until stage01_blocks.py was
    # written to build it from the paper's own drawn Markush table (SI Fig 1)
    # -- and because the stages run in order, stage03's weaker version then
    # silently overwrote stage01's every single time. Two writers for one file
    # is a defect that only shows up as "the file is not deterministic", which
    # is how it was found here.
    print(f"built {len([r for r in rows if r['canonical_SMILES']])} structures,"
          f" {len(rows)} rows")
    print(f"HRMS checked: {sum(1 for r in rows if r['hrms_formula'])}")
    for p in probs:
        print("  PROBLEM", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
