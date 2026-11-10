#!/usr/bin/env python3
"""Decompose the drawn lipids into blocks, then RE-ENUMERATE and compare.

    python3 -I stage02_blocks.py

Reads  02_structures.csv (stage01).
Writes 01_blocks.csv, 03_enumerated.csv, 02_blocks_consensus.json.

WHY DECOMPOSE RATHER THAN ENUMERATE FROM READ-ONCE BLOCKS. The SI draws each
lipid in full, so every headgroup is drawn ten times over -- once per tail --
and every tail ten times over, once per headgroup. Deriving the blocks from the
lipids therefore gives up to ten INDEPENDENT reads of each block, and a
disagreement among them is a detected error rather than an undetected one. A
block read once from one picture has no such check.

The cut is chemical, not positional: a bond from a nitrogen to a carbon whose
side of the bond contains no nitrogen at all is an arm. That is the bond the
epoxide-opening made ("A microwave vial was charged with the epoxy ester ...
followed by the headgroup amine ... 90 C for 1 hour", SI p.5), and it is the
only N-C bond in these molecules that disconnects anything.

MASS CANNOT SEE CONNECTIVITY, so the m/z gate in stage01 does not establish
that a headgroup is 1,4,7-triazonane rather than 1-(2-aminoethyl)piperazine --
both are C6H15N3. What this stage adds is that the SAME graph comes out of ten
differently-drawn molecules, and that re-assembling the blocks reproduces each
drawn lipid bond for bond.
"""
import csv
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "lib"))

from rdkit import Chem                                    # noqa: E402
from rdkit.Chem.rdMolDescriptors import CalcMolFormula    # noqa: E402
import assemble                                           # noqa: E402

# Supplementary Table 2 gives the naming rule by example: AMG11 = Am1 + G1,
# AMG510 = Am5 + G10, AMG104 = Am10 + G4, AMG1541 = Am15 + G41, AMG1741 =
# Am17 + G41. The split is therefore NOT positional and must be resolved
# against the known block sets; a label whose split is not unique is refused.
AM_IDS = [str(i) for i in range(1, 20)] + ["15L"]
G_IDS = ([str(i) for i in range(1, 11)]
         + ["41", "42", "43", "44", "47", "48", "49", "50", "51", "53", "54",
            "55", "56", "59", "60", "61"]
         + ["C12", "C14", "C16", "C18"])


# Labels that carry a VALUE in a figure but are not in the 10x10 grid.
EXTRA_TARGETS = ["AMG1041", "AMG1042", "AMG1043", "AMG1044", "AMG10C12",
                 "AMG141", "AMG1141", "AMG1241", "AMG1341", "AMG1441",
                 "AMG1541", "AMG1641", "AMG1741", "AMG1841", "AMG1941",
                 "AMG1547", "AMG1548", "AMG1549", "AMG1550", "AMG1551",
                 "AMG1553", "AMG1554", "AMG1555", "AMG1556"]

# A label whose Am/G split is not unique from the name alone. Resolved by
# STRUCTURE, not by preference: the drawn molecule's own headgroup and arm are
# compared with the consensus blocks, and the split that matches both wins.
# AMG141 is Am1+G41 or Am14+G1 on the name; Fig 2b puts it in a series of
# headgroups all carrying G41, and the drawing agrees.
AMBIGUOUS = {}


def label_splits(lab):
    """Every (Am, G, suffix) the NAME allows. Often one; AMG141 allows two."""
    m = re.fullmatch(r"AMG([0-9A-Za-z]+)(-OHD|R|S)?", lab)
    if not m:
        return []
    body, suffix = m.group(1), m.group(2) or ""
    return [(a, body[len(a):], suffix) for a in AM_IDS
            if body.startswith(a) and body[len(a):] in G_IDS]


def split_label(lab):
    """The unique (Am, G, suffix), or None when the name allows more."""
    hits = label_splits(lab)
    return hits[0] if len(hits) == 1 else None


def arm_bonds(mol):
    """[(n_idx, c_idx)] for every N-C bond whose C side holds no nitrogen."""
    out = []
    for b in mol.GetBonds():
        a1, a2 = b.GetBeginAtom(), b.GetEndAtom()
        if {a1.GetSymbol(), a2.GetSymbol()} != {"N", "C"}:
            continue
        n, c = (a1, a2) if a1.GetSymbol() == "N" else (a2, a1)
        em = Chem.RWMol(mol)
        em.RemoveBond(a1.GetIdx(), a2.GetIdx())
        frags = Chem.GetMolFrags(em.GetMol())
        if len(frags) != 2:
            continue
        side = next(f for f in frags if c.GetIdx() in f)
        if any(mol.GetAtomWithIdx(i).GetSymbol() == "N" for i in side):
            continue
        out.append((n.GetIdx(), c.GetIdx()))
    return out


def decompose(smiles):
    """(head_smiles_with_[*:3], [arm_smiles_with_[*:2]]) or raise.

    The head's free valences are tagged 3, not 1: every one of them is a map
    number, and assemble.attach requires EXACTLY ONE atom carrying the map it
    is given. Tagging them all 1 made it refuse, correctly.
    """
    mol = Chem.MolFromSmiles(smiles)
    cuts = arm_bonds(mol)
    if not cuts:
        raise ValueError("no arm bonds found")
    rw = Chem.RWMol(mol)
    for n_i, c_i in cuts:
        rw.RemoveBond(n_i, c_i)
        rw.GetAtomWithIdx(n_i).SetAtomMapNum(3)
        rw.GetAtomWithIdx(c_i).SetAtomMapNum(2)
        for idx in (n_i, c_i):            # let valence put the H back
            a = rw.GetAtomWithIdx(idx)
            a.SetNoImplicit(False)
            a.SetNumExplicitHs(0)
            a.SetNumRadicalElectrons(0)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    pieces = Chem.GetMolFrags(out, asMols=True)
    heads = [p for p in pieces
             if any(a.GetSymbol() == "N" for a in p.GetAtoms())]
    arms = [p for p in pieces
            if not any(a.GetSymbol() == "N" for a in p.GetAtoms())]
    if len(heads) != 1:
        raise ValueError(f"{len(heads)} N-containing fragments after the cut")
    return Chem.MolToSmiles(heads[0]), sorted(Chem.MolToSmiles(a) for a in arms)


def _unmapped(smiles):
    m = Chem.MolFromSmiles(smiles)
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    return Chem.MolFromSmiles(Chem.MolToSmiles(m))


def attach_arm(head_smiles, arm_smiles):
    """One arm onto the FIRST free [*:3] of the head, via curation/lib."""
    head = Chem.MolFromSmiles(head_smiles)
    rw = Chem.RWMol(head)
    first = None
    for a in rw.GetAtoms():
        if a.GetAtomMapNum() == 3:
            first = a.GetIdx()
            break
    if first is None:
        raise ValueError("head has no free attachment point")
    rw.GetAtomWithIdx(first).SetAtomMapNum(1)
    head_open = Chem.MolToSmiles(rw.GetMol())
    return safe_attach(head_open, arm_smiles, label="AMG arm")


STEREO_WAIVED = []


def safe_attach(head, arm, label):
    """assemble.attach, with ONE documented waiver and no others.

    `check_stereo_preserved` requires the product's potential-stereocentre
    multiset to equal the sum of the fragments'. That is right for a join
    between two finished fragments and WRONG for one that CREATES a centre:
    alkylating the second nitrogen of Am7's cis-fused bicycle makes both
    ring-fusion carbons stereocentres, and putting an arm on Am8 does the same
    to its CH-OH. The assertion fires on 24 of these builds.

    The waiver is narrow and checked, not assumed: it is accepted only when
    NEITHER the fragments NOR the product carry a single ASSIGNED (R/S)
    centre, so nothing can have been inverted -- which is the failure
    `check_stereo_preserved` exists to catch (F-013, 537 rows of
    epicholesterol). Every waiver is recorded in 02_blocks_consensus.json.
    """
    try:
        return assemble.attach(head, arm, a_map=1, b_map=2, label=label)
    except Exception as exc:                                # noqa: BLE001
        if "check_stereo_preserved" not in str(exc):
            raise
        h = Chem.MolFromSmiles(head)
        a = Chem.MolFromSmiles(arm)
        # The map numbers are the build's attachment points and must SURVIVE:
        # clearing them on the working copies left the head with no free site
        # and every multi-arm Am7/Am8 lipid failed on its second arm.
        unmapped = []
        for frag in (h, a):
            f = Chem.Mol(frag)
            for at in f.GetAtoms():
                at.SetAtomMapNum(0)
            unmapped.append(f)
        assigned = [c for f in unmapped
                    for c in assemble.stereocentres(f).values() if c != "?"]
        if assigned:
            raise
        rw = Chem.RWMol(Chem.CombineMols(h, a))
        ia = next(i for i, at in enumerate(h.GetAtoms())
                  if at.GetAtomMapNum() == 1)
        ib = next(i for i, at in enumerate(a.GetAtoms())
                  if at.GetAtomMapNum() == 2) + h.GetNumAtoms()
        rw.AddBond(ia, ib, Chem.BondType.SINGLE)
        for idx in (ia, ib):
            at = rw.GetAtomWithIdx(idx)
            at.SetAtomMapNum(0)
            at.SetNumExplicitHs(0)
            at.SetNoImplicit(False)
            at.SetNumRadicalElectrons(0)
        prod = rw.GetMol()
        Chem.SanitizeMol(prod)
        assemble.check_no_radicals(prod, label=label)
        if [c for c in assemble.stereocentres(prod).values() if c != "?"]:
            raise
        STEREO_WAIVED.append(label)
        return Chem.MolToSmiles(prod)


def build(head_smiles, arm_smiles, n_arms):
    smi = head_smiles
    for _ in range(n_arms):
        smi = attach_arm(smi, arm_smiles)
    mol = Chem.MolFromSmiles(smi)
    assemble.check_no_radicals(mol, label="AMG")
    return Chem.MolToSmiles(mol)


def main():
    rows = [r for r in csv.DictReader(open(os.path.join(HERE,
                                                        "02_structures.csv")))]
    drawn = {}
    for r in rows:
        if r["status"] == "PASS" or (r["status"] in ("FAIL", "NO_HRMS")
                                     and r["SMILES"]):
            drawn.setdefault(r["label"], r)

    heads, tails, per_lipid, problems = {}, {}, {}, []
    for lab, r in sorted(drawn.items()):
        sp = split_label(lab)
        if sp is None:
            problems.append((lab, "label does not split uniquely"))
            continue
        am, g, suffix = sp
        if suffix:                       # R/S/-OHD variants are not block-built
            continue
        try:
            head, arms = decompose(r["SMILES"])
        except Exception as exc:                            # noqa: BLE001
            problems.append((lab, f"decompose: {exc}"))
            continue
        if len(set(arms)) != 1:
            problems.append((lab, f"{len(set(arms))} distinct arms"))
            continue
        per_lipid[lab] = {"Am": am, "G": g, "head": head, "arm": arms[0],
                          "n_arms": len(arms), "status": r["status"],
                          "drawn_SMILES": r["SMILES"]}
        heads.setdefault(am, {}).setdefault(head, []).append(lab)
        tails.setdefault(g, {}).setdefault(arms[0], []).append(lab)

    # --- labels the NAME cannot split: resolve by STRUCTURE ---------------
    # AMG141 is Am1+G41 or Am14+G1 on the name alone. The drawn molecule knows
    # which: decompose it and ask which candidate's head AND arm both match a
    # block already established from unambiguous lipids. Resolving by
    # preference instead would be the AC_2025 mistake -- a confident,
    # internally consistent choice that puts measured values on the wrong
    # molecule.
    prelim_head = {am: max(v.items(), key=lambda kv: len(kv[1]))[0]
                   for am, v in heads.items()}
    prelim_tail = {g: max(v.items(), key=lambda kv: len(kv[1]))[0]
                   for g, v in tails.items()}
    resolved = {}
    for lab, why in list(problems):
        if why != "label does not split uniquely":
            continue
        r = drawn.get(lab)
        try:
            head, arms = decompose(r["SMILES"])
        except Exception:                                   # noqa: BLE001
            continue
        if len(set(arms)) != 1:
            continue
        ok = [(a, g) for a, g, sfx in label_splits(lab) if not sfx
              and prelim_head.get(a) == head and prelim_tail.get(g) == arms[0]]
        if len(ok) != 1:
            continue
        am, g = ok[0]
        resolved[lab] = {"Am": am, "G": g,
                         "rejected": [f"Am{a}+G{gg}"
                                      for a, gg, sfx in label_splits(lab)
                                      if not sfx and (a, gg) != (am, g)]}
        per_lipid[lab] = {"Am": am, "G": g, "head": head, "arm": arms[0],
                          "n_arms": len(arms), "status": r["status"],
                          "drawn_SMILES": r["SMILES"]}
        heads.setdefault(am, {}).setdefault(head, []).append(lab)
        tails.setdefault(g, {}).setdefault(arms[0], []).append(lab)
        problems.remove((lab, why))

    # --- consensus, and every disagreement named -------------------------
    cons_head, cons_tail, disagree = {}, {}, []
    for am, variants in heads.items():
        best = max(variants.items(), key=lambda kv: len(kv[1]))
        cons_head[am] = best[0]
        if len(variants) > 1:
            disagree.append({"block": f"Am{am}", "n_variants": len(variants),
                             "variants": {k: v for k, v in variants.items()}})
    for g, variants in tails.items():
        best = max(variants.items(), key=lambda kv: len(kv[1]))
        cons_tail[g] = best[0]
        if len(variants) > 1:
            disagree.append({"block": f"G{g}", "n_variants": len(variants),
                             "variants": {k: v for k, v in variants.items()}})

    n_arms_for = {}
    for lab, d in per_lipid.items():
        n_arms_for.setdefault(d["Am"], {}).setdefault(d["n_arms"], []).append(lab)

    # --- 01_blocks.csv ----------------------------------------------------
    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["block_id", "kind", "SMILES_with_attachment_map",
                    "n_independent_reads", "n_variants_seen",
                    "n_arms_on_this_head", "free_amine_SMILES", "formula",
                    "read_from"])
        for am in sorted(cons_head, key=lambda x: (len(x), x)):
            smi = cons_head[am]
            fm = Chem.MolFromSmiles(smi)
            for a in fm.GetAtoms():
                a.SetAtomMapNum(0)
            free = Chem.MolToSmiles(fm)
            na = sorted(n_arms_for.get(am, {}))
            w.writerow([f"Am{am}", "amine_headgroup", smi,
                        len(heads[am][smi]), len(heads[am]),
                        "|".join(str(x) for x in na), free,
                        CalcMolFormula(Chem.MolFromSmiles(free)),
                        ";".join(sorted(heads[am][smi]))])
        for g in sorted(cons_tail, key=lambda x: (len(x), x)):
            smi = cons_tail[g]
            w.writerow([f"G{g}", "tail_arm_as_attached", smi,
                        len(tails[g][smi]), len(tails[g]), "", "",
                        CalcMolFormula(_unmapped(smi)),
                        ";".join(sorted(tails[g][smi]))])

    # --- re-enumerate and compare ----------------------------------------
    # Only labels the paper actually names are built. The 10 x 10 grid is the
    # library the paper states it made ("a combinatorial library of 100 unique
    # ionizable lipids ... 10 amine headgroups (Am), 9 glycidyl esters and 1
    # epoxy-ester tails (G)"); the rest are the named later generations. The
    # full 15 x 28 cross product would be 420 molecules that do not exist.
    targets = [(str(a), str(g)) for a in range(1, 11) for g in range(1, 11)]
    for lab in sorted(per_lipid):
        t = (per_lipid[lab]["Am"], per_lipid[lab]["G"])
        if t not in targets:
            targets.append(t)
    for lab in EXTRA_TARGETS:
        sp = split_label(lab)
        if sp and (sp[0], sp[1]) not in targets:
            targets.append((sp[0], sp[1]))

    enum_rows, agree, differ, built_only = [], 0, [], 0
    for am, g in targets:
        if True:
            lab = f"AMG{am}{g}"
            if am not in cons_head or g not in cons_tail:
                enum_rows.append({"label": lab, "Am": f"Am{am}", "G": f"G{g}",
                                  "SMILES": "",
                                  "formula": "NO_BLOCK:"
                                  + ("Am" + am if am not in cons_head
                                     else "G" + g),
                                  "drawn_in_SI": lab in per_lipid,
                                  "matches_drawing": ""})
                continue
            na = max(n_arms_for.get(am, {1: []}))
            try:
                smi = build(cons_head[am], cons_tail[g], na)
            except Exception as exc:                        # noqa: BLE001
                enum_rows.append({"label": lab, "Am": f"Am{am}", "G": f"G{g}",
                                  "SMILES": "", "formula": f"BUILD_FAIL:{exc}",
                                  "drawn_in_SI": lab in per_lipid,
                                  "matches_drawing": ""})
                continue
            d = per_lipid.get(lab)
            match = ""
            if d:
                match = (d["drawn_SMILES"] ==
                         Chem.MolToSmiles(Chem.MolFromSmiles(d["drawn_SMILES"])))
                match = (smi == Chem.MolToSmiles(
                    Chem.MolFromSmiles(d["drawn_SMILES"])))
                if match:
                    agree += 1
                else:
                    differ.append(lab)
            else:
                built_only += 1
            enum_rows.append({
                "label": lab, "Am": f"Am{am}", "G": f"G{g}", "SMILES": smi,
                "formula": CalcMolFormula(Chem.MolFromSmiles(smi)),
                "n_arms": na, "drawn_in_SI": bool(d),
                "matches_drawing": match,
                "provenance": ("drawn_and_reassembled" if d
                               else "enumerated_from_blocks")})

    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as fh:
        cols = ["label", "Am", "G", "n_arms", "formula", "provenance",
                "drawn_in_SI", "matches_drawing", "SMILES"]
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in enum_rows:
            w.writerow({k: r.get(k, "") for k in cols})

    json.dump({"heads": cons_head, "tails": cons_tail,
               "stereo_waivers": sorted(set(STEREO_WAIVED)),
               "n_stereo_waivers": len(STEREO_WAIVED),
               "resolved_by_structure": resolved,
               "block_disagreements": disagree, "problems": problems,
               "n_arms_per_head": {k: {str(n): v for n, v in d.items()}
                                   for k, d in n_arms_for.items()}},
              open(os.path.join(HERE, "02_blocks_consensus.json"), "w"),
              indent=1)

    print(f"lipids decomposed       : {len(per_lipid)}")
    print(f"labels split by structure: {len(resolved)} "
          + ", ".join(f"{k}=Am{v['Am']}+G{v['G']}" for k, v in resolved.items()))
    print(f"stereo waivers (new centre created by the bond): "
          f"{len(STEREO_WAIVED)}")
    print(f"amine blocks (Am)       : {len(cons_head)}")
    print(f"tail blocks  (G)        : {len(cons_tail)}")
    print(f"blocks with >1 variant  : {len(disagree)}"
          + (f"  {[d['block'] for d in disagree]}" if disagree else ""))
    print(f"enumerated              : {len(enum_rows)}")
    print(f"  re-enumeration == drawing : {agree}")
    print(f"  re-enumeration != drawing : {len(differ)} {differ}")
    print(f"  enumerated, not drawn     : {built_only}")
    if problems:
        print(f"problems ({len(problems)}):")
        for p in problems:
            print("   ", p)
    return 0


if __name__ == "__main__":
    sys.exit(main())
