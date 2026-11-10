#!/usr/bin/env python3
"""Stage 03 -- build the 29 blocks and enumerate the 100 AID-lipids.

EVERY CHAIN LENGTH IN THIS FILE IS READ FROM 01/02's MEASUREMENT JSON. No
integer describing a chain is typed here. What IS written here is TOPOLOGY --
which atom is bonded to which -- read off Fig. 1a, because that is the thing a
pixel count cannot give you and a mass cannot check (CLAUDE.md section 5,
"mass cannot see connectivity").

THE REACTION, from Fig. 1a and Supplementary Fig. 1a:

    R-NH2  +  2-iminothiolane  ->  R-NH-C(=NH)-CH2CH2CH2-SH   (fast, amidine)
    R'-SH  +  CH2=CH-CO2R''    ->  R'-S-CH2CH2-CO2R''         (thiol-Michael)

so an AID-lipid arm is  -NH-C(=NH)-(CH2)3-S-(CH2)2-C(=O)O-alkyl, and the two
CH2 runs in it are measured in stage 02 off the one product the paper draws.

HOW MANY ARMS. One per PRIMARY amine. Supplementary Fig. 1a draws 19T-O12 with
an arm on each of amine 19's two primary nitrogens and none on its tertiary
N-methyl, and prints the mass spectrum that confirms it (07_qc.json P1b). The
cases that drawing does NOT settle are in DECISIONS.md: amines with a SECONDARY
NH as well (18, 24), and amines with three or four primary nitrogens (20, 22,
25), where the printed stoichiometry is sub-stoichiometric for exhaustive
substitution. Those rows are flagged on the row, not in this docstring.

Assembly goes through curation/lib/assemble.attach, which adds a bond between
two atom-MAPPED atoms and deletes nothing -- so the F-013 mechanism
(RemoveAtom beside a stereocentre) cannot occur here at all. The one
stereocentre in the corpus of blocks is amine 6's C=C, and
check_stereo_preserved runs on every build.

Run: python3 -I 03_enumerate.py
"""
import csv
import json
import os
import sys

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import assemble                                            # noqa: E402
import emit                                                # noqa: E402
from errors import StructureError                          # noqa: E402

CHAINS = os.path.join(HERE, "01_chain_measurements.json")
LINKER = os.path.join(HERE, "02_si_linker_measurements.json")
BLOCKS_CSV = os.path.join(HERE, "01_blocks.csv")
ENUM_CSV = os.path.join(HERE, "03_enumerated.csv")

PRIMARY_AMINE = Chem.MolFromSmarts("[NX3;H2;!$(N=*)]")


def load():
    with open(CHAINS) as fh:
        c = json.load(fh)["measurements"]
    with open(LINKER) as fh:
        k = json.load(fh)
    return c, k


def n(c, key, field="heavy_atoms_in_chain"):
    """One measured integer, by name. Raises if stage 01 did not measure it."""
    if key not in c:
        raise KeyError(f"03_enumerate: stage 01 has no measurement {key!r}.")
    v = c[key][field]
    if not isinstance(v, int):
        raise TypeError(f"03_enumerate: {key}.{field} is {v!r}, not an int.")
    return v


def build_blocks(c, k):
    """The 25 amines and 4 acrylates. Topology typed; every length measured."""
    A = {}

    # --- 1-5: n-alkylamines ------------------------------------------------
    for i in (1, 2, 3, 4, 5):
        A[i] = ("C" * n(c, f"amine0{i}_chain") + "N",
                f"n-alkylamine, C{n(c, f'amine0{i}_chain')}")

    # --- 6: mono-unsaturated alkylamine ------------------------------------
    n6 = n(c, "amine06_chain")
    d6 = c["amine06_double_bond"]["double_bond_position_from_amine_end"]
    geom = c["amine06_double_bond"]["geometry"]
    if geom not in ("Z (cis)", "E (trans)"):
        raise StructureError(f"03_enumerate: amine 6 geometry {geom!r}.")
    bond = "/C=C\\" if geom.startswith("Z") else "/C=C/"
    A[6] = ("C" * (n6 - d6 - 1) + bond + "C" * (d6 - 1) + "N",
            f"mono-unsaturated alkylamine, C{n6}, delta-{d6}, {geom}")

    # --- 7-8: fatty hydrazides --------------------------------------------
    for i, key in ((7, "amine07_acyl"), (8, "amine08_acyl")):
        m = n(c, key)
        A[i] = ("C" * (m - 1) + "C(=O)NN", f"fatty acid hydrazide, C{m} acyl")

    # --- 9: 2-alkyl-1-alkylamine ------------------------------------------
    # Stage 01 reports each branch's carbon count with the shared junction
    # carbon removed AND the label end removed where the branch stops on a
    # heteroatom's text -- the aminomethyl branch does, the two alkyl branches
    # do not. Taking `vertices` here instead would silently add a carbon to
    # the head.
    br = sorted((b["carbons_excluding_junction"]
                 for b in c["amine09_branched"]["branches"]), reverse=True)
    long_arm, short_arm, head = br
    if head != 1:
        raise StructureError(
            f"03_enumerate: amine 9's aminomethyl branch measures {head} "
            "carbon(s), expected 1.")
    A[9] = ("C" * long_arm + "C(" + "C" * short_arm + ")CN",
            f"2-alkyl-1-alkylamine: a {long_arm}C and a {short_arm}C branch "
            "on the carbon next to CH2NH2")

    # --- 10: alkyl glycinate ----------------------------------------------
    g = n(c, "amine10_glycyl")
    if g != 2:
        raise StructureError(
            f"03_enumerate: amine 10's amino-acyl unit measures {g} heavy "
            "atoms, expected 2 (CH2 + C=O, i.e. glycine).")
    if n(c, "amine10_ester_link") != 1:
        raise StructureError(
            "03_enumerate: amine 10's ester O is not bonded straight to the "
            "carbonyl carbon.")
    t10 = n(c, "amine10_ester_alkyl")
    A[10] = ("NCC(=O)O" + "C" * t10, f"glycine C{t10} alkyl ester")

    # --- 11-14: para-substituted anilines ----------------------------------
    for i, key in ((11, "amine11_alkyl"), (12, "amine12_alkyl"),
                   (13, "amine13_alkyl")):
        pat = c[f"amine{i}_ring_substitution"]["pattern"]
        if pat != "para":
            raise StructureError(f"03_enumerate: amine {i} ring is {pat}.")
        A[i] = ("C" * n(c, key) + "c1ccc(N)cc1",
                f"4-alkylaniline, C{n(c, key)} chain, measured {pat}")
    if c["amine14_ring_substitution"]["pattern"] != "para":
        raise StructureError("03_enumerate: amine 14 ring is not para.")
    if c["amine14_ring"]["ring_vertices"] != 6:
        raise StructureError("03_enumerate: amine 14 ring is not 6-membered.")
    t14 = n(c, "amine14_alkoxy")
    A[14] = ("C" * t14 + "Oc1ccc(N)cc1",
             f"4-alkoxyaniline, C{t14} chain, measured para")

    # --- 15-17: alpha,omega-diaminoalkanes ---------------------------------
    for i in (15, 16, 17):
        m = n(c, f"amine{i}_linker")
        A[i] = ("N" + "C" * m + "N", f"alpha,omega-diaminoalkane, C{m}")

    # --- 18: linear triamine with a secondary NH ---------------------------
    a, b = n(c, "amine18_linkerA"), n(c, "amine18_linkerB")
    A[18] = ("N" + "C" * a + "N" + "C" * b + "N",
             f"linear triamine, C{a}/C{b} linkers, central NH secondary")

    # --- 19, 21: N-methyl bis(aminoalkyl)amines ----------------------------
    for i in (19, 21):
        if not c[f"amine{i}_Nmethyl"]["is_methyl"]:
            raise StructureError(f"03_enumerate: amine {i}'s stub is labelled.")
        a, b = n(c, f"amine{i}_linkerA"), n(c, f"amine{i}_linkerB")
        A[i] = ("CN(" + "C" * a + "N)" + "C" * b + "N",
                f"N-methyl bis(aminoalkyl)amine, C{a}/C{b} arms")

    # --- 20, 22: tris(aminoalkyl)amines ------------------------------------
    for i in (20, 22):
        a = n(c, f"amine{i}_linkerA")
        b = n(c, f"amine{i}_linkerB")
        d = n(c, f"amine{i}_linkerC")
        A[i] = ("N" + "C" * a + "N(" + "C" * b + "N)" + "C" * d + "N",
                f"tris(aminoalkyl)amine, C{a}/C{b}/C{d} arms")

    # --- 23, 24: piperazine-cored polyamines -------------------------------
    for i in (23, 24):
        ra, rb = n(c, f"amine{i}_ringA"), n(c, f"amine{i}_ringB")
        if (ra, rb) != (2, 2):
            raise StructureError(
                f"03_enumerate: amine {i}'s ring halves measure {ra}/{rb} "
                "carbons; piperazine needs 2 and 2.")
    a23, b23 = n(c, "amine23_linkerA"), n(c, "amine23_linkerB")
    A[23] = ("N" + "C" * a23 + "N1CCN(" + "C" * b23 + "N)CC1",
             f"1,4-bis(aminoalkyl)piperazine, C{a23}/C{b23} arms")
    a24 = n(c, "amine24_linkerA")
    b24 = n(c, "amine24_linkerB")
    d24 = n(c, "amine24_linkerC")
    A[24] = ("N" + "C" * a24 + "N" + "C" * b24 + "N1CCN(" + "C" * d24 + "N)CC1",
             f"piperazine with a C{a24}/C{b24} aminoalkyl-amino arm and a "
             f"C{d24} aminoalkyl arm; the mid-chain NH is secondary")

    # --- 25: tetra-armed diaminoalkane core --------------------------------
    arms = [n(c, f"amine25_arm{s}") for s in "ABCD"]
    core = n(c, "amine25_core")
    if len(set(arms)) != 1:
        raise StructureError(f"03_enumerate: amine 25's arms measure {arms}.")
    q = arms[0]
    A[25] = ("N" + "C" * q + "N(" + "C" * q + "N)" + "C" * core + "N("
             + "C" * q + "N)" + "C" * q + "N",
             f"N,N,N',N'-tetrakis(aminoalkyl) C{core} diamine; four C{q} arms")

    # --- acrylates ---------------------------------------------------------
    T = {}
    for tok in (9, 12, 14, 16):
        m = n(c, f"acrylate_O{tok}_tail")
        chk = c[f"acrylate_O{tok}_tail"]["name_token_check"]
        if not chk["matches"]:
            raise StructureError(
                f"03_enumerate: acrylate O{tok}'s name token does not equal "
                f"its measured chain length ({m}). {chk['note']}")
        T[f"O{tok}"] = ("C=CC(=O)O" + "C" * m,
                        f"n-alkyl acrylate, C{m} ester", m)
    return A, T


def build_lipid(amine_smiles, tail_alkyl, k, label):
    x = k["linker"]["CH2_between_amidine_C_and_S"]
    y = k["linker"]["CH2_between_S_and_carbonyl"]
    # The amidine carbon is written [CH:1], not [C:1]. A bracket atom carries
    # noImplicit, so [C:1] with a double bond to N and one C neighbour is
    # valence 3 and RDKit gives it an unpaired electron -- the product then
    # parses, canonicalises and has self-consistent descriptors while being a
    # radical. assemble.check_no_radicals caught it here on the first build.
    traut_open = "N=[CH:1]" + "C" * x + "[SH:3]"
    acceptor = "[CH2:2]" + "C" * (y - 1) + "C(=O)O" + tail_alkyl
    arm = assemble.attach(traut_open, acceptor, a_map=3, b_map=2,
                          label=f"{label}:thiol-Michael")

    cur = amine_smiles
    n_arms = len(Chem.MolFromSmiles(cur).GetSubstructMatches(PRIMARY_AMINE))
    if n_arms == 0:
        raise StructureError(f"03_enumerate[{label}]: no primary amine.")
    for step in range(n_arms):
        m = Chem.MolFromSmiles(cur)
        hits = m.GetSubstructMatches(PRIMARY_AMINE)
        if len(hits) != n_arms - step:
            raise StructureError(
                f"03_enumerate[{label}]: after {step} arm(s) the molecule has "
                f"{len(hits)} primary amine(s), expected {n_arms - step}. An "
                "arm was installed somewhere other than a primary nitrogen.")
        for a in m.GetAtoms():
            a.SetAtomMapNum(0)
        m.GetAtomWithIdx(hits[0][0]).SetAtomMapNum(5)
        cur = assemble.attach(Chem.MolToSmiles(m), arm, a_map=5, b_map=1,
                              label=f"{label}:amidine{step + 1}")
    return cur, n_arms, arm


def descriptors(smi):
    m = Chem.MolFromSmiles(smi)
    return {"Formula": rdMolDescriptors.CalcMolFormula(m),
            "MolWt": round(Descriptors.MolWt(m), 3),
            "ExactMolWt": round(Descriptors.ExactMolWt(m), 5),
            "InChIKey": Chem.MolToInchiKey(m)}


def main():
    c, k = load()
    A, T = build_blocks(c, k)

    # ---------------- 01_blocks.csv ---------------------------------------
    rows = []
    for i in sorted(A):
        smi, desc = A[i]
        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise StructureError(f"03_enumerate: amine {i} SMILES {smi!r}.")
        rows.append({"block_type": "amine", "block_id": str(i),
                     "SMILES": Chem.MolToSmiles(m), "description": desc,
                     "n_primary_amines":
                         len(m.GetSubstructMatches(PRIMARY_AMINE)),
                     **descriptors(smi)})
    for tok in ("O9", "O12", "O14", "O16"):
        smi, desc, _m = T[tok]
        m = Chem.MolFromSmiles(smi)
        rows.append({"block_type": "acrylate", "block_id": tok,
                     "SMILES": Chem.MolToSmiles(m), "description": desc,
                     "n_primary_amines": 0, **descriptors(smi)})
    with open(BLOCKS_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {BLOCKS_CSV}: {len(rows)} blocks")

    # ---------------- 03_enumerated.csv -----------------------------------
    recs, order = {}, []
    for i in sorted(A):
        for tok in ("O9", "O12", "O14", "O16"):
            name = f"{i}T-{tok}"
            tail_alkyl = "C" * T[tok][2]
            smi, n_arms, _arm = build_lipid(A[i][0], tail_alkyl, k, name)
            d = descriptors(smi)
            recs[name] = {"lipid_name": name, "amine": i, "tail": tok,
                          "amine_SMILES": A[i][0], "n_arms": n_arms,
                          "IL_SMILES": smi, **d}
            order.append(name)

    if len(order) != len(A) * len(T):
        raise StructureError(
            f"03_enumerate: enumerated {len(order)}, expected "
            f"{len(A)} x {len(T)}.")

    # ---------------- flags, written ON THE ROW ---------------------------
    # Every one of these is derived from the item's own record. emit_rows
    # refuses a literal, which is error class 29 made unwritable.
    SECONDARY_NH = {18, 24}
    MANY_PRIMARY = {20, 22, 25}

    # Amine 19 is the one polyamine whose product the paper DRAWS and whose
    # mass it PRINTS, so it must not carry the same caveat as the amines whose
    # arm count is carried over from it. Labelling the best-evidenced row like
    # the extrapolated ones is the inverted-confidence failure (F-011/F-025)
    # pointing the other way.
    DEMONSTRATED = {19}

    def confidence(item, rec):
        a = rec["amine"]
        if rec["n_arms"] == 1 or a in DEMONSTRATED:
            return "high"
        if a in MANY_PRIMARY:
            return "low"
        return "medium"

    def note(item, rec):
        a, k_ = rec["amine"], rec["n_arms"]
        if k_ == 1:
            return ""
        bits = [f"amine {a} carries {k_} primary nitrogens and is drawn here "
                f"with an arm on each"]
        if a == 19:
            bits.append("which is the substitution pattern Supplementary "
                        "Fig. 1a draws and its mass spectrum confirms "
                        "(401.129 [M+2H]2+)")
        else:
            bits.append("a pattern shown directly only for amine 19 "
                        "(Supplementary Fig. 1a) and carried over here")
        if a in SECONDARY_NH:
            bits.append("its secondary NH is left unsubstituted, which the "
                        "paper nowhere states (DECISIONS.md D-2)")
        if a in MANY_PRIMARY:
            bits.append("and the printed stoichiometry, 3 equiv Traut's "
                        "reagent and 2.5 equiv acrylate per 1 equiv amine, is "
                        "below what exhaustive substitution of its nitrogens "
                        "needs (DECISIONS.md D-3)")
        return "; ".join(bits) + "."

    q = emit.emit_rows(
        order, recs,
        confidence=confidence, note=note,
        structure_assumed=lambda i, r: ("1" if r["amine"] in SECONDARY_NH
                                        or r["amine"] in MANY_PRIMARY else "0"),
        positional_ambiguity=lambda i, r: ("1" if r["amine"] in SECONDARY_NH
                                           else "0"),
        value_is_bound=lambda i, r: "",
        value_normalization_unit=lambda i, r: "",
        label="MW_2024")

    cols = ["lipid_name", "amine", "tail", "n_arms", "amine_SMILES",
            "IL_SMILES", "Formula", "MolWt", "ExactMolWt", "InChIKey"]
    with open(ENUM_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols + emit.columns())
        w.writeheader()
        for name, extra in zip(order, q):
            w.writerow({**{c_: recs[name][c_] for c_ in cols}, **extra})
    print(f"wrote {ENUM_CSV}: {len(order)} lipids, "
          f"{emit.caveat_count(q)} carrying a caveat, "
          f"{len({recs[o]['InChIKey'] for o in order})} distinct InChIKeys")


if __name__ == "__main__":
    main()
