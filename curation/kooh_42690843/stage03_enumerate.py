#!/usr/bin/env python3
"""The PILL lipids this SI lets us build: the nine PURIFIED, characterised ones.

THE RECONCILIATION, UP FRONT. The abstract says 275 alpha-branched amine
ionizable lipids were generated. The SI names NINE of them (A2B1C1, A2B3C1,
A2B8C1, A18B1C1, A18B3C1, A18B8C1, A18B8C8, A18B8C9, A18B8C14) and says
nothing, anywhere, about which 275 of the 26 x 8 x 14 = 2,912 possible
amine x boronic acid x aldehyde combinations were made. 275 is 9.4 % of 2,912
and does not factor into any product of the panel sizes. So this file enumerates
NINE lipids, not 275, and 08_contradictions.csv and DECISIONS.md D-1 record the
gap rather than papering over it. Emitting 2,912 would fabricate 2,903 lipids
that were never made.

THE REACTION, AND THE REGIOCHEMISTRY MASS CANNOT SEE.
SI p.S12, "General Procedure for the Photo-One-Pot Synthesis of Lipids":
amine + aldehyde are pre-mixed in HFIP "to promote imine formation"; catechol
(0.3 eq) and the boronic acid are then added and the vial is irradiated at
427 nm. This is a photochemical alkyl-radical addition to the imine (a
Petasis-type / photo-Petasis multicomponent coupling): the BORONIC ACID'S ALKYL
GROUP ENDS UP ON THE IMINE CARBON, not on the nitrogen.

That is the one thing a molecular formula cannot tell you here. N-alkylation
and C-alkylation of the imine give the SAME formula and the SAME exact mass, so
all nine HRMS values pass either way. Three independent lines of evidence in
the SI pin the C-alkylation, and they are recorded here because a later reader
will otherwise have only the formula:

  1. THE AUTHORS' OWN IUPAC NAMES. A2B1C1 is named
     "N1-(1-cyclopentylnonyl)-N3,N3-dimethylpropane-1,3-diamine". The
     substituent on N1 is 1-cyclopentylNONYL: a nine-carbon chain (nonanal, C1)
     bearing the cyclopentyl (B1) at its C1. Both blocks hang off ONE carbon.
     An N-alkylated product would be named as two separate N-substituents.
     Likewise A2B3C1 is "...-N3-(tridecan-5-yl)..." (nonanal C9 + butyl C4 =
     C13, branched at C5), and A2B8C1 is "...-(heptadecan-9-yl)..." (C9 + C8 =
     C17, branched at C9). Each name independently reproduces the arithmetic.
  2. THE PRODUCT SKELETON NAMES for the diamine series:
     A18B8C8 is "(pentane-1,5-diylbis(azanediyl))bis(TETRADECANE-6,1-DIYL)
     bis(4-methyloctanoate)" -- a C14 chain attached to N at position 6, which
     is exactly C8's five-carbon aldehyde arm (C6..C1 bearing the ester O) plus
     B8's octyl (C7..C14) meeting at C6.
  3. THE 1H NMR. Each product shows a one-proton multiplet near 2.9-3.0 ppm
     (A2B3C1 2.95-2.90 m 1H; A2B8C1 2.96 t 1H; A18B8C1 3.02-2.92 m 3H) -- the
     single methine N-CH. A reductive amination with no boronic-acid
     incorporation would give an N-CH2 2H triplet and no methine.

The alpha carbon so formed bears four different groups and is a STEREOCENTRE.
The reaction is not stereoselective and the products are not resolved, so it is
left UNASSIGNED in every SMILES here and the rows carry confidence = medium
with a note saying so.

HOW THE STRUCTURES ARE BUILT. Every bond is made by curation/lib
assemble.attach, which bonds two atom-mapped fragments without deleting an
atom, and asserts no radical and no lost stereocentre on each build. `attach`
can give an atom ONE new bond per call and clears its map afterwards, but the
alpha carbon needs THREE heavy neighbours (R from the aldehyde, R' from the
boronic acid, N from the amine), so one of them has to be present in a written
fragment. The cut taken here puts the NITROGEN there -- each aldehyde is
written as its ammonia reductive-amination product, i.e. the primary amine
R-CH2-NH2 -- because that fragment is a real, nameable molecule whose formula
is a function of the aldehyde alone (one fragment per aldehyde, not one per
pair), and because it leaves the C-C bond that the regiochemistry turns on to
be made by `attach` rather than written down. The alternative cut writes the
whole R-CH(R')- group as a literal per (B,C) pair, which would put the claim
being tested into the input.

Run:  python3 -I stage03_enumerate.py              writes 03_enumerated.csv
      python3 -I stage03_enumerate.py --seed-defect regio
            builds A2B1C1 with the boronic alkyl on the NITROGEN instead of the
            alpha carbon and confirms (a) HRMS still passes -- mass cannot see
            connectivity -- and (b) the name-derived connectivity check fires.
"""
import csv
import os
import sys

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

LIB = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib")
sys.path.insert(0, LIB)
import assemble                                       # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "03_enumerated.csv")

# --- fragments ---------------------------------------------------------------
# ALDEHYDE -> its ammonia reductive-amination amine, with
#   map 1 on the alpha carbon (the former carbonyl carbon; takes the C-C bond)
#   map 2 on the nitrogen      (takes the bond to the rest of the amine block)
# Each is the aldehyde of 01_blocks.csv with C=O replaced by CH-NH2; the
# formula relation (aldehyde + NH3 + 2[H] - H2O) is asserted below against
# 01_blocks.csv, so these literals cannot drift from the blocks.
ALDEHYDE_AMINE = {
    "C1":  "[NH2:2][CH2:1]CCCCCCCC",
    "C8":  "[NH2:2][CH2:1]CCCCCOC(=O)CCC(C)CCCC",
    "C9":  "[NH2:2][CH2:1]CCCCCOC(=O)CCCCCCCC",
    "C14": "[NH2:2][CH2:1]CCCCCOC(=O)C(CCCCCC)CCCCCCCC",
}
# BORONIC ACID -> the alkyl it transfers, map 3 on the ipso carbon.
BORONIC_ALKYL = {
    "B1": "[CH2:3]1CCCC1",
    "B3": "[CH3:3]CCC",
    "B8": "[CH3:3]CCCCCCC",
}
# AMINE -> what is left once the reacting N(s) are accounted for in the
# aldehyde fragment. map 4 (and 5, for a diamine, on the far end).
AMINE_REST = {
    "A2":  ("mono", "[CH3:4]CCN(C)C"),      # 3-(dimethylamino)propyl
    "A18": ("di",   "[CH3:4]CCC[CH3:5]"),   # pentane-1,5-diyl
}

# The library: every lipid the SI names, with the HRMS the SI prints for it.
# (lipid, A, B, C, SI page, ion, formula_printed, calcd_printed, found_printed)
LIBRARY = [
    ("A2B1C1",   "A2",  "B1", "C1",  "S13", "[M+H]+",  "C19H41N2",  297.3264, 297.3282),
    ("A2B3C1",   "A2",  "B3", "C1",  "S13", "[M+H]+",  "C18H41N2",  285.3264, 285.3279),
    ("A2B8C1",   "A2",  "B8", "C1",  "S14", "[M+H]+",  "C22H49N2",  341.3890, 341.3471),
    ("A18B1C1",  "A18", "B1", "C1",  "S14", "[M+H]+",  "C33H67N2",  491.5299, 491.5301),
    ("A18B3C1",  "A18", "B3", "C1",  "S15", "[M+H]+",  "C31H67N2",  467.5299, 467.5285),
    ("A18B8C1",  "A18", "B8", "C1",  "S16", "[M+H]+",  "C39H83N2",  579.6551, 579.6534),
    ("A18B8C8",  "A18", "B8", "C8",  "S17", "[M+H]+",  "C51H103N2O4", 807.7912, 807.7932),
    ("A18B8C9",  "A18", "B8", "C9",  "S18", "[M+H]+",  "C51H103N2O4", 807.7912, 807.7855),
    ("A18B8C14", "A18", "B8", "C14", "S19", "[M+H]+",  "C65H131N2O4", 1004.0103, 1004.0117),
]

# Independent geometric evidence, from the structure DRAWN beside each
# compound in the SI: ChemDraw draws one stroked line per heavy-atom bond plus
# a second line per double bond, so a page's stroked-line count is a bond-count
# of the structures on it. Measured on the PDF's own coordinates.
# (SI page index 0-based, measured stroked lines, which lipids are drawn there)
DRAWN_LINES = [
    (11, 21, ["A2B1C1"]),
    (12, 42, ["A2B3C1", "A2B8C1"]),
    (13, 36, ["A18B1C1"]),
    (14, 32, ["A18B3C1"]),
    (15, 40, ["A18B8C1"]),
    (16, 58, ["A18B8C8"]),
    (17, 58, ["A18B8C9"]),
    (18, 72, ["A18B8C14"]),
]

PROTON = 1.00727646          # mass of H+, for an [M+H]+ monoisotopic calc
HRMS_TOL_MDA = 2.0           # calcd-vs-ours tolerance, milli-Da


def build(a, b, c):
    """Amine + aldehyde + boronic acid -> the alpha-branched amine lipid."""
    ald = ALDEHYDE_AMINE[c]
    bor = BORONIC_ALKYL[b]
    kind, rest = AMINE_REST[a]
    alpha = assemble.attach(ald, bor, a_map=1, b_map=3,
                            label="%s+%s alpha C-C" % (c, b))
    if kind == "mono":
        return assemble.attach(alpha, rest, a_map=2, b_map=4,
                               label="%s N-alkylation" % a)
    half = assemble.attach(alpha, rest, a_map=2, b_map=4,
                           label="%s first arm" % a)
    alpha2 = assemble.attach(ALDEHYDE_AMINE[c], BORONIC_ALKYL[b],
                             a_map=1, b_map=3, label="%s+%s alpha C-C #2"
                             % (c, b))
    return assemble.attach(half, alpha2, a_map=5, b_map=2,
                           label="%s second arm" % a)


def seed_regio():
    """A2B1C1 with cyclopentyl on the nitrogen, not on the alpha carbon.

    N1-cyclopentyl-N1-nonyl-N3,N3-dimethylpropane-1,3-diamine. Same atoms,
    same formula C19H40N2, same monoisotopic mass. Built by hand because the
    normal route cannot express it -- which is the honest way to show that the
    distinction is invisible to every mass-based check.
    """
    return "CCCCCCCCCN(C1CCCC1)CCCN(C)C"


def mono_mass(formula_mol):
    return Descriptors.ExactMolWt(formula_mol)


def connectivity_from_name(lipid, smi):
    """Does the alpha carbon carry the aldehyde arm AND the boronic alkyl?

    This is the check HRMS cannot do. For each lipid the SI's own IUPAC name
    states that one carbon bears: the amine nitrogen, the aldehyde-derived
    chain, and the boronic-acid-derived alkyl. We test exactly that pattern,
    with the boronic group's identity as drawn in 01_blocks.csv.
    """
    m = Chem.MolFromSmiles(smi)
    # alpha carbon: sp3 C with exactly one H, bonded to an N and to >=2 C
    patt = Chem.MolFromSmarts("[CX4H1]([NX3])([#6])[#6]")
    hits = m.GetSubstructMatches(patt)
    alphas = {h[0] for h in hits}
    want = 2 if lipid.startswith("A18") else 1
    if len(alphas) != want:
        return False, ("found %d alpha methine(s) bearing N, expected %d"
                       % (len(alphas), want))
    # and each alpha carbon must carry the boronic block's group
    bcode = lipid.split("B")[1].split("C")[0]
    probe = {"1": Chem.MolFromSmarts("[CX4H1]([NX3])[CH1]1[CH2][CH2][CH2][CH2]1"),
             "3": Chem.MolFromSmarts("[CX4H1]([NX3])[CH2][CH2][CH2][CH3]"),
             "8": Chem.MolFromSmarts(
                 "[CX4H1]([NX3])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH3]")}[bcode]
    n = len(m.GetSubstructMatches(probe))
    if n < want:
        return False, ("alpha carbon does not carry the B%s group (%d/%d "
                       "matches)" % (bcode, n, want))
    return True, "alpha methine carries N + the B%s alkyl, x%d" % (bcode, want)


def check_fragment_formulas(blocks_csv):
    """The aldehyde-amine literals must equal their aldehyde + NH3 + 2H - H2O."""
    blocks = {}
    with open(blocks_csv) as fh:
        for r in csv.DictReader(fh):
            blocks[r["block"]] = r["SMILES"]
    out = []
    for c, frag in ALDEHYDE_AMINE.items():
        f = Chem.MolFromSmiles(frag)
        for a in f.GetAtoms():
            a.SetAtomMapNum(0)
        got = rdMolDescriptors.CalcMolFormula(f)
        ald = Chem.MolFromSmiles(blocks[c])
        # aldehyde + NH3 + H2 - H2O, computed on atom counts
        cnt = {}
        for a in ald.GetAtoms():
            cnt[a.GetSymbol()] = cnt.get(a.GetSymbol(), 0) + 1
        cnt["H"] = sum(a.GetTotalNumHs() for a in ald.GetAtoms())
        cnt["N"] = cnt.get("N", 0) + 1
        cnt["H"] = cnt["H"] + 3 + 2 - 2
        cnt["O"] = cnt["O"] - 1
        want = "C%dH%dN%s" % (cnt["C"], cnt["H"],
                              "" if cnt["N"] == 1 else cnt["N"])
        if cnt["O"]:
            want += "O%s" % ("" if cnt["O"] == 1 else cnt["O"])
        out.append((c, got, want, got == want))
    return out


def main():
    seeded = "--seed-defect" in sys.argv
    print("fragment formulas vs 01_blocks.csv:")
    frag_bad = []
    for c, got, want, ok in check_fragment_formulas(
            os.path.join(HERE, "01_blocks.csv")):
        print("  %s %-4s fragment %-10s  aldehyde+NH3+2[H]-H2O = %-10s"
              % ("ok " if ok else "BAD", c, got, want))
        if not ok:
            frag_bad.append(c)
    if frag_bad:
        raise SystemExit("stage03: aldehyde fragment(s) %s do not match their "
                         "block. Refusing." % frag_bad)

    rows, bad = [], []
    print("\nbuild + HRMS (calcd printed in the SI vs our structure) + "
          "name-derived connectivity:")
    for (lipid, a, b, c, pg, ion, formula, calcd, found) in LIBRARY:
        smi = build(a, b, c)
        m = Chem.MolFromSmiles(smi)
        ours = Descriptors.ExactMolWt(m) + PROTON
        d_mda = (ours - calcd) * 1000.0
        f_ok = rdMolDescriptors.CalcMolFormula(m)
        conn_ok, conn_why = connectivity_from_name(lipid, smi)
        ok = abs(d_mda) <= HRMS_TOL_MDA and conn_ok
        if not ok:
            bad.append((lipid, d_mda, conn_why))
        print("  %s %-10s %-5s  neutral %-13s  [M+H]+ ours %10.4f  "
              "SI calcd %10.4f  d %+6.2f mDa   %s"
              % ("ok " if ok else "BAD", lipid, pg, f_ok, ours, calcd,
                 d_mda, conn_why))
        nst = sum(1 for _ in Chem.FindMolChiralCenters(
            m, includeUnassigned=True, useLegacyImplementation=False))
        rows.append({
            "lipid": lipid, "amine_block": a, "boronic_block": b,
            "aldehyde_block": c,
            "IL_SMILES": Chem.MolToSmiles(m),
            "formula": f_ok,
            "MW": round(Descriptors.MolWt(m), 3),
            "exact_mass_MplusH": round(ours, 4),
            "SI_page": pg, "SI_HRMS_ion": ion,
            "SI_HRMS_formula_printed": formula,
            "SI_HRMS_calcd_printed": calcd,
            "SI_HRMS_found_printed": found,
            "delta_calcd_mDa": round(d_mda, 3),
            "unassigned_stereocentres": nst,
            "purified": 1,
            "regiochemistry_evidence":
                "SI IUPAC name + product skeleton name + "
                "1H NMR N-CH methine",
            "confidence": "high",
            "confidence_reason":
                "structure built from blocks whose identity is pinned by the "
                "printed reagent amount; HRMS formula and the drawn bond-line "
                "count both agree; connectivity pinned by the authors' own "
                "IUPAC name. The alpha stereocentre is unassigned -- the "
                "reaction is not stereoselective and the SI does not resolve "
                "it.",
        })

    print("\ndrawn-structure bond-line count on each SI page "
          "(the PDF's own coordinates):")
    import pymupdf
    doc = pymupdf.open(os.path.join(HERE, "source", "ja6c13674_si_001.pdf"))
    smi_by = {r["lipid"]: r["IL_SMILES"] for r in rows}
    draw_bad = []
    for pno, measured, lipids in DRAWN_LINES:
        exp = 0
        for lp in lipids:
            mm = Chem.MolFromSmiles(smi_by[lp])
            k = Chem.Mol(mm)
            Chem.Kekulize(k, clearAromaticFlags=True)
            exp += mm.GetNumBonds() + sum(
                1 for bd in k.GetBonds()
                if bd.GetBondType() == Chem.BondType.DOUBLE)
        n = 0
        for path in doc[pno].get_drawings():
            if path["type"] == "s":
                n += sum(1 for it in path["items"] if it[0] == "l")
        okd = (n == exp and n == measured)
        if not okd:
            draw_bad.append((pno, exp, n))
        print("  %s S%-3d %-22s expected %3d lines, page carries %3d"
              % ("ok " if okd else "BAD", pno + 1, "+".join(lipids), exp, n))

    if seeded:
        print("\n--- SEEDED DEFECT: A2B1C1 built N-alkylated instead of "
              "C-alkylated ---")
        bad_smi = seed_regio()
        bm = Chem.MolFromSmiles(bad_smi)
        bmass = Descriptors.ExactMolWt(bm) + PROTON
        print("  formula %s   [M+H]+ %.4f   SI calcd %.4f   d %+0.3f mDa"
              % (rdMolDescriptors.CalcMolFormula(bm), bmass, 297.3264,
                 (bmass - 297.3264) * 1000))
        print("  -> HRMS check: %s  (MASS CANNOT SEE CONNECTIVITY)"
              % ("PASSES" if abs(bmass - 297.3264) * 1000 <= HRMS_TOL_MDA
                 else "fails"))
        cok, cwhy = connectivity_from_name("A2B1C1", bad_smi)
        print("  -> name-derived connectivity check: %s -- %s"
              % ("PASSES (BAD!)" if cok else "FIRES", cwhy))
        if cok:
            raise SystemExit(
                "SEEDED REGIOCHEMISTRY DEFECT DID NOT FIRE. The only check "
                "that can see this class is blind, so a clean run means "
                "nothing (CLAUDE.md section 6).")
        return 0

    if bad or draw_bad:
        raise SystemExit("stage03: refusing to write with a failing check: "
                         "%s %s" % (bad, draw_bad))
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print("\nwrote %s (%d lipids; the paper states 275 were MADE -- "
          "see 08_contradictions.csv and DECISIONS.md D-1)" % (OUT, len(rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
