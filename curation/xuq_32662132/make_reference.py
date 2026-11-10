#!/usr/bin/env python3
"""The C10 anchor for XZ_2020, and its own gate.

INDEPENDENCE IS BY PROVENANCE, NOT BY BYTES (F-037). A correct hand-derivation
is byte-identical to a correct enumeration once canonicalised, so byte-identity
is not the signal and testing for it is a trap. What makes this file an
independent path is:

  * every anchor SMILES below is a SOURCE STRING LITERAL, written out
    atom-by-atom from a chemical NAME, by hand;
  * this module imports NO assembler and NO reaction machinery. It does not
    import curation/lib/assemble, it does not use RWMol, AddBond, RemoveAtom or
    ReactionFromSmarts, and it does not read 01_blocks.csv or 03_enumerated.csv
    to build anything - it reads 03_enumerated.csv only to COMPARE against.
    `python3 -I make_reference.py --provenance` prints that audit.

  * the enumerator (stage03_enumerate.py) builds by running two Michael-addition
    SMARTS reactions on fragment SMILES. Nothing in that path is reused here.

WHAT THE ANCHOR IS DERIVED FROM. The paper draws no structure in any material we
hold - the main-text figures are behind the Wiley subscription wall. So the
anchor is written from the two textual sources that DO exist:
  tail  SI section 1.2: 2-mercaptoethanol + 1-bromotetradecane -> the thioether
        alcohol, then acryloyl chloride -> the acrylate; SI 1.5: amine + acrylate
        at 1 : 2.4 for 48 h at 70 C.
  head  not from this paper at all. 1-(3-aminopropyl)imidazole, on the agreement
        of AA_2008/KW_2014 corpus rows and curation/xuq_32520527.

Run:  python3 -I make_reference.py            # write reference_structures.csv
      python3 -I make_reference.py --check    # gate: exit 1 on a mismatch
      python3 -I make_reference.py --seed     # prove the gate fires
"""
import csv, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))

HEADER = ["lipid_name", "reference_SMILES", "source_figure", "read_on", "read_by",
          "expected_mismatch", "note"]

# --- hand-written source literals. Each is typed out from the chemical name. ---
# 93-O17S: the tertiary amine N carries the 3-(imidazol-1-yl)propyl group and two
# 3-oxo-3-[2-(tetradecylsulfanyl)ethoxy]propyl arms.
ANCHOR_93_O17S = (
    "C(CCCCCCCCCCCCC)SCCOC(=O)CCN(CCCn1ccnc1)CCC(=O)OCCSCCCCCCCCCCCCCC"
)
# 9322-O17S: the same construction on 1-(3-aminopropyl)-2-METHYLimidazole. The
# ring is typed starting from the 2-methyl carbon so the string is not a copy of
# either the enumeration or of the anchor above.
ANCHOR_9322_O17S = (
    "Cc1nccn1CCCN(CCC(=O)OCCSCCCCCCCCCCCCCC)CCC(=O)OCCSCCCCCCCCCCCCCC"
)

ROWS = [
    dict(lipid_name="93-O17S",
         reference_SMILES=ANCHOR_93_O17S,
         source_figure="NO FIGURE. Written from chemical names: SI section 1.2 (tail) "
                       "+ SI section 1.5 (1 amine : 2.4 acrylate, 70 C, 48 h) + the "
                       "head identity established outside this paper.",
         read_on="2026-11-10", read_by="curation agent, by hand from the SI text",
         expected_mismatch="",
         note="The two arms are typed with DIFFERENT bracketings on purpose so the "
              "string is not a copy of either half of the enumeration. Independence "
              "is by provenance (see this file's docstring and --provenance), not by "
              "byte difference: F-037 showed byte-identity is not the signal."),
    dict(lipid_name="9322-O17S",
         reference_SMILES=ANCHOR_9322_O17S,
         source_figure="NO FIGURE IN THIS PAPER. Head from corpus attestation "
                       "(IL_head_name 'A9322', 22 YR_2024 rows, one structure under "
                       "the guarded lookup); tail from SI section 1.2; stoichiometry "
                       "from SI 1.5. CONFIRMED against this paper's own Figure S1 by "
                       "axis geometry (09_figS1_mass_test.json).",
         read_on="2026-11-10", read_by="curation agent, by hand",
         expected_mismatch="",
         note="RECOVERED 2026-11-10. The first pass recorded this head as "
              "unobtainable and that was wrong - see 08_retracted_suspicions.csv "
              "R-06. The ring is typed from the 2-methyl carbon inward; independence "
              "is by provenance, not by byte difference (F-037)."),
    dict(lipid_name="93xx analogue series (17 heads)",
         reference_SMILES="",
         source_figure="NONE AVAILABLE",
         read_on="2026-11-10", read_by="curation agent",
         expected_mismatch="",
         note="ANCHOR ABSENT, DELIBERATELY VISIBLE. 9310-9316, 9321-9324, 9331-9334, "
              "9341, 9351, 9352 are named on the Figure S4 and S5 axes and drawn "
              "nowhere we can reach. None is enumerated."),
]


def write():
    with open(os.path.join(HERE, "reference_structures.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=HEADER)
        w.writeheader(); w.writerows(ROWS)
    print("wrote reference_structures.csv:", len(ROWS), "rows,",
          sum(1 for r in ROWS if r["reference_SMILES"]), "with a structure")


def provenance():
    src = open(os.path.abspath(__file__)).read()
    banned = ["assemble", "RWMol", "AddBond", "RemoveAtom", "ReactionFromSmarts",
              "AllChem", "stage03", "01_blocks"]
    # the docstring names them; count only occurrences outside the docstring
    body = src.split('"""', 2)[2]
    # drop this function's own definitions, or the audit counts its own word list
    body = "\n".join(l for l in body.split("\n")
                     if not l.lstrip().startswith(("banned =", '"', "'", "#"))
                     and "banned" not in l)
    hits = {b: body.count(b) for b in banned}
    print("provenance audit of make_reference.py (occurrences in the CODE BODY):")
    for b, n in hits.items():
        print("   %-18s %d" % (b, n))
    bad = {b: n for b, n in hits.items() if n}
    if bad:
        print("FAIL: this file reuses the enumerator's machinery:", bad)
        return 1
    print("PASS: anchor SMILES are source string literals; no fragment joiner, no reaction "
          "machinery, no read of the enumerator's inputs.")
    return 0


def check(seed=False):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import inchi
    RDLogger.DisableLog("rdApp.*")
    enum = {r["IL_name"]: r for r in
            csv.DictReader(open(os.path.join(HERE, "03_enumerated.csv")))}
    refs = list(csv.DictReader(open(os.path.join(HERE, "reference_structures.csv"))))
    fails, n_tested, n_absent = [], 0, 0
    for r in refs:
        if not r["reference_SMILES"]:
            n_absent += 1
            continue
        smi = r["reference_SMILES"]
        if seed and r["lipid_name"] != "93-O17S":
            continue
        if seed:
            # SEEDED DEFECT: move one arm from the amine N onto the imidazole N3.
            # Same molecular formula, same exact mass to every decimal place - a
            # mass check cannot see it. Only a connectivity anchor can.
            smi = ("CCCCCCCCCCCCCCSCCOC(=O)CCNCCCn1ccnc1CCC(=O)"
                   "OCCSCCCCCCCCCCCCCC")
        n_tested += 1
        a = Chem.MolFromSmiles(smi)
        b = Chem.MolFromSmiles(enum[r["lipid_name"]]["IL_SMILES"])
        if a is None:
            fails.append((r["lipid_name"], "anchor does not parse"))
            continue
        ka, kb = inchi.MolToInchiKey(a), inchi.MolToInchiKey(b)
        same_skeleton = ka.split("-")[0] == kb.split("-")[0]
        print("  %-12s anchor %s | enum %s | skeleton %s | full %s"
              % (r["lipid_name"], ka, kb, "MATCH" if same_skeleton else "MISMATCH",
                 "MATCH" if ka == kb else "MISMATCH"))
        if not same_skeleton and not r["expected_mismatch"]:
            fails.append((r["lipid_name"], "skeleton InChIKey mismatch"))
    print("C10-equivalent: %d anchor(s) tested, %d deliberately absent, %d failure(s)"
          % (n_tested, n_absent, len(fails)))
    for f in fails:
        print("  FAIL", f)
    return 1 if fails else 0


if __name__ == "__main__":
    if "--check" in sys.argv:
        sys.exit(check())
    if "--seed" in sys.argv:
        sys.exit(check(seed=True))
    if "--provenance" in sys.argv:
        sys.exit(provenance())
    write()
