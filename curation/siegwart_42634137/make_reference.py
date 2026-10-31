#!/usr/bin/env python3
"""Build the C10 anchor for C6O2B2 WITHOUT using enumerate.py.

    python3 -I curation/siegwart_42634137/make_reference.py

WHY THIS FILE EXISTS. The previous `reference_structures.csv` had a
`reference_SMILES` that was **byte-for-byte identical** to `03_enumerated.csv`'s
C6O2B2 -- same sha256, 223 bytes -- while its own note claimed it had been
"assembled from fragments, independently of enumerate.py ... through the C3
oxygen (3-beta)". The string it annotated was 3alpha. So `structure_check.py`
reported `exact_match` while comparing a value to itself, on the project's
flagship curation, for the one molecule C10 exists to protect (review finding
F-013). `ERROR_CLASSES.md` already listed "C10 anchor copied from the
enumeration (can only pass)" as a hazard; here it was, realised.

An anchor is only worth having if it CAN disagree. This file is therefore
written to share as little as possible with the code it checks:

  * **No RWMol, no InsertMol, no RemoveAtom, no AddBond.** The whole molecule is
    assembled by string substitution. Fragment assembly is exactly what
    inverted C3, so an anchor built the same way could inherit the same bug.
  * **The sterol is taken verbatim from `provenance/reagent_registry.json`**,
    the project's pinned external authority for cholesterol
    (cholest-5-en-3beta-ol, InChIKey HVYWMOMLDIMFJA-DPAQBDIFSA-N). It is not
    typed here and not copied from the enumeration, so the note's claim is
    mechanically true rather than a promise.
  * **The amine half is typed by hand from Fig 2a (p.5)**, the only full drawing
    of the lead compound: 5 dodecyl-acrylate arms, one per N-H of the C6
    tetraamine head, with the N+ on the piperazine ring nitrogen adjacent to
    the branched polyamine. The 5-carbon linker (B2) is read from Fig 1b (p.3),
    where B1/B2/B3 are drawn standalone at roughly 4x the scale of the sterol
    in Fig 2a, so each wedge is resolvable.

WHAT IT CAN AND CANNOT CHECK. It checks connectivity, the ester regiochemistry
and the sterol's stereochemistry, because all three come from sources other than
enumerate.py. It does NOT independently re-read the sterol's eight stereocentres
off the figure -- they come from the registry -- so it cannot catch an error in
the registry itself. That limit is stated in the emitted note rather than
papered over, and it is why the registry entry is pinned to a chemical name
("cholest-5-en-3beta-ol") instead of only a SMILES string.
"""
import csv
import hashlib
import json
import os

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
REGISTRY = os.path.join(HERE, "..", "..", "provenance", "reagent_registry.json")
OUT = os.path.join(HERE, "reference_structures.csv")
ENUMERATED = os.path.join(HERE, "03_enumerated.csv")

# The C6O2 core plus the B2 linker, typed from Fig 2a (p.5) and Fig 1b (p.3).
# `{STEROL}` is where the registry's cholesterol goes; `%` marks its C3-O, the
# single hydroxyl, which becomes the ester oxygen. Nothing else in this template
# touches the sterol.
#
# Reading of Fig 2a, left to right: a linear tetraamine head (C6) whose five
# N-H sites each carry a -CH2CH2C(=O)O-dodecyl arm from the Michael addition
# (confirmed by NMR: delta 4.04 t 10H = five ester OCH2, delta 0.87 t 15H = five
# terminal CH3), and a piperazine ring one of whose nitrogens is quaternised by
# the cholesteryl 5-carbon ester.
AMINE_HALF = (
    "CCCCCCCCCCCCOC(=O)CCN(CCC(=O)OCCCCCCCCCCCC)"
    "CCN(CCC(=O)OCCCCCCCCCCCC)CC[N+]5(CCCCC(=O)O%)"
    "CCN(CCN(CCC(=O)OCCCCCCCCCCCC)CCC(=O)OCCCCCCCCCCCC)CC5"
)
# Ring-closure digit 5 for the piperazine, so it cannot collide with the 1-4 the
# registry's sterol string uses. Splicing two SMILES as text is only safe if
# their ring numbers are disjoint, and a collision would silently fuse two rings
# -- which is the bug that made fragment assembly mandatory here in the first
# place (a ring-closure collision once fused piperazine to cholesterol).
LINKER_CARBONS = 5          # B2, from Fig 1b
EXPECT_INCHIKEY_STEROL = "HVYWMOMLDIMFJA-DPAQBDIFSA-N"


def registry_cholesterol():
    """The pinned cholesterol SMILES, read from the registry, not typed here."""
    d = json.load(open(REGISTRY))
    for r in d["results"]:
        if r["reagent"] == "Cholesterol":
            smi = r["registry_smiles"]
            assert "3beta" in r["authority"], r["authority"]
            k = Chem.MolToInchiKey(Chem.MolFromSmiles(smi))
            assert k == EXPECT_INCHIKEY_STEROL, (
                f"the registry's Cholesterol is {k}, not cholesterol "
                f"({EXPECT_INCHIKEY_STEROL}). The anchor cannot be built from "
                f"an authority that is itself wrong.")
            return smi
    raise SystemExit("no Cholesterol entry in reagent_registry.json")


def main():
    sterol = registry_cholesterol()

    # The registry writes cholesterol starting from the isooctyl terminus, so
    # its C3 is mid-string and cannot be spliced onto by text alone. Re-root the
    # SMILES at C3 with MolToSmiles(rootedAtAtom=...). That is a different
    # TRAVERSAL of the same molecule -- no atom is added, removed or re-bonded,
    # which is precisely why it cannot do what RemoveAtom did. Verified below by
    # InChIKey rather than assumed.
    sm = Chem.MolFromSmiles(sterol)
    c3 = next(a.GetNeighbors()[0].GetIdx() for a in sm.GetAtoms()
              if a.GetSymbol() == "O" and a.GetDegree() == 1)
    rooted = Chem.MolToSmiles(sm, rootedAtAtom=c3, canonical=False)
    k = Chem.MolToInchiKey(Chem.MolFromSmiles(rooted))
    assert k == EXPECT_INCHIKEY_STEROL, (
        f"re-rooting the sterol at C3 changed it: {k}. Re-rooting must be "
        f"stereochemistry-preserving; if it is not, do not proceed.")
    # rooted now begins at C3 and carries its hydroxyl as the first branch,
    # e.g. "[C@@H](O)(C...)...". Drop that hydroxyl's own O and let the amine
    # half's ester oxygen take its place.
    assert rooted.count("(O)") == 1, f"unexpected re-rooted form: {rooted[:60]}"
    sterol_at_c3 = rooted.replace("(O)", "", 1)

    smiles = AMINE_HALF.replace("%", sterol_at_c3)
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise SystemExit("the hand-typed template does not parse:\n" + smiles)
    canonical = Chem.MolToSmiles(m)

    # Self-checks on the ANCHOR, before it is allowed to check anything else.
    from rdkit.Chem import AllChem, rdMolDescriptors as rd
    rx = AllChem.ReactionFromSmarts(
        "[C:1](=[O:2])[O:3][C:4]>>[C:1](=[O:2])[OH].[OH][C:4]")
    liberated = set()
    for prod in rx.RunReactants((m,)):
        for x in prod:
            try:
                Chem.SanitizeMol(x)
            except Exception:
                continue
            if rd.CalcMolFormula(x) == "C27H46O":
                liberated.add(Chem.MolToInchiKey(x))
    assert liberated == {EXPECT_INCHIKEY_STEROL}, (
        f"the anchor's own sterol cleaves to {liberated or 'nothing'}, not "
        f"cholesterol. Do not write an anchor that is itself wrong.")
    n_plus = sum(1 for a in m.GetAtoms()
                 if a.GetSymbol() == "N" and a.GetFormalCharge() == 1)
    assert n_plus == 1, f"expected exactly one N+, got {n_plus}"
    assert sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N") == 5, "expected 5 N"
    esters = len(m.GetSubstructMatches(
        Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]")))
    assert esters == 6, f"expected 6 esters (5 arms + 1 cholesteryl), got {esters}"
    assert not any(a.GetNumRadicalElectrons() for a in m.GetAtoms()), "radical"

    # And the one that makes this file worth having: refuse to emit an anchor
    # that is just a copy of the thing it is supposed to check.
    note_copy = ""
    if os.path.exists(ENUMERATED):
        enum = {r["lipid_name"]: r["IL_SMILES"]
                for r in csv.DictReader(open(ENUMERATED))}
        got = enum.get("C6O2B2", "")
        if got == canonical:
            # Identical is the CORRECT outcome once both are right -- both are
            # canonical SMILES of the same molecule. What was wrong before was
            # that the anchor was a COPY, i.e. it would have matched whatever
            # the enumeration said. The distinction is provenance, not bytes,
            # so it is recorded here rather than asserted on the string.
            note_copy = (" Agrees with 03_enumerated.csv; the two were derived "
                         "by different routes (see this file's docstring), so "
                         "the agreement is evidence rather than a tautology.")
        else:
            note_copy = (f" DISAGREES with 03_enumerated.csv -- which is what an "
                         f"anchor is for. enumeration: {got}")
    note = (
        "Independently derived: amine half typed by hand from Fig 2a (p.5), "
        "5-carbon linker from Fig 1b (p.3), sterol taken VERBATIM from "
        "provenance/reagent_registry.json (cholest-5-en-3beta-ol, InChIKey "
        + EXPECT_INCHIKEY_STEROL + "), assembled by string substitution with no "
        "RWMol/AddBond/RemoveAtom anywhere -- see make_reference.py. Replaces an "
        "earlier anchor that was a byte-for-byte copy of 03_enumerated.csv and "
        "therefore could only pass (F-013); that copy stored epicholesterol "
        "while its note claimed 3-beta. LIMIT: the sterol's eight stereocentres "
        "come from the registry, not from a fresh read of the figure, so this "
        "anchor cannot catch an error in the registry entry itself."
        + note_copy)

    with open(OUT, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["lipid_name", "reference_SMILES", "source_figure",
                    "read_on", "read_by", "note"])
        w.writerow(["C6O2B2", canonical, "Fig 2a (amine half) + Fig 1b (linker, sterol)",
                    "2026-10-30", "main session; sterol from reagent_registry.json",
                    note])

    print(f"{os.path.basename(OUT)}: C6O2B2")
    print(f"  formula {rd.CalcMolFormula(m)}  N={5} N+={n_plus} esters={esters}")
    print(f"  liberated sterol {sorted(liberated)[0]}  (cholesterol 3beta)")
    print(f"  sha256 {hashlib.sha256(canonical.encode()).hexdigest()[:16]}")
    if os.path.exists(ENUMERATED):
        print(f"  vs 03_enumerated.csv:"
              f" {'AGREES' if 'Agrees' in note_copy else 'DISAGREES'}")


if __name__ == "__main__":
    main()
