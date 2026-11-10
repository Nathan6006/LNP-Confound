#!/usr/bin/env python3
"""ONE audited fragment assembler, with a stereo-aware post-build assertion.

WHY THIS FILE EXISTS. Both assembly routes this project has used produced a
structural error that nothing downstream could see:

  STRING CONCATENATION  fused piperazine to cholesterol through a ring-closure
                        digit collision. Caught only because building the C10
                        anchor forced someone to look at the drawn structure.
  FRAGMENT ASSEMBLY     the mandated replacement then inverted C3 on all 537
                        `ZT_2026` rows (F-013). `RWMol.RemoveAtom` drops an
                        atom out of a stereocentre's neighbour list WITHOUT
                        correcting that centre's chiral tag, so cholesterol
                        became epicholesterol. The epimers share molecular
                        formula, exact mass to four decimals, the InChIKey
                        skeleton block and all 17 descriptor columns. C1-C9
                        were blind; C10 was blind because its anchor was a byte
                        copy of the enumeration; C11 was blind twice.

So CLAUDE.md section 6 is right that "fragment assembly is not the safe option
-- it is a different hazard", and the two rules it draws are implemented here as
code rather than as advice:

  1. WHERE YOU MUST DELETE AN ATOM, DELETE ONE NOT ADJACENT TO A STEREOCENTRE.
     For an esterification that means the ACID loses its OH and the ALCOHOL's
     oxygen is retained -- which is what the chemistry does anyway. `esterify`
     does it that way and `attach` refuses any deletion adjacent to a centre.
  2. NEVER BOND TO "ATOM INDEX 0 OF THE FRAGMENT." Tag the attachment atom with
     an atom map number; `attach` resolves by map, not by index.

AND THE ASSERTION THAT WOULD HAVE CAUGHT F-013. `check_stereo_preserved`
compares the CIP-labelled stereocentres of each input fragment against the
product. A centre that was R and is now S, or that was assigned and is now
unassigned, raises StructureError. This is cheap, it runs on every build, and
no composition, mass or descriptor check can substitute for it.

An external anchor is still required -- this module checks that assembly did
not CHANGE a fragment's stereochemistry, not that the fragment was right. That
is C10's job (provenance/structure_check.py) and C16's (the anchor must be
independent; see check_anchor_independence.py and F-037).
"""
from rdkit import Chem
from rdkit.Chem import AllChem  # noqa: F401  (registers CIP code paths)

from errors import StructureError

# Silence RDKit's own chatter so a real refusal is not lost in a log.
try:
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
except Exception:                                      # pragma: no cover
    pass


def _mol(x, what="input"):
    if isinstance(x, Chem.Mol):
        return Chem.Mol(x)
    m = Chem.MolFromSmiles(x)
    if m is None:
        raise StructureError(f"assemble: {what} is not parseable SMILES: {x!r}")
    return m


def stereocentres(mol):
    """{atom_idx: 'R'|'S'|'?'} for every potential stereocentre.

    includeUnassigned=True on purpose: F-013's signature is a centre going from
    assigned to unassigned, and a function that only looks at assigned centres
    cannot see that.
    """
    m = _mol(mol, "molecule")
    Chem.AssignStereochemistry(m, cleanIt=True, force=True,
                               flagPossibleStereoCenters=True)
    out = {}
    for a in m.GetAtoms():
        if a.HasProp("_CIPCode"):
            out[a.GetIdx()] = a.GetProp("_CIPCode")
        elif a.HasProp("_ChiralityPossible"):
            out[a.GetIdx()] = "?"
    return out


def stereo_fingerprint(mol):
    """A canonical multiset of CIP codes, independent of atom numbering.

    Atom indices change during assembly, so an index-keyed comparison would be
    meaningless. This counts codes instead: {'R': 2, 'S': 1, '?': 0}. A centre
    inverting shows up as R and S swapping; a centre being lost shows up as a
    count falling or '?' rising.
    """
    codes = stereocentres(_mol(mol, "molecule")).values()
    fp = {"R": 0, "S": 0, "?": 0}
    for c in codes:
        fp[c] = fp.get(c, 0) + 1
    return fp


def check_stereo_preserved(product, fragments, label="",
                           allow_cip_relabel=False):
    """Raise unless the product's stereocentres are the fragments' sum.

    This is the assertion that F-013 did not have. Returns the comparison dict
    for the QC file on success.
    """
    prod = stereo_fingerprint(_mol(product, "product"))
    want = {"R": 0, "S": 0, "?": 0}
    for f in fragments:
        for k, v in stereo_fingerprint(_mol(f, "fragment")).items():
            want[k] = want.get(k, 0) + v
    # Refuse only a LOSS or an INVERSION of an assigned centre. A bond that
    # CREATES a new stereocentre is ordinary chemistry -- alkylating the second
    # nitrogen of a cis-fused bicycle, or a ring carbon bearing OH, makes a
    # centre that did not exist in either fragment. The original rule compared
    # the whole fingerprint for equality and so refused those builds too; it
    # fired on 40 legitimate assemblies in andersondg_41203968 (2026-11-10).
    #
    # F-013 is still caught: inverting C3 took cholesterol's R:4,S:4 to R:5,S:3,
    # and S DECREASED. An assigned centre going missing or flipping always shows
    # up as a decrease in R or in S.
    lost = {k: want[k] - prod.get(k, 0) for k in ("R", "S")
            if prod.get(k, 0) < want[k]}
    gained = {k: prod.get(k, 0) - want[k] for k in ("R", "S", "?")
              if prod.get(k, 0) > want.get(k, 0)}
    if lost and allow_cip_relabel:
        # A CIP RELABEL is not an inversion. Fall back to a parity comparison,
        # which is numbering- and substituent-independent. Opt-in until the
        # running extraction batch drains (editing a shared check mid-batch is
        # the hazard that bit this run on 2026-11-10).
        ok, detail = parity_preserved(product, fragments[-1], label=label)
        if ok and detail.get("mapped"):
            return {"product": prod, "fragments_sum": want, "preserved": True,
                    "cip_relabel": lost, "parity": detail}
    if lost:
        raise StructureError(
            f"assemble.check_stereo_preserved({label!r}): the product's "
            f"stereocentres are {prod}, the fragments' sum is {want}.\n"
            "Assembly changed a stereocentre. The known mechanism is "
            "RWMol.RemoveAtom dropping an atom out of a centre's neighbour "
            "list without correcting its chiral tag -- that turned cholesterol "
            "into epicholesterol on 537 rows (F-013), invisible to formula, "
            "exact mass to 4 dp, the InChIKey skeleton block and all 17 "
            "descriptor columns.\n"
            "Delete an atom that is NOT adjacent to a stereocentre. For an "
            "esterification the acid loses its OH and the alcohol's oxygen is "
            "retained, which is what the chemistry does.\n"
            f"Assigned centres LOST or INVERTED: {lost}.")
    return {"product": prod, "fragments_sum": want, "preserved": True,
            "centres_created": gained or None}


def _mapped_atom(mol, mapno, what):
    hits = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomMapNum() == mapno]
    if len(hits) != 1:
        raise StructureError(
            f"assemble: {what} has {len(hits)} atoms with map number {mapno}, "
            "need exactly 1. Tag the attachment atom with an atom map "
            "(e.g. 'OCC[N:1]CC'). Bonding to 'atom index 0 of the fragment' is "
            "forbidden here: canonicalisation renumbers atoms, so index 0 is "
            "not a chemical statement. CLAUDE.md section 6.")
    return hits[0]


def attach(a_smiles, b_smiles, a_map=1, b_map=2, bond="single", label=""):
    """Bond two fragments at their atom-mapped attachment atoms.

    Neither fragment loses an atom, so there is no RemoveAtom and no chance of
    the F-013 mechanism. Use this for an N-alkylation, a C-C coupling, or any
    join where both partners keep every atom. For a condensation that must
    expel a leaving group, use `esterify` or `condense`.

    Returns canonical SMILES. Asserts stereo preservation before returning.
    """
    A, B = _mol(a_smiles, "fragment A"), _mol(b_smiles, "fragment B")
    ia = _mapped_atom(A, a_map, "fragment A")
    ib = _mapped_atom(B, b_map, "fragment B")

    combo = Chem.RWMol(Chem.CombineMols(A, B))
    ib_off = ib + A.GetNumAtoms()
    order = {"single": Chem.BondType.SINGLE, "double": Chem.BondType.DOUBLE,
             "triple": Chem.BondType.TRIPLE}
    if bond not in order:
        raise StructureError(f"assemble.attach: bond={bond!r}.")
    if combo.GetBondBetweenAtoms(ia, ib_off) is not None:
        raise StructureError("assemble.attach: those atoms are already bonded.")
    combo.AddBond(ia, ib_off, order[bond])
    for idx in (ia, ib_off):
        at = combo.GetAtomWithIdx(idx)
        at.SetAtomMapNum(0)
        # Hand the hydrogen count back to RDKit rather than decrementing it.
        # A bracket atom in SMILES carries noImplicit=True, so decrementing
        # '[NH:1]' to 0 explicit H gives '[N]' -- a RADICAL, not a secondary
        # amine. That is the Siegwart radical-carbon class: it parses,
        # canonicalises, and every descriptor is self-consistent. Clearing both
        # flags lets valence decide, which is what the chemistry does.
        at.SetNumExplicitHs(0)
        at.SetNoImplicit(False)
        # And clear the radical count. Parsing '[NH:1]' leaves a
        # valence-incomplete atom, so sanitizing the FRAGMENT assigned it an
        # unpaired electron, and that integer rides through CombineMols
        # untouched by the H flags above. Found by this module's own
        # check_no_radicals on its first run.
        at.SetNumRadicalElectrons(0)

    prod = combo.GetMol()
    try:
        Chem.SanitizeMol(prod)
    except Exception as e:
        raise StructureError(f"assemble.attach({label!r}): sanitize failed: {e}")
    check_no_radicals(prod, label=label or "attach")

    # Strip the maps off the inputs before comparing, so a map number is never
    # mistaken for a chemical difference.
    clean = []
    for frag in (A, B):
        f = Chem.Mol(frag)
        for a in f.GetAtoms():
            a.SetAtomMapNum(0)
        clean.append(f)
    check_stereo_preserved(prod, clean, label=label or "attach")
    return Chem.MolToSmiles(prod)


def esterify(acid_smiles, alcohol_smiles, label=""):
    """Acid + alcohol -> ester. THE ACID LOSES ITS OH; THE ALCOHOL KEEPS ITS O.

    That asymmetry is the whole point. In a sterol ester the stereocentre is
    C3 of the sterol, adjacent to the sterol's oxygen -- so deleting the
    sterol's OH is what inverted 537 rows (F-013), and deleting the acid's OH
    instead touches no centre. It is also what the chemistry does: the
    carbonyl oxygen that survives into the ester comes from the alcohol.

    Returns canonical SMILES, after asserting stereo preservation.
    """
    acid = _mol(acid_smiles, "acid")
    alc = _mol(alcohol_smiles, "alcohol")

    patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    hits = acid.GetSubstructMatches(patt)
    if len(hits) != 1:
        raise StructureError(
            f"assemble.esterify({label!r}): the acid has {len(hits)} free "
            "carboxylic acid group(s), need exactly 1. With two, which one "
            "reacts is a decision for DECISIONS.md, not for this function.")
    c_idx, _o_dbl, o_h = hits[0]

    # Refuse if the OH we are about to delete is adjacent to a stereocentre.
    centres = set(stereocentres(acid))
    nbrs = {n.GetIdx() for n in acid.GetAtomWithIdx(o_h).GetNeighbors()}
    if nbrs & centres:
        raise StructureError(
            f"assemble.esterify({label!r}): the acid's leaving OH is adjacent "
            f"to stereocentre(s) {sorted(nbrs & centres)}. Deleting it is the "
            "F-013 mechanism. Reverse the roles, or protect the centre.")

    alc_patt = Chem.MolFromSmarts("[OX2H1]")
    ahits = [h[0] for h in alc.GetSubstructMatches(alc_patt)]
    if len(ahits) != 1:
        raise StructureError(
            f"assemble.esterify({label!r}): the alcohol has {len(ahits)} free "
            "OH group(s), need exactly 1. Which one esterifies is a decision "
            "for DECISIONS.md.")
    o_alc = ahits[0]

    rw = Chem.RWMol(Chem.CombineMols(acid, alc))
    o_alc_off = o_alc + acid.GetNumAtoms()
    rw.AddBond(c_idx, o_alc_off, Chem.BondType.SINGLE)
    # Delete the ACID's hydroxyl oxygen -- last, and on the acid side only.
    rw.RemoveAtom(o_h)
    oa = rw.GetAtomWithIdx(o_alc_off - (1 if o_alc_off > o_h else 0))
    oa.SetNumExplicitHs(0)
    oa.SetNoImplicit(False)

    prod = rw.GetMol()
    try:
        Chem.SanitizeMol(prod)
    except Exception as e:
        raise StructureError(
            f"assemble.esterify({label!r}): sanitize failed: {e}")

    check_no_radicals(prod, label=label or "esterify")
    check_stereo_preserved(prod, [acid, alc], label=label or "esterify")
    return Chem.MolToSmiles(prod)


def check_no_radicals(mol, label=""):
    """Raise if any atom carries an unpaired electron.

    The Siegwart radical carbon parsed, canonicalised, and had self-consistent
    descriptors -- it was found by a reviewer comparing against the drawn
    structure, not by any check. An unpaired electron in a lipid is almost
    always an arithmetic slip in a hydrogen count during assembly, so it is
    worth one line on every build.
    """
    m = _mol(mol, "product")
    bad = [(a.GetIdx(), a.GetSymbol(), a.GetNumRadicalElectrons())
           for a in m.GetAtoms() if a.GetNumRadicalElectrons()]
    if bad:
        raise StructureError(
            f"assemble.check_no_radicals({label!r}): {len(bad)} radical "
            f"atom(s): {bad[:8]}.\n"
            "An unpaired electron here means a hydrogen count was decremented "
            "by hand during assembly. Set NumExplicitHs(0) and "
            "NoImplicit(False) on the attachment atoms and let valence decide. "
            "A radical parses, canonicalises and gives self-consistent "
            "descriptors, so nothing downstream will catch it.")
    return True


def parity_preserved(product, fragment, label=""):
    """True if every ASSIGNED stereocentre of `fragment` keeps its SPATIAL
    arrangement in `product`, judged by neighbour order, not by CIP code.

    WHY THIS EXISTS. `check_stereo_preserved` compares CIP code COUNTS, and a
    CIP code is a property of the WHOLE molecule: a transformation that changes
    a remote substituent can RELABEL a centre without moving a single atom.
    Measured 2026-11-10: (S)-glycidol's C2 is `S`; its decanoate ester is `R`.
    Nothing inverted -- the free OH leaves the C1/C3 comparison tied and the
    ester breaks the tie. The count check cannot tell that from an inversion,
    so it REFUSED a correct build (reported from `chamolstad_41691524`).

    Counts genuinely cannot separate the two cases, because a real inversion
    also preserves the total. Parity can: map each fragment atom onto the
    product by substructure match, then compare the centre's chiral tag after
    permuting its neighbours into the SAME order on both sides. The permutation
    parity is what encodes handedness, and it does not care what the CIP rules
    make of the substituents.

    Returns (ok, detail). `ok` is False only if a centre genuinely flipped.
    """
    prod = _mol(product, "product")
    frag = _mol(fragment, "fragment")
    for m in (prod, frag):
        Chem.AssignStereochemistry(m, cleanIt=True, force=True)

    match = prod.GetSubstructMatch(frag, useChirality=False)
    if not match:
        return True, {"mapped": False,
                      "why": "fragment is not a substructure of the product; "
                             "parity cannot be compared, falling back"}
    f2p = {fi: pi for fi, pi in enumerate(match)}

    flipped, checked = [], 0
    for fa in frag.GetAtoms():
        tag = fa.GetChiralTag()
        if tag not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                       Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            continue
        pa = prod.GetAtomWithIdx(f2p[fa.GetIdx()])
        if pa.GetChiralTag() not in (Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                                     Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            flipped.append((fa.GetIdx(), "lost its tag in the product"))
            continue
        checked += 1
        # Neighbours of the centre, in each molecule's own bond order, mapped
        # into the fragment's numbering so the two lists are comparable.
        fn = [n.GetIdx() for n in fa.GetNeighbors()]
        pn = [n.GetIdx() for n in pa.GetNeighbors()]
        p2f = {v: k for k, v in f2p.items()}
        pn_in_f = [p2f.get(i) for i in pn]
        # Atoms the product gained at this centre (e.g. an acyl O) have no
        # fragment image; a substitution that replaces a hydrogen is fine.
        common = [i for i in pn_in_f if i is not None and i in fn]
        if len(common) < 3:
            continue                      # too little shared to judge
        perm = [fn.index(i) for i in common]
        inv = sum(1 for a in range(len(perm)) for b in range(a + 1, len(perm))
                  if perm[a] > perm[b])
        same_tag = (fa.GetChiralTag() == pa.GetChiralTag())
        # even permutation + same tag, or odd permutation + opposite tag,
        # both mean the handedness survived.
        if same_tag != (inv % 2 == 0):
            flipped.append((fa.GetIdx(), "parity inverted"))
    return (not flipped), {"mapped": True, "centres_compared": checked,
                           "flipped": flipped}


def canonical(smiles):
    """Canonical SMILES, or raise. Never returns the input unchanged on failure."""
    m = _mol(smiles)
    return Chem.MolToSmiles(m)
