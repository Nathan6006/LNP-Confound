#!/usr/bin/env python3
"""The peptoid monomer pool, as drawn in main-article Figure 1B.

WHAT A PEPTOID MONOMER IS HERE. Submonomer synthesis (main-article Experimental
Methods, and Figure 1B): the resin-bound amine is bromoacetylated, then a
PRIMARY amine displaces the bromide. So every residue is an N-substituted
glycine --

        -N(R)-CH2-C(=O)-

-- and the "monomer" named by a three-letter code is the substituent R, which
is whatever the displacing amine carried. The chain is built C-terminus to
N-terminus on Rink amide resin, so after cleavage the C-terminus is a PRIMARY
CARBOXAMIDE (-C(=O)NH2) and the N-terminus is the last amine added, still
bearing its hydrogen.

THE ONE EXCEPTION, and it is load-bearing. Dea (diethanolamine) and Eae
(2-(ethylamino)ethanol) are SECONDARY amines, not primary. Displacing with a
secondary amine gives a backbone nitrogen that already carries two arms and NO
hydrogen, and which therefore cannot be acylated again -- so such a residue can
only ever be the N-terminal one. Getting this wrong is a 15.011 Da (one NH)
error, and it is exactly the error the calculated-exact-mass check caught on
this curation's first run: peptoid 227 came out at 1077.95 against the paper's
printed 1062.94 until the N-terminal hydrogen was removed for Dea. Figure 1B
marks Dea with an asterisk; this is what the asterisk is about.

WHERE THE SIDE-CHAIN STRUCTURES COME FROM. Figure 1B draws them. It names none
of them, so the drawing is the only definition. The figure raster
(nn4c05513_0001.jpg, 4183x2463 px, 600 dpi) is the single main-article figure
PMC served; it was read visually.

WHERE THE CHAIN LENGTHS COME FROM -- three independent routes, not one:
  1. the three-letter code itself (Hex/Oct/Dec/Dod). NOT TRUSTED on its own:
     CLAUDE.md section 5, "name tokens are not chain lengths", KS_2024's C24 is
     two C4 tails.
  2. SI Table S4's printed "Total Lipid Carbons" column, over 34 peptoids. This
     is an arithmetic constraint on the whole monomer set at once and it is
     exactly determined -- see stage03's check carbon_arithmetic.
  3. the PDF's OWN stroked bond coordinates for Dod and Ehx, measured with
     curation/lib chains on SI page 20's Figure S7 legend. 13 and 9 merged
     vertices, attachment point included, so 12 and 8 carbons.
All three agree. Route 2 is the one that actually pins the set; route 3 is the
only one that touches geometry.

WHAT IS DELIBERATELY ABSENT. The 39-monomer pool of Figure 1B is NOT fully
transcribed here, because this curation enumerates only peptoids 222-266 (SI
Table S9) and those use 18 monomers. Monomers that appear only in the 221-member
screening library (Table S1) are listed with resolved=False and no SMILES, for
two reasons: that library has no extractable value (00_value_sources.json), and
Figure 1B gives the SAME three-letter code 'Mpp' to two different structures, so
a blanket transcription would be guessing on at least one of them.
"""

# --- side chains R, written as SMILES fragments attached at their FIRST atom
#     to the backbone nitrogen. Read off Figure 1B.
#
# `n_carbons` is the count used by the Table S4 arithmetic check. It is recorded
# here so the check can FAIL; it is not where the count was decided.
#
# `amine_order`: "primary" (normal residue, N-terminal copy keeps one H) or
# "secondary" (N-terminal only, no H, two arms on the backbone nitrogen).

MONOMERS = {
    # ---- lipid monomers (Figure 1B, green panel) -------------------------
    "Hex": dict(sidechain="CCCCCC",                  n_carbons=6,  amine_order="primary",
                kind="lipid", drawn_as="C6H13", branched=False, unsaturated=False,
                note="n-hexyl"),
    "Oct": dict(sidechain="CCCCCCCC",                n_carbons=8,  amine_order="primary",
                kind="lipid", drawn_as="C8H17", branched=False, unsaturated=False,
                note="n-octyl"),
    "Dec": dict(sidechain="CCCCCCCCCC",              n_carbons=10, amine_order="primary",
                kind="lipid", drawn_as="C10H21", branched=False, unsaturated=False,
                note="n-decyl; n-decylamine is in the SI Materials list"),
    "Dod": dict(sidechain="CCCCCCCCCCCC",            n_carbons=12, amine_order="primary",
                kind="lipid", drawn_as="C12H23", branched=False, unsaturated=False,
                note="n-dodecyl. Figure 1B prints the formula C12H23, which is "
                     "one hydrogen short of dodecyl (C12H25); SI Figure S7's "
                     "legend prints 'Dod (C12)' and draws 12 vertices. See "
                     "08_contradictions.csv."),
    "Ole": dict(sidechain="CCCCCCCC/C=C\\CCCCCCCC",  n_carbons=18, amine_order="primary",
                kind="lipid", drawn_as="Oleyl", branched=False, unsaturated=True,
                note="(Z)-octadec-9-en-1-yl. Figure 1B writes the word 'Oleyl' "
                     "and draws no bond-line structure, so the 9Z position and "
                     "geometry are the standard meaning of the word, not a "
                     "reading. Flagged structure_assumed on every row using it."),
    "Ehx": dict(sidechain="CC(CC)CCCC",              n_carbons=8,  amine_order="primary",
                kind="lipid", drawn_as="2-ethylhexyl", branched=True, unsaturated=False,
                note="2-ethylhexyl, drawn in full in Figure 1B and again in SI "
                     "Figure S7's legend. RACEMIC at C2: 2-ethylhexylamine is "
                     "sold and used as the racemate (SI Materials: VWR) and the "
                     "paper draws no wedge. Stored without a stereocentre."),

    # ---- headgroup monomers used in 222-233 (Figure 1B, cyan + salmon) ---
    "Aet": dict(sidechain="CCN",                     n_carbons=2,  amine_order="primary",
                kind="cationic", drawn_as="H2N-CH2CH2-", branched=False, unsaturated=False,
                note="2-aminoethyl"),
    "Apr": dict(sidechain="CCCN",                    n_carbons=3,  amine_order="primary",
                kind="cationic", drawn_as="H2N-CH2CH2CH2-", branched=False, unsaturated=False,
                note="3-aminopropyl"),
    "Nhe": dict(sidechain="CCO",                     n_carbons=2,  amine_order="primary",
                kind="hydrophilic", drawn_as="HO-CH2CH2-", branched=False, unsaturated=False,
                note="2-hydroxyethyl"),
    "Nhp": dict(sidechain="CCCO",                    n_carbons=3,  amine_order="primary",
                kind="hydrophilic", drawn_as="HO-CH2CH2CH2-", branched=False, unsaturated=False,
                note="3-hydroxypropyl"),
    "Apd": dict(sidechain="C(CO)CO",                 n_carbons=3,  amine_order="primary",
                kind="hydrophilic", drawn_as="HO-CH2-CH(-)-CH2-OH", branched=True, unsaturated=False,
                note="1,3-dihydroxypropan-2-yl, from 2-amino-1,3-propanediol "
                     "(SI Materials). Attaches through the CENTRAL carbon -- "
                     "Figure 1B draws the squiggle there, and SI Figure S7 "
                     "draws it again in the full structures of 236 and 266."),
    "Dea": dict(sidechain=None, arms=("CCO", "CCO"), n_carbons=4,  amine_order="secondary",
                kind="hydrophilic", drawn_as="(HO-CH2CH2)2N-", branched=False, unsaturated=False,
                note="diethanolamine. SECONDARY amine: the backbone nitrogen "
                     "itself bears both hydroxyethyl arms and no hydrogen. "
                     "Figure 1B marks it with an asterisk and draws the squiggle "
                     "ON THE NITROGEN."),
    "Eae": dict(sidechain=None, arms=("CC", "CCO"),  n_carbons=4,  amine_order="secondary",
                kind="hydrophilic", drawn_as="HO-CH2CH2-N(C2H5)-", branched=False, unsaturated=False,
                note="2-(ethylamino)ethanol (SI Materials, TCI). SECONDARY "
                     "amine, same treatment as Dea. NOT DRAWN in Figure 1B's "
                     "pool -- the code appears only in SI Table S9. The "
                     "structure is the Materials-list reagent, and the "
                     "calculated exact mass confirms the composition."),
    "Nme": dict(sidechain="CCOC",                    n_carbons=3,  amine_order="primary",
                kind="hydrophilic", drawn_as="CH3-O-CH2CH2-", branched=False, unsaturated=False,
                note="2-methoxyethyl. SMILES reads attachment-C, then O, then "
                     "methyl: C-O-C."),
    "Nmp": dict(sidechain="CCCOC",                   n_carbons=4,  amine_order="primary",
                kind="hydrophilic", drawn_as="CH3-O-CH2CH2CH2-", branched=False, unsaturated=False,
                note="3-methoxypropyl (3-methoxypropan-1-amine, SI Materials)."),
    "Mip": dict(sidechain="CCCn1ccnc1C",             n_carbons=7,  amine_order="primary",
                kind="cationic", drawn_as="2-methylimidazol-1-yl-propyl", branched=False, unsaturated=False,
                note="3-(2-methyl-1H-imidazol-1-yl)propyl."),
    "Mop": dict(sidechain="CCCN1CCOCC1",             n_carbons=7,  amine_order="primary",
                kind="cationic", drawn_as="morpholin-4-yl-propyl", branched=False, unsaturated=False,
                note="3-morpholinopropyl (3-morpholinopropylamine, SI Materials)."),

    # ---- used in Table S9 but NOT drawn in Figure 1B's pool --------------
    "Dme": dict(sidechain=None,                      n_carbons=4,  amine_order="primary",
                kind="cationic", drawn_as=None, branched=False, unsaturated=False,
                note="UNRESOLVED CONNECTIVITY. The code appears in SI Table S9 "
                     "(peptoid 231) and nowhere in Figure 1B's drawn pool. The "
                     "printed calculated exact mass pins the COMPOSITION at "
                     "C4H10N, but C4H10N is satisfied by 2-(dimethylamino)ethyl, "
                     "3-(methylamino)propyl, 2-(methylamino)propyl and others, "
                     "and mass cannot see connectivity. Figure 1B's 'Mae' is "
                     "drawn as methylaminopropyl and has the SAME formula, which "
                     "is precisely why the mass does not decide. NOT ENUMERATED; "
                     "peptoid 231 is excluded from 03_enumerated.csv."),
}

# Monomers that appear in SI Table S1's 221-member screening library but are not
# used by peptoids 222-266. Listed so their absence is visible, not silent.
NOT_ENUMERATED_POOL = [
    "Mae", "Eed", "Mpl", "Mpp", "Mpi", "Imp", "Pip", "Bpp", "Spm", "Tae",
    "Met", "Pet", "Smb", "Mpe", "Fur", "Frm", "Cym", "Iso", "Mme", "Cme",
    "Cpr", "Ode", "Frm",
]

# Backbone, stated once.
C_TERMINUS = "C(N)=O"   # primary carboxamide, from Rink amide resin
BACKBONE_UNIT_FORMULA = "C2H2NO"   # -N-CH2-C(=O)- contributed by every residue
