#!/usr/bin/env python3
"""The CBB / CBBZ library, transcribed from the SI's own IUPAC names.

SOURCE OF CONNECTIVITY: `source/41467_2026_69077_MOESM1_ESM.pdf`, pp 18-50,
"Synthesis of CBB intermediates and lipidoids". Every compound there carries a
full IUPAC name, 1H/13C NMR, and "HRMS calculated for <formula>, <m/z>". The
NAME is the connectivity; the HRMS formula is an independent atom-count check
that is BLIND TO CONNECTIVITY (it cannot see which nitrogen an arm sits on).

This file holds only the paper's own library logic -- scaffold templates and
tail fragments. Every structural OPERATION is done by curation/lib/assemble.py.

WHERE THE RING METHYL SITS. Three phrasings in the SI describe the same ring:

  CBB21  "...methylene)bis(3-METHYL-4,1-phenylene))bis(ethylazanediyl)..."
  CBBZ1  "...methylene)bis(2-METHYL-4,1-phenylene))bis(ethylazanediyl)..."
  CBB6   "...(4-(ethyl(...)amino)-2-METHYLphenyl)methyl)-3-METHYLphenyl)(ethyl)amino..."

The two divalent "4,1-phenylene" forms disagree with each other under any one
locant convention. CBB6's two MONOVALENT prefixes do not: in a monovalent
substituent prefix position 1 IS the attachment point, so "(4-(amino)-2-
methylphenyl)" attached to the central methine puts the methyl ORTHO to that
methine, and "-3-methylphenyl)(ethyl)amino" -- a phenyl cited as a substituent
ON the amine nitrogen -- puts the methyl META to the nitrogen. Same ring, and
they agree. The synthesis corroborates it: both series descend from I2,
"3-((ethyl(m-tolyl)amino)methyl)benzenesulfonic acid", whose N-alkyl-m-toluidine
is attacked PARA TO THE NITROGEN (C4), which is ORTHO to the C3 methyl.

So: methyl ortho to the central methine, meta to the nitrogen, in BOTH series.
The CBBZ "2-methyl" locant is a nomenclature inconsistency, not a second
structure; recorded in 08_contradictions.csv.

INDEPENDENT CORROBORATION OF THE WHOLE DYE CORE, not used to build it:
C47H49N3O7S2 (the free acid of Brilliant Blue G) is 832.00 g/mol, so the sodium
salt is 854.01 -- and the SI weighs out "commercial brilliant G250 (854 mg, 1.0
mmol)". That is a round number the authors USED and never printed as a formula.
"""

# ---------------------------------------------------------------- scaffolds
# Ring-closure digits are kept globally distinct: 1 quinoid, 2/5 toluidine,
# 3/6 benzenesulfonyl, 4 the aniline head. `{a}` and `{b}` are the two
# sulfur substituents, which is what makes an unsymmetric arm pair expressible
# (CBB5, CBB16, CBB18 all need it).
_DYE = ("{head}N=C1C=CC(=C("
        "c2ccc(N(CC)Cc3cccc(c3)S(=O)(=O){a})cc2C)"
        "c5ccc(N(CC)Cc6cccc(c6)S(=O)(=O){b})cc5C"
        ")C=C1")

HEAD_ETHOXY   = "CCOc4ccc(cc4)"            # 4-ethoxyphenylimino  (CBB series)
HEAD_PHENOL   = "Oc4ccc(cc4)"              # 4-hydroxyphenylimino (CBBZ click handle)
HEAD_PROPYNYL = "C#CCOc4ccc(cc4)"          # 4-(prop-2-yn-1-yloxy)phenylimino (CBBZ0)

# sulfur substituents
SO3H = "O"                                               # -SO2-OH, free sulfonic acid
K1_A = "NCCCN(C)CCC[NH2:2]"      # -SO2-NH-(CH2)3-N(CH3)-(CH2)3-NH2, arm a
K1_B = "NCCCN(C)CCC[NH2:3]"      #                                    arm b
K2_A = "NCCC[N+](C)(C)CCC[NH2:2]"   # the N-methylated (quaternary) linker
K2_B = "NCCC[N+](C)(C)CCC[NH2:3]"
# triethylene-glycol sulfonic ESTER arm, terminal OH palmitoylated (CBB15/16)
TEG_PALMITATE = "OCCOCCOCCOC(=O)CCCCCCCCCCCCCCC"


def dye(head, a, b):
    return _DYE.format(head=head, a=a, b=b)


# ------------------------------------------------------------------- tails
# Acyl fragments in ALDEHYDE form: the mapped carbon carries one explicit H,
# which assemble.attach hands back to RDKit on bonding, so R-CH=O + H2N-  ->
# R-C(=O)-NH-.  No atom is removed anywhere, so the F-013 mechanism cannot
# arise. Named exactly as the SI names the acid that was used.
ACYL = {
 "butyryl":        ("[CH:1](=O)CCC", "sodium butyrate / butyric acid, C4:0"),
 "lauroyl":        ("[CH:1](=O)CCCCCCCCCCC", "lauric acid, C12:0"),
 "palmitoyl":      ("[CH:1](=O)CCCCCCCCCCCCCCC", "palmitoyl chloride, C16:0"),
 "triacontanoyl":  ("[CH:1](=O)" + "C" * 29, "melissic acid, C30:0"),
 "15-hexadecynoyl":("[CH:1](=O)" + "C" * 13 + "C#C",
                    "alkynyl palmitic acid = hexadec-15-ynoic acid"),
 "linoleoyl":      ("[CH:1](=O)CCCCCCC/C=C\\C/C=C\\CCCCC",
                    "linoleic acid, (9Z,12Z)-octadeca-9,12-dienoic"),
 "dha":            ("[CH:1](=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC",
                    "docosahexaenoic acid, (4Z,7Z,10Z,13Z,16Z,19Z)-docosa-"
                    "4,7,10,13,16,19-hexaenoic"),
 "retinoyl":       ("[CH:1](=O)/C=C(C)/C=C/C=C(C)/C=C/C7=C(C)CCCC7(C)C",
                    "all-trans-retinoic acid, named verbatim in the SI as the "
                    "reagent; the SI's own locants for the amide "
                    "('2,6-dimethyl-9-(2,3,3-trimethylcyclohex-1-en-1-yl)') are "
                    "a mis-numbering of retinoic acid's 3,7-dimethyl-9-(2,6,6-"
                    "trimethylcyclohex-1-en-1-yl), and the HRMS formula "
                    "C101H136N9O7S2+ fixes the composition as retinoyl"),
 "trifluoroacetyl":("[CH:1](=O)C(F)(F)F", "trifluoroacetic acid, from the boc "
                    "deprotection step"),
 "formyl_for_carbamate": ("[CH:1](=O)O",
                    "formic acid. Bonded to the core's terminal N this gives a "
                    "CARBAMIC ACID, which assemble.esterify then esterifies "
                    "with cholesterol -- so cholesterol KEEPS its own oxygen "
                    "and its C3 is never touched (the F-013 rule, and what "
                    "cholesteryl chloroformate + amine actually does)."),
}

# The pentapeptide. SEQUENCE IS CONTESTED -- see DECISIONS.md D-1. Residues are
# written C-terminus-first so the free COOH is the one that acylates the core.
CPP5 = {
 "KLVPM": ("[CH:1](=O)[C@@H](CCSC)NC(=O)[C@@H]1CCCN1C(=O)[C@@H](C(C)C)"
           "NC(=O)[C@@H](CC(C)C)NC(=O)[C@@H](N)CCCCN",
           "as written in the SI synthesis heading: 'boc-protected CPP5 "
           "peptide (KLVPM)'"),
 "KLPVM": ("[CH:1](=O)[C@@H](CCSC)NC(=O)[C@@H](C(C)C)NC(=O)[C@@H]1CCCN1"
           "C(=O)[C@@H](CC(C)C)NC(=O)[C@@H](N)CCCCN",
           "as written in the main text: 'a cell-penetrating peptide (CPP5, "
           "sequence KLPVM)'"),
}

# ---------------------------------------- CBBZ azide tails, post-click
# Written as the 1,4-disubstituted triazole that the CuAAC forms: N1 carries the
# azide's carbon, C4 carries the methylene that was the propargyl CH2 and is
# mapped for attachment to the dye's phenolic oxygen. Ring written from C4:
# C4-C5-N1-N2-N3, so `c7cn(R)nn7` is unambiguously 1,4.
_DODECYL      = "CCCCCCCCCCCC"
_HEXADECYL    = "CCCCCCCCCCCCCCCC"
_HEPTADECAN9  = "C(CCCCCCCC)CCCCCCCC"
_UNDECYL      = "CCCCCCCCCCC"


def _triazole(r):
    return "[CH2:2]c7cn(" + r + ")nn7"


def _bis_propionate(spacer, alkoxy):
    """N(spacer)(CH2CH2CO2R)2 -- the '3,3'-(R-azanediyl)dipropionate' motif."""
    return (spacer + "N(CCC(=O)O" + alkoxy + ")CCC(=O)O" + alkoxy)


AZIDE_TAIL = {
 "Z1":  (_bis_propionate("CCC", _DODECYL),
         "I6, didodecyl 3,3'-((3-azidopropyl)azanediyl)dipropionate"),
 "Z2":  (_bis_propionate("CCOCCOCC", _DODECYL),
         "I8, didodecyl 3,3'-((2-(2-(2-azidoethoxy)ethoxy)ethyl)azanediyl)dipropionate"),
 "Z3":  (_bis_propionate("CCC", _HEXADECYL),
         "I10, dihexadecyl 3,3'-((3-azidopropyl)azanediyl)dipropionate"),
 "Z4":  (_bis_propionate("CCC", _HEPTADECAN9),
         "I12, di(heptadecan-9-yl) 3,3'-((3-azidopropyl)azanediyl)dipropionate"),
 "Z5":  ("CCN(CCCCCCCC(=O)O" + _HEPTADECAN9 + ")CCCCCC(=O)O" + _UNDECYL,
         "I15, heptadecan-9-yl 8-((2-azidoethyl)(6-oxo-6-(undecyloxy)hexyl)"
         "amino)octanoate"),
 "Z6":  (_bis_propionate("CCOCCOCCOCCOCC", _DODECYL),
         "I17, dodecyl 1-azido-15-(3-(dodecyloxy)-3-oxopropyl)-3,6,9,12-"
         "tetraoxa-15-azaoctadecan-18-oate"),
 "Z7":  (_bis_propionate("CCOCCOCCOCCOCC", _HEXADECYL),
         "hexadecyl 1-azido-15-(3-(hexadecyloxy)-3-oxopropyl)-3,6,9,12-"
         "tetraoxa-15-azaoctadecan-18-oate"),
 "Z8":  (_bis_propionate("CCOCCOCC", _HEXADECYL),
         "I20, dihexadecyl 3,3'-((2-(2-(2-azidoethoxy)ethoxy)ethyl)azanediyl)"
         "dipropionate"),
 "Z9":  (_bis_propionate("CCOCCOCCOCCOCCOCCOCC", _HEPTADECAN9),
         "heptadecan-9-yl 1-azido-21-(3-(heptadecan-9-yloxy)-3-oxopropyl)-"
         "3,6,9,12,15,18-hexaoxa-21-azatetracosan-24-oate"),
 "Z10": (_bis_propionate("CCOCCOCCOCCOCCOCCOCC", _DODECYL),
         "I24, dodecyl 1-azido-21-(3-(dodecyloxy)-3-oxopropyl)-3,6,9,12,15,18-"
         "hexaoxa-21-azatetracosan-24-oate"),
}

# -------------------------------------------------- one whole-molecule template
# CBB9 is a MACROCYCLE: the two terminal amines are bridged by one oxalyl, so
# closing it needs two bonds and assemble.attach makes one. It has no
# combinatorial tail, so the whole molecule is transcribed the same way the
# scaffolds above are, then canonicalised and radical-checked by assemble.
# SI name: "...-7,28-dithia-4,8,12,16,19,23,27,31-octaaza-1,3(1,4),6,29(1,3)-
# tetrabenzenacyclohentriacontaphane-17,18-dione 7,7,28,28-tetraoxide (9)";
# "17,18-dione" on adjacent skeletal atoms is the oxalamide bridge.
# The macrocyclic bond is ring-closure digit 8, opened ON THE CENTRAL METHINE
# and closed on the second toluidine ring. Toluidine ring 5 is written from its
# NITROGEN, so its methyl must sit at the position ADJACENT TO the closure (the
# methyl is ortho to the central methine, meta to the nitrogen) -- `c5cc(C)c8cc5`,
# not `c5ccc8cc5C`, which would put it ortho to the nitrogen instead. The first
# version of this string left ring 8 out entirely; it parsed, canonicalised and
# gave self-consistent descriptors, and only the HRMS atom count caught it
# (C63H83 against the SI's C63H81, i.e. one bond short).
CBB9_MACROCYCLE = (
    "CCOc4ccc(cc4)N=C1C=CC(=C8"
    "c2ccc(N(CC)Cc3cccc(c3)S(=O)(=O)NCCCN(C)CCCNC(=O)C(=O)NCCCN(C)CCCNS(=O)(=O)"
    "c6cccc(c6)CN(CC)c5cc(C)c8cc5"
    ")cc2C)C=C1"
)

# HRMS as the SI prints it: compound -> (formula string, charge, calculated m/z)
# Transcribed by hand from the PDF text layer, on purpose: a formula that is
# re-derived from our own SMILES would compare a value to itself.
HRMS = {
 "CBB1":  ("C117H174N9O9S2", 3, 637.7620),
 "CBB2":  ("C119H178N9O9S2", 3, 647.1057),
 "CBB3":  ("C105H144N9O7S2", 1, 1707.0625),
 "CBB4":  ("C107H150N9O7S2", 3, 579.0361),
 "CBB5":  ("C76H97N6O7S2",   1, 1269.6855),
 "CBB6":  ("C93H144N9O7S2",  1, 1563.0625),
 "CBB7":  ("C95H150N9O7S2",  3, 531.3706),
 "CBB8":  ("C93H136N9O7S2",  1, 1554.9999),
 "CBB9":  ("C63H82N9O7S2",   1, 1140.5773),
 "CBB10": ("C121H200N9O7S2", 1, 1955.5007),
 "CBB11": ("C97H144N9O7S2",  1, 1611.0625),
 "CBB12": ("C85H130N9O7S2",  1, 1450.9373),
 "CBB13": ("C69H96N9O7S2",   1, 1226.6869),
 "CBB14": ("C71H102N9O7S2",  3, 418.9109),
 "CBB15": ("C91H134N3O15S2", 1, 1572.9251),
 "CBB16": ("C69H92N3O11S2",  1, 1202.6168),
 "CBB17": ("C65H82F6N9O7S3", 1, 1278.5678),
 "CBB18": ("C90H131F3N15O11S3", 1, 1750.9262),
 "CBB19": ("C115H182N21O15S4", 3, 741.7664),
 "CBB20": ("C101H136N9O7S2", 1, 1650.9999),
 "CBB21": ("C61H84N9O5S2",   1, 1086.6031),
 "CBBZ1": ("C81H112N7O11S2", 1, 1422.7856),
 "CBBZ2": ("C84H118N7O13S2", 1, 1496.8224),
 "CBBZ3": ("C89H128N7O11S2", 1, 1534.9108),
 "CBBZ4": ("C91H132N7O11S2", 1, 1562.9421),
 "CBBZ5": ("C92H134N7O11S2", 1, 1576.9578),
 "CBBZ6": ("C88H126N7O15S2", 1, 1584.8748),
 "CBBZ7": ("C96H142N7O15S2", 1, 1697.0000),
 "CBBZ8": ("C92H134N7O13S2", 1, 1608.9476),
 "CBBZ9": ("C102H154N7O17S2", 1, 1813.0837),
 "CBBZ10":("C92H134N7O17S2", 1, 1672.9273),
 # CBB22 is commercial Brilliant Blue G and is not re-characterised in the SI;
 # CBBZ0 is reported with 1H/13C NMR only, no HRMS. Both absent on purpose.
}


# ---------------------------------------------------------------------------
# THE PAPER'S OWN BLOCK TABLE, read off Supplementary Figure 1 (SI p2, rendered
# at 400 dpi and LOOKED AT -- renders/si_p2.png, renders/crop_table.png,
# renders/crop_linkers.png, renders/crop_core.png).
#
# This is a SECOND, INDEPENDENT statement of every assignment: the IUPAC names
# on pp 18-50 say it in words, this figure says it as a Markush table and a
# drawing. stage07_qc.py P10 compares the two compound by compound. The drawing
# also settles the regiochemistry argument in P8 by eye: in both toluidine
# rings the methyl sits on the carbon ADJACENT to the central methine and three
# bonds from the N(ethyl)(benzyl), which is the ortho/meta assignment we built.
#
# The table covers compounds 1-21 only. CBB22 is absent from it because CBB22 is
# the purchased dye, not a synthesised member.
DRAWN_LIPIDS = {
 "L1": "Cholesteryl Formate", "L2": "Docosahexaenoate", "L3": "Palmitate",
 "L4": "Alkynyl palmitate", "L5": "Melissate", "L6": "Linoleate",
 "L7": "Laurate", "L8": "Butyrate", "L9": "Trifluoroacetate",
 "L10": "Retinoate", "L11": "CPP5", "L12": "(-SO2-)OH",
 "H": "none (free primary amine)", "-COCO-": "oxalyl bridge, intramolecular",
}
DRAWN_LINKERS = {
 "K1": "H2N-(CH2)3-N(CH3)-(CH2)3-NH2",
 "K2": "H2N-(CH2)3-N(+)(CH3)2-(CH2)3-NH2",
 "K3": "HO-CH2CH2-O-CH2CH2-O-CH2CH2-OH (triethylene glycol)",
}
DRAWN_TABLE = {
 "CBB1":  ("L1", "L1", "K1"),   "CBB2":  ("L1", "L1", "K2"),
 "CBB3":  ("L2", "L2", "K1"),   "CBB4":  ("L2", "L2", "K2"),
 "CBB5":  ("L2", "L12", "K1"),  "CBB6":  ("L3", "L3", "K1"),
 "CBB7":  ("L3", "L3", "K2"),   "CBB8":  ("L4", "L4", "K1"),
 "CBB9":  ("-COCO-", None, "K1"), "CBB10": ("L5", "L5", "K1"),
 "CBB11": ("L6", "L6", "K1"),   "CBB12": ("L7", "L7", "K1"),
 "CBB13": ("L8", "L8", "K1"),   "CBB14": ("L8", "L8", "K2"),
 "CBB15": ("L3", "L3", "K3"),   "CBB16": ("L3", "L12", "K3"),
 "CBB17": ("L9", "L9", "K1"),   "CBB18": ("L9", "L11", "K1"),
 "CBB19": ("L11", "L11", "K1"), "CBB20": ("L10", "L10", "K1"),
 "CBB21": ("H", "H", "K1"),
}
# how the figure's lipid codes map onto the tail ids stage03 derives from the
# IUPAC names. Written out so P10 compares two readings, not one reading twice.
LIPID_TO_TAIL = {
 "L1": "cholesteryl carbamate", "L2": "dha", "L3": "palmitoyl",
 "L4": "15-hexadecynoyl", "L5": "triacontanoyl", "L6": "linoleoyl",
 "L7": "lauroyl", "L8": "butyryl", "L9": "trifluoroacetyl",
 "L10": "retinoyl", "L11": "CPP5", "L12": "free sulfonic acid",
 "H": "none", "-COCO-": "oxalyl bridge (intramolecular)",
}
LINKER_TO_SCAFFOLD = {"K1": "bis-sulfonamide, K1 linker",
                      "K2": "bis-sulfonamide, K2 quaternised linker",
                      "K3": "bis(TEG sulfonic ester)"}


# ---------------------------------------------------------------------------
# THE CBBZ BLOCK TABLE, read off Supplementary Figure 3 (SI p4, 400 dpi,
# renders/crop_fig3.png). "CBB Structure = (R1, R2, n1, n2, X)".
#
# The drawn generic arm is  N-CH2-(CH2)n -C(=O)-O-R , so n=1 is a propionate
# and n=6 an octanoate; the drawn azide side is  N3-CH2CH2-[X]-N< .
#
# Z10 DISAGREES WITH THE TEXT and the text wins -- see 08_contradictions.csv.
DRAWN_TABLE_Z = {
 "Z1":  ("L13", "L13", 1, 1, "K4"),  "Z2":  ("L13", "L13", 1, 1, "K5"),
 "Z3":  ("L14", "L14", 1, 1, "K4"),  "Z4":  ("L15", "L15", 1, 1, "K4"),
 "Z5":  ("L15", "L16", 6, 4, "N.A."), "Z6": ("L13", "L13", 1, 1, "K6"),
 "Z7":  ("L14", "L14", 1, 1, "K6"),  "Z8":  ("L14", "L14", 1, 1, "K5"),
 "Z9":  ("L15", "L15", 1, 1, "K7"),  "Z10": ("L14", "L14", 1, 1, "K7"),
}
DRAWN_R = {"L13": _DODECYL, "L14": _HEXADECYL, "L15": _HEPTADECAN9,
           "L16": _UNDECYL}
DRAWN_R_NAME = {"L13": "Dodecyl", "L14": "Hexadecyl",
                "L15": "Heptadecan-9-yl", "L16": "Undecyl"}
# X as drawn: K4 is a CH2, K5/K6/K7 are (PEG)2/4/6, N.A. is nothing.
DRAWN_X = {"K4": "CCC", "K5": "CC" + "OCC" * 2, "K6": "CC" + "OCC" * 4,
           "K7": "CC" + "OCC" * 6, "N.A.": "CC"}
DRAWN_X_NAME = {"K4": "CH2", "K5": "(PEG)2", "K6": "(PEG)4", "K7": "(PEG)6",
                "N.A.": "none"}
# Z10's R block is the one cell where SI Fig 3 and the SI text disagree. The
# text is taken, so the exception is DECLARED here rather than left to a
# tolerance: the product's own IUPAC name says "21-(3-(dodecyloxy)-3-oxopropyl)",
# the reagent is named "dodecyl 1-azido-21-(3-(dodecyloxy)-3-oxopropyl)-...-
# oate", and the HRMS formula C92H134N7O17S2+ AND its printed m/z 1672.9273
# both reproduce on dodecyl to 1e-4. Hexadecyl would be C100H150N7O17S2+.
DRAWN_Z_EXPECTED_MISMATCH = {
 "Z10": "SI Fig 3 draws R = L14 (hexadecyl); the SI text and the HRMS both "
        "say dodecyl (L13). We emit dodecyl.",
}
