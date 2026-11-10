#!/usr/bin/env python3
"""Stage 01 -- the building blocks this paper's library is made from.

PMID 23651337, Wang/Sun/Alberti/Xu, ACS Synth Biol 2012;1(9).
"A combinatorial library of unsaturated lipidoids for efficient intracellular
gene delivery."

WHAT IS IN HAND. Only the 3-page Supporting Information. The main article is
behind the ACS subscription wall and has no PMC deposit; the checkable
negatives are in 00_value_sources.json -> acquisition.unobtainable.

CONSEQUENCE FOR THIS STAGE. The SI DRAWS NO STRUCTURE AND NAMES NO AMINE.
The only chemical identity obtainable is the TAIL reagent, named in the
published abstract. The 16 amine heads -- and therefore all 16 lipidoids --
are NOT enumerable. This stage records that boundary; it does not cross it.

Run:  python3 -I stage01_blocks.py
"""
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

# The published abstract, verbatim from the Europe PMC REST core record
# (EXT_ID:23651337 AND SRC:MED, resultType=core, fetched 2026-11-10). This is
# the publisher's own text and it is the ONLY primary statement of the
# chemistry available to this curation.
ABSTRACT = (
    "A combinatorial library of unsaturated lipidoids was synthesized through "
    "the Michael addition of amines to oleyl acrylamide. Their capability in "
    "facilitating in vitro gene delivery was evaluated by transfecting HeLa "
    "cells with EGFP-encoding plasmid DNA and mRNA. The preliminary screening "
    "results indicated that lipidoids with unsaturated oleyl tails are "
    "superior transfection agents compared to saturated lipidoids with "
    "n-octadecyl tails under the same conditions. The different transfection "
    "abilities of the unsaturated and saturated lipidioids were ascribed to "
    "the large, tightly packed lipoplexes of saturated lipidoids. The "
    "potential applications of the library of lipidoids were further expanded "
    "by looking at their ability to transfect fibroblasts as well as "
    "different cancerous cell lines."
)

# ---------------------------------------------------------------------------
# The tail reagents.
#
# "oleyl acrylamide" is a chemical NAME, and a name is not a structure until
# something external pins it (CLAUDE.md 3). What pins it here:
#   * "acrylamide", not "acrylate" -- the acryloyl group is on NITROGEN. The
#     Xu lab writes both in its other papers and distinguishes them by name
#     (curation/xuq_32662132/01_blocks.csv: O16S-S is ...OCC..., N16S-S is
#     ...NCC..., built from 2-mercaptoethanol vs cysteamine respectively).
#   * "oleyl" = (9Z)-octadec-9-en-1-yl. A single cis C=C at C9 of an
#     18-carbon chain. This is the IUPAC-recognised trivial name; no part of
#     it is inferred from this paper.
#   * "unsaturated" in the title and "unsaturated oleyl tails" in the abstract
#     confirm exactly one degree of unsaturation is intended in the tail.
#
# A SMILES FRAGMENT HAS A DIRECTION THAT MASS CANNOT SEE. The direction here is
# fixed by the name, not by the formula: the acryloyl carbonyl is bonded to N,
# N is bonded to C1 of the chain, and the C=C is at C9-C10 counting from the
# N-bearing end. Writing the chain from the methyl end would put the alkene at
# C9 from the WRONG end (= C10 from N) with an IDENTICAL molecular formula and
# an IDENTICAL exact mass. 07_qc.json seeds precisely that defect.
# ---------------------------------------------------------------------------
BLOCKS = [
    {
        "block_id": "oleyl_acrylamide",
        "block_type": "acrylamide_tail",
        "chemical_name": "N-[(9Z)-octadec-9-en-1-yl]prop-2-enamide "
                         "(= oleyl acrylamide = N-oleylacrylamide)",
        "SMILES": "DETERMINABLE -- built in 03_enumerated.csv "
                  "(tail_SMILES) by stage03 through curation/lib/assemble. "
                  "Deliberately not restated here: two copies of a structure "
                  "is how they drift apart.",
        "n_michael_acceptors": 1,
        "evidence": "NAMED in the published abstract: 'the Michael addition of "
                    "amines to oleyl acrylamide'. Also 'unsaturated oleyl "
                    "tails' in the same abstract and 'unsaturated lipidoids' "
                    "in the title.",
        "confidence": "high",
        "what_is_drawn": "NOTHING. The SI contains three figures (two TEM "
                         "micrographs, one line plot, one bar chart) and no "
                         "chemical structure of any kind. The main-text "
                         "figure that would draw the library is behind the "
                         "ACS subscription wall.",
        "note": "The one block in this paper whose structure is determinable. "
                "Derived from its chemical name, not from a drawing.",
    },
    {
        "block_id": "octadecyl_acrylamide",
        "block_type": "acrylamide_tail",
        "chemical_name": "N-octadecylprop-2-enamide "
                         "(= n-octadecyl acrylamide)",
        "SMILES": "DETERMINABLE -- see 03_enumerated.csv, as above.",
        "n_michael_acceptors": 1,
        "evidence": "PARTIALLY NAMED. The abstract names the comparator only "
                    "as 'saturated lipidoids with n-octadecyl tails'. It does "
                    "NOT print the phrase 'octadecyl acrylamide'. The "
                    "acrylamide linkage is carried over by parallel "
                    "construction from the unsaturated series, which the "
                    "abstract says was run 'under the same conditions'.",
        "confidence": "medium",
        "what_is_drawn": "NOTHING -- as above.",
        "note": "ASSUMED FEATURE: that the saturated comparator is the "
                "acrylAMIDE and not an acrylATE. n-octadecyl is the exact "
                "saturated analogue of oleyl (both C18, differing only in the "
                "C9-C10 double bond), which is what makes the paper's "
                "unsaturated-vs-saturated comparison a one-variable one; but "
                "no printed sentence states the linkage for this member.",
    },
]

# ---------------------------------------------------------------------------
# The amine heads. THE GAP.
# ---------------------------------------------------------------------------
AMINE_GAP = {
    "block_type": "amine_head",
    "how_many_exist": 16,
    "how_many_enumerable": 0,
    "evidence_for_16": "Figure S3's x-axis carries 16 major tick marks, each "
                       "at the centre of one bar, labelled 1..16 ('Lipidoid'). "
                       "Measured, not counted by eye: stage05_values.py "
                       "x_tick_vs_bar_centre.",
    "why_zero_enumerable": "The SI names no amine, draws no amine, prints no "
                           "formula, no HRMS and no NMR. There is no amine "
                           "catalogue page, no synthesis section and no "
                           "experimental section in the 3 pages we have.",
    "name_token_warning": "The labels 1..16 are a LOCAL SERIES INDEX. They are "
                          "NOT the Xu/Anderson catalogue amine numbers, which "
                          "run 75-400 (see curation/xuq_32520527/01_blocks.csv "
                          "and curation/xuq_32438474/01_blocks.csv: 75, 76, 77, "
                          "78, 80, 81, 82, 87, 90, 93, 304, 306, 400). The "
                          "sibling curation/xuq_24574196 states the same for "
                          "its own paper: \"The paper's local label '1' is NOT "
                          "catalogue amine 1.\" Reading these 1..16 as "
                          "catalogue numbers would assign 16 measured values "
                          "to 16 wrong molecules.",
    "consequence": "A lipidoid in this library is (amine head + n x oleyl "
                   "acrylamide). With the head unknown, the lipidoid is "
                   "unknown. 03_enumerated.csv therefore contains ZERO "
                   "lipidoids. Nothing is guessed.",
}

# ---------------------------------------------------------------------------
# Lipidoid 16's two forms: 16I and 16II. Named in the SI, not resolvable.
# ---------------------------------------------------------------------------
SERIES_AMBIGUITY = {
    "tokens": ["16I", "16II"],
    "where": "Figure S2 caption (p. S2): 'Ethidium bromide (EB) exclusion "
             "assay of 16I/DNA and 16II/DNA'. Figure S1 caption (p. S2): 'TEM "
             "images of 16II/DNA (A) and 16II/mRNA complexes (B)'.",
    "what_is_certain": "Lipidoid 16 exists in two forms distinguished by a "
                       "roman numeral, and Figure S2's own legend names the "
                       "two plotted series 'Unsaturated Lipidoid' and "
                       "'Saturated Lipidoid' -- so the numeral distinguishes "
                       "the oleyl series from the n-octadecyl series.",
    "what_is_NOT_determinable": "WHICH numeral is which. The caption lists "
                                "'16I ... and 16II', the legend lists "
                                "'Unsaturated ... Saturated'; nothing states "
                                "that the two orders correspond, and the SI "
                                "defines neither numeral. See DECISIONS.md D-2.",
    "a_tempting_inference_that_is_NOT_made": "Figure S1 shows 16II with both "
        "DNA and mRNA, which looks like the lead compound, and the abstract "
        "says the unsaturated lipidoids are superior -- so 16II 'should' be "
        "unsaturated. But the same abstract ascribes the difference to 'the "
        "large, tightly packed lipoplexes of SATURATED lipidoids', which is a "
        "reason to show a TEM of the saturated one. The two readings point "
        "opposite ways and the SI settles neither. Not assigned.",
}


def main():
    rows = []
    for b in BLOCKS:
        rows.append(b)
    out = os.path.join(HERE, "01_blocks.csv")
    cols = ["block_id", "block_type", "chemical_name", "SMILES",
            "n_michael_acceptors", "evidence", "confidence", "what_is_drawn",
            "note"]
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})

    qc = {
        "abstract_verbatim": ABSTRACT,
        "abstract_source": "Europe PMC REST, search?query=EXT_ID:23651337 AND "
                           "SRC:MED&resultType=core, fetched 2026-11-10. "
                           "abstractText length = 797 characters.",
        "blocks_emitted": len(rows),
        "tails_enumerable": 2,
        "amine_heads": AMINE_GAP,
        "series_numeral_ambiguity": SERIES_AMBIGUITY,
        "library_size_arithmetic": {
            "distinct_structures_determinable": 0,
            "distinct_structures_that_exist": "NOT DETERMINABLE from the SI. "
                "16 unsaturated lipidoids are individually measured in Figure "
                "S3, and at least one saturated counterpart (lipidoid 16's "
                "other form) exists. Whether the saturated series has 16 "
                "members or only the one shown in Figure S2 is not stated "
                "anywhere in the SI.",
            "formulations": "NOT RECOVERABLE. The SI states no lipidoid:DNA "
                "weight ratio, no N:P, no helper-lipid composition, no dose "
                "and no particle concentration. It has no experimental "
                "section.",
            "individually_measured_in_the_file_we_have": 16,
            "individually_measured_on_a_DELIVERY_endpoint": 0,
        },
    }
    with open(os.path.join(HERE, "01_blocks_qc.json"), "w") as fh:
        json.dump(qc, fh, indent=2)

    print(f"01_blocks.csv: {len(rows)} blocks (both tails; ZERO amine heads)")
    print(f"amine heads existing / enumerable: "
          f"{AMINE_GAP['how_many_exist']} / {AMINE_GAP['how_many_enumerable']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
