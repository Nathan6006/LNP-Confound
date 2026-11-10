#!/usr/bin/env python3
"""Add the 7 quality columns and fix two IDs, producing the merge-ready table.

WHY THIS IS A SEPARATE STAGE. `07_NG_2025_rows.csv` was built against the
72-column schema, before the quality/provenance columns existed. Rather than
edit another session's stages 01-10, this reads their output and writes
`08_NG_2025_merge.csv` in the current 79-column schema, so their pipeline stays
intact and every flag added here is auditable in one place.

Column names come from `provenance/schema_quality.json`, never restated here --
that single-source rule is why these columns exist at all (F-001).

THE FLAGS, and the evidence for each. Settled from the Tier-2 review
(`provenance/review/round_2026-11-07/report_NG_2025.md`), with each claim
re-checked against the primary source before being acted on.

  confidence = medium on every crude lipid (321 rows)
      The paper, verbatim (Methods, main text line 128): "the solvents were
      evaporated to obtain **310 different crude lipid compounds**. These crude
      compounds were subsequently utilized to prepare LNPs for in vitro
      screening." Only C-a16 was made at scale, column-purified and
      characterised. So every emitted structure except the commercial benchmark
      is the NOMINAL major product of a one-pot reaction that was never purified
      and never characterised. The curation recorded this nowhere, and it is the
      single largest honest caveat on the slice. I verified the sentence in the
      PDF text layer rather than taking the reviewer's quotation.

  confidence = low on the 28 b3/b4 lipids
      Amines b3 and b4 are drawn TWO mutually exclusive ways in the same paper:
      as an aminal (H2N-CH2-NR2) in the four amine-list panels and the four
      dialkyl-ketone product panels, and as a hydrazine (N-N direct) in the
      pyrrole and acetophenone product panels. The two differ by CH2 = 14 Da and
      the paper gives no mass, formula or NMR for any b3/b4 lipid, so nothing in
      it can decide. The aminal is emitted; the alternative is kept in
      `07_flags.csv`.

  positional_ambiguity = 1 on every crude lipid (321 rows)
      A Mannich on a phenol gives ortho AND para products, both ortho positions
      of a 4-alkylphenol are open, and nothing was purified. The paper draws one
      isomer and we store that isomer; a consumer should know it is a convention
      and not a measurement.

  structure_assumed = 1 on the 28 b3/b4 lipids and the 22 tail-D rows
      b3/b4: the assumed feature is WHICH of two drawings of the same reagent.
      tail D: emitted as linear 4-n-nonylphenol as drawn, but commercial
      "4-nonylphenol" is predominantly a BRANCHED technical isomer mixture --
      this column's own documented example. Tail C deliberately gets NO flag:
      its C-a16 NMR shows delta 0.90 (t, J = 6.8 Hz, 12H), a linear terminal
      CH3, which rules out the tert-octyl isomer. Worth stating because it is
      the obvious next question.

  value_is_bound = blank on all 333
      Deliberate, and a real strength of this curation: every in-vitro value is
      the publisher's own Source Data, not a digitised colormap, and the low
      tail has no repeated clamp value (minimum 0.845098 = log10 7, all distinct).
      Nothing here is censored or floored.

WHAT IS NOT FLAGGED, and why -- this is the review's most important finding.
The folder's README asserts that "15 of 52" 1:2:1 dialkyl-ketone lipids are
drawn as the open-chain mono-adduct, and names them. The reviewer rendered NINE
of the fifteen named lipids from the native Supplementary Fig 8/9 rasters and
found **every one drawn as the closed 4-oxopiperidine ring**, including 6 of 6
on tail I. So the count does not reproduce. All 104 are enumerated as rings
either way -- no structure is affected -- but a `confidence = low` resting on
that count would be a flag carrying a false reason, which is worse than no flag
(the F-011 lesson). The 104 therefore carry `medium` from the crude-product rule
and nothing else, and the README's count is retracted pending the mechanical
hole-count the reviewer describes.

    python3 -I curation/ng_2025/stage11_quality.py
"""
import csv
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", ".."))
SRC = os.path.join(HERE, "07_NG_2025_rows.csv")
OUT = os.path.join(HERE, "08_NG_2025_merge.csv")
SCHEMA = os.path.join(ROOT, "provenance", "schema_quality.json")

CRUDE = ("Structure is the nominal product of an unpurified crude Mannich "
         "reaction: the paper evaporated solvent to obtain 310 crude lipid "
         "compounds and screened those directly, characterising only C-a16. "
         "Value is the authors' own Source Data, not a digitisation.")
B34 = ("Amine b3/b4 is drawn as an aminal (H2N-CH2-NR2) in the amine-list and "
       "dialkyl-ketone panels and as a hydrazine (N-N) in the pyrrole and "
       "acetophenone panels; the two differ by CH2 and the paper gives no mass, "
       "formula or NMR for any b3/b4 lipid. The aminal is emitted; the "
       "alternative is in 07_flags.csv.")
TAIL_D = ("Phenolic tail D is emitted as linear 4-n-nonylphenol as drawn; "
          "commercial 4-nonylphenol is predominantly a branched technical "
          "isomer mixture.")
OB10 = ("The Fig. 1 Source Data make this the highest value in the pyrrole "
        "panel (4.07e4); Supplementary Fig 17b re-measures it ~24x lower while "
        "five companion lipids on the same panel agree within 1.9x.")
VIVO = ("Dose 2.5 ug is derived from the printed 0.125 mg/kg and the ~20 g body "
        "weight stated in Methods.")
VIVO_NORM = (
    "z-scored SEPARATELY for liver (intravenous, Fig. 1m) and muscle "
    "(intramuscular injection site, Fig. 1p): different route, different "
    "tissue, different imaging session. The panels' absolute scales carry "
    "printed multipliers and so the offset is recoverable in principle, but it "
    "is not recorded in this curation, and pooling would assert that a liver "
    "IV signal and a muscle IM injection-site signal are directly comparable. "
    "Liver-vs-muscle Spearman over the 11 shared lipids is 0.982. NOTE that "
    "Experiment_value_within_model is therefore NOT what its own definition "
    "says for these rows -- it re-standardises the union of two already "
    "flattened blocks (see F-028).")
MC3 = ("Commercial benchmark, formulated with DSPC rather than DOPE exactly as "
       "Methods state. The matching Fig. 1p benchmark is deliberately held back: "
       "that panel's axis labels it DLin-MC3-DMA while the Source Data column "
       "for the same three numbers is headed SM-102, and those are different "
       "molecules.")


def main():
    quality = [c["name"] for c in json.load(open(SCHEMA))["columns"]]
    rows = list(csv.DictReader(open(SRC)))
    hdr = list(rows[0])
    for c in quality:
        if c not in hdr:
            hdr.append(c)

    n_b34 = n_d = n_crude = 0
    for r in rows:
        name = r["IL_name"].strip()
        vivo = r["Model"] == "in_vivo"
        is_mc3 = "MC3" in name

        # E4: Formulation_ID must start with its own Experiment_ID. 23 in-vivo
        # rows carried the NG_2025 prefix under Experiment_ID NG_2025_2; only
        # 292 of 30,665 corpus rows break that invariant and all of them are one
        # known upstream oddity.
        eid, fid = r["Experiment_ID"], r["Formulation_ID"]
        if not fid.startswith(eid + "_"):
            r["Formulation_ID"] = fid.replace("NG_2025_", eid + "_", 1)

        r["value_normalization_unit"] = (
            "experiment_x_model_type" if vivo else "experiment")

        notes = []
        if is_mc3:
            r["confidence"] = "high"
            notes.append(MC3)
        else:
            n_crude += 1
            r["positional_ambiguity"] = "1"
            notes.append(CRUDE)
            if re.search(r"-b(3|4)$", name):
                r["confidence"] = "low"
                r["structure_assumed"] = "1"
                notes.append(B34)
                n_b34 += 1
            else:
                r["confidence"] = "medium"
            if name.startswith("D-"):
                r["structure_assumed"] = "1"
                notes.append(TAIL_D)
                n_d += 1
            if name == "O-b10":
                notes.append(OB10)
            if vivo:
                notes.append(VIVO)
        if vivo:
            notes.append(VIVO_NORM)
        r["curation_note"] = " ".join(notes)

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in hdr})

    import collections
    conf = collections.Counter((r["Model"], r["confidence"]) for r in rows)
    print(f"{os.path.basename(OUT)}: {len(rows)} rows x {len(hdr)} cols")
    print(f"  crude (medium/low + positional_ambiguity): {n_crude}")
    print(f"  b3/b4 low + structure_assumed: {n_b34}   tail-D structure_assumed: {n_d}")
    for k, v in sorted(conf.items()):
        print(f"  {k[0]:9} confidence={k[1] or '(blank)':7} {v}")
    bad = [r["Formulation_ID"] for r in rows
           if not r["Formulation_ID"].startswith(r["Experiment_ID"] + "_")]
    print(f"  Formulation_ID prefix violations after fix: {len(bad)}")
    assert not bad, bad
    assert all(not r.get("value_is_bound") for r in rows), "nothing here is censored"


if __name__ == "__main__":
    main()
