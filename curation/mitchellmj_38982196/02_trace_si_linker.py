#!/usr/bin/env python3
"""Stage 02 -- MEASURE the T-MCR linker from the one product the paper draws.

Supplementary Fig. 1a (SI p.3) draws 19T-O12 in full. That drawing is the only
place in the paper where the amidine-thioether-ester linker appears as a
structure rather than as "R1" in a scheme, and it is also what settles how many
arms a POLYAMINE receives -- it shows BOTH primary amines of amine 19 carrying
an arm while the tertiary N-methyl is untouched, under printed stoichiometry
"19 (1 equiv.), Traut's reagent (3 equiv.), O12 (2.5 equiv.)".

Two counts matter and both are measured here, not assumed:
  * CH2 between the amidine carbon and the sulfur  (expected 3, from the
    5-membered ring of Traut's reagent opening)
  * CH2 between the sulfur and the ester carbonyl  (expected 2, from the
    thiol-Michael addition to an acrylate)

The arithmetic cross-check is in 06_checks.py: these two counts enter both MS
anchors, and a linker one CH2 longer or shorter would move each printed mass by
14.

Run: python3 -I 02_trace_si_linker.py
"""
import importlib.util
import json
import os
import sys

import numpy as np
from PIL import Image
from scipy import ndimage

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import chains                                              # noqa: E402
from errors import DetectionError                          # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "trace1", os.path.join(HERE, "01_trace_chains.py"))
T = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(T)

PNG = os.path.join(HERE, "source", "si_fig1_p3_400dpi.png")
SWEEP = (4.0, 5.0, 6.0, 7.0)
OUT = os.path.join(HERE, "02_si_linker_measurements.json")

# Seeds on the drawn 19T-O12, at 400 dpi. Each must be ink or the script
# refuses: a seed that has drifted must not quietly measure its neighbour.
M = [
    dict(id="upper_arm_amidine_to_S", seed=(1576, 535), minus=2,
         ends=("label", "label"),
         role="upper arm: from the amidine carbon's HN label across to the S "
              "label. Vertices minus the two label positions is the amidine "
              "carbon plus the CH2 run"),
    dict(id="upper_arm_S_to_carbonylO", seed=(1729, 531), minus=2,
         ends=("label", "heteroatom"),
         role="upper arm: from the S label through to the carbonyl oxygen. "
              "Three heavy atoms means CH2-CH2-C(=O) and nothing else, i.e. "
              "the thiol added to the acrylate's TERMINAL carbon. The "
              "carbonyl O is drawn as an outline letter and is not in the "
              "glyph map, so its end is checked by adjacent ink"),
    dict(id="upper_arm_tail", seed=(1834, 534), minus=1,
         ends=("label", "methyl"), name_token=12,
         role="upper arm: the ester O-alkyl tail of O12"),
    dict(id="lower_arm_amidine_to_S", seed=(1576, 658), minus=2,
         ends=("label", "label"), role="lower arm, amidine carbon to S"),
    dict(id="lower_arm_S_to_carbonylO", seed=(1728, 661), minus=2,
         ends=("label", "heteroatom"),
         role="lower arm, S through to the carbonyl oxygen"),
    dict(id="lower_arm_tail", seed=(1834, 658), minus=1,
         ends=("label", "methyl"), name_token=12,
         role="lower arm: the ester O-alkyl tail of O12"),
    dict(id="core_upper_N_to_NH", seed=(1500, 585), minus=2,
         ends=("label", "label"),
         role="the drawn product's own amine-19 core: tertiary N to the upper "
              "arm's amidine NH. Measured here, on the SI page, so the C10 "
              "anchor in make_reference.py can be built WITHOUT reusing any "
              "integer from stage 01's reading of the article figure"),
    dict(id="core_lower_N_to_NH", seed=(1501, 612), minus=2,
         ends=("label", "label"),
         role="the drawn product's core: tertiary N to the lower arm's "
              "amidine NH"),
    dict(id="core_N_methyl", seed=(1463, 598), minus=1,
         ends=("methyl", "heteroatom"),
         role="the drawn product's core: the unlabelled bond on the tertiary "
              "N, i.e. its methyl"),
    dict(id="reactant_O12_tail", seed=(1014, 592), minus=1,
         ends=("label", "methyl"), name_token=12,
         role="the O12 acrylate REACTANT's tail, drawn separately on the same "
              "page -- an independent instance of the same chain"),
]


def main():
    g = np.array(Image.open(PNG).convert("L"))
    ink = g < 128
    lab, n = ndimage.label(ink, structure=np.ones((3, 3)))
    objs = ndimage.find_objects(lab)
    boxes = T.glyph_boxes(lab, objs)

    out = {"_what": "Stage 02: the T-MCR linker, measured off the only product "
                    "structure the paper draws in full (Supplementary Fig. 1a, "
                    "19T-O12).",
           "figure": os.path.basename(PNG), "dpi": 400,
           "second_method": "NOT an arc-length ruler. On this page the only "
                            "chains long enough to calibrate one are the three "
                            "tails being measured, which would not be an "
                            "independent check. The checks used instead are "
                            "all external to the trace: the product's two arms "
                            "are drawn separately and must agree atom for "
                            "atom; the acrylate REACTANT's tail is a third "
                            "independent drawing of the same chain; the same "
                            "O12 tail is measured again in stage 01 off a "
                            "different figure at a different scale; and both "
                            "printed masses are reproduced in 07_qc.json.",
           "printed_stoichiometry": "Supplementary Fig. 1a, verbatim: "
                                    "'19 (1 equiv.)  Traut's reagent (3 equiv.)"
                                    "  O12 (2.5 equiv.)'",
           "both_primary_amines_substituted": True,
           "both_primary_amines_basis":
               "The drawn product carries one amidine-thioether-ester arm on "
               "EACH of amine 19's two primary nitrogens, and none on its "
               "tertiary N-methyl. Confirmed numerically by the mass spectrum "
               "printed in the same figure, 401.129 [M+2H]2+ -- see "
               "07_qc.json check P1b.",
           "measurements": {}}

    for m in M:
        cid = T.component_at(m["seed"], ink, lab)
        S, nb = T.skeleton(cid, lab, objs)
        if m["ends"][1] == "junction":
            path, juncs = T.path_to_junction(S, nb)
            checkable = [(path[0], m["ends"][0])]
        else:
            path, ends, juncs = T.longest_path(S, nb)
            if (path[0][0], path[0][1]) > (path[-1][0], path[-1][1]):
                path = path[::-1]
            checkable = [(path[0], m["ends"][0]), (path[-1], m["ends"][1])]

        kinds = []
        for pt, declared in checkable:
            got, d = T.endpoint_kind(pt, boxes, lab, cid)
            # "heteroatom" accepts either test. What the +/-1 rule needs is
            # that the end is an ATOM POSITION and not a bare methyl; whether
            # that atom's letter is its own component or fused into the bond
            # is a drawing accident, and here it differs between the product's
            # two arms (the upper arm's carbonyl O is 29 px from the glyph
            # map, the lower arm's is 12). "methyl" stays discriminating.
            ok = (got == declared or
                  (declared == "heteroatom"
                   and got in ("label", "label_fused")))
            if not ok:
                raise DetectionError(
                    f"02_trace_si_linker[{m['id']}]: endpoint {pt} measures "
                    f"{got!r} (nearest glyph {d:.1f} px), declared "
                    f"{declared!r}. Refusing.")
            kinds.append([got, round(d, 1)])

        tight = {e: chains.raster_chain_length(path, e) for e in (3.0, 4.0, 5.0, 6.0)}
        wide = {e: chains.raster_chain_length(path, e) for e in SWEEP}
        sweep = {"tolerance_3_to_6": {str(k): v for k, v in tight.items()},
                 "tolerance_4_to_7": {str(k): v for k, v in wide.items()}}
        uniq = sorted(set(wide.values()))
        if len(uniq) != 1:
            raise DetectionError(
                f"02_trace_si_linker[{m['id']}]: not stable even over the "
                f"widened sweep: {wide}. Refusing.")
        v = uniq[0]
        resolved = None if len(set(tight.values())) == 1 else (
            "a 3.0 px Douglas-Peucker tolerance adds one spurious vertex to "
            "every LINKER component on this 400 dpi render, and to none of "
            "the three tails, because the linker's steeper bonds leave a "
            "coarser skeleton staircase. The sweep is therefore run at "
            "4.0-7.0 px for this page. That widening is a no-op on every "
            "measurement that had already converged -- asserted below -- and "
            "the value it leaves is confirmed externally by both printed "
            "masses (07_qc.json P1 and P1b).")
        rec = {"role": m["role"], "component": cid, "vertices": v,
               "dp_sweep": sweep, "dp_resolved_by": resolved, "endpoints_measured": kinds,
               "minus": m["minus"], "heavy_atoms": v - m["minus"],
               "count_carbons": chains.count_carbons(
                   raster=v, name_token=m.get("name_token"), label=m["id"])}
        if m.get("name_token") is not None:
            rec["name_token_check"] = chains.check_name_token(
                m["name_token"], v - m["minus"], label=m["id"])
        out["measurements"][m["id"]] = rec

    ms = out["measurements"]
    # The widened sweep must change nothing that had already converged.
    for k, r in ms.items():
        t = set(r["dp_sweep"]["tolerance_3_to_6"].values())
        w = set(r["dp_sweep"]["tolerance_4_to_7"].values())
        if len(t) == 1 and t != w:
            raise DetectionError(
                f"02_trace_si_linker[{k}]: widening the tolerance sweep "
                f"CHANGED a converged measurement ({t} -> {w}). The widening "
                "would then be choosing the answer. Refusing.")
    tails = [k for k in ms if k.endswith("_tail")]
    if len({ms[k]["heavy_atoms"] for k in tails}) != 1:
        raise DetectionError(
            "02_trace_si_linker: the three independently drawn O12 tails on "
            f"this page do not agree: "
            f"{ {k: ms[k]['heavy_atoms'] for k in tails} }.")
    out["tail_agreement"] = {
        "n_independent_drawings_on_this_page": len(tails),
        "heavy_atoms_each": ms[tails[0]]["heavy_atoms"],
        "note": "Two product arms plus the reactant acrylate, drawn "
                "separately, all give the same tail length; stage 01 measures "
                "the same O12 off the article figure at a different scale."}

    for a, b in (("core_upper_N_to_NH", "core_lower_N_to_NH"),
                 ("upper_arm_amidine_to_S", "lower_arm_amidine_to_S"),
                 ("upper_arm_S_to_carbonylO", "lower_arm_S_to_carbonylO"),
                 ("upper_arm_tail", "lower_arm_tail")):
        if ms[a]["heavy_atoms"] != ms[b]["heavy_atoms"]:
            raise DetectionError(
                f"02_trace_si_linker: the molecule's two arms measure "
                f"differently ({a}={ms[a]['heavy_atoms']}, "
                f"{b}={ms[b]['heavy_atoms']}). They are drawn as the same arm.")
    out["core_from_this_drawing"] = {
        "N_to_amidine_NH_carbons": ms["core_upper_N_to_NH"]["heavy_atoms"],
        "n_arms_on_the_tertiary_N": 2,
        "methyl_carbons": ms["core_N_methyl"]["heavy_atoms"],
        "note": "Amine 19 re-measured from the PRODUCT drawing on the SI "
                "page, independently of stage 01's reading of Fig. 1a. "
                "make_reference.py builds the C10 anchor from these numbers "
                "and from the linker below, so the anchor shares neither the "
                "figure, the scale, the calibration nor the assembly code "
                "path with the enumerator."}
    out["linker"] = {
        "amidine_C_plus_CH2_run": ms["upper_arm_amidine_to_S"]["heavy_atoms"],
        "CH2_between_amidine_C_and_S":
            ms["upper_arm_amidine_to_S"]["heavy_atoms"] - 1,
        "S_to_carbonylO_heavy_atoms":
            ms["upper_arm_S_to_carbonylO"]["heavy_atoms"],
        "CH2_between_S_and_carbonyl":
            ms["upper_arm_S_to_carbonylO"]["heavy_atoms"] - 1,
        "note": "Traut's reagent is 2-iminothiolane, a five-membered ring: "
                "ring-opening by a primary amine gives an amidine carbon "
                "followed by three CH2 and a free thiol. The thiol then adds "
                "to the acrylate's terminal carbon, putting two CH2 between S "
                "and the ester carbonyl. Both counts are measured above."}
    with open(OUT, "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"wrote {OUT}")
    print(json.dumps(out["linker"], indent=1))
    for k, v in ms.items():
        print(f"  {k:26s} heavy_atoms={v['heavy_atoms']}")


if __name__ == "__main__":
    main()
