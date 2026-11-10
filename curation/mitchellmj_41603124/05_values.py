#!/usr/bin/env python3
"""Stage 05 - every VALUE this paper reports, with an admissibility verdict.

One row per (panel x quantity x group) block. `corpus_admissible` is yes/no
with a reason on EVERY row, and the counts sum to the total.

The corpus endpoint vocabulary is delivery-only. Measured from data/ on
2026-11-10: Experiment_method over 1,803 in-vivo rows takes exactly four
values -- luminescence_normalized (836), protein_adundance_normalized (766),
editing_efficiency_normalized (141), luminescence_relative_to_Spikevax (60).
Particle biodistribution, cellular uptake, viability, size, PDI, encapsulation
efficiency, cytokines, body weight and tumour volume are none of these.
"""
import csv, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "05_values.csv")

POOL = ("one pool of 133 b-mRNA LNPs co-injected per mouse; barcode pool "
        "narrowed from 200 to 134 (Supp Table 2)")

R = []


def add(**kw):
    R.append(kw)


# ---------------------------------------------------------------- the screen
for strain, fig in (("C57BL/6", "Figure 2a + Supplementary Figure 14a"),
                    ("APOE-/-", "Figure 2b + Supplementary Figure 14b")):
    add(panel=fig,
        quantity="Mean normalized intracellular b-mRNA accumulation",
        printed_units="log(1 + Mean normalized accumulation), continuous "
                      "colour ramp, scale ticks at 1, 2, 3",
        readout_class="cargo_accumulation_by_RNA_sequencing",
        cargo="barcoded EGFP mRNA (b-mRNA)",
        biological_system="%s mouse" % strain,
        organ_cell="16 sorted populations: blood B/T/NK/monocyte; lung "
                   "endothelial/epithelial/leukocyte; liver endothelial/"
                   "hepatocyte/Kupffer/leukocyte; spleen B/T/NK/monocyte/DC",
        route="intravenous", dose="not stated per LNP (pooled)",
        timepoint="short, 'relatively short time following treatment'; "
                  "no hours given",
        n_values=120 * 16,
        attributable="YES - Supp Fig 14 labels x = 'Amine or thiol' 1..20 and "
                     "y = 'Dialkyl maleate' 4/6.2/8/8i/9.2/18, so every cell "
                     "maps to an IL_name. Figure 2a/2b alone would NOT be: "
                     "its x axis is 'Barcode number' with NO tick labels.",
        pooled="yes", pool_identity=POOL, pool_size=133,
        compositional="YES",
        corpus_admissible="no",
        reason="ENDPOINT. The quantity is how much barcoded mRNA is PRESENT in "
               "the sorted cells, measured by sequencing RNA, and the authors "
               "sampled deliberately before peak translation ('we opted to "
               "analyze b-mRNA accumulation at a relatively short time ... "
               "peak mRNA translation occurs relatively late'). It is cargo "
               "accumulation, not expression, and maps to none of the four "
               "Experiment_method values. The corpus's only barcoded "
               "experiment, LR_2023 (407 of 1,803 in-vivo rows), carries "
               "Cargo_type = peptide_barcode and Experiment_method = "
               "protein_adundance_normalized -- a TRANSLATED readout. Pooling "
               "is not the objection; the readout is.")

add(panel="Figure 2c, 2d, 2e",
    quantity="Pairwise Pearson r between sorted populations",
    printed_units="r_p, 0.00-1.00 colour ramp",
    readout_class="derived_statistic", cargo="b-mRNA",
    biological_system="C57BL/6 and APOE-/-",
    organ_cell="population x population", route="intravenous",
    dose="n/a", timepoint="n/a", n_values=0,
    attributable="n/a - a statistic over the screen, not a per-LNP value",
    pooled="yes", pool_identity=POOL, pool_size=133, compositional="n/a",
    corpus_admissible="no",
    reason="Not a per-LNP delivery measurement; a correlation matrix derived "
           "from the screen that is itself inadmissible.")

# ---------------------------------------------------- the validation subset
LNPS = ["LNP 9 (2D8)", "LNP 37 (7D4)", "LNP 94 (17D8i)", "LNP 96 (17D9.2)",
        "LNP 97 (18D4)", "LNP 112 (20D8i)", "SM-102"]
for organ in ("liver", "spleen"):
    add(panel="Figure 3b",
        quantity="Total %s luminescence" % organ,
        printed_units="p/s, log axis 10^5..10^12 (x) / 10^4..10^9 (y)",
        readout_class="luminescence", cargo="NanoLuc mRNA",
        biological_system="C57BL/6 mouse", organ_cell=organ,
        route="intravenous", dose="0.1 mg/kg encapsulated mRNA",
        timepoint="6 h", n_values=7,
        attributable="YES - 7 groups named in the panel's own legend: "
                     + "; ".join(LNPS),
        pooled="no - single-plex validation injections", pool_identity="",
        pool_size="", compositional="no",
        corpus_admissible="no",
        reason="ENDPOINT IS ADMISSIBLE; STRUCTURE AND PRECISION ARE NOT. Of "
               "the 7 groups only 18D4 has a structure this paper determines "
               "(MS-confirmed, 45.8x margin) and only SM-102 comes from the "
               "registry; 2D8 and 7D4 are undeterminable (amines 2 and 7 and "
               "maleate 8 are drawn nowhere), and 17D9.2/20D8i/17D8i are "
               "drawn but not readable to a structure. data/ has 0 blank "
               "IL_SMILES in 30,998 rows, so 5 of 7 cannot be emitted -- and "
               "the two undeterminable ones are the two LOWEST liver values, "
               "so dropping them truncates the bottom of the range that every "
               "z-score in the experiment is computed over. See DECISIONS D-2.")

add(panel="Figure 3d, 3h",
    quantity="DiR+ cells", printed_units="% of population",
    readout_class="particle_accumulation", cargo="DiR-labelled EGFP mRNA LNP",
    biological_system="C57BL/6 mouse",
    organ_cell="liver endothelial/hepatocyte/Kupffer/leukocyte (d); spleen "
               "B/T/NK/macrophage/DC (h)",
    route="intravenous", dose="0.6 mg/kg encapsulated mRNA", timepoint="12 h",
    n_values=7 * 9,
    attributable="YES - same 7 groups, by legend colour",
    pooled="no", pool_identity="", pool_size="", compositional="no",
    corpus_admissible="no",
    reason="ENDPOINT. DiR is a lipophilic membrane dye: this is where the "
           "PARTICLE went, i.e. cellular uptake / biodistribution, which the "
           "corpus vocabulary excludes.")

for organ, pops, n in (("liver", "endothelial/hepatocyte/Kupffer/leukocyte", 4),
                       ("spleen", "B/T/NK/macrophage/DC", 5)):
    add(panel="Figure 3f" if organ == "liver" else "Figure 3j",
        quantity="EGFP+ cells in %s" % organ, printed_units="% of population",
        readout_class="protein_abundance", cargo="EGFP mRNA",
        biological_system="C57BL/6 mouse", organ_cell=pops,
        route="intravenous", dose="0.6 mg/kg encapsulated mRNA",
        timepoint="12 h", n_values=7 * n,
        attributable="YES - same 7 groups, by legend colour",
        pooled="no", pool_identity="", pool_size="", compositional="no",
        corpus_admissible="no",
        reason="ENDPOINT IS ADMISSIBLE (positivity rate for a translated "
               "reporter = protein abundance). Blocked by the SAME structure "
               "gap as Figure 3b: 5 of the 7 ionizable lipids have no "
               "determinable structure. See DECISIONS D-2.")

for p, q, why in (
    ("Figure 4b, 4c", "Tumour volume; survival",
     "Not a delivery measurement - therapeutic efficacy."),
    ("Figure 4d-4g", "Reciprocal endpoint anti-OVA IgG / IgG1 / IgG2c titre",
     "Not a delivery measurement - humoral immune response."),
    ("Figure 5e, 5g, 5i", "ApoE and beta2-GPI adsorption to single LNPs",
     "Not a delivery measurement - protein corona composition."),
    ("Supplementary Figure 1", "Hydrodynamic diameter of C12-494 LNPs",
     "Not a delivery measurement - particle size."),
    ("Supplementary Figure 10", "LNP diameter; b-mRNA concentration; "
     "entrapment efficiency",
     "Not a delivery measurement - size and encapsulation efficiency."),
    ("Supplementary Figure 15, 16", "Body weight; tumour growth curves",
     "Not a delivery measurement - tolerability and efficacy."),
    ("Supplementary Figure 17", "Protein adsorption to SORT LNPs",
     "Not a delivery measurement - protein corona."),
    ("Figure 1f, 1g", "b-mRNA pool-narrowing enrichment/depletion",
     "Not an LNP measurement at all - it characterises the BARCODES, in a "
     "single formulation (C12-494), to pick which barcodes to use."),
):
    add(panel=p, quantity=q, printed_units="as printed",
        readout_class="not_delivery", cargo="", biological_system="",
        organ_cell="", route="", dose="", timepoint="", n_values="",
        attributable="", pooled="", pool_identity="", pool_size="",
        compositional="", corpus_admissible="no", reason=why)


def main():
    cols = ["panel", "quantity", "printed_units", "readout_class", "cargo",
            "biological_system", "organ_cell", "route", "dose", "timepoint",
            "n_values", "attributable", "pooled", "pool_identity", "pool_size",
            "compositional", "corpus_admissible", "reason"]
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in R:
            w.writerow({c: r.get(c, "") for c in cols})
    yes = sum(1 for r in R if r["corpus_admissible"] == "yes")
    no = sum(1 for r in R if r["corpus_admissible"] == "no")
    assert yes + no == len(R), "verdict counts do not sum"
    print(json.dumps({"value_blocks": len(R), "corpus_admissible_yes": yes,
                      "corpus_admissible_no": no,
                      "counts_sum": yes + no == len(R)}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
