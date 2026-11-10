"""Reconstruct the paper's OWN definition of 'normalized delivery' from its
worked example (SI Appendix Fig S1b, p.2) and check it cell by cell.

WHY THIS IS THE DECISIVE STAGE. The scope question for this paper is whether a
pooled barcode readout is comparable in kind to the corpus's per-LNP values.
That cannot be settled by reading the prose; it is settled by the arithmetic.
Fig S1b prints three tables -- raw counts, normalised counts, and counts
normalised to input -- and prints `Total (%) = 100` under BOTH normalised
tables. If that reproduces, the readout is a COMPOSITIONAL share of a pool:
every LNP's value is a function of every other member of its pool.

The brief's suggested free check ("normalised delivery across all members of
one pool must sum to 1 or to the stated normalisation") is exactly this, and
the paper drew the answer without anybody recomputing it.

SEEDING. The agreement count must MOVE when a defect is seeded. Seeds that
merely raise are no-ops; this file's own first chain-check seeds were.
"""
import sys, json

SEED = sys.argv[1] if len(sys.argv) > 1 else "none"

# --- transcribed from SI Fig S1b, table 1 (raw counts) ---------------------
COLS = ["Lung M1", "Lung M2", "Lung M3", "Kidney M1", "Kidney M2", "Kidney M3"]
RAW = {           # lnp: [6 sample counts] + input
    "LNP1": ([100,  80, 200, 300, 200, 250], 100),
    "LNP2": ([ 50,  45, 110, 100,  60,  70], 120),
    "LNP3": ([120, 105, 250, 150,  90, 110], 110),
}
PRINTED_TOTALS = [270, 230, 560, 550, 350, 430]
PRINTED_INPUT_TOTAL = 330
# --- table 2 (normalised counts %) and table 3 (normalised to input %) -----
PRINTED_T2 = {"LNP1": [37, 35, 36, 55, 57, 58], "LNP2": [19, 20, 20, 18, 17, 16],
              "LNP3": [44, 46, 45, 27, 26, 26]}
PRINTED_T2_INPUT = {"LNP1": 30, "LNP2": 36, "LNP3": 33}
PRINTED_T3 = {"LNP1": [40, 38, 39, 58, 60, 61], "LNP2": [17, 18, 18, 16, 15, 14],
              "LNP3": [44, 45, 44, 26, 25, 24]}

if SEED == "wrong_raw_count":
    RAW["LNP2"] = ([50, 45, 110, 100, 60, 999], 120)
elif SEED == "drop_input_normalisation":
    for k in RAW:
        RAW[k] = (RAW[k][0], 100)            # pretend every input share is equal
elif SEED == "no_renormalisation":
    pass                                      # handled below

lnps = list(RAW)

# step 1: column totals
tot = [sum(RAW[l][0][c] for l in lnps) for c in range(6)]
itot = sum(RAW[l][1] for l in lnps)

# step 2: share of the sample
t2 = {l: [100.0 * RAW[l][0][c] / tot[c] for c in range(6)] for l in lnps}
t2_in = {l: 100.0 * RAW[l][1] / itot for l in lnps}

# step 3: divide by the input share, then RENORMALISE so the pool sums to 100
ratio = {l: [t2[l][c] / t2_in[l] for c in range(6)] for l in lnps}
rsum = [sum(ratio[l][c] for l in lnps) for c in range(6)]
if SEED == "no_renormalisation":
    t3 = {l: [100.0 * ratio[l][c] for c in range(6)] for l in lnps}
else:
    t3 = {l: [100.0 * ratio[l][c] / rsum[c] for c in range(6)] for l in lnps}

def agree(printed, computed, tol=0.5):
    hit = miss = 0; detail = []
    for l in lnps:
        for c in range(6):
            p, q = printed[l][c], computed[l][c]
            if abs(p - q) <= tol: hit += 1
            else:
                miss += 1
                detail.append({"lnp": l, "col": COLS[c], "printed": p,
                               "computed": round(q, 2)})
    return hit, miss, detail

h1 = sum(1 for c in range(6) if tot[c] == PRINTED_TOTALS[c])
h2, m2, d2 = agree(PRINTED_T2, t2)
h3, m3, d3 = agree(PRINTED_T3, t3)
hin = sum(1 for l in lnps if abs(t2_in[l] - PRINTED_T2_INPUT[l]) <= 0.5)

out = {
 "seed": SEED,
 "step1_column_totals": {"AGREEMENT": "%d/6" % h1, "computed": tot,
                         "printed": PRINTED_TOTALS,
                         "input_total_computed": itot,
                         "input_total_printed": PRINTED_INPUT_TOTAL,
                         "input_total_agrees": itot == PRINTED_INPUT_TOTAL},
 "step2_share_of_sample": {"AGREEMENT": "%d/18" % h2, "misses": d2,
                           "input_column_AGREEMENT": "%d/3" % hin},
 "step3_normalised_to_input_then_RENORMALISED":
     {"AGREEMENT": "%d/18" % h3, "misses": d3},
 "THE_FREE_CHECK_the_authors_drew_but_did_not_compute": {
     "claim_printed_under_table_2": "Total (%) = 100 in every column",
     "computed_table2_column_sums": [round(sum(t2[l][c] for l in lnps), 6) for c in range(6)],
     "claim_printed_under_table_3": "Total (%) = 100 in every column",
     "computed_table3_column_sums": [round(sum(t3[l][c] for l in lnps), 6) for c in range(6)],
     "max_abs_residual_from_100":
        round(max(abs(sum(t[l][c] for l in lnps) - 100.0)
                  for t in (t2, t3) for c in range(6)), 12)},
 "FORMULA_ESTABLISHED":
    "normalized_delivery_i = [ (count_i / SUM_j count_j) / (input_i / SUM_j input_j) ] "
    "/ SUM_k[ same ] x 100%.  It is a COMPOSITIONAL share: it sums to 100% over "
    "the members of the pool, so every LNP's value is a function of which other "
    "LNPs were co-injected with it.",
}
print(json.dumps(out, indent=1))
