# Stage 8 — Negative Controls and Robustness

The split-half control is a validity GATE: a failure here means the upstream Stage 3/5 machinery is producing signal out of nothing and everything preceding this stage would need to be distrusted. It is reported plainly, not softened.

## Split-half negative control (GATE) — `LNPDB_vitro`

- Largest source by unique lipids: **KZ_2016** (1,470 lipids)

**Separability collapse test**: random 50/50 lipid-level pseudo-source split -> balanced accuracy = **0.496** (chance = 0.500, tolerance +/-0.1). **PASS**

**Concordance-approaches-ceiling test**: within-source (KZ_2016) replicate split-half Spearman = -0.178, but n=14 replicate-group pairs is below Stage 5's own reliability threshold (20) for trusting a concordance estimate. **INCONCLUSIVE (LOW-CONFIDENCE, n too small)** — not scored as a pass or fail, consistent with how Stage 5 treats any concordance built on this few replicate pairs.

## Featurization sensitivity — `LNPDB_vitro`

- Headline (Morgan radius=2, 2048-bit): balanced accuracy = 0.967
- Alternate (Morgan radius=3, 1024-bit): balanced accuracy = **0.971** (chance = 0.034)
- Difference: 0.004 — qualitatively robust to featurization choice.

## Seed sensitivity — `LNPDB_vitro`

- Separability balanced accuracy across seeds [1, 2, 3, 4, 5]: [0.967, 0.968, 0.966, 0.967, 0.966] -> mean=0.967, std=0.001, range=[0.966, 0.968]
- Pooled concordance point estimate (seed-invariant, Spearman is deterministic given data): 0.083
- Bootstrap CI width across seeds [1, 2, 3, 4, 5]: [0.269, 0.251, 0.255, 0.262, 0.261] -> mean=0.260, std=0.006

## Split-half negative control (GATE) — `LNPDB_vivo`

- Largest source by unique lipids: **LR_2023** (383 lipids)

**Separability collapse test**: random 50/50 lipid-level pseudo-source split -> balanced accuracy = **0.490** (chance = 0.500, tolerance +/-0.1). **PASS**

**Concordance-approaches-ceiling test**: within-source (LR_2023) replicate split-half Spearman = 0.717, but n=9 replicate-group pairs is below Stage 5's own reliability threshold (20) for trusting a concordance estimate. **INCONCLUSIVE (LOW-CONFIDENCE, n too small)** — not scored as a pass or fail, consistent with how Stage 5 treats any concordance built on this few replicate pairs.

## Featurization sensitivity — `LNPDB_vivo`

- Headline (Morgan radius=2, 2048-bit): balanced accuracy = 0.982
- Alternate (Morgan radius=3, 1024-bit): balanced accuracy = **0.982** (chance = 0.077)
- Difference: 0.000 — qualitatively robust to featurization choice.

## Seed sensitivity — `LNPDB_vivo`

- Separability balanced accuracy across seeds [1, 2, 3, 4, 5]: [0.988, 0.989, 0.983, 0.981, 0.982] -> mean=0.985, std=0.003, range=[0.981, 0.989]
- Pooled concordance point estimate (seed-invariant, Spearman is deterministic given data): -0.123
- Bootstrap CI width across seeds [1, 2, 3, 4, 5]: [1.514, 1.396, 1.419, 1.45, 1.424] -> mean=1.440, std=0.040

## Gate summary

| dataset | gate | result |
|---|---|---|
| LNPDB_vitro | separability_collapse | PASS |
| LNPDB_vitro | concordance_ceiling | INCONCLUSIVE |
| LNPDB_vivo | separability_collapse | PASS |
| LNPDB_vivo | concordance_ceiling | INCONCLUSIVE |

**Overall gate status: PASS**
