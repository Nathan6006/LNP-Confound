# Stage 4 — Source-Exclusivity and Excipient Invariance

## `LNPDB_vitro`

### Ionizable-lipid exclusivity

- Distinct lipids: 11,128
- Appear in exactly 1 source: 10,897 (97.9%)
- Distribution of #sources-per-lipid: {1: 10897, 2: 230, 4: 1}

### Component-level exclusivity (DB-provided head/linker/tail1/tail2 decomposition)

Using the pre-existing `IL_head/linker/tail*_SMILES` decomposition columns already in LNPDB, rather than a fresh BRICS/fragmentation pass — these are the DB curators' chemically-informed split, which is more reliable than an automatic fragmentation run cold on already-decomposed lipids. Coverage (non-missing) varies by component and is reported for each.

- **head** (`IL_head_SMILES`): coverage=67.9% (12,699/18,700 rows), 385 distinct components, 77.7% single-source
- **linker** (`IL_linker_SMILES`): coverage=14.9% (2,786/18,700 rows), 24 distinct components, 83.3% single-source
- **tail1** (`IL_tail1_SMILES`): coverage=54.1% (10,114/18,700 rows), 181 distinct components, 78.5% single-source
- **tail2** (`IL_tail2_SMILES`): coverage=17.4% (3,248/18,700 rows), 11 distinct components, 90.9% single-source


### Excipient invariance

- **helper_lipid** (`HL_name`): 7 distinct identities across 18,131 rows, top-1 = 48.5%, top-3 = 90.9%, Herfindahl = 0.339 — top values: {'DOPE': 8786, 'DSPC': 5261, 'DOTAP': 2432, 'MDOA': 572, 'DDAB': 360}
- **peg_lipid** (`PEG_name`): 9 distinct identities across 18,516 rows, top-1 = 64.6%, top-3 = 99.6%, Herfindahl = 0.496 — top values: {'DMG-PEG2000': 11954, 'DMPE-PEG2000': 4978, 'Unknown': 1512, 'ALC-0159': 12, 'DMG-C-PEG2000': 12}
- **cholesterol** (`CHL_name`): 13 distinct identities across 18,700 rows, top-1 = 99.9%, top-3 = 99.9%, Herfindahl = 0.997 — top values: {'Cholesterol': 18675, 'Beta-Sitosterol': 5, 'Campesterol': 4, 'Fucosterol': 4, 'Stigmastanol': 4}

Interpretation: if the excipient identities are highly concentrated (few distinct HL/PEG/CHL choices shared across nearly all sources — high top-1/top-3 fraction, high Herfindahl), the formulation backbone is nearly invariant across the corpus and **cannot itself explain** cross-lab disagreement in Stage 5 (there is little formulation variance to correlate with the disagreement). If instead formulation is a strong source signature, it becomes an alternative explanation for Stage 3 separability that competes with lipid-chemistry exclusivity.

### Formulation-source association (Cramer's V, NMI)

- Formulation key = HL/CHL/PEG identity + molar ratios (rounded to 0.1)
- Distinct formulation keys: 2,030
- Cramer's V(source, formulation) = **0.754** (chi2 p=0.00e+00)
- Normalized mutual information(source, formulation) = **0.639**

Figure: `figures/stage4_exclusivity_LNPDB_vitro.png` (+ .svg)

## `LNPDB_vivo`

### Ionizable-lipid exclusivity

- Distinct lipids: 1,434
- Appear in exactly 1 source: 1,432 (99.9%)
- Distribution of #sources-per-lipid: {1: 1432, 3: 1, 4: 1}

### Component-level exclusivity (DB-provided head/linker/tail1/tail2 decomposition)

Using the pre-existing `IL_head/linker/tail*_SMILES` decomposition columns already in LNPDB, rather than a fresh BRICS/fragmentation pass — these are the DB curators' chemically-informed split, which is more reliable than an automatic fragmentation run cold on already-decomposed lipids. Coverage (non-missing) varies by component and is reported for each.

- **head** (`IL_head_SMILES`): coverage=72.3% (1,242/1,717 rows), 62 distinct components, 66.1% single-source
- **linker** (`IL_linker_SMILES`): coverage=30.4% (522/1,717 rows), 24 distinct components, 100.0% single-source
- **tail1** (`IL_tail1_SMILES`): coverage=76.8% (1,319/1,717 rows), 74 distinct components, 93.2% single-source
- **tail2** (`IL_tail2_SMILES`): coverage=39.1% (672/1,717 rows), 38 distinct components, 97.4% single-source


### Excipient invariance

- **helper_lipid** (`HL_name`): 3 distinct identities across 1,717 rows, top-1 = 61.4%, top-3 = 100.0%, Herfindahl = 0.510 — top values: {'DOPE': 1055, 'DSPC': 624, 'DOTAP': 38}
- **peg_lipid** (`PEG_name`): 10 distinct identities across 1,717 rows, top-1 = 49.0%, top-3 = 99.6%, Herfindahl = 0.383 — top values: {'DMG-PEG2000': 841, 'DMPE-PEG2000': 583, 'Unknown': 286, 'DMG-PEG5000': 1, 'DPG-PEG2000': 1}
- **cholesterol** (`CHL_name`): 4 distinct identities across 1,717 rows, top-1 = 97.7%, top-3 = 99.9%, Herfindahl = 0.956 — top values: {'Cholesterol': 1678, 'n-butyl lithocholate': 35, '3-methylpentyl lithocholate': 2, 'iso-Pentyl lithocholate': 2}

Interpretation: if the excipient identities are highly concentrated (few distinct HL/PEG/CHL choices shared across nearly all sources — high top-1/top-3 fraction, high Herfindahl), the formulation backbone is nearly invariant across the corpus and **cannot itself explain** cross-lab disagreement in Stage 5 (there is little formulation variance to correlate with the disagreement). If instead formulation is a strong source signature, it becomes an alternative explanation for Stage 3 separability that competes with lipid-chemistry exclusivity.

### Formulation-source association (Cramer's V, NMI)

- Formulation key = HL/CHL/PEG identity + molar ratios (rounded to 0.1)
- Distinct formulation keys: 60
- Cramer's V(source, formulation) = **0.761** (chi2 p=0.00e+00)
- Normalized mutual information(source, formulation) = **0.778**

Figure: `figures/stage4_exclusivity_LNPDB_vivo.png` (+ .svg)
