# Stage 1 — Curation, Canonicalization, and the Frozen Analysis Set

RDKit version: **2025.03.6**

All exclusion criteria below were fixed before any downstream metric (separability, concordance, protocol comparison) was computed. Nothing here was chosen to make a result look better or worse.

## 1. Raw row counts

- `LNPDB_vitro.csv`: 19,028 rows
- `LNPDB_vivo.csv`: 1,717 rows
- `molea.csv`: 1,152 rows

## 2. SMILES canonicalization

Canonicalized with RDKit `Chem.MolFromSmiles` → `Chem.MolToSmiles(canonical=True)`. Parse failures on the primary ionizable-lipid column are excluded (Rule E1 below); failures on formulation SMILES columns are logged but do not drop rows (formulation identity is not the primary key).

| column | n present | n parse-failed | % failed |
|---|---|---|---|
| LNPDB.IL_SMILES | 20,745 | 0 | 0.0% |
| molea.smiles | 1,152 | 0 | 0.0% |
| LNPDB.IL_protonated_SMILES | 18,767 | 0 | 0.0% |
| LNPDB.HL_SMILES | 20,176 | 0 | 0.0% |
| LNPDB.CHL_SMILES | 20,745 | 0 | 0.0% |
| LNPDB.PEG_SMILES | 18,763 | 0 | 0.0% |
| LNPDB.fifthcomponent_SMILES | 38 | 0 | 0.0% |


## 3. Frozen exclusion rules applied (in order)

| rule | dataset | reason | rows affected |
|---|---|---|---|
| E1_il_smiles_parse_failure | LNPDB | IL_SMILES did not parse as a valid molecule under RDKit | 0 |
| E1_il_smiles_parse_failure | molea | smiles did not parse as a valid molecule under RDKit | 0 |
| E2_missing_delivery_readout | LNPDB | Experiment_value (delivery readout) is missing | 328 |
| E2_missing_delivery_readout | molea | Both mTP_(Log2)_Hepatocyte and mTP_(Log2)_Chondrocyte are missing | 0 |
| E3_exact_duplicate_row | LNPDB | Row is a byte-for-byte content duplicate of an earlier row (all fields incl. readout value identical) | 0 |
| E3_exact_duplicate_row | molea | Row is a byte-for-byte content duplicate of an earlier row | 0 |


## 4. Post-exclusion row counts

- LNPDB in-vitro: 19,028 → **18,700** (328 excluded)
- LNPDB in-vivo: 1,717 → **1,717** (0 excluded)
- molea: 1,152 → **1,152** (0 excluded)

## 5. Source-label and endpoint definitions

- **LNPDB**: source label = `Experiment_ID` (native, 1:1 with `Publication_PMID`). Endpoint = the row itself; no further split needed since `Experiment_ID` + `Experiment_method` + `Model_type` already uniquely identify an assay context.
- **MOLEA**: no source/publication column exists in the raw file. **DECISION (frozen, matches CLAUDE.md §6, 2026-07-23/24):** treat MOLEA as **one source, two endpoints** — `source_label='MOLEA'`, `endpoint in {hepatocyte, chondrocyte}`. This is a single combinatorial screen (`Markush code`) from one lab/library measuring two cell types, not two independent publications, so splitting it into two 'sources' would misrepresent its provenance. `mTP_Difference` is a derived selectivity contrast between the two endpoints, not an independently poolable readout, and is carried through unchanged for reference but not treated as a third endpoint.

## 6. Ionizable lipid vs. full formulation

Downstream chemical analyses (source separability, shared-lipid concordance) use `canonical_IL_SMILES` only. Formulation identity (`HL_*`, `CHL_*`, `PEG_*`, `fifthcomponent_*`, each canonicalized where present) is retained as separate columns for the Stage 4 excipient-invariance analysis and is never mixed into the lipid identity key.

## 7. Endpoint-comparable environments

An **environment** = a unique combination of (in-vitro/in-vivo, `Experiment_method` [readout modality], `Model_type` [cell line or animal strain — the biological system], `Route_of_administration`). An environment is **comparable/analyzable** if it has **≥ 20 rows** (frozen threshold, chosen for basic statistical stability, not to hit a target count). This is stricter than counting raw sources: two sources using the identical readout/cell-line/route combination are the same environment, and a source split across multiple readouts contributes multiple environments.

**LNPDB_vitro**: 34 raw sources → 14 distinct environments, of which **13 are comparable** (≥20 rows).

| environment | n_rows | comparable |
|---|---|---|
| LNPDB_vitro | luminescence_normalized | HeLa | in_vitro | 5,160 | yes |
| LNPDB_vitro | luminescence_normalized | DC2.4 | in_vitro | 3,681 | yes |
| LNPDB_vitro | luminescence_normalized | A549 | in_vitro | 2,569 | yes |
| LNPDB_vitro | luminescence_discretized_normalized | HeLa | in_vitro | 2,295 | yes |
| LNPDB_vitro | luminescence_normalized | HepG2 | in_vitro | 2,000 | yes |
| LNPDB_vitro | luminescence_normalized | RAW264.7 | in_vitro | 1,200 | yes |
| LNPDB_vitro | luminescence_discretized_normalized | IGROV1 | in_vitro | 763 | yes |
| LNPDB_vitro | luminescence_normalized | IGROV1 | in_vitro | 439 | yes |
| LNPDB_vitro | luminescence_normalized | BeWo_b30 | in_vitro | 260 | yes |
| LNPDB_vitro | luminescence_normalized | HEK293T | in_vitro | 144 | yes |
| LNPDB_vitro | LRP6_knockdown_normalized | HEK293T | in_vitro | 112 | yes |
| LNPDB_vitro | luminescence_normalized | BMDC | in_vitro | 36 | yes |
| LNPDB_vitro | luminescence_normalized | HBEC_ALI | in_vitro | 29 | yes |
| LNPDB_vitro | luminescence_normalized | BMDM | in_vitro | 12 | no |


**LNPDB_vivo**: 17 raw sources → 9 distinct environments, of which **8 are comparable** (≥20 rows).

| environment | n_rows | comparable |
|---|---|---|
| LNPDB_vivo | protein_adundance_normalized | Mouse_B6 | intravenous | 498 | yes |
| LNPDB_vivo | luminescence_normalized | Mouse_B6 | intravenous | 473 | yes |
| LNPDB_vivo | protein_adundance_normalized | Mouse_BALBc | intramuscular | 268 | yes |
| LNPDB_vivo | editing_efficiency_normalized | Mouse_Ai14 | intravenous | 141 | yes |
| LNPDB_vivo | luminescence_normalized | Mouse_B6 | intramuscular | 127 | yes |
| LNPDB_vivo | luminescence_normalized | Mouse_BALBc | intramuscular | 91 | yes |
| LNPDB_vivo | luminescence_relative_to_Spikevax | Mouse_ICR | intradermal | 60 | yes |
| LNPDB_vivo | luminescence_normalized | Mouse_B6 | intratracheal | 49 | yes |
| LNPDB_vivo | luminescence_normalized | Mouse_CD1 | intravenous | 10 | no |


## 8. MOLEA structure notes

- Distinct `Markush code` values: 1,152 (of 1,152 rows; 0 duplicate codes)
- Distinct canonical SMILES: 1,152 (0 rows share a canonical SMILES with another row)

## 9. Frozen outputs

- `results/curated.parquet`: 20,417 rows (LNPDB vitro+vivo, post-exclusion)
- `results/curated_molea.parquet`: 1,152 rows
- `results/exclusions.json`, `results/environments.csv`
