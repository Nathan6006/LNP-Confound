# Review packet — NG_2025

Generated 2026-11-09 by `provenance/review/brief.py` from the live `data/` and `curation/` trees. **Every number below was read out of the repository at generation time.** If the packet and the repository disagree, the packet is the bug — say so.

- **Paper:** Gong et al., Nat Biomed Eng 2025, 9, 2181-2195 - Mannich-reaction antioxidant ionizable lipids, 310 in vitro on MC38 + 23 in vivo (liver IV, muscle IM)
- **Curation folder:** `curation/ng_2025/`
- **Experiment_ID(s) in the corpus:** NG_2025, NG_2025_2

## 1. Corpus slice as it stands now

### in vitro — `data/data_vitro.csv`

```json
{
 "rows": 310,
 "unique_IL_SMILES": 310,
 "unique_IL_name": 310,
 "Experiment_value": {
  "n": 310,
  "min": -2.288,
  "max": 2.9436,
  "mean": "-3.226e-09",
  "pop_sd": "0.999999983983"
 },
 "Model_type": {
  "MC38": 310
 },
 "Model": {
  "in_vitro": 310
 },
 "Model_target": {
  "in_vitro": 310
 },
 "Cargo": {
  "mRNA": 310
 },
 "Cargo_type": {
  "FLuc": 310
 },
 "Experiment_method": {
  "luminescence_normalized": 310
 },
 "Route_of_administration": {
  "in_vitro": 310
 },
 "HL_name": {
  "DOPE": 310
 },
 "CHL_name": {
  "Cholesterol": 310
 },
 "PEG_name": {
  "DMG-PEG2000": 310
 },
 "confidence": {
  "medium": 282,
  "low": 28
 },
 "molratios": {
  "IL_molratio": [
   "50"
  ],
  "HL_molratio": [
   "10"
  ],
  "CHL_molratio": [
   "38.5"
  ],
  "PEG_molratio": [
   "1.5"
  ],
  "_sum_over_all_five": "100.0-100.0"
 },
 "flag_columns": {
  "value_normalization_unit": {
   "experiment": 310
  },
  "confidence": {
   "medium": 282,
   "low": 28
  },
  "positional_ambiguity": {
   "1": 310
  },
  "structure_assumed": {
   "1": 44
  },
  "curation_note": {
   "Structure is the nominal product of an unpurified crude Mannich reaction: the paper evaporated solvent to obtain 310 crude lipid compounds and screened those directly, characterising only C-a16. Value is the authors' own Source Data, not a digitisation.": 265,
   "Structure is the nominal product of an unpurified crude Mannich reaction: the paper evaporated solvent to obtain 310 crude lipid compounds and screened those directly, characterising only C-a16. Value is the authors' own Source Data, not a digitisation. Amine b3/b4 is drawn as an aminal (H2N-CH2-NR2) in the amine-list and dialkyl-ketone panels and as a hydrazine (N-N) in the pyrrole and acetophenone panels; the two differ by CH2 and the paper gives no mass, formula or NMR for any b3/b4 lipid. The aminal is emitted; the alternative is in 07_flags.csv.": 28,
   "Structure is the nominal product of an unpurified crude Mannich reaction: the paper evaporated solvent to obtain 310 crude lipid compounds and screened those directly, characterising only C-a16. Value is the authors' own Source Data, not a digitisation. Phenolic tail D is emitted as linear 4-n-nonylphenol as drawn; commercial 4-nonylphenol is predominantly a branched technical isomer mixture.": 16,
   "Structure is the nominal product of an unpurified crude Mannich reaction: the paper evaporated solvent to obtain 310 crude lipid compounds and screened those directly, characterising only C-a16. Value is the authors' own Source Data, not a digitisation. The Fig. 1 Source Data make this the highest value in the pyrrole panel (4.07e4); Supplementary Fig 17b re-measures it ~24x lower while five companion lipids on the same panel agree within 1.9x.": 1
  }
 }
}
```

### in vivo — `data/data_vivo.csv`

```json
{
 "rows": 23,
 "unique_IL_SMILES": 12,
 "unique_IL_name": 12,
 "Experiment_value": {
  "n": 23,
  "min": -2.293,
  "max": 1.4737,
  "mean": "4.348e-08",
  "pop_sd": "1.000000129333"
 },
 "Model_type": {
  "Mouse_B6": 23
 },
 "Model": {
  "in_vivo": 23
 },
 "Model_target": {
  "liver": 12,
  "muscle": 11
 },
 "Cargo": {
  "mRNA": 23
 },
 "Cargo_type": {
  "FLuc": 23
 },
 "Experiment_method": {
  "luminescence_normalized": 23
 },
 "Route_of_administration": {
  "intravenous": 12,
  "intramuscular": 11
 },
 "HL_name": {
  "DOPE": 22,
  "DSPC": 1
 },
 "CHL_name": {
  "Cholesterol": 23
 },
 "PEG_name": {
  "DMG-PEG2000": 23
 },
 "confidence": {
  "medium": 22,
  "high": 1
 },
 "molratios": {
  "IL_molratio": [
   "50"
  ],
  "HL_molratio": [
   "10"
  ],
  "CHL_molratio": [
   "38.5"
  ],
  "PEG_molratio": [
   "1.5"
  ],
  "_sum_over_all_five": "100.0-100.0"
 },
 "flag_columns": {
  "value_normalization_unit": {
   "experiment_x_model_type": 23
  },
  "confidence": {
   "medium": 22,
   "high": 1
  },
  "positional_ambiguity": {
   "1": 22
  },
  "structure_assumed": {
   "1": 6
  },
  "curation_note": {
   "Structure is the nominal product of an unpurified crude Mannich reaction: the paper evaporated solvent to obtain 310 crude lipid compounds and screened those directly, characterising only C-a16. Value is the authors' own Source Data, not a digitisation. Dose 2.5 ug is derived from the printed 0.125 mg/kg and the ~20 g body weight stated in Methods. z-scored SEPARATELY for liver (intravenous, Fig. 1m) and muscle (intramuscular injection site, Fig. 1p): different route, different tissue, different imaging session. The panels' absolute scales carry printed multipliers and so the offset is recoverable in principle, but it is not recorded in this curation, and pooling would assert that a liver IV signal and a muscle IM injection-site signal are directly comparable. Liver-vs-muscle Spearman over the 11 shared lipids is 0.982. NOTE that Experiment_value_within_model is therefore NOT what its own definition says for these rows -- it re-standardises the union of two already flattened blocks (see F-028).": 16,
   "Structure is the nominal product of an unpurified crude Mannich reaction: the paper evaporated solvent to obtain 310 crude lipid compounds and screened those directly, characterising only C-a16. Value is the authors' own Source Data, not a digitisation. Phenolic tail D is emitted as linear 4-n-nonylphenol as drawn; commercial 4-nonylphenol is predominantly a branched technical isomer mixture. Dose 2.5 ug is derived from the printed 0.125 mg/kg and the ~20 g body weight stated in Methods. z-scored SEPARATELY for liver (intravenous, Fig. 1m) and muscle (intramuscular injection site, Fig. 1p): different route, different tissue, different imaging session. The panels' absolute scales carry printed multipliers and so the offset is recoverable in principle, but it is not recorded in this curation, and pooling would assert that a liver IV signal and a muscle IM injection-site signal are directly comparable. Liver-vs-muscle Spearman over the 11 shared lipids is 0.982. NOTE that Experiment_value_within_model is therefore NOT what its own definition says for these rows -- it re-standardises the union of two already flattened blocks (see F-028).": 6,
   "Commercial benchmark, formulated with DSPC rather than DOPE exactly as Methods state. The matching Fig. 1p benchmark is deliberately held back: that panel's axis labels it DLin-MC3-DMA while the Source Data column for the same three numbers is headed SM-102, and those are different molecules. z-scored SEPARATELY for liver (intravenous, Fig. 1m) and muscle (intramuscular injection site, Fig. 1p): different route, different tissue, different imaging session. The panels' absolute scales carry printed multipliers and so the offset is recoverable in principle, but it is not recorded in this curation, and pooling would assert that a liver IV signal and a muscle IM injection-site signal are directly comparable. Liver-vs-muscle Spearman over the 11 shared lipids is 0.982. NOTE that Experiment_value_within_model is therefore NOT what its own definition says for these rows -- it re-standardises the union of two already flattened blocks (see F-028).": 1
  }
 }
}
```

### Six rows drawn at random (seed 11)

```json
[
 {
  "Experiment_ID": "NG_2025",
  "IL_name": "L-b-b13",
  "IL_SMILES": "CCCCCCCCOC(=O)C12CN(CCCN3CCCC3)CC(C(=O)OCCCCCCCC)(CN(CCCN3CCCC3)C1)C2=O",
  "Experiment_value": "-0.125818",
  "Model_type": "MC38",
  "Model": "in_vitro",
  "IL_molratio": "50",
  "HL_name": "DOPE",
  "HL_molratio": "10",
  "CHL_molratio": "38.5",
  "PEG_molratio": "1.5",
  "fifthcomponent_name": "",
  "fifthcomponent_molratio": "",
  "Publication_PMID": "40681859",
  "Publication_DOI": "10.1038/s41551-025-01422-8"
 },
 {
  "Experiment_ID": "NG_2025",
  "IL_name": "P-b6",
  "IL_SMILES": "CN1CCN(CCCN(CCC(=O)c2ccccc2)CCC(=O)c2ccccc2)CC1",
  "Experiment_value": "-0.34621",
  "Model_type": "MC38",
  "Model": "in_vitro",
  "IL_molratio": "50",
  "HL_name": "DOPE",
  "HL_molratio": "10",
  "CHL_molratio": "38.5",
  "PEG_molratio": "1.5",
  "fifthcomponent_name": "",
  "fifthcomponent_molratio": "",
  "Publication_PMID": "40681859",
  "Publication_DOI": "10.1038/s41551-025-01422-8"
 },
 {
  "Experiment_ID": "NG_2025",
  "IL_name": "M-b3",
  "IL_SMILES": "CN(C)CN(Cc1cccn1C)Cc1cccn1C",
  "Experiment_value": "-0.98887",
  "Model_type": "MC38",
  "Model": "in_vitro",
  "IL_molratio": "50",
  "HL_name": "DOPE",
  "HL_molratio": "10",
  "CHL_molratio": "38.5",
  "PEG_molratio": "1.5",
  "fifthcomponent_name": "",
  "fifthcomponent_molratio": "",
  "Publication_PMID": "40681859",
  "Publication_DOI": "10.1038/s41551-025-01422-8"
 },
 {
  "Experiment_ID": "NG_2025",
  "IL_name": "N-b10",
  "IL_SMILES": "CCN(CC)CCN(Cc1cccn1CC)Cc1cccn1CC",
  "Experiment_value": "-0.29267",
  "Model_type": "MC38",
  "Model": "in_vitro",
  "IL_molratio": "50",
  "HL_name": "DOPE",
  "HL_molratio": "10",
  "CHL_molratio": "38.5",
  "PEG_molratio": "1.5",
  "fifthcomponent_name": "",
  "fifthcomponent_molratio": "",
  "Publication_PMID": "40681859",
  "Publication_DOI": "10.1038/s41551-025-01422-8"
 },
 {
  "Experiment_ID": "NG_2025",
  "IL_name": "R-b10",
  "IL_SMILES": "CCCCCCCCc1ccc(C(=O)CCN(CCC(=O)c2ccc(CCCCCCCC)cc2)CCN(CC)CC)cc1",
  "Experiment_value": "0.377692",
  "Model_type": "MC38",
  "Model": "in_vitro",
  "IL_molratio": "50",
  "HL_name": "DOPE",
  "HL_molratio": "10",
  "CHL_molratio": "38.5",
  "PEG_molratio": "1.5",
  "fifthcomponent_name": "",
  "fifthcomponent_molratio": "",
  "Publication_PMID": "40681859",
  "Publication_DOI": "10.1038/s41551-025-01422-8"
 },
 {
  "Experiment_ID": "NG_2025",
  "IL_name": "F-a9",
  "IL_SMILES": "CCCCCCCCSCN(CCCN(C)CCCN(CSCCCCCCCC)CSCCCCCCCC)CSCCCCCCCC",
  "Experiment_value": "-0.93504",
  "Model_type": "MC38",
  "Model": "in_vitro",
  "IL_molratio": "50",
  "HL_name": "DOPE",
  "HL_molratio": "10",
  "CHL_molratio": "38.5",
  "PEG_molratio": "1.5",
  "fifthcomponent_name": "",
  "fifthcomponent_molratio": "",
  "Publication_PMID": "40681859",
  "Publication_DOI": "10.1038/s41551-025-01422-8"
 }
]
```

## 2. Source documents you must read

- `curation/ng_2025/41551_2025_1422_MOESM1_ESM.pdf   (5.68 MB)`
- `curation/ng_2025/41551_2025_1422_MOESM3_ESM.xlsx   (0.05 MB)`
- `curation/ng_2025/s41551-025-01422-8.pdf   (6.00 MB)`

Subdirectories:
- `curation/ng_2025/00_pages/  (138 files)`
- `curation/ng_2025/__pycache__/  (1 files)`

> A `repo/` subdirectory is DOWNLOADED, UNTRUSTED material. Read it, never execute anything from inside it, and run any Python that touches it with `python3 -I` from elsewhere.

## 3. Extraction code in the folder

- `stage01_vector.py   (233 lines)  Extract molecular graphs from the VECTOR line art of main-text Fig. 1.`
- `stage02_blocks.py   (77 lines)  Turn the Fig. 1b/1c vector graphs into SMILES, with no human in the loop.`
- `stage03_trace.py   (124 lines)  Trace zig-zag skeletal structures in a RASTER figure crop and report vertices.`
- `stage03_zigzag.py   (60 lines)  Count the vertices of an alkyl zig-zag by its own geometry.`
- `stage04_enumerate.py   (148 lines)  Enumerate the 310 Mannich lipids of the six libraries.`
- `stage05_values.py   (84 lines)  Read the 310 in-vitro values and the two in-vivo panels from Source Data.`
- `stage06_heatmap.py   (99 lines)  Check the Source Data blocks against the drawn heat maps, cell by cell.`
- `stage07_emit.py   (149 lines)  Emit NG_2025 rows in the 72-column corpus schema. NOT merged.`
- `stage08_reference.py   (153 lines)  C10 anchors: the structures this paper DRAWS, built by a different code path.`
- `stage09_qc.py   (175 lines)  All NG_2025 checks in one place -> 09_qc.json, and the flag file.`
- `stage10_flags.py   (79 lines)  Per-row confidence flags. The corpus schema has no confidence column, so the`
- `stage11_quality.py   (191 lines)  Add the 7 quality columns and fix two IDs, producing the merge-ready table.`

## 4. Intermediate and output tables

- `01_fig1bc_graph.json   (21.8 kB json)`
- `01_fig1d_graph.json   (5.3 kB json)`
- `02_blocks_hand.csv   33 rows x 3 cols
        cols: block_id, SMILES, source`
- `02_vector_vs_hand.csv   20 rows x 7 cols
        cols: block_id, n_atoms_drawn, x, y, vector_SMILES, hand_SMILES, agree`
- `03_enumerated.csv   310 rows x 8 cols
        cols: lipid_name, heatmap_label, library, tail, amine, n_tail_arms, IL_SMILES, source_figure`
- `03_ketone_p09.json   (6.2 kB json)`
- `03_thiol_tails.json   (2.2 kB json)`
- `04_values.json   (41.4 kB json)`
- `06_heatmap_check.json   (1.0 kB json)`
- `07_NG_2025_rows.csv   333 rows x 72 cols
        cols: Index, LNP_ID, Experiment_ID, Formulation_ID, IL_name, IL_SMILES, IL_protonated_SMILES, IL_head_name, IL_head_SMILES, IL_linker_name, IL_linker_SMILES, IL_tail1_name, IL_tail1_SMILES, IL_tail2_name, IL_tail2_SMILES, IL_tail3_name, IL_tail3_SMILES, IL_tail4_name ...`
- `07_flags.csv   133 rows x 5 cols
        cols: lipid_name, confidence, issue, detail, alternative_amine_SMILES`
- `07_held_back.csv   1 rows x 3 cols
        cols: panel, label, reason`
- `08_NG_2025_merge.csv   333 rows x 79 cols
        cols: Index, LNP_ID, Experiment_ID, Formulation_ID, IL_name, IL_SMILES, IL_protonated_SMILES, IL_head_name, IL_head_SMILES, IL_linker_name, IL_linker_SMILES, IL_tail1_name, IL_tail1_SMILES, IL_tail2_name, IL_tail2_SMILES, IL_tail3_name, IL_tail3_SMILES, IL_tail4_name ...`
- `09_qc.json   (5.0 kB json)`
- `NG_2025_blocks.csv   55 rows x 16 cols
        cols: block_id, block_class, name_in_paper, SMILES, n_reactive_sites, reactive_site_type, MW_reported, formula_reported, HRMS_reported, CAS_or_vendor, source_figure, source_table, source_page, read_method, read_confidence, ambiguity_note`
- `NG_2025_checks.csv   25 rows x 9 cols
        cols: check_id, check_description, what_it_would_catch, method, expected, observed, result, severity, note`
- `NG_2025_conditions.csv   9 rows x 45 cols
        cols: condition_id, block_ids, naming_rule, degree_of_substitution, reaction_description, library_size_stated, cargo_type, cargo_identity, formulation_molar_ratios_verbatim, IL_mol_pct, helper_lipid_identity, helper_mol_pct, sterol_identity, sterol_mol_pct, peg_lipid_identity, peg_mol_pct, IL_to_nucleic_acid_wt_ratio, N_to_P_ratio ...`
- `reference_structures.csv   6 rows x 7 cols
        cols: lipid_name, reference_SMILES, source_figure, read_on, read_by, expected_mismatch, note`

- `README.md   (11.9 kB)`

## 5. C10 drawn-structure anchors for this folder

```json
[
 {
  "paper": "ng_2025",
  "lipid_name": "C-a16",
  "source_figure": "Supplementary Fig 3a; Supplementary Fig 19a (MS inset)",
  "read_on": "2026-10-29",
  "note": "four 2-hydroxy-5-octylbenzyl arms on the two lysine side-chain NH2 of cyclo(Lys-Lys); the diketopiperazine N-H stay intact",
  "reference_inchikey": "AIOBALMNQABGKF-UHFFFAOYSA-N",
  "enumerated_inchikey": "AIOBALMNQABGKF-UHFFFAOYSA-N",
  "skeleton_match": true,
  "stereo_match": true,
  "status": "exact_match",
  "anchor_may_be_a_copy": "stored reference is byte-identical to the enumeration and already canonical, and this folder has no make_reference.py showing a separate derivation. An anchor copied from the enumeration can only pass (F-013). Byte-equality is not proof of a copy -- two correct routes give the same canonical SMILES -- so this is a prompt to confirm, not a failure."
 },
 {
  "paper": "ng_2025",
  "lipid_name": "M-b1",
  "source_figure": "Supplementary Fig 14a",
  "read_on": "2026-10-29",
  "note": "N-methyl, two (1-methylpyrrol-2-yl)methyl arms",
  "reference_inchikey": "XETWGVUCRKHXPQ-UHFFFAOYSA-N",
  "enumerated_inchikey": "XETWGVUCRKHXPQ-UHFFFAOYSA-N",
  "skeleton_match": true,
  "stereo_match": true,
  "status": "exact_match",
  "anchor_may_be_a_copy": "stored reference is byte-identical to the enumeration and already canonical, and this folder has no make_reference.py showing a separate derivation. An anchor copied from the enumeration can only pass (F-013). Byte-equality is not proof of a copy -- two correct routes give the same canonical SMILES -- so this is a prompt to confirm, not a failure."
 },
 {
  "paper": "ng_2025",
  "lipid_name": "P-b1",
  "source_figure": "Supplementary Fig 16a",
  "read_on": "2026-10-29",
  "note": "N-methyl, two 3-oxo-3-phenylpropyl arms (the acetophenone methyl becomes the CH2 beta to nitrogen)",
  "reference_inchikey": "FCEMMQWGASQRTL-UHFFFAOYSA-N",
  "enumerated_inchikey": "FCEMMQWGASQRTL-UHFFFAOYSA-N",
  "skeleton_match": true,
  "stereo_match": true,
  "status": "exact_match",
  "anchor_may_be_a_copy": "stored reference is byte-identical to the enumeration and already canonical, and this folder has no make_reference.py showing a separate derivation. An anchor copied from the enumeration can only pass (F-013). Byte-equality is not proof of a copy -- two correct routes give the same canonical SMILES -- so this is a prompt to confirm, not a failure."
 },
 {
  "paper": "ng_2025",
  "lipid_name": "I-a-b1",
  "source_figure": "Supplementary Fig 8a",
  "read_on": "2026-10-29",
  "note": "diethyl 1-methyl-4-oxopiperidine-3,5-dicarboxylate: the 1:2:1 double Mannich closes a ring, it is not an acyclic bis-adduct",
  "reference_inchikey": "FMZCQFLLEXZUJM-UHFFFAOYSA-N",
  "enumerated_inchikey": "FMZCQFLLEXZUJM-UHFFFAOYSA-N",
  "skeleton_match": true,
  "stereo_match": true,
  "status": "exact_match",
  "anchor_may_be_a_copy": "stored reference is byte-identical to the enumeration and already canonical, and this folder has no make_reference.py showing a separate derivation. An anchor copied from the enumeration can only pass (F-013). Byte-equality is not proof of a copy -- two correct routes give the same canonical SMILES -- so this is a prompt to confirm, not a failure."
 },
 {
  "paper": "ng_2025",
  "lipid_name": "I-b-b1",
  "source_figure": "Supplementary Fig 11a",
  "read_on": "2026-10-29",
  "note": "diethyl 3,7-dimethyl-9-oxo-3,7-diazabicyclo[3.3.1]nonane-1,5-dicarboxylate (bispidinone): the 1:4:2 ratio closes a second ring",
  "reference_inchikey": "CKUKECUYDJFXEN-UHFFFAOYSA-N",
  "enumerated_inchikey": "CKUKECUYDJFXEN-UHFFFAOYSA-N",
  "skeleton_match": true,
  "stereo_match": true,
  "status": "exact_match",
  "anchor_may_be_a_copy": "stored reference is byte-identical to the enumeration and already canonical, and this folder has no make_reference.py showing a separate derivation. An anchor copied from the enumeration can only pass (F-013). Byte-equality is not proof of a copy -- two correct routes give the same canonical SMILES -- so this is a prompt to confirm, not a failure."
 },
 {
  "paper": "ng_2025",
  "lipid_name": "B-a1",
  "source_figure": "Fig. 1d",
  "read_on": "2026-10-29",
  "note": "DERIVED at run time from the page vector bond coordinates in 01_fig1d_graph.json (stage01_vector.py, glyph merge radius 3.0): glyph signatures give the atoms, geometry gives the chain lengths. No human typed this string, so it cannot have been copied from the enumeration",
  "reference_inchikey": "YOZSZGYTICTJMX-UHFFFAOYSA-N",
  "enumerated_inchikey": "YOZSZGYTICTJMX-UHFFFAOYSA-N",
  "skeleton_match": true,
  "stereo_match": true,
  "status": "exact_match",
  "anchor_may_be_a_copy": "stored reference is byte-identical to the enumeration and already canonical, and this folder has no make_reference.py showing a separate derivation. An anchor copied from the enumeration can only pass (F-013). Byte-equality is not proof of a copy -- two correct routes give the same canonical SMILES -- so this is a prompt to confirm, not a failure."
 }
]
```

> An anchor can only be trusted if it was assembled by a DIFFERENT code path than the enumerator. Check that, do not assume it.

## 6. C11 named commercial reagents appearing in this slice

```json
{
 "DOPE (HL_name)": 310,
 "Cholesterol (CHL_name)": 310
}
```

> C11 checks these by InChIKey. It CANNOT check a name that is not in the registry, and it cannot tell you whether the paper's own shorthand means the lipid the registry assumes. `18PG` is the precedent: it was read as symmetric 18:0/18:0 across 360 rows and is in fact 18:0/18:1.

## 7. Folders that are NOT yours

- `curation/sd_2026/` — queued, not started - PDFs only
- `curation/jl_2026/` — PAUSED pending the authors' reply about source data; the SI heatmaps are 192 dpi and the printed digits are 5.1x3.5 px, which we established is not extractable
- `curation/liu_iphos/` — blind-validation exercise against the existing SL_2021 rows, not a merge - reviewed under corpus target CORPUS-C
- `curation/amino_lipid_fg_2026/` — blind-validation exercise against existing FG_2026 rows, not a merge - reviewed under CORPUS-C

## 8. Error-class registry (verbatim)

The second table, *Known, not yet automatable*, is where your time is worth most: Tier-1 already covers the first table on every build.

# Error-class registry

The ratchet. Every bug any reviewer finds gets a row here and, wherever possible,
a permanent automated check. **This file only grows.** A class with a check can
never silently recur; a class without one is where review budget should go.

The point of keeping it: reviewers are expensive and find a class *once*; an
automated invariant is free and prevents it *forever*. Discovering a class is
worth paying for. Re-discovering one is not.

## Covered by an automated check

| # | Error class | Found | Check | Where |
|---|---|---|---|---|
| 1 | SMILES that does not parse | early | C1 | `audit.py` |
| 2 | Non-canonical SMILES (same molecule, two strings) | early | C1 | `audit.py` |
| 3 | Disconnected fragments / stray counter-ion | early | C1 | `audit.py` |
| 4 | **Radical electron from a deleted substituent** | 2026-10-29, Siegwart | C1 | `audit.py` + hard guard in both `quaternize()` implementations |
| 5 | Descriptor columns stale vs their own SMILES | 2026-10-28 | C8 | `audit.py` |
| 6 | Mass non-monotonic within a homologous series | 2026-10-28 | C2 | `audit.py` |
| 7 | Duplicate SMILES under two names in one experiment | 2026-10-28 | C3 | `audit.py` |
| 8 | Exact duplicate rows | 2026-10-28 | C4 | `audit.py` |
| 9 | Formulation ratios not summing to 100 | 2026-10-28 | C5 | `audit.py` |
| 10 | Non-finite or absurd `Experiment_value` | 2026-10-28 | C6 | `audit.py` |
| 11 | Field inconsistency; blank required fields | 2026-10-28 | C7 | `audit.py` |
| 12 | Lost E/Z on a stereogenic double bond | 2026-10-28 | C9 | `audit.py` |
| 13 | z-score convention drift (ddof, post-hoc row drops) | 2026-10-28 | rebuild | `build_dataset.py` |
| 14 | **Attachment at the wrong atom** (acyl on a side-chain methyl, not C3) | 2026-10-29, Siegwart | **C10** | `structure_check.py` |
| 15 | **Substituent on the wrong nitrogen** (regiochemistry) | 2026-10-28 Lin, 2026-10-29 Siegwart | **C10** | `structure_check.py` |
| 16 | Heatmap geometry off by a row/column (coloured header band) | 2026-10-28, twice | colour-distance gate | each `digitize_*.py`; >~5 means geometry, not colour |
| 17 | Wrong publication identifier (PMID ≠ DOI ≠ link) | 2026-10-29 | manifest build | `build_manifest.py` resolves every PMID/DOI against NCBI |
| 18 | **The PAPER draws a known benchmark wrong** (ALC-0315 head one carbon short) | 2026-10-29, LipidAI | **C10** `expected_mismatch` | `structure_check.py` — the override gets its own status so it stays visible instead of living in a prose note |
| 19 | **Chain length miscounted by eye from a figure** | standing risk, acted on 2026-10-29 | geometric vertex trace | `curation/xd_2026_lipidai/stage02_trace.py` — skeletonise, trace, RDP-simplify, count; never count by eye |
| 20 | **Log axis anchored on the panel box instead of its tick marks** (10³ tick sat 47 px above the bottom rule; every value high by up to 1.76×) | 2026-10-29, `ql_2026` Fig 3c | per-curation | reproduce every printed summary statistic, never a subset |
| 21 | **Provenance text asserts one thing, the code does another** (docs said the heatmap value was carried; the code took the conflicting bar-chart value) | 2026-10-29, `ql_2026` | per-curation | a check that inverts the stored value back to its source |
| 22 | **Blocks standardised separately when they should be pooled**, deleting a real between-block difference | 2026-10-29, `ql_2026` | per-curation | assert the between-block gap SURVIVES; per-block mean 0/sd 1 is a tautology |
| 23 | **Colormap floor/ceiling degeneracy left unflagged** (113 of 320 cells in a band where identical RGB spans ~980 RLU) | 2026-10-29, `ql_2026` Fig 2b | per-curation | flag floored and saturated cells symmetrically |
| 24 | **A C10 anchor documented as machine-derived but actually a typed literal** (the note and the README both said the B-a1 reference was read from vector bond coordinates; it was a string in the script, and nothing in any stage read the graph file) | 2026-10-30, `ng_2025` | per-curation | derive the anchor IN the script from the stored graph, so the provenance claim is executed rather than asserted; a claim no code runs is the same failure as copying the anchor from the enumeration |
| 25 | **Per-item drawings contradicting the paper's own reaction scheme** (about a third of the 104 dialkyl-ketone lipids are drawn as the open-chain mono-adduct, and every one of those also gives the amine nitrogen an extra methyl it cannot have) | 2026-10-30, `ng_2025` | per-curation | classify EVERY drawn member against the scheme and report the count; a spot-check gave 5 where the true figure was 15, understating the ambiguity six-fold |
| 26 | **A reconstructed rule validated on a subset but reported library-wide** (the >5x-MC3 hit threshold reproduces the 11 in-vivo lipids among the 64 phenols, but over all 310 it also catches a thiol lipid that was never taken into mice) | 2026-10-30, `ng_2025` | per-curation | run the reconstructed rule over the FULL set, and state the scope in the claim itself |
| 27 | **A validation calibrated on the quantity it validates** (`digitize_fig2fgh.py` fitted its pixels-per-percent scale by least squares onto the Fig 2e hit rates the panel was checking; the measured tick spacing was computed and never used, and the file had no `assert`) | 2026-10-30, `LL_2026` | **none possible in `audit.py`** — it is a property of a curation's validation code, not of the corpus | review rule: the pixel-to-unit scale must derive ONLY from axis furniture (ticks, labelled gridlines, printed anchors), never from the comparison target. Signature to grep for: a fit whose x-vector is the expected answer |
| 28 | **A measurement failure disguised as a measurement** (`else 0.0` fallback recorded a bar that exists as 0.0 because a `min_area` floor dropped its 450 px² component — and that fabricated zero became the set's single largest "disagreement") | 2026-10-30, `LL_2026` Fig 2h EP6 | per-curation | never let a detector miss fall through to a numeric literal: record `None` and fail, so an absent measurement cannot be averaged, deltaed or reported as agreement |

| 29 | **A hard-coded confidence literal discarding the per-item assessment** (`build_new_curation.py` read `lipid_name -> IL_SMILES` from the enumeration and then wrote `"high"` for every lipid it found, throwing away the `confidence` and `quaternization_site_assumed` columns beside it — so the 24 lipids whose record reads *"ASSUMED: no ring N and the paper draws no structure"* were labelled the project's highest confidence) | 2026-10-30, `ZT_2026` | **C13** | `audit.py` — plus: never write a confidence a per-item record already states; carry it. Inverting a disclosed assumption into an assurance is worse than publishing nothing, because the only surviving provenance then asserts the opposite of the truth |
| 30 | **A component looked up in the live corpus BY NAME, resolved by file order** (`comp.setdefault(name, smiles)` = first row wins, over a corpus written in alphabetical experiment order; 292 names map to more than one structure. `DLin-MC3-DMA` was *ambiguous* — 2,190 rows of the correct ester plus 203 dioxolane rows from `AC_2025` — so the lookup silently returned whichever experiment sorted first, and `MC3-CB1/2/3` were built on KC2) | 2026-10-30, `ZT_2026` | per-curation | the lookup must RAISE on an ambiguous name instead of picking one (`load_corpus_components` now returns the ambiguity set and `smi()` refuses). Regenerating after a `NAME_PATCHES` change fixes the symptom by luck of row order, not the cause; nothing tracks which artefacts were built from a corpus lookup |
| 31 | **Fragment assembly silently inverting a chiral tag** (`RWMol.RemoveAtom` deleted cholesterol's hydroxyl oxygen out of the middle of C3's neighbour list without correcting C3's tag, so all 537 `ZT_2026` rows stored epicholesterol, 3α. The epimers share molecular formula, exact mass to four decimals, the InChIKey SKELETON block and all 17 descriptor columns, so C1–C9 were blind; C10 was blind because its anchor was a byte copy of this code's own output; C11 was blind twice over, comparing skeleton blocks only and inspecting `CHL_*` columns that are empty on every one of these rows) | 2026-10-30, `ZT_2026` | **C14** | `audit.py` — cleave every ester inside `IL_SMILES` corpus-wide and pin the liberated sterol's FULL InChIKey against `reagent_registry.json`. Also: do the condensation the way the chemistry does (the ACID loses its OH, the alcohol's oxygen is retained), so no atom adjacent to a stereocentre is ever deleted |
| 32 | **A bar-top detector fooled by replicate markers** (replicate dots straddling the top edge occlude the MIDDLE of the bar width, so a "first row where >60% of the width is fill" rule fires several pixels INSIDE the bar; the bias is one-sided and non-constant, so it does not cancel in a z-score) | 2026-10-30, independently in `QL_2026` Fig 2c/3c and `XD_2026` Fig 5C **on the same day** | per-curation | measure from the FLANK columns, which a centred marker cannot reach, and require the left and right flank to agree to ≤1–2 px. Note the limit: every printed hit rate in the `QL_2026` panels still reproduces on the corrected values, so the paper's own statistics cannot discriminate this one — it has to be caught geometrically |
| 33 | **A check whose findings are printed but never gate** (`audit.py` discarded `main()`'s return value; C1–C9 printed and the exit code came only from the C10/C11 delegates, so eight simultaneous seeded violations exited 0. Of 43 mutations, 30 appeared in the report and 5 failed a command. `build_dataset.py` never imported `sys` at all) | 2026-10-30, corpus-wide | **the gate itself** | `audit_baseline.json` + `test_checks.py` — each metric has an accepted value WITH its reason, departure in EITHER direction fails (a vanished known condition is how a check that stopped looking announces itself), and a seeded-mutation suite proves the gate fires. The obvious fix is wrong on its own: failing on anything printed fails forever on conditions we accept on purpose, which is how it became print-only |
| 34 | **A quality flag produced, documented, and then dropped at the schema boundary** (1,137 rows carried a convention rather than a measurement; `CLAUDE.md` claimed the project flagged them "in the data (`positional_ambiguity`, `*_assumed`, `confidence`)" and none of those columns existed anywhere. The loss was one stage EARLIER than the merge — no file `build_dataset.py` reads had a flag column — so no change to the builder could have preserved them) | 2026-10-30, four of five curations | **C13** | `schema_quality.json` as the single source of the column names and their vocabulary; `audit.py` enforces the vocabulary and counts rows flagged without a note. A flag that lives only in a per-folder QC JSON is not a flag, it is a memo |

| 35 | **A mask change silently invalidating every window tuned against the old mask** (switching Fig 3c's bar detector from a blue-hue mask to a luminance mask made the panel box's own top rule visible at y = 2301-2314, and the search window started at 2311 — inside it — so every bar would have returned the same wrong top) | 2026-10-30, `QL_2026` Fig 3c | per-curation | assert the first scanned row is not already a full-width rule. Generally: a mask defines what counts as ink, so every threshold, window and offset calibrated under the old one has to be re-derived, not inherited |
| 36 | **A calibration measured mostly from our own detector error** (the `QL_2026` cross-library offset was computed from four benchmark LNPs plotted in both panels — and three of the four were dot-contaminated, so the 1.1116 "offset between the two screens" was largely the bar-top bug measuring itself. Correcting the detector did not sharpen the offset, it removed it: 1.0071 ± 0.0247, not distinguishable from 1.0) | 2026-10-30, `QL_2026` | per-curation | a quantity derived from a small set of shared anchors inherits every bias those anchors carry; check the anchors for the defect under investigation BEFORE trusting a constant derived from them |
| 37 | **A check made blind to a class by the normalization it performs before testing** (915 `IL_name` cells carried a trailing space from upstream. C7's required-field check called `.strip()` on each value *before* asking whether it was blank, so a padded value passed as present and the check could not see padding at all. `"A1I5R22C8 "` and `"A1I5R22C8"` are distinct dict keys carrying the *same* SMILES, so this is a silent join hazard, not a cosmetic one — and the condition was DISCLOSED in `CLAUDE.md` and worked around with two defensive `.strip()` calls in `build_dataset.py` for the corpus's whole life) | 2026-11-07, `BL_2023` + `ZH_2023` | **C7 `padded_cells`** | strip at the single emit choke point, and test for padding on EVERY cell rather than normalizing it away inside the test. Generally: a check that canonicalises its input before comparing has deliberately given up the ability to detect non-canonical input, so whatever it normalises away needs its own check. A known-and-worked-around defect is worse than an unknown one, because the workarounds are evidence everybody saw it and nobody owned it |
| 38 | **A tracked artifact whose input lives in a session scratchpad** (`build_manifest.py` read its NCBI cross-check cache from `.../claude-501/<old-session-uuid>/scratchpad/esum.json`. Scratchpads are session-scoped and get cleared, so the script raised `FileNotFoundError` and `data/experiment_manifest.csv` — committed to git — could not be regenerated from anything in the repository. The cache was also silently STALE: it predated the `QL_2026` and `XD_2026` merges, so the manifest recorded "PMID … not checked against NCBI" for both) | 2026-11-07, second instance — `lit_review/build_data.py` had it on 2026-10-30 | the rule, plus `fetch_esummary.py` | **no script in this repo may read an input from a session scratchpad**; anything a tracked output depends on is written into the repo together with the script that fetches it. Derive paths from `__file__`, not from an absolute path under `/Users/`. Two instances make this a missing rule rather than an oversight |
| 39 | **A number asserted in the project's own documentation that nothing measures** (`CLAUDE.md` stated the F-001 flag count as 1,152 and told the reader to re-measure it with `C13.rows_carrying_a_flag`, which returns 1,652 because it counts `confidence = high`; the figure matching the prose is 1,157. It also gave the chance baseline as `1/31` vitro — the S from before a July re-trim, matching neither the live corpus (42) nor the frozen run (34) — and described a trailing-space defect as affecting "several" rows when it was 915 cells. The count had moved four times, ≥266 → 1,137 → 1,152 → 1,157, each correction by hand) | 2026-11-07, corpus-wide (F-022/F-024/F-026) | **C15** | `provenance/check_docs.py` + `docs_claims.json` — every documented number paired with an expression that measures it, each tagged with the section that claims it. Deliberately does NOT parse the Markdown: the expected value is transcribed by hand, so changing a doc number means changing both and the pair either agrees or the build fails. **A quantity with a prose definition and no implementation gets re-implemented, slightly differently, by every pass that touches it.** Narrower cousin of class 21, which is a curation's provenance text contradicting its own code |
| 40 | **A source stating the same quantity in two parameterisations where only ONE is row-aligned** (`WT_2026`'s URS file gives each formulation as four design parameters AND four molar percents. On that one file of seven they describe different formulations: 0 of 129 design rows agree with their own percent row, while 126 derive exactly one OTHER row's percents, injectively. The design block's `PEG_PEG+Chol` is strictly monotone over all 129 rows with zero violations -- impossible under uniform random sampling, and the signature of a spreadsheet sort applied to a contiguous column selection. `IL_to_nucleicacid_chargeratio` was taken from `NP_ratio`, which travels with the sorted block, so it is wrong on 129 corpus rows: median error 6.86 against a file range of 4.01-25.54. C4, C5 and C6 are all blind because every row's percents still sum to exactly 100) | 2026-11-07, `WT_2026` | per-curation, automatable | where a source gives both parameterisations, DERIVE one from the other and require **row-wise** agreement. The positive signature of the defect is a cross-row injective match with ZERO self-matches. Note what does not work: a range check passes identically on misassigned values because a permutation preserves the multiset, which is why this curation's `validate.py` recorded `PASS N/P range 3.33-37.12` on 129 wrong rows |
| 41 | **A derived "within-model" column silently pooling normalisation units the source says are incomparable** (`Experiment_value_within_model` re-standardises inside (`Experiment_ID` x `Model_type`), but `WT_2026` is normalised per SOURCE FILE and one `Model_type` spans several files -- `Ramos` alone covers 959 rows across five. So the derived column straddles five scales the authors' own `dataset_info.txt` forbids comparing, and `value_normalization_unit` does not describe it: that column declares the unit of `Experiment_value` only. A consumer reaching for the "within model" column to remove an assay offset gets exactly the cross-file mixing the `source_file` flag exists to prevent) | 2026-11-07, `WT_2026` | **open** | a derived column needs its own declared unit, or must refuse to compute where the base unit is not the one it assumes. Generally: adding a re-normalised convenience column creates a SECOND normalisation contract, and a flag that describes only the first one is not a flag on the second |
| 42 | **A quantity reported to a precision the measurement cannot reach, and dominated by the reading convention used to measure it** (the `QL_2026` cross-library offset. It went 1.1116 -> 1.0071 as the bar-top detector was repaired and an independent re-measurement returned 0.9929, the RECIPROCAL -- two honest measurements differing by more than the correction either applies. 0.71% is 0.2 of one native pixel. And both panels are read at the TOP EDGE of the bar's outline stroke while the value sits at the stroke CENTRE, so the pen-width difference between the panels (10 px vs 6 px top rule) biases the ratio +1.8-3.1% on its own -- 2.5-4.4x the offset measured through it) | 2026-11-09, `QL_2026` | per-curation | before quoting a ratio of two digitized quantities, compute what ONE pixel is worth and what the reading convention contributes. If either exceeds the effect, report the comparison as an AGREEMENT INTERVAL and apply no correction. The check keeps its value -- it would have caught a real 1.5x offset -- while the number does not |
| 43 | **A tie in a lookup table resolved by index order, biasing every value in the degenerate band one way** (`digitize_fig2b.py` matched each heatmap cell to the nearest colorbar row with `argmin`. Near the floor several adjacent rows carry IDENTICAL RGB, so every distance ties and `argmin` returns the FIRST tie = smallest y = LARGEST value. Band widths grow toward the floor, so the bias is one-sided AND non-uniform -- median +12.8%, up to +17.9% over 186 cells -- and a non-uniform one-sided bias does not cancel in a z-score) | 2026-11-09, `QL_2026` Fig 2b | per-curation | take the MIDPOINT of the tied run, which is what the colour actually says. Keep the band EDGE only where the band is open on one side and the row is flagged as a bound, and check the flag's threshold separates the open band from the closed ones -- here a single `< 4.0` colour test straddled both and put `value_is_bound = upper` on 65 two-sided measurements |
| 44 | **A known class repaired only in the file where it was found** (error class 30 -- component lookup by name over the corpus, first row wins -- was hardened in `build_new_curation.py` in round 1 and left live in THREE more emitters, one of which resolves `DLin-MC3-DMA`, the very name whose ambiguity put 419 measurements on the wrong molecule. It is unambiguous today only because an unrelated name repair landed, i.e. fixed by luck of row order rather than by a guard. The same round also found a retracted count surviving in the one file the retraction never reached, and a repaired anchor re-copying the half its own bug had not touched) | 2026-11-09, corpus-wide (F-037/F-044) | **shared module + C16** | when a class is fixed, GREP THE CORPUS FOR THE CONSTRUCT, not for the symptom, and put the remedy in one importable place so the next caller inherits it. A per-file repair of a corpus-wide class is a sample, not a fix |
| 45 | **A check keyed on a column that a rename moved, defaulting to a number instead of raising** (three of four `ZT_2026` Fig 2c marginal checks were keyed on `*_molpct` against keys that are raw design levels; the 2026-10-28 rename moved those levels to `*_molratio`, every lookup missed, and `got.get(k, 0)` scored a literal 0. The QC file recorded the three failures verbatim while the summary beside it asserted ALL FOUR EXACT, for twelve days) | 2026-11-09, `ZT_2026` | per-curation | index with `got[k]` so a miss raises, and assert the key column's actual contents in the message. Classes 28 and 35 compounding: a default turned a miss into a datum, and a rename invalidated everything keyed on the old name |
| 46 | **A reagent whose structure already exists, independently curated, elsewhere in the corpus — and was never consulted** (`ZT_2026` read the `OF-2` core off a figure by eye and wrote all four tails one carbon short, C17 for C18. `LL_2025` had carried the same reagent as `OF-02` all along, at the correct C84H152N4O6 / UMSCHHWDFWCONK-RGLFHLPNSA-N. The error was internally consistent — formula, descriptors and canonical SMILES all agreed with each other — and this was the one lipid whose own `confidence_basis` records that the apex chain counter could NOT be used on it, because its alkenes break the zig-zag, so a by-eye count was the only evidence it ever had) | 2026-11-07, `ZT_2026` `OF-2` (F-031) | **not yet automated** | for every reagent name appearing in more than one experiment, compare skeleton InChIKeys and report disagreements. Needs a small alias table (`OF-2` vs `OF-02`) and is otherwise mechanical. The general point is stronger than the check: **the corpus is its own best external anchor, and a second independent curation of the same reagent is the cheapest C10 there is** — cheaper than reading the figure again, and available for free whenever two papers share a commercial lipid |
| 47 | **Two reagent names swapped** (MC3 ↔ KC2, identical formula and exact mass) | 2026-10-29, `AC_2025` | **C11** | `audit.py` — named-reagent registry | *(was numbered 18; renumbered 2026-11-09, see the note on numbering)*
| 48 | **Vendor shorthand expanded to the wrong lipid** (`18PG` read as symmetric 18:0/18:0; it is 18:0/18:1 SOPG) | 2026-10-29, `YZ_2022`/`YZ_2024` | **C11** | `audit.py` — named-reagent registry | *(was numbered 19; renumbered 2026-11-09, see the note on numbering)*

## A note on the numbering above

**Two numbers were used twice, from 2026-10-29 until 2026-11-09.** Classes 18 and
19 each named two different things:

| was | is now | class |
|---|---|---|
| 18 | **18** | the PAPER draws a known benchmark wrong (ALC-0315) |
| 18 | **47** | two reagent NAMES swapped (MC3 / KC2) |
| 19 | **19** | chain length miscounted by eye from a figure |
| 19 | **48** | vendor shorthand expanded to the wrong lipid (`18PG`) |

That is not cosmetic, and it had already done its damage: two separate documents
cite "error class 19" meaning **different classes** — the round-2 ledger means
the by-eye chain count, and `reagent_registry.py` meant the `18PG` expansion.
A ledger whose identifiers are not unique cannot be cited, and the whole point
of this file is that a class costs tokens to discover once and ~0 to prevent
afterwards — which only holds if a reference resolves.

Renumbered 2026-11-09, keeping the numbers that were already cited correctly and
moving the later pair to the end so no existing citation silently repoints. If
you add a class, take `max + 1` and check the number is unused first.


Classes **18** and **19** each appear twice: the figure-drawing/vertex-trace pair
and the reagent-name/vendor-shorthand pair were added on the same day by
different passes. **Still duplicated, deliberately.** Round 1's reports are
written against these numbers, so renumbering now would make every reference in
`review/round_2026-10-30/` point at the wrong row. Renumber when round 2 opens,
in one pass, and retire the old numbers rather than reusing them.

**This file is edited concurrently by more than one session**, which is how the
duplication happened and how classes 27-28 below were first written as 24-25,
colliding with 24-26 added from another session minutes earlier. Before
appending, re-read the table and take the next free number; do not assume the
highest number you last saw is still the highest.

## Known, not yet automatable

These need eyes on a figure or judgement about a paper's text. This is where a
reviewer's time is actually worth spending.

| Error class | Why no check | Mitigation |
|---|---|---|
| Name token misread (`C24` = two C4 tails, not 24 carbons) | needs the library's own naming convention | C2 restricted to conventions verified against source metadata |
| Value misread from a figure | needs the figure | overlay PNG per extraction; spot-check ≥10 cells |
| Experimental metadata contradicting itself between caption and Methods | needs both read in full | record BOTH and flag UNRESOLVED (e.g. Siegwart 20 ng/n=3 vs 15 ng/n=5) |
| Commercial reagent that is an isomer mixture | no single true structure | flag `confidence = medium` and say so |
| Regiochemistry a paper never pins | genuinely unknowable | flag `positional_ambiguity` / `*_assumed`, confidence `low` |
| **A C10 anchor copied from the enumeration** (can only pass) | once both sides are canonicalised a copied reference is byte-identical to an honestly re-derived one, so a string comparison would fail the curations that did it right | assemble the reference by a different code path — fragment assembly, or just a different atom ordering — and say which in `read_by` / `note` |
| **De-duplicating reference strings on an assumed identifier format** | LipidAI appended a per-row counter to the DOI *sometimes*, and to a PMCID or a page range otherwise; cutting only at `doi:` over-counts badly. The class is real; the NUMBERS first recorded here were not. "299 against a true 47" was a reviewer claim propagated unchecked; the folder's own `06_lipiddb_qc.json` records 56/50 and the paper states 47, and nine defensible de-duplication rules give 50-116. | cut at the first identifier of ANY kind, cross-check against the count the paper states itself, and **report the rule with the number** - a de-duplicated count is meaningless without it |
| **A structure with no external anchor at all** | commercial lipids supplied without name, formula, HRMS or synthesis (LipidAI's E1–E9) | emit with `confidence = medium`, record an empty `reference_SMILES` row so the absence is visible, and never upgrade without an external structure |

## What the composition checks cannot see

Worth stating plainly, because it is the reason three structural bugs survived
every check that existed at the time:

- **Mass and molecular formula are blind to connectivity.** HRMS matched 4/4 to
  ±0.0001 Da on a molecule with the acyl group on the wrong end of the steroid.
- Deleting a hydroxyl left the **oxygen count unchanged**; only hydrogen moved,
  by one.
- Moving an arm between nitrogens changes **nothing** about the formula.

C10 exists specifically to close this gap, and it is the only Tier-1 check that
compares against an external ground truth rather than internal consistency.

**2026-10-29 added a second instance of exactly this.** `AC_2025` held the
DLin-KC2-DMA dioxolane under the name `DLin-MC3-DMA` and the DLin-MC3-DMA ester
under the name `KC2`. Both are C43H79NO2 at **641.6111 Da** — the same formula
and the same exact mass to four decimals. Every composition check passed over
419 rows for as long as they were in the corpus.

It also taught a second lesson, the harder one. Finding that a name and a
structure disagree does **not** tell you which of the two is wrong, and the two
repairs are not equivalent: rewriting the structure keeps the value with the
name, rewriting the name keeps the value with the structure, and they assign
every measurement to opposite molecules. The first fix here rewrote the
structures and was wrong. Only the upstream source settled it — COMET's master
table stores SMILES and labels with **no names at all**, so the values belong to
the structures and the names were a downstream annotation (built, per its own
comment, by matching molecular weight, which cannot tell these two apart).

**So: when a name and a structure disagree, go to the source and find out which
one the measurement was recorded against. Do not repair the cheaper half.** An
internal-consistency check cannot answer this — it is satisfied by either
repair, including the wrong one.

## What a figure's own printed summaries are for

The `ql_2026` Fig 3c axis error (class 20) is the cleanest case yet of a bug
that **no internal check could see** and that the paper itself gave away three
times over. The values were self-consistent, the bars were detected correctly,
the z-scores were perfect — and every number was wrong by up to 1.76×.

What caught it: the panel's own printed hit rates, and a sentence in the
Results. What let it survive a first pass: `validate.py` asserted two of the
three printed rates and **excused the third** as "a boundary case" — and the
excused one was the only one that discriminated.

The rule this earns: **when a paper prints a statistic derived from a figure you
are digitizing, reproduce ALL of it, and never excuse the part that disagrees.**
A disagreeing summary statistic is the single most informative thing available,
because it is computed from the data you are trying to recover and is
independent of how you recover it. Excusing one is how a calibration error gets
through a passing test suite.

C10 could not have caught it either: C10 only looks at lipids a *curation folder*
enumerated, and these are catalogue reagents that arrive named. That gap is what
C11 closes — a registry of named commercial reagents whose structure is fixed and
public, checked by InChIKey skeleton wherever that name appears in the corpus.
A name like `DSPC` or `SM-102` means one molecule; the corpus should be able to
prove it means that molecule everywhere, and a shorthand like `18PG` should be
expanded once, in one place, with the paper that defines it cited next to it.

