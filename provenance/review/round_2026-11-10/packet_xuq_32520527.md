# Review packet — `xuq_32520527`

Generated from live files. Queue says **33** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 32520527.

**NO SCHEMA-VALID ROW FILE FOUND** (needs `IL_SMILES` and a column starting `Experiment_value`). That is itself a finding: this folder cannot be merged as it stands.

Files that look row-like but lack one of those columns:
  - `01_blocks.csv` — 9 cols, IL_SMILES=False, value=[]
  - `03_enumerated.csv` — 15 cols, IL_SMILES=True, value=[]
  - `05_values.csv` — 23 cols, IL_SMILES=False, value=['Experiment_value_z_of_log10', 'Experiment_value_within_model']
  - `08_contradictions.csv` — 7 cols, IL_SMILES=False, value=[]
  - `08_retracted_suspicions.csv` — 4 cols, IL_SMILES=False, value=[]
  - `reference_structures.csv` — 7 cols, IL_SMILES=False, value=[]
