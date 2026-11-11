# Review packet — `xuq_32520527`

Generated from live files. Queue says **33** emitted rows; status `extracted`; proposed `Experiment_ID` **none**; PMID 32520527.

**NO FILE MATCHED THE DISCOVERY RULE** (a column in ('IL_SMILES', 'smiles', 'SMILES', 'canonical_SMILES', 'il_smiles', 'emitted_SMILES', 'product_SMILES', 'parent_SMILES') plus one starting `Experiment_value`).

**This is a statement about the discovery rule, NOT about the data.** An earlier version of this generator reported the same condition as an affirmative finding -- "this folder cannot be merged" -- and was wrong for 2 of 4 folders, writing off ~200 verified rows. Check the files below yourself before concluding anything.

Files that look row-like but lack one of those columns:
  - `01_blocks.csv` — 9 cols, IL_SMILES=False, value=[]
  - `03_enumerated.csv` — 15 cols, IL_SMILES=True, value=[]
  - `05_values.csv` — 23 cols, IL_SMILES=False, value=['Experiment_value_z_of_log10', 'Experiment_value_within_model']
  - `08_contradictions.csv` — 7 cols, IL_SMILES=False, value=[]
  - `08_retracted_suspicions.csv` — 4 cols, IL_SMILES=False, value=[]
  - `reference_structures.csv` — 7 cols, IL_SMILES=False, value=[]
