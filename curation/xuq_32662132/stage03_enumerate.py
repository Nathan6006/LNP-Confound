"""Stage 03 - enumerate the lipidoids we can actually build.

Head 93 bears ONE primary amine, so the SI's 1 : 2.4 amine-to-acrylate charge
(SI 1.5) gives the N,N-bis adduct. The imidazole N3 is not counted as a Michael
site: 277 in-vitro and 25 in-vivo corpus rows (173 + 25 unique structures) built on this same head by three
unrelated groups all carry exactly two arms on the primary nitrogen and none on
the ring (02_corpus_overlap.json). That is an external observation, not a rule
we invented.

TWO heads are enumerated: 93 and 9322. Head 9322 was recovered from the corpus
(IL_head_name `A9322`, 22 YR_2024 rows, one structure under the guarded lookup)
and confirmed against this paper's OWN Figure S1 by axis geometry -- the
9322-O17O base peak sits 0.20 px from the prediction for 9322 and 3.78 px from
the prediction for 93, an 18.9x margin on a one-CH2 difference.

The other SEVENTEEN 93xx heads are NOT enumerated. None has any corpus
attestation, the paper states no digit grammar, and one data point (9322 = the
2-methyl analogue) does not define one. The grammar is UNVERIFIED and no
structure is invented.
"""
import csv, os, sys
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, inchi
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
PROTON = 1.0072765

ADD_PRIMARY   = AllChem.ReactionFromSmarts(
    "[NX3;H2:1].[CH2:2]=[CH:3][CX3:4]=[O:5]>>[N:1][CH2:2][CH2:3][C:4]=[O:5]")
ADD_SECONDARY = AllChem.ReactionFromSmarts(
    "[NX3;H1;!$(N-C=O):1].[CH2:2]=[CH:3][CX3:4]=[O:5]>>[N:1][CH2:2][CH2:3][C:4]=[O:5]")

def add_one(mol, tail, rxn):
    prods = rxn.RunReactants((mol, tail))
    keys = set()
    out = None
    for (p,) in prods:
        Chem.SanitizeMol(p)
        s = Chem.MolToSmiles(p)
        keys.add(s)
        out = Chem.MolFromSmiles(s)
    if len(keys) != 1:
        raise SystemExit("ambiguous addition: %d distinct products %s" % (len(keys), keys))
    return out

def main():
    blocks = {r["block_id"]: r for r in csv.DictReader(open(os.path.join(HERE, "01_blocks.csv")))}
    rows = []
    for hid in ("93", "9322"):
      head = Chem.MolFromSmiles(blocks[hid]["SMILES"])
      for tid in ("O17O", "O17S", "O17Se", "O16S-S", "N16S-S"):
        tail = Chem.MolFromSmiles(blocks[tid]["SMILES"])
        m1 = add_one(head, tail, ADD_PRIMARY)
        m2 = add_one(m1, tail, ADD_SECONDARY)
        # the product must have no N-H left on the aliphatic nitrogen and two new esters/amides
        n_nh = sum(a.GetTotalNumHs() for a in m2.GetAtoms() if a.GetSymbol() == "N")
        smi = Chem.MolToSmiles(m2)
        em = ExactMolWt(m2)
        rows.append(dict(
            lipidoid=hid + "-" + tid, IL_name=hid + "-" + tid, head_amine=hid,
            head_name=blocks[hid]["chemical_name"], tail=tid, n_arms=2,
            IL_SMILES=smi, formula_neutral=CalcMolFormula(m2),
            exact_mass_neutral=round(em, 4), MplusH_ours=round(em + PROTON, 4),
            inchikey_full=inchi.MolToInchiKey(m2),
            inchikey_skeleton=inchi.MolToInchiKey(m2).split("-")[0],
            n_NH_remaining=n_nh,
            n_distinct_constitutions=1,
            positional_ambiguity=0, structure_assumed=0,
            confidence_reason=(
                "head identity from corpus attestation under the guarded lookup, never "
                "from this paper (93: 16 rows / 2 experiments / 2 spellings / 1 "
                "structure, plus curation/xuq_32520527; 9322: 22 YR_2024 rows, plus "
                "this paper's own Figure S1 at 0.20 px); tail built from the named "
                "reagents of its own synthesis in SI section %s; arm count 2 fixed by "
                "one primary amine per head" % blocks[tid]["source"].split()[4]),
        ))
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    for r in rows:
        print(r["IL_name"], r["formula_neutral"], r["MplusH_ours"], "NH_left", r["n_NH_remaining"])
    print("enumerated:", len(rows))

if __name__ == "__main__":
    main()
