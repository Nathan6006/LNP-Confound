"""Stage 01/03 -- build the eight ionizable lipids of Tilstra/Couture-Senecal
et al., JACS 2023 (PMID 36652629) from named reagents, and the three benchmark
lipids from the project reagent registry.

HOW THE STRUCTURES ARE PINNED, in descending order of strength:

 1. The SI's own synthesis prose names every reagent and every step
    (pages 9-17). Nothing here is inferred from a name token.
 2. The SI prints a MOLECULAR FORMULA and an m/z for each final lipid. The
    formula pins the number of tails: four, i.e. both terminal primary amines
    exhaustively dialkylated. 07_qc.json checks ours against theirs.
 3. The paper DRAWS each lipid (Figs S1a, S2a, S3a, S4a, S7b, and a vector
    line drawing above each synthesis step). That is the C10 anchor material.

MASS CANNOT SEE CONNECTIVITY. The formula would be identical if a tail sat on
the amide nitrogen instead of the terminal amine, or if three tails sat on one
nitrogen and one on the other. What rules those out is the chemistry the SI
describes -- a large excess of epoxide or alkyl halide onto a primary amine,
which exhaustively dialkylates it -- together with the drawn structures.

ASSEMBLY. Every bond is made by curation/lib assemble.attach between two
ATOM-MAPPED attachment atoms. No string surgery on a canonical SMILES, no
RemoveAtom, no "atom index 0". Every fragment below is a real, valence-complete
molecule, so no fragment carries a radical into a product.
"""
import os, sys, csv, json

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lib"))
import assemble, reagents                                    # noqa: E402
from rdkit import Chem                                        # noqa: E402
from rdkit.Chem import rdMolDescriptors, Descriptors          # noqa: E402
from rdkit import RDLogger                                    # noqa: E402
RDLogger.DisableLog("rdApp.*")


# ---------------------------------------------------------------- remap -----
def remap(smiles, smarts, mapno, label=""):
    """Put atom map `mapno` on the ONE atom matching `smarts`. Raise otherwise.

    This is the only glue this folder adds to assemble.attach, and it exists
    because attach clears the map on an atom once it is bonded: a nitrogen that
    must be alkylated twice has to be re-tagged between the two calls. It is a
    MAP operation on a molecule object -- it never edits SMILES text, and it
    refuses when the target is not unique, so it cannot quietly tag the wrong
    nitrogen.
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError(f"remap({label}): unparseable {smiles!r}")
    hits = m.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    if len(hits) != 1:
        raise ValueError(
            f"remap({label}): SMARTS {smarts!r} matched {len(hits)} atoms in "
            f"{smiles!r}, need exactly 1. Refusing to pick one.")
    m.GetAtomWithIdx(hits[0][0]).SetAtomMapNum(mapno)
    return Chem.MolToSmiles(m)


# SMARTS for a secondary ALIPHATIC amine -- explicitly not an amide N, because
# after the first arm is attached the molecule contains both.
SEC_AMINE = "[NX3;H1;!$(N[CX3]=[OX1])]"

# ------------------------------------------------------------- fragments ----
# Each is a real molecule, named as the SI names it. The atom map marks the
# atom that forms the next bond; on attachment RDKit recomputes its hydrogens.

CORES = {   # the primary amine that takes the double aza-Michael (SI p9, p13, p15)
    "C3":  ("CCC[NH2:1]",   "n-propylamine"),
    "OC2": ("OCC[NH2:1]",   "2-amino-1-ethanol (used TBS-protected, then deprotected)"),
    "OC3": ("OCCC[NH2:1]",  "3-amino-1-propanol (used TBS-protected, then deprotected)"),
}

LINKERS = {  # the diamine that opens the methyl ester to the amide (SI p10, p12)
    "K2": ("[CH3:5]C[NH2:6]",  "ethylenediamine", 2),
    "K3": ("[CH3:5]CC[NH2:6]", "1,3-diaminopropane", 3),
}

TAILS = {   # four of these go on, two per terminal nitrogen
    "E10": ("[CH3:9]C(O)CCCCCCCC",            "1,2-epoxydecane",        10, "epoxide ring-opening"),
    "E14": ("[CH3:9]C(O)CCCCCCCCCCCC",        "1,2-epoxytetradecane",   14, "epoxide ring-opening"),
    "E15": ("[CH3:9]C(O)CCCCCCCCCCCCC",       "1,2-epoxypentadecane",   15, "epoxide ring-opening"),
    "A10": ("[CH3:9]CCCCCCCCC",               "1-bromodecane",          10, "N-alkylation with DIPEA"),
    "A14": ("[CH3:9]CCCCCCCCCCCCC",           "tetradecanal",           14, "reductive amination (STAB/AcOH)"),
}

# The propanamide arm contributed by methyl acrylate. Map 1 is the carbon that
# bonds to the core nitrogen; map 2 is the amide nitrogen, which bonds to the
# linker's first carbon. 3-aminopropanamide, valence-complete as written.
ARM = "[CH3:1]CC(=O)[NH2:2]"

# The eight lipids the paper actually made. 3 cores x 2 linkers x 5 tails = 30
# possible; 8 exist. This is an ITERATIVE design, not a combinatorial library.
LIPIDS = [
    ("C3-K2-E14",  "C3",  "K2", "E14", "Fig S1a; SI p10"),
    ("C3-K2-A14",  "C3",  "K2", "A14", "Fig S1a; SI p11"),
    ("C3-K2-E15",  "C3",  "K2", "E15", "SI p12 (structure drawn p11 bottom)"),
    ("C3-K2-E10",  "C3",  "K2", "E10", "Fig S2a; SI p12"),
    ("C3-K3-E10",  "C3",  "K3", "E10", "Fig S3a; SI p13"),
    ("OC2-K3-E10", "OC2", "K3", "E10", "Fig S4a; SI p14-15"),
    ("OC3-K3-E10", "OC3", "K3", "E10", "Fig S4a; SI p16"),
    ("OC2-K2-A10", "OC2", "K2", "A10", "Fig S7b; SI p16-17"),
]

# Printed by the SI for each final lipid: formula, the ion, and the m/z given
# in the "calc'd" and "found" columns. Transcribed by hand from the text layer.
SI_MS = {
    "C3-K2-E14":  dict(formula="C69H141N5O6",  ion="[M+H]+", calc=1137.92, found=1137.09, method="ESI+",       page=10),
    "C3-K2-A14":  dict(formula="C69H141N5O2",  ion="[M+H]+", calc=1073.11, found=1073.11, method="MALDI-TOF+", page=11),
    "C3-K2-E15":  dict(formula="C73H149N5O6",  ion="[M+H]+", calc=1193.15, found=1193.2,  method="MALDI-TOF+", page=12),
    "C3-K2-E10":  dict(formula="C53H109N5O6",  ion="[M]+",   calc=912.5,   found=912.8,   method="MALDI-TOF+", page=12),
    "C3-K3-E10":  dict(formula="C55H113N5O6",  ion="[M+H]+", calc=940.87,  found=940.87,  method="ESI+",       page=13),
    "OC2-K3-E10": dict(formula="C54H111N5O7",  ion="[M+H]+", calc=942.85,  found=942.85,  method="ESI+",       page=15),
    "OC3-K3-E10": dict(formula="C55H113N5O7",  ion="[M+H]+", calc=956.86,  found=956.88,  method="DART+",      page=16),
    "OC2-K2-A10": dict(formula="C52H107N5O3",  ion="[M+H]+", calc=850.84,  found=850.84,  method="ESI+",       page=17),
}

# Benchmarks. Taken from the project registry by their FULL registered names --
# reagents.smi refuses an ambiguous name, which is the class-30 guard. The
# paper's Fig S7a axis says "MC3"; Table S2 rows 25-27 say "DLin-MC3-DMA", so
# the expansion is the paper's own, not ours.
BENCHMARKS = [
    ("SM-102",       "SM-102",       "Table S2 IDs 1, 9, 13, 17, 19, 20, 21; drawn in Fig S1a"),
    ("DLin-MC3-DMA", "DLin-MC3-DMA", "Table S2 IDs 25-27; axis label 'MC3' in Fig S7a"),
    ("ALC-0315",     "ALC-0315",     "Fig S11 pKa table only (6.27). NOT in Table S2 and in no flux panel."),
]


def build_tailed_linker(lk, tl):
    """linker diamine, one nitrogen exhaustively dialkylated with two tails."""
    frag, _, _ = LINKERS[lk]
    tail = TAILS[tl][0]
    s = assemble.attach(frag, tail, a_map=6, b_map=9, label=f"{lk}+{tl} #1")
    s = remap(s, SEC_AMINE, 6, label=f"{lk}+{tl} re-tag N")
    s = assemble.attach(s, tail, a_map=6, b_map=9, label=f"{lk}+{tl} #2")
    return s


def build(core, lk, tl):
    tailed = build_tailed_linker(lk, tl)
    # map 5 (the linker's open carbon) survived both tail attachments.
    # map 1 on the arm (the carbon that bonds to the core N) survives this.
    arm = assemble.attach(ARM, tailed, a_map=2, b_map=5, label=f"arm-{lk}-{tl}")
    s = assemble.attach(CORES[core][0], arm, a_map=1, b_map=1, label=f"{core} arm1")
    s = remap(s, SEC_AMINE, 1, label=f"{core} re-tag core N")
    s = assemble.attach(s, arm, a_map=1, b_map=1, label=f"{core} arm2")
    return s


def props(smi):
    m = Chem.MolFromSmiles(smi)
    return dict(
        formula=rdMolDescriptors.CalcMolFormula(m),
        exact_mass=round(Descriptors.ExactMolWt(m), 4),
        average_mw=round(Descriptors.MolWt(m), 2),
        n_unspecified_stereocentres=len(
            Chem.FindMolChiralCenters(m, includeUnassigned=True,
                                      useLegacyImplementation=False)),
    )


def main():
    blocks, rows = [], []
    for k, (smi, name) in CORES.items():
        blocks.append(dict(block_type="core_amine", block_id=k, reagent=name,
                           fragment_SMILES=smi,
                           source="SI Supplementary Methods, pages 9/13/15"))
    for k, (smi, name, n) in LINKERS.items():
        blocks.append(dict(block_type="linker_diamine", block_id=k, reagent=name,
                           fragment_SMILES=smi,
                           source=f"SI pages 10/12; {n}-carbon bridge, token K{n}"))
    for k, (smi, name, n, how) in TAILS.items():
        blocks.append(dict(block_type="tail", block_id=k, reagent=name,
                           fragment_SMILES=smi,
                           source=f"SI; {n} carbons; {how}"))

    drawn = {a: e for a, _, _, _, e in LIPIDS}
    for nm, core, lk, tl, _where in LIPIDS:
        smi = build(core, lk, tl)
        p = props(smi)
        rows.append(dict(
            IL_name=nm, IL_SMILES=smi, core=core, linker=lk, tail=tl,
            n_tails=4, tail_carbons=TAILS[tl][2],
            formula_ours=p["formula"], formula_SI=SI_MS[nm]["formula"],
            exact_mass_ours=p["exact_mass"], average_mw_ours=p["average_mw"],
            SI_ion=SI_MS[nm]["ion"], SI_calc_mz=SI_MS[nm]["calc"],
            SI_found_mz=SI_MS[nm]["found"], SI_method=SI_MS[nm]["method"],
            SI_page=SI_MS[nm]["page"],
            n_unspecified_stereocentres=p["n_unspecified_stereocentres"],
            origin="enumerated from SI synthesis",
            drawn_at=drawn[nm],
        ))

    for disp, regname, where in BENCHMARKS:
        smi = reagents.smi(regname)
        p = props(smi)
        rows.append(dict(
            IL_name=disp, IL_SMILES=smi, core="", linker="", tail="",
            n_tails="", tail_carbons="",
            formula_ours=p["formula"], formula_SI="",
            exact_mass_ours=p["exact_mass"], average_mw_ours=p["average_mw"],
            SI_ion="", SI_calc_mz="", SI_found_mz="", SI_method="", SI_page="",
            n_unspecified_stereocentres=p["n_unspecified_stereocentres"],
            origin="provenance/reagent_registry (C11), not synthesised here",
            drawn_at=where,
        ))

    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(blocks[0]))
        w.writeheader(); w.writerows(blocks)
    with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

    bad = [r for r in rows if r["formula_SI"] and
           r["formula_ours"] != r["formula_SI"]]
    print(f"blocks        : {len(blocks)}")
    print(f"lipids built  : {len(LIPIDS)} novel + {len(BENCHMARKS)} benchmark "
          f"= {len(rows)}")
    print(f"formula ours == SI : {len(LIPIDS) - len(bad)}/{len(LIPIDS)}")
    for r in rows:
        print(f"  {r['IL_name']:12s} {r['formula_ours']:16s} "
              f"SI {r['formula_SI'] or '-':16s} "
              f"exact {r['exact_mass_ours']:9.4f} avg {r['average_mw_ours']:8.2f} "
              f"stereo {r['n_unspecified_stereocentres']}")
    if bad:
        print("\nFORMULA MISMATCH:", [b["IL_name"] for b in bad])
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
