#!/usr/bin/env python3
"""Build the only ionizable-lipid structure this paper discloses: H9T6.

H9  = tris(2-aminoethyl)amine (TREN)  -- SI p.5, named verbatim in the H9T6 prep.
T6  = (S)-glycidyl decanoate          -- derived three independent ways, see 07_qc.json.

Nothing here is hand-written as a product. Every bond is formed through
curation/lib/assemble, which asserts stereo preservation and refuses radicals.
The one operation assemble does not cover -- opening the epoxide -- is done by
removing a BOND (never an atom) at the carbon that is NOT the stereocentre, and
the stereocentre's neighbour ordering is asserted unchanged across it.
"""
import sys, os, csv, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "lib"))
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import assemble

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}

# ---------------------------------------------------------------- the tail
# The SI NAMES the reagent: "(S)-(-)-Glycidyl decanoate" (SI p.4). So the
# stereo assignment is derived FROM THE PUBLISHED NAME, by asking RDKit's CIP
# implementation which of the two candidates is S -- not by carrying a chiral
# tag through a bond-forming step. That also sidesteps a genuine false positive
# in assemble.esterify's stereo guard: esterifying (S)-glycidol RE-LABELS C2
# from S to R with no spatial change, because the ester breaks the C1/C3 CIP
# tie that glycidol's free OH leaves open. check_stereo_preserved compares CIP
# COUNTS, so it reads that relabel as an inversion and refuses. Reported in
# 08_contradictions.csv; curation/lib is NOT edited here.
DECANOIC = "CCCCCCCCCC(=O)O"
assert rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(DECANOIC)) == "C10H20O2"

# constitution: built by the library, on ACHIRAL glycidol, so nothing can be lost
FLAT_TAIL = assemble.esterify(DECANOIC, "OCC1CO1", label="glycidyl decanoate (flat)")
assert not assemble.stereocentres(FLAT_TAIL) or \
       set(assemble.stereocentres(FLAT_TAIL).values()) == {"?"}, FLAT_TAIL
out["T6_epoxide_flat"] = FLAT_TAIL

# stereo: enumerate both tags on that constitution and select by CIP code
_ring_c = [i for i, c in assemble.stereocentres(FLAT_TAIL).items()]
assert len(_ring_c) == 1, _ring_c
_cands = {}
for tag in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
    _m = Chem.RWMol(Chem.MolFromSmiles(FLAT_TAIL))
    _m.GetAtomWithIdx(_ring_c[0]).SetChiralTag(tag)
    _mm = _m.GetMol(); Chem.SanitizeMol(_mm)
    Chem.AssignStereochemistry(_mm, cleanIt=True, force=True)
    _code = [a.GetProp("_CIPCode") for a in _mm.GetAtoms() if a.HasProp("_CIPCode")]
    assert len(_code) == 1, _code
    _cands[_code[0]] = Chem.MolToSmiles(_mm)
assert set(_cands) == {"R", "S"}, _cands
TAIL_EPOXIDE = _cands["S"]           # "(S)-(-)-Glycidyl decanoate", SI p.4
out["T6_epoxide"] = TAIL_EPOXIDE
out["T2_epoxide_R_enantiomer"] = _cands["R"]
assert rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(TAIL_EPOXIDE)) == "C13H24O3"

# ------------------------------------------------------- open the epoxide
def open_epoxide(epoxide_smiles):
    """Break the O-CH2 bond of the oxirane. Removes NO atom.

    The stereocentre is the SUBSTITUTED ring carbon; we never touch a bond at
    it, and we assert its neighbour ordering and chiral tag are byte-identical
    before and after. That is the F-013 guard made mechanical: the bug there
    was RemoveAtom silently editing a stereocentre's neighbour list.
    """
    m = Chem.MolFromSmiles(epoxide_smiles)
    ring = Chem.MolFromSmarts("[CX4H2;r3][OX2;r3][CX4H1;r3]")
    hits = m.GetSubstructMatches(ring)
    if len(hits) != 1:
        raise SystemExit(f"open_epoxide: {len(hits)} oxirane matches, need 1")
    c_term, o, c_sub = hits[0]
    before = (m.GetAtomWithIdx(c_sub).GetChiralTag(),
              [n.GetIdx() for n in m.GetAtomWithIdx(c_sub).GetNeighbors()])
    rw = Chem.RWMol(m)
    rw.RemoveBond(c_term, o)
    for idx in (c_term, o):
        a = rw.GetAtomWithIdx(idx)
        a.SetNumExplicitHs(0); a.SetNoImplicit(False); a.SetNumRadicalElectrons(0)
    rw.GetAtomWithIdx(c_term).SetAtomMapNum(2)      # the attachment carbon
    prod = rw.GetMol()
    Chem.SanitizeMol(prod)
    after = (prod.GetAtomWithIdx(c_sub).GetChiralTag(),
             [n.GetIdx() for n in prod.GetAtomWithIdx(c_sub).GetNeighbors()])
    if before != after:
        raise SystemExit(f"open_epoxide: stereocentre disturbed {before} -> {after}")
    assemble.check_no_radicals(prod, label="opened arm")
    return Chem.MolToSmiles(prod), before == after

ARM_MAPPED, _ = open_epoxide(TAIL_EPOXIDE)
out["T6_arm_mapped"] = ARM_MAPPED
arm_plain = Chem.MolToSmiles(Chem.MolFromSmiles(ARM_MAPPED.replace(":2", "")))
assert rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(arm_plain)) == "C13H26O3", arm_plain

# --------------------------------------------------------------- the head
TREN = "NCCN(CCN)CCN"
assert rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(TREN)) == "C6H18N4"
out["H9_TREN"] = Chem.MolToSmiles(Chem.MolFromSmiles(TREN))

# ------------------------------------- add arms, one N-H at a time, n = 1..6
def add_arm(core_smiles):
    """Alkylate the FIRST remaining N-H on a primary or secondary amine."""
    m = Chem.MolFromSmiles(core_smiles)
    patt = Chem.MolFromSmarts("[NX3;H1,H2;!$(N-C=O)]")
    hits = [h[0] for h in m.GetSubstructMatches(patt)]
    if not hits:
        raise SystemExit("add_arm: no N-H left")
    rw = Chem.RWMol(m)
    rw.GetAtomWithIdx(hits[0]).SetAtomMapNum(1)
    return assemble.attach(Chem.MolToSmiles(rw.GetMol()), ARM_MAPPED,
                           a_map=1, b_map=2, label="N-alkylation")

series, cur = {}, TREN
for n in range(1, 7):
    cur = add_arm(cur)
    mol = Chem.MolFromSmiles(cur)
    series[n] = {
        "n_tails": n, "smiles": cur,
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "exact_mass": Descriptors.ExactMolWt(mol),
        "mw": Descriptors.MolWt(mol),
        "n_stereocentres": len(assemble.stereocentres(mol)),
    }
out["series"] = series

H9T6 = series[6]["smiles"]
assert series[6]["n_stereocentres"] == 6, series[6]
assemble.check_no_radicals(H9T6, label="H9T6")

# -------- the glycidol-capped 6-armed species (for the 1362.11 LC-MS peak)
glyc_arm, _ = open_epoxide("OCC1CO1")  # cap species: mass only, stereo irrelevant
capped = assemble.attach(
    Chem.MolToSmiles(
        (lambda m: (m.GetAtomWithIdx(
            [h[0] for h in m.GetSubstructMatches(
                Chem.MolFromSmarts("[NX3;H1,H2;!$(N-C=O)]"))][0]).SetAtomMapNum(1), m)[1])(
        Chem.RWMol(Chem.MolFromSmiles(series[5]["smiles"])))),
    glyc_arm, a_map=1, b_map=2, label="glycidol cap")
cm = Chem.MolFromSmiles(capped)
out["capped_5tail_plus_glycidol"] = {
    "smiles": capped, "formula": rdMolDescriptors.CalcMolFormula(cm),
    "exact_mass": Descriptors.ExactMolWt(cm)}

with open(os.path.join(HERE, "02_enumerate.json"), "w") as f:
    json.dump(out, f, indent=1)

# ------------------------------------------------------------- 01_blocks.csv
with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["block_role","block_id","name","SMILES","source","note"])
    w.writerow(["head","H9","tris(2-aminoethyl)amine (TREN)",out["H9_TREN"],
                "SI p.5, named verbatim in the H9T6 preparation",
                "the ONLY one of the library's 12 head amines disclosed anywhere on disk"])
    w.writerow(["tail","T6","(S)-glycidyl decanoate",TAIL_EPOXIDE,
                "SI p.4 epoxide-tail list + p.5 H9T6 NMR/LC-MS",
                "chain length C10 derived 3 ways (NMR integral, LC-MS ladder, R/S index pairing); see 07_qc.json"])
    for nm, smi, src in [
        ("T_(S)-glycidyl octanoate", "(S)", "SI p.4"),
        ("T_(S)-glycidyl dodecanoate", "(S)", "SI p.4"),
        ("T_(S)-glycidyl tetradecanoate", "(S)", "SI p.4"),
        ("T_(R)-glycidyl octanoate", "(R)", "SI p.4"),
        ("T_(R)-glycidyl decanoate", "(R)", "SI p.4"),
        ("T_(R)-glycidyl dodecanoate", "(R)", "SI p.4"),
        ("T_(R)-glycidyl tetradecanoate", "(R)", "SI p.4")]:
        w.writerow(["tail_unindexed","", nm, "", src,
                    "NMR-characterised in the SI but NEVER mapped to a T index; "
                    "the T1..T8 numbering is not printed anywhere on disk"])

# --------------------------------------------------------- 03_enumerated.csv
with open(os.path.join(HERE, "03_enumerated.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["IL_name","IL_SMILES","formula","exact_mass","MW","n_tails",
                "n_stereocentres","emitted_in_05_values","confidence",
                "structure_assumed","positional_ambiguity","curation_note"])
    s6 = series[6]
    w.writerow(["H9T6", H9T6, s6["formula"], f"{s6['exact_mass']:.4f}",
                f"{s6['mw']:.2f}", 6, 6, "yes", "medium", 0, 0,
                "Designed hexa-substituted product of tris(2-aminoethyl)amine with "
                "(S)-glycidyl decanoate; 1H NMR of the MPLC-purified material integrates "
                "to 6 tails (18H CH3, 12H Hbeta, 12H Halpha). HPLC-ELSD reports the "
                "4-tailed species as 55.26% of total area, so the isolated material is a "
                "mixture and this structure is its designed, NMR-supported major-by-NMR "
                "component, not a single pure compound."])
    for n in (4, 5):
        s = series[n]
        w.writerow([f"H9T6_{n}tail_LCMS_species", s["smiles"], s["formula"],
                    f"{s['exact_mass']:.4f}", f"{s['mw']:.2f}", n, n, "no", "low", 1, 1,
                    f"Partially substituted co-eluting species detected by LC-MS "
                    f"([M+H]+ computed {s['exact_mass']+1.00728:.2f}). Which of the three "
                    f"primary nitrogens carries how many arms is NOT determined by any "
                    f"evidence on disk -- the constitution is a positional isomer mixture. "
                    f"NOT emitted."])
print(json.dumps({k: v for k, v in out.items() if k != "series"}, indent=1))
for n, s in series.items():
    print(n, s["formula"], f"{s['exact_mass']:.4f}", f"[M+H]+ {s['exact_mass']+1.00728:.2f}")
print("capped", out["capped_5tail_plus_glycidol"]["formula"],
      f"[M+H]+ {out['capped_5tail_plus_glycidol']['exact_mass']+1.00728:.2f}")
print("H9T6 =", H9T6)
print("InChIKey", Chem.MolToInchiKey(Chem.MolFromSmiles(H9T6)))
