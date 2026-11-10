"""Stage 01 - building blocks for Zhao et al. 2020 (PMID 32662132).

Every tail here is built from the SI's OWN Experimental Procedures (tier 4,
PDF text layer), reagent by named reagent. Nothing is read off a figure.

The head is NOT drawn anywhere in the SI. It is taken from three independent
external records that agree (see 01_blocks_qc.json), not guessed.
"""
import csv, json, os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from rdkit import Chem, RDLogger
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))

# --------------------------------------------------------------- the heads
# Read from the corpus through the GUARDED lookup in curation/corpus_components.py,
# which RAISES on a name that carries more than one structure rather than picking
# by row order (error class 30, F-014/F-044). Both names resolve to exactly one
# structure each, and `93` / `A93` -- two spellings in two experiments -- resolve
# to the SAME structure, which is itself a cross-experiment agreement.
def _heads_from_corpus():
    import csv as _csv
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
    from corpus_components import load, smi
    rows = []
    for pool in ("vitro", "vivo"):
        pth = os.path.join(HERE, "..", "..", "data", "data_%s.csv" % pool)
        rows += list(_csv.DictReader(open(pth)))
    comp, amb = load(rows, prefixes=("IL_head",))
    return {n: smi(comp, amb, n) for n in ("93", "A93", "A9322")}


HEAD = {
    "block_id": "93",
    "block_type": "amine_head",
    "SMILES": "NCCCn1ccnc1",
    "chemical_name": "1-(3-aminopropyl)imidazole",
    "n_michael_slots": 2,
    "source": "NOT drawn in this paper's SI. Catalogue number 93 of the Xu/Anderson "
              "aliphatic-amine set, resolved from three independent external records "
              "that agree exactly (01_blocks_qc.json head_external_attestation).",
    "what_is_drawn": "nothing - the SI draws no amine head. The main-text figure that "
                     "draws them is behind the Wiley subscription wall.",
}

HEAD_9322 = {
    "block_id": "9322",
    "block_type": "amine_head",
    "SMILES": "Cc1nccn1CCCN",
    "chemical_name": "1-(3-aminopropyl)-2-methylimidazole",
    "n_michael_slots": 2,
    "source": "NOT drawn in this paper's SI. Resolved from the corpus through the "
              "guarded lookup: IL_head_name 'A9322' carries ONE structure across 22 "
              "YR_2024 rows. Independently CONFIRMED against this paper's own "
              "Figure S1 ESI-MS by axis geometry: the 9322-O17O base peak sits "
              "0.20 px from the prediction for this head and 3.78 px from the "
              "prediction for the plain-imidazole head 93 (09_figS1_mass_test.json).",
    "what_is_drawn": "nothing - the SI draws no amine head.",
}

# ------------------------------------------------- tails, from the SI text
# Each `route` is the SI section that builds it and the named reagents it uses.
TAILS = [
    dict(block_id="O17O", SMILES="C=CC(=O)OCCOCCCCCCCCCCCCCC",
         chemical_name="2-(tetradecyloxy)ethyl prop-2-enoate",
         si_section="1.1",
         route="ethylene glycol + NaH + 1-bromotetradecane/KI -> 2-(tetradecyloxy)ethanol; "
               "then acryloyl chloride/TEA -> acrylate",
         chain_after_linker="CCOCCCCCCCCCCCCCC", prefix="O", linker="ester"),
    dict(block_id="O17S", SMILES="C=CC(=O)OCCSCCCCCCCCCCCCCC",
         chemical_name="2-(tetradecylsulfanyl)ethyl prop-2-enoate",
         si_section="1.2",
         route="2-mercaptoethanol + 1-bromotetradecane/K2CO3 -> 2-(tetradecylthio)ethanol; "
               "then acryloyl chloride (as for O17O)",
         chain_after_linker="CCSCCCCCCCCCCCCCC", prefix="O", linker="ester"),
    dict(block_id="O17Se", SMILES="C=CC(=O)OCC[Se]CCCCCCCCCCCCCC",
         chemical_name="2-(tetradecylselanyl)ethyl prop-2-enoate",
         si_section="1.3",
         route="KSeCN + 2-bromoethanol -> selenocyanate; NaBH4 reduction; "
               "+ 1-bromotetradecane -> 2-(tetradecylselanyl)ethanol; then acryloyl chloride",
         chain_after_linker="CC[Se]CCCCCCCCCCCCCC", prefix="O", linker="ester"),
    dict(block_id="O16S-S", SMILES="C=CC(=O)OCCSSCCCCCCCCCCCC",
         chemical_name="2-(dodecyldisulfanyl)ethyl prop-2-enoate",
         si_section="1.4",
         route="1-dodecanol + 2,2'-dipyridyl disulfide -> dodecyl 2-pyridyl disulfide; "
               "+ 2-mercaptoethanol -> 2-(dodecyldisulfanyl)ethanol; then acryloyl chloride",
         chain_after_linker="CCSSCCCCCCCCCCCC", prefix="O", linker="ester"),
    dict(block_id="N16S-S", SMILES="C=CC(=O)NCCSSCCCCCCCCCCCC",
         chemical_name="N-[2-(dodecyldisulfanyl)ethyl]prop-2-enamide",
         si_section="1.4",
         route="as O16S-S but with cysteamine in place of 2-mercaptoethanol, so the "
               "acryloyl chloride acylates N and the linker is an amide",
         chain_after_linker="CCSSCCCCCCCCCCCC", prefix="N", linker="amide"),
]

# Named but NOT specified anywhere we can reach. These are NOT enumerated.
UNSPECIFIED = {
    "O18S-S": "named in the Figure S4 and S5 legends; the SI synthesises O16S-S and "
              "N16S-S only. The chain-atom naming rule (01_blocks_qc.json) predicts "
              "2-(tetradecyldisulfanyl)ethyl acrylate, but a prediction from a rule is "
              "not a structure read from a source. FLAGGED, NOT BUILT.",
    "O17C": "named in the Figure S1, S4 and S5 legends; synthesised nowhere in the SI. "
            "The rule fixes the chain at 17 heavy atoms but NOT what 'C' designates - "
            "a plain C17 alkyl chain, a carbonate, or a carbamate all fit the name. "
            "Linkers are formula-blind and the SI has no 13C spectrum and no legible "
            "HRMS, so THE TEST IS UNAVAILABLE. FLAGGED, NOT BUILT.",
    "EC18": "epoxide tail, attributed to reference [3] (Sun et al. Bioconjug. Chem. 2012). "
            "Not drawn, not synthesised here. FLAGGED, NOT BUILT.",
    "EC16": "epoxide tail, same provenance as EC18. FLAGGED, NOT BUILT.",
    "O16B": "appears ONLY as a Figure S1 panel label. Probably O16S-S under a second "
            "name ('B' for bisulfide); the chain count 16 matches. Unprovable here "
            "because Figure S1's m/z labels are 5-7 px tall and illegible. See "
            "08_contradictions.csv C-02. FLAGGED, NOT BUILT.",
    "N16B": "same as O16B, against N16S-S. FLAGGED, NOT BUILT.",
}

def chain_heavy_atoms(smi_fragment):
    """Heavy-atom count of the chain hanging off the ester O / amide N."""
    m = Chem.MolFromSmiles(smi_fragment)
    if m is None:
        raise ValueError("unparseable chain fragment: " + smi_fragment)
    return m.GetNumHeavyAtoms()

def main():
    rows, qc = [], {}
    corpus_heads = _heads_from_corpus()
    attest = {
        "93":    {"corpus_name": "93",    "corpus_SMILES": corpus_heads["93"]},
        "A93":   {"corpus_name": "A93",   "corpus_SMILES": corpus_heads["A93"]},
        "9322":  {"corpus_name": "A9322", "corpus_SMILES": corpus_heads["A9322"]},
    }
    for H in (HEAD, HEAD_9322):
        m = Chem.MolFromSmiles(H["SMILES"])
        if m is None:
            raise SystemExit("head SMILES does not parse: " + H["block_id"])
        key = H["block_id"]
        want = Chem.MolToSmiles(Chem.MolFromSmiles(attest[key]["corpus_SMILES"]))
        if Chem.MolToSmiles(m) != want:
            raise SystemExit("head %s disagrees with the corpus attestation %s"
                             % (key, want))
        rows.append(dict(block_id=H["block_id"], block_type="amine_head",
                         SMILES=Chem.MolToSmiles(m), formula=CalcMolFormula(m),
                         exact_mass=round(ExactMolWt(m), 4),
                         chemical_name=H["chemical_name"],
                         n_michael_slots=H["n_michael_slots"],
                         source=H["source"], what_is_drawn=H["what_is_drawn"]))

    name_rule = []
    for t in TAILS:
        mt = Chem.MolFromSmiles(t["SMILES"])
        if mt is None:
            raise SystemExit("tail does not parse: " + t["block_id"])
        n = chain_heavy_atoms(t["chain_after_linker"])
        # the number embedded in the block id, e.g. O16S-S -> 16
        digits = "".join(c for c in t["block_id"] if c.isdigit())
        name_rule.append(dict(tail=t["block_id"], name_token=int(digits),
                              measured_chain_heavy_atoms=n, agrees=(int(digits) == n),
                              prefix=t["prefix"], linker_from_route=t["linker"],
                              prefix_agrees=((t["prefix"] == "O") == (t["linker"] == "ester"))))
        rows.append(dict(block_id=t["block_id"], block_type="acrylate_tail",
                         SMILES=Chem.MolToSmiles(mt), formula=CalcMolFormula(mt),
                         exact_mass=round(ExactMolWt(mt), 4),
                         chemical_name=t["chemical_name"], n_michael_slots=0,
                         source="SI Experimental Procedures section %s (PDF text layer): %s"
                                % (t["si_section"], t["route"]),
                         what_is_drawn="nothing is drawn; the structure is read from the "
                                       "named reagents of its own synthesis"))

    with open(os.path.join(HERE, "01_blocks.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    qc["head_external_attestation"] = {
        "_what": "the heads are not drawn in this paper. Each is resolved from the "
                 "corpus through corpus_components.smi, which raises on an ambiguous "
                 "name instead of picking by row order.",
        "93": {"corpus_IL_head_name": ["93 (AA_2008, 11 rows)", "A93 (KW_2014, 5 rows)"],
               "corpus_SMILES": corpus_heads["93"],
               "two_spellings_one_structure": corpus_heads["93"] == corpus_heads["A93"],
               "third_record": "curation/xuq_32520527 read head 93 off a drawn figure "
                               "in a third Xu-lab paper and got the same structure"},
        "9322": {"corpus_IL_head_name": ["A9322 (YR_2024, 22 rows)"],
                 "corpus_SMILES": corpus_heads["A9322"],
                 "is_the_2_methyl_analogue_of_93": True,
                 "mass_confirmation": "09_figS1_mass_test.json"},
        "the_other_17_heads": {
            "names": ["9310", "9311", "9312", "9313", "9315", "9316", "9321", "9323",
                      "9324", "9331", "9332", "9333", "9334", "9341", "9351", "9352"],
            "corpus_attestation": "NONE. Zero rows in either pool carry any of these "
                                  "tokens in IL_name or IL_head_name.",
            "verdict": "NOT ENUMERABLE. 9322 = the 2-methyl analogue of 93 is ONE data "
                       "point; it does not define a digit grammar. The paper states no "
                       "grammar, and nothing here maps 9313 or 9331 to a substitution "
                       "pattern. THE GRAMMAR IS UNVERIFIED and no structure is invented."},
    }
    qc["n_blocks"] = len(rows)
    qc["n_heads"] = 2
    qc["n_tails_built"] = len(TAILS)
    qc["n_tails_named_but_not_built"] = len(UNSPECIFIED)
    qc["tails_not_built"] = UNSPECIFIED
    qc["name_token_rule"] = {
        "claim": "for every tail the SI actually synthesises, the number in its name "
                 "equals the heavy-atom count of the chain hanging off the ester oxygen "
                 "or amide nitrogen, and the leading letter is O for an ester and N for "
                 "an amide. Zero free parameters.",
        "per_tail": name_rule,
        "n_agree": sum(1 for r in name_rule if r["agrees"]),
        "n_tested": len(name_rule),
        "prefix_n_agree": sum(1 for r in name_rule if r["prefix_agrees"]),
        "why_it_matters": "the run brief's trap 5 says a name token is not a chain length. "
                          "Here it is, and that is MEASURED on five tails rather than "
                          "assumed. The rule is what makes O18S-S's and O17C's chain "
                          "lengths predictable - and it is still not enough to build them, "
                          "because it does not fix what 'C' designates.",
    }
    with open(os.path.join(HERE, "01_blocks_qc.json"), "w") as f:
        json.dump(qc, f, indent=1)
    print("blocks:", len(rows), "| name-rule", qc["name_token_rule"]["n_agree"], "/",
          qc["name_token_rule"]["n_tested"], "| prefix",
          qc["name_token_rule"]["prefix_n_agree"], "/", len(name_rule))

if __name__ == "__main__":
    main()
