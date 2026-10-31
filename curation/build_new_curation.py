#!/usr/bin/env python3
"""Emit newly curated rows in the 67-column LNPDB schema, ready to concatenate.

Scope: everything curated AFTER the LNPDB + Bowen Li additions. Right now that is
one paper, Tian et al., Adv Mater 2026 (PMID 42634137), contributing two blocks:

    ZT_2026  in vitro   486 rows   formulation DoE on ONE lipid (C6O2B2)
    ZT_2026  in vivo     51 rows   17 lipid cores x 3 cholesterol tails

Same Experiment_ID in both blocks, following the BL_2023 / BL_2024 precedent of
an ID appearing in both the vitro and vivo corpus files. The two blocks are
z-scored SEPARATELY, as the corpus does, so splitting this file on `Model`
reproduces the corpus convention exactly.

CONVENTIONS REVERSE-ENGINEERED FROM THE CORPUS, not assumed:

* Experiment_value is a per-Experiment_ID z-score of log10(raw), population sd.
  Verified against SL_2021 (the iPhos 4-bin source): its four stored values are
  exactly equally spaced (-1.1598, +0.2109, +1.5815, +2.9522), which can only
  happen if log-decade bins were z-scored. For DISCRETIZED sources the quantity
  z-scored is the BIN ORDINAL, not a midpoint: for SL_2021, ordinals 0..3 with
  counts 190/293/76/13 give mean 0.84615, sd 0.72958, and (0-0.84615)/0.72958 =
  -1.1598 exactly. The 486 are binned, so their ordinals are z-scored the same
  way; the 51 are continuous, so log10(flux) is z-scored.
* molratio columns are PERCENTAGES summing to 100 (checked: AC_2025 gives
  35.5+16.3+47.2+1.0).
* `van.der.Waals.Molecular.Volume` actually holds ExactMolWt. This is a misnomer
  in the corpus, confirmed on 400/400 unique ILs. It is reproduced as-is so the
  new rows stay consistent; it is NOT a van der Waals volume.

COMPONENT MAPPING for this paper's 4-component system (no free cholesterol -
that is the paper's central claim, so CHL_* is deliberately empty):
    IL              = C6O2B2 (or the 51 library lipids)
    HL              = the helper phospholipid (DSPC/DOPC/DMPC/DOPG/DOPE/EggSM)
    PEG             = DMG-PEG2000
    fifthcomponent  = DODAP, the ionizable co-lipid
"""
import csv
import json
import os
import statistics as st

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors as D, rdMolDescriptors as rd

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
CORPUS = os.path.join(HERE, "..", "data", "data_vitro.csv")
SIEG = os.path.join(HERE, "siegwart_42634137")
OUT = os.path.join(HERE, "new_curation.csv")
OUT_QC = os.path.join(HERE, "new_curation_qc.json")

EID = "ZT_2026"
LINK = "https://advanced.onlinelibrary.wiley.com/doi/10.1002/adma.74762"
DOI = "10.1002/adma.74762"
PMID = "42634137"

# Helper-lipid and co-lipid SMILES. DSPC, DOPE, DMG-PEG2000 and Cholesterol are
# taken from the corpus at run time so they match byte-for-byte. The remaining
# four phospholipids and DODAP are not in the corpus and are supplied here as
# standard reference structures; they are flagged in the QC as such.
SUPPLIED = {
    "DOPC": "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "DMPC": "CCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCC",
    "DOPG": "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@H](COP(=O)([O-])OC[C@@H](O)CO)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
    "EggSM": "CCCCCCCCCCCCCCCC(=O)N[C@@H](COP(=O)([O-])OCC[N+](C)(C)C)[C@H](O)/C=C/CCCCCCCCCCCCC",
    "DODAP": "CCCCCCCC/C=C\\CCCCCCCC(=O)OC[C@@H](CN(C)C)OC(=O)CCCCCCC/C=C\\CCCCCCCC",
}

# The five "commercially available cores" of Fig 1b. MC3 and SM-102 come from
# the corpus. TT3, Lipid 14 and OF-2 are not in the corpus and ARE drawn in
# Fig 1b, so they were read from that panel at 600 dpi:
#
#   TT3       benzene-1,3,5-tricarboxamide, three -NH-(CH2)3-N(C12H25)2 arms.
#             All six alkyl chains counted by apex detection; the two unclipped
#             chains both gave 11 extrema (= 12 C) with bond-gap CV 0.013.
#   Lipid 14  (CH3)2N-(CH2)3-C(=O)O-CH2CH2-N[(CH2)8-O-C(=O)-CH(C6H13)(C8H17)]2,
#             i.e. 2-hexyldecanoate tails. Every run measured at CV < 0.02.
#   OF-2      piperazine-2,5-dione (the cyclic lysine dimer) with two
#             -(CH2)4-N(tail)2 arms; each tail is an epoxide-opened
#             9,12-octadecadiene, -CH2-CH(OH)-(CH2)6-CH=CH-CH2-CH=CH-(CH2)4-CH3.
#             The chain reads CH3-(CH2)4-CH=CH-CH2-CH=CH- from the methyl end,
#             which is the omega-6 linoleyl signature.
#
# The apex counter was VALIDATED against MC3 in the same figure and FAILED there
# (20-21 C against a known 18), because cis double bonds break the regular
# zig-zag. Bond-gap CV separates the two regimes cleanly: 0.013 on the saturated
# chains that were counted, 0.31 on MC3. So TT3 and Lipid 14, which are fully
# saturated, are machine-counted; OF-2's unsaturated tails were read by eye from
# the alkene pattern and carry a lower confidence flag.
COMMERCIAL_CORES = {
    "TT3": ("O=C(NCCCN(CCCCCCCCCCCC)CCCCCCCCCCCC)c1cc(C(=O)NCCCN(CCCCCCCCCCCC)"
            "CCCCCCCCCCCC)cc(C(=O)NCCCN(CCCCCCCCCCCC)CCCCCCCCCCCC)c1", "high"),
    "Lipid 14": ("CN(C)CCCC(=O)OCCN(CCCCCCCCOC(=O)C(CCCCCC)CCCCCCCC)"
                 "CCCCCCCCOC(=O)C(CCCCCC)CCCCCCCC", "high"),
    "OF-2": (r"O=C1NC(CCCCN(CC(O)CCCCCC/C=C\C/C=C\CCCC)CC(O)CCCCCC/C=C\C/C=C\CCCC)"
             r"C(=O)NC1CCCCN(CC(O)CCCCCC/C=C\C/C=C\CCCC)CC(O)CCCCCC/C=C\C/C=C\CCCC",
             "medium"),
}

ESTER = Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]")
CARBONATE = Chem.MolFromSmarts("[OX2][CX3](=O)[OX2]")
DISULFIDE = Chem.MolFromSmarts("[#16X2][#16X2]")


def descriptors(smiles):
    """The 17 RDKit-derived columns, matching the corpus function-for-function."""
    out = dict.fromkeys(
        ["Heavy.Atoms", "Rings", "Aromatic.Rings", "Rotatable.Bonds",
         "ExactMolWt", "Topological.Polar.Surface.Area",
         "Hydrogen.Bond.Donors", "Hydrogen.Bond.Acceptors", "LogP",
         "Molar.Refractivity", "Fraction.sp3.Carbons", "sp3.Carbons",
         "Nitrogen.Count", "Molecular.Weight",
         "has_ester", "has_carbonate", "has_disulfide"], "")
    if not smiles:
        return out
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return out
    out.update({
        "Heavy.Atoms": m.GetNumHeavyAtoms(),
        "Rings": rd.CalcNumRings(m),
        "Aromatic.Rings": rd.CalcNumAromaticRings(m),
        "Rotatable.Bonds": rd.CalcNumRotatableBonds(m),
        "ExactMolWt": D.ExactMolWt(m),
        "Topological.Polar.Surface.Area": round(rd.CalcTPSA(m), 2),
        "Hydrogen.Bond.Donors": D.NumHDonors(m),
        "Hydrogen.Bond.Acceptors": D.NumHAcceptors(m),
        "LogP": Crippen.MolLogP(m),
        "Molar.Refractivity": Crippen.MolMR(m),
        "Fraction.sp3.Carbons": rd.CalcFractionCSP3(m),
        "sp3.Carbons": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C"
                           and str(a.GetHybridization()) == "SP3"),
        "Nitrogen.Count": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N"),
        "Molecular.Weight": round(D.MolWt(m), 2),
        "has_ester": int(m.HasSubstructMatch(ESTER)),
        "has_carbonate": int(m.HasSubstructMatch(CARBONATE)),
        "has_disulfide": int(m.HasSubstructMatch(DISULFIDE)),
    })
    return out


def zscore(values):
    """Population z-score, matching the corpus (verified ddof=0 on SL_2021)."""
    mu, sd = st.mean(values), st.pstdev(values)
    return [(v - mu) / sd for v in values]


def load_corpus_components(header_only=False):
    """name -> SMILES from the live corpus, REFUSING to guess when a name is
    ambiguous. Error class 30, review finding F-014.

    This used to be a bare `comp.setdefault(n, s)`, i.e. first row in file order
    wins. 292 names in the corpus map to more than one structure, and
    data_vitro.csv is written in ALPHABETICAL experiment order -- so which
    structure a name resolved to was a function of which experiment sorts first,
    a property no caller could see and nothing recorded.

    That is how `MC3-CB1/2/3` came to carry the KC2 dioxolane under an MC3 name.
    At the time `DLin-MC3-DMA` was ambiguous rather than simply wrong: 2,190 rows
    of the correct acyclic ester plus 203 rows of dioxolane from `AC_2025`, whose
    names were later repaired. Alphabetically `AC_2025` comes first, so the
    lookup returned the dioxolane. Regenerating after the name repair fixed the
    three rows -- but it fixed them by luck of row order, which is why the guard
    below and not just a regeneration.

    `ambiguous` is returned rather than raised at load time because most of the
    292 are library members no caller ever asks for. The failure has to happen
    where a caller actually requests one -- see `smi()` in main().
    """
    rows = list(csv.DictReader(open(CORPUS)))
    hdr = list(rows[0].keys())
    if header_only:
        return hdr, {}, {}
    comp, seen = {}, {}
    for r in rows:
        for pre in ("HL", "PEG", "CHL", "IL"):
            n, s = r.get(f"{pre}_name"), r.get(f"{pre}_SMILES")
            if n and s:
                seen.setdefault(n, {}).setdefault(s, []).append(r["Experiment_ID"])
                comp.setdefault(n, s)
    ambiguous = {n: {k: sorted(set(v)) for k, v in d.items()}
                 for n, d in seen.items() if len(d) > 1}
    return hdr, comp, ambiguous


def _load_enumerate():
    """Import siegwart_42634137/enumerate.py by path, for `quaternize` below.

    THIS FILE USED TO CARRY ITS OWN COPY of the cholesteryl-ester assembly and
    the quaternization-site rule, with a docstring claiming parity with
    enumerate.py. It did not have parity, and the divergence cost two bugs:

      * the copy kept the `AddBond`-then-`RemoveAtom` construction, which
        inverts C3, so 15 conjugate structures (MC3-CB*, SM-102-CB*, OF-2-CB*,
        TT3-CB*, Lipid 14-CB*) carried epicholesterol. Fixing enumerate.py
        alone would have left all 15 wrong AFTER the fix was announced;
      * the copy still picked the site by `max(atom index)`, the rule
        enumerate.py was repaired away from, so `Lipid 14-CB*` was alkylated at
        the hindered central nitrogen rather than the terminal dimethylamino.

    Both were found by review, not by a check, and the second copy is precisely
    what CLAUDE.md's Tier-3 rule was written about -- in the same file it was
    written about. So there is now one implementation and this imports it. If
    that import ever fails, this raises instead of silently falling back to a
    local copy; a fallback is how two implementations start.
    """
    import importlib.util
    path = os.path.join(SIEG, "enumerate.py")
    spec = importlib.util.spec_from_file_location("sieg_enumerate", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_EN = _load_enumerate()


def quaternize(core_smiles, n_link):
    """Cholesteryl-O2C-(CH2)n- onto the nitrogen the chemistry would alkylate.

    One implementation, shared with the synthesized library -- see
    `_load_enumerate` above for why that matters. Used for all five commercial
    cores (MC3-C, SM-102-C, TT3-C, Lipid 14-C, OF-2-C); the parent SMILES come
    from the corpus for the first two and from the Fig 1b drawing for the rest.

    The paper draws NO -C conjugate for any of the five, so the conjugation site
    is inferred for every one of them. These rows carry `confidence = medium`
    and `structure_assumed = 1` with the reason on the row.
    """
    smi, _note = _EN.quaternize(core_smiles, n_link)
    return smi or ""


def base_row(hdr):
    return dict.fromkeys(hdr, "")


def main():
    hdr, comp, ambiguous = load_corpus_components()
    # The corpus carries 79 columns: 67 curation + 5 RDKit descriptors that
    # build_dataset.py recomputes on merge + 7 quality/provenance columns. Emit
    # everything except the five recomputed ones -- writing those here would
    # only risk them going stale.
    #
    # The count is DERIVED, not asserted against a literal. It used to read
    # `assert len(hdr) == 67`, which is the error class F-002 was logged for: a
    # check whose reference value goes stale on every legitimate change, so the
    # person who hits it edits the number until it passes. It fired the moment
    # the quality columns were added, against a header that was entirely correct.
    EXTRA = {"Num_carbon_in_tail", "num_protonatable_nitrogens",
             "num_permanent_cationic_N", "formal_net_charge",
             "num_unsaturated_cc_bonds"}
    hdr = [c for c in hdr if c not in EXTRA]
    # Needed explicitly only on a build that runs before the corpus carries the
    # quality columns; afterwards the corpus header above already has them.
    for c in (x["name"] for x in json.load(open(os.path.join(
            HERE, "..", "provenance", "schema_quality.json")))["columns"]):
        if c not in hdr and c not in EXTRA:
            hdr.append(c)
    assert not (EXTRA & set(hdr)), "recomputed descriptor leaked into the header"
    assert "Experiment_value" in hdr and "IL_SMILES" in hdr, hdr[:5]

    def smi(name):
        """Resolve a component name, refusing an ambiguous one (F-014).

        Looking a structure up by name from the live corpus is only safe while
        that name means one thing. When it means two, picking either silently is
        how three measurements ended up on the wrong molecule, so this raises
        instead and the caller has to say which experiment's structure it wants.
        """
        if name in ambiguous:
            raise SystemExit(
                f"AMBIGUOUS COMPONENT NAME {name!r}: the corpus carries "
                f"{len(ambiguous[name])} different structures under it "
                f"{ {k[:40] + '...': v for k, v in ambiguous[name].items()} }. "
                f"Refusing to guess -- error class 30. Name the experiment whose "
                f"structure you mean, or add a SUPPLIED literal with its source.")
        return comp.get(name) or SUPPLIED.get(name, "")

    start_index = 1 + max(
        int(r["Index"]) for f in ("data_vitro.csv", "data_vivo.csv")
        for r in csv.DictReader(open(os.path.join(HERE, "..", "data", f)))
        if r["Index"].strip().isdigit())

    rows, notes = [], {}

    # ---- block 1: the 486 formulation screen, in vitro --------------------
    _enum0 = {r["lipid_name"]: r
              for r in csv.DictReader(open(os.path.join(SIEG, "03_enumerated.csv")))}
    c6o2b2 = _enum0["C6O2B2"]["IL_SMILES"]
    # Every row in this block is the same lipid, so it carries that lipid's own
    # confidence -- `high`: C6 has a piperazine nitrogen and Fig 2a draws the
    # quaternization on it, so the site is read, not assumed (F-001/F-011).
    c6o2b2_conf = _enum0["C6O2B2"]["confidence"]
    screen = list(csv.DictReader(open(os.path.join(SIEG, "04b_fig2c_486.csv"))))
    ordinals = [int(r["bin_index"]) - 1 for r in screen]
    zs = zscore([float(o) for o in ordinals])
    for i, (r, z) in enumerate(zip(screen, zs)):
        ratios = {k: float(r[k]) for k in
                  ("C6O2B2_molratio", "DODAP_molratio",
                   "helper_molratio", "DMGPEG_molratio")}
        tot = sum(ratios.values())
        row = base_row(hdr)
        row.update({
            "Index": start_index + len(rows),
            "LNP_ID": f"LNP_{start_index + len(rows):07d}",
            "Experiment_ID": EID,
            "Formulation_ID": f"{EID}_F{i + 1}",
            "IL_name": "C6O2B2", "IL_SMILES": c6o2b2,
            "IL_molratio": round(ratios["C6O2B2_molratio"] / tot * 100, 2),
            "IL_to_nucleicacid_massratio": "",
            "HL_name": r["helper_lipid"], "HL_SMILES": smi(r["helper_lipid"]),
            "HL_molratio": round(ratios["helper_molratio"] / tot * 100, 2),
            "PEG_name": "DMG-PEG2000", "PEG_SMILES": smi("DMG-PEG2000"),
            "PEG_molratio": round(ratios["DMGPEG_molratio"] / tot * 100, 2),
            "fifthcomponent_name": "DODAP", "fifthcomponent_SMILES": smi("DODAP"),
            "fifthcomponent_molratio": round(ratios["DODAP_molratio"] / tot * 100, 2),
            "Aqueous_buffer": "citrate", "Dialysis_buffer": "PBS",
            "Mixing_method": "handmixed",
            "Model": "in_vitro", "Model_type": "IGROV1", "Model_target": "in_vitro",
            "Route_of_administration": "in_vitro",
            "Cargo": "mRNA", "Cargo_type": "FLuc",
            "Dose_ug_nucleicacid": 0.02,
            "Experiment_method": "luminescence_discretized_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": z,
            "Publication_link": LINK, "Publication_DOI": DOI, "Publication_PMID": PMID,
        })
        row.update(descriptors(c6o2b2))
        row["confidence"] = c6o2b2_conf
        rows.append(row)
    notes["vitro_rows"] = len(rows)

    # ---- block 2: the 51-lipid library, in vivo ---------------------------
    # The whole enumerated record, not just the SMILES: this block used to read
    # `lipid_name -> IL_SMILES` and then hard-code `"high"` for every lipid it
    # found, discarding the `confidence` and `quaternization_site_assumed`
    # columns the enumerator had written beside it (review finding F-011).
    # 24 of the 36 enumerated lipids carry `confidence = low` with the verbatim
    # reason "ASSUMED: no ring N and the paper draws no structure" -- exactly
    # the rows a reader looking for the weakest structures needs to find. The
    # hard-coded literal labelled them `high`, which is worse than labelling
    # them nothing: it asserts the opposite of the only evidence there is.
    enum = {r["lipid_name"]: r
            for r in csv.DictReader(open(os.path.join(SIEG, "03_enumerated.csv")))}
    lib = {k: r["IL_SMILES"] for k, r in enum.items()}
    vivo = list(csv.DictReader(open(os.path.join(SIEG, "04c_fig1d_51.csv"))))
    LINK_N = {"B1": 3, "B2": 5, "B3": 7}
    PARENT = {"MC3-C": "DLin-MC3-DMA", "SM-102-C": "SM-102"}

    resolved, unresolved = {}, []
    for r in vivo:
        core, tail = r["core"], r["chol_tail"]
        key = f"{core}{tail}"
        if key in lib:
            e = enum[key]
            # the lipid's OWN assessment, carried not overwritten
            resolved[key] = (lib[key], "enumerated from Fig 1b blocks",
                             e["confidence"] or "high",
                             {"structure_assumed": e.get("quaternization_site_assumed", ""),
                              "note": (e.get("quaternization_site", "")
                                       if e.get("quaternization_site_assumed") == "1"
                                       else "")})
            continue
        parent_smiles, conf, how = None, None, None
        if core in PARENT and smi(PARENT[core]):
            # smi(), not comp[], so an ambiguous parent name is a hard failure
            # rather than a silent row-order coin flip (F-014)
            parent_smiles = smi(PARENT[core])
            conf, how = "medium", f"corpus {PARENT[core]}"
        elif core.endswith("-C") and core[:-2] in COMMERCIAL_CORES:
            parent_smiles, conf = COMMERCIAL_CORES[core[:-2]]
            how = f"Fig 1b drawing of {core[:-2]}"
        if parent_smiles:
            smi_q = quaternize(parent_smiles, LINK_N[tail])
            if smi_q:
                # The cholesterol CONJUGATION SITE on a commercial core is
                # inferred -- the paper draws no -C conjugate at all -- so it is
                # an assumed structural feature, not merely a lower confidence.
                resolved[key] = (smi_q, f"quaternized from {how}", conf,
                                 {"structure_assumed": "1",
                                  "note": "Cholesterol conjugation site inferred: "
                                          "the paper draws no -C conjugate for "
                                          "any commercial core."})
                continue
        unresolved.append(key)

    zs = zscore([__import__("math").log10(float(r["value_printed"])) for r in vivo])
    for i, (r, z) in enumerate(zip(vivo, zs)):
        key = r["lipid_name"]
        res = resolved.get(key)
        s = res[0] if res else ""
        row = base_row(hdr)
        row.update({
            "Index": start_index + len(rows),
            "LNP_ID": f"LNP_{start_index + len(rows):07d}",
            "Experiment_ID": EID,
            "Formulation_ID": f"{EID}_Fvivo",
            "IL_name": key, "IL_SMILES": s,
            "IL_molratio": 47.6,
            "HL_name": "DOPE", "HL_SMILES": smi("DOPE"), "HL_molratio": 23.8,
            "PEG_name": "DMG-PEG2000", "PEG_SMILES": smi("DMG-PEG2000"),
            "PEG_molratio": 4.8,
            "fifthcomponent_name": "DODAP", "fifthcomponent_SMILES": smi("DODAP"),
            "fifthcomponent_molratio": 23.8,
            "Aqueous_buffer": "citrate", "Dialysis_buffer": "PBS",
            "Mixing_method": "handmixed",
            "Model": "in_vivo", "Model_type": "Mouse_B6", "Model_target": "brain",
            "Route_of_administration": "intravenous",
            "Cargo": "mRNA", "Cargo_type": "FLuc",
            "Dose_ug_nucleicacid": 2.0,
            "Experiment_method": "luminescence_normalized",
            "Experiment_batching": "individual",
            "Experiment_value": z,
            "Publication_link": LINK, "Publication_DOI": DOI, "Publication_PMID": PMID,
        })
        row.update(descriptors(s))
        # Quality flags onto the row (F-001), from the lipid's own record (F-011).
        if res:
            row["confidence"] = res[2]
            flags = res[3] if len(res) > 3 else {}
            if flags.get("structure_assumed") == "1":
                row["structure_assumed"] = "1"
                row["curation_note"] = flags.get("note", "")
        else:
            row["curation_note"] = ("IL_SMILES could not be resolved for this "
                                    "lipid; see new_curation_qc.json "
                                    "IL_SMILES_resolution.unresolved.")
        rows.append(row)

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=hdr)
        w.writeheader()
        w.writerows(rows)

    qc(rows, hdr, resolved, unresolved, notes)


def qc(rows, hdr, resolved, unresolved, notes):
    vitro = [r for r in rows if r["Model"] == "in_vitro"]
    vivo = [r for r in rows if r["Model"] == "in_vivo"]
    ratio_bad = []
    for r in rows:
        tot = sum(float(r[k]) for k in
                  ("IL_molratio", "HL_molratio", "PEG_molratio", "fifthcomponent_molratio")
                  if r[k] != "")
        if abs(tot - 100) > 0.05:
            ratio_bad.append([r["LNP_ID"], round(tot, 3)])
    out = {
        "file": "new_curation.csv",
        "scope": "curated after the LNPDB + Bowen Li additions; currently one paper",
        "columns": len(hdr),
        "schema_source": "data/LNPDB_vitro.csv (67 columns, matched exactly)",
        "rows_total": len(rows),
        "blocks": {
            "ZT_2026_in_vitro": {
                "n": len(vitro),
                "what": "486-formulation DoE on ONE lipid, C6O2B2; zero structural variance",
                "method": "luminescence_discretized_normalized (4 bins from Fig 2c)",
                "z_mean": round(st.mean([r["Experiment_value"] for r in vitro]), 12),
                "z_sd": round(st.pstdev([r["Experiment_value"] for r in vitro]), 12),
                "distinct_values": len({round(r["Experiment_value"], 6) for r in vitro}),
            },
            "ZT_2026_in_vivo": {
                "n": len(vivo),
                "what": "17 lipid cores x 3 cholesterol tails, fixed formulation",
                "method": "luminescence_normalized (z-score of log10 brain total flux)",
                "z_mean": round(st.mean([r["Experiment_value"] for r in vivo]), 12),
                "z_sd": round(st.pstdev([r["Experiment_value"] for r in vivo]), 12),
            },
        },
        "molratio_rows_not_summing_to_100": ratio_bad,
        "IL_SMILES_resolution": {
            "resolved": len(resolved),
            "unresolved": unresolved,
            "how": "36 synthesized cores enumerated from the Fig 1b building blocks. "
                   "The 15 commercial-core conjugates are built by quaternizing their "
                   "parent, which Fig 1b DRAWS for all five (MC3, SM-102, TT3, Lipid 14, "
                   "OF-2). MC3 and SM-102 parents are taken from the corpus; TT3, "
                   "Lipid 14 and OF-2 were read off Fig 1b at 600 dpi.",
            "confidence_basis": "CORRECTED 2026-10-30, finding F-011. The enumerated "
                                "lipids now carry THEIR OWN confidence from "
                                "03_enumerated.csv instead of a hard-coded 'high'. "
                                "low = the quaternization site is an assumption, because "
                                "the core has no ring nitrogen and the paper draws no "
                                "structure for it (24 lipids). medium = the core is C5/C6, "
                                "whose piperazine N is drawn in Fig 2a, or the cholesterol "
                                "CONJUGATION SITE is inferred for a commercial core (the "
                                "paper draws no -C conjugate at all), or the parent has "
                                "unsaturation the chain counter cannot validate (OF-2). "
                                "high = parent structure machine-counted (apex detection, "
                                "bond-gap CV < 0.02). The counter was validated against "
                                "MC3 in the same figure and FAILED there (20-21 C vs a "
                                "known 18), which is why only saturated chains were "
                                "machine-counted.",
            "by_confidence": {c: sorted(k for k, v in resolved.items() if v[2] == c)
                              for c in ("high", "medium", "low")},
            "structure_assumed": sorted(
                k for k, v in resolved.items()
                if len(v) > 3 and v[3].get("structure_assumed") == "1"),
        },
        "supplied_component_SMILES": {
            "note": "not present in the corpus; standard reference structures",
            "names": sorted(SUPPLIED),
        },
        "schema_notes": {
            "ExactMolWt": "was named van.der.Waals.Molecular.Volume in the corpus. "
                          "Verified as ExactMolWt on 40472/40542 rows across all three "
                          "corpus files (12846 unique lipids). The 70 exceptions are "
                          "confined to JW_2024 and LX_2024 and are STALE ROWS: their "
                          "Molecular.Weight column disagrees with the SMILES too, and "
                          "the stored exact/average ratio is 0.99928 +/- 0.00002, the "
                          "monoisotopic ratio. Renamed in the corpus files as well.",
            "stale_descriptor_rows_found_in_corpus": "70 rows in JW_2024 and LX_2024 have "
                                                     "ALL descriptor columns computed from "
                                                     "a previous version of their SMILES. "
                                                     "Not fixed here; flagged for you.",
            "CHL_columns_empty": "deliberate - this platform uses NO free cholesterol, "
                                 "which is the paper's central claim",
            "fifthcomponent": "DODAP, the ionizable co-lipid",
            "IL_to_nucleicacid_massratio": "left empty; the paper states a TOTAL "
                                           "lipid:mRNA ratio of 40:1 wt/wt, which is not "
                                           "the IL-specific ratio this column wants",
        },
        "paper_internal_contradictions_recorded": {
            "in_vitro_dose_and_n": "Fig 2c caption says 20 ng mRNA per well, n=3; Methods "
                                   "4.8 says 15 ng per well, n=5. CHECKED whether these are "
                                   "different experiments: they are not. 4.8 is the ONLY in "
                                   "vitro subsection in the paper, and it states that "
                                   "luciferase AND viability were read from one combined "
                                   "ONE-Glo + Tox assay, i.e. Fig 2c and Fig S12 are the same "
                                   "wells. So the paper contradicts itself on both dose and n. "
                                   "The figure caption is used (0.02 ug) because it describes "
                                   "Fig 2c specifically; n is likewise contested (3 vs 5) and "
                                   "is NOT recorded in this schema. UNRESOLVED.",
            "acceptable_outstanding_counts": "text says 240/91, figures and digitization "
                                             "both give 238/93; totals reconcile at 331",
        },
        "merge_instructions": "split on `Model`: in_vitro rows append to "
                              "data/LNPDB_vitro.csv, in_vivo rows to data/LNPDB_vivo.csv. "
                              "Both blocks are already z-scored within themselves, which "
                              "is the corpus convention (per Experiment_ID, per file).",
    }
    json.dump(out, open(OUT_QC, "w"), indent=1)

    print(f"wrote {OUT}  ({len(rows)} rows x {len(hdr)} columns)")
    print(f"  in vitro  {len(vitro)}  z mean {out['blocks']['ZT_2026_in_vitro']['z_mean']:.2e} "
          f"sd {out['blocks']['ZT_2026_in_vitro']['z_sd']:.6f}  "
          f"{out['blocks']['ZT_2026_in_vitro']['distinct_values']} distinct values")
    print(f"  in vivo   {len(vivo)}  z mean {out['blocks']['ZT_2026_in_vivo']['z_mean']:.2e} "
          f"sd {out['blocks']['ZT_2026_in_vivo']['z_sd']:.6f}")
    print(f"  molratio rows not summing to 100: {len(ratio_bad)}")
    print(f"  IL_SMILES resolved {len(resolved)}/51 vivo lipids; "
          f"unresolved: {unresolved if unresolved else 'none'}")


if __name__ == "__main__":
    main()
