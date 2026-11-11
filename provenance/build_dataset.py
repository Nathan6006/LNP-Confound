#!/usr/bin/env python3
"""Build data/data_vitro.csv and data/data_vivo.csv, the finalized corpus.

LINEAGE
    original LNPDB (vitro)
      + the ML-work change log (see CHANGELOG.md)
      + the Bowen Li lab additions
      - three deletions REVERSED at the user's instruction
      = data_vitro.csv

SOURCES, and which one wins where:

  P = LNP_ML/new_data/LNPDB_vitro_del_processed.csv  (23194 rows, 83 cols)
      The post-change-log dataset. AUTHORITATIVE for IL_SMILES and
      Experiment_value, because it carries the manual corrections: iPhos
      switched from the charged form to neutral, the Lee_Unsat/SL_2020 E/Z
      fixes, and the LX_2024 SMILES fixes. All verified present here and absent
      from V (see validate()).
      It LACKS 21 curation columns (Model, Route_of_administration,
      Experiment_batching, the IL head/linker/tail decomposition,
      fifthcomponent_*, IL_protonated_SMILES).

  V = provenance/superseded/LNPDB_vitro.csv  (19028 rows, 67 cols)
      The curation schema. AUTHORITATIVE for those 21 columns, and the ONLY
      source for the three restored experiments.

ROW ASSEMBLY
  * every P row except YX_2024            21994
  * ALL of YX_2024, taken from V           2400   (restores 1200 RAW264.7)
  * ZC_2023 from V                          100   (restored)
  * XH_2024 from V                           30   (restored)
                                          ------
                                           24524

  Taking YX_2024 wholly from V is safe and was checked, not assumed: V's
  z-score over all 2400 rows is bit-identical to P's over its 1200 (affine fit
  slope 1.000000, intercept 0.000000, max residual 7e-15), i.e. P never
  re-z-scored after dropping RAW264.7. So the restored experiment is already
  correctly scaled. The 20 YX_2024 lipids whose SMILES P corrected are mapped
  onto the restored RAW rows through a V-canonical -> P-canonical table; RAW's
  SMILES set is a subset of HeLa's, so every restored row is covered.

DELETIONS KEPT (the user's call, recorded here so a reviewer can see it)
  * ZL_2022 (144) and YX_2023 (248): z-scored WITHOUT taking the log of RLU,
    so not comparable with the rest of the corpus. 392 rows.
  * 77 out-of-distribution cell-line rows: BMDC (36) + BMDM (12) inside
    LM_2019, HBEC_ALI (29) inside JW_2024.
  * 328 KL_2010 rows. NOT a change-log item - these simply have no
    Experiment_value in V. Reported separately so it is not mistaken for one.

COLUMNS: the 67 curation columns, then the 5 extra RDKit descriptors from the
change log. One-hot columns are deliberately NOT stored - see CHANGELOG.md.

Every descriptor column is RECOMPUTED from the canonical IL_SMILES rather than
carried over, which repairs the 70 stale rows in JW_2024/LX_2024 whose
descriptors were computed from a superseded structure.
"""
import csv
import json
import math
import io
import os
import statistics as st
import sys
from collections import Counter, defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors as D, rdMolDescriptors as rd

RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
P_PATH = os.path.join(ROOT, "..", "LNP_ML", "new_data", "LNPDB_vitro_del_processed.csv")
V_PATH = os.path.join(HERE, "superseded", "LNPDB_vitro.csv")
VIVO_PATH = os.path.join(HERE, "superseded", "LNPDB_vivo.csv")
OUT_VITRO = os.path.join(ROOT, "data", "data_vitro.csv")
OUT_VIVO = os.path.join(ROOT, "data", "data_vivo.csv")
OUT_QC = os.path.join(HERE, "dataset_validation.json")

# Curated papers appended after the corpus is assembled. Each is already in the
# 67-column schema and already z-scored within its own experiment blocks, so it
# is NOT renormalized again - only the 5 extra descriptors are computed and the
# rows are routed to vitro/vivo by their Model column.
MERGE = [os.path.join(HERE, "..", "curation", "new_curation.csv"),
         # Gong et al., Nat Biomed Eng 2025 (PMID 40910568-era record; DOI
         # 10.1038/s41551-025-01422-8). 333 rows: NG_2025, 310 in vitro on MC38
         # -- a Model_type new to the corpus -- and NG_2025_2, 23 in vivo in two
         # endpoint blocks, liver after IV and muscle after IM.
         #
         # Merged 2026-11-09 after a Tier-2 adversarial review (MERGE WITH
         # FIXES; report in provenance/review/round_2026-11-07/report_NG_2025.md).
         # The reviewer re-derived both z-vectors from the publisher's Source
         # Data workbook -- in vitro to max |delta| 5e-7, both in vivo blocks
         # EXACTLY -- and independently pinned the workbook-to-lipid mapping by
         # reducing each heat-map panel's vector fill colours to one principal
         # component and testing every column permutation (identity wins at
         # |rho| 0.998-0.9999 against a best competitor of 0.64-0.91). No
         # emitted structure and no emitted value was found wrong.
         #
         # It merges the output of stage11_quality.py, not stage07's, because
         # the folder was built against the 72-column schema: that stage adds
         # the 7 quality columns and repairs the Formulation_ID prefix on the 23
         # in-vivo rows. One row is deliberately held back by the curation (the
         # Fig. 1p benchmark, whose axis label and Source Data header name
         # different molecules).
         os.path.join(HERE, "..", "curation", "ng_2025",
                      "08_NG_2025_merge.csv"),
         os.path.join(HERE, "..", "curation", "lin_aromatic_2026",
                      "10_LL_2026_rows.csv"),
         os.path.join(HERE, "..", "curation", "lin_aromatic_2026",
                      "12_LL_2026_2_rows.csv"),
         # FALCON: all 2,544 author-reported measurements under one WT_2026
         # block, spanning seven cell lines - five of them new Model_type
         # categories (Ramos, THP-1, INS-1E, NIH-3T3, C2C12). An earlier version
         # merged only the subset on cell lines already in the corpus; that
         # restriction was withdrawn as indefensible - a cell line being new is
         # not a reason to discard a measurement, and Model_type is a feature,
         # not an inclusion criterion. See that folder's make_merge_subset.py.
         os.path.join(HERE, "..", "curation", "wht_2026_falcon",
                      "05_WT_2026_merge.csv"),
         # XD_2026: the 12-lipid in-vivo screen Duan et al. ran themselves
         # (Fig 5). The 485-row LipidDB in the same GitHub repo is deliberately
         # NOT here - it is a compilation across the 47 studies the paper cites,
         # so one Experiment_ID would fabricate a lab. See that folder's
         # stage06_lipiddb.py and 06_lipiddb_qc.json.
         os.path.join(HERE, "..", "curation", "xd_2026_lipidai",
                      "07_XD_2026_rows.csv"),
         # QL_2026: Liu et al., Adv Healthc Mater 2026 - 344 thiosuccinimide
         # lipids, HeLa, one pooled z after a 1.112x scalar offset tied the
         # six-tail library to the main library through their four shared
         # benchmark LNPs. The four benchmarks themselves (MC3, SM102, ALC-0315,
         # Lipid 5) are NOT emitted: they are calibration, not library members.
         os.path.join(HERE, "..", "curation", "ql_2026",
                      "11_QL_2026_rows.csv")]

# ------------------------------------------------------- publication metadata
# Each entry was verified against NCBI E-utilities before being applied; none is
# a guess. Surfaced by an independent manifest audit on 2026-10-29.
#
# SB_2024: the PMID is WRONG UPSTREAM - it is present in Original_LNPDB.csv and
#   points at "Multisociety Consensus ... Endovascular Therapy of Acute Ischemic
#   Stroke", Int J Stroke 2018. The experiment's own DOI resolves to PMID
#   38779592 (Bae SH, Bioact Mater 2024, trehalose-glycolipid LNP vaccine),
#   which matches the SB/2024 initials-and-year convention. DOI and PMID now agree.
# XH_2025: bioRxiv DOI carried a version suffix (...640222v1); the registered
#   DOI has no suffix, so a doi.org or [doi] lookup on the stored value fails.
# JW_2024: PMID and DOI agree with each other and with NCBI (Jiang AY,
#   Nat Nanotechnol 2024), but Publication_link pointed at an unrelated article
#   (s41587-024-02490-y). Link repointed at the DOI.
# MO_2026: Publication_link was not a URL at all - it read
#   "molea.csv (user-supplied HTS, ...)". Replaced with the DOI resolver.
METADATA_PATCHES = {
    "SB_2024": {"Publication_PMID": "38779592"},
    "XH_2025": {"Publication_DOI": "10.1101/2025.02.25.640222"},
    "JW_2024": {"Publication_link": "https://www.nature.com/articles/s41565-023-01548-3"},
    "MO_2026": {"Publication_link": "https://doi.org/10.1038/s41587-026-03109-0"},
}

# Experiments whose Experiment_value is NOT one z-score per experiment, declared
# explicitly so audit.py C12 can tell a known condition from a new one. Error
# class 22: standardizing each (experiment x cell type) block separately forces
# the between-assay offset to identically zero. It was found four times by hand
# -- ql_2026, WT_2026, and these two -- and never once by a check.
#
# WT_2026 is repaired at source (one z per source file, see that folder's
# build_rows.py) and declares its own unit on every row. These two cannot be
# repaired: both arrive as z-scores in the upstream P/V files with no raw signal
# anywhere in the repository, so the deleted offset is not recoverable from what
# we hold. It would take re-curating each paper's own supplementary data. So the
# honest move is to state it on the row rather than leave a reader to infer it.
ZSCORE_UNIT_DECLARED = {
    "YX_2024": ("experiment_x_model_type",
                "HeLa and RAW264.7 were z-scored separately upstream, so the "
                "between-assay offset is identically zero by construction and "
                "is NOT recoverable: this arrives already standardized in "
                "LNPDB_vitro.csv with no raw RLU in the repository. All 1,200 "
                "lipids appear in both assays, so cross-assay RANKING is "
                "usable (Spearman 0.461) but absolute potency between the two "
                "assays is not. See F-016."),
    "LZ_2025": ("experiment_x_model_type",
                "A549 and HepG2 were z-scored separately upstream; as YX_2024, "
                "the offset is not recoverable from the data we hold. See F-016."),

    # ---- F-080, added 2026-11-10. ----------------------------------------
    # The three entries below are standardized per BLOCK on a key that is NOT
    # Model_type, which is why C12 could not see them: it partitioned on
    # Model_type, and NO in-vivo experiment has more than one. C12 therefore
    # examined 0 of 20 in-vivo experiments while reporting clean.
    #
    # Per-block standardization is not an inference here. Three independent
    # lines, the third new:
    #   1. Each block sits at mean 0 / population sd 1 to the 9-dp ROUNDING
    #      FLOOR -- ZH_2023 liver |mean| 4.3e-11 against a predicted floor of
    #      2.3e-11, spleen 3.7e-11 against 2.3e-11, BL_2023 muscle (n=4) 2.8e-17
    #      with its four stored values summing to exactly 0 at 9 dp.
    #   2. The alternative -- that the raw blocks happened to share a mean AND
    #      an sd -- has a probability. Over the 160 lipids ZH_2023 measures in
    #      both organs the per-lipid liver-spleen difference has sd 1.1966, so
    #      the sample mean difference has SE 0.0946 even if the true offset were
    #      zero; P(|mean| < 8e-11) is about 7e-10, and the sd agreement is a
    #      second, independent coincidence of about 5e-10.
    #   3. ZH_2023's values are DISCRETE, and each organ lies on its own
    #      arithmetic lattice -- liver on 0.005823944 (127 unique values, every
    #      gap an exact integer multiple, max residual 1.6e-7), spleen on
    #      0.007970627 -- while the union of the two lies on no common lattice
    #      (residual 0.490). Three continuous controls return "no lattice", so
    #      the test is not manufacturing structure. Two different affine maps
    #      were applied. The ratio 1.3686 is what a single global z would have
    #      left between the two stored sds instead of 1.0000000000.
    #
    # NOT RESCALED, deliberately. Our values are byte-identical to
    # Original_LNPDB.csv (max |sorted(ours) - sorted(upstream)| = 0.000e+00) and
    # README.md's list of the eleven experiments we rescaled contains neither.
    # The per-block z is UPSTREAM, both papers are closed-access, and no raw
    # signal exists in the repository. Same posture as YX_2024/LZ_2025 above.
    #
    # NOT SPLIT into separate Experiment_IDs either, and that is the sharpest
    # reason to keep them whole: 160 of ZH_2023's lipids and all 4 of BL_2023's
    # muscle lipids appear in both of their own blocks, so splitting would
    # manufacture 164 FAKE cross-source shared lipids inside single papers and
    # feed them into Result 3 -- a pool that has only 2 real ones. It would also
    # contradict policy P4(c).
    "ZH_2023": ("experiment_x_model_target",
                "Liver and spleen were z-scored SEPARATELY upstream, so the "
                "between-organ offset is identically zero by construction and "
                "is not recoverable: this arrives already standardized in "
                "Original_LNPDB.csv with no raw luminescence in the repository. "
                "Model_target is the ONLY column that reproduces the partition "
                "(route, dose, cargo, method and mixing are all single-valued "
                "here). 160 lipids are measured in both organs, so cross-organ "
                "RANKING is usable (Spearman 0.241) but absolute potency "
                "between organs is not. See F-080.",
                "Model_target"),
    "BL_2023": ("experiment_x_model_target",
                "The two in-vivo arms were z-scored SEPARATELY upstream and the "
                "offset between them is not recoverable. The block is THE "
                "FIGURE ARM, not the organ as such: five columns reproduce the "
                "same partition exactly -- Model_target, "
                "Route_of_administration, Dose_ug_nucleicacid, Dialysis_buffer "
                "and Mixing_method (lung_epithelium/intratracheal/5ug/PBS/"
                "microfluidics against muscle/intramuscular/10ug/blank/"
                "handmixed) -- so the token names one of five co-varying keys "
                "and asserts no more than that. This experiment's 720 in-vitro "
                "rows are a THIRD separately standardized block. See F-080.",
                "Model_target"),
    "NG_2025_2": ("experiment_x_model_target",
                  "Liver (intravenous, Fig. 1m) and muscle (intramuscular "
                  "injection site, Fig. 1p) were z-scored separately by this "
                  "curation. The token was 'experiment_x_model_type' until "
                  "2026-11-10, which was FALSE on these rows: Model_type is "
                  "the single value Mouse_B6 and the block key is Model_target. "
                  "The curation_note had said so all along while the enum said "
                  "otherwise. See F-080.",
                  "Model_target"),
}

# The block key for each declared exception, where it is not Model_type. Used to
# BLANK Experiment_value_within_model rather than compute a meaningless one: that
# column re-standardizes inside (experiment x Model_type), and on these rows
# Model_type is constant, so it collapses to the identity -- byte-identical to
# Experiment_value on all 325 ZH_2023 rows. That satisfies the letter of its
# definition while defeating its stated purpose ("use this where a between-assay
# offset would be a nuisance"), because here the assay block is the organ and
# that offset is already gone. Re-standardizing the union of two already-
# flattened blocks is error class 41, and a blank is the honest value: the
# schema's own rule is that blank means "not stated", not "clean".
ZSCORE_BLOCK_KEY = {eid: t[2] for eid, t in ZSCORE_UNIT_DECLARED.items()
                    if len(t) > 2}

# A DECLARATION MUST NOT LEAK ACROSS POOLS. BL_2023 spans both files: its 53
# in-vivo rows are standardized per organ, and its 720 in-vitro rows are a
# THIRD, separately standardized block whose key is the pool, not the organ.
# Keying the declaration on Experiment_ID alone put `experiment_x_model_target`
# on all 773 -- a false statement about the 720 in-vitro rows, introduced by
# this very repair and caught by reading its own output. Those rows keep
# `experiment`, which already means experiment-x-pool for the nine experiments
# that span both files (see schema_quality.json).
ZSCORE_UNIT_POOL = {"ZH_2023": "in_vivo", "BL_2023": "in_vivo",
                    "NG_2025_2": "in_vivo"}

RESTORE = {"ZC_2023", "XH_2024"}
V_WINS = {"YX_2024"} | RESTORE          # experiments taken wholly from V
EXTRA = ["Num_carbon_in_tail", "num_protonatable_nitrogens",
         "num_permanent_cationic_N", "formal_net_charge",
         "num_unsaturated_cc_bonds"]

# Quality / provenance columns, 2026-10-30 (review findings F-001 and F-016).
# The names and the allowed values live in schema_quality.json and are read, not
# restated: the whole reason these columns exist is that per-curation flags were
# restated in six places, drifted, and then died one stage before the merge.
sys.path.insert(0, HERE)                 # so `import premerge_gate` resolves
import premerge_gate                      # the row-level merge gate

# Experiment_ID -> (folder, path) of whoever claimed it first, so a tag reused
# by a DIFFERENT paper is caught at merge time (policy P4b). Two files in one
# folder sharing a tag is normal -- a screen can span both pools.
_SEEN_EIDS = {}

QUALITY_SCHEMA = json.load(open(os.path.join(HERE, "schema_quality.json")))
QUALITY = [c["name"] for c in QUALITY_SCHEMA["columns"]]
SMILES_COLS = ["IL_SMILES", "IL_protonated_SMILES", "IL_head_SMILES",
               "IL_linker_SMILES", "IL_tail1_SMILES", "IL_tail2_SMILES",
               "IL_tail3_SMILES", "IL_tail4_SMILES", "HL_SMILES",
               "CHL_SMILES", "PEG_SMILES", "fifthcomponent_SMILES"]

AMIDE = Chem.MolFromSmarts("[NX3][CX3]=[OX1]")
ESTER = Chem.MolFromSmarts("[CX3](=O)[OX2H0][#6]")
CARBONATE = Chem.MolFromSmarts("[OX2][CX3](=O)[OX2]")
DISULFIDE = Chem.MolFromSmarts("[#16X2][#16X2]")

# ---------------------------------------------------------------------------
# Structure patches applied after canonicalization. Each one fixes a lipid whose
# stored SMILES collided with a sibling's, which is always a transcription bug:
# two names that differ in tail length cannot share a structure.
#
# SL_2020 / 4A3-SC6 ------------------------------------------------------------
# Held the SC8 (octyl) structure, identical to 4A3-SC8. The SL_2020 tail series
# is 6-carbon (2C/3C/4C/5N name the double-bond position), so SC6 is the
# saturated hexyl control and SC8 the octyl one. Every other SC6/SC8 pair in
# SL_2020 is already correct and differs by exactly 2 CH2 per tail, which fixes
# the expected mass: 4A3 has 4 tails, so SC6 = SC8 - 8 CH2 = 1467.21 - 112.22 =
# 1354.99. The patched structure hits 1354.99 exactly.
#
# KZ_2016 / *-SC3 (42 lipids) --------------------------------------------------
# SC3 is **n-PROPYL** (C3). Corrected 2026-11-10, review finding F-077.
#
# Every SC3 originally held the SC2 (isopentyl, C5) structure. The 2026-10-28
# patch replaced it with **neopentyl**, also C5, on this reasoning: "6A3-SC2 and
# 6A3-SC3 share MolWt 1876.692 but are recorded with DIFFERENT tail chain
# lengths, 4 and 3, so SC3 is a more branched isomer of the same C5 formula."
#
# THAT MolWt WAS THE BUG'S OWN OUTPUT. `individual_metadata.csv`'s MolWt is a
# derived column computed from the PRE-FIX isopentyl SMILES. Measured over all
# 1,502 upstream rows: 42 disagree with the MolWt recomputed from their own
# `main_data.csv` SMILES, all 42 are SC3, and every deficit is exactly
# -28.054 x Num_tails (2 CH2 per tail). The shared 1876.692 is therefore not
# evidence of isomerism, it is the fingerprint of the defect that
# `Zhou_dendrimer/note.txt` -- "fixed 42 smiles (SC3 tails were being
# represented as SC2)" -- says was repaired, and `temp.py` beside it names the
# intent in code: wrong_pattern "SCCC(C)C" -> correct_pattern "SCCC".
# The patch read the stale half of that sentence and discarded the live half.
#
# Five independent witnesses, none of which is the other's output:
#   1. `IL_tail1_SMILES` = `CCCS` (propane-1-thiol) on 42/42 rows, present in
#      `Original_LNPDB.csv` and so predating both temp.py and our patch.
#   2. The homologue mass ladder. Repaired, SC3 = 360.2083 for head 1A1 =
#      SC4 - 14.0157 exactly, making SC3..SC18 an unbroken +1 CH2 ladder with
#      token = carbon count. As stored, SC3 was isobaric with BOTH SC2 and SC5
#      and heavier than SC4 -- three anomalies at once.
#   3. `main_data.csv` carries n-propyl on 42/42 after temp.py ran.
#   4. Figure 1b of PMID 40644303 (the 2025 Siegwart paper) DRAWS every thiol in
#      this library; SC3 is an unbranched zigzag one bond shorter than SC4, and
#      no quaternary carbon appears anywhere in the SC series.
#   5. The 2016 SI anchors the nomenclature by HRMS and NMR: "1-octanethiol" for
#      SC8, 1-tetradecanethiol for SC14, and SAR bins written "SC5-8"/"SC9-12".
#      Five SI HRMS formulas match the corpus enumeration 5/5.
# Against them stood one stale derived number. SC1/SC2 are the two BRANCHED
# specials (iso-C4, iso-C5, isobaric with SC4/SC5) and SC19 is 2-phenylethane-
# thiol, an index not a count; SC3..SC18 are linear.
#
# THE REPAIR IS SMARTS-BASED AND COVERS EVERY IL SMILES COLUMN, for two reasons
# the 2026-10-28 patch got wrong:
#   * It edited only `IL_SMILES`. `IL_protonated_SMILES` still held the ORIGINAL
#     isopentyl structure on 42/42 rows -- byte-identical to the same head's SC2
#     -- so the corpus stated this one tail three different ways in a single
#     row. No composition check could see it: isopentyl and neopentyl are
#     isomers, so formula and exact mass agree.
#   * A literal string replace is traversal-dependent. The isopentyl pattern
#     `CC(C)CCS` appears in only 32 of the 42 protonated strings and `SCCC(C)C`
#     in 37; the other rows write the same tail a different way. A SMARTS query
#     matches 42/42 in both columns. After the repair the protonated form sits
#     exactly one proton above the neutral (+1.00728 Da, a proton less an
#     electron) for the first time.
# ZH_2023 / A4I18-1R2C70 (2 vivo rows) ----------------------------------------
#
# RESOLVED 2026-11-09, review finding F-021 -- found by C14 on its first run, in
# an experiment no reviewer had been assigned.
#
# The upstream structure's sterol cleaves to **C26H44O**, a 19-NOR sterol
# missing the C10 angular methyl, with no stereochemistry at all. Cholesterol is
# C27H46O. That is a connectivity error, not a stereochemistry one, so it is a
# different class from F-013's epimer.
#
# Settled from the paper, which we did not hold until the user supplied it. The
# library is a Ugi four-component reaction, so the `C` token is the CARBOXYLIC
# ACID component -- and SI "Section D. Detailed List of Ugi Reactants Used in
# the Study" (p. 34) draws C70 in full as **cholesteryl hemisuccinate**:
# HOOC-CH2CH2-C(=O)-O-cholesteryl, with both angular methyls present, the C8/
# C14/C17 hydrogens marked, the iso-octyl side chain, and the C3-O on a BOLD
# WEDGE (3-beta).
#
# Worth recording how nearly this went wrong in the other direction: Figure 2a's
# `-C` menu lists only C18, C18-1, C18-2, C18-3 and C23, the SI text layer
# mentions C70 exactly once (a data row in Table S1), and the word "cholesterol"
# appears in this paper only as a FORMULATION EXCIPIENT. On that evidence the
# reasonable conclusion was that no sterol belonged in this lipid at all. The
# structure is drawn in an image with no text layer, in an appendix section the
# table of contents names but the body never cross-references.
#
# Rebuilt rather than hand-edited: the corrected string is the esterification of
# the Ugi acid with the registry's cholesterol, run as a reaction so the ACID
# loses its OH and the alcohol's oxygen is retained -- nothing adjacent to a
# stereocentre is deleted, which is the F-013 hazard. Verified with audit.py's
# OWN C14 reaction object: the corrected lipid liberates C27H46O /
# HVYWMOMLDIMFJA-DPAQBDIFSA-N, matching reagent_registry.json exactly, where the
# upstream one liberated C26H44O / DZKYBYMQFBHTTM-UHFFFAOYSA-N.
#
# Formula C62H111N3O4 -> C63H113N3O4 (the restored C19 methyl).
# A COMPONENT patch, not a substring PATCH: the whole IL structure is
# replaced, and COMPONENT_PATCHES strips the name -- which matters here,
# because ZH_2023 is one of the two experiments whose IL_name arrives from
# upstream with a trailing space (F-022), so an exact-match key would
# silently never fire.
ZH_C70 = ("CCCCCCCC/C=C\\CCCCCCCCNC(=O)C(C(CC)CCCC)N(CCCN(C)C)C(=O)CCC(=O)O[C@H]1CC[C@@]2(C)C(=CC[C@H]3[C@@H]4CC[C@H]([C@H](C)CCCC(C)C)[C@@]4(C)CC[C@@H]32)C1")

# SMARTS patches. Keyed like PATCHES, but each rule is
# (query SMARTS, replacement SMILES, index of the replacement's attachment atom)
# and is applied to EVERY column that carries the ionizable lipid, so the
# neutral and protonated strings cannot disagree. Use these, not PATCHES, for
# anything where the stored traversal varies between rows (F-077).
SMARTS_PATCHES = {
    # neopentylthio (what the 2026-10-28 patch wrote) AND isopentylthio (the
    # original upstream bug, still live in IL_protonated_SMILES) both become
    # n-propylthio. Listing both is deliberate: the two columns were wrong in
    # two DIFFERENT ways, and a rule that only knew about one would have fixed
    # half the row and reported success.
    ("KZ_2016", "*-SC3"): [
        ("[#16X2][CH2][CX4]([CH3])([CH3])[CH3]", "CCCS", 3),   # neopentyl -> propyl
        ("[#16X2][CH2][CH2][CH1]([CH3])[CH3]", "CCCS", 3),     # isopentyl -> propyl
    ],
    # SL_2020 / 4A3-SC6: IL_SMILES was patched to hexyl in 2026-10-28 but
    # IL_protonated_SMILES still held the SC8 octyl structure, a +113.13 Da
    # disagreement with its own row -- visible to a mass check, and nothing was
    # checking. Octylthio -> hexylthio.
    ("SL_2020", "4A3-SC6"): [
        ("[#16X2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH3]", "CCCCCCS", 6),
    ],
}

PATCHES = {
}

# ---------------------------------------------------------------------------
# Whole-value component replacements, found 2026-10-29 while curating
# wht_2026_falcon. Unlike PATCHES these are not substring edits and they are not
# all on IL_SMILES, so they get their own table keyed by
# (Experiment_ID, column, component name).
#
# YZ_2022 / YZ_2024 / 18PG (360 rows) ------------------------------------------
# Stored as a symmetric 18:0/18:0 lipid (DSPG). Avanti's "18PG" is not
# symmetric. The paper that introduced this shorthand into the corpus - YZ_2022,
# Nat Commun 2022, PMC9310361 - defines it in its own words: "the anionic
# 1,2-dimyristoyl-sn-glycero-3-phosphate (14PA) and 1-STEAROYL-2-OLEOYL-sn-
# glycero-3-phospho-(1'-rac-glycerol) (18PG)". That is SOPG: sn-1 stearoyl 18:0,
# sn-2 oleoyl 18:1 cis-9. Formula goes C42H82O10P- (777.5651) -> C42H80O10P-
# (775.5495). 14PA and DDAB in the same panel were checked and are correct.
# YZ_2024 is not open access, but it reuses the platform exactly - the two
# experiments share all 360 formulations, the same single ionizable lipid and
# the same six helper lipids - so the same reagent is meant. If YZ_2024 is ever
# found to differ, drop its key here; nothing else depends on it.
#
# AC_2025 / DLin-MC3-DMA and KC2 (419 rows) ------------------------------------
# The two names carry each other's structures. By definition DLin-MC3-DMA is
# (6Z,9Z,28Z,31Z)-heptatriaconta-6,9,28,31-tetraen-19-yl 4-(dimethylamino)-
# butanoate, an ACYCLIC ESTER, and DLin-KC2-DMA is 2,2-dilinoleyl-4-(2-dimethyl-
# aminoethyl)-1,3-DIOXOLANE, a ring with no ester. AC_2025 stores the dioxolane
# under "DLin-MC3-DMA" (rings=1, ester=False) and the ester under "KC2"
# (rings=0, ester=True). Both are C43H79NO2 at exactly 641.6111 Da, which is
# why every mass, formula and descriptor check in the corpus passed over them:
# only connectivity distinguishes the two. The names are kept and the structures
# put back on the right name, which also restores agreement with the 2,194 MC3
# rows in the rest of the corpus.
MC3 = "CCCCC/C=C\\C/C=C\\CCCCCCCCC(CCCCCCCC/C=C\\C/C=C\\CCCCC)OC(=O)CCCN(C)C"
KC2 = "CCCCC/C=C\\C/C=C\\CCCCCCCCC1(CCCCCCCC/C=C\\C/C=C\\CCCCC)OCC(CCN(C)C)O1"
SOPG = ("CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC(O)CO)"
        "OC(=O)CCCCCCC/C=C\\CCCCCCCC")

COMPONENT_PATCHES = {
    ("ZH_2023", "IL_SMILES", "A4I18-1R2C70"): ZH_C70,
    ("YZ_2022", "HL_SMILES", "18PG"): SOPG,
    ("YZ_2024", "HL_SMILES", "18PG"): SOPG,
}

# ---------------------------------------------------------------------------
# NAME patches: the structure is right and the label on it is wrong.
#
# AC_2025 / MC3 <-> KC2, revised 2026-10-29 after checking the upstream source.
# The first fix here rewrote the STRUCTURES to match the names. That was
# backwards, and it silently moved every measured value onto the other molecule.
#
# github.com/alvinchangw/COMET is the source. Its master table
# (experiments/data_json/in_house_lnp_data_overall_new_full_without_pbae_
# NPratios_updated_09222023_npratios.json, 3,119 entries) stores each LNP as
# component SMILES plus labels and carries NO LIPID NAMES AT ALL. The names
# reached the corpus from LNP_ML/data_files/COMET_LANCE/convert_json_to_csv.py,
# whose IL_SMILES_TO_NAME dict is annotated "verified via RDKit MW" -- and MW
# cannot separate MC3 from KC2, which are both C43H79NO2 at 641.6111 Da. That
# dict maps the dioxolane to "DLin-MC3-DMA" and the acyclic ester to "KC2",
# which is the wrong way round.
#
# So the measured values belong to the STRUCTURES, which came straight from the
# source, and only the labels need correcting. Patching by structure rather than
# by name keeps every (structure, value) pair exactly as COMET recorded it.
# built lazily: canon() is defined further down
NAME_PATCHES = {("AC_2025", "IL"): {MC3: "DLin-MC3-DMA", KC2: "KC2"}}


def apply_name_patches(rows):
    """Relabel a component whose STRUCTURE is trusted and whose name is not.
    Keyed by canonical SMILES so the value a row carries never moves."""
    log = []
    for r in rows:
        for (eid, pre), mapping in NAME_PATCHES.items():
            if r["Experiment_ID"] != eid:
                continue
            want = {canon(k): v for k, v in mapping.items()}.get(
                r.get(f"{pre}_SMILES", ""))
            if want and (r.get(f"{pre}_name") or "").strip() != want:
                log.append({"experiment": eid, "column": f"{pre}_name",
                            "before": r[f"{pre}_name"], "after": want,
                            "smiles": r[f"{pre}_SMILES"]})
                r[f"{pre}_name"] = want
    return log


def apply_component_patches(rows):
    """Replace a component's whole SMILES where the stored one is the wrong
    molecule. Runs after canonicalization and before descriptors, so descriptors
    are computed from the corrected structure."""
    log = []
    for r in rows:
        for (eid, col, name), smi in COMPONENT_PATCHES.items():
            if r["Experiment_ID"] != eid:
                continue
            pre = col.split("_")[0]
            if (r.get(f"{pre}_name") or "").strip() != name:
                continue
            new = canon(smi)
            if r.get(col) != new:
                log.append({"experiment": eid, "column": col, "component": name,
                            "before": r[col], "after": new})
                r[col] = new
    return log


def apply_patches(rows):
    """Apply PATCHES and SMARTS_PATCHES, returning what changed.

    Two mechanisms, and the SMARTS one exists because the literal one was not
    enough (F-077):

      * PATCHES is literal string surgery on `IL_SMILES` only. Anchored '^'
        rules match at the string start, where the first tail is written
        chain-first. It is kept for rules where the traversal is known fixed.
      * SMARTS_PATCHES matches a substructure and so is independent of how the
        row happens to write it, and it is applied to EVERY column carrying the
        ionizable lipid. The 2026-10-28 patch edited `IL_SMILES` alone, leaving
        `IL_protonated_SMILES` holding the original wrong tail on 42/42
        `KZ_2016` SC3 rows and the pre-patch octyl tail on `SL_2020`/`4A3-SC6`.
        A literal replace would also have missed 10 of those 42, because
        `CC(C)CCS` appears in only 32 of the protonated strings.

    A patch that fixes one column and leaves a sibling column disagreeing is
    worse than no patch: it makes the row internally inconsistent while looking
    repaired, and no formula or mass check can see it when the wrong tail is an
    isomer of the right one.
    """
    log = []
    # Only the columns that carry the IONIZABLE LIPID itself. Helper lipid,
    # cholesterol and PEG columns are untouched -- a tail patch has no business
    # rewriting DSPC.
    IL_COLS = ["IL_SMILES", "IL_protonated_SMILES", "IL_head_SMILES",
               "IL_linker_SMILES", "IL_tail1_SMILES", "IL_tail2_SMILES",
               "IL_tail3_SMILES", "IL_tail4_SMILES"]
    for r in rows:
        key = None
        if (r["Experiment_ID"], r["IL_name"]) in PATCHES:
            key = (r["Experiment_ID"], r["IL_name"])
        elif (r["Experiment_ID"], "*-SC3") in PATCHES and r["IL_name"].endswith("-SC3"):
            key = (r["Experiment_ID"], "*-SC3")
        if key:
            before = r["IL_SMILES"]
            s = before
            for old, new in PATCHES[key]:
                if old.startswith("^"):
                    stem = old[1:]
                    if s.startswith(stem):
                        s = new + s[len(stem):]
                else:
                    s = s.replace(old, new)
            if s != before:
                r["IL_SMILES"] = canon(s)
                log.append({"experiment": r["Experiment_ID"],
                            "lipid": r["IL_name"], "column": "IL_SMILES",
                            "how": "literal",
                            "before": before, "after": r["IL_SMILES"]})

        skey = None
        if (r["Experiment_ID"], r["IL_name"]) in SMARTS_PATCHES:
            skey = (r["Experiment_ID"], r["IL_name"])
        elif ((r["Experiment_ID"], "*-SC3") in SMARTS_PATCHES
              and r["IL_name"].endswith("-SC3")):
            skey = (r["Experiment_ID"], "*-SC3")
        if not skey:
            continue
        for col in IL_COLS:
            before = r.get(col)
            if not before:
                continue
            m = Chem.MolFromSmiles(before)
            if m is None:
                continue
            hit = 0
            for smarts, repl, cp in SMARTS_PATCHES[skey]:
                q = _smarts(smarts)
                if q is None or not m.HasSubstructMatch(q):
                    continue
                n = len(m.GetSubstructMatches(q))
                m = Chem.ReplaceSubstructs(
                    m, q, Chem.MolFromSmiles(repl), replaceAll=True,
                    replacementConnectionPoint=cp)[0]
                # ReplaceSubstructs can leave the molecule needing a re-parse
                # before another query will match it correctly.
                m = Chem.MolFromSmiles(Chem.MolToSmiles(m))
                if m is None:
                    raise SystemExit(
                        f"SMARTS patch {skey} broke {r['IL_name']} {col}: the "
                        f"replacement did not re-parse. Refusing to emit a "
                        f"structure nobody has looked at.")
                hit += n
            if hit:
                after = canon(Chem.MolToSmiles(m))
                if after != before:
                    r[col] = after
                    log.append({"experiment": r["Experiment_ID"],
                                "lipid": r["IL_name"], "column": col,
                                "how": f"smarts x{hit}",
                                "before": before, "after": after})
    return log


_smarts_cache = {}


def _smarts(p):
    if p not in _smarts_cache:
        _smarts_cache[p] = Chem.MolFromSmarts(p)
    return _smarts_cache[p]


_canon_cache = {}


def canon(s):
    if not s or not s.strip():
        return ""
    if s not in _canon_cache:
        m = Chem.MolFromSmiles(s)
        _canon_cache[s] = Chem.MolToSmiles(m) if m else s
    return _canon_cache[s]


def n_protonatable(m):
    """Neutral non-amide amines, plus pyridine-type aromatic N.

    The aromatic clause is not cosmetic: imidazole-headed lipids carry one
    basic (pyridine-type) N and one non-basic pyrrole-type NH, and omitting it
    mis-scored every one of them. Reproduces the change-log column on 1200/1200
    unique structures.
    """
    am = {t[0] for t in m.GetSubstructMatches(AMIDE)}
    n = 0
    for a in m.GetAtoms():
        if a.GetSymbol() != "N":
            continue
        hv = sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1)
        if a.GetIsAromatic():
            if a.GetTotalNumHs() == 0 and hv == 2:
                n += 1
        elif a.GetFormalCharge() == 0 and a.GetIdx() not in am and hv < 4:
            n += 1
    return n


def longest_chain(m):
    """Longest path through acyclic carbons; reproduces Num_carbon_in_tail."""
    idx = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "C" and not a.IsInRing()]
    S = set(idx)
    adj = {i: [nb.GetIdx() for nb in m.GetAtomWithIdx(i).GetNeighbors() if nb.GetIdx() in S]
           for i in idx}
    best = 0
    for start in idx:
        stack = [(start, frozenset([start]))]
        while stack:
            u, vis = stack.pop()
            if len(vis) > best:
                best = len(vis)
            for v in adj[u]:
                if v not in vis:
                    stack.append((v, vis | {v}))
    return best


_desc_cache = {}


def descriptors(smiles):
    keys = ["Heavy.Atoms", "Rings", "Aromatic.Rings", "Rotatable.Bonds", "ExactMolWt",
            "Topological.Polar.Surface.Area", "Hydrogen.Bond.Donors",
            "Hydrogen.Bond.Acceptors", "LogP", "Molar.Refractivity",
            "Fraction.sp3.Carbons", "sp3.Carbons", "Nitrogen.Count",
            "Molecular.Weight", "has_ester", "has_carbonate", "has_disulfide"] + EXTRA
    if not smiles:
        return dict.fromkeys(keys, "")
    if smiles in _desc_cache:
        return _desc_cache[smiles]
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        out = dict.fromkeys(keys, "")
    else:
        out = {
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
            "sp3.Carbons": sum(1 for a in m.GetAtoms()
                               if a.GetSymbol() == "C" and str(a.GetHybridization()) == "SP3"),
            "Nitrogen.Count": sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N"),
            "Molecular.Weight": round(D.MolWt(m), 2),
            "has_ester": int(m.HasSubstructMatch(ESTER)),
            "has_carbonate": int(m.HasSubstructMatch(CARBONATE)),
            "has_disulfide": int(m.HasSubstructMatch(DISULFIDE)),
            "Num_carbon_in_tail": longest_chain(m),
            "num_protonatable_nitrogens": n_protonatable(m),
            "num_permanent_cationic_N": sum(
                1 for a in m.GetAtoms() if a.GetSymbol() == "N" and a.GetFormalCharge() == 1
                and sum(1 for nb in a.GetNeighbors() if nb.GetAtomicNum() > 1) == 4),
            "formal_net_charge": sum(a.GetFormalCharge() for a in m.GetAtoms()),
            "num_unsaturated_cc_bonds": sum(
                1 for b in m.GetBonds()
                if b.GetBondType() == Chem.BondType.DOUBLE and not b.GetIsAromatic()
                and b.GetBeginAtom().GetSymbol() == "C" and b.GetEndAtom().GetSymbol() == "C"),
        }
    _desc_cache[smiles] = out
    return out


def main():
    P = list(csv.DictReader(open(P_PATH)))
    V = list(csv.DictReader(open(V_PATH)))
    hdr67 = list(V[0].keys())
    out_cols = hdr67 + EXTRA + QUALITY
    v_by_id = {r["LNP_ID"]: r for r in V}

    # SMILES corrections P made, keyed by V's canonical form, so they can be
    # carried onto restored rows that never reached P.
    fix = {}
    for r in P:
        v = v_by_id.get(r["LNP_ID"])
        if v:
            a, b = canon(v["IL_SMILES"]), canon(r["IL_SMILES"])
            if a != b:
                fix[a] = b

    rows, provenance = [], Counter()
    for r in P:
        if r["Experiment_ID"] in V_WINS:
            continue
        v = v_by_id.get(r["LNP_ID"], {})
        row = {c: v.get(c, "") for c in hdr67}
        row.update({c: r[c] for c in hdr67 if c in r and r[c] != ""})
        row["IL_SMILES"] = r["IL_SMILES"]
        row["Experiment_value"] = r["Experiment_value"]
        if not row.get("Model"):
            row["Model"] = "in_vitro"
        if not row.get("Route_of_administration"):
            row["Route_of_administration"] = "in_vitro"
        rows.append(row)
        provenance["processed" if v else "processed_only(BowenLi)"] += 1

    for eid in sorted(V_WINS):
        for v in V:
            if v["Experiment_ID"] != eid:
                continue
            if not v["Experiment_value"].strip():
                continue
            row = {c: v.get(c, "") for c in hdr67}
            c0 = canon(row["IL_SMILES"])
            row["IL_SMILES"] = fix.get(c0, c0)
            rows.append(row)
            provenance[f"restored_from_V:{eid}"] += 1

    for row in rows:
        for c in SMILES_COLS:
            if row.get(c):
                row[c] = canon(row[c])
    patch_log = apply_patches(rows)
    component_log = apply_component_patches(rows)
    name_log = apply_name_patches(rows)
    for row in rows:
        row.update(descriptors(row["IL_SMILES"]))

    vitro_fix = renormalize_ddof(rows)
    rows.sort(key=lambda r: (r["Experiment_ID"], int(r["Index"]) if str(r["Index"]).isdigit() else 0))
    vivo = list(csv.DictReader(open(VIVO_PATH)))
    vivo_fix = renormalize_ddof(vivo)
    for row in vivo:
        for c in SMILES_COLS:
            if row.get(c):
                row[c] = canon(row[c])
    # ORDER CORRECTED 2026-11-09 (F-047). This path used to compute descriptors
    # INSIDE the canonicalisation loop above, i.e. BEFORE the component patches
    # -- the opposite of the vitro path twenty lines up, which patches first and
    # then computes. So any component patch touching `IL_SMILES` on a vivo row
    # left all 14 descriptor columns describing the PRE-patch molecule.
    #
    # It had never bitten because the only two COMPONENT_PATCHES were on
    # `HL_SMILES` (the 18PG expansion), which does not feed the IL descriptors.
    # F-021's cholesteryl-hemisuccinate repair is the first one on `IL_SMILES`
    # in vivo, and it surfaced immediately: `Molecular.Weight` stayed at 962.59
    # for a structure whose formula had changed to C63H113N3O4.
    #
    # C8 caught it, which is the system working -- but C8 compares descriptors to
    # the SMILES, so it can only ever catch this AFTER a patch has shipped wrong
    # numbers in an intermediate. The ordering is the actual fix.
    component_log += apply_component_patches(vivo)
    name_log += apply_name_patches(vivo)
    for row in vivo:
        row.update(descriptors(row["IL_SMILES"]))
    meta_fixed = apply_metadata_patches(rows) + apply_metadata_patches(vivo)
    merged = merge_curations(rows, vivo, out_cols)
    quality = {"vitro": add_quality_columns(rows),
               "vivo": add_quality_columns(vivo)}
    stripped = {"vitro": strip_cells(rows), "vivo": strip_cells(vivo)}
    rows.sort(key=lambda r: (r["Experiment_ID"],
                             int(r["Index"]) if str(r["Index"]).isdigit() else 0))
    vivo.sort(key=lambda r: (r["Experiment_ID"],
                             int(r["Index"]) if str(r["Index"]).isdigit() else 0))
    write(OUT_VITRO, out_cols, rows)
    write(OUT_VIVO, out_cols, vivo)

    qc = validate(rows, vivo, P, V, out_cols, provenance, fix)
    qc["cells_whitespace_stripped"] = stripped
    qc["structure_patches"] = patch_summary(rows, patch_log)
    qc["name_patches"] = {
        "n_rows_relabelled": len(name_log),
        "by_change": {f'{d["experiment"]}: {d["before"]} -> {d["after"]}': sum(
            1 for x in name_log if x["before"] == d["before"]
            and x["experiment"] == d["experiment"]) for d in name_log},
    }
    qc["component_patches"] = {
        "n_rows_changed": len(component_log),
        "by_component": {f"{d['experiment']}/{d['component']}": sum(
            1 for x in component_log if x["experiment"] == d["experiment"]
            and x["component"] == d["component"]) for d in component_log},
        "examples": {f"{d['experiment']}/{d['component']}":
                     {"before": d["before"], "after": d["after"]}
                     for d in component_log},
    }
    qc["merged_curations"] = merged
    qc["metadata_fields_patched"] = meta_fixed
    qc["renormalized_to_population_sd"] = {"vitro": vitro_fix, "vivo": vivo_fix}
    qc["quality_columns"] = quality
    json.dump(qc, open(OUT_QC, "w"), indent=1)
    report(qc)


def renormalize_ddof(rows):
    """Put every experiment on the corpus's population-sd (ddof=0) convention.

    Most experiments sit at population sd exactly 1.0000, but seven do not:
    SS_2018, LL_2025, SX_2025, BL_2023_3CR and KL_2010 were z-scored with the
    SAMPLE sd (ddof=1), and SP_2020 and LZ_2025 are short by two rows' worth of
    variance, consistent with rows being dropped AFTER normalization. The
    deviations run from 0.1% to 3.4%.

    Any experiment already centred at zero is rescaled to population sd 1. That
    is a pure affine map within an experiment, so it preserves every ordering
    and every relative spacing - a per-experiment z-score carries no more
    information than that - while making the convention uniform across the
    corpus. In a dataset built to study cross-laboratory comparability, a silent
    normalization difference is exactly the wrong thing to leave in.

    Experiments that are NOT z-scores at all are deliberately left alone: forcing
    them to mean 0 / sd 1 would disguise the problem rather than surface it. See
    the Spikevax note in validate().
    """
    by = defaultdict(list)
    for r in rows:
        try:
            by[r["Experiment_ID"]].append(r)
        except KeyError:
            pass
    fixed = {}
    for eid, rs in by.items():
        v = [float(r["Experiment_value"]) for r in rs]
        n = len(v)
        if n < 3:
            continue
        pop = st.pstdev(v)
        if abs(st.mean(v)) < 1e-6 and pop > 0 and abs(pop - 1) > 1e-9:
            k = 1.0 / pop
            for r in rs:
                r["Experiment_value"] = float(r["Experiment_value"]) * k
            fixed[eid] = {"n": n, "scale": round(k, 6),
                          "old_pop_sd": round(pop, 6), "new_pop_sd": 1.0}
    return fixed


def patch_summary(rows, log):
    by = defaultdict(list)
    for r in rows:
        by[(r["Experiment_ID"], r["IL_SMILES"])].append(r["IL_name"])
    dup = {}
    for (eid, smi), names in by.items():
        u = sorted(set(names))
        if len(u) > 1:
            dup.setdefault(eid, []).append(u)
    counts = Counter(e["experiment"] for e in log)
    return {"lipids_patched": {k: v for k, v in counts.items()},
            "n_rows_patched": len(log),
            "examples": log[:2],
            "duplicate_SMILES_remaining_by_experiment": dup or "none"}


def apply_metadata_patches(rows):
    """Correct verified publication-identifier errors. Returns rows changed."""
    n = 0
    for r in rows:
        patch = METADATA_PATCHES.get(r["Experiment_ID"])
        if not patch:
            continue
        for k, v in patch.items():
            if r.get(k) != v:
                r[k] = v
                n += 1
    return n


def open_decisions(path):
    """OPEN entries in the curation folder's DECISIONS.md, if it has one.

    THE MERGE REFUSES A FOLDER WITH AN OPEN DECISION. This is what makes a long
    unattended curation run safe: the pipeline is good at detecting that a
    source is silent, and it is not entitled to resolve it.

    The case for the gate is `AC_2025`, where the first repair was confident,
    internally consistent, and would have assigned **419 measured values to the
    wrong molecules** -- the NAMES were wrong and the structures were right,
    and nothing in the data said which half to fix. F-004 is the same shape:
    three mutually inconsistent `N:P` candidates (7.199 / 16.766 / 6.0), a
    reconstruction that was implemented, committed and then reverted, and a
    final decision to blank the column that was the user's to make and not the
    pipeline's.

    A folder with no DECISIONS.md is not gated -- the five curations merged
    before 2026-11-09 predate the convention, and retro-fitting a file to a
    merged folder would gate history rather than new work.
    """
    folder = os.path.dirname(os.path.abspath(path))
    dec = os.path.join(folder, "DECISIONS.md")
    if not os.path.exists(dec):
        return None
    open_ids = []
    current = None
    in_fence = False
    for line in io.open(dec, encoding="utf-8"):
        t = line.strip()
        # SKIP FENCED CODE BLOCKS. The scaffold's own DECISIONS.md documents
        # the entry format inside a ``` fence, and that example registered as
        # a real OPEN entry ("D-1  <one line>") -- so EVERY scaffolded folder
        # was permanently blocked by a phantom. A gate that always fires
        # trains people to ignore it, which is worse than no gate.
        if t.startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        # Strip markdown emphasis before matching. A folder wrote
        # `**STATUS: OPEN**` and this gate -- which only matched a bare
        # `STATUS:` -- reported ZERO open decisions for a folder with five.
        # That is a FALSE NEGATIVE in the one check standing between an
        # unresolved question and the corpus, and it is far worse than the
        # phantom entry fixed alongside it. Found 2026-11-10.
        bare = t.lstrip("*_# 	")
        if t.startswith("### "):
            current = t[4:].strip().lstrip("*_ ").rstrip("*_ ")
        elif bare.upper().startswith("STATUS:") and current:
            if "OPEN" in bare.upper().split("STATUS:", 1)[1]:
                open_ids.append(current)
            current = None
    return open_ids

def merge_curations(vitro, vivo, cols):
    """Append per-paper curations, routing each row by its Model column.

    These files are produced by curation/build_new_curation.py in the same
    67-column schema and are already z-scored within their own experiment
    blocks, so renormalize_ddof is deliberately not re-run over them. Index and
    LNP_ID are reassigned above the current maximum so a rebuild can never
    collide with the corpus.
    """
    out = {}
    nxt = 1 + max(int(r["Index"]) for r in vitro + vivo
                  if str(r["Index"]).isdigit())
    for path in MERGE:
        if not os.path.exists(path):
            out[os.path.basename(path)] = "NOT FOUND"
            continue
        opened = open_decisions(path)
        if opened:
            raise SystemExit(
                f"REFUSING TO MERGE {os.path.relpath(path)}: its folder's "
                f"DECISIONS.md has {len(opened)} OPEN entry(ies):\n  "
                + "\n  ".join(opened)
                + "\n\nResolve each one -- set STATUS to 'RESOLVED <date> by "
                  "<who>' and record the DECISION line -- before merging. The "
                  "pipeline detects that a source is silent; it does not get "
                  "to decide what the source meant. AC_2025's first repair was "
                  "confident, internally consistent, and would have put 419 "
                  "measured values on the wrong molecules.")
        # THE PRE-MERGE GATE (provenance/premerge_gate.py). open_decisions()
        # above gates the FOLDER's unresolved questions; nothing gated the
        # ROWS, and the routing line below sent any row with a blank `Model`
        # to the in-vivo pool. `lij_42130331` has 300 such rows, which would
        # have become 16.6% of data_vivo.csv in silence.
        violations = premerge_gate.check(path, known_experiments=_SEEN_EIDS,
                                         in_merge=True)
        if violations:
            raise SystemExit(
                f"REFUSING TO MERGE {os.path.relpath(path)}: it fails the "
                f"pre-merge gate on {len(violations)} count(s):\n  "
                + "\n  ".join(violations)
                + "\n\nFix the emitted file, not this gate.")
        added = Counter()
        for r in csv.DictReader(open(path)):
            row = {c: r.get(c, "") for c in cols}
            for c in SMILES_COLS:
                if row.get(c):
                    row[c] = canon(row[c])
            row.update(descriptors(row["IL_SMILES"]))
            row["Index"] = nxt
            row["LNP_ID"] = f"LNP_{nxt:07d}"
            nxt += 1
            # NO ELSE-BRANCH FALLTHROUGH. This line used to read
            #     (vitro if r["Model"] == "in_vitro" else vivo).append(row)
            # so a blank, misspelled or absent Model did not fail -- it
            # silently became an in-vivo row. The gate above now catches that
            # before we get here, and this raise is the second line of defence
            # so the two cannot drift apart.
            m = r.get("Model")
            if m not in ("in_vitro", "in_vivo"):
                raise SystemExit(
                    f"{os.path.relpath(path)}: Model={m!r} is neither "
                    f"'in_vitro' nor 'in_vivo'. A row that does not say which "
                    f"pool it belongs in has not been curated.")
            (vitro if m == "in_vitro" else vivo).append(row)
            _SEEN_EIDS.setdefault(
                r["Experiment_ID"],
                (os.path.basename(os.path.dirname(os.path.abspath(path))),
                 os.path.relpath(path)))
            added[f'{r["Experiment_ID"]}/{r["Model"]}'] += 1
        out[os.path.basename(path)] = dict(added)
    return out


def add_quality_columns(rows):
    """Populate the two normalization columns; validate the five flag columns.

    `Experiment_value_within_model` is computed here rather than in each
    curation because it is exactly derivable from what is already stored: a
    per-block z-score of an affine function of the log signal is the per-block
    z-score of the signal itself, so re-standardizing Experiment_value inside
    each (Experiment_ID x Model_type) block gives the same numbers as going back
    to the raw data would. The reverse is NOT true, which is the whole of F-016:
    once the per-cell-type z has been taken, the between-assay offset is gone
    and no amount of arithmetic brings it back.

    `value_normalization_unit` is stated on every row that has a z-score, from
    ZSCORE_UNIT_DECLARED where an experiment is a known exception and from the
    curation's own emitted value where it set one. Rows whose values are not
    z-scores at all (the vivo blocks floored at 0) are left blank in both
    columns rather than given a meaningless renormalization.

    The five flag columns are NOT invented here -- only checked. A flag has to
    come from the curation that knows the evidence; anything this function could
    infer would be a guess wearing a provenance column's name.
    """
    allowed = {c["name"]: c["allowed"] for c in QUALITY_SCHEMA["columns"]}
    by_exp = defaultdict(list)
    for r in rows:
        by_exp[r["Experiment_ID"]].append(r)

    log = {"experiments_on_zscore_scale": 0, "experiments_not_zscores": [],
           "within_model_populated": 0, "declared_units": {}, "bad_values": []}

    for eid, block in sorted(by_exp.items()):
        vals = []
        for r in block:
            try:
                vals.append(float(r["Experiment_value"]))
            except ValueError:
                vals.append(None)
        good = [v for v in vals if v is not None and math.isfinite(v)]
        on_z = (len(good) > 2 and abs(st.mean(good)) < 5e-3
                and abs(st.pstdev(good) - 1) < 0.05)
        if not on_z:
            log["experiments_not_zscores"].append(eid)
            continue
        log["experiments_on_zscore_scale"] += 1

        declared = ZSCORE_UNIT_DECLARED.get(eid)
        pool = ZSCORE_UNIT_POOL.get(eid)
        for r in block:
            # A DECLARED EXCEPTION OVERRIDES THE CURATION'S OWN TOKEN. It used
            # to only fill a BLANK, which let NG_2025_2 keep
            # `experiment_x_model_type` -- a token that is literally false on
            # those rows, since Model_type is the single value Mouse_B6 and the
            # block key is Model_target. Its curation_note had said "z-scored
            # SEPARATELY for liver ... and muscle" all along while the enum said
            # otherwise, and C12's declaration test is exact set equality
            # against one literal, so it read the false token as a pass. Where
            # this table and a curation disagree, this table is the adjudicated
            # answer and the override is logged.
            if pool and r.get("Model") != pool:
                if not r.get("value_normalization_unit"):
                    r["value_normalization_unit"] = "experiment"
                continue
            if declared:
                was = r.get("value_normalization_unit")
                if was and was != declared[0]:
                    log.setdefault("units_overridden", {}).setdefault(
                        f"{eid}: {was} -> {declared[0]}", 0)
                    log["units_overridden"][
                        f"{eid}: {was} -> {declared[0]}"] += 1
                r["value_normalization_unit"] = declared[0]
            elif not r.get("value_normalization_unit"):
                r["value_normalization_unit"] = "experiment"
        if declared:
            log["declared_units"][eid] = declared[0]
            for r in block:
                if pool and r.get("Model") != pool:
                    continue
                if not r.get("curation_note"):
                    r["curation_note"] = declared[1]

        # F-080: where the experiment's standardization block is NOT Model_type,
        # Experiment_value_within_model would re-standardize the union of blocks
        # whose offsets are already deleted -- error class 41, and on ZH_2023 it
        # comes out byte-identical to Experiment_value on all 325 rows. Leave it
        # blank: "not stated" is true, "identical to Experiment_value" is a
        # number that invites a consumer to use it for exactly the purpose the
        # schema says it is for, and which it cannot serve here.
        if eid in ZSCORE_BLOCK_KEY:
            pool = ZSCORE_UNIT_POOL.get(eid)
            hit = [r for r in block if not pool or r.get("Model") == pool]
            for r in hit:
                r["Experiment_value_within_model"] = ""
            log.setdefault("within_model_blanked", {})[eid] = len(hit)
            if len(hit) == len(block):
                continue
            block = [r for r in block if r not in hit]
            vals = []
            for r in block:
                try:
                    vals.append(float(r["Experiment_value"]))
                except ValueError:
                    vals.append(None)

        cells = defaultdict(list)
        for r, v in zip(block, vals):
            if v is not None and math.isfinite(v):
                cells[r["Model_type"]].append((r, v))
        for cell, pairs in cells.items():
            v = [x[1] for x in pairs]
            if len(v) < 2:
                continue
            mu, sd = st.mean(v), st.pstdev(v)
            if sd == 0:
                continue
            for r, x in pairs:
                r["Experiment_value_within_model"] = round((x - mu) / sd, 10)
                log["within_model_populated"] += 1

    for r in rows:
        for col, ok in allowed.items():
            if ok is None:
                continue
            if str(r.get(col, "")) not in ok:
                log["bad_values"].append((r["LNP_ID"], col, r.get(col)))
    return log


def strip_cells(rows):
    """Strip surrounding whitespace from every cell, in place. Returns the count.

    915 `IL_name` cells arrived from upstream padded with a trailing space --
    720 `BL_2023` vitro, 150 `ZH_2023` + 45 `BL_2023` vivo -- and in `ZH_2023`
    the padding tracks the organ exactly (every liver row padded, every spleen
    row not), which is the fingerprint of two source tables concatenated. No
    cell anywhere else in the corpus carries padding.

    The padding is a join hazard that no check could see: `"A1I5R22C8 "` and
    `"A1I5R22C8"` are different dict keys carrying the *same* SMILES, while C7
    only ever tested for a blank -- and it stripped before testing, so a padded
    value passed as present. Two comparisons in this file already called
    `.strip()` defensively, which is the workaround and not the fix (F-022).

    Called from main() BEFORE write() and validate(), not inside write(). Doing
    it inside write() emitted stripped rows while validate() went on reading the
    padded ones in memory, so the QC JSON named `A55-T3 ` as a duplicate-SMILES
    partner and the emitted corpus contained no such lipid -- the exact
    docs-disagree-with-data split this project keeps hitting.

    Verified safe before applying: stripping merges no two distinct names inside
    any experiment, maps no name onto a second SMILES, and leaves the duplicate
    `(Experiment_ID, IL_name, Model_target, Model_type)` key count unchanged at
    40 -- the padded/unpadded pairs in `ZH_2023` are one lipid measured in two
    organs, which is two legitimate rows.
    """
    n = 0
    for r in rows:
        for k, v in r.items():
            if isinstance(v, str) and v != v.strip():
                r[k] = v.strip()
                n += 1
    return n


def write(path, cols, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})


def validate(rows, vivo, P, V, cols, provenance, fix):
    out = {"columns": len(cols), "vitro_rows": len(rows), "vivo_rows": len(vivo)}
    out["row_provenance"] = dict(provenance)

    exp = Counter(r["Experiment_ID"] for r in rows)
    out["experiments"] = len(exp)
    out["rows_per_experiment"] = dict(sorted(exp.items()))

    # restorations and retained deletions
    pexp = Counter(r["Experiment_ID"] for r in P)
    out["restored"] = {
        "ZC_2023": exp.get("ZC_2023", 0),
        "XH_2024": exp.get("XH_2024", 0),
        "YX_2024_RAW264.7": sum(1 for r in rows
                                if r["Experiment_ID"] == "YX_2024"
                                and r["Model_type"] == "RAW264.7"),
        "expected": {"ZC_2023": 100, "XH_2024": 30, "YX_2024_RAW264.7": 1200},
    }
    out["deletions_kept"] = {
        "ZL_2022_and_YX_2023": {"present": exp.get("ZL_2022", 0) + exp.get("YX_2023", 0),
                                "expected": 0,
                                "why": "z-scored without log(RLU); not comparable"},
        "ood_cell_lines": {
            "present": Counter(r["Model_type"] for r in rows),
            "expected_absent": ["BMDC", "BMDM", "HBEC_ALI"],
            "found": [k for k in ("BMDC", "BMDM", "HBEC_ALI")
                      if any(r["Model_type"] == k for r in rows)],
        },
        "KL_2010_blank_values": {
            "note": "NOT a change-log deletion: 328 V rows have no Experiment_value",
            "kept": exp.get("KL_2010", 0), "expected": 50},
    }

    # z-score integrity, per experiment
    bad_z = []
    by = defaultdict(list)
    for r in rows:
        try:
            by[r["Experiment_ID"]].append(float(r["Experiment_value"]))
        except ValueError:
            bad_z.append(r["LNP_ID"])
    zstat = {}
    for k, v in by.items():
        mu, sd = st.mean(v), st.pstdev(v)
        zstat[k] = [round(mu, 6), round(sd, 6), len(v)]
    off = {k: v for k, v in zstat.items() if abs(v[0]) > 1e-3 or abs(v[1] - 1) > 0.05}
    out["zscore_per_experiment_mean_sd_n"] = zstat
    out["zscore_experiments_off_spec"] = off
    out["rows_with_unparseable_value"] = len(bad_z)

    # chemistry integrity
    bad_smiles = [r["LNP_ID"] for r in rows if not Chem.MolFromSmiles(r["IL_SMILES"])]
    out["unparseable_IL_SMILES"] = len(bad_smiles)
    out["canonical"] = all(r["IL_SMILES"] == canon(r["IL_SMILES"]) for r in rows[:4000])

    # change-log fixes still present
    sl21 = [r for r in rows if r["Experiment_ID"] == "SL_2021"]
    charges = Counter(sum(a.GetFormalCharge() for a in Chem.MolFromSmiles(r["IL_SMILES"]).GetAtoms())
                      for r in sl21)
    sl20 = [r for r in rows if r["Experiment_ID"] == "SL_2020"]
    name_by_smiles = defaultdict(set)
    for r in sl20:
        name_by_smiles[r["IL_SMILES"]].add(r["IL_name"])
    dup = {k: sorted(v) for k, v in name_by_smiles.items() if len(v) > 1}
    out["changelog_fixes_verified"] = {
        "SL_2021_iPhos_neutral": {"formal_charge_histogram": {str(k): v for k, v in charges.items()},
                                  "pass": set(charges) <= {0, 1}},
        "SL_2020_EZ_trans_present": {"n_trans_motifs": sum(r["IL_SMILES"].count("/C=C/") for r in sl20),
                                     "pass": any("/C=C/" in r["IL_SMILES"] for r in sl20)},
        "SL_2020_duplicate_SMILES": {"remaining": {",".join(v): k[:60] for k, v in dup.items()},
                                     "n_remaining": len(dup),
                                     "note": "28 T/C collisions in V are resolved; any listed "
                                             "here were already present in the processed file"},
    }
    out["SMILES_corrections_carried_from_processed"] = len(fix)

    # descriptors recomputed: how many rows changed vs what P/V stored
    # Compare against V, the curation file, which is where stale descriptors
    # live: 70 rows in JW_2024/LX_2024 had EVERY descriptor computed from a
    # superseded SMILES, and more rows drift because P corrected structures
    # without V being updated. P's own descriptors were already current.
    by_id = {r["LNP_ID"]: r for r in rows}
    stale = 0
    for v in V:
        n = by_id.get(v["LNP_ID"])
        if not n or not v.get("Molecular.Weight"):
            continue
        try:
            if abs(float(v["Molecular.Weight"]) - float(n["Molecular.Weight"])) > 0.05:
                stale += 1
        except ValueError:
            pass
    out["descriptor_rows_corrected_vs_curation_file"] = stale
    out["vivo_experiments"] = dict(sorted(Counter(r["Experiment_ID"] for r in vivo).items()))

    nz = {}
    byv = defaultdict(list)
    for r in vivo:
        byv[r["Experiment_ID"]].append(float(r["Experiment_value"]))
    for k, v in byv.items():
        if abs(st.mean(v)) > 1e-3 or abs(st.pstdev(v) - 1) > 0.05:
            nz[k] = {"n": len(v), "mean": round(st.mean(v), 4),
                     "sd": round(st.pstdev(v), 4),
                     "range": [round(min(v), 3), round(max(v), 3)],
                     "method": sorted({r["Experiment_method"] for r in vivo
                                       if r["Experiment_ID"] == k})}
    out["vivo_experiments_not_on_the_zscore_scale"] = {
        "found": nz,
        "verdict": "NOT FIXED and NOT DROPPED - flagged for the user. These are "
                   "ratios to a benchmark (luminescence_relative_to_Spikevax), "
                   "floored at 0, not z-scores. They are the vivo analogue of the "
                   "ZL_2022 / YX_2023 problem that got those two dropped from the "
                   "vitro corpus, and the same decision should probably apply.",
    }
    return out


def report(q):
    print(f"data_vitro.csv  {q['vitro_rows']} rows x {q['columns']} cols, "
          f"{q['experiments']} experiments")
    print(f"data_vivo.csv   {q['vivo_rows']} rows x {q['columns']} cols")
    print(f"  row provenance: {q['row_provenance']}")
    r = q["restored"]
    print(f"  restored: ZC_2023 {r['ZC_2023']}/100, XH_2024 {r['XH_2024']}/30, "
          f"YX_2024 RAW {r['YX_2024_RAW264.7']}/1200")
    d = q["deletions_kept"]
    print(f"  kept deleted: ZL_2022+YX_2023 present={d['ZL_2022_and_YX_2023']['present']} (want 0); "
          f"OOD lines found={d['ood_cell_lines']['found']} (want [])")
    print(f"  z-score off spec: {list(q['zscore_experiments_off_spec'])or 'none'}")
    print(f"  unparseable IL_SMILES: {q['unparseable_IL_SMILES']}   "
          f"all canonical: {q['canonical']}")
    f = q["changelog_fixes_verified"]
    print(f"  iPhos neutral: {f['SL_2021_iPhos_neutral']['pass']} "
          f"{f['SL_2021_iPhos_neutral']['formal_charge_histogram']}")
    print(f"  SL_2020 trans motifs: {f['SL_2020_EZ_trans_present']['n_trans_motifs']}")
    sp = q["structure_patches"]
    print(f"  structure patches: {sp['lipids_patched']} ({sp['n_rows_patched']} rows)")
    print(f"  duplicate SMILES remaining: {sp['duplicate_SMILES_remaining_by_experiment']}")
    print(f"  descriptor rows corrected vs curation file: "
          f"{q['descriptor_rows_corrected_vs_curation_file']}")
    rn = q.get("renormalized_to_population_sd", {})
    print(f"  renormalized to pop sd 1 -- vitro: {list(rn.get('vitro', {}))}")
    print(f"                              vivo:  {list(rn.get('vivo', {}))}")
    nz = q["vivo_experiments_not_on_the_zscore_scale"]["found"]
    print(f"  vivo NOT on z-score scale (flagged, untouched): {list(nz)}")
    print(f"  merged curations: {q.get('merged_curations')}")
    print(f"  publication-metadata fields patched: {q.get('metadata_fields_patched')} "
          f"({', '.join(METADATA_PATCHES)})")


if __name__ == "__main__":
    main()
