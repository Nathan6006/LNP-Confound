#!/usr/bin/env python3
"""
PROVENANCE: the 20 PubMed queries behind this file live in `provenance/` -
queries_batch1.txt (14) + queries_batch2.txt (6), verbatim, plus the raw and
de-duplicated PMID lists, the corpus map used for duplicate checking, and the XML
parsers. See provenance/README.md for search parameters and the sweep's known blind
spots. Those files were recovered from a temporary scratchpad on 2026-10-27; do not
delete them, they are the only record of the exact search.
"""

"""Merge auto-extracted PubMed metadata with manual triage annotations -> data.json.

Auto metadata comes from the literature sweep (recs3.json, produced by the
scratchpad parse); manual annotations below are hand-assigned after reading each
abstract. Run with -I from outside the scratchpad.
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data.json")

# The sweep's parsed output originally lived in a session scratchpad, which is
# ephemeral - by 2026-10-30 it was gone and this script could not run at all.
# `provenance/recs3_archived.json` is that metadata recovered from the data.json
# it had already produced, and `provenance/corpus_pmids_at_sweep.json` is the
# have-map, archived at sweep time. The recovery is faithful for every field the
# script reads EXCEPT `maxlib`, which in data.json is already the manual `lib`
# wherever one was assigned; `maxlib` is only a fallback for rows with no manual
# library size, so no output row changes. Verified by rebuilding and diffing:
# the only rows that moved were the ones whose ADDED/status entries changed.
SCRATCH = "/private/tmp/claude-501/-Users-nathanliu-Documents-LNP-confound/ebd1fdd9-2e2e-4284-8eb1-27c940f1a93a/scratchpad/lit/"
ARCHIVE_RECS = os.path.join(HERE, "provenance", "recs3_archived.json")
ARCHIVE_HAVE = os.path.join(HERE, "provenance", "corpus_pmids_at_sweep.json")

# ---------------------------------------------------------------- manual triage
# kind:  ivl = in-vitro ionizable-lipid structure library
#        ivv = in-vivo barcoded / pooled screen
#        form = formulation-space screen (few ILs, vary ratios/helpers)
#        comp = compiled/curated multi-study dataset
#        meth = platform/method paper (little usable SAR data)
# status: new | added | have | dup (same data already in corpus under another PMID)
#   added = curated BY US and merged into data/data_vitro.csv or data_vivo.csv.
#           Distinct from "have", which means it was already in LNPDB upstream.
# prio:   1 = add first, 2 = worth adding, 3 = marginal, 0 = do not add

# ------------------------------------------------------------------ ADDED BY US
# Papers curated in curation/<folder>/ and merged into the canonical corpus.
# `rows` is what actually landed, which is not always what the abstract claims:
# Lin advertises "over 900" but that is 840 (Fig 2e) + 147 (Fig 2a) = 987 across
# two screens run at DIFFERENT formulations, so they carry separate IDs.
ADDED = {
 "42634137": dict(expid="ZT_2026", folder="siegwart_42634137",
   rows_vitro=486, rows_vivo=51, added="2026-10-28",
   what="486-formulation DoE on one lipid (C6O2B2) in IGROV1 + 51 lipids x 3 "
        "cholesterol tails measured as in vivo brain flux",
   validation="7 independent checks; 4 marginal bar charts reproduced exactly"),
 "41689530": dict(expid="LL_2026,LL_2026_2", folder="lin_aromatic_2026",
   rows_vitro=987, rows_vivo=0, added="2026-10-28",
   what="840 aromatic ionizable lipids (Fig 2e) + 147 initial screen (Fig 2a), "
        "both IGROV1 luciferase",
   validation="HRMS 4/4 exact to 0.0001 Da; 48/48 printed hit rates reproduced. NOTE 2026-10-30 (F-006): that hit-rate check was weaker than stated until this date -- its pixel scale was fitted to the rates it validated, so a uniform scale error was invisible. Re-anchored on axis ticks; the 48 values themselves hold, max |delta| 1.0 pp"),
 # Preprint: no PMID, so this key never matches a tracker row. Recorded here so
 # the ADDED registry stays the single list of what we curated.
 "doi:10.64898/2026.06.13.731774": dict(expid="WT_2026", folder="wht_2026_falcon",
   rows_vitro=2544, rows_vivo=0, added="2026-10-29",
   what="FALCON closed-loop formulation optimization (Mao lab). 2,544 "
        "measurements over 1,218 formulations from the authors' own repo, "
        "7 cell lines, all merged under one WT_2026 block",
   validation="17 independent checks; all 5 published initial-library sizes "
              "(Suppl Table 9) reproduced; surfaced and fixed 2 corpus errors"),
 "42592089": dict(expid="XD_2026", folder="xd_2026_lipidai",
   rows_vitro=0, rows_vivo=12, added="2026-10-29",
   what="the paper's OWN in-vivo screen: E1-E9 plus ALC-0315/SM-102/MC3, "
        "BALB/c IV, whole-body IVIS luciferase at 6 h (Fig 5). The repo's "
        "485-row LipidDB is deliberately NOT merged - it is a compilation "
        "across the 47 studies the paper cites, so one Experiment_ID would "
        "fabricate a lab; extracted standalone in the curation folder",
   validation="structures traced geometrically from Fig 5A, not counted by eye; "
              "SM-102 and MC3 reproduce the corpus exactly while ALC-0315 is "
              "DRAWN one carbon short (benchmarks emitted with corpus SMILES, "
              "E1-E9 flagged confidence=medium); Fig 5C digitisation agrees with "
              "the independent Fig 5D class matrix 12/12"),
 "41145968": dict(expid="QL_2026", folder="ql_2026",
   rows_vitro=344, rows_vivo=0, added="2026-10-29",
   what="344 thiosuccinimide lipids in HeLa from a bispecific-linker library: "
        "8 amine cores x 8 NHS-ester/maleimide linkers x 5 linear tails (320, "
        "Fig 2b/2c) plus a 24-member six-tail library on the T3A core (Fig 3c). "
        "The 4 benchmark LNPs on the same panels are calibration, not library "
        "members, and are not emitted",
   validation="21 checks; all 3 of Fig 3d's printed hit rates, all 8 of Fig "
              "3e's per-linker rates and the '17 lipids >1e5' sentence "
              "reproduced; 29 of 32 Table S1/S2 masses agree as printed, and "
              "Table S2 holds the tail at 2C10 so all 8 linkers are pinned "
              "independently of tail length; the two libraries are tied by a "
              "1.1116x offset from their 4 shared benchmarks and z-scored as "
              "ONE pool, which preserves the paper's ~37x library difference"),
}

M = {
 # ---- P1 new in-vitro IL libraries -------------------------------------------
 "41689530": dict(kind="ivl", status="new", prio=1, lib=900, lab="Shuai Liu (Zhejiang Univ.)",
   readout="in vitro luciferase + in vivo IVIS",
   note="Largest brand-new IL library found. 900+ aromatic ionizable lipids, systematic aromatic-group SAR. Benchmarked against SM-102, so anchorable. Lab already partly in corpus (SW_2024, YY_2023) but this library is entirely new chemistry."),
 "38888972": dict(kind="ivl", status="new", prio=1, lib=580, lab="Kathryn Whitehead (CMU)",
   readout="in vitro + in vivo luciferase, cell-level flow",
   note="580 ionizable lipidoids, branched-tail series. Whitehead lab is in your corpus only via KW_2014 (2014). Pre-dates LNPDB's Aug-2025 cutoff yet was NOT curated - a genuine omission, not a new paper. Only ~40 worked in vivo, so in vitro is the usable layer."),
 "42634137": dict(kind="ivl", status="new", prio=1, lib=486, lab="Daniel Siegwart (UTSW)",
   readout="brain endothelial transfection; in vitro screen then in vivo",
   note="486 LNPs from a modular cholesterol-based lipid library. Siegwart is in corpus only via SL_2020/SL_2021. New lipid class (esterified cholesterol conjugates) = fully disjoint chemistry."),
 "42806111": dict(kind="ivl", status="new", prio=1, lib=384, lab="Bowen Li (Toronto)",
   readout="5.7-kb ABE-NanoLuc reporter (NEW endpoint)",
   note="384-lipid combinatorial library screened with a large (5.7 kb) cargo - a readout axis nothing in your corpus has. You already hold 4 datasets from this lab; ask them directly rather than curating from SI."),
 "41145968": dict(kind="ivl", status="new", prio=1, lib=344, lab="Jinghong Li (Tsinghua)",
   readout="in vitro luciferase + in vivo organ tropism",
   note="344 thiosuccinimide lipids across 2 libraries built on NHS-ester-maleimide bispecific linkers. Explicitly a LINKER-space library - orthogonal to the head/tail-space libraries that dominate your corpus."),
 "42130331": dict(kind="ivl", status="new", prio=1, lib=300, lab="Jinghong Li (Tsinghua)",
   readout="in vitro + spleen-tropic in vivo",
   note="300 beta-hydroxy thioether ILs. Same group as the thiosuccinimide library; the two together are ~644 new lipids from one lab, which also makes them a within-lab cross-library concordance test."),
 "42690843": dict(kind="ivl", status="new", prio=1, lib=275, lab="Heebeom Koo (Catholic Univ. Korea)",
   readout="purification-free in vitro screen, then in vivo",
   note="PILL platform: 275 alpha-branched amine ILs from a visible-light multicomponent reaction. Fully independent of the Anderson/Li lineage - the single most lineage-diverse addition available. MC3 benchmark included."),
 "40681859": dict(kind="ivl", status="new", prio=1, lib=None, lab="Michael Mitchell (Penn)",
   readout="in vitro + in vivo luciferase",
   note="Mannich-reaction chemistry, SIX combinatorial libraries screened in vitro and in vivo. Total size not stated in abstract - needs SI check, but 6 libraries implies several hundred. Published Dec 2025, just past LNPDB's cutoff."),
 "41851281": dict(kind="ivl", status="new", prio=1, lib=None, lab="Yao-Xin Lin (NCNST, Beijing)",
   readout="in vitro + organ-selective in vivo",
   note="Library spanning amino heads x degradable linkers x hydrophobic tails, with MD-derived 3D conformations converted to 2D density images for ML. ML-trained means per-lipid data must be tabulated. Chinese Academy of Sciences - lineage-independent."),
 "41203968": dict(kind="ivl", status="new", prio=2, lib=None, lab="Daniel Anderson (MIT)",
   readout="vaccine immunogenicity + in vivo expression",
   note="Degradable cyclic amino alcohol ILs via sequential combinatorial chemistry. Anderson lab, so same lineage as most of your corpus, but post-cutoff and not curated."),
 "42150080": dict(kind="ivl", status="new", prio=2, lib=None, lab="Daniel Siegwart (UTSW)",
   readout="hiPSC-derived cardiomyocyte transfection (NEW cell type)",
   note="First human iPSC-cardiomyocyte LNP screen. Library size not in abstract. Value is the cell type - your corpus has 12 cell lines, none cardiac, none hiPSC-derived."),
 "38982196": dict(kind="ivl", status="new", prio=2, lib=100, lab="Michael Mitchell (Penn)",
   readout="in vitro + in vivo luciferase",
   note="~100 amidine-incorporated degradable lipids, fast/facile synthesis. Nature Chemistry 2024, pre-cutoff, not in LNPDB - another genuine omission."),
 "29558113": dict(kind="ivl", status="new", prio=2, lib=288, lab="Pavel Levkin (KIT)",
   readout="in vitro transfection (plasmid/siRNA era)",
   note="288 lipids by one-pot thiolactone ring-opening. Older (2018) and non-mRNA, but a real 288-member library from a lab with zero corpus representation. Good for testing whether the ~0 cross-lab concordance holds in an even more distant chemical space."),
 "23909648": dict(kind="ivl", status="new", prio=3, lib=200, lab="Pavel Levkin (KIT)",
   readout="in vitro transfection, HTS",
   note="200 alkyl amine lipids, combinatorial + HTS. Early (2013), non-mRNA. Same lab as above."),
 "39105751": dict(kind="ivl", status="new", prio=2, lib=200, lab="Colin McKinlay",
   readout="in vitro mRNA transfection",
   note="~200-member peptoid-based nanoparticle library. Peptoid scaffold is structurally unlike every IL in your corpus."),
 # Demoted to P3 2026-10-26: verified n is 24, not 407 ("YK-407" is a lipid name).
 "40856486": dict(kind="ivl", status="new", prio=3, lib=407, lab="Gang Song",
   readout="spleen-targeted, deep-learning model",
   note="407 formulations with a deep-learning predictor for spleen targeting. DL-trained so the data must be tabulated somewhere. FLAGGED: the abstract says 'formulations', not lipids, so this may be a ratio sweep like the Mao-lab papers that were struck - count the distinct ILs in the SI before committing to it."),
 "40028722": dict(kind="ivl", status="new", prio=2, lib=216, lab="Ying Chau (HKUST)",
   readout="in vitro, immunomodulatory readouts",
   note="216 LNPs with ML-assisted design for dendritic-cell immunomodulation."),
 "39264000": dict(kind="ivl", status="new", prio=2, lib=213, lab="Jae-Hwan Nam (Catholic Univ. Korea)",
   readout="mRNA vaccine delivery, in vitro + in vivo",
   note="213-member rational-design series. Same group as 42446186 - pair them for a within-lab replicate."),
 "40500862": dict(kind="ivl", status="new", prio=2, lib=161, lab="Shuai Liu (Zhejiang Univ.)",
   readout="splenic NK-cell mRNA delivery",
   note="161 carbonate-bearing ILs. Same lab as the 900-lipid aromatic library - curate both and you get a within-lab, cross-library concordance pair for free."),
 "41691524": dict(kind="ivl", status="new", prio=3, lib=96, lab="Hyun-Jung Cha-Molstad (KRIBB)",
   readout="vaccine, IL + mRNA element co-engineering",
   note="96 LNPs, but the design varies mRNA elements as well as lipid - confounded by construct, harder to pool cleanly."),
 "42446186": dict(kind="ivl", status="new", prio=3, lib=56, lab="Jae-Hwan Nam (Catholic Univ. Korea)",
   readout="tail-modified ILs, safety + delivery",
   note="56 tail-modified ILs. Small; main value is as a within-lab replicate of 39264000."),
 "40689592": dict(kind="ivl", status="new", prio=3, lib=43, lab="Filip Du Prez (Ghent)",
   readout="in vitro transfection vs MC3 + in vivo",
   note="43 trialkyl galloyl amidoamine lipids, chromatography-free synthesis. Small but MC3-benchmarked and a completely independent European lab."),
 "42585110": dict(kind="ivl", status="new", prio=3, lib=40, lab="Shuai Liu (Zhejiang Univ.)",
   readout="spleen-selective antigen expression",
   note="40 retinoic-acid-derived ILs. Third Liu-lab library in this list."),
 "41775665": dict(kind="ivl", status="new", prio=3, lib=46, lab="Ming Wang",
   readout="mass-spectrometry-assisted HT discovery",
   note="46 lung-specific lipids found via an MS-assisted screen. Novel assay modality; small n."),
 "42489054": dict(kind="ivl", status="new", prio=3, lib=45, lab="Jing Zhang",
   readout="PEG-lipid diversification",
   note="45 PEG-lipids, not ionizable lipids. Only useful if you ever model the excipient scaffold, which Result 2 says is near-invariant."),
 "41617015": dict(kind="ivl", status="new", prio=3, lib=40, lab="Xin Zhang",
   readout="nose-to-brain, ML-guided",
   note="40 formulations, ML-guided optimisation. Small."),
 "42592089": dict(kind="comp", status="new", prio=1, lib=485, lab="Xiangrong Song (Sichuan Univ.)",
   readout="mRNA expression, ML framework",
   note="ML framework for rapid IL screening; library size not in abstract. Needs SI check - ML papers almost always ship the training table."),

 # ---- already in corpus but with NO Publication_PMID (matched by DOI only) ----
 # These three fell through both the PMID map and the title map because the corpus
 # identifies them only via Publication_link. Verified by DOI in this session.
 "41922838": dict(kind="ivl", status="have", prio=0, lib=960, lab="Bowen Li (Toronto)",
   readout="BEAS-2B in vitro luciferase", expid="FG_2026",
   note="ALREADY IN CORPUS as FG_2026 (960 rows, 960 lipids, BEAS-2B). Matched by doi 10.1038/s41563-026-02555-0 in Publication_link; the corpus row carries no Publication_PMID. 960 amino-acid-derived ionizable lipids, CHCha/CHGly series."),
 "41742414": dict(kind="ivl", status="have", prio=0, lib=1700, lab="Bowen Li (Toronto)",
   readout="HeLa (1,920 rows) + HBE (1,781 rows) in vitro luciferase", expid="YX_2026,YX_2026_SDL",
   note="ALREADY IN CORPUS as YX_2026 + YX_2026_SDL (3,701 rows combined). Matched by doi 10.1016/j.cell.2026.01.012; no Publication_PMID in the corpus. LUMI-lab self-driving platform, >1,700 LNPs autonomously synthesised and screened; 'SDL' = the self-driving-lab HBE arm."),
 "42050330": dict(kind="ivl", status="have", prio=0, lib=1152, lab="Bowen Li (Toronto)",
   readout="mouse primary hepatocyte + chondrocyte, mTP log2", expid="MO_2026",
   note="ALREADY IN CORPUS as MO_2026, and separately as data/molea.csv. This is MOLEA (doi 10.1038/s41587-026-03109-0). The corpus Publication_link literally reads 'molea.csv (user-supplied HTS, 1152-lipid screen)' with no PMID. NOTE: MO_2026 and molea.csv are the same 1,152 lipids - pooling both would double-count, which is why CLAUDE.md section 6 keeps molea separate."),

 # ---- corpus duplicates / corrections ----------------------------------------
 "42190011": dict(kind="ivl", status="dup", prio=0, lib=500, lab="Michael Mitchell (Penn)",
   readout="in vitro luciferase (HepG2)", expid="XH_2025",
   note="CORRECTION to my earlier recommendation: this is ALREADY in your corpus as XH_2025, curated from the bioRxiv preprint (PMID 40060499, doi 10.1101/2025.02.25.640222). The PNAS 2026 article is the same 500-lipid dataset. Do not re-add. Worth re-checking the SI against XH_2025's 500 rows in case the final version added lipids."),
 "40060499": dict(kind="ivl", status="have", prio=0, lib=500, lab="Michael Mitchell (Penn)",
   readout="in vitro luciferase (HepG2)", expid="XH_2025",
   note="The bioRxiv preprint LNPDB actually curated. Published version = PNAS 2026 (PMID 42190011)."),

 # ---- P2 in-vivo barcoded screens --------------------------------------------
 "29489381": dict(kind="ivv", status="new", prio=1, lib=281, lab="James Dahlman (Georgia Tech)",
   readout="DNA-barcode delivery, in vitro AND in vivo, same 281 LNPs",
   note="HIGHEST conceptual value in this whole list. 281 LNPs measured BOTH in vitro and in vivo by one lab, explicitly asking 'does in vitro predict in vivo?' - the answer was largely no. That is a same-lab, same-chemistry, cross-context concordance rung exactly like your molea hepatocyte-vs-chondrocyte point, and it is already published, so it is a citation rather than a curation job."),
 "35132167": dict(kind="ivv", status="new", prio=1, lib=89, lab="James Dahlman (Georgia Tech)",
   readout="barcoded delivery, mouse vs NHP",
   note="Species-dependent delivery: same 89 nanoparticles, different species, one lab. Another ready-made rung on the concordance ladder. Small n but the design is the point."),
 "40571690": dict(kind="ivv", status="new", prio=2, lib=45, lab="James Dahlman (Georgia Tech)",
   readout="barcoded aVHH mRNA, mouse + end-of-life NHP",
   note="45 LNP chemistries dosed in BOTH mice and NHPs; expression ranked differently between species. Too small to train on, ideal as a cross-species analogue of Result 3."),
 "30275336": dict(kind="ivv", status="new", prio=2, lib=250, lab="James Dahlman (Georgia Tech)",
   readout="functional mRNA delivery, barcoded (FIND)",
   note="250-LNP in vivo functional screen. Pre-cutoff and not in LNPDB."),
 "30216729": dict(kind="ivv", status="new", prio=3, lib=226, lab="James Dahlman (Georgia Tech)",
   readout="barcoded biodistribution",
   note="226 LNPs, receptor-retargeting study."),
 "41183194": dict(kind="ivv", status="new", prio=3, lib=109, lab="James Dahlman (Georgia Tech)",
   readout="barcoded, spleen vs liver",
   note="109 glycolipid nanoparticles, spleen-targeting without charge changes."),
 "41603124": dict(kind="ivv", status="new", prio=3, lib=122, lab="Michael Mitchell (Penn)",
   readout="barcoded mRNA (b-mRNA), in vivo only",
   note="122 LNPs (133-member library incl. 120 new plug-and-play ILs). VERIFIED data availability: 'available from the corresponding author upon reasonable request', SI is a single PDF. Low priority purely because the data is not deposited."),
 "41555650": dict(kind="ivv", status="new", prio=3, lib=None, lab="Michael Mitchell (Penn)",
   readout="barcoded mRNA, brain",
   note="Barcoded in vivo screen varying IL structure, IL mol%, and PEG mol%. Library size not in abstract."),
 "41269090": dict(kind="ivv", status="new", prio=3, lib=150, lab="J. Brandon Dixon (Georgia Tech)",
   readout="lymphatic endothelial cells, intradermal",
   note="150 LNPs over 3 screening rounds. These are formulation variants more than distinct ILs."),
 "36701517": dict(kind="ivv", status="new", prio=3, lib=137, lab="James Dahlman (Georgia Tech)",
   readout="barcoded + transcriptional response, lung",
   note="137 LNPs with paired transcriptomics."),
 "37028451": dict(kind="ivv", status="new", prio=3, lib=94, lab="James Dahlman (Georgia Tech)",
   readout="barcoded, high-throughput screens",
   note="94 LNPs, preferential delivery study."),
 "35671473": dict(kind="ivv", status="new", prio=3, lib=98, lab="James Dahlman (Georgia Tech)",
   readout="universal barcoding, ApoE-independent",
   note="98 LNPs."),
 "33758781": dict(kind="ivv", status="new", prio=3, lib=109, lab="James Dahlman (Georgia Tech)",
   readout="barcoded, liver immune cells",
   note="109 LNPs with constrained phospholipids."),
 "30394729": dict(kind="ivv", status="new", prio=3, lib=100, lab="James Dahlman (Georgia Tech)",
   readout="in vivo directed evolution, bone marrow",
   note="~100 LNPs."),

 # ---- P3 formulation-space screens -------------------------------------------
 "42220328": dict(kind="form", status="new", prio=2, lib=1944, lab="Hai-Quan Mao (Johns Hopkins)",
   readout="in vitro mRNA delivery + cellular immune response",
   note="1,944 formulations - the largest raw count in this list - but only EIGHT imidazolium cores; the rest is helper lipid x PEG density. Huge for formulation modelling, near-useless for IL-structure SAR. Include only if you add a formulation-space arm."),
 "38583366": dict(kind="form", status="new", prio=2, lib=1000, lab="Hai-Quan Mao (Johns Hopkins)",
   readout="liver gene editing, intraduodenal",
   note=">1,000 formulations, stepwise down-selection. Same caveat: formulation space, few ILs. Mao lab is in your corpus via YZ_2022/YZ_2024, both of which have exactly 1 unique IL - consistent with this lab's formulation-first style."),
 "41453807": dict(kind="meth", status="new", prio=3, lib=96, lab="David Issadore (Penn)",
   readout="in vitro transfection + physicochemical",
   note="Platform paper: parallelised microfluidics at ~1,000 formulations/hour, but only 96 were actually characterised. Interesting as a future data firehose, not as data today."),
 "42569578": dict(kind="form", status="new", prio=3, lib=864, lab="Ning Zhang",
   readout="siRNA knockdown",
   note="864 candidate formulations varying IL, helper, molar ratio and lipid concentration, with ML. siRNA endpoint; formulation-dominated."),
 "36133441": dict(kind="form", status="new", prio=0, lib=54, lab="Chun-Wan Yen (Genentech)",
   readout="oligonucleotide knockdown",
   note="54 PEGylated lipids. Excipient-space only."),

 # ---- P4 compiled datasets / benchmarks (positioning) ------------------------
 "42763081": dict(kind="comp", status="new", prio=1, lib=613, lab="Wang Sheng",
   readout="storage stability (particle size, EE)",
   note="MUST-CITE FOR YOUR PAPER. 613 LNP formulations pooled from 18 published studies. They ran nested, formulation-GROUPED cross-validation to kill repeated-measures leakage and found per-study ML models had 'no reliable predictive validity in any of the six adequately powered studies'. Different endpoint (storage stability, not transfection) but it is independent, published confirmation of your exact methodological thesis - honest grouped CV collapses pooled-literature LNP predictions."),
 "40267508": dict(kind="comp", status="new", prio=1, lib=6454, lab="Arezoo Ardekani (Purdue)",
   readout="activity + cell viability, QSAR",
   note="6,454 LNP formulations curated across 21 independent studies - a second, independent pooled corpus to LNPDB. Two uses: (a) cross-check your source list for publications LNPDB missed; (b) it is a live instance of pooled-literature modelling you can evaluate against, which is a much cleaner T11 target than re-running someone's architecture."),
 "41605942": dict(kind="comp", status="have", prio=0, lib=19528, lab="Evan Collins / Wonpil Im (Anderson lineage)",
   readout="n/a - the database itself", expid="(whole corpus)",
   note="This is LNPDB itself: 19,528 LNPs, 12,845 unique ILs, 42 publications + BroadPharm, counts stated 'as of August 2025'. Reconciles exactly to your 19,797-row Original_LNPDB.csv. Cite as the provenance of your corpus and as the source of the Aug-2025 cutoff."),
 "39658727": dict(kind="comp", status="have", prio=0, lib=9000, lab="Jacob Witten / Daniel Anderson (MIT)",
   readout=">9,000 LNP activity measurements",
   note="LION/LiON. Your local ~/LNP_ML repo is this project's codebase (github.com/jswitten/LNP_ML). Its >9,000 measurements are Anderson-lab historical screens, which are already inside LNPDB - so this is not a new source, it is the origin of much of your corpus."),
 "39554995": dict(kind="comp", status="new", prio=3, lib=7674, lab="Yuhua Hu",
   readout="patent-derived structures (no transfection values)",
   note="7,674 compounds profiled from siRNA-LNP patents. Structures without activity data, so it cannot feed Result 3 - but it is a large source of IL SMILES for chemical-space coverage arguments."),

 # ---- rejected ---------------------------------------------------------------
 "42916086": dict(kind="meth", status="new", prio=0, lib=None, lab="Justin Hanes (Johns Hopkins)",
   readout="hepatic delivery, PEG-lipid pairing",
   note="REJECTED: the '2000' my extractor caught was PEG2000, not a library size. This is a focused PEG/ligand-pairing study, not a library screen."),
}

# ---------------------------------------------------------------- second batch
# Found by mining the per-lipid `References` column of LipidDB
# (github.com/yuhanchenzz/LipidAI, the public dataset behind PMID 42592089) and by
# re-scanning the sweep for library-language papers that never state a size in the
# abstract - which is exactly why both LNPDB and my first pass missed them.
# Library sizes here are UNVERIFIED unless stated; abstracts do not give them.
M2 = {
 # ---- Virgil Percec (Penn): one-component ionizable amphiphilic Janus dendrimers
 # Demoted from priority 1 on 2026-10-26: IAJDs are not LNPs, so this cannot sit in
 # the "add first" tier even though it is the largest absent-lab contribution to LipidDB.
 "37376020": dict(kind="ivl", status="new", prio=2, lib=None, lab="Virgil Percec (Penn)",
   readout="in vivo organ-targeted mRNA, luciferase",
   note="'Screening Libraries to Discover Molecular Design Principles for Targeted Delivery of mRNA with One-Component IAJDs'. Percec's lab contributed 124 lipids to LipidDB - the single largest contribution from any lab absent from your corpus. CAVEAT, and it is a real one: IAJDs are sequence-defined Janus dendrimers used as ONE-COMPONENT carriers - no helper lipid, no cholesterol, no PEG - so these are not LNPs, the formulation columns would all be null, and calling them 'ionizable lipids' is a stretch. They do carry per-compound mRNA luciferase values, so IL_SMILES and Experiment_value would populate. Maximally disjoint support, which cuts both ways: good for the extrapolation argument, awkward for pooling. Your call, not mine."),
 "34324336": dict(kind="ivl", status="new", prio=2, lib=None, lab="Virgil Percec (Penn)",
   readout="in vivo mRNA delivery", note="Founding IAJD paper (JACS 2021). Start of a ~12-paper series (2021-2026) spanning JACS, Nat Commun, Sci Adv and Biomacromolecules, none of it curated anywhere in your files. Treat the series as one candidate programme rather than 12 separate adds."),
 "38306714": dict(kind="ivl", status="new", prio=2, lib=None, lab="Virgil Percec (Penn)",
   readout="total + targeted organ delivery",
   note="Constitutional-isomerism IAJD study (JACS 2024); contributed 13 lipids to LipidDB. Isomer series means near-identical structures with different activities - a useful stress test for fingerprint-based models."),

 # ---- Qiaobing Xu (Tufts): lipidoid libraries, ZERO corpus representation -----
 "32662132": dict(kind="ivl", status="new", prio=1, lib=None, lab="Qiaobing Xu (Tufts)",
   readout="ex vivo + in vivo mRNA to primary T lymphocytes",
   note="Imidazole-based lipidoids. Screened a broad library then a focused imidazole library with explicit structure-activity correlation. Qiaobing Xu is a top-5 LNP lab with NOT ONE source in your corpus - the most conspicuous single omission the whole search turned up. Primary T cells are also a cell type you do not have."),
 "32438474": dict(kind="ivl", status="new", prio=2, lib=None, lab="Qiaobing Xu (Tufts)",
   readout="protein + mRNA, in vitro and in vivo (muscle/lung/spleen/brain)",
   note="Cholesteryl-based biodegradable disulfide lipidoid library. Four organs via four routes."),
 "32520527": dict(kind="ivl", status="new", prio=2, lib=None, lab="Qiaobing Xu (Tufts)",
   readout="intracellular mRNA delivery, in vitro",
   note="Cyclic benzylidene acetal pH-responsive lipidoid library. The paper itself notes the lab had already built reduction- and photo-responsive libraries - i.e. there is a whole family of these."),
 "22148515": dict(kind="ivl", status="new", prio=2, lib=None, lab="Qiaobing Xu (Tufts)",
   readout="in vitro plasmid DNA transfection vs Lipofectamine 2000",
   note="Combinatorial lipidoid library from amines x acrylates/acrylamides/epoxides, with an explicit SAR section. CAVEAT: pDNA not mRNA, and benchmarked against Lipofectamine rather than MC3 - endpoint-comparability work needed before pooling."),
 "23651337": dict(kind="ivl", status="new", prio=3, lib=None, lab="Qiaobing Xu (Tufts)",
   readout="HeLa EGFP, plasmid DNA and mRNA",
   note="Unsaturated lipidoid library via Michael addition. Saturated-vs-unsaturated contrast on matched backbones. HeLa + EGFP, so partially comparable to your discretized-luminescence sources."),
 "24574196": dict(kind="ivl", status="new", prio=3, lib=None, lab="Qiaobing Xu (Tufts)",
   readout="siRNA knockdown, in vitro",
   note="Bioreducible disulfide lipidoid library, paired against non-bioreducible analogues."),

 # ---- Dan Peer (Tel Aviv) ------------------------------------------------------
 "31999380": dict(kind="ivl", status="new", prio=1, lib=None, lab="Dan Peer (Tel Aviv)",
   readout="siRNA gene silencing in leukocytes, in vitro + in vivo",
   note="'A Combinatorial Library of Lipid Nanoparticles for RNA Delivery to Leukocytes'. Novel ionizable amino lipids varying linker moiety (hydrazine / hydroxylamine / ethanolamine). Peer is another major lab with no corpus presence. Leukocytes are a hard-to-transfect cell type absent from your 12 cell lines."),
 "37092557": dict(kind="ivl", status="new", prio=1, lib=None, lab="Dan Peer (Tel Aviv)",
   readout="luciferase / mCherry / Cre mRNA, cell-type-specific",
   note="'A Combinatorial Library of LNPs for Cell Type-Specific mRNA Delivery'. Library varies hydrophobic tails AND linkers; identifies Lipid 23 (liver) vs Lipid 16 (CD11b+ cells). Contributed 7 lipids to LipidDB. Directly relevant to Result 2 (excipient/linker invariance)."),

 # ---- other labs surfaced by the LipidDB reference mining ---------------------
 "30637934": dict(kind="ivl", status="new", prio=1, lib=None, lab="Kathryn Whitehead (CMU)",
   readout="in vivo mRNA luciferase + apparent pKa",
   note="The 2019 branched-tail precursor to the 580-lipidoid 2024 paper (also on this list). Pairing the two gives a within-lab, cross-library concordance test on related chemistry - the same trick available with the Liu and Jinghong Li pairs. 2019, pre-cutoff, uncurated."),
 "36652629": dict(kind="ivl", status="new", prio=2, lib=None, lab="Omar Khan (Toronto)",
   readout="intramuscular mRNA delivery",
   note="'Iterative Design of Ionizable Lipids for Intramuscular mRNA Delivery' (JACS 2023); 7 lipids in LipidDB. Iterative design implies several rounds of per-lipid data. Intramuscular is a route your in-vitro corpus cannot speak to."),
 # PMID corrected 2026-10-26: the hand-entered 29307505 was a medical-cannabis paper.
 # Sabnis et al. Mol Ther 2018 is 29653760, and it is in PMC (PMC5986714).
 "29653760": dict(kind="ivl", status="new", prio=2, lib=None, lab="Kerry Benenato (Moderna)",
   readout="in vivo expression + tolerability",
   note="Sabnis et al., the SM-102 discovery paper - an INDUSTRIAL lipid series, a provenance class entirely missing from your corpus (10 lipids in LipidDB, plus Cornebise et al. for the follow-up amino lipid). Industry screens are a distinct environment and a strong test of the provenance-recognition claim. Library size not stated in the abstract - needs an SI check."),
 "26469188": dict(kind="form", status="new", prio=3, lib=200, lab="Daniel Anderson (MIT)",
   readout="in vivo mRNA, fractional factorial design",
   note="Kauffman et al. 2015: formulation optimisation by fractional factorial / definitive screening design. Anderson lab, pre-cutoff, uncurated. Formulation space, not IL space."),

 # ---- historical cationic-lipid libraries (pre-lipidoid era) ------------------
 "9465047": dict(kind="ivl", status="new", prio=3, lib=None, lab="Ronald Zuckermann (Chiron)",
   readout="in vitro plasmid transfection",
   note="PNAS 1998, 'A combinatorial approach to the discovery of efficient cationic peptoid reagents for gene delivery'. The oldest combinatorial transfection library found. Peptoid/lipitoid chemistry, pDNA endpoint - almost certainly not poolable, but it establishes that combinatorial delivery screening is 27 years old, which is useful framing for the Introduction."),
 "12210986": dict(kind="ivl", status="new", prio=3, lib=None, lab="Ulrich Massing (Freiburg)",
   readout="in vitro transfection, HTS",
   note="Chembiochem 2002, combinatorial synthesis of cationic lipids + HTS; companion to the 2001 J Biomol Screen method paper. Same framing value as above."),

 # ---- additional 2026 library papers the first pass missed --------------------
 "42610355": dict(kind="ivl", status="new", prio=2, lib=None, lab="Jee-Yeon Lee",
   readout="mRNA delivery to T cells at reduced lipid mass",
   note="'A double-amine ionizable lipid library' (Biomater Sci 2026). Two backbone nitrogens; T-cell endpoint. Size not stated in abstract."),
 "42886223": dict(kind="ivl", status="new", prio=3, lib=None, lab="Ming Zhao",
   readout="EGFP disruption, CRISPR-Cas9 RNP",
   note="HTS of an IL library for Cas9 RNP delivery. Only 4 authors and no size given - likely a modest library. RNP cargo is unlike anything in your corpus."),
 "40433897": dict(kind="ivl", status="new", prio=3, lib=None, lab="Zhen Zhang",
   readout="iterative in vitro screening",
   note="Vitamin-E-based functional LNPs, iterative screening (ACS Nano 2025). Size not stated."),
 "36839799": dict(kind="ivl", status="new", prio=3, lib=None, lab="Peng George Wang",
   readout="in vitro transfection",
   note="'Development of a Library of Disulfide Bond-Containing Cationic Lipids' (Pharmaceutics 2023). Open-access journal, so SI is likely reachable."),
 "39375194": dict(kind="form", status="new", prio=2, lib=None, lab="Hai-Quan Mao (Johns Hopkins)",
   readout="plasmid DNA LNP, ML-derived design features",
   note="ML on plasmid-DNA LNP design features (ACS Nano 2024). pDNA cargo; Mao-lab formulation-space style. ML-trained so a table should exist."),
 # Promoted to P1 2026-10-26: full text confirms 584 experimentally made ILs.
 "39664425": dict(kind="ivl", status="new", prio=1, lib=None, lab="Xiangrong Song (Sichuan Univ.)",
   readout="mRNA expression",
   note="'Pioneering integration of combinatorial chemistry and machine learning' (Acta Pharm Sin B 2024) - the predecessor to PMID 42592089, whose LipidAI dataset IS public on GitHub. Check whether this earlier paper's data is in the same repo."),
}

# ------------------------------------------------- access / data-availability grade
# A = machine-readable data file reachable right now
# B = SI document reachable, manual extraction needed
# C = nothing reachable from here (NOT a rejection - CWRU access may well work)
# D = confirmed unobtainable
ACCESS_OVERRIDE = {
 # User emailed the authors and was refused. This is the only hard rejection.
 "40267508": ("D", "CONFIRMED UNOBTAINABLE: data was curated but never published, and the "
                   "authors declined when asked directly (user's own correspondence, Oct 2026). "
                   "Keep only as evidence that pooled-literature LNP modelling is being done "
                   "without releasing data - do not plan to curate it."),
 "42592089": ("A", "Data and code public at github.com/yuhanchenzz/LipidAI (MIT licence): "
                   "'LipidDB-Ionizable lipid Library.xlsx' = 485 lipids with SMILES, ECFP, a "
                   "'Total expression' value and a per-lipid References column, plus "
                   "descriptors_1107.csv. VERIFIED by download. Note it is itself a "
                   "literature-compiled corpus, not new experiments - so it is a coverage "
                   "cross-check on LNPDB, and its reference column is a ready-made map of "
                   "sources (210 of its 218 cited publications are absent from your corpus)."),
 "41603124": ("C", "Data availability statement VERIFIED: 'available from the corresponding "
                   "author upon reasonable request'; SI is a single PDF. High refusal risk."),
 "42634137": ("C", "Data availability statement VERIFIED: 'upon reasonable request'. SI is one "
                   "PDF (ADMA-38-e74762-s001.pdf). Downgraded from priority 1 on this basis."),
 "39264000": ("C", "Data availability VERIFIED: 'upon reasonable request'. SI is a .docx."),
 "42446186": ("C", "Data availability VERIFIED: 'upon reasonable request'."),
 "42569578": ("B", "Mixed signal: statement says 'Data will be made available on request', but "
                   "an mmc2.xlsx IS deposited in the open archive. Check that spreadsheet first."),
 "40571690": ("A", "41587_2025_2711_MOESM3_ESM.xlsx deposited and reachable. Nature Portfolio "
                   "source-data policy. Only 45 LNPs, but mouse AND NHP for the same chemistries."),
 "39554995": ("A", "Four spreadsheets (mmc2-mmc5.xlsx) deposited; statement says all data are in "
                   "the text or supplementary materials. 7,674 patent-derived structures - "
                   "structures only, no activity values."),
 "40028722": ("B", "Data availability VERIFIED: included in the Supplementary Information (.docx)."),
}

# publications with no PMID (preprints / arXiv) -------------------------------
EXTRA = [
 dict(pmid=None, ident="arXiv:2609.17721", title="Decoding Extrahepatic Targeting of Lipid Nanoparticles with Interpretable Machine Learning",
      journal="arXiv", year="2026", first="", last="", doi="10.48550/arXiv.2609.17721",
      url="https://arxiv.org/abs/2609.17721", kind="comp", status="new", prio=1, lib=476,
      lab="(preprint)", readout="IVIS biodistribution, in vivo",
      note="476 intravenous LNP formulations curated from 81 separate studies, classifiers reported at ROC-AUC 0.839-0.874, with no cross-study or leave-one-study-out validation discussed. This is precisely the practice your paper critiques, published Sep 2026. Strongest available evidence that the failure mode is live and current, not historical."),
 dict(pmid=None, ident="arXiv:2507.03209", title="A Machine Learning Benchmarking Framework for Lipid Nanoparticle Transfection Efficiency Prediction",
      journal="arXiv", year="2025", first="", last="", doi="10.48550/arXiv.2507.03209",
      url="https://arxiv.org/abs/2507.03209", kind="comp", status="new", prio=2, lib=1100,
      lab="(preprint)", readout="HeLa transfection",
      note="Benchmarking framework built on 1,100 unique ILs taken from a SINGLE source (Xu et al. HeLa data, i.e. your YX_2024/AGILE). Claims to 'support assessment of model generalization' while evaluating within one lab's data. A second live instance of the gap, and a direct foil for your protocol-comparison figure."),
 dict(pmid=None, ident="BroadPharm", title="BroadPharm commercial ionizable lipid catalogue",
      journal="(vendor)", year="2024", first="", last="", doi="",
      url="https://www.broadpharm.com/", kind="comp", status="have", prio=0, lib=269,
      lab="(vendor)", readout="none - structures only", expid="BroadPharm_2024",
      note="269 ionizable lipids, structures only, no delivery values. Already present in your full_data_del+phys.csv as BroadPharm_2024 and correctly excluded from both delivery corpora."),
]

# ------------------------------------------------- what varies / what is measured
# Two axes kept separate from priority, so a row is never excluded for being a
# different KIND of experiment - it is labelled and left in. Defaults are
# varies="IL structure", measures="expression"; only the exceptions are listed.
#
# varies:   "IL structure"  the library sweeps the ionizable lipid
#           "formulation"   few ILs, many helper/chol/PEG/ratio combinations
#                           (the AC_2025 "comet_lance" shape: 1,978 rows, 7 ILs)
#           "both"          IL structure and formulation swept together - IL effect is
#                           confounded with formulation unless the design is crossed
# measures: "expression"           luciferase / GFP / protein abundance
#           "knockdown"            siRNA target reduction
#           "editing"              indel or base-edit efficiency
#           "barcode accumulation" deep-sequenced barcode abundance in tissue. Pooled
#                                  and competitive: all LNPs share one animal, so values
#                                  are within-pool relative, unlike one-LNP-per-well.
#                                  Treat as its own environment.
#           "functional barcode"   barcoded readout that requires translation
#                                  (Cre->tdTomato, aVHH protein), so expression-like
TAGS = {
    # pooled barcode-accumulation screens
    "30216729": (None, "barcode accumulation"),
    "36701517": (None, "barcode accumulation"),
    "41183194": (None, "barcode accumulation"),
    "33758781": (None, "barcode accumulation"),
    "30394729": (None, "barcode accumulation"),
    "35671473": (None, "barcode accumulation"),
    "37028451": (None, "barcode accumulation"),
    "29489381": (None, "barcode accumulation"),
    "35132167": (None, "barcode accumulation"),
    "41603124": (None, "barcode accumulation"),
    "41269090": ("formulation", "barcode accumulation"),
    "41555650": ("both", "barcode accumulation"),
    # barcoded but translation-dependent, so expression-like
    "30275336": (None, "functional barcode"),
    "40571690": (None, "functional barcode"),
    # formulation / ratio sweeps
    "42220328": ("formulation", None),
    "38583366": ("formulation", None),
    "39375194": ("formulation", None),
    "26469188": ("formulation", None),
    "42569578": ("formulation", "knockdown"),
    "41453807": ("formulation", None),
    # Reclassified once the real n was known: 24 bis-hydroxyethylamine-derived lipids
    # were synthesised, so this sweeps IL structure after all, not formulation.
    "41691524": ("both", None),          # mRNA elements swept alongside the lipid
    # knockdown endpoints
    "31999380": (None, "knockdown"),
    "24574196": (None, "knockdown"),
    "36133441": ("formulation", "knockdown"),
    # editing endpoints
    "42150080": (None, "editing"),
    "41203968": (None, "expression"),
}

# ------------------------------------------------- verified library sizes (2026-10-26)
# Re-derived for all 61 candidates. `lib` previously came from a regex over abstracts and
# had been hand-checked only for priority 1-2 rows; two values were plain wrong (see
# below). Each entry is (n, what_n_counts, how_verified).
#   abstract  = stated in the PubMed abstract, read in context
#   fulltext  = stated in the Europe PMC full text
#   download  = counted in a file actually retrieved
#   None      = not stated in any text reachable from here; needs the SI
# "what_n_counts" matters: a paper's synthesis count, screening count and formulation
# count routinely differ, and for pooling it is the number actually ASSAYED that counts.
SIZES = {
    # --- ionizable-lipid libraries, size stated in the abstract
    "41689530": (900, "ionizable lipids ('over 900')", "abstract"),
    "38888972": (580, "ionizable lipidoids synthesised and tested", "abstract"),
    "42634137": (486, "LNPs screened", "abstract"),
    "42806111": (384, "lipids in the combinatorial library", "abstract"),
    "41145968": (344, "thiosuccinimide lipids across 2 libraries", "abstract"),
    "42130331": (300, "distinct ionizable lipids", "abstract"),
    "42690843": (275, "a-branched amine ionizable lipids generated", "abstract"),
    "29558113": (288, "structurally diverse lipidoids", "abstract"),
    "40028722": (216, "LNP formulations screened", "abstract"),
    "39264000": (213, "LNPs analysed", "abstract"),
    "39105751": (200, "peptoid molecules ('>200') clustered", "abstract"),
    "40500862": (161, "carbonate-bearing ionizable lipids", "abstract"),
    "38982196": (100, "chemically diverse AID-lipids", "abstract"),
    "23909648": (200, "lipidoids synthesised and screened", "abstract"),
    "41691524": (96, "member biodegradable ionizable lipid library", "abstract"),
    "42446186": (56, "member tail-modified lipid library", "abstract"),
    "34324336": (54, "sequence-defined IAJDs across 6 libraries", "abstract"),
    "40689592": (43, "trialkyl galloyl amidoamine lipids", "abstract"),
    "42585110": (40, "member retinoic-acid-derived library", "abstract"),
    "38306714": (14, "constitutional isomeric IAJDs", "abstract"),
    # --- in vivo / barcoded screens, size stated in the abstract
    "29489381": (281, "LNPs, measured both in vitro and in vivo", "abstract"),
    "30275336": (250, "LNPs ('>250') in the FIND screen", "abstract"),
    "30216729": (226, "LNPs", "abstract"),
    "41269090": (150, "LNPs over 3 screening rounds", "abstract"),
    "36701517": (137, "LNPs, across 19 cell types", "abstract"),
    "41603124": (122, "LNPs in the b-mRNA screen", "abstract"),
    "41183194": (109, "chemically distinct LNPs", "abstract"),
    "33758781": (109, "LNPs with constrained phospholipids", "abstract"),
    "30394729": (100, "different nanoparticles ('>100')", "abstract"),
    "35671473": (98, "different LNPs, across 11 cell types", "abstract"),
    "37028451": (94, "chemically distinct nanoparticles", "abstract"),
    "35132167": (89, "lipid nanoparticles, mouse vs humanised/NHP", "abstract"),
    "40571690": (45, "LNP-aVHH mRNA chemistries pooled", "abstract"),
    # --- formulation sweeps, size stated in the abstract
    "42220328": (1944, "formulations over only 8 imidazolium cores", "abstract"),
    "38583366": (1000, "formulations ('over 1000'), few distinct ILs", "abstract"),
    "41453807": (96, "formulations actually characterised (1000/hr is throughput)", "abstract"),
    # --- resolved from Europe PMC full text
    "39664425": (584, "ionizable lipids built on the 4CR platform (the 40,000-lipid "
                      "library in the same paper is VIRTUAL, screened in silico only)",
                 "fulltext"),
    "37376020": (78, "nsSS IAJD molecules across two libraries", "fulltext"),
    "37092557": (24, "ionizable lipids (approximate - the text numbers them to 23 and "
                     "excludes 2 before in vivo; confirm in the SI)", "fulltext"),
    "36839799": (16, "ionizable lipids containing disulfide bonds", "fulltext"),
    # --- counted in a file I downloaded
    "42592089": (485, "lipids in LipidDB-Ionizable lipid Library.xlsx, compiled from 47 "
                      "published studies - a literature corpus, not new experiments",
                 "download"),
    # --- CORRECTIONS: both of these were wrong in the previous build
    "40856486": (24, "candidates synthesised. The '407' I had was the lipid NAME 'YK-407', "
                     "not a library size - exactly the false positive I had flagged as a "
                     "risk for this row", "abstract"),
    "42569578": (203, "formulations experimentally characterised, of which only 25 were "
                      "taken forward with siRNA. The '864' I had is the GENERATED library, "
                      "most of it never made - so the assayable n is 203 physicochemical "
                      "and 25 biological", "abstract"),
    # --- size not stated in any text reachable from here: needs the SI
    "40681859": (None, "six combinatorial Mannich libraries, total not stated", None),
    "41851281": (None, "amino heads x degradable linkers x tails, total not stated", None),
    "32662132": (None, "a broad lipidoid library then a focused imidazole library", None),
    "31999380": (None, "combinatorial library, size not stated", None),
    "30637934": (None, "branched-tail library plus a secondary library", None),
    "41203968": (None, "sequential combinatorial chemistry, total not stated", None),
    "42150080": (None, "'chemically diverse LNP library', total not stated", None),
    "32438474": (None, "cholesteryl disulfide lipidoid library, size not stated", None),
    "32520527": (None, "cyclic benzylidene acetal library, size not stated", None),
    "22148515": (None, "amines x acrylates/acrylamides/epoxides, total not stated", None),
    "36652629": (None, "iterative design over several rounds, total not stated", None),
    "42610355": (None, "double-amine library, size not stated", None),
    "39375194": (None, "pDNA LNP formulation space, size not stated", None),
    "23651337": (None, "unsaturated lipidoid library, size not stated", None),
    "24574196": (None, "bioreducible lipidoid library, size not stated", None),
    "40433897": (None, "iterative vitamin-E screening, size not stated", None),
    "26469188": (None, "fractional-factorial design; the '200' I had was the lipid name "
                       "C12-200", None),
    "41555650": (None, "barcoded brain screen; the '306' I had was the lipid name C14-306",
                 None),
}

GRADE_NOTE = {
    "A": "data file reachable now",
    "B": "SI document reachable; manual extraction",
    "C": "not reachable from here - try CWRU library access",
}

# ------------------------------------------------- unauthenticated access probe
# 2026-10-26, probe_access.py over all 46 candidates. Two routes, no credentials:
# download each Europe PMC supplementary archive and LIST it, and look for a preprint
# carrying the same SI. This corrected two flaws in the first grading pass, both of
# which had wrongly told the user to skip papers:
#   * SI presence was inferred from HTTP status. Many archives return 200 while being
#     ~300-byte EMPTY zips, and many non-empty ones hold only figure images.
#   * Elsevier and Wiley were skipped entirely - Elsevier names SI files mmc1.pdf with
#     no "supplement" in the path, and Wiley's downloadSupplement endpoint has no file
#     extension, so both failed the old link filter.
# Result: 13 of 46 archives have real contents, 10 of those are a usable SI document,
# and exactly ONE holds a spreadsheet - which turned out not to be transfection data.
ACCESS_PROBE = {
    "42634137": ("B", "SI IS FREE despite the verified 'upon reasonable request' statement: "
                      "ADMA-38-e74762-s001.pdf, 12.8 MB, in the Europe PMC archive. Upgraded "
                      "from C - no email needed. Tables are in the PDF, so extraction is "
                      "manual."),
    "30275336": ("B", "SI free: pnas.1811276115.sapp.pdf, 27.3 MB. Largest SI on the list."),
    "40028722": ("B", "SI free as a WORD file: IDRD_A_2465909_SM6840.docx, 4.8 MB. One of only "
                      "two .docx SIs found - Word tables convert far more cleanly than PDF, so "
                      "this is among the best extraction prospects here."),
    "39264000": ("B", "SI free as a WORD file: SMLL-21-2405618-s001.docx, 2.6 MB. Upgraded from "
                      "C. The other good extraction prospect."),
    "39105751": ("B", "SI free: nn4c05513_si_001.pdf, 1.6 MB."),
    "41269090": ("B", "SI free: nn5c12080_si_001.pdf, 2.1 MB."),
    "36701517": ("B", "SI free: nl2c04479_si_001.pdf, 12.5 MB."),
    "41183194": ("B", "SI free: pnas.2409569122.sapp.pdf, 1.8 MB."),
    "33758781": ("B", "SI free: BTM2-5-e10161-s001.pdf, 3.2 MB."),
    "41603124": ("B", "PREPRINT ROUTES AROUND THE PAYWALL: bioRxiv 10.1101/2025.06.09.658602 "
                      "(Hamilton et al., Jun 2025) serves its SI free as a 1.24 MB PDF, so the "
                      "published version's 'data upon reasonable request' statement does not "
                      "bind. Upgraded from C. Published SI (ADMA-38-e14370-s001.pdf, 2.1 MB) is "
                      "also in the Europe PMC archive."),
    "39375194": ("B", "PREPRINT: bioRxiv 10.1101/2023.12.07.570602 (title match 0.93), SI free "
                      "as a 2.15 MB PDF. Upgraded from C."),
    "42569578": ("B", "The ONLY machine-readable SI file the probe found anywhere - mmc2.xlsx - "
                      "and I opened it to check rather than trusting the extension. It is NOT "
                      "transfection data: ~48 rows of size / PDI / zeta potential / Cy3-siRNA "
                      "fluorescence over five COMMERCIAL lipids (ALC-0315, SM-102, MC-3, DSPC, "
                      "CEE-E12). Not the 203 formulations the paper describes, and no knockdown "
                      "values. mmc1.docx (2.7 MB) may hold more."),
    "39664425": ("C", "Europe PMC archive has contents but they are figure images only (gr1.jpg, "
                      "gr1.gif). The 584 lipids are not in any reachable file."),
    "42592089": ("A", "Data public at github.com/yuhanchenzz/LipidAI (MIT), verified by "
                      "download: 485 lipids with SMILES, ECFP, a Total expression value and a "
                      "per-lipid References column. NOTE the Europe PMC archive for this paper "
                      "is 40 files of figure images only - the GitHub repo is the data, and it "
                      "is a literature-compiled corpus, so deduplicate against your own before "
                      "use."),
}

# Grade D = confirmed unobtainable. Dropped from data.json entirely so the site
# carries only A-C, but kept listed here so a later sweep cannot silently
# re-propose them - the builder re-drops anything in this set.
REJECTED = {
    "40267508": "Ardekani: data curated but never published, and the authors declined when "
                "asked directly (user's own correspondence, Oct 2026).",
}

# ------------------------------------------------- inclusion triage (2026-10-26)
# A row is a candidate DATA SOURCE only if it clears all three tests:
#   1. per-LNP quantified delivery of a NUCLEIC ACID (mRNA / siRNA / pDNA - all three
#      are present in the corpus: Cargo = mRNA 14,195, siRNA 3,641, pDNA 1,192),
#      read out as expression, knockdown or editing;
#   2. not already in the corpus (re-verified: 0 of 74 match by PMID, DOI or
#      >=85% title similarity).
# Anything that fails is removed from the candidate list. Rows that fail test 1 but
# still matter to the manuscript become CITE_ONLY instead of vanishing.
#
# 2026-10-26, corrected after user review. I had added a third test - "enough distinct
# ILs to carry IL-structure signal" - and struck formulation sweeps and barcode screens
# under it. Both strikes were wrong:
#   * Formulation sweeps are already first-class in this corpus. AC_2025 ("comet_lance")
#     is 1,978 rows on SEVEN unique ILs (283 rows per IL) and is the largest in-vitro
#     source by that ratio; YZ_2022 and YZ_2024 are 1,080 rows on ONE IL each. A test
#     that excludes new formulation sweeps would exclude the corpus's own biggest ones.
#   * Barcodes are nucleic acids and barcode screens are nucleic-acid delivery. The
#     corpus already pools five different Experiment_method quantities, and the paper's
#     whole claim is that quantities do not rank-correlate across environments - so a
#     deliberately different readout is evidence, not contamination.
# Both classes are reinstated. Instead of excluding them, every row now carries
# `varies` (what the library sweeps) and `measures` (what the assay quantifies), so the
# heterogeneity is explicit and filterable at analysis time rather than hidden by a
# curation decision. The one caveat worth carrying: pooled barcode screens co-administer
# all LNPs in one animal, so each value is competitive and within-pool, unlike
# one-LNP-per-well - that is a reason to model them as their own environment, not to drop
# them.

STRUCK = {
    # --- cargo is not a nucleic acid, or there is no transfection readout at all
    "41617015": "Not an LNP paper: nasal-spray olfactory deposition, ML-guided. Swept in by "
                "keyword overlap only. No lipid library, no transfection readout.",
    "41775665": "Cargo is PROTEIN, not a nucleic acid ('In Vivo Protein Delivery'). "
                "Not poolable with an mRNA/siRNA/pDNA expression corpus.",
    "42886223": "Cargo is Cas9 ribonucleoprotein (protein + gRNA complex), read out as EGFP "
                "disruption. Neither the cargo class nor the editing readout is comparable to "
                "the corpus's expression values. Reversible if an editing arm is ever added.",
    "39554995": "Structures only - 7,674 patent compounds with NO transfection values. "
                "Cannot feed any delivery analysis. (Still a large source of IL SMILES if a "
                "chemical-space-coverage figure ever needs one.)",
    # --- the ionizable lipid is held constant; only the excipient varies
    # These are not formulation sweeps in the AC_2025 sense (which sweeps ratios around a
    # few ILs and still carries IL identity). Here the swept variable is the PEG-lipid
    # itself, so there is no ionizable-lipid axis at all. Reinstate if a PEG-lipid arm is
    # ever added - the corpus does have PEG_name / PEG_SMILES / PEG_molratio columns.
    "36133441": "54 PEGylated lipids with the ionizable lipid held constant: no IL axis.",
    "42489054": "45 PEG-lipids, ionizable lipid held constant. Same reason as 36133441.",
    "42916086": "Not a library at all. The '2000' my size extractor caught was PEG2000. "
                "Focused PEG/ligand-pairing study.",
}

# Fails the data-source tests but earns a place in the manuscript.
CITE_ONLY = {
    "42763081": "Endpoint is STORAGE STABILITY (particle size, encapsulation efficiency), not "
                "transfection - so it is not addable data. But it is the single most important "
                "citation available to you: 613 formulations from 18 studies, nested "
                "formulation-GROUPED CV, and 'no reliable predictive validity in any of the six "
                "adequately powered studies'. Independent published confirmation of your thesis "
                "on a different endpoint.",
    "9465047":  "PNAS 1998 cationic peptoids, pDNA: pre-LNP chemistry, not poolable. Cite for "
                "the Introduction - combinatorial transfection screening is 27 years old, so "
                "the field's pooling problem is not new.",
    "12210986": "Chembiochem 2002 combinatorial cationic lipids + HTS. Same framing value.",
}

# Preprint/arXiv rows are keyed by ident, not PMID.
CITE_ONLY_IDENT = {
    "arXiv:2507.03209": "Built on 1,100 unique ILs taken from a SINGLE source - Xu et al.'s HeLa "
                        "data, which is already your YX_2024/AGILE. So it is not new data by "
                        "definition; it is a foil. Claims to 'support assessment of model "
                        "generalization' while evaluating inside one lab's data.",
    "arXiv:2609.17721": "476 formulations compiled from 81 studies, IVIS biodistribution rather "
                        "than transfection, ROC-AUC 0.839-0.874 with no cross-study or "
                        "leave-one-study-out validation. The live failure mode, dated Sep 2026.",
}


def load_access():
    """Merge the automated SI probe with hand-verified overrides."""
    here = os.path.dirname(os.path.abspath(__file__))
    si = {}
    p = os.path.join(here, "vet_si.json")
    if os.path.exists(p):
        si = json.load(open(p))
    out = {}
    for pmid, r in si.items():
        out[pmid] = (r.get("grade", "C"), r.get("note", ""))
    for pmid, (g, n) in ACCESS_OVERRIDE.items():
        out[pmid] = (g, n)
    # the unauthenticated probe is the most recent evidence, so it wins
    for pmid, (g, n) in ACCESS_PROBE.items():
        out[pmid] = (g, n)
    return out


def main():
    src = os.path.join(SCRATCH, "recs3.json")
    have_p = os.path.join(SCRATCH, "have.json")
    if not os.path.exists(src):
        src, have_p = ARCHIVE_RECS, ARCHIVE_HAVE
    if not os.path.exists(src):
        sys.exit(f"missing {src}")
    auto = {r["pmid"]: r for r in json.load(open(src))}
    have = json.load(open(have_p))
    extra_meta_p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "meta_batch2.json")
    if os.path.exists(extra_meta_p):
        for pmid, m in json.load(open(extra_meta_p)).items():
            auto.setdefault(pmid, {}).update(m)
    access = load_access()
    M.update(M2)

    rows = []
    for pmid, man in M.items():
        a = auto.get(pmid, {})
        expid = man.get("expid") or ",".join(have.get(pmid, [])) or ""
        status = man["status"]
        if expid and status == "new":
            status = "have"  # safety net: never advertise something already curated
        rows.append(dict(
            pmid=pmid, ident=f"PMID {pmid}",
            title=a.get("title", ""), journal=a.get("journal", ""), year=a.get("year", ""),
            first=a.get("first", ""), last=a.get("last", ""), doi=a.get("doi", ""),
            url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            lib=man.get("lib") if man.get("lib") is not None else (a.get("maxlib") or None),
            lab=man["lab"], readout=man["readout"], kind=man["kind"],
            status=status, prio=man["prio"], expid=expid, note=man["note"],
            grade=access.get(pmid, ("C", "not assessed"))[0],
            grade_note=access.get(pmid, ("C", "not assessed"))[1],
            abstract=(a.get("abstract", "") or "")[:1500],
        ))

    # every corpus PMID that is not already a manual row -> listed as "have"
    for pmid, eids in sorted(have.items()):
        if pmid in M:
            continue
        a = auto.get(pmid, {})
        rows.append(dict(
            pmid=pmid, ident=f"PMID {pmid}", title=a.get("title", "(not retrieved in sweep)"),
            journal=a.get("journal", ""), year=a.get("year", ""), first=a.get("first", ""),
            last=a.get("last", ""), doi=a.get("doi", ""),
            url=f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            lib=a.get("maxlib") or None, lab="", readout="", kind="ivl", status="have",
            prio=0, expid=",".join(eids), grade="-", grade_note="",
            note="Already curated in your corpus.", abstract=(a.get("abstract", "") or "")[:1500]))

    for e in EXTRA:
        e.setdefault("grade", "A" if "arxiv" in (e.get("url") or "") else "-")
        e.setdefault("grade_note", "open preprint, full text free" if e["grade"] == "A" else "")
    rows.extend(EXTRA)

    # ---- inclusion triage -----------------------------------------------------
    # 1. confirmed unobtainable -> gone entirely
    kept = []
    for r in rows:
        if r["grade"] == "D" or r["pmid"] in REJECTED:
            print(f"DROP unobtainable  {r['ident']}: {r['title'][:58]}")
        else:
            kept.append(r)
    rows = kept

    # 2. fails the nucleic-acid-transfection / meaningful-addition tests -> gone,
    #    with the reason preserved in STRUCK so a later sweep re-strikes it
    kept, n_struck = [], 0
    for r in rows:
        why = STRUCK.get(r["pmid"] or "")
        if why and r["status"] == "new":
            print(f"STRIKE {r['ident']}: {why[:70]}")
            n_struck += 1
        else:
            kept.append(r)
    rows = kept

    # 3. not addable data, but belongs in the manuscript -> status 'cite'
    n_cite = 0
    for r in rows:
        why = CITE_ONLY.get(r["pmid"] or "") or CITE_ONLY_IDENT.get(r["ident"])
        if why and r["status"] == "new":
            r["status"], r["prio"] = "cite", 0
            r["note"] = why
            n_cite += 1

    # 4. label what varies and what is measured, rather than excluding on either
    for r in rows:
        v, m = TAGS.get(r["pmid"] or "", (None, None))
        r["varies"] = v or "IL structure"
        r["measures"] = m or "expression"

    # 5. overwrite the regex-derived library size with the re-verified one
    for r in rows:
        if r["pmid"] in SIZES:
            n, counts, how = SIZES[r["pmid"]]
            r["lib"], r["lib_counts"], r["lib_src"] = n, counts, how
        else:
            r.setdefault("lib_counts", "")
            r.setdefault("lib_src", None)

    # 6. drop candidates whose VERIFIED size is under 100 unique LNPs. Together these
    #    were 15 papers holding 867 LNPs - under 8% of the volume - so the per-paper
    #    curation cost is not repaid. Rows with an UNKNOWN size are kept: not stated is
    #    not the same as small, and several of them are large libraries whose abstracts
    #    simply never give a total.
    MIN_LNPS = 100
    # Exceptions kept as citations: too small to train on, but each is a published
    # cross-species concordance rung, which is the paper's own argument.
    KEEP_AS_CITE = {
        "40571690": "45 chemistries dosed in BOTH mice and NHPs, with expression ranking "
                    "differently between species. Too small to train on - cite it as a "
                    "cross-species analogue of Result 3.",
        "35132167": "89 nanoparticles across mice and humanised-liver/NHP models in one lab. "
                    "Same role: a ready-made cross-species concordance rung.",
    }
    kept, n_small = [], 0
    for r in rows:
        if r["status"] == "new" and r["lib"] and r["lib"] < MIN_LNPS:
            if r["pmid"] in KEEP_AS_CITE:
                r["status"], r["prio"] = "cite", 0
                r["note"] = KEEP_AS_CITE[r["pmid"]]
                kept.append(r)
                continue
            print(f"DROP <{MIN_LNPS} LNPs  {r['pmid']}  n={r['lib']:>4}  {r['lab'][:24]:24} "
                  f"{r['title'][:44]}")
            n_small += 1
        else:
            kept.append(r)
    rows = kept
    print(f"dropped {n_small} candidates under {MIN_LNPS} unique LNPs "
          f"({len(KEEP_AS_CITE)} of them kept as citations instead)")

    rows.sort(key=lambda r: (r["prio"] == 0, r["prio"], -(r["lib"] or 0)))
    from collections import Counter
    # mark what we have since curated and merged
    for r in rows:
        a = ADDED.get(r["pmid"])
        if a:
            r["status"], r["prio"] = "added", 0
            r["expid"] = a["expid"]
            r["added_on"] = a["added"]
            r["added_rows_vitro"] = a["rows_vitro"]
            r["added_rows_vivo"] = a["rows_vivo"]
            r["added_folder"] = a["folder"]
            r["added_what"] = a["what"]
            r["added_validation"] = a["validation"]
            if a["rows_vitro"] or a["rows_vivo"]:
                r["lib"] = a["rows_vitro"] or a["rows_vivo"]

    nv = [r for r in rows if r["status"] == "new"]
    print(f"\ntriage: {n_struck} struck, {n_cite} reclassified cite-only")
    print("varies:  ", dict(Counter(r["varies"] for r in nv)))
    print("measures:", dict(Counter(r["measures"] for r in nv)))

    meta = dict(
        generated="2026-10-25",
        swept=1996,
        queries=20,
        cutoff="LNPDB curation cutoff: August 2025",
        corpus_pmids=len(have),
        grade_legend=GRADE_NOTE,
        struck=n_struck,
        rejected=len(REJECTED),
    )
    json.dump(dict(meta=meta, rows=rows), open(OUT, "w"), indent=1)
    print(f"wrote {OUT}: {len(rows)} rows")
    from collections import Counter
    print(Counter((r["status"], r["prio"]) for r in rows))


if __name__ == "__main__":
    main()
