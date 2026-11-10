# Source documents — what is NOT in this repository

Every paper and supplement this corpus is built from is excluded by
`.gitignore` (`*.pdf`, `*.xlsx`, `*.docx`) because none of it is ours to
redistribute. This file is the index: for each curation folder it gives the
paper, its DOI, and the **exact filenames** the scripts expect, so a fresh
clone can be made runnable by downloading each one from the publisher and
dropping it in at the path shown.

Page renders (`curation/*/00_pages/`) and vendored upstream clones
(`curation/*/repo/`) are also excluded and are **regenerated**, not downloaded.

## Curated and merged into `data/`

| folder | `Experiment_ID` | paper | DOI | files expected |
|---|---|---|---|---|
| `lin_aromatic_2026` | `LL_2026`, `LL_2026_2` | Lin J-H, …, Liu Q. *Aromatic Ionizable Lipids Enhance mRNA Delivery via Coordinating Multiple Mechanisms.* JACS 2026, 148, 6960–6969 | `10.1021/jacs.5c16971` | `ja5c16971.pdf`, `ja5c16971_si_001.pdf` |
| `ql_2026` | `QL_2026` | Liu Q, …, Li W. *Synthesis and Screening of Bispecific Linker-based Combinatorial Libraries Identified Thiosuccinimide Lipids for mRNA Delivery.* Adv Healthc Mater 2025 | `10.1002/adhm.202502225` | `Adv Healthcare Materials - 2025 - Liu - Synthesis and Screening of Bispecific Linker‐based Combinatorial Libraries.pdf` |
| `siegwart_42634137` | `ZT_2026` | Tian Z, …, Siegwart DJ. *Esterified Cholesterol Conjugates Enable LNP-Mediated mRNA Delivery to the Blood–Brain Barrier Following Intravenous Administration.* Adv Mater 2026 | `10.1002/adma.74762` | `Advanced Materials - 2026 - Tian - Esterified Cholesterol Conjugates Enable LNP‐Mediated mRNA Delivery to the Blood Brain.pdf`, `adma74762-sup-0001-suppmat.pdf` |
| `xd_2026_lipidai` | `XD_2026` | Duan X, …, Song Y. *A novel machine learning framework-based rapid screening of ionizable lipids in LNPs for highly-efficient mRNA expression.* Acta Pharm Sin B 2026 | `10.1016/j.apsb.2025.10.040` | `1-s2.0-S2211383525007221-main.pdf`; `repo/` = `github.com/yuhanchenzz/LipidAI` |
| `wht_2026_falcon` | `WT_2026` | Toh W-H, Cheng L, …, Mao H-Q. *FALCON: Closed-Loop Multi-Objective Optimization of Lipid Nanoparticles for Cell-Selective mRNA Delivery.* bioRxiv 2026 (preprint, no PMID) | `10.64898/2026.06.13.731774` | `2026.06.13.731774v1.full.pdf`, `media-1.pdf` (SI, 62 pp); `repo/` = `github.com/MaoResearchGroup/falcon-lnp-optimization` |
| `ng_2025` | — (in progress) | Gong N, Kim D, …, Mitchell MJ. *Mannich reaction-based combinatorial libraries identify antioxidant ionizable lipids for mRNA delivery with reduced immunogenicity.* Nat Biomed Eng 2025, 9, 2181–2195 | `10.1038/s41551-025-01422-8` | `s41551-025-01422-8.pdf`, `41551_2025_1422_MOESM1_ESM.pdf` |

## Downloaded, not yet curated

| folder | paper | DOI | files expected |
|---|---|---|---|
| `amino_lipid_fg_2026` | Gong F, …, Li W. *Amino acid-derived ionizable lipids enable inhaled base editing for therapeutic gene correction in the lung.* Nat Mater 2026 | `10.1038/s41563-026-02555-0` | `s41563-026-02555-0.pdf`, `41563_2026_2555_MOESM1_ESM.pdf` |
| `jl_2026` | Lin J, Zhu Y, …, Mao H-Q. *Imidazolium lipid-based nanoparticles enable effective mRNA delivery and cellular immune response.* Materials Today 2026, 95, 103265 | — (Elsevier PII `S1369702126001112`) | `1-s2.0-S1369702126001112-main.pdf`, `…-mmc1.pdf` |
| `liu_iphos` | Liu S, …, Siegwart DJ. *Membrane-destabilizing ionizable phospholipids for organ-selective mRNA delivery and CRISPR–Cas gene editing.* Nat Mater 2021 | `10.1038/s41563-020-00886-0` | `s41563-020-00886-0.pdf`, `41563_2020_886_MOESM1_ESM.pdf` |
| `sd_2026` | Dong S, …, Li W. *Lipid nanoparticles optimized for large RNA cargo and tissue targeting enhance in vivo genome editing.* Nat Biotechnol 2026 | `10.1038/s41587-026-03298-8` | `s41587-026-03298-8.pdf`, `41587_2026_3298_MOESM1_ESM.pdf` (45 MB) |

## Needed but not held — open finding F-021

| what | paper | DOI / PMID |
|---|---|---|
| `ZH_2023`'s sterol conjugate `A4I18-1R2C70` cleaves to a C26H44O 19-nor sterol, not cholesterol. The corpus inherited these 325 in-vivo rows from upstream LNPDB; the source paper has never been held locally, so the correct structure cannot be settled. | He Z, Le Z, …, Chen Y. *A Multidimensional Approach to Modulating Ionizable Lipids for High-Performing and Organ-Selective mRNA Delivery.* Angew Chem Int Ed 2023, 62, e202310401 | `10.1002/anie.202310401` / PMID `37661193` |

## Upstream dataset

`data/Original_LNPDB.csv` and `provenance/superseded/LNPDB_official_download_20261028.csv`
are the pristine upstream LNPDB download (both tracked; they differ, and
`provenance/CHANGELOG.md` records how). Everything else under
`provenance/superseded/` is a pre-rebuild snapshot and is excluded — check out
an earlier commit instead.

## Staged by the 2026-11-10 autonomous run — NOT merged into `data/`

Every folder below still carries OPEN decisions, so `build_dataset.py`
refuses it by design. `rows` is what the curation proposes, not what is in
the corpus. Files are the **exact names on disk**; all are gitignored, so a
fresh clone must re-download each from the publisher. Acquisition routes,
URLs and hashes are in each folder's `PROVENANCE.md` and
`00_value_sources.json`; `provenance/BROWSER_ACQUISITION.md` records which
publishers serve anonymously and which need a browser.

| folder | proposed `Experiment_ID` | rows | PMID | DOI | files expected |
|---|---|---|---|---|---|
| `maohq_39375194` | `LC_2024` | 5112 | 39375194 | 10.1021/acsnano.4c07615 | `nn4c07615_si_001.pdf`, `github/1080_LNP_OHE.csv`, `github/Compiled_validation.xlsx` |
| `siegwart_40644303` | `—` | 505 | 40644303 | 10.1021/acsbiomaterials.5c00463 | `PMC12869450.fulltext.txt`, `PMC12869450.jats.xml`, `ab5c00463_si_001.pdf`, `figures/fig1.jpg`, `figures/fig2.jpg`, `figures/fig3.jpg` … |
| `lij_42130331` | `—` | 300 | 42130331 | 10.1021/acsnano.6c00706 | `nn6c00706_si_001.pdf` |
| `peer_39912611` | `GN_2025` | 144 | 39912611 | 10.1021/acsnano.4c18636 | `PMC11841047.jats.xml`, `nn4c18636_si_001.pdf`, `si_figs/sf01.jpeg`, `si_figs/sf03.jpeg`, `si_figs/sf04.png`, `si_figs/sf05.jpeg` … |
| `linyx_41851281` | `—` | 140 | 41851281 | 10.1038/s41551-026-01640-8 | `41551_2026_1640_MOESM10_ESM.xlsx`, `41551_2026_1640_MOESM11_ESM.xlsx`, `41551_2026_1640_MOESM12_ESM.xlsx`, `41551_2026_1640_MOESM13_ESM.xlsx`, `41551_2026_1640_MOESM14_ESM.xlsx`, `41551_2026_1640_MOESM15_ESM.xlsx` … |
| `andersondg_41203968` | `—` | 132 | 41203968 | 10.1038/s41565-025-02044-6 | `41565_2025_2044_MOESM10_ESM.xlsx`, `41565_2025_2044_MOESM11_ESM.xlsx`, `41565_2025_2044_MOESM12_ESM.xlsx`, `41565_2025_2044_MOESM1_ESM.pdf`, `41565_2025_2044_MOESM2_ESM.pdf`, `41565_2025_2044_MOESM3_ESM.xlsx` … |
| `namjh_39264000` | `—` | 108 | 39264000 | 10.1002/smll.202405618 | `PMC11855220.jats.xml`, `SMLL-21-2405618.pdf`, `smll202405618-sup-0001-SuppMat.docx`, `si_media/figS1_structures.tif`, `si_media/figS2_ID_RLuc.tif`, `si_media/figS3_size_pdi.tif` … |
| `mitchellmj_38982196` | `—` | 100 | 38982196 | 10.1038/s41557-024-01557-2 | `41557_2024_1557_MOESM1_ESM.pdf`, `41557_2024_1557_MOESM2_ESM.pdf`, `41557_2024_1557_MOESM3_ESM.xlsx`, `41557_2024_1557_MOESM4_ESM.xlsx`, `41557_2024_1557_MOESM5_ESM.xlsx`, `41557_2024_1557_MOESM6_ESM.xlsx` … |
| `xuq_32438474` | `—` | 90 | 32438474 | 10.1002/anie.202004994 | `PMC7679290.jats.xml`, `anie202004994-sup-0001-misc_information.pdf`, `pmc/NIHMS1624779-supplement-supp_info.pdf`, `pmc/nihms-1624779-f0001.jpg`, `pmc/nihms-1624779-f0002.jpg`, `pmc/nihms-1624779-f0003.jpg` … |
| `chamolstad_41691524` | `—` | 77 | 41691524 | 10.1021/acsnano.5c19036 | `nn5c19036_si_001.xlsx`, `nn5c19036_si_002.pdf` |
| `zhangz_40433897` | `—` | 75 | 40433897 | 10.1021/acsnano.5c01378 | `nn5c01378_si_001.pdf` |
| `palczewski_41654520` | `—` | 64 | 41654520 | 10.1038/s41467-026-69077-w | `41467_2026_69077_MOESM1_ESM.pdf`, `41467_2026_69077_MOESM2_ESM.pdf`, `41467_2026_69077_MOESM3_ESM.xlsx`, `41467_2026_69077_MOESM4_ESM.pdf`, `41467_2026_69077_MOESM5_ESM.pdf`, `41467_2026_69077_MOESM6_ESM.xlsx` … |
| `whiteheadka_38888972` | `—` | 49 | 38888972 | 10.1002/adhm.202400225 | `PMC11368637.jats.xml`, `adhm202400225-sup-0001-suppmat.docx`, `nihms-2006250-f0001.jpg`, `nihms-2006250-f0002.jpg`, `nihms-2006250-f0003.jpg`, `nihms-2006250-f0004.jpg` … |
| `mitchellmj_41555650` | `EH_2026` | 48 | 41555650 | 10.1021/acsnano.5c19138 | `nn5c19138_si_001.pdf` |
| `dahlmanje_36701517` | `—` | 37 | 36701517 | 10.1021/acs.nanolett.2c04479 | `nl2c04479_si_001.pdf` |
| `andersondg_26469188` | `KK_2015` | 36 | 26469188 | 10.1021/acs.nanolett.5b02497 | `nl5b02497_si_001.pdf` |
| `xuq_32520527` | `—` | 33 | 32520527 | 10.1021/acs.bioconjchem.0c00295 | `PMC9069324.jats.xml`, `bc0c00295_si_001.pdf`, `nihms-1799810-f0001.jpg`, `nihms-1799810-f0002.jpg`, `nihms-1799810-f0003.jpg`, `nihms-1799810-f0004.jpg` … |
| `levkinpa_23909648` | `LL_2013` | 29 | 23909648 | 10.1021/bc400158w | `bc400158w_si_001.pdf` |
| `dixonjb_41269090` | `EM_2025` | 28 | 41269090 | 10.1021/acsnano.5c12080 | `nn5c12080_si_001.pdf`, `pmc_oa/PMC12874650.jats.xml`, `pmc_oa/nn5c12080_0001.jpg`, `pmc_oa/nn5c12080_0002.jpg`, `pmc_oa/nn5c12080_0004.jpg` |
| `sd_2026` | `—` | 18 | 42806111 | 10.1038/s41587-026-03298-8 | `41587_2026_3298_MOESM2_ESM.pdf`, `41587_2026_3298_MOESM3_ESM.pdf`, `41587_2026_3298_MOESM4_ESM.pdf` |
| `khanof_36652629` | `—` | 18 | 36652629 | 10.1021/jacs.2c10670 | `ja2c10670_si_001.pdf` |
| `li_39883839` | `—` | 13 | 39883839 | 10.1073/pnas.2409572122 | `PMC11804478.jats.xml`, `pnas.2409572122.sapp.pdf`, `fig/pnas.2409572122fig01.jpg`, `fig/pnas.2409572122fig02.jpg`, `fig/pnas.2409572122fig03.jpg`, `fig/pnas.2409572122fig04.jpg` … |
| `xuq_32662132` | `—` | 13 | 32662132 | 10.1002/anie.202008082 | `anie202008082-sup-0001-misc_information.pdf`, `figS1_page4_xref95.png`, `figS2_page5_xref101.png`, `figS3_page6_xref106.png`, `figS4_page7_xref111.png`, `figS5_page8_xref116.png` |
| `whiteheadka_30637934` | `KW_2019` | 11 | 30637934 | 10.1002/smll.201805097 | `smll201805097-sup-0001-S1.pdf` |
| `dahlmanje_41183194` | `KG_2025` | 10 | 41183194 | 10.1073/pnas.2409569122 | `pnas.2409569122.sapp.pdf`, `pnas.2409569122fig01.jpg`, `pnas.2409569122fig02.jpg`, `pnas.2409569122fig03.jpg`, `pnas.2409569122fig04.jpg`, `pnas.2409569122fig05.jpg` |
| `mckinlaycj_39105751` | `—` | 5 | 39105751 | 10.1021/acsnano.4c05513 | `PMC11342374.jats.xml`, `nn4c05513_si_001.pdf` |
| `levkinpa_29558113` | `—` | 5 | 29558113 | 10.1021/acs.bioconjchem.8b00007 | `bc8b00007_si_001.pdf` |
| `lius_40500862` | `—` | 4 | 40500862 | 10.1021/jacs.5c01389 | `ja5c01389_si_001.pdf` |
| `goulart_41138890` | `MR_2025` | 3 | 41138890 | 10.1016/j.ijpharm.2025.126278 | `mmc1.docx`, `figs/gr1_lrg.jpg`, `figs/gr3_lrg.jpg`, `figs/gr4_lrg.jpg`, `figs/gr6_lrg.jpg` |
| `songx_39664425` | `—` | 0 | 39664425 | 10.1016/j.apsb.2024.08.032 | `PMC11628822.html`, `europepmc_PMC11628822_fullTextXML.xml`, `europepmc_core_39664425.json`, `gr1.jpg` |
| `pabst_37602460` | `—` | 0 | 37602460 | 10.1021/acsami.3c05481 | `PMC10472336.jats.xml`, `am3c05481_si_001.pdf`, `figs/am3c05481_0002.jpg`, `figs/am3c05481_0003.jpg`, `figs/am3c05481_0004.jpg`, `figs/am3c05481_0005.jpg` … |
| `siegwartdj_42150080` | `—` | 0 | 42150080 | 10.1073/pnas.2528477123 | `PMC13214034.jats.xml`, `pnas.2528477123.sapp.pdf` |
| `liu_41162324` | `—` | 0 | 41162324 | 10.1021/jacs.5c13982 | `ja5c13982_si_001.pdf` |
| `kooh_42690843` | `—` | 0 | 42690843 | 10.1021/jacs.6c13674 | `ja6c13674_si_001.pdf` |
| `peerd_31999380` | `—` | 0 | 31999380 | 10.1002/adma.201906128 | `adma201906128-sup-0001-SuppMat.pdf` |
| `xuq_23651337` | `—` | 0 | 23651337 | 10.1021/sb300023h | `sb300023h_si_001.pdf` |
| `xuq_24574196` | `—` | 0 | 24574196 | 10.1002/adhm.201400039 | `adhm201400039-sup-0001-S1.pdf` |
| `dahlmanje_29489381` | `—` | 0 | 29489381 | 10.1021/acs.nanolett.8b00432 | `PMC6054134.jats.xml`, `nl8b00432_si_001.pdf` |
| `dahlmanje_30216729` | `CS_2018` | 0 | 30216729 | 10.1021/acs.nanolett.8b03149 | `nl8b03149_si_001.pdf` |
| `dahlmanje_30275336` | `—` | 0 | 30275336 | 10.1073/pnas.1811276115 | `PMC6196543.jats.xml`, `pnas.1811276115.sapp.pdf` |
| `zhangn_42569578` | `MY_2026` | 0 | 42569578 | 10.1016/j.mtbio.2026.103483 | `PMC13450591.jats.xml`, `ga1.jpg`, `gr1.jpg`, `gr2.jpg`, `gr3.jpg`, `gr4.jpg` … |
| `mitchellmj_41603124` | `AH_2026` | 0 | 41603124 | 10.1002/adma.202514370 | `PMC12957869.jats.xml`, `adma72049-sup-0001-SuppMat.pdf`, `figs/ADMA-38-e14370-g001.jpg`, `figs/ADMA-38-e14370-g002.jpg`, `figs/ADMA-38-e14370-g003.jpg`, `figs/ADMA-38-e14370-g004.jpg` … |
| `dahlmanje_33758781` | `ZG_2020` | 0 | 33758781 | 10.1002/btm2.10161 | `article_epmc.pdf`, `btm210161-sup-0001-AppendixS1.pdf`, `epmc_suppl/suppl.zip`, `epmc_suppl/extracted/BTM2-5-e10161-g001.jpg`, `epmc_suppl/extracted/BTM2-5-e10161-g002.jpg`, `epmc_suppl/extracted/BTM2-5-e10161-g003.jpg` … |
| `dahlmanje_30394729` | `CS_2018` | 0 | 30394729 | 10.1021/jacs.8b08976 | `ja8b08976_si_001.pdf` |
