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
