# Four papers that need a human download

Everything else in the queue was acquired without a login. These four were not,
and the reason differs for each — none of them is simply "needs a subscription".

Save each file **exactly as named** into the folder shown. The filenames are the
publisher's own; the extraction scripts and `curation/SOURCES.md` expect them.
`source/` is gitignored, so nothing you drop in will be committed.

---

## 1. `curation/cheng_41656866/source/`

**Wang Q, Lin Y, Xiao J, Xu K, Luo Z, Ren H, Liu F, Jia L, Wei T, Cheng Q.**
*Optimizing Peptide Ionizable Lipids Enables Efficient and Low-Toxicity mRNA
Delivery for In Vivo Prime Editing and Protein Replacement Therapy.*
**Advanced Materials 38 (2026)** · DOI `10.1002/adma.202522552` · PMID 41656866

- Page: <https://advanced.onlinelibrary.wiley.com/doi/10.1002/adma.202522552>
- Open **Supporting Information**, download **`adma72475-sup-0001-SuppMat.pdf`** (**44.3 MB**)

**No login needed** — the page loads anonymously; I read the title and the SI
link off it. The automation failed only on the 44 MB transfer.

---

## 2. `curation/liu_42585110/source/`

**Jiang K, Shao B, Li R, Liang T, You C, Wang L, Cen J, Tan J, Zhang Y, Chen M,
Liu S.** *Retinoic-Acid-Derived Ionizable Lipids Enable Spleen-Selective Antigen
Expression for mRNA Cancer Vaccination.*
**Angew. Chem. Int. Ed. 65 (2026)** · DOI `10.1002/anie.3823572` · PMID 42585110

- Page: <https://onlinelibrary.wiley.com/doi/10.1002/anie.3823572>
- Open **Supporting Information**, download **`anie74215-sup-0001-SuppMat.docx`** (**13.5 MB**)

**No login needed**, same as above. 13.5 MB also defeated the bridge, so the
size threshold is lower than I assumed — a 4 MB file went through earlier.

---

## 3. `curation/leejy_42610355/source/`

**Sriram V, Phatangare V, Lee JY.** *A double-amine ionizable lipid library
enables highly efficient mRNA delivery to T cells at reduced lipid mass.*
**Biomaterials Science 14 (2026)** · DOI `10.1039/d6bm00568c` · PMID 42610355

- Page: <https://pubs.rsc.org/en/content/articlelanding/2026/bm/d6bm00568c>
- Download **both** the article PDF **and** the ESI (supplementary) PDF.

**This one is Cloudflare, not a paywall** — it is **CC BY**, so it is free to
read; RSC simply refuses the automated browser (403 on the landing page, the
PDF route and the suppdata route, and it did not clear after a wait). I stopped
rather than work around a bot check. **It is also the only one of the four whose
licence passes §8**, so it is the most valuable to get.

---

## 4. `curation/xuq_22148515/source/`

**Sun S, Wang M, Knupp SA, Soto-Feliciano Y, Hu X, Kaplan DL, Langer R,
Anderson DG, Xu Q.** *Combinatorial library of lipidoids for in vitro DNA
delivery.* **Bioconjugate Chemistry 23 (2012)** · DOI `10.1021/bc200572w` ·
PMID 22148515

- Download the **article PDF**. There is no SI to get: the PMC record
  (`PMC3261308`) declares `hasSuppl=N` and its JATS has **0 `<table-wrap>` and
  0 `<supplementary-material>`**, so **all the values are in its 7 figures**,
  which PMC's asset route will not serve.

**Lowest priority, and it may not be obtainable.** ACS is **not covered by the
institutional session in use** — the browser identifies as Phillips Exeter
Academy, measured earlier in this run and recorded in
`BROWSER_ACQUISITION.md`. If it is not reachable, say so and it stays
`pending_acquisition`; that is a fine outcome.

---

## After you drop the files

Nothing else is needed — tell me and I will launch the extractions. Each folder
already exists with an empty `source/`. If a download lands with a browser-added
suffix (`(1)`, `.download`), rename it to the exact name above, because
`00_value_sources.json` records the filename as part of the provenance.
