# DECISIONS — curation/xuq_24574196

Wang, Alberti, Varone, Pouli, Georgakoudi, Xu. *Enhanced Intracellular siRNA
Delivery Using Bioreducible Lipid-Like Nanoparticles.* Adv. Healthcare Mater.
2014. PMID 24574196, DOI 10.1002/adhm.201400039.

`build_dataset.py` refuses a folder with an `OPEN` entry. Every entry below is
deliberately left open: **this folder emits zero rows and is not merge-ready**,
and the gate should keep saying so until a human decides otherwise.

### D-1
Amine 2 is never drawn, named chemically, or given a formula. SI p.2 says only
"amine 1 or 2 were mixed with different acrylates". Every lipidoid built on it —
an unknown fraction of the paper's library — is therefore unenumerable from this
SI. The Xu lab's amine catalogue numbers recur across its papers, so amine 2 may
be recoverable from a sibling paper or from the main text's Figure 1, but it is
not recoverable from anything in hand. Guessing it would invent library members.
STATUS: OPEN

### D-2
`O14B` and `O18B` are named on SI p.2 and given no formula, no NMR and no
drawing. Their S-alkyls are extrapolated from the single pinned B-member
(`O16B`, dodecyl) by the chain-heavy-atom convention, which is a one-point fit
corroborated only by its consistency with `O16` and with the paper's own
"similar chemical structures" claim. They are emitted with
`confidence = medium`, `structure_assumed = 1` and a note. Whether a
two-member extrapolation belongs in the corpus at all is a decision for a human.
The primary source for these acrylates is Ryu, Roy, Ventura & Thayumanavan,
*Langmuir* 2010, **26**, 7086, which was NOT acquired and would settle it.
STATUS: OPEN

### D-3
The main article is closed access: no PMC deposit (NCBI `elink` returns
`pubmed_pmc_refs` only, no self-link), not in Europe PMC, `oa_status: closed`
in OpenAlex with no repository copy, and a Cloudflare interstitial on the Wiley
page to a plain client. Browser acquisition was outside this run's scope and no
bypass was attempted. **Every corpus-admissible value in this paper is in that
article.** Whether to acquire it through the institutional route, or to leave
this folder as a structures-only deposit, is a human decision.
STATUS: OPEN

### D-4
The N/P ratio of the Fig S3 viability experiment is stated as 2:1 in the
Methods and 5:1 in the caption (contradiction X-01). It does not affect anything
emitted — the panel is inadmissible — but if the main article is acquired, the
same document's unduplicated numbers inherit the doubt this raises.
STATUS: OPEN

### D-5
This SI draws no chemical structure anywhere, so the C10 anchor rests on the
printed ¹H NMR **assignment list** rather than on a drawn structure. That is
weaker ground truth than C10 normally has: an assignment list states
connectivity in prose that the authors wrote, where a drawing is something they
had to render. It is recorded as such in `reference_structures.csv` rather than
presented as a figure read. Whether that clears C10's bar is a human decision.
STATUS: OPEN
