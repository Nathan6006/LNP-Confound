05_values.csv IS INTENTIONALLY EMPTY -- header, zero data rows.

Not "we could not digitise it well enough". There is nothing in the material
held to digitise. Every transfection measurement in this paper is in the main
article (Figures 2-4), which this folder does not have; see 00_acquisition.json
"NOT_held" for the three acquisition routes tried and what each returned, and
00_value_sources.json for the per-panel tier assessment.

The claim is checkable in one command:

    python3 -I stage00_scan_for_values.py

It scans all 49 SI pages for transfection-readout vocabulary and for the
library's own member labels (AxByCz), and reports any page carrying both. One
page does -- S29, the Figure S7 caption -- and the script adjudicates it inline:
that figure is Cy5-mRNA CELLULAR UPTAKE by flow cytometry for a single lipid,
prints no number, and is not a transfection readout. An unadjudicated candidate
page makes the script exit non-zero. Run it with --seed to watch it fire on a
synthetic screening page.

CONSEQUENCE FOR THE CORPUS: zero rows are emittable from this folder today.
The nine structures in 03_enumerated.csv are correct and fully checked, and
07_rows_quality.csv already carries their quality columns, so the folder is one
value source away from being mergeable. It is not mergeable now, and
DECISIONS.md D-2 is open on exactly that.
