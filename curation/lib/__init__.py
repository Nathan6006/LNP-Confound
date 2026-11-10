"""curation/lib -- one audited implementation per hazard class.

Import these modules by path, not as a package, so a per-paper script can be
run from its own folder with `python3 -I`:

    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "lib"))
    import axis, bars, heatmap, chains, assemble, reagents, emit

WHAT BELONGS HERE. A function belongs in this package when the SAME MISTAKE has
been made in two places, or when getting it wrong is invisible downstream. What
does not belong here is a paper's own geometry constants, panel bounds or
library logic -- those stay in the curation folder, where a reviewer reads them
next to the figure.

WHAT EVERY FUNCTION HERE DOES. It refuses. A detector that cannot find its
anchor raises; it never returns a plausible number. Error class 28 is a detector
miss falling through to a literal, and F-033 is what it cost: an origin left as
three hard-coded literals after the fit that set it was removed, biasing 48 of
48 measurements one way. If a caller wants a default it writes the default at
the call site, where a reviewer can see it.

THE RULE THAT GOVERNS CHANGES HERE. CLAUDE.md section 6: "A NEW CHECK IS NOT
DONE UNTIL IT HAS BEEN MADE TO FAIL ON PURPOSE." Every invariant in this package
has a seeded mutation in provenance/test_lib.py, and each was confirmed to fire
before the module was trusted. Three incidents paid for that rule: a re-copied
C10 anchor that could only pass, three ZT_2026 marginal checks that scored a
literal 0 for twelve days while the summary beside them said "ALL FOUR EXACT",
and two versions of C16 that reported "all derived" while examining nothing.
Every one printed a clean line.

MIGRATION. Do not retro-fit a MERGED curation onto this package. Rewriting its
scripts changes nothing in data/ and destroys the audit trail its findings were
adjudicated against. Migrate a merged folder only when a finding sends you into
that file anyway, and then as its own commit.
"""
__all__ = ["axis", "bars", "heatmap", "chains", "assemble", "reagents",
           "emit", "errors"]
