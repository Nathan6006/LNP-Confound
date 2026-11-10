#!/usr/bin/env python3
"""One emit seed, in its own process. Prints the refusal and exits non-zero."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "lib"))
import emit
from errors import EmitError

mode = sys.argv[1]
items = ["LC_2024|1|RLU_B16", "LC_2024|181|RLU_B16"]
records = {items[0]: {"hl": "18PG"}, items[1]: {"hl": "DSPC"}}
conf = (lambda i, r: "medium" if r["hl"] == "18PG" else "high")
note = (lambda i, r: "18PG is emitted as dioleoyl PG; see DECISIONS.md D-2.")

if mode == "literal_confidence":
    conf = "high"
elif mode == "caveat_without_note":
    note = lambda i, r: ""
elif mode == "bad_vocabulary":
    conf = lambda i, r: "very high"
elif mode == "missing_record":
    records = {items[1]: {"hl": "DSPC"}}

try:
    rows = emit.emit_rows(items, records, confidence=conf, note=note,
                          positional_ambiguity=lambda i, r: "0",
                          structure_assumed=lambda i, r: "0",
                          value_is_bound=lambda i, r: "",
                          value_normalization_unit=lambda i, r: "experiment",
                          label="seed:" + mode)
except EmitError as e:
    print("REFUSED:", str(e).splitlines()[0])
    sys.exit(1)
print("emitted %d rows, caveats %d" % (len(rows), emit.caveat_count(rows)))
sys.exit(0)
