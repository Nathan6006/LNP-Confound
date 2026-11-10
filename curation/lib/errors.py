#!/usr/bin/env python3
"""The exception types curation/lib raises.

Every one of these exists so a detector can REFUSE instead of returning a
plausible number. That posture is the whole point of the package: error class 28
is a detector miss falling through to a hard-coded literal, and F-033 is what it
costs -- an origin left as three literals after the fit that set it was removed,
biasing 48 of 48 measurements one way.

A caller that wants a default must write the default itself, at the call site,
where a reviewer can see it.
"""


class CurationLibError(Exception):
    """Base. Catch this to catch any refusal from curation/lib."""


class CalibrationError(CurationLibError):
    """An axis could not be calibrated from its own furniture.

    Raised rather than returning a fallback scale. Error classes 20 and 27.
    """


class DetectionError(CurationLibError):
    """A feature the caller asked for was not found, or was found ambiguously.

    Raised rather than returning the best candidate. Error classes 28, 32, 35.
    """


class ValidationTargetError(CurationLibError):
    """A function was handed the quantity it exists to validate.

    Error class 27: `digitize_fig2fgh.py` fitted its pixels-per-percent by least
    squares onto the very hit rates the panel was there to validate, so a uniform
    error in the heatmap read was absorbed and invisible, and all 48 printed
    statistics passed on values that were ~1.1% short.

    The signature to grep for is a fit whose x-vector is the expected answer.
    """


class StructureError(CurationLibError):
    """A built structure failed its own post-build assertion.

    Error classes 14, 15, 31. F-013: `RWMol.RemoveAtom` dropped cholesterol's
    hydroxyl out of C3's neighbour list without correcting C3's chiral tag, and
    537 rows stored epicholesterol. The epimers share molecular formula, exact
    mass to four decimals, the InChIKey skeleton block and all 17 descriptor
    columns, so nothing downstream could see it.
    """


class AmbiguityError(CurationLibError):
    """A name resolves to more than one structure.

    Error classes 30 and 46. Never resolved by row order; see
    curation/corpus_components.py.
    """


class EmitError(CurationLibError):
    """A row was emitted with a flag that is not traceable to its own record.

    Error classes 29 and 34. F-011: a hard-coded "high" relabelled the 24 lipids
    whose enumeration record reads "ASSUMED: no ring N and the paper draws no
    structure" as the project's highest confidence, inverting the labelling so
    the most uncertain rows carried the strongest claim.
    """
