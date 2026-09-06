"""Library of common compounds with optional melting/boiling points and UV-Vis peaks.

Import named species and interpolate the live object into a reaction string::

    from ChemCompute.compounds import water, h_plus, oh_minus, mno4
    Reaction.from_string(f"{water().token} > {h_plus().token} & {oh_minus().token}")

Each name is a function that **returns** the shared library compound (the same
object every time). Copy before mutating ``spectrum``, ``mp``, or ``bp``.

``mp`` and ``bp`` are 1 atm values in **kelvin** (the same scale as ``env.T``).
Ions and species that sublime or decompose omit them.

UV-Vis points use wavelength in metres and ε in M^-1 m^-1 (literature
M^-1 cm^-1 × 100). With the default ``uvvis_spectrum(..., path_length=0.01)``
that is a 1 cm cuvette. Spectra are piecewise-linear envelopes: several
(wavelength, ε) points are connected, with ε → 0 on the wings, so a scan
across wavelengths looks like a UV-Vis trace rather than a single spike.
Peaks are omitted when a clear aqueous (or neat) envelope was not available.

Look up a species by formula with ``get("H2O")``.
"""

from __future__ import annotations

from .._general import Compound
from ._complexes import *
from ._complexes import __all__ as _COMPLEXES
from ._ions import *
from ._ions import __all__ as _IONS
from ._molecules import *
from ._molecules import __all__ as _MOLECULES
from ._solids import *
from ._solids import __all__ as _SOLIDS

__all__ = [
    "all_compounds",
    "get",
    *_MOLECULES,
    *_IONS,
    *_COMPLEXES,
    *_SOLIDS,
]


def all_compounds():
    """Unique library compounds (aliases such as ``water`` / ``h2o`` collapse)."""
    seen = {}
    for name in __all__:
        obj = globals()[name]
        if callable(obj) and name not in ("all_compounds", "get"):
            obj = obj()
        if isinstance(obj, Compound) and id(obj) not in seen:
            seen[id(obj)] = obj
    return tuple(seen.values())


def get(formula: str) -> Compound:
    """Return the library compound with this formula.

    Raises
    ------
    KeyError
        If the formula is missing or used by more than one distinct object.
    """
    matches = [c for c in all_compounds() if c.formula == formula]
    if not matches:
        raise KeyError(f"no library compound with formula {formula!r}")
    first = matches[0]
    if any(item is not first for item in matches):
        raise KeyError(
            f"formula {formula!r} is used by more than one library compound"
        )
    return first
