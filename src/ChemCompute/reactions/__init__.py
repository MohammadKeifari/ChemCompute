"""Library of common aqueous reactions: acid–base, precipitation, complexes, redox.

Each name is a function that **returns** the shared :class:`~ChemCompute.Reaction`
(the same object every time), with ``K`` or ``infinite_K`` set and every
concentration 0. Put amounts on the environment::

    from ChemCompute import Enviroment
    from ChemCompute.reactions import water_kw, acetic_acid, fescn_kf

    env = Enviroment(
        water_kw(),
        acetic_acid(),
        concentrations={"CH3COOH": 0.10},
    )

Library compound functions are reused via ``.token`` when the species already
exists with the right phase. Molecular acids stored as a gas/liquid/solid in
``compounds`` are written ``HF.aq`` so they remain in Q.
"""

from __future__ import annotations

from .._general import Reaction
from ._acid_base import *
from ._acid_base import __all__ as _ACID_BASE
from ._complexes import *
from ._complexes import __all__ as _COMPLEXES
from ._precipitation import *
from ._precipitation import __all__ as _PRECIPITATION
from ._redox import *
from ._redox import __all__ as _REDOX

__all__ = [
    "all_reactions",
    *_ACID_BASE,
    *_PRECIPITATION,
    *_COMPLEXES,
    *_REDOX,
]


def all_reactions():
    """Unique library reactions (one object per name)."""
    seen = {}
    for name in __all__:
        obj = globals()[name]
        if callable(obj) and name != "all_reactions":
            obj = obj()
        if isinstance(obj, Reaction) and id(obj) not in seen:
            seen[id(obj)] = obj
    return tuple(seen.values())
