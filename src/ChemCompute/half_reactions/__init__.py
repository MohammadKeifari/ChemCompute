"""Library of common aqueous half-reactions vs SHE at 25 °C.

Each name is a function that **returns** the shared :class:`~ChemCompute.HalfReaction`
(the same object every time). ``E0`` is the standard reduction potential in volts.
Every concentration starts at 0; put amounts on the environment::

    from ChemCompute import Enviroment
    from ChemCompute.half_reactions import hydrogen, oxygen, iron_iii

    env = Enviroment(
        hydrogen(),
        oxygen(),
        iron_iii(),
        concentrations={"Fe+3": 0.01, "Fe+2": 0.001},
        buffer=["H+"],
    )

Library compound functions are reused via ``.token`` when the species already
exists with the right phase. Metal solids are ``M.s``. Liquid water is omitted
from Q. Aqueous peroxide and sulfur dioxide are written ``H2O2.aq`` / ``SO2.aq``
so they stay in Q.
"""

from __future__ import annotations

from .._half_reaction import HalfReaction
from ._complexes import *
from ._complexes import __all__ as _COMPLEXES
from ._halogens import *
from ._halogens import __all__ as _HALOGENS
from ._metals import *
from ._metals import __all__ as _METALS
from ._oxoanions import *
from ._oxoanions import __all__ as _OXOANIONS
from ._water import *
from ._water import __all__ as _WATER

__all__ = [
    "all_half_reactions",
    *_WATER,
    *_HALOGENS,
    *_METALS,
    *_OXOANIONS,
    *_COMPLEXES,
]


def all_half_reactions():
    """Unique library half-reactions (aliases such as ``she`` / ``hydrogen`` collapse)."""
    seen = {}
    for name in __all__:
        obj = globals()[name]
        if callable(obj) and name != "all_half_reactions":
            obj = obj()
        if isinstance(obj, HalfReaction) and id(obj) not in seen:
            seen[id(obj)] = obj
    return tuple(seen.values())
