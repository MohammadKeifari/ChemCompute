"""Library environments that couple related library reactions and half-reactions.

Each name is a factory that **returns** a new :class:`~ChemCompute.Enviroment`.
The environment uses the library objects returned by those functions (the same
``water()``, ``water_kw()``, ``agcl()``, ``iron_iii()``, …). Concentrations live
on the environment and start at 0.

    from ChemCompute.environments import phosphoric_acid, silver_chloride, water_limits

    env = phosphoric_acid(concentrations={"H3PO4": 0.10})
    agcl = silver_chloride(concentrations={"Cl-": 0.10, "Ag+": 1e-3})
    window = water_limits()
"""

from __future__ import annotations

from .._general import Enviroment
from ._acids import *
from ._acids import __all__ as _ACIDS
from ._complexes import *
from ._complexes import __all__ as _COMPLEXES
from ._precipitation import *
from ._precipitation import __all__ as _PRECIPITATION
from ._redox import *
from ._redox import __all__ as _REDOX

__all__ = [
    "all_environments",
    *_ACIDS,
    *_COMPLEXES,
    *_PRECIPITATION,
    *_REDOX,
]


def all_environments():
    """Return one environment for every library factory."""
    built = []
    for name in __all__:
        if name == "all_environments":
            continue
        factory = globals()[name]
        if callable(factory):
            env = factory()
            if isinstance(env, Enviroment):
                built.append(env)
    return tuple(built)
