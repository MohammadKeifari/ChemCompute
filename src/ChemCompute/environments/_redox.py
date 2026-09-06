"""Coupled half-reaction environments (library couples used as-is)."""

from .. import half_reactions, reactions
from ._builders import assemble


def water_limits(concentrations=None, *, T=298, volume=1.0):
    """O2/H2O and H+/H2 water window, with Kw."""
    return assemble(
        reactions.water_kw(),
        half_reactions.hydrogen(),
        half_reactions.oxygen(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def iron_couple(concentrations=None, *, T=298, volume=1.0):
    """Fe+3/Fe+2, with Kw."""
    return assemble(
        reactions.water_kw(),
        half_reactions.iron_iii(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def daniel_cell(concentrations=None, *, T=298, volume=1.0):
    """Cu+2/Cu and Zn+2/Zn (two-electrode shared Eh)."""
    return assemble(
        half_reactions.copper(),
        half_reactions.zinc(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


__all__ = [
    "daniel_cell",
    "iron_couple",
    "water_limits",
]
