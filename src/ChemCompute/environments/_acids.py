"""Coupled polyprotic acid–base environments (all Ka steps plus Kw)."""

from .. import reactions
from ._builders import assemble


def phosphoric_acid(concentrations=None, *, T=298, volume=1.0):
    """H3PO4 ⇌ H2PO4- ⇌ HPO4-2 ⇌ PO4-3, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.phosphoric_acid_1(),
        reactions.phosphoric_acid_2(),
        reactions.phosphoric_acid_3(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def carbonic_acid(concentrations=None, *, T=298, volume=1.0):
    """H2CO3 ⇌ HCO3- ⇌ CO3-2, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.carbonic_acid_1(),
        reactions.carbonic_acid_2(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def hydrogen_sulfide(concentrations=None, *, T=298, volume=1.0):
    """H2S ⇌ HS- ⇌ S-2, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.hydrogen_sulfide_1(),
        reactions.hydrogen_sulfide_2(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def sulfurous_acid(concentrations=None, *, T=298, volume=1.0):
    """H2SO3 ⇌ HSO3- ⇌ SO3-2, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.sulfurous_acid_1(),
        reactions.sulfurous_acid_2(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def sulfuric_acid(concentrations=None, *, T=298, volume=1.0):
    """H2SO4 → HSO4- ⇌ SO4-2, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.sulfuric_acid_1(),
        reactions.sulfuric_acid_2(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def ammonia(concentrations=None, *, T=298, volume=1.0):
    """NH4+ ⇌ NH3 + H+, with Kw (do not also add ammonia Kb)."""
    return assemble(
        reactions.water_kw(),
        reactions.ammonium(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


__all__ = [
    "ammonia",
    "carbonic_acid",
    "hydrogen_sulfide",
    "phosphoric_acid",
    "sulfuric_acid",
    "sulfurous_acid",
]
