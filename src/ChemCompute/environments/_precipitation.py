"""Precipitation environments: Ksp plus, where it matters, soluble complexes.

``silver_chloride`` couples AgCl(s) dissolution to AgCl2- formation (excess chloride).
``silver_chloride_ammonia`` is AgCl(s) dissolving as Ag(NH3)2+.
Hydroxide precipitates are coupled to ammine complexes the same way.
"""

from .. import reactions
from ._builders import assemble


def silver_chloride(concentrations=None, *, T=298, volume=1.0):
    """AgCl(s) ⇌ Ag+ + Cl- and Ag+ + 2 Cl- ⇌ AgCl2-, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.agcl_ksp(),
        reactions.agcl2_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def silver_chloride_ammonia(concentrations=None, *, T=298, volume=1.0):
    """AgCl(s) with Ag(NH3)2+ and NH4+/NH3, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.ammonium(),
        reactions.agcl_ksp(),
        reactions.ag_nh3_2_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def silver_bromide(concentrations=None, *, T=298, volume=1.0):
    """AgBr(s) ⇌ Ag+ + Br-, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.agbr_ksp(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def silver_iodide(concentrations=None, *, T=298, volume=1.0):
    """AgI(s) ⇌ Ag+ + I-, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.agi_ksp(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def calcium_fluoride(concentrations=None, *, T=298, volume=1.0):
    """CaF2(s) ⇌ Ca+2 + 2 F-, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.caf2_ksp(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def barium_sulfate(concentrations=None, *, T=298, volume=1.0):
    """BaSO4(s) ⇌ Ba+2 + SO4-2, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.baso4_ksp(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def calcium_carbonate(concentrations=None, *, T=298, volume=1.0):
    """CaCO3(s) coupled to H2CO3 / HCO3- / CO3-2 and Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.carbonic_acid_1(),
        reactions.carbonic_acid_2(),
        reactions.caco3_ksp(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def iron_hydroxide(concentrations=None, *, T=298, volume=1.0):
    """Fe(OH)3(s) ⇌ Fe+3 + 3 OH-, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.feoh3_ksp(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def copper_hydroxide_ammine(concentrations=None, *, T=298, volume=1.0):
    """Cu(OH)2(s) with Cu(NH3)4+2 and NH4+/NH3, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.ammonium(),
        reactions.cuoh2_ksp(),
        reactions.cu_nh3_4_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def zinc_hydroxide_ammine(concentrations=None, *, T=298, volume=1.0):
    """Zn(OH)2(s) with Zn(NH3)4+2 and NH4+/NH3, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.ammonium(),
        reactions.znoh2_ksp(),
        reactions.zn_nh3_4_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def iron_thiocyanate_hydroxide(concentrations=None, *, T=298, volume=1.0):
    """Fe(OH)3(s) competing with FeSCN+2, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.feoh3_ksp(),
        reactions.fescn_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


__all__ = [
    "barium_sulfate",
    "calcium_carbonate",
    "calcium_fluoride",
    "copper_hydroxide_ammine",
    "iron_hydroxide",
    "iron_thiocyanate_hydroxide",
    "silver_bromide",
    "silver_chloride",
    "silver_chloride_ammonia",
    "silver_iodide",
    "zinc_hydroxide_ammine",
]
