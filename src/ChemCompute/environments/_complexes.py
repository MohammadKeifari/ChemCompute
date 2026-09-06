"""Complex-formation environments, including ammine systems coupled to NH4+/NH3."""

from .. import reactions
from ._builders import assemble


def fescn(concentrations=None, *, T=298, volume=1.0):
    """Fe+3 + SCN- ⇌ FeSCN+2, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.fescn_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def ferroin(concentrations=None, *, T=298, volume=1.0):
    """Fe+2 + 3 phen ⇌ Fe(phen)3+2, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.ferroin_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def triiodide(concentrations=None, *, T=298, volume=1.0):
    """I2 + I- ⇌ I3-, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.i3_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def copper_ammine(concentrations=None, *, T=298, volume=1.0):
    """Cu(NH3)4+2 formation coupled to NH4+/NH3 and Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.ammonium(),
        reactions.cu_nh3_4_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def nickel_ammine(concentrations=None, *, T=298, volume=1.0):
    """Ni(NH3)6+2 formation coupled to NH4+/NH3 and Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.ammonium(),
        reactions.ni_nh3_6_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def silver_ammine(concentrations=None, *, T=298, volume=1.0):
    """Ag(NH3)2+ formation coupled to NH4+/NH3 and Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.ammonium(),
        reactions.ag_nh3_2_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def zinc_ammine(concentrations=None, *, T=298, volume=1.0):
    """Zn(NH3)4+2 formation coupled to NH4+/NH3 and Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.ammonium(),
        reactions.zn_nh3_4_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def ferrocyanide(concentrations=None, *, T=298, volume=1.0):
    """Fe+2 + 6 CN- ⇌ Fe(CN)6-4, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.fe_cn6_4_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


def ferricyanide(concentrations=None, *, T=298, volume=1.0):
    """Fe+3 + 6 CN- ⇌ Fe(CN)6-3, with Kw."""
    return assemble(
        reactions.water_kw(),
        reactions.fe_cn6_3_kf(),
        concentrations=concentrations,
        T=T,
        volume=volume,
    )


__all__ = [
    "copper_ammine",
    "ferricyanide",
    "ferrocyanide",
    "ferroin",
    "fescn",
    "nickel_ammine",
    "silver_ammine",
    "triiodide",
    "zinc_ammine",
]
