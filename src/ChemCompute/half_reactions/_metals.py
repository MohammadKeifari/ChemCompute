"""Simple metal ion / metal and adjacent-oxidation-state couples vs SHE (25 °C).

Metal solids are written ``M.s`` (not library compounds). Aqueous ions reuse
library species via ``.token``.
"""

from ..compounds import (
    ag_plus,
    al3,
    cd2,
    co_2plus,
    co_3plus,
    cr2,
    cr3,
    cu2,
    cu_plus,
    fe2,
    fe3,
    h_plus,
    hg2,
    mn2,
    ni2,
    pb2,
    sn2,
    sn4,
    v2,
    v3,
    vo2,
    water,
    zn2,
)
from ._builders import hr, library


# Fe+3 + e- ⇌ Fe+2.
@library
def iron_iii():
    return hr(
        f"{fe3().token} & @e = {fe2().token}",
        E0=0.771,
        name="Fe+3/Fe+2",
    )


# Fe+2 + 2 e- ⇌ Fe(s).
@library
def iron():
    return hr(
        f"{fe2().token} & 2_@e = Fe.s",
        E0=-0.447,
        name="Fe+2/Fe",
    )


# Cu+2 + 2 e- ⇌ Cu(s).
@library
def copper():
    return hr(
        f"{cu2().token} & 2_@e = Cu.s",
        E0=0.342,
        name="Cu+2/Cu",
    )


# Cu+ + e- ⇌ Cu(s).
@library
def copper_i():
    return hr(
        f"{cu_plus().token} & @e = Cu.s",
        E0=0.521,
        name="Cu+/Cu",
    )


# Cu+2 + e- ⇌ Cu+.
@library
def copper_ii_i():
    return hr(
        f"{cu2().token} & @e = {cu_plus().token}",
        E0=0.153,
        name="Cu+2/Cu+",
    )


# Ag+ + e- ⇌ Ag(s).
@library
def silver():
    return hr(
        f"{ag_plus().token} & @e = Ag.s",
        E0=0.799,
        name="Ag+/Ag",
    )


# Zn+2 + 2 e- ⇌ Zn(s).
@library
def zinc():
    return hr(
        f"{zn2().token} & 2_@e = Zn.s",
        E0=-0.762,
        name="Zn+2/Zn",
    )


# Ni+2 + 2 e- ⇌ Ni(s).
@library
def nickel():
    return hr(
        f"{ni2().token} & 2_@e = Ni.s",
        E0=-0.257,
        name="Ni+2/Ni",
    )


# Pb+2 + 2 e- ⇌ Pb(s).
@library
def lead():
    return hr(
        f"{pb2().token} & 2_@e = Pb.s",
        E0=-0.126,
        name="Pb+2/Pb",
    )


# Sn+2 + 2 e- ⇌ Sn(s).
@library
def tin():
    return hr(
        f"{sn2().token} & 2_@e = Sn.s",
        E0=-0.138,
        name="Sn+2/Sn",
    )


# Sn+4 + 2 e- ⇌ Sn+2.
@library
def tin_iv():
    return hr(
        f"{sn4().token} & 2_@e = {sn2().token}",
        E0=0.151,
        name="Sn+4/Sn+2",
    )


# Cd+2 + 2 e- ⇌ Cd(s).
@library
def cadmium():
    return hr(
        f"{cd2().token} & 2_@e = Cd.s",
        E0=-0.403,
        name="Cd+2/Cd",
    )


# Co+2 + 2 e- ⇌ Co(s).
@library
def cobalt():
    return hr(
        f"{co_2plus().token} & 2_@e = Co.s",
        E0=-0.280,
        name="Co+2/Co",
    )


# Co+3 + e- ⇌ Co+2.
@library
def cobalt_iii():
    return hr(
        f"{co_3plus().token} & @e = {co_2plus().token}",
        E0=1.82,
        name="Co+3/Co+2",
    )


# Al+3 + 3 e- ⇌ Al(s).
@library
def aluminum():
    return hr(
        f"{al3().token} & 3_@e = Al.s",
        E0=-1.662,
        name="Al+3/Al",
    )


# Cr+3 + 3 e- ⇌ Cr(s).
@library
def chromium():
    return hr(
        f"{cr3().token} & 3_@e = Cr.s",
        E0=-0.744,
        name="Cr+3/Cr",
    )


# Cr+3 + e- ⇌ Cr+2.
@library
def chromium_ii():
    return hr(
        f"{cr3().token} & @e = {cr2().token}",
        E0=-0.407,
        name="Cr+3/Cr+2",
    )


# Mn+2 + 2 e- ⇌ Mn(s).
@library
def manganese():
    return hr(
        f"{mn2().token} & 2_@e = Mn.s",
        E0=-1.185,
        name="Mn+2/Mn",
    )


# Hg+2 + 2 e- ⇌ Hg(l).
@library
def mercury():
    return hr(
        f"{hg2().token} & 2_@e = Hg.l",
        E0=0.851,
        name="Hg+2/Hg",
    )


# VO+2 + 2 H+ + e- ⇌ V+3 + H2O(l).
@library
def vanadyl():
    return hr(
        f"{vo2().token} & 2_{h_plus().token} & @e = {v3().token} & {water().token}",
        E0=0.337,
        name="VO+2/V+3",
    )


# V+3 + e- ⇌ V+2.
@library
def vanadium_iii():
    return hr(
        f"{v3().token} & @e = {v2().token}",
        E0=-0.255,
        name="V+3/V+2",
    )


__all__ = [
    "aluminum",
    "cadmium",
    "chromium",
    "chromium_ii",
    "cobalt",
    "cobalt_iii",
    "copper",
    "copper_i",
    "copper_ii_i",
    "iron",
    "iron_iii",
    "lead",
    "manganese",
    "mercury",
    "nickel",
    "silver",
    "tin",
    "tin_iv",
    "vanadium_iii",
    "vanadyl",
    "zinc",
]
