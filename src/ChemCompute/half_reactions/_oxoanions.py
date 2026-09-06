"""Oxoanion and elemental-nonmetal couples vs SHE (25 °C)."""

from ..compounds import (
    cr2o7,
    cr3,
    h_plus,
    hno2,
    mn2,
    mno4,
    no,
    no3_minus,
    s2o3_2minus,
    so4_2minus,
    water,
)
from ._builders import hr, library


# MnO4- + 8 H+ + 5 e- ⇌ Mn+2 + 4 H2O(l).
@library
def permanganate():
    return hr(
        f"{mno4().token} & 8_{h_plus().token} & 5_@e = {mn2().token} & 4_{water().token}",
        E0=1.507,
        name="MnO4-/Mn+2",
    )


# MnO4- + 4 H+ + 3 e- ⇌ MnO2(s) + 2 H2O(l).
@library
def permanganate_dioxide():
    return hr(
        f"{mno4().token} & 4_{h_plus().token} & 3_@e = MnO2.s & 2_{water().token}",
        E0=1.679,
        name="MnO4-/MnO2",
    )


# MnO2(s) + 4 H+ + 2 e- ⇌ Mn+2 + 2 H2O(l).
@library
def manganese_dioxide():
    return hr(
        f"MnO2.s & 4_{h_plus().token} & 2_@e = {mn2().token} & 2_{water().token}",
        E0=1.224,
        name="MnO2/Mn+2",
    )


# Cr2O7-2 + 14 H+ + 6 e- ⇌ 2 Cr+3 + 7 H2O(l).
@library
def dichromate():
    return hr(
        f"{cr2o7().token} & 14_{h_plus().token} & 6_@e = 2_{cr3().token} & 7_{water().token}",
        E0=1.33,
        name="Cr2O7-2/Cr+3",
    )


# NO3- + 4 H+ + 3 e- ⇌ NO(g) + 2 H2O(l).
@library
def nitrate():
    return hr(
        f"{no3_minus().token} & 4_{h_plus().token} & 3_@e = {no().token} & 2_{water().token}",
        E0=0.957,
        name="NO3-/NO",
    )


# NO3- + 3 H+ + 2 e- ⇌ HNO2(aq) + H2O(l).
@library
def nitrate_nitrous():
    return hr(
        f"{no3_minus().token} & 3_{h_plus().token} & 2_@e = {hno2().token} & {water().token}",
        E0=0.934,
        name="NO3-/HNO2",
    )


# SO4-2 + 4 H+ + 2 e- ⇌ SO2(aq) + 2 H2O(l).
@library
def sulfate():
    return hr(
        f"{so4_2minus().token} & 4_{h_plus().token} & 2_@e = SO2.aq & 2_{water().token}",
        E0=0.17,
        name="SO4-2/SO2",
    )


# S(s) + 2 H+ + 2 e- ⇌ H2S(aq).
@library
def sulfur():
    return hr(
        f"S.s & 2_{h_plus().token} & 2_@e = H2S.aq",
        E0=0.142,
        name="S/H2S",
    )


# S4O6-2 + 2 e- ⇌ 2 S2O3-2.
@library
def tetrathionate():
    return hr(
        f"S4O6-2 & 2_@e = 2_{s2o3_2minus().token}",
        E0=0.08,
        name="S4O6-2/S2O3-2",
    )


__all__ = [
    "dichromate",
    "manganese_dioxide",
    "nitrate",
    "nitrate_nitrous",
    "permanganate",
    "permanganate_dioxide",
    "sulfate",
    "sulfur",
    "tetrathionate",
]
