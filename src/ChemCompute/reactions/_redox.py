"""Strongly driven aqueous redox titrations. Concentrations start at 0.

These use ``infinite_K`` (analytical completeness), not a tabulated Kc.
Water is the library liquid so it drops out of Q.
"""

from ..compounds import (
    cr2o7,
    cr3,
    cu2,
    fe2,
    fe3,
    h_plus,
    i_minus,
    mn2,
    mno4,
    water,
    zn2,
)
from ._builders import library, rxn


# MnO4- + 8 H+ + 5 Fe+2 → Mn+2 + 5 Fe+3 + 4 H2O(l).
@library
def permanganate_iron():
    return rxn(
        f"{mno4().token} & 8_{h_plus().token} & 5_{fe2().token} > {mn2().token} & 5_{fe3().token} & 4_{water().token}",
        infinite_K=True,
    )


# Cr2O7-2 + 14 H+ + 6 Fe+2 → 2 Cr+3 + 6 Fe+3 + 7 H2O(l).
@library
def dichromate_iron():
    return rxn(
        f"{cr2o7().token} & 14_{h_plus().token} & 6_{fe2().token} > 2_{cr3().token} & 6_{fe3().token} & 7_{water().token}",
        infinite_K=True,
    )


# Cu+2 + Zn(s) → Cu(s) + Zn+2.
@library
def copper_zinc():
    return rxn(
        f"{cu2().token} & Zn.s > Cu.s & {zn2().token}",
        infinite_K=True,
    )


# I2(aq) + 2 S2O3-2 → 2 I- + S4O6-2.
@library
def iodine_thiosulfate():
    return rxn(
        f"I2.aq & 2_S2O3-2 > 2_{i_minus().token} & S4O6-2",
        infinite_K=True,
    )


__all__ = [
    "copper_zinc",
    "dichromate_iron",
    "iodine_thiosulfate",
    "permanganate_iron",
]
