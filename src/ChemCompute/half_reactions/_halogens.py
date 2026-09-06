"""Halogen / halide couples vs SHE (25 °C)."""

from ..compounds import (
    br2,
    br_minus,
    cl2,
    cl_minus,
    f2,
    f_minus,
    h_plus,
    hocl,
    i2,
    i3,
    i_minus,
    water,
)
from ._builders import hr, library


# F2(g) + 2 e- ⇌ 2 F-.
@library
def fluorine():
    return hr(
        f"{f2().token} & 2_@e = 2_{f_minus().token}",
        E0=2.866,
        name="F2/F-",
    )


# Cl2(g) + 2 e- ⇌ 2 Cl-.
@library
def chlorine():
    return hr(
        f"{cl2().token} & 2_@e = 2_{cl_minus().token}",
        E0=1.358,
        name="Cl2/Cl-",
    )


# Br2(l) + 2 e- ⇌ 2 Br-.
@library
def bromine():
    return hr(
        f"{br2().token} & 2_@e = 2_{br_minus().token}",
        E0=1.066,
        name="Br2/Br-",
    )


# I2(s) + 2 e- ⇌ 2 I-.
@library
def iodine():
    return hr(
        f"{i2().token} & 2_@e = 2_{i_minus().token}",
        E0=0.535,
        name="I2/I-",
    )


# I3- + 2 e- ⇌ 3 I-.
@library
def triiodide():
    return hr(
        f"{i3().token} & 2_@e = 3_{i_minus().token}",
        E0=0.536,
        name="I3-/I-",
    )


# HOCl(aq) + H+ + 2 e- ⇌ Cl- + H2O(l).
@library
def hypochlorite():
    return hr(
        f"{hocl().token} & {h_plus().token} & 2_@e = {cl_minus().token} & {water().token}",
        E0=1.482,
        name="HOCl/Cl-",
    )


__all__ = [
    "bromine",
    "chlorine",
    "fluorine",
    "hypochlorite",
    "iodine",
    "triiodide",
]
