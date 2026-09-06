"""Water, hydrogen, oxygen, and peroxide couples vs SHE (25 °C)."""

from ..compounds import h2, h_plus, o2, oh_minus, water
from ._builders import hr, library


# 2 H+ + 2 e- ⇌ H2(g). E° = 0 V by definition of SHE.
@library
def hydrogen():
    return hr(
        f"2_{h_plus().token} & 2_@e = {h2().token}",
        E0=0.000,
        name="H+/H2",
    )


she = hydrogen


# O2(g) + 4 H+ + 4 e- ⇌ 2 H2O(l).
@library
def oxygen():
    return hr(
        f"{o2().token} & 4_{h_plus().token} & 4_@e = 2_{water().token}",
        E0=1.229,
        name="O2/H2O",
    )


# O2(g) + 2 H2O(l) + 4 e- ⇌ 4 OH-.
@library
def oxygen_hydroxide():
    return hr(
        f"{o2().token} & 2_{water().token} & 4_@e = 4_{oh_minus().token}",
        E0=0.401,
        name="O2/OH-",
    )


# O2(g) + 2 H+ + 2 e- ⇌ H2O2(aq). Library H2O2 is the neat liquid; keep aqueous in Q.
@library
def hydrogen_peroxide():
    return hr(
        f"{o2().token} & 2_{h_plus().token} & 2_@e = H2O2.aq",
        E0=0.695,
        name="O2/H2O2",
    )


# H2O2(aq) + 2 H+ + 2 e- ⇌ 2 H2O(l).
@library
def peroxide_water():
    return hr(
        f"H2O2.aq & 2_{h_plus().token} & 2_@e = 2_{water().token}",
        E0=1.776,
        name="H2O2/H2O",
    )


__all__ = [
    "hydrogen",
    "hydrogen_peroxide",
    "oxygen",
    "oxygen_hydroxide",
    "peroxide_water",
    "she",
]
