"""Common aqueous complex-formation equilibria (overall β / Kf, 25 °C). Concentrations start at 0."""

from ..compounds import (
    ag_nh3_2_plus,
    ag_plus,
    agcl2,
    cl_minus,
    cn_minus,
    co_3plus,
    co_nh3_6_3plus,
    cu2,
    cu_nh3_4_2plus,
    en,
    fe2,
    fe3,
    fe_cn6_3minus,
    fe_cn6_4minus,
    ferroin,
    fescn,
    i3,
    i_minus,
    ni2,
    ni_en3_2plus,
    ni_nh3_6_2plus,
    phen,
    scn_minus,
    zn2,
    zn_nh3_4_2plus,
)
from ._builders import library, rxn


# Fe+3 + SCN- ⇌ FeSCN+2. Kf = 10^3.22.
@library
def fescn_kf():
    return rxn(
        f"{fe3().token} & {scn_minus().token} > {fescn().token}",
        K=10 ** 3.22,
    )


# Cu+2 + 4 NH3(aq) ⇌ Cu(NH3)4+2. β4 = 2.1e13.
@library
def cu_nh3_4_kf():
    return rxn(
        f"{cu2().token} & 4_NH3.aq > {cu_nh3_4_2plus().token}",
        K=2.1e13,
    )


# Ni+2 + 6 NH3(aq) ⇌ Ni(NH3)6+2. β6 = 5.5e8.
@library
def ni_nh3_6_kf():
    return rxn(
        f"{ni2().token} & 6_NH3.aq > {ni_nh3_6_2plus().token}",
        K=5.5e8,
    )


# Ag+ + 2 NH3(aq) ⇌ Ag(NH3)2+. β2 = 1.7e7.
@library
def ag_nh3_2_kf():
    return rxn(
        f"{ag_plus().token} & 2_NH3.aq > {ag_nh3_2_plus().token}",
        K=1.7e7,
    )


# Ag+ + 2 Cl- ⇌ AgCl2-. β2 = 1.8e5.
@library
def agcl2_kf():
    return rxn(
        f"{ag_plus().token} & 2_{cl_minus().token} > {agcl2().token}",
        K=1.8e5,
    )


# Zn+2 + 4 NH3(aq) ⇌ Zn(NH3)4+2. β4 = 2.9e9.
@library
def zn_nh3_4_kf():
    return rxn(
        f"{zn2().token} & 4_NH3.aq > {zn_nh3_4_2plus().token}",
        K=2.9e9,
    )


# Co+3 + 6 NH3(aq) ⇌ Co(NH3)6+3. β6 = 2.3e33.
@library
def co_nh3_6_kf():
    return rxn(
        f"{co_3plus().token} & 6_NH3.aq > {co_nh3_6_3plus().token}",
        K=2.3e33,
    )


# Fe+2 + 6 CN- ⇌ Fe(CN)6-4. β6 = 1e35.
@library
def fe_cn6_4_kf():
    return rxn(
        f"{fe2().token} & 6_{cn_minus().token} > {fe_cn6_4minus().token}",
        K=1.0e35,
    )


# Fe+3 + 6 CN- ⇌ Fe(CN)6-3. β6 = 1e42.
@library
def fe_cn6_3_kf():
    return rxn(
        f"{fe3().token} & 6_{cn_minus().token} > {fe_cn6_3minus().token}",
        K=1.0e42,
    )


# I2(aq) + I- ⇌ I3-. K = 710.
@library
def i3_kf():
    return rxn(
        f"I2.aq & {i_minus().token} > {i3().token}",
        K=710.0,
    )


# Fe+2 + 3 phen ⇌ Fe(phen)3+2. β3 = 2e21.
@library
def ferroin_kf():
    return rxn(
        f"{fe2().token} & 3_{phen().token} > {ferroin().token}",
        K=2.0e21,
    )


# Ni+2 + 3 en ⇌ Ni(en)3+2. β3 = 2e18.
@library
def ni_en3_kf():
    return rxn(
        f"{ni2().token} & 3_{en().token} > {ni_en3_2plus().token}",
        K=2.0e18,
    )


__all__ = [
    "ag_nh3_2_kf",
    "agcl2_kf",
    "co_nh3_6_kf",
    "cu_nh3_4_kf",
    "fe_cn6_3_kf",
    "fe_cn6_4_kf",
    "ferroin_kf",
    "fescn_kf",
    "i3_kf",
    "ni_en3_kf",
    "ni_nh3_6_kf",
    "zn_nh3_4_kf",
]
