"""Common Ksp dissolution equilibria (25 °C). Solid is omitted from Q. Concentrations start at 0."""

from ..compounds import (
    ag2cro4,
    ag2s,
    ag2so4,
    ag_plus,
    agbr,
    agcl,
    agi,
    al3,
    aloh3,
    ba2,
    baco3,
    baf2,
    baso4,
    ca2,
    caf2,
    caco3,
    caoh2,
    caso4,
    cd2,
    cds,
    cl_minus,
    co3_2minus,
    cro4,
    cu2,
    cuoh2,
    cus,
    f_minus,
    fe2,
    fe3,
    feoh2,
    feoh3,
    fes,
    hg2,
    hgs,
    i_minus,
    mg2,
    mgco3,
    mgoh2,
    mn2,
    mnoh2,
    ni2,
    nioh2,
    oh_minus,
    pb2,
    pbcl2,
    pbi2,
    pbs,
    pbso4,
    s_2minus,
    so4_2minus,
    sr2,
    srf2,
    srso4,
    zn2,
    znoh2,
    zns,
)
from ._builders import library, rxn


# AgCl(s) ⇌ Ag+ + Cl-. Ksp = 1.8e-10.
@library
def agcl_ksp():
    return rxn(
        f"{agcl().token} > {ag_plus().token} & {cl_minus().token}",
        K=1.8e-10,
    )


# AgBr(s) ⇌ Ag+ + Br-. Ksp = 5.0e-13.
@library
def agbr_ksp():
    return rxn(
        f"{agbr().token} > {ag_plus().token} & Br-",
        K=5.0e-13,
    )


# AgI(s) ⇌ Ag+ + I-. Ksp = 8.3e-17.
@library
def agi_ksp():
    return rxn(
        f"{agi().token} > {ag_plus().token} & {i_minus().token}",
        K=8.3e-17,
    )


# Ag2CrO4(s) ⇌ 2 Ag+ + CrO4-2. Ksp = 1.2e-12.
@library
def ag2cro4_ksp():
    return rxn(
        f"{ag2cro4().token} > 2_{ag_plus().token} & {cro4().token}",
        K=1.2e-12,
    )


# Ag2S(s) ⇌ 2 Ag+ + S-2. Ksp = 6e-51.
@library
def ag2s_ksp():
    return rxn(
        f"{ag2s().token} > 2_{ag_plus().token} & {s_2minus().token}",
        K=6.0e-51,
    )


# Ag2SO4(s) ⇌ 2 Ag+ + SO4-2. Ksp = 1.2e-5.
@library
def ag2so4_ksp():
    return rxn(
        f"{ag2so4().token} > 2_{ag_plus().token} & {so4_2minus().token}",
        K=1.2e-5,
    )


# BaSO4(s) ⇌ Ba+2 + SO4-2. Ksp = 1.1e-10.
@library
def baso4_ksp():
    return rxn(
        f"{baso4().token} > {ba2().token} & {so4_2minus().token}",
        K=1.1e-10,
    )


# BaCO3(s) ⇌ Ba+2 + CO3-2. Ksp = 2.6e-9.
@library
def baco3_ksp():
    return rxn(
        f"{baco3().token} > {ba2().token} & {co3_2minus().token}",
        K=2.6e-9,
    )


# BaF2(s) ⇌ Ba+2 + 2 F-. Ksp = 1.8e-7.
@library
def baf2_ksp():
    return rxn(
        f"{baf2().token} > {ba2().token} & 2_{f_minus().token}",
        K=1.8e-7,
    )


# CaF2(s) ⇌ Ca+2 + 2 F-. Ksp = 3.9e-11.
@library
def caf2_ksp():
    return rxn(
        f"{caf2().token} > {ca2().token} & 2_{f_minus().token}",
        K=3.9e-11,
    )


# CaCO3(s) ⇌ Ca+2 + CO3-2. Ksp = 4.5e-9.
@library
def caco3_ksp():
    return rxn(
        f"{caco3().token} > {ca2().token} & {co3_2minus().token}",
        K=4.5e-9,
    )


# CaSO4(s) ⇌ Ca+2 + SO4-2. Ksp = 2.4e-5.
@library
def caso4_ksp():
    return rxn(
        f"{caso4().token} > {ca2().token} & {so4_2minus().token}",
        K=2.4e-5,
    )


# Ca(OH)2(s) ⇌ Ca+2 + 2 OH-. Ksp = 5.5e-6.
@library
def caoh2_ksp():
    return rxn(
        f"{caoh2().token} > {ca2().token} & 2_{oh_minus().token}",
        K=5.5e-6,
    )


# Mg(OH)2(s) ⇌ Mg+2 + 2 OH-. Ksp = 5.6e-12.
@library
def mgoh2_ksp():
    return rxn(
        f"{mgoh2().token} > {mg2().token} & 2_{oh_minus().token}",
        K=5.6e-12,
    )


# MgCO3(s) ⇌ Mg+2 + CO3-2. Ksp = 3.5e-8.
@library
def mgco3_ksp():
    return rxn(
        f"{mgco3().token} > {mg2().token} & {co3_2minus().token}",
        K=3.5e-8,
    )


# Fe(OH)3(s) ⇌ Fe+3 + 3 OH-. Ksp = 1e-38.
@library
def feoh3_ksp():
    return rxn(
        f"{feoh3().token} > {fe3().token} & 3_{oh_minus().token}",
        K=1.0e-38,
    )


# Fe(OH)2(s) ⇌ Fe+2 + 2 OH-. Ksp = 4.9e-17.
@library
def feoh2_ksp():
    return rxn(
        f"{feoh2().token} > {fe2().token} & 2_{oh_minus().token}",
        K=4.9e-17,
    )


# Al(OH)3(s) ⇌ Al+3 + 3 OH-. Ksp = 3e-34.
@library
def aloh3_ksp():
    return rxn(
        f"{aloh3().token} > {al3().token} & 3_{oh_minus().token}",
        K=3.0e-34,
    )


# Zn(OH)2(s) ⇌ Zn+2 + 2 OH-. Ksp = 3e-17.
@library
def znoh2_ksp():
    return rxn(
        f"{znoh2().token} > {zn2().token} & 2_{oh_minus().token}",
        K=3.0e-17,
    )


# Cu(OH)2(s) ⇌ Cu+2 + 2 OH-. Ksp = 2.2e-20.
@library
def cuoh2_ksp():
    return rxn(
        f"{cuoh2().token} > {cu2().token} & 2_{oh_minus().token}",
        K=2.2e-20,
    )


# Ni(OH)2(s) ⇌ Ni+2 + 2 OH-. Ksp = 5.5e-16.
@library
def nioh2_ksp():
    return rxn(
        f"{nioh2().token} > {ni2().token} & 2_{oh_minus().token}",
        K=5.5e-16,
    )


# Mn(OH)2(s) ⇌ Mn+2 + 2 OH-. Ksp = 1.6e-13.
@library
def mnoh2_ksp():
    return rxn(
        f"{mnoh2().token} > {mn2().token} & 2_{oh_minus().token}",
        K=1.6e-13,
    )


# PbS(s) ⇌ Pb+2 + S-2. Ksp = 3e-28.
@library
def pbs_ksp():
    return rxn(
        f"{pbs().token} > {pb2().token} & {s_2minus().token}",
        K=3.0e-28,
    )


# PbI2(s) ⇌ Pb+2 + 2 I-. Ksp = 7.1e-9.
@library
def pbi2_ksp():
    return rxn(
        f"{pbi2().token} > {pb2().token} & 2_{i_minus().token}",
        K=7.1e-9,
    )


# PbCl2(s) ⇌ Pb+2 + 2 Cl-. Ksp = 1.7e-5.
@library
def pbcl2_ksp():
    return rxn(
        f"{pbcl2().token} > {pb2().token} & 2_{cl_minus().token}",
        K=1.7e-5,
    )


# PbSO4(s) ⇌ Pb+2 + SO4-2. Ksp = 1.6e-8.
@library
def pbso4_ksp():
    return rxn(
        f"{pbso4().token} > {pb2().token} & {so4_2minus().token}",
        K=1.6e-8,
    )


# CuS(s) ⇌ Cu+2 + S-2. Ksp = 6e-37.
@library
def cus_ksp():
    return rxn(
        f"{cus().token} > {cu2().token} & {s_2minus().token}",
        K=6.0e-37,
    )


# ZnS(s) ⇌ Zn+2 + S-2. Ksp = 2e-25.
@library
def zns_ksp():
    return rxn(
        f"{zns().token} > {zn2().token} & {s_2minus().token}",
        K=2.0e-25,
    )


# CdS(s) ⇌ Cd+2 + S-2. Ksp = 8e-28.
@library
def cds_ksp():
    return rxn(
        f"{cds().token} > {cd2().token} & {s_2minus().token}",
        K=8.0e-28,
    )


# FeS(s) ⇌ Fe+2 + S-2. Ksp = 6e-19.
@library
def fes_ksp():
    return rxn(
        f"{fes().token} > {fe2().token} & {s_2minus().token}",
        K=6.0e-19,
    )


# HgS(s) ⇌ Hg+2 + S-2. Ksp = 4e-53.
@library
def hgs_ksp():
    return rxn(
        f"{hgs().token} > {hg2().token} & {s_2minus().token}",
        K=4.0e-53,
    )


# SrF2(s) ⇌ Sr+2 + 2 F-. Ksp = 4.3e-9.
@library
def srf2_ksp():
    return rxn(
        f"{srf2().token} > {sr2().token} & 2_{f_minus().token}",
        K=4.3e-9,
    )


# SrSO4(s) ⇌ Sr+2 + SO4-2. Ksp = 3.4e-7.
@library
def srso4_ksp():
    return rxn(
        f"{srso4().token} > {sr2().token} & {so4_2minus().token}",
        K=3.4e-7,
    )


__all__ = [
    "ag2cro4_ksp",
    "ag2s_ksp",
    "ag2so4_ksp",
    "agbr_ksp",
    "agcl_ksp",
    "agi_ksp",
    "aloh3_ksp",
    "baco3_ksp",
    "baf2_ksp",
    "baso4_ksp",
    "caco3_ksp",
    "caf2_ksp",
    "caoh2_ksp",
    "caso4_ksp",
    "cds_ksp",
    "cuoh2_ksp",
    "cus_ksp",
    "feoh2_ksp",
    "feoh3_ksp",
    "fes_ksp",
    "hgs_ksp",
    "mgco3_ksp",
    "mgoh2_ksp",
    "mnoh2_ksp",
    "nioh2_ksp",
    "pbcl2_ksp",
    "pbi2_ksp",
    "pbs_ksp",
    "pbso4_ksp",
    "srf2_ksp",
    "srso4_ksp",
    "znoh2_ksp",
    "zns_ksp",
]
