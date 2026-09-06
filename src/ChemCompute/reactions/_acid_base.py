"""Common aqueous acid–base equilibria (25 °C). Concentrations start at 0.

Molecular acids that exist in the compound library as a gas, liquid, or solid
are written with an ``.aq`` formula so they stay in Q. Water is the library
liquid, so Kw is [H+][OH-]. Strong acids use ``infinite_K``.
"""

from ..compounds import (
    ch3coo_minus,
    clo_minus,
    co3_2minus,
    h2co3,
    h2po4_minus,
    h_plus,
    hco3_minus,
    hno2,
    hocl,
    hpo4_2minus,
    hs_minus,
    hso3_minus,
    hso4_minus,
    nh4,
    no2_minus,
    oh_minus,
    po4_3minus,
    s_2minus,
    so3_2minus,
    so4_2minus,
    water,
)
from ._builders import library, rxn


# H2O(l) ⇌ H+ + OH-. Kw = 1.0e-14.
@library
def water_kw():
    return rxn(
        f"{water().token} > {h_plus().token} & {oh_minus().token}",
        K=1.0e-14,
    )


# HCl(aq) → H+ + Cl- (strong).
@library
def hydrochloric_acid():
    return rxn("HCl.aq > H+ & Cl-", infinite_K=True)


# HBr(aq) → H+ + Br- (strong).
@library
def hydrobromic_acid():
    return rxn("HBr.aq > H+ & Br-", infinite_K=True)


# HI(aq) → H+ + I- (strong).
@library
def hydroiodic_acid():
    return rxn("HI.aq > H+ & I-", infinite_K=True)


# HNO3(aq) → H+ + NO3- (strong).
@library
def nitric_acid():
    return rxn("HNO3.aq > H+ & NO3-", infinite_K=True)


# HClO4(aq) → H+ + ClO4- (strong).
@library
def perchloric_acid():
    return rxn("HClO4.aq > H+ & ClO4-", infinite_K=True)


# H2SO4(aq) → H+ + HSO4- (strong first proton).
@library
def sulfuric_acid_1():
    return rxn("H2SO4.aq > H+ & HSO4-", infinite_K=True)


# HSO4- ⇌ H+ + SO4-2. Ka2 = 1.2e-2.
@library
def sulfuric_acid_2():
    return rxn(
        f"{hso4_minus().token} > {h_plus().token} & {so4_2minus().token}",
        K=1.2e-2,
    )


# HF(aq) ⇌ H+ + F-. Ka = 6.6e-4.
@library
def hydrogen_fluoride():
    return rxn("HF.aq > H+ & F-", K=6.6e-4)


# HNO2 ⇌ H+ + NO2-. Ka = 4.5e-4.
@library
def nitrous_acid():
    return rxn(
        f"{hno2().token} > {h_plus().token} & {no2_minus().token}",
        K=4.5e-4,
    )


# HCOOH(aq) ⇌ H+ + HCOO-. Ka = 1.8e-4.
@library
def formic_acid():
    return rxn("HCOOH.aq > H+ & HCOO-", K=1.8e-4)


# CH3COOH(aq) ⇌ H+ + CH3COO-. Ka = 1.8e-5.
@library
def acetic_acid():
    return rxn(
        f"CH3COOH.aq > {h_plus().token} & {ch3coo_minus().token}",
        K=1.8e-5,
    )


# H2CO3 ⇌ H+ + HCO3-. Ka1 = 4.45e-7.
@library
def carbonic_acid_1():
    return rxn(
        f"{h2co3().token} > {h_plus().token} & {hco3_minus().token}",
        K=4.45e-7,
    )


# HCO3- ⇌ H+ + CO3-2. Ka2 = 4.69e-11.
@library
def carbonic_acid_2():
    return rxn(
        f"{hco3_minus().token} > {h_plus().token} & {co3_2minus().token}",
        K=4.69e-11,
    )


# H2S(aq) ⇌ H+ + HS-. Ka1 = 8.9e-8.
@library
def hydrogen_sulfide_1():
    return rxn("H2S.aq > H+ & HS-", K=8.9e-8)


# HS- ⇌ H+ + S-2. Ka2 = 1.2e-13.
@library
def hydrogen_sulfide_2():
    return rxn(
        f"{hs_minus().token} > {h_plus().token} & {s_2minus().token}",
        K=1.2e-13,
    )


# H2SO3(aq) ⇌ H+ + HSO3-. Ka1 = 1.4e-2.
@library
def sulfurous_acid_1():
    return rxn("H2SO3.aq > H+ & HSO3-", K=1.4e-2)


# HSO3- ⇌ H+ + SO3-2. Ka2 = 6.3e-8.
@library
def sulfurous_acid_2():
    return rxn(
        f"{hso3_minus().token} > {h_plus().token} & {so3_2minus().token}",
        K=6.3e-8,
    )


# H3PO4(aq) ⇌ H+ + H2PO4-. Ka1 = 7.1e-3.
@library
def phosphoric_acid_1():
    return rxn("H3PO4.aq > H+ & H2PO4-", K=7.1e-3)


# H2PO4- ⇌ H+ + HPO4-2. Ka2 = 6.3e-8.
@library
def phosphoric_acid_2():
    return rxn(
        f"{h2po4_minus().token} > {h_plus().token} & {hpo4_2minus().token}",
        K=6.3e-8,
    )


# HPO4-2 ⇌ H+ + PO4-3. Ka3 = 4.2e-13.
@library
def phosphoric_acid_3():
    return rxn(
        f"{hpo4_2minus().token} > {h_plus().token} & {po4_3minus().token}",
        K=4.2e-13,
    )


# HOCl ⇌ H+ + ClO-. Ka = 3.5e-8.
@library
def hypochlorous_acid():
    return rxn(
        f"{hocl().token} > {h_plus().token} & {clo_minus().token}",
        K=3.5e-8,
    )


# HCN(aq) ⇌ H+ + CN-. Ka = 6.2e-10.
@library
def hydrogen_cyanide():
    return rxn("HCN.aq > H+ & CN-", K=6.2e-10)


# NH4+ ⇌ NH3(aq) + H+. Ka = 5.6e-10.
@library
def ammonium():
    return rxn(
        f"{nh4().token} > NH3.aq & {h_plus().token}",
        K=5.6e-10,
    )


# NH3(aq) + H2O(l) ⇌ NH4+ + OH-. Kb = 1.8e-5.
@library
def ammonia():
    return rxn(
        f"NH3.aq & {water().token} > {nh4().token} & {oh_minus().token}",
        K=1.8e-5,
    )


# C6H5OH(aq) ⇌ H+ + C6H5O-. Ka = 1.3e-10.
@library
def phenol():
    return rxn("C6H5OH.aq > H+ & C6H5O-", K=1.3e-10)


__all__ = [
    "acetic_acid",
    "ammonia",
    "ammonium",
    "carbonic_acid_1",
    "carbonic_acid_2",
    "formic_acid",
    "hydrobromic_acid",
    "hydrochloric_acid",
    "hydrogen_cyanide",
    "hydrogen_fluoride",
    "hydrogen_sulfide_1",
    "hydrogen_sulfide_2",
    "hydroiodic_acid",
    "hypochlorous_acid",
    "nitric_acid",
    "nitrous_acid",
    "perchloric_acid",
    "phenol",
    "phosphoric_acid_1",
    "phosphoric_acid_2",
    "phosphoric_acid_3",
    "sulfuric_acid_1",
    "sulfuric_acid_2",
    "sulfurous_acid_1",
    "sulfurous_acid_2",
    "water_kw",
]
