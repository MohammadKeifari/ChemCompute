"""Aqueous ions. No mp/bp. Colored aqua ions and oxoanions carry UV-Vis envelopes."""

from ._builders import compound, library, vis_bands, vis_nm

# --- hydrogen / hydroxide ---
@library
def h_plus():
    return compound("H+", "aq", charge=1)
@library
def h3o():
    return compound("H3O+", "aq", charge=1)
@library
def oh_minus():
    return compound("OH-", "aq", charge=-1)

# --- alkali / alkaline earth ---
@library
def li_plus():
    return compound("Li+", "aq", charge=1)
@library
def na_plus():
    return compound("Na+", "aq", charge=1)
@library
def k_plus():
    return compound("K+", "aq", charge=1)
@library
def rb_plus():
    return compound("Rb+", "aq", charge=1)
@library
def cs_plus():
    return compound("Cs+", "aq", charge=1)
@library
def be2():
    return compound("Be+2", "aq", charge=2)
@library
def mg2():
    return compound("Mg+2", "aq", charge=2)
@library
def ca2():
    return compound("Ca+2", "aq", charge=2)
@library
def sr2():
    return compound("Sr+2", "aq", charge=2)
@library
def ba2():
    return compound("Ba+2", "aq", charge=2)

# --- p-block / post-transition ---
@library
def al3():
    return compound("Al+3", "aq", charge=3)
@library
def sn2():
    return compound("Sn+2", "aq", charge=2)
@library
def sn4():
    return compound("Sn+4", "aq", charge=4)
@library
def pb2():
    return compound("Pb+2", "aq", charge=2)
@library
def bi3():
    return compound("Bi+3", "aq", charge=3)
@library
def ag_plus():
    return compound("Ag+", "aq", charge=1)
@library
def cd2():
    return compound("Cd+2", "aq", charge=2)
@library
def hg2():
    return compound("Hg+2", "aq", charge=2)
@library
def zn2():
    return compound("Zn+2", "aq", charge=2)

# --- d-block aqua ions (ligand-field envelopes around published λmax/ε) ---
@library
def ti3():
    return compound("Ti+3", "aq", charge=3, spectrum=vis_bands((500, 6.0)))
@library
def v2():
    return compound("V+2", "aq", charge=2, spectrum=vis_bands((560, 4.0), (850, 2.0)))
@library
def v3():
    return compound("V+3", "aq", charge=3, spectrum=vis_bands((400, 8.0), (580, 6.0)))
@library
def vo2():
    return compound("VO+2", "aq", charge=2, spectrum=vis_bands((760, 17.0), fwhm_cm=4000.0))
@library
def cr2():
    return compound("Cr+2", "aq", charge=2, spectrum=vis_bands((710, 5.0)))
@library
def cr3():
    return compound("Cr+3", "aq", charge=3, spectrum=vis_bands((407, 15.0), (575, 13.0)))
@library
def mn2():
    return compound("Mn+2", "aq", charge=2)
@library
def fe2():
    return compound("Fe+2", "aq", charge=2, spectrum=vis_bands((1000, 1.6), fwhm_cm=3000.0))
@library
def fe3():
    return compound("Fe+3", "aq", charge=3)
@library
def co_2plus():
    return compound("Co+2", "aq", charge=2, spectrum=vis_bands((510, 4.8)))
@library
def co_3plus():
    return compound("Co+3", "aq", charge=3)
@library
def ni2():
    return compound("Ni+2", "aq", charge=2, spectrum=vis_bands((395, 5.0), (658, 2.0)))
@library
def cu_plus():
    return compound("Cu+", "aq", charge=1)
@library
def cu2():
    return compound("Cu+2", "aq", charge=2, spectrum=vis_bands((800, 12.0), fwhm_cm=4000.0))

# --- halides / simple anions ---
@library
def f_minus():
    return compound("F-", "aq", charge=-1)
@library
def cl_minus():
    return compound("Cl-", "aq", charge=-1)
@library
def br_minus():
    return compound("Br-", "aq", charge=-1)
@library
def i_minus():
    return compound("I-", "aq", charge=-1)
@library
def hs_minus():
    return compound("HS-", "aq", charge=-1)
@library
def s_2minus():
    return compound("S-2", "aq", charge=-2)
@library
def cn_minus():
    return compound("CN-", "aq", charge=-1)
@library
def scn_minus():
    return compound("SCN-", "aq", charge=-1)
@library
def ocn_minus():
    return compound("OCN-", "aq", charge=-1)
@library
def n3_minus():
    return compound("N3-", "aq", charge=-1)
@library
def nh4():
    return compound("NH4+", "aq", charge=1)

# --- oxoanions ---
@library
def no2_minus():
    return compound("NO2-", "aq", charge=-1)
@library
def no3_minus():
    return compound("NO3-", "aq", charge=-1)
@library
def so3_2minus():
    return compound("SO3-2", "aq", charge=-2)
@library
def hso3_minus():
    return compound("HSO3-", "aq", charge=-1)
@library
def so4_2minus():
    return compound("SO4-2", "aq", charge=-2)
@library
def hso4_minus():
    return compound("HSO4-", "aq", charge=-1)
@library
def s2o3_2minus():
    return compound("S2O3-2", "aq", charge=-2)
@library
def co3_2minus():
    return compound("CO3-2", "aq", charge=-2)
@library
def hco3_minus():
    return compound("HCO3-", "aq", charge=-1)
@library
def po4_3minus():
    return compound("PO4-3", "aq", charge=-3)
@library
def hpo4_2minus():
    return compound("HPO4-2", "aq", charge=-2)
@library
def h2po4_minus():
    return compound("H2PO4-", "aq", charge=-1)
@library
def clo_minus():
    return compound("ClO-", "aq", charge=-1)
@library
def clo2_minus():
    return compound("ClO2-", "aq", charge=-1)
@library
def clo3_minus():
    return compound("ClO3-", "aq", charge=-1)
@library
def clo4_minus():
    return compound("ClO4-", "aq", charge=-1)
@library
def ch3coo_minus():
    return compound("CH3COO-", "aq", charge=-1)
@library
def c2o4_2minus():
    return compound("C2O4-2", "aq", charge=-2)

# colored oxoanions / polyhalides — connected envelopes, λmax/ε from analytical tables
# MnO4- visible vibronic structure plus the near-UV charge-transfer band
@library
def mno4():
    return compound(
    "MnO4-",
    "aq",
    charge=-1,
    spectrum=vis_nm(
        (200, 0.0),
        (210, 8000.0),
        (221, 20000.0),
        (240, 3500.0),
        (270, 500.0),
        (290, 900.0),
        (310, 1800.0),
        (340, 700.0),
        (380, 350.0),
        (430, 280.0),
        (470, 700.0),
        (490, 1250.0),
        (505, 2050.0),
        (525, 2450.0),
        (545, 2180.0),
        (565, 1250.0),
        (580, 350.0),
        (610, 40.0),
        (650, 0.0),
    ),
)
@library
def cro4():
    return compound(
    "CrO4-2",
    "aq",
    charge=-2,
    spectrum=vis_bands((273, 3620.0, 4000.0), (372, 4830.0, 2800.0)),
)
@library
def cr2o7():
    return compound(
    "Cr2O7-2",
    "aq",
    charge=-2,
    spectrum=vis_bands((257, 2300.0, 4000.0), (350, 1570.0, 3800.0), (440, 370.0, 3000.0)),
)
@library
def i3():
    return compound(
    "I3-",
    "aq",
    charge=-1,
    spectrum=vis_bands((287, 40000.0, 2500.0), (353, 26400.0, 2800.0)),
)

__all__ = [
    "ag_plus",
    "al3",
    "ba2",
    "be2",
    "bi3",
    "br_minus",
    "c2o4_2minus",
    "ca2",
    "cd2",
    "ch3coo_minus",
    "cl_minus",
    "clo2_minus",
    "clo3_minus",
    "clo4_minus",
    "clo_minus",
    "cn_minus",
    "co_2plus",
    "co_3plus",
    "co3_2minus",
    "cr2",
    "cr2o7",
    "cr3",
    "cro4",
    "cs_plus",
    "cu2",
    "cu_plus",
    "f_minus",
    "fe2",
    "fe3",
    "h2po4_minus",
    "h3o",
    "h_plus",
    "hco3_minus",
    "hg2",
    "hpo4_2minus",
    "hs_minus",
    "hso3_minus",
    "hso4_minus",
    "i3",
    "i_minus",
    "k_plus",
    "li_plus",
    "mg2",
    "mn2",
    "mno4",
    "n3_minus",
    "na_plus",
    "nh4",
    "ni2",
    "no2_minus",
    "no3_minus",
    "ocn_minus",
    "oh_minus",
    "pb2",
    "po4_3minus",
    "rb_plus",
    "s2o3_2minus",
    "s_2minus",
    "scn_minus",
    "sn2",
    "sn4",
    "so3_2minus",
    "so4_2minus",
    "sr2",
    "ti3",
    "v2",
    "v3",
    "vo2",
    "zn2",
]
