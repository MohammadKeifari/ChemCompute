"""Aqueous ions. No mp/bp. Colored aqua ions and oxoanions carry UV-Vis envelopes."""

from ._builders import compound, vis_bands, vis_nm

# --- hydrogen / hydroxide ---
h_plus = compound("H+", "aq", charge=1)
h3o = compound("H3O+", "aq", charge=1)
oh_minus = compound("OH-", "aq", charge=-1)

# --- alkali / alkaline earth ---
li_plus = compound("Li+", "aq", charge=1)
na_plus = compound("Na+", "aq", charge=1)
k_plus = compound("K+", "aq", charge=1)
rb_plus = compound("Rb+", "aq", charge=1)
cs_plus = compound("Cs+", "aq", charge=1)
be2 = compound("Be+2", "aq", charge=2)
mg2 = compound("Mg+2", "aq", charge=2)
ca2 = compound("Ca+2", "aq", charge=2)
sr2 = compound("Sr+2", "aq", charge=2)
ba2 = compound("Ba+2", "aq", charge=2)

# --- p-block / post-transition ---
al3 = compound("Al+3", "aq", charge=3)
sn2 = compound("Sn+2", "aq", charge=2)
sn4 = compound("Sn+4", "aq", charge=4)
pb2 = compound("Pb+2", "aq", charge=2)
bi3 = compound("Bi+3", "aq", charge=3)
ag_plus = compound("Ag+", "aq", charge=1)
cd2 = compound("Cd+2", "aq", charge=2)
hg2 = compound("Hg+2", "aq", charge=2)
zn2 = compound("Zn+2", "aq", charge=2)

# --- d-block aqua ions (ligand-field envelopes around published λmax/ε) ---
ti3 = compound("Ti+3", "aq", charge=3, spectrum=vis_bands((500, 6.0)))
v2 = compound("V+2", "aq", charge=2, spectrum=vis_bands((560, 4.0), (850, 2.0)))
v3 = compound("V+3", "aq", charge=3, spectrum=vis_bands((400, 8.0), (580, 6.0)))
vo2 = compound("VO+2", "aq", charge=2, spectrum=vis_bands((760, 17.0), fwhm_cm=4000.0))
cr2 = compound("Cr+2", "aq", charge=2, spectrum=vis_bands((710, 5.0)))
cr3 = compound("Cr+3", "aq", charge=3, spectrum=vis_bands((407, 15.0), (575, 13.0)))
mn2 = compound("Mn+2", "aq", charge=2)
fe2 = compound("Fe+2", "aq", charge=2, spectrum=vis_bands((1000, 1.6), fwhm_cm=3000.0))
fe3 = compound("Fe+3", "aq", charge=3)
co_2plus = compound("Co+2", "aq", charge=2, spectrum=vis_bands((510, 4.8)))
co_3plus = compound("Co+3", "aq", charge=3)
ni2 = compound("Ni+2", "aq", charge=2, spectrum=vis_bands((395, 5.0), (658, 2.0)))
cu_plus = compound("Cu+", "aq", charge=1)
cu2 = compound("Cu+2", "aq", charge=2, spectrum=vis_bands((800, 12.0), fwhm_cm=4000.0))

# --- halides / simple anions ---
f_minus = compound("F-", "aq", charge=-1)
cl_minus = compound("Cl-", "aq", charge=-1)
br_minus = compound("Br-", "aq", charge=-1)
i_minus = compound("I-", "aq", charge=-1)
hs_minus = compound("HS-", "aq", charge=-1)
s_2minus = compound("S-2", "aq", charge=-2)
cn_minus = compound("CN-", "aq", charge=-1)
scn_minus = compound("SCN-", "aq", charge=-1)
ocn_minus = compound("OCN-", "aq", charge=-1)
n3_minus = compound("N3-", "aq", charge=-1)
nh4 = compound("NH4+", "aq", charge=1)

# --- oxoanions ---
no2_minus = compound("NO2-", "aq", charge=-1)
no3_minus = compound("NO3-", "aq", charge=-1)
so3_2minus = compound("SO3-2", "aq", charge=-2)
hso3_minus = compound("HSO3-", "aq", charge=-1)
so4_2minus = compound("SO4-2", "aq", charge=-2)
hso4_minus = compound("HSO4-", "aq", charge=-1)
s2o3_2minus = compound("S2O3-2", "aq", charge=-2)
co3_2minus = compound("CO3-2", "aq", charge=-2)
hco3_minus = compound("HCO3-", "aq", charge=-1)
po4_3minus = compound("PO4-3", "aq", charge=-3)
hpo4_2minus = compound("HPO4-2", "aq", charge=-2)
h2po4_minus = compound("H2PO4-", "aq", charge=-1)
clo_minus = compound("ClO-", "aq", charge=-1)
clo2_minus = compound("ClO2-", "aq", charge=-1)
clo3_minus = compound("ClO3-", "aq", charge=-1)
clo4_minus = compound("ClO4-", "aq", charge=-1)
ch3coo_minus = compound("CH3COO-", "aq", charge=-1)
c2o4_2minus = compound("C2O4-2", "aq", charge=-2)

# colored oxoanions / polyhalides — connected envelopes, λmax/ε from analytical tables
# MnO4- visible vibronic structure plus the near-UV charge-transfer band
mno4 = compound(
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
cro4 = compound(
    "CrO4-2",
    "aq",
    charge=-2,
    spectrum=vis_bands((273, 3620.0, 4000.0), (372, 4830.0, 2800.0)),
)
cr2o7 = compound(
    "Cr2O7-2",
    "aq",
    charge=-2,
    spectrum=vis_bands((257, 2300.0, 4000.0), (350, 1570.0, 3800.0), (440, 370.0, 3000.0)),
)
i3 = compound(
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
