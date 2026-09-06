"""Neutral molecules: diatomics through common solvents, acids, and salts.

Melting and boiling points are 1 atm values in kelvin. Species that sublime
or decompose instead of a clean 1 atm melt/boil omit mp/bp.
"""

from ._builders import compound

# --- diatomics ---
h2 = compound("H2", "g", mp=13.99, bp=20.271)
n2 = compound("N2", "g", mp=63.15, bp=77.355)
o2 = compound("O2", "g", mp=54.36, bp=90.188)
f2 = compound("F2", "g", mp=53.48, bp=85.03)
cl2 = compound("Cl2", "g", mp=171.6, bp=239.11)
br2 = compound("Br2", "l", mp=265.8, bp=332.0)
i2 = compound("I2", "s", mp=386.85, bp=457.4)
co = compound("CO", "g", mp=68.15, bp=81.65)
no = compound("NO", "g", mp=109.5, bp=121.4)
hf = compound("HF", "g", mp=189.60, bp=292.67)
hcl = compound("HCl", "g", mp=158.97, bp=188.10)
hbr = compound("HBr", "g", mp=186.0, bp=206.4)
hi = compound("HI", "g", mp=222.35, bp=237.6)

# --- triatomics / small inorganics ---
h2o = compound("H2O", "l", mp=273.15, bp=373.15)
water = h2o
co2 = compound("CO2", "g")
n2o = compound("N2O", "g", mp=182.33, bp=184.67)
no2 = compound("NO2", "g", mp=261.88, bp=294.0)
n2o4 = compound("N2O4", "g", mp=261.9, bp=294.3)
o3 = compound("O3", "g", mp=80.7, bp=161.3)
h2s = compound("H2S", "g", mp=187.66, bp=212.84)
so2 = compound("SO2", "g", mp=197.65, bp=263.13)
so3 = compound("SO3", "l", mp=289.95, bp=317.90)
hcn = compound("HCN", "l", mp=259.86, bp=298.85)
cs2 = compound("CS2", "l", mp=161.58, bp=319.22)
ocs = compound("OCS", "g", mp=134.3, bp=222.9)
h2o2 = compound("H2O2", "l", mp=272.72, bp=423.35)
nh3 = compound("NH3", "g", mp=195.40, bp=239.82)
ph3 = compound("PH3", "g", mp=139.4, bp=185.5)

# --- hydrides / C1–C2 organics ---
ch4 = compound("CH4", "g", mp=90.67, bp=111.66)
c2h6 = compound("C2H6", "g", mp=90.35, bp=184.55)
c2h4 = compound("C2H4", "g", mp=104.0, bp=169.4)
c2h2 = compound("C2H2", "g")
c3h8 = compound("C3H8", "g", mp=85.5, bp=231.1)
c3h6 = compound("C3H6", "g", mp=87.9, bp=225.5)
ch3oh = compound("CH3OH", "l", mp=175.47, bp=337.85)
methanol = ch3oh
c2h5oh = compound("C2H5OH", "l", mp=159.05, bp=351.52)
ethanol = c2h5oh
isopropanol = compound("C3H7OH", "l", mp=184.65, bp=355.4)
hcho = compound("HCHO", "g", mp=181.15, bp=254.05)
ch3cho = compound("CH3CHO", "l", mp=150.15, bp=293.35)
hcooh = compound("HCOOH", "l", mp=281.55, bp=373.95)
ch3cooh = compound("CH3COOH", "l", mp=289.81, bp=391.2)
acetic_acid = ch3cooh
acetone = compound("CH3COCH3", "l", mp=178.5, bp=329.22)

# --- common solvents ---
benzene = compound("C6H6", "l", mp=278.68, bp=353.2)
toluene = compound("C7H8", "l", mp=178.2, bp=383.8)
hexane = compound("C6H14", "l", mp=177.8, bp=341.9)
chloroform = compound("CHCl3", "l", mp=209.6, bp=334.3)
ccl4 = compound("CCl4", "l", mp=250.2, bp=349.8)
dichloromethane = compound("CH2Cl2", "l", mp=176.5, bp=313.0)
diethyl_ether = compound("C2H5OC2H5", "l", mp=156.85, bp=307.6)
ethyl_acetate = compound("CH3COOC2H5", "l", mp=189.6, bp=350.3)
acetonitrile = compound("CH3CN", "l", mp=227.4, bp=354.8)
thf = compound("C4H8O", "l", mp=164.7, bp=339.0)
dmf = compound("HCON(CH3)2", "l", mp=212.7, bp=426.0)
dmso = compound("CH3SOCH3", "l", mp=291.7, bp=462.0)
phenol = compound("C6H5OH", "s", mp=314.0, bp=454.9)
aniline = compound("C6H5NH2", "l", mp=267.1, bp=457.3)
pyridine = compound("C5H5N", "l", mp=231.6, bp=388.4)
n2h4 = compound("N2H4", "l", mp=275.15, bp=386.7)
urea = compound("CO(NH2)2", "s", mp=406.0)
glucose = compound("C6H12O6", "s", mp=419.0)
hocl = compound("HOCl", "aq")
bf3 = compound("BF3", "g", mp=146.8, bp=173.2)
pcl3 = compound("PCl3", "l", mp=161.2, bp=349.3)
pcl5 = compound("PCl5", "s")
sicl4 = compound("SiCl4", "l", mp=204.3, bp=330.8)

# --- oxoacids ---
hno3 = compound("HNO3", "l", mp=231.55, bp=356.15)
nitric_acid = hno3
h2so4 = compound("H2SO4", "l", mp=283.46, bp=610.0)
sulfuric_acid = h2so4
h3po4 = compound("H3PO4", "s", mp=315.5)
phosphoric_acid = h3po4

# --- molecular elements and common salts / oxides ---
p4 = compound("P4", "s", mp=317.3, bp=553.7)
s8 = compound("S8", "s", mp=388.36, bp=717.8)
nacl = compound("NaCl", "s", mp=1073.8, bp=1738.0)
kcl = compound("KCl", "s", mp=1043.0, bp=1693.0)
nabr = compound("NaBr", "s", mp=1020.0, bp=1663.0)
kbr = compound("KBr", "s", mp=1007.0, bp=1656.0)
nai = compound("NaI", "s", mp=933.0, bp=1577.0)
ki = compound("KI", "s", mp=954.0, bp=1593.0)
naf = compound("NaF", "s", mp=1266.0, bp=1978.0)
naoh = compound("NaOH", "s", mp=596.0, bp=1661.0)
koh = compound("KOH", "s", mp=633.0, bp=1600.0)
cao = compound("CaO", "s", mp=2886.0, bp=3123.0)
mgo = compound("MgO", "s", mp=3125.0, bp=3873.0)
al2o3 = compound("Al2O3", "s", mp=2345.0, bp=3250.0)
sio2 = compound("SiO2", "s", mp=1986.0, bp=3223.0)
fe2o3 = compound("Fe2O3", "s", mp=1812.0)
caco3 = compound("CaCO3", "s")
cacl2 = compound("CaCl2", "s", mp=1045.0, bp=1935.0)
mgcl2 = compound("MgCl2", "s", mp=987.0, bp=1685.0)
fecl2 = compound("FeCl2", "s", mp=950.0, bp=1296.0)
fecl3 = compound("FeCl3", "s", mp=577.0, bp=589.0)
cucl2 = compound("CuCl2", "s", mp=771.0, bp=1266.0)
zncl2 = compound("ZnCl2", "s", mp=563.0, bp=1005.0)
agno3 = compound("AgNO3", "s", mp=485.0)
kno3 = compound("KNO3", "s", mp=607.0)
na2so4 = compound("Na2SO4", "s", mp=1157.0)
baso4 = compound("BaSO4", "s", mp=1853.0)
cuso4 = compound("CuSO4", "s")
nh4cl = compound("NH4Cl", "s")

__all__ = [
    "acetic_acid",
    "acetone",
    "acetonitrile",
    "agno3",
    "al2o3",
    "aniline",
    "baso4",
    "benzene",
    "bf3",
    "br2",
    "c2h2",
    "c2h4",
    "c2h5oh",
    "c2h6",
    "c3h6",
    "c3h8",
    "caco3",
    "cacl2",
    "cao",
    "ccl4",
    "ch3cho",
    "ch3cooh",
    "ch3oh",
    "ch4",
    "cl2",
    "chloroform",
    "co",
    "co2",
    "cs2",
    "cucl2",
    "cuso4",
    "dichloromethane",
    "diethyl_ether",
    "dmf",
    "dmso",
    "ethanol",
    "ethyl_acetate",
    "f2",
    "fe2o3",
    "fecl2",
    "fecl3",
    "glucose",
    "h2",
    "h2o",
    "h2o2",
    "h2s",
    "h2so4",
    "h3po4",
    "hbr",
    "hcho",
    "hcl",
    "hcn",
    "hcooh",
    "hexane",
    "hf",
    "hi",
    "hno3",
    "hocl",
    "i2",
    "isopropanol",
    "kbr",
    "kcl",
    "ki",
    "kno3",
    "koh",
    "methanol",
    "mgcl2",
    "mgo",
    "n2",
    "n2h4",
    "n2o",
    "n2o4",
    "na2so4",
    "nabr",
    "nacl",
    "naf",
    "nai",
    "naoh",
    "nh3",
    "nh4cl",
    "nitric_acid",
    "no",
    "no2",
    "o2",
    "o3",
    "ocs",
    "p4",
    "pcl3",
    "pcl5",
    "ph3",
    "phenol",
    "phosphoric_acid",
    "pyridine",
    "s8",
    "sicl4",
    "sio2",
    "so2",
    "so3",
    "sulfuric_acid",
    "thf",
    "toluene",
    "urea",
    "water",
    "zncl2",
]
