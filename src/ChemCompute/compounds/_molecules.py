"""Neutral molecules: diatomics through common solvents, acids, and salts.

Melting and boiling points are 1 atm values in kelvin. Species that sublime
or decompose instead of a clean 1 atm melt/boil omit mp/bp.
"""

from ._builders import compound, library

# --- diatomics ---
@library
def h2():
    return compound("H2", "g", mp=13.99, bp=20.271)
@library
def n2():
    return compound("N2", "g", mp=63.15, bp=77.355)
@library
def o2():
    return compound("O2", "g", mp=54.36, bp=90.188)
@library
def f2():
    return compound("F2", "g", mp=53.48, bp=85.03)
@library
def cl2():
    return compound("Cl2", "g", mp=171.6, bp=239.11)
@library
def br2():
    return compound("Br2", "l", mp=265.8, bp=332.0)
@library
def i2():
    return compound("I2", "s", mp=386.85, bp=457.4)
@library
def co():
    return compound("CO", "g", mp=68.15, bp=81.65)
@library
def no():
    return compound("NO", "g", mp=109.5, bp=121.4)
@library
def hf():
    return compound("HF", "g", mp=189.60, bp=292.67)
@library
def hcl():
    return compound("HCl", "g", mp=158.97, bp=188.10)
@library
def hbr():
    return compound("HBr", "g", mp=186.0, bp=206.4)
@library
def hi():
    return compound("HI", "g", mp=222.35, bp=237.6)

# --- triatomics / small inorganics ---
@library
def h2o():
    return compound("H2O", "l", mp=273.15, bp=373.15)
water = h2o
@library
def co2():
    return compound("CO2", "g")
@library
def n2o():
    return compound("N2O", "g", mp=182.33, bp=184.67)
@library
def no2():
    return compound("NO2", "g", mp=261.88, bp=294.0)
@library
def n2o4():
    return compound("N2O4", "g", mp=261.9, bp=294.3)
@library
def o3():
    return compound("O3", "g", mp=80.7, bp=161.3)
@library
def h2s():
    return compound("H2S", "g", mp=187.66, bp=212.84)
@library
def so2():
    return compound("SO2", "g", mp=197.65, bp=263.13)
@library
def so3():
    return compound("SO3", "l", mp=289.95, bp=317.90)
@library
def hcn():
    return compound("HCN", "l", mp=259.86, bp=298.85)
@library
def cs2():
    return compound("CS2", "l", mp=161.58, bp=319.22)
@library
def ocs():
    return compound("OCS", "g", mp=134.3, bp=222.9)
@library
def h2o2():
    return compound("H2O2", "l", mp=272.72, bp=423.35)
@library
def nh3():
    return compound("NH3", "g", mp=195.40, bp=239.82)
@library
def ph3():
    return compound("PH3", "g", mp=139.4, bp=185.5)

# --- hydrides / C1–C2 organics ---
@library
def ch4():
    return compound("CH4", "g", mp=90.67, bp=111.66)
@library
def c2h6():
    return compound("C2H6", "g", mp=90.35, bp=184.55)
@library
def c2h4():
    return compound("C2H4", "g", mp=104.0, bp=169.4)
@library
def c2h2():
    return compound("C2H2", "g")
@library
def c3h8():
    return compound("C3H8", "g", mp=85.5, bp=231.1)
@library
def c3h6():
    return compound("C3H6", "g", mp=87.9, bp=225.5)
@library
def ch3oh():
    return compound("CH3OH", "l", mp=175.47, bp=337.85)
methanol = ch3oh
@library
def c2h5oh():
    return compound("C2H5OH", "l", mp=159.05, bp=351.52)
ethanol = c2h5oh
@library
def isopropanol():
    return compound("C3H7OH", "l", mp=184.65, bp=355.4)
@library
def hcho():
    return compound("HCHO", "g", mp=181.15, bp=254.05)
@library
def ch3cho():
    return compound("CH3CHO", "l", mp=150.15, bp=293.35)
@library
def hcooh():
    return compound("HCOOH", "l", mp=281.55, bp=373.95)
@library
def ch3cooh():
    return compound("CH3COOH", "l", mp=289.81, bp=391.2)
acetic_acid = ch3cooh
@library
def acetone():
    return compound("CH3COCH3", "l", mp=178.5, bp=329.22)

# --- common solvents ---
@library
def benzene():
    return compound("C6H6", "l", mp=278.68, bp=353.2)
@library
def toluene():
    return compound("C7H8", "l", mp=178.2, bp=383.8)
@library
def hexane():
    return compound("C6H14", "l", mp=177.8, bp=341.9)
@library
def chloroform():
    return compound("CHCl3", "l", mp=209.6, bp=334.3)
@library
def ccl4():
    return compound("CCl4", "l", mp=250.2, bp=349.8)
@library
def dichloromethane():
    return compound("CH2Cl2", "l", mp=176.5, bp=313.0)
@library
def diethyl_ether():
    return compound("C2H5OC2H5", "l", mp=156.85, bp=307.6)
@library
def ethyl_acetate():
    return compound("CH3COOC2H5", "l", mp=189.6, bp=350.3)
@library
def acetonitrile():
    return compound("CH3CN", "l", mp=227.4, bp=354.8)
@library
def thf():
    return compound("C4H8O", "l", mp=164.7, bp=339.0)
@library
def dmf():
    return compound("HCON(CH3)2", "l", mp=212.7, bp=426.0)
@library
def dmso():
    return compound("CH3SOCH3", "l", mp=291.7, bp=462.0)
@library
def phenol():
    return compound("C6H5OH", "s", mp=314.0, bp=454.9)
@library
def aniline():
    return compound("C6H5NH2", "l", mp=267.1, bp=457.3)
@library
def pyridine():
    return compound("C5H5N", "l", mp=231.6, bp=388.4)
@library
def n2h4():
    return compound("N2H4", "l", mp=275.15, bp=386.7)
@library
def urea():
    return compound("CO(NH2)2", "s", mp=406.0)
@library
def glucose():
    return compound("C6H12O6", "s", mp=419.0)
@library
def hocl():
    return compound("HOCl", "aq")
@library
def h2co3():
    return compound("H2CO3", "aq")
@library
def hno2():
    return compound("HNO2", "aq")
@library
def en():
    return compound("en", "aq")
@library
def phen():
    return compound("phen", "aq")
@library
def bf3():
    return compound("BF3", "g", mp=146.8, bp=173.2)
@library
def pcl3():
    return compound("PCl3", "l", mp=161.2, bp=349.3)
@library
def pcl5():
    return compound("PCl5", "s")
@library
def sicl4():
    return compound("SiCl4", "l", mp=204.3, bp=330.8)

# --- oxoacids ---
@library
def hno3():
    return compound("HNO3", "l", mp=231.55, bp=356.15)
nitric_acid = hno3
@library
def h2so4():
    return compound("H2SO4", "l", mp=283.46, bp=610.0)
sulfuric_acid = h2so4
@library
def h3po4():
    return compound("H3PO4", "s", mp=315.5)
phosphoric_acid = h3po4

# --- molecular elements and common salts / oxides ---
@library
def p4():
    return compound("P4", "s", mp=317.3, bp=553.7)
@library
def s8():
    return compound("S8", "s", mp=388.36, bp=717.8)
@library
def nacl():
    return compound("NaCl", "s", mp=1073.8, bp=1738.0)
@library
def kcl():
    return compound("KCl", "s", mp=1043.0, bp=1693.0)
@library
def nabr():
    return compound("NaBr", "s", mp=1020.0, bp=1663.0)
@library
def kbr():
    return compound("KBr", "s", mp=1007.0, bp=1656.0)
@library
def nai():
    return compound("NaI", "s", mp=933.0, bp=1577.0)
@library
def ki():
    return compound("KI", "s", mp=954.0, bp=1593.0)
@library
def naf():
    return compound("NaF", "s", mp=1266.0, bp=1978.0)
@library
def naoh():
    return compound("NaOH", "s", mp=596.0, bp=1661.0)
@library
def koh():
    return compound("KOH", "s", mp=633.0, bp=1600.0)
@library
def cao():
    return compound("CaO", "s", mp=2886.0, bp=3123.0)
@library
def mgo():
    return compound("MgO", "s", mp=3125.0, bp=3873.0)
@library
def al2o3():
    return compound("Al2O3", "s", mp=2345.0, bp=3250.0)
@library
def sio2():
    return compound("SiO2", "s", mp=1986.0, bp=3223.0)
@library
def fe2o3():
    return compound("Fe2O3", "s", mp=1812.0)
@library
def caco3():
    return compound("CaCO3", "s")
@library
def cacl2():
    return compound("CaCl2", "s", mp=1045.0, bp=1935.0)
@library
def mgcl2():
    return compound("MgCl2", "s", mp=987.0, bp=1685.0)
@library
def fecl2():
    return compound("FeCl2", "s", mp=950.0, bp=1296.0)
@library
def fecl3():
    return compound("FeCl3", "s", mp=577.0, bp=589.0)
@library
def cucl2():
    return compound("CuCl2", "s", mp=771.0, bp=1266.0)
@library
def zncl2():
    return compound("ZnCl2", "s", mp=563.0, bp=1005.0)
@library
def agno3():
    return compound("AgNO3", "s", mp=485.0)
@library
def kno3():
    return compound("KNO3", "s", mp=607.0)
@library
def na2so4():
    return compound("Na2SO4", "s", mp=1157.0)
@library
def baso4():
    return compound("BaSO4", "s", mp=1853.0)
@library
def cuso4():
    return compound("CuSO4", "s")
@library
def nh4cl():
    return compound("NH4Cl", "s")

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
    "en",
    "ethanol",
    "ethyl_acetate",
    "f2",
    "fe2o3",
    "fecl2",
    "fecl3",
    "glucose",
    "h2",
    "h2co3",
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
    "hno2",
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
    "phen",
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
