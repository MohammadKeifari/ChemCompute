"""Sparingly soluble salts used in precipitation (Ksp) reactions. Phase is solid."""

from ._builders import compound, library

@library
def agcl():
    return compound("AgCl", "s", mp=728.0, bp=1823.0)
@library
def agbr():
    return compound("AgBr", "s", mp=705.0)
@library
def agi():
    return compound("AgI", "s", mp=831.0)
@library
def ag2cro4():
    return compound("Ag2CrO4", "s")
@library
def ag2s():
    return compound("Ag2S", "s", mp=1098.0)
@library
def ag2so4():
    return compound("Ag2SO4", "s")
@library
def baco3():
    return compound("BaCO3", "s")
@library
def baf2():
    return compound("BaF2", "s", mp=1641.0)
@library
def caf2():
    return compound("CaF2", "s", mp=1691.0)
@library
def caso4():
    return compound("CaSO4", "s")
@library
def caoh2():
    return compound("Ca(OH)2", "s")
@library
def mgoh2():
    return compound("Mg(OH)2", "s")
@library
def mgco3():
    return compound("MgCO3", "s")
@library
def feoh3():
    return compound("Fe(OH)3", "s")
@library
def feoh2():
    return compound("Fe(OH)2", "s")
@library
def aloh3():
    return compound("Al(OH)3", "s")
@library
def znoh2():
    return compound("Zn(OH)2", "s")
@library
def cuoh2():
    return compound("Cu(OH)2", "s")
@library
def nioh2():
    return compound("Ni(OH)2", "s")
@library
def mnoh2():
    return compound("Mn(OH)2", "s")
@library
def pbs():
    return compound("PbS", "s", mp=1391.0)
@library
def pbi2():
    return compound("PbI2", "s")
@library
def pbcl2():
    return compound("PbCl2", "s", mp=774.0)
@library
def pbso4():
    return compound("PbSO4", "s")
@library
def cus():
    return compound("CuS", "s")
@library
def zns():
    return compound("ZnS", "s")
@library
def cds():
    return compound("CdS", "s")
@library
def fes():
    return compound("FeS", "s")
@library
def hgs():
    return compound("HgS", "s")
@library
def srf2():
    return compound("SrF2", "s", mp=1750.0)
@library
def srso4():
    return compound("SrSO4", "s")

__all__ = [
    "ag2cro4",
    "ag2s",
    "ag2so4",
    "agbr",
    "agcl",
    "agi",
    "aloh3",
    "baco3",
    "baf2",
    "caf2",
    "caoh2",
    "caso4",
    "cds",
    "cuoh2",
    "cus",
    "feoh2",
    "feoh3",
    "fes",
    "hgs",
    "mgco3",
    "mgoh2",
    "mnoh2",
    "nioh2",
    "pbcl2",
    "pbi2",
    "pbs",
    "pbso4",
    "srf2",
    "srso4",
    "znoh2",
    "zns",
]
