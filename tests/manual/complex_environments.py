"""Multi-reaction environments for manual validation (env16–env20)."""

from __future__ import annotations

from ChemCompute import Enviroment, Reaction, XS

T = 298


def _rxn(reaction_str: str, concentrations: dict, **kwargs) -> Reaction:
    return Reaction.from_string(reaction_str, concentrations=concentrations, T=T, **kwargs)


def build_env16() -> Enviroment:
    """HF fluoride speciation in 0.06 M HCl with excess CaF2(s) and H2O(l)."""
    return Enviroment(
        _rxn("HF.aq & F-.aq > HF2-.aq", {}, K=0.1),
        _rxn("2_HF.aq > H2F2.aq", {}, K=0.5),
        _rxn("HF.aq > H+ & F-", {}, K=10 ** (-2.93)),
        _rxn(
            "H2O.l > H+ & OH-",
            {"H2O": XS(55.5), "H+": 1e-14 / 0.06},
            K=1e-14,
        ),
        _rxn("CaF2.s > Ca+2 & 2_F-", {"CaF2": XS(10.0)}, K=5e-9),
        _rxn("HCl.aq > H+ & Cl-", {"HCl": 0.06}, K=1.0, infinite_K=True),
        T=T,
    )


def build_env17() -> Enviroment:
    """Carbonate system with excess CaCO3(s), 0.01 M HCl (infinite K)."""
    return Enviroment(
        _rxn("CaCO3.s > Ca+2 & CO3-2", {"CaCO3": XS(5.0)}, K=4.7e-9),
        _rxn("H+ & CO3-2 > HCO3-", {}, K=1 / (4.7e-11)),
        _rxn(
            "H2O.l > H+ & OH-",
            {"H2O": XS(55.5), "OH-": 1e-14 / 0.01},
            K=1e-14,
        ),
        _rxn("HCl.aq > H+ & Cl-", {"HCl": 0.01}, K=1.0, infinite_K=True),
        T=T,
    )


def build_env18() -> Enviroment:
    """Ammonia buffer, AgCl(s) precipitation, and HCl dissociation (infinite K)."""
    return Enviroment(
        _rxn("NH4+ > H+ & NH3", {"NH4+": 0.1, "NH3": 0.01}, K=5.6e-10),
        _rxn("Ag+ & Cl- > AgCl.s", {"Ag+": 0.01, "Cl-": 0.05}, K=1.0, infinite_K=True),
        _rxn("HCl.aq > H+ & Cl-", {"HCl": 0.05}, K=1.0, infinite_K=True),
        T=T,
    )


def build_env19() -> Enviroment:
    """Phosphate buffer with CaHPO4(s) precipitation and excess solid."""
    return Enviroment(
        _rxn("H3PO4 > H+ & H2PO4-", {"H3PO4": 0.01}, K=7.5e-3),
        _rxn("H2PO4- > H+ & HPO4-2", {}, K=6.2e-8),
        _rxn(
            "Ca+2 & HPO4-2 > CaHPO4.s",
            {"Ca+2": 0.02, "CaHPO4": XS(1.0)},
            K=1e5,
        ),
        _rxn("H2O.l > H+ & OH-", {"H2O": XS(55.5), "H+": 1e-7, "OH-": 1e-7}, K=1e-14),
        T=T,
    )


def build_env20() -> Enviroment:
    """Cu-EDTA complexation (infinite K), NaCl dissociation, water autoionization."""
    return Enviroment(
        _rxn(
            "Cu+2 & EDTA-4 > CuEDTA-2",
            {"Cu+2": 0.005, "EDTA-4": 0.01},
            K=1.0,
            infinite_K=True,
        ),
        _rxn("NaCl.aq > Na+ & Cl-", {"NaCl": 0.1}, K=1.0, infinite_K=True),
        _rxn("H2O.l > H+ & OH-", {"H2O": XS(55.5), "H+": 1e-7, "OH-": 1e-7}, K=1e-14),
        T=T,
    )


# Reference concentrations from `python tests/manual/generate_expected.py`
# (Newton, tol=1e-10, min_concentration=1e-20, max_iter=8000).
# Water: Kw = [H+][OH-] = 1e-14. Excess s/l phases: activity = 1, omitted from Q and K.
ENV16_EXPECTED = [
    0.021509986486924065,
    0.000664600379739335,
    1.4295545187397702e-06,
    0.0002313397593338279,
    0.03802590443998585,
    55.5,
    2.629786233216643e-13,
    10.0,
    0.011320062747184268,
    0.0,
    0.06,
]
ENV17_EXPECTED = [
    5.0,
    0.0099994801759049,
    4.700244330031744e-07,
    9.99849037952633e-07,
    0.009999010151471897,
    55.5,
    1.0001509849098483e-08,
    0.0,
    0.01,
]
ENV18_EXPECTED = [
    0.10999999846000008,
    0.04000000153999992,
    1.5399999192533542e-09,
    0.0,
    0.09,
    0.01,
    0.0,
]
ENV19_EXPECTED = [
    1.0000727707418222e-05,
    3.754732296543353e-06,
    0.019976246475594338,
    0.00032985767923505126,
    0.03031610488253681,
    1.0,
    55.5,
    2.663305719343177e-09,
]
ENV20_EXPECTED = [0.0, 0.005, 0.005, 0.0, 0.1, 0.1, 55.5, 1e-07, 1e-07]

REFERENCE_EQUILIBRIUM_KWARGS = {
    "method": "newton",
    "loss": "log_quotient",
    "tol": 1e-10,
    "max_iter": 8000,
    "min_concentration": 1e-20,
}

COMPLEX_EQUILIBRIUM_KWARGS = {name: dict(REFERENCE_EQUILIBRIUM_KWARGS) for name in (
    "env16", "env17", "env18", "env19", "env20"
)}
