"""Multi-reaction environments for manual validation (env16–env20)."""

from __future__ import annotations

from ChemCompute import Compound, Enviroment, Reaction

T = 298


def _aq(name: str, **kwargs) -> Compound:
    return Compound(name, phase_point_list=[{"phase": "aq", "temperature": T}], **kwargs)


def _solid(name: str, *, excess: bool = False) -> Compound:
    return Compound(name, phase_point_list=[{"phase": "s", "temperature": T}], excess=excess)


def _liquid(name: str, *, excess: bool = False) -> Compound:
    return Compound(name, phase_point_list=[{"phase": "l", "temperature": T}], excess=excess)


def build_env16() -> Enviroment:
    """HF fluoride speciation in 0.06 M HCl with excess CaF2(s) and H2O(l)."""
    hf = _aq("HF")
    fm = _aq("F-")
    hf2 = _aq("HF2-")
    h2f2 = _aq("H2F2")
    hp = _aq("H+")
    ohm = _aq("OH-")
    ca2 = _aq("Ca+2")
    hcl = _aq("HCl")
    clm = _aq("Cl-")
    h2o = _liquid("H2O", excess=True)
    caf2 = _solid("CaF2", excess=True)

    return Enviroment(
        Reaction(
            [
                {"stoichiometric_coefficient": 1, "compound": hf, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": fm, "rate_dependency": 1},
            ],
            [{"stoichiometric_coefficient": 1, "compound": hf2, "rate_dependency": 1}],
            [0.0, 0.0],
            [0.0],
            K=0.1,
        ),
        Reaction(
            [{"stoichiometric_coefficient": 2, "compound": hf, "rate_dependency": 2}],
            [{"stoichiometric_coefficient": 1, "compound": h2f2, "rate_dependency": 1}],
            [0.0],
            [0.0],
            K=0.5,
        ),
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": hf, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": fm, "rate_dependency": 1},
            ],
            [0.0],
            [0.0, 0.0],
            K=10 ** (-2.93),
        ),
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": h2o, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": ohm, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
            ],
            [55.5],
            [1e-14 / 0.06, 0.0],
            K=1e-14,
        ),
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": caf2, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": ca2, "rate_dependency": 1},
                {"stoichiometric_coefficient": 2, "compound": fm, "rate_dependency": 2},
            ],
            [10.0],
            [0.0, 0.0],
            K=5e-9,
        ),
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": hcl, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": clm, "rate_dependency": 1},
            ],
            [0.06],
            [0.0, 0.0],
            infinite_K=True,
        ),
        T=T,
    )


def build_env17() -> Enviroment:
    """Carbonate system with excess CaCO3(s), 0.01 M HCl (infinite K)."""
    caco3 = _solid("CaCO3", excess=True)
    ca2 = _aq("Ca+2")
    co3 = _aq("CO3-2")
    hco3 = _aq("HCO3-")
    hp = _aq("H+")
    ohm = _aq("OH-")
    h2o = _liquid("H2O", excess=True)
    hcl = _aq("HCl")
    clm = _aq("Cl-")

    return Enviroment(
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": caco3, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": ca2, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": co3, "rate_dependency": 1},
            ],
            [5.0],
            [0.0, 0.0],
            K=4.7e-9,
        ),
        Reaction(
            [
                {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": co3, "rate_dependency": 1},
            ],
            [{"stoichiometric_coefficient": 1, "compound": hco3, "rate_dependency": 1}],
            [0.0, 0.0],
            [0.0],
            K=1 / (4.7e-11),
        ),
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": h2o, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": ohm, "rate_dependency": 1},
            ],
            [55.5],
            [0.0, 1e-14 / 0.01],
            K=1e-14,
        ),
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": hcl, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": clm, "rate_dependency": 1},
            ],
            [0.01],
            [0.0, 0.0],
            infinite_K=True,
        ),
        T=T,
    )


def build_env18() -> Enviroment:
    """Ammonia buffer, AgCl(s) precipitation, and HCl dissociation (infinite K)."""
    nh4 = _aq("NH4+")
    hp = _aq("H+")
    nh3 = _aq("NH3")
    ag = _aq("Ag+")
    clm = _aq("Cl-")
    agcl = _solid("AgCl")
    hcl = _aq("HCl")

    return Enviroment(
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": nh4, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": nh3, "rate_dependency": 1},
            ],
            [0.1],
            [0.0, 0.01],
            K=5.6e-10,
        ),
        Reaction(
            [
                {"stoichiometric_coefficient": 1, "compound": ag, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": clm, "rate_dependency": 1},
            ],
            [{"stoichiometric_coefficient": 1, "compound": agcl, "rate_dependency": 1}],
            [0.01, 0.05],
            [0.0],
            infinite_K=True,
        ),
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": hcl, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": clm, "rate_dependency": 1},
            ],
            [0.05],
            [0.0, 0.0],
            infinite_K=True,
        ),
        T=T,
    )


def build_env19() -> Enviroment:
    """Phosphate buffer with CaHPO4(s) precipitation and excess solid."""
    h3po4 = _aq("H3PO4")
    h2po4 = _aq("H2PO4-")
    hpo4 = _aq("HPO4-2")
    ca2 = _aq("Ca+2")
    cahpo4 = _solid("CaHPO4", excess=True)
    hp = _aq("H+")
    ohm = _aq("OH-")
    h2o = _liquid("H2O", excess=True)

    return Enviroment(
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": h3po4, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": h2po4, "rate_dependency": 1},
            ],
            [0.01],
            [0.0, 0.0],
            K=7.5e-3,
        ),
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": h2po4, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": hpo4, "rate_dependency": 1},
            ],
            [0.0],
            [0.0, 0.0],
            K=6.2e-8,
        ),
        Reaction(
            [
                {"stoichiometric_coefficient": 1, "compound": ca2, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": hpo4, "rate_dependency": 1},
            ],
            [{"stoichiometric_coefficient": 1, "compound": cahpo4, "rate_dependency": 1}],
            [0.02, 0.0],
            [1.0],
            K=1e5,
        ),
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": h2o, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": ohm, "rate_dependency": 1},
            ],
            [55.5],
            [1e-7, 1e-7],
            K=1e-14,
        ),
        T=T,
    )


def build_env20() -> Enviroment:
    """Cu-EDTA complexation (infinite K), NaCl dissociation, water autoionization."""
    cu = _aq("Cu+2")
    edta = _aq("EDTA-4")
    cuedta = _aq("CuEDTA-2")
    ohm = _aq("OH-")
    nacl = _aq("NaCl")
    na = _aq("Na+")
    clm = _aq("Cl-")
    hp = _aq("H+")
    h2o = _liquid("H2O", excess=True)

    return Enviroment(
        Reaction(
            [
                {"stoichiometric_coefficient": 1, "compound": cu, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": edta, "rate_dependency": 1},
            ],
            [{"stoichiometric_coefficient": 1, "compound": cuedta, "rate_dependency": 1}],
            [0.005, 0.01],
            [0.0],
            infinite_K=True,
        ),
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": nacl, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": na, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": clm, "rate_dependency": 1},
            ],
            [0.1],
            [0.0, 0.0],
            infinite_K=True,
        ),
        Reaction(
            [{"stoichiometric_coefficient": 1, "compound": h2o, "rate_dependency": 1}],
            [
                {"stoichiometric_coefficient": 1, "compound": hp, "rate_dependency": 1},
                {"stoichiometric_coefficient": 1, "compound": ohm, "rate_dependency": 1},
            ],
            [55.5],
            [1e-7, 1e-7],
            K=1e-14,
        ),
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
