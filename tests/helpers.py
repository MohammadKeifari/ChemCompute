"""Shared builders for environment tests."""

from ChemCompute import Compound, Enviroment, Reaction


def aq(formula, **kwargs):
    return Compound(formula, phase_point_list=[{"phase": "aq", "temperature": 298}], **kwargs)


def simple_rxn(initial_a=1.0, initial_b=0.0, K=2.0):
    a = aq("A")
    b = aq("B")
    return Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": a, "rate_dependency": 1}],
        products=[{"stoichiometric_coefficient": 1, "compound": b, "rate_dependency": 1}],
        reactants_concentration=[initial_a],
        products_concentration=[initial_b],
        K=K,
    )


def water_env(h_plus=1e-7, oh_minus=1e-7, volume=0.1):
    h = aq("H+", charge=1)
    oh = aq("OH-", charge=-1)
    water = Compound("H2O", excess=True, phase_point_list=[{"phase": "l", "temperature": 298}])
    rxn = Reaction(
        reactants=[
            {"stoichiometric_coefficient": 1, "compound": h, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": oh, "rate_dependency": 1},
        ],
        products=[{"stoichiometric_coefficient": 1, "compound": water, "rate_dependency": 0}],
        reactants_concentration=[h_plus, oh_minus],
        products_concentration=[0.0],
        K=1e-14,
    )
    return Enviroment(
        rxn,
        concentrations={"H+": h_plus, "OH-": oh_minus},
        volume=volume,
    )


def weak_acid_env(*, h_plus=1e-7, ha=0.1, a_minus=0.0, buffer=None):
    ha_c = aq("HA")
    h_c = aq("H+", charge=1)
    a_c = aq("A-", charge=-1)
    rxn = Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": ha_c, "rate_dependency": 1}],
        products=[
            {"stoichiometric_coefficient": 1, "compound": h_c, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": a_c, "rate_dependency": 1},
        ],
        reactants_concentration=[ha],
        products_concentration=[h_plus, a_minus],
        K=1e-5,
        kf=1.0,
        kb=1e5,
    )
    return Enviroment(
        rxn,
        concentrations={"H+": h_plus, "HA": ha, "A-": a_minus},
        buffer=buffer,
    )


def ammonia_buffer_env():
    nh4 = aq("NH4+", charge=1)
    nh3 = aq("NH3")
    h = aq("H+", charge=1)
    oh = aq("OH-", charge=-1)
    water = Compound("H2O", excess=True)
    ka = 5.6e-10
    rxn1 = Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": nh4, "rate_dependency": 1}],
        products=[
            {"stoichiometric_coefficient": 1, "compound": h, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": nh3, "rate_dependency": 1},
        ],
        reactants_concentration=[0.05],
        products_concentration=[1e-9, 0.05],
        K=ka,
    )
    rxn2 = Reaction(
        reactants=[
            {"stoichiometric_coefficient": 1, "compound": h, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": oh, "rate_dependency": 1},
        ],
        products=[{"stoichiometric_coefficient": 1, "compound": water, "rate_dependency": 0}],
        reactants_concentration=[1e-9, 1e-7],
        products_concentration=[0.0],
        K=1e-14,
    )
    return Enviroment(rxn1, rxn2)
