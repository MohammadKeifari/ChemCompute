"""Premade enzyme-kinetics Enviroment factories."""

from __future__ import annotations

from .._general import Compound, Enviroment, Reaction


def _mm_reaction(substrate, product, Vmax, Km, rate_law="michaelis_menten", **extra):
    s = Compound(substrate, phase_point_list=[{"temperature": 298, "phase": "aq"}])
    p = Compound(product, phase_point_list=[{"temperature": 298, "phase": "aq"}])
    rxn = Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": s, "rate_dependency": 1}],
        products=[{"stoichiometric_coefficient": 1, "compound": p, "rate_dependency": 1}],
        reactants_concentration=[1.0],
        products_concentration=[0.0],
        kf=0.0,
        kb=0.0,
    )
    rxn.rate_law = rate_law
    rxn.bio_params = {"substrate": substrate, "Vmax": Vmax, "Km": Km, **extra}
    return rxn


def single_substrate_mm(
    substrate="S",
    product="P",
    *,
    s0=1.0,
    p0=0.0,
    Vmax=1e-6,
    Km=1e-4,
    T=298,
):
    """Michaelis-Menten: v = Vmax[S]/(Km+[S])."""
    rxn = _mm_reaction(substrate, product, Vmax, Km)
    rxn.compounds[0]["concentration"] = s0
    rxn.compounds[1]["concentration"] = p0
    return Enviroment(rxn, T=T)


def competitive_inhibition(
    substrate="S",
    product="P",
    inhibitor="I",
    *,
    s0=1.0,
    p0=0.0,
    i0=0.1,
    Vmax=1e-6,
    Km=1e-4,
    Ki=1e-4,
    T=298,
):
    """Competitive inhibition: apparent Km increases."""
    i = Compound(inhibitor, phase_point_list=[{"temperature": 298, "phase": "aq"}])
    rxn = _mm_reaction(
        substrate,
        product,
        Vmax,
        Km,
        rate_law="mm_competitive",
        inhibitor=inhibitor,
        Ki=Ki,
    )
    rxn.compounds[0]["concentration"] = s0
    rxn.compounds[1]["concentration"] = p0
    env = Enviroment(rxn, T=T)
    env.compounds_concentration.append({"compound": i, "concentration": i0})
    env.compounds.append(i)
    return env


def uncompetitive_inhibition(
    substrate="S",
    product="P",
    inhibitor="I",
    *,
    s0=1.0,
    p0=0.0,
    i0=0.1,
    Vmax=1e-6,
    Km=1e-4,
    Ki=1e-4,
    T=298,
):
    """Uncompetitive inhibition."""
    i = Compound(inhibitor, phase_point_list=[{"temperature": 298, "phase": "aq"}])
    rxn = _mm_reaction(
        substrate,
        product,
        Vmax,
        Km,
        rate_law="mm_uncompetitive",
        inhibitor=inhibitor,
        Ki=Ki,
    )
    rxn.compounds[0]["concentration"] = s0
    rxn.compounds[1]["concentration"] = p0
    env = Enviroment(rxn, T=T)
    env.compounds_concentration.append({"compound": i, "concentration": i0})
    env.compounds.append(i)
    return env


def noncompetitive_inhibition(
    substrate="S",
    product="P",
    inhibitor="I",
    *,
    s0=1.0,
    p0=0.0,
    i0=0.1,
    Vmax=1e-6,
    Km=1e-4,
    Ki=1e-4,
    T=298,
):
    """Noncompetitive inhibition."""
    i = Compound(inhibitor, phase_point_list=[{"temperature": 298, "phase": "aq"}])
    rxn = _mm_reaction(
        substrate,
        product,
        Vmax,
        Km,
        rate_law="mm_noncompetitive",
        inhibitor=inhibitor,
        Ki=Ki,
    )
    rxn.compounds[0]["concentration"] = s0
    rxn.compounds[1]["concentration"] = p0
    env = Enviroment(rxn, T=T)
    env.compounds_concentration.append({"compound": i, "concentration": i0})
    env.compounds.append(i)
    return env


def mixed_inhibition(
    substrate="S",
    product="P",
    inhibitor="I",
    *,
    s0=1.0,
    p0=0.0,
    i0=0.1,
    Vmax=1e-6,
    Km=1e-4,
    Ki=1e-4,
    alpha=1.0,
    alpha_prime=1.0,
    T=298,
):
    """Mixed inhibition with alpha and alpha_prime factors."""
    i = Compound(inhibitor, phase_point_list=[{"temperature": 298, "phase": "aq"}])
    rxn = _mm_reaction(
        substrate,
        product,
        Vmax,
        Km,
        rate_law="mm_mixed",
        inhibitor=inhibitor,
        Ki=Ki,
        alpha=alpha,
        alpha_prime=alpha_prime,
    )
    rxn.compounds[0]["concentration"] = s0
    rxn.compounds[1]["concentration"] = p0
    env = Enviroment(rxn, T=T)
    env.compounds_concentration.append({"compound": i, "concentration": i0})
    env.compounds.append(i)
    return env


def sequential_pathway(
    *,
    s0=1.0,
    a0=0.0,
    b0=0.0,
    p0=0.0,
    Vmax1=1e-6,
    Km1=1e-4,
    Vmax2=8e-7,
    Km2=1e-4,
    T=298,
):
    """Sequential S -> A -> B -> P with two MM steps."""
    s = Compound("S", phase_point_list=[{"temperature": 298, "phase": "aq"}])
    a = Compound("A", phase_point_list=[{"temperature": 298, "phase": "aq"}])
    b = Compound("B", phase_point_list=[{"temperature": 298, "phase": "aq"}])
    p = Compound("P", phase_point_list=[{"temperature": 298, "phase": "aq"}])

    rxn1 = Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": s, "rate_dependency": 1}],
        products=[{"stoichiometric_coefficient": 1, "compound": a, "rate_dependency": 1}],
        reactants_concentration=[s0],
        products_concentration=[a0],
        kf=0.0,
        kb=0.0,
    )
    rxn1.rate_law = "michaelis_menten"
    rxn1.bio_params = {"substrate": "S", "Vmax": Vmax1, "Km": Km1}

    rxn2 = Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": a, "rate_dependency": 1}],
        products=[{"stoichiometric_coefficient": 1, "compound": b, "rate_dependency": 1}],
        reactants_concentration=[0.0],
        products_concentration=[b0],
        kf=0.0,
        kb=0.0,
    )
    rxn2.rate_law = "michaelis_menten"
    rxn2.bio_params = {"substrate": "A", "Vmax": Vmax2, "Km": Km2}

    rxn3 = Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": b, "rate_dependency": 1}],
        products=[{"stoichiometric_coefficient": 1, "compound": p, "rate_dependency": 1}],
        reactants_concentration=[0.0],
        products_concentration=[p0],
        kf=0.0,
        kb=0.0,
    )
    rxn3.rate_law = "michaelis_menten"
    rxn3.bio_params = {"substrate": "B", "Vmax": Vmax2, "Km": Km2}

    return Enviroment(rxn1, rxn2, rxn3, T=T)
