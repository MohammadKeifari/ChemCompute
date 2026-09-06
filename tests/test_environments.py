from ChemCompute.compounds import agcl, agcl2, water
from ChemCompute.environments import (
    all_environments,
    phosphoric_acid,
    silver_chloride,
    silver_chloride_ammonia,
    water_limits,
)
from ChemCompute.reactions import water_kw


def _formulas(env):
    return {compound.formula for compound in env.compounds}


def _all_zero(env):
    return all(value == 0.0 for value in env.concentrations)


def test_phosphoric_acid_couples_all_deprotonations_and_kw():
    env = phosphoric_acid()
    assert len(env.reactions) == 4
    assert _formulas(env) >= {"H3PO4", "H2PO4-", "HPO4-2", "PO4-3", "H+", "OH-", "H2O"}
    assert _all_zero(env)
    assert any(rxn.reactants[0]["compound"] is water() for rxn in env.reactions)


def test_phosphoric_acid_accepts_concentrations():
    env = phosphoric_acid(concentrations={"H3PO4": 0.10})
    assert env.concentrations_dict["H3PO4"] == 0.10
    assert phosphoric_acid().concentrations_dict["H3PO4"] == 0.0


def test_each_call_is_an_independent_environment():
    first = phosphoric_acid(concentrations={"H3PO4": 0.10})
    second = phosphoric_acid()
    assert first is not second
    assert first.reactions[0] is water_kw()
    assert second.reactions[0] is water_kw()
    assert second.concentrations_dict["H3PO4"] == 0.0


def test_silver_chloride_couples_ksp_and_chloro_complex():
    env = silver_chloride()
    assert agcl() in env.compounds
    assert agcl2() in env.compounds
    assert _formulas(env) >= {"AgCl", "Ag+", "Cl-", "AgCl2-"}
    assert _all_zero(env)


def test_silver_chloride_ammonia_has_ammine_and_solid():
    env = silver_chloride_ammonia()
    assert _formulas(env) >= {"AgCl", "Ag+", "Cl-", "Ag(NH3)2+", "NH4+", "NH3"}


def test_all_environments_build_at_zero():
    built = all_environments()
    assert len(built) >= 15
    assert all(_all_zero(env) for env in built)
    assert all(len(env.reactions) >= 2 for env in built)


def test_water_limits_uses_library_half_reactions():
    env = water_limits()
    assert len(env.half_reactions) == 2
    assert _all_zero(env)

