from ChemCompute import Enviroment
from ChemCompute.compounds import agcl, fescn, water
from ChemCompute.reactions import (
    acetic_acid,
    agcl_ksp,
    all_reactions,
    fescn_kf,
    hydrochloric_acid,
    water_kw,
)


def _all_zero(rxn):
    return all(species["concentration"] == 0.0 for species in rxn.compounds)


def test_library_reactions_start_at_zero_concentration():
    assert _all_zero(water_kw())
    assert _all_zero(acetic_acid())
    assert _all_zero(agcl_ksp())
    assert _all_zero(fescn_kf())


def test_water_kw_reuses_library_liquid_water():
    assert water_kw().K == 1.0e-14
    assert water_kw().reactants[0]["compound"] is water()
    assert water().phase(298.0) == "l"


def test_acetic_acid_stays_aqueous_in_q():
    ha = acetic_acid().reactants[0]["compound"]
    assert ha.formula == "CH3COOH"
    assert ha.phase(298.0) == "aq"
    assert acetic_acid().K == 1.8e-5


def test_agcl_ksp_uses_library_solid():
    assert agcl_ksp().reactants[0]["compound"] is agcl()
    assert agcl().phase(298.0) == "s"
    assert agcl_ksp().K == 1.8e-10


def test_fescn_uses_library_complex():
    assert fescn_kf().products[0]["compound"] is fescn()
    assert fescn_kf().K == 10 ** 3.22


def test_strong_acid_is_irreversible():
    assert hydrochloric_acid().infinite_K is True
    assert _all_zero(hydrochloric_acid())


def test_environment_overrides_do_not_mutate_library():
    env = Enviroment(water_kw(), acetic_acid(), concentrations={"CH3COOH": 0.10})
    assert env.concentrations_dict["CH3COOH"] == 0.10
    assert _all_zero(acetic_acid())
    assert _all_zero(water_kw())


def test_all_reactions_are_zero_concentration_objects():
    built = all_reactions()
    assert len(built) >= 40
    assert all(_all_zero(rxn) for rxn in built)


def test_repeated_calls_return_the_same_reaction():
    assert water_kw() is water_kw()
