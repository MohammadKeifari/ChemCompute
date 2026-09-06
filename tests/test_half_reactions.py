from ChemCompute import Enviroment, HalfReaction
from ChemCompute.compounds import fe2, fe3, h2, h_plus, i2, i_minus, water
from ChemCompute.environments import daniel_cell, iron_couple, water_limits
from ChemCompute.half_reactions import (
    all_half_reactions,
    hydrogen,
    iodine,
    iron_iii,
    oxygen,
    permanganate,
    she,
)


def _all_zero(hr):
    return all(
        value == 0.0
        for value in hr.oxidized_concentration + hr.reduced_concentration
    )


def test_library_half_reactions_start_at_zero():
    assert _all_zero(hydrogen())
    assert _all_zero(iron_iii())
    assert _all_zero(permanganate())
    assert _all_zero(iodine())


def test_repeated_calls_return_the_same_half_reaction():
    assert hydrogen() is hydrogen()
    assert iron_iii() is iron_iii()
    assert she is hydrogen
    assert she() is hydrogen()


def test_hydrogen_is_the_she():
    assert hydrogen().E0_SHE == 0.0
    assert hydrogen().n_electrons == 2.0
    assert hydrogen().oxidized[0]["compound"] is h_plus()
    assert hydrogen().reduced[0]["compound"] is h2()


def test_oxygen_reuses_library_water():
    assert oxygen().E0_SHE == 1.229
    assert oxygen().n_electrons == 4.0
    waters = [
        entry["compound"]
        for entry in oxygen().reduced
        if entry["compound"].formula == "H2O"
    ]
    assert waters[0] is water()
    assert water().phase(298.0) == "l"


def test_iron_iii_reuses_library_ions():
    assert iron_iii().E0_SHE == 0.771
    assert iron_iii().n_electrons == 1.0
    assert iron_iii().oxidized[0]["compound"] is fe3()
    assert iron_iii().reduced[0]["compound"] is fe2()


def test_permanganate_electron_count():
    assert permanganate().n_electrons == 5.0
    assert permanganate().E0_SHE == 1.507


def test_iodine_uses_library_solid():
    assert iodine().oxidized[0]["compound"] is i2()
    assert i2().phase(298.0) == "s"
    assert iodine().reduced[0]["compound"] is i_minus()


def test_environment_overrides_do_not_mutate_library():
    env = Enviroment(
        iron_iii(),
        concentrations={"Fe+3": 0.01, "Fe+2": 0.001},
    )
    assert env.concentrations_dict["Fe+3"] == 0.01
    assert _all_zero(iron_iii())
    assert env.half_reactions[0] is iron_iii()


def test_all_half_reactions_are_zero_concentration_objects():
    built = all_half_reactions()
    assert len(built) >= 40
    assert all(isinstance(hr, HalfReaction) for hr in built)
    assert all(_all_zero(hr) for hr in built)
    assert all(hr.n_electrons > 0 for hr in built)


def test_water_limits_environment():
    env = water_limits()
    assert env.half_reactions[0] is hydrogen()
    assert oxygen() in env.half_reactions
    assert {compound.formula for compound in env.compounds} >= {"H+", "OH-", "H2O", "H2", "O2"}
    assert all(value == 0.0 for value in env.concentrations)


def test_iron_environment_accepts_concentrations():
    env = iron_couple(concentrations={"Fe+3": 0.10})
    assert env.concentrations_dict["Fe+3"] == 0.10
    assert iron_couple().concentrations_dict["Fe+3"] == 0.0
    assert env.half_reactions[0] is iron_iii()


def test_daniel_cell_couples_copper_and_zinc():
    env = daniel_cell()
    assert {compound.formula for compound in env.compounds} >= {"Cu+2", "Cu", "Zn+2", "Zn"}
    assert len(env.half_reactions) == 2
