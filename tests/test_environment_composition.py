"""Tests for Environment Composition API."""

import math

import numpy as np

from ChemCompute import Compound, Enviroment, Reaction, ScaledEnviroment


def _simple_rxn(initial_a=1.0, initial_b=0.0, K=2.0):
    a = Compound("A", phase_point_list=[{"phase": "aq", "temperature": 298}])
    b = Compound("B", phase_point_list=[{"phase": "aq", "temperature": 298}])
    return Reaction(
        reactants=[{"stoichiometric_coefficient": 1, "compound": a, "rate_dependency": 1}],
        products=[{"stoichiometric_coefficient": 1, "compound": b, "rate_dependency": 1}],
        reactants_concentration=[initial_a],
        products_concentration=[initial_b],
        K=K,
    )


class TestFromCompounds:
    def test_compounds_only_no_reactions(self):
        env = Enviroment.from_compounds({"Na+": 0.1, "Cl-": 0.1}, volume=1.0)
        assert env.reactions == []
        assert np.isclose(env.concentrations_dict["Na+"], 0.1)
        assert np.isclose(env.concentrations_dict["Cl-"], 0.1)
        assert env.volume == 1.0

    def test_equilibrium_no_reactions(self):
        env = Enviroment.from_compounds({"A": 0.5}, volume=1.0)
        result = env.equilibrium(return_details=True)
        assert result.stop_reason == "no_reactions"
        assert result.concentrations == [0.5]


class TestConcentrationOverride:
    def test_override_beats_reaction_sum(self):
        rxn = _simple_rxn(initial_a=1.0, initial_b=0.0)
        env = Enviroment(rxn, concentrations={"A": 0.25})
        assert np.isclose(env.concentrations_dict["A"], 0.25)

    def test_override_adds_new_species(self):
        rxn = _simple_rxn()
        extra = Compound("X", phase_point_list=[{"phase": "aq", "temperature": 298}])
        env = Enviroment(rxn, concentrations={extra: 0.03})
        assert "X" in env.concentrations_dict
        assert np.isclose(env.concentrations_dict["X"], 0.03)


class TestStandaloneReaction:
    def test_equilibrium_matches_env_wrapper(self):
        rxn = _simple_rxn()
        direct = rxn.equilibrium(method="newton", tol=1e-10)
        wrapped = Enviroment(rxn).equilibrium(method="newton", tol=1e-10)
        assert np.allclose(direct, wrapped)

    def test_kinetics_runs(self):
        rxn = Reaction.from_string_simple_syntax(
            "A > B",
            [1.0, 0.0],
            K=2.0,
            kf=0.5,
            kb=0.25,
        )
        checkpoints = rxn.kinetics(time=1.0, accuracy=0.1, plot=False)
        assert checkpoints[-1][0] < 1.0


class TestEnvironmentMixing:
    def test_add_equal_volumes(self):
        envA = Enviroment.from_compounds({"A": 1.0}, volume=1.0)
        envB = Enviroment.from_compounds({"B": 0.1}, volume=1.0)
        envC = envA + envB
        assert envC.volume == 2.0
        assert np.isclose(envC.concentrations_dict["A"], 0.5)
        assert np.isclose(envC.concentrations_dict["B"], 0.05)

    def test_coefficient_mix(self):
        envA = Enviroment.from_compounds({"A": 1.0}, volume=1.0)
        envB = Enviroment.from_compounds({"B": 0.1}, volume=1.0)
        envD = 0.5 * envA + 4 * envB
        assert np.isclose(envD.volume, 4.5)
        assert np.isclose(envD.concentrations_dict["A"], 0.5 / 4.5)
        assert np.isclose(envD.concentrations_dict["B"], 0.4 / 4.5)

    def test_combine_matches_operator(self):
        envA = Enviroment.from_compounds({"A": 1.0}, volume=1.0)
        envB = Enviroment.from_compounds({"B": 0.1}, volume=1.0)
        env_op = 0.5 * envA + 4 * envB
        env_fn = Enviroment.combine((0.5, envA), (4.0, envB))
        assert np.isclose(env_op.volume, env_fn.volume)
        assert env_op.concentrations_dict == env_fn.concentrations_dict

    def test_add_compounds_slug(self):
        envB = Enviroment.from_compounds({"B": 0.1}, volume=1.0)
        envC = envB.add_compounds({"A": 1.0}, volume=1.0)
        assert envC.volume == 2.0
        assert np.isclose(envC.concentrations_dict["A"], 0.5)
        assert np.isclose(envC.concentrations_dict["B"], 0.05)

    def test_add_compounds_with_coefficient(self):
        envB = Enviroment.from_compounds({"B": 0.1}, volume=1.0)
        envE = envB.add_compounds({"A": 1.0}, volume=1.0, coefficient=4.0)
        assert np.isclose(envE.volume, 5.0)

    def test_temperature_mismatch_raises(self):
        envA = Enviroment.from_compounds({"A": 1.0}, T=298)
        envB = Enviroment.from_compounds({"B": 0.1}, T=310)
        try:
            envA + envB
            assert False, "expected ValueError"
        except ValueError:
            pass

    def test_invalid_coefficient_raises(self):
        envA = Enviroment.from_compounds({"A": 1.0})
        try:
            0 * envA + envA
            assert False, "expected ValueError"
        except ValueError:
            pass

    def test_scaled_type(self):
        envA = Enviroment.from_compounds({"A": 1.0})
        assert isinstance(0.5 * envA, ScaledEnviroment)
