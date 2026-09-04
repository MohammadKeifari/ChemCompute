import numpy as np
import pytest

from ChemCompute import Compound, Enviroment, EquilibriumResult, Reaction


@pytest.fixture
def simple_equilibrium_environment():
    A = Compound("A")
    B = Compound("B")
    reactants = [{"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1}]
    products = [{"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}]
    rxn = Reaction(reactants, products, [1.0], [0.0], K=2.0, kf=0.5, kb=0.25)
    return Enviroment(rxn, T=298)


@pytest.fixture
def simple_kinetic_environment():
    rxn = Reaction.from_string_simple_syntax(
        "A > B",
        concentrations=[1.0, 0.0],
        K=2.0,
        kf=0.5,
        kb=0.25,
    )
    return Enviroment(rxn, T=298)


def test_env_equilibrium_default(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(max_iter=1000, tol=1e-8)
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


@pytest.mark.parametrize("method", ["bgd", "sgd", "newton"])
def test_env_equilibrium_methods(simple_equilibrium_environment, method):
    env = simple_equilibrium_environment
    result = env.equilibrium(method=method, max_iter=1000, tol=1e-8)
    assert len(result) == 2
    assert np.isclose(result[0] + result[1], 1.0, atol=0.01)


@pytest.mark.parametrize("loss", ["log_quotient", "quotient_error", "log_huber"])
def test_env_equilibrium_loss_functions(simple_equilibrium_environment, loss):
    env = simple_equilibrium_environment
    result = env.equilibrium(method="bgd", loss=loss, max_iter=1000, tol=1e-6)
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


def test_env_equilibrium_invalid_method(simple_equilibrium_environment):
    with pytest.raises(ValueError, match="Invalid method"):
        simple_equilibrium_environment.equilibrium(method="invalid")


def test_env_equilibrium_invalid_loss(simple_equilibrium_environment):
    with pytest.raises(ValueError, match="Invalid loss"):
        simple_equilibrium_environment.equilibrium(loss="invalid")


def test_compound_labels(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    assert env.compound_labels == ["A", "B"]
    assert env.compound_labels == [c.formula for c in env.compounds]


def test_equilibrium_return_details(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(method="newton", tol=1e-10, return_details=True)

    assert isinstance(result, EquilibriumResult)
    assert result.compounds == ["A", "B"]
    assert len(result.concentrations) == 2
    assert len(result.reaction_extents) == 1
    assert len(result.reaction_extent_percent) == 1
    assert len(result.reaction_quotient_ratio) == 1
    assert result.stop_reason in {"residual_tol", "reaction_extent_limit", "max_iter"}
    assert result.iterations >= 1


def test_reaction_extent_percent_near_equilibrium(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(method="newton", tol=1e-12, return_details=True)

    assert result.max_reaction_extent_percent < 1e-4
    assert np.allclose(result.reaction_quotient_ratio, [1.0], rtol=1e-3, atol=1e-3)


def test_reaction_extent_error_limit_stopping(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(
        method="bgd",
        reaction_extent_error_limit=0.01,
        max_iter=5000,
        return_details=True,
    )

    assert result.stop_reason == "reaction_extent_limit"
    assert result.max_reaction_extent_percent <= 0.01
    assert result.converged is True


def test_reaction_extent_error_limit_ignores_tol(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    loose_tol_result = env.equilibrium(
        method="bgd",
        tol=1e30,
        reaction_extent_error_limit=0.01,
        max_iter=5000,
        return_details=True,
    )
    assert loose_tol_result.stop_reason == "reaction_extent_limit"
    assert loose_tol_result.max_reaction_extent_percent <= 0.01


def test_equilibrium_criterion_met_with_tol(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(method="newton", tol=1e-10, return_details=True)

    assert result.criterion_type == "residual_tol"
    assert result.criterion_met is True
    assert result.criterion_value < result.criterion_limit


def test_equilibrium_criterion_met_with_extent_limit(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(
        method="bgd",
        reaction_extent_error_limit=0.01,
        max_iter=5000,
        return_details=True,
    )

    assert result.criterion_type == "reaction_extent"
    assert result.criterion_met is True
    assert result.criterion_value <= result.criterion_limit


def test_equilibrium_criterion_not_met_at_max_iter(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    result = env.equilibrium(method="bgd", tol=1e-30, max_iter=1, return_details=True)

    assert result.stop_reason == "max_iter"
    assert result.criterion_type == "residual_tol"
    assert result.criterion_met is False


def test_equilibrium_last_result_without_return_details(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    concentrations = env.equilibrium(method="newton", tol=1e-10)

    assert isinstance(concentrations, list)
    last = env.last_equilibrium_result
    assert isinstance(last, EquilibriumResult)
    assert last.concentrations == concentrations
    assert last.stop_reason in {"residual_tol", "reaction_extent_limit", "max_iter"}
    assert last.iterations >= 1
    assert len(last.q_over_k) == 1
    assert np.isclose(last.q_over_k[0], 1.0, rtol=1e-3, atol=1e-3)


def test_apply_equilibrium_updates_concentrations(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    assert env.concentrations == [1.0, 0.0]

    result = env.apply_equilibrium(method="newton", tol=1e-10)

    assert isinstance(result, EquilibriumResult)
    assert env.concentrations == result.concentrations
    assert np.isclose(env.concentrations[0], 1.0 / 3.0, rtol=1e-4)
    assert np.isclose(env.concentrations[1], 2.0 / 3.0, rtol=1e-4)
    assert env.last_equilibrium_result is result
    assert result.iterations >= 1


def test_env_kinetics_basic(simple_kinetic_environment):
    env = simple_kinetic_environment
    result = env.kinetics(time=1.0, accuracy=1e-2, plot=False)
    assert len(result) >= 1
    final = result[-1]
    assert len(final) == 2
    assert all(c >= 0 for c in final)


def test_env_kinetics_mass_conservation(simple_kinetic_environment):
    env = simple_kinetic_environment
    initial_total = sum(env.concentrations)
    result = env.kinetics(time=2.0, accuracy=1e-2, plot=False)
    final_total = sum(result[-1])
    assert np.isclose(final_total, initial_total, atol=0.05)


def test_env_equilibrium_phase_exclusion():
    A_g = Compound("A", phase_point_list=[{"phase": "g", "temperature": 298}])
    B_l = Compound("B", phase_point_list=[{"phase": "l", "temperature": 298}])
    C_s = Compound("C", phase_point_list=[{"phase": "s", "temperature": 298}])
    rxn = Reaction(
        [{"stoichiometric_coefficient": 1, "compound": A_g, "rate_dependency": 1}],
        [
            {"stoichiometric_coefficient": 1, "compound": B_l, "rate_dependency": 1},
            {"stoichiometric_coefficient": 1, "compound": C_s, "rate_dependency": 1},
        ],
        [1.0],
        [0.0, 0.0],
        K=5.0,
        kf=0.3,
        kb=0.06,
    )
    env = Enviroment(rxn, T=298)
    result = env.equilibrium(max_iter=1000, tol=1e-8)
    assert len(result) == 3
    assert np.isclose(result[0], 0.2, atol=0.02)
