import numpy as np
import pytest

from ChemCompute import Compound, Enviroment, Reaction


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


def test_env_equilibrium_concentration_error_limit(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    reference = env.equilibrium(method="bgd", tol=1e-8, max_iter=5000)
    result = env.equilibrium(
        method="bgd",
        tol=1.0,
        concentration_error_limit=1e-8,
        max_iter=5000,
    )
    assert len(result) == len(reference)
    np.testing.assert_allclose(result, reference, rtol=1e-4, atol=1e-6)


def test_env_equilibrium_concentration_error_limit_ignores_tol(simple_equilibrium_environment):
    env = simple_equilibrium_environment
    reference = env.equilibrium(method="bgd", tol=1e-8, max_iter=5000)
    result = env.equilibrium(
        method="bgd",
        tol=1e30,
        concentration_error_limit=1e-8,
        max_iter=5000,
    )
    np.testing.assert_allclose(result, reference, rtol=1e-4, atol=1e-6)


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
