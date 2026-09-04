import warnings

import numpy as np
import pytest

from ChemCompute import Compound, Enviroment, Reaction
from ChemCompute.Kinetic import KineticalCalculator
from ChemCompute.Thermodynamic import EquilibriumCalculator


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
def test_env_equilibrium_methods_match_calculator(simple_equilibrium_environment, method):
    env = simple_equilibrium_environment

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        calc = EquilibriumCalculator(method_of_calculation=method)
        calc.fit(env)
        expected = calc.calculate(max_iter=1000, tol=1e-8)

    actual = env.equilibrium(method=method, max_iter=1000, tol=1e-8)
    assert len(actual) == len(expected)
    np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=1e-8)


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
    loose = env.equilibrium(
        method="bgd",
        max_iter=5000,
        tol=1e-20,
        concentration_error_limit=1.0,
    )
    tight = env.equilibrium(
        method="bgd",
        max_iter=5000,
        tol=1e-12,
        concentration_error_limit=None,
    )
    assert isinstance(loose, list)
    assert len(loose) == len(tight)


def test_env_kinetics_matches_calculator(simple_kinetic_environment):
    env = simple_kinetic_environment

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        kc = KineticalCalculator(accuracy=1e-2)
        kc.fit(env)
        expected = kc.calculate(time=1.0, plot=False)

    actual = env.kinetics(time=1.0, accuracy=1e-2, plot=False)
    assert len(actual) == len(expected)
    np.testing.assert_allclose(actual[-1], expected[-1], rtol=1e-5, atol=1e-8)


def test_equilibrium_calculator_deprecation_warning():
    with pytest.warns(DeprecationWarning, match="EquilibriumCalculator is deprecated"):
        EquilibriumCalculator()


def test_kinetical_calculator_deprecation_warning():
    with pytest.warns(DeprecationWarning, match="KineticalCalculator is deprecated"):
        KineticalCalculator()
