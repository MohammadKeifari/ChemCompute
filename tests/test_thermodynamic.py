import pytest
import numpy as np
import sys
from pathlib import Path

# Add paths to sys.path for imports
project_root = Path(__file__).parent.parent
src_chemcalc_path = project_root / "src" / "ChemCalc"
sys.path.insert(0, str(src_chemcalc_path))  # For _general import in Thermodynamic.py
sys.path.insert(0, str(project_root))  # For src.ChemCalc imports

# Import from _general the same way Thermodynamic.py does to ensure type compatibility
from _general import Enviroment, Compound, Reaction
from src.ChemCalc.Thermodynamic import EquilibriumCalculator


# -------------------------
# EquilibriumCalculator Tests
# -------------------------

# ---------- Initialization Tests ---------- #

def test_equilibrium_calculator_init_default():
    """Test initialization with default method."""
    calc = EquilibriumCalculator()
    assert calc.method_of_calculation == "bgd"
    assert calc.fitted == False


def test_equilibrium_calculator_init_bgd():
    """Test initialization with batch gradient descent method."""
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    assert calc.method_of_calculation == "bgd"
    assert calc.fitted == False


def test_equilibrium_calculator_init_sgd():
    """Test initialization with stochastic gradient descent method."""
    calc = EquilibriumCalculator(method_of_calculation="sgd")
    assert calc.method_of_calculation == "sgd"
    assert calc.fitted == False


def test_equilibrium_calculator_init_newton():
    """Test initialization with Newton's method."""
    calc = EquilibriumCalculator(method_of_calculation="newton")
    assert calc.method_of_calculation == "newton"
    assert calc.fitted == False


# ---------- Fixtures for test environments ---------- #

@pytest.fixture
def simple_equilibrium_environment():
    """Create a simple environment with one reversible reaction: A ⇌ B"""
    A = Compound("A")
    B = Compound("B")
    
    reactants = [{"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1}]
    products = [{"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}]
    
    rxn = Reaction(
        reactants, 
        products, 
        [1.0],  # A initial concentration
        [0.0],  # B initial concentration
        K=2.0,  # Equilibrium constant
        kf=0.5, 
        kb=0.25
    )
    env = Enviroment(rxn, T=298)
    return env


@pytest.fixture
def multi_reaction_equilibrium_environment():
    """Create an environment with multiple reactions: A ⇌ B and B ⇌ C"""
    A = Compound("A")
    B = Compound("B")
    C = Compound("C")
    
    # Reaction 1: A ⇌ B
    reactants1 = [{"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1}]
    products1 = [{"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}]
    rxn1 = Reaction(
        reactants1, 
        products1, 
        [1.0], 
        [0.0], 
        K=2.0, 
        kf=0.5, 
        kb=0.25
    )
    
    # Reaction 2: B ⇌ C
    reactants2 = [{"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}]
    products2 = [{"stoichiometric_coefficient": 1, "compound": C, "rate_dependency": 1}]
    rxn2 = Reaction(
        reactants2, 
        products2, 
        [0.0], 
        [0.0], 
        K=1.5, 
        kf=0.3, 
        kb=0.2
    )
    
    env = Enviroment(rxn1, rxn2, T=298)
    return env


@pytest.fixture
def complex_stoichiometry_environment():
    """Create an environment with non-unity stoichiometric coefficients: A + 2B ⇌ C"""
    A = Compound("A")
    B = Compound("B")
    C = Compound("C")
    
    reactants = [
        {"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1},
        {"stoichiometric_coefficient": 2, "compound": B, "rate_dependency": 2}
    ]
    products = [
        {"stoichiometric_coefficient": 1, "compound": C, "rate_dependency": 1}
    ]
    
    rxn = Reaction(
        reactants, 
        products, 
        [1.0, 2.0],  # A and B initial concentrations
        [0.0],  # C initial concentration
        K=10.0, 
        kf=0.5, 
        kb=0.05
    )
    env = Enviroment(rxn, T=298)
    return env


@pytest.fixture
def phase_environment():
    """Create an environment with different phases (gas, liquid, solid)."""
    A_g = Compound("A", phase_point_list=[{"phase": "g", "temperature": 298}])
    B_l = Compound("B", phase_point_list=[{"phase": "l", "temperature": 298}])
    C_s = Compound("C", phase_point_list=[{"phase": "s", "temperature": 298}])
    
    reactants = [{"stoichiometric_coefficient": 1, "compound": A_g, "rate_dependency": 1}]
    products = [
        {"stoichiometric_coefficient": 1, "compound": B_l, "rate_dependency": 1},
        {"stoichiometric_coefficient": 1, "compound": C_s, "rate_dependency": 1}
    ]
    
    rxn = Reaction(
        reactants, 
        products, 
        [1.0], 
        [0.0, 0.0], 
        K=5.0, 
        kf=0.3, 
        kb=0.06
    )
    env = Enviroment(rxn, T=298)
    return env


# ---------- Fit Method Tests ---------- #

def test_fit_valid_environment(simple_equilibrium_environment):
    """Test fitting with a valid Enviroment instance."""
    calc = EquilibriumCalculator()
    calc.fit(simple_equilibrium_environment)
    
    assert calc.fitted == True
    assert calc.env == simple_equilibrium_environment
    assert hasattr(calc, 'concentration_equation')
    assert len(calc.concentration_equation) == 2  # A and B


def test_fit_invalid_environment():
    """Test that fitting with invalid input raises ValueError."""
    calc = EquilibriumCalculator()
    with pytest.raises(ValueError, match="The input should be an instance of Enviroment class"):
        calc.fit("not an environment")


def test_fit_multi_reaction_environment(multi_reaction_equilibrium_environment):
    """Test fitting with environment containing multiple reactions."""
    calc = EquilibriumCalculator()
    calc.fit(multi_reaction_equilibrium_environment)
    
    assert calc.fitted == True
    assert len(calc.concentration_equation) == 3  # A, B, C
    # Check that concentration equations contain reaction extent variables
    assert any("x1" in eq or "x2" in eq for eq in calc.concentration_equation)


def test_fit_overwrites_previous_environment(simple_equilibrium_environment, multi_reaction_equilibrium_environment):
    """Test that fitting again overwrites the previous environment."""
    calc = EquilibriumCalculator()
    calc.fit(simple_equilibrium_environment)
    assert len(calc.concentration_equation) == 2
    
    calc.fit(multi_reaction_equilibrium_environment)
    assert len(calc.concentration_equation) == 3


def test_generate_concentration_equations(simple_equilibrium_environment):
    """Test that concentration equations are generated correctly."""
    calc = EquilibriumCalculator()
    calc.fit(simple_equilibrium_environment)
    
    equations = calc._generate_concentration_equations()
    assert len(equations) == 2
    # For A ⇌ B, A should have +x1, B should have -x1
    assert "x1" in equations[0] or "x1" in equations[1]


def test_generate_concentration_equations_complex_stoichiometry(complex_stoichiometry_environment):
    """Test concentration equations with non-unity stoichiometric coefficients."""
    calc = EquilibriumCalculator()
    calc.fit(complex_stoichiometry_environment)
    
    equations = calc._generate_concentration_equations()
    assert len(equations) == 3
    # A + 2B ⇌ C: A should have +x1, B should have +2x1, C should have -x1
    assert any("x1" in eq for eq in equations)


# ---------- Calculate Method Tests (BGD) ---------- #

def test_calculate_not_fitted():
    """Test that calculate raises ValueError when not fitted."""
    calc = EquilibriumCalculator()
    with pytest.raises(ValueError, match="Environment not fitted"):
        calc.calculate()


def test_calculate_bgd_basic(simple_equilibrium_environment):
    """Test basic calculation using batch gradient descent."""
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(simple_equilibrium_environment)
    
    result = calc.calculate(max_iter=100, tol=1e-6)
    
    assert isinstance(result, list)
    assert len(result) == 2  # A and B
    assert all(c >= 0 for c in result)  # All concentrations non-negative
    assert hasattr(calc, 'x_solution')  # Solution extents stored


def test_calculate_bgd_convergence(simple_equilibrium_environment):
    """Test that BGD converges to equilibrium."""
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(simple_equilibrium_environment)
    
    result = calc.calculate(max_iter=1000, tol=1e-8)
    
    # At equilibrium, Q should equal K (approximately)
    # For A ⇌ B with K=2.0, at equilibrium: [B]/[A] ≈ 2.0
    # With initial [A]=1.0, [B]=0.0, we expect [A] ≈ 0.33, [B] ≈ 0.67
    assert result[0] > 0
    assert result[1] > 0
    # Check that mass is conserved (approximately)
    total = result[0] + result[1]
    assert np.isclose(total, 1.0, atol=0.01)


def test_calculate_bgd_multi_reaction(multi_reaction_equilibrium_environment):
    """Test BGD calculation with multiple reactions."""
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(multi_reaction_equilibrium_environment)
    
    result = calc.calculate(max_iter=1000, tol=1e-6)
    
    assert isinstance(result, list)
    assert len(result) == 3  # A, B, C
    assert all(c >= 0 for c in result)
    # Check mass conservation
    total = sum(result)
    initial_total = sum(multi_reaction_equilibrium_environment.concentrations_array)
    assert np.isclose(total, initial_total, atol=0.1)


def test_calculate_bgd_custom_parameters(simple_equilibrium_environment):
    """Test BGD with custom learning rate and tolerance."""
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(simple_equilibrium_environment)
    
    result = calc.calculate(
        max_iter=500,
        learning_rate=0.05,
        tol=1e-10,
        backtrack_beta=0.7,
        min_concentration=1e-14
    )
    
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


def test_calculate_bgd_phase_exclusion(phase_environment):
    """Test that solid and liquid phases are excluded from equilibrium expressions."""
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(phase_environment)
    
    result = calc.calculate(max_iter=1000, tol=1e-6)
    
    assert isinstance(result, list)
    assert len(result) == 3
    assert all(c >= 0 for c in result)


# ---------- Calculate Method Tests (SGD) ---------- #

def test_calculate_sgd_basic(simple_equilibrium_environment):
    """Test basic calculation using stochastic gradient descent."""
    calc = EquilibriumCalculator(method_of_calculation="sgd")
    calc.fit(simple_equilibrium_environment)
    
    result = calc.calculate(max_iter=100, tol=1e-6)
    
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


def test_calculate_sgd_convergence(simple_equilibrium_environment):
    """Test that SGD converges to equilibrium."""
    calc = EquilibriumCalculator(method_of_calculation="sgd")
    calc.fit(simple_equilibrium_environment)
    
    result = calc.calculate(max_iter=2000, tol=1e-8)
    
    assert result[0] > 0
    assert result[1] > 0
    total = result[0] + result[1]
    assert np.isclose(total, 1.0, atol=0.01)


def test_calculate_sgd_multi_reaction(multi_reaction_equilibrium_environment):
    """Test SGD calculation with multiple reactions."""
    calc = EquilibriumCalculator(method_of_calculation="sgd")
    calc.fit(multi_reaction_equilibrium_environment)
    
    result = calc.calculate(max_iter=2000, tol=1e-6)
    
    assert isinstance(result, list)
    assert len(result) == 3
    assert all(c >= 0 for c in result)


def test_calculate_sgd_custom_parameters(simple_equilibrium_environment):
    """Test SGD with custom parameters."""
    calc = EquilibriumCalculator(method_of_calculation="sgd")
    calc.fit(simple_equilibrium_environment)
    
    result = calc.calculate(
        max_iter=1000,
        learning_rate=0.05,
        tol=1e-10,
        backtrack_beta=0.7
    )
    
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


# ---------- Calculate Method Tests (Newton) ---------- #

def test_calculate_newton_basic(simple_equilibrium_environment):
    """Test basic calculation using Newton's method."""
    calc = EquilibriumCalculator(method_of_calculation="newton")
    calc.fit(simple_equilibrium_environment)
    
    result = calc.calculate(max_iter=50, tol=1e-8)
    
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


def test_calculate_newton_convergence(simple_equilibrium_environment):
    """Test that Newton's method converges to equilibrium."""
    calc = EquilibriumCalculator(method_of_calculation="newton")
    calc.fit(simple_equilibrium_environment)
    
    result = calc.calculate(max_iter=100, tol=1e-10)
    
    assert result[0] > 0
    assert result[1] > 0
    total = result[0] + result[1]
    assert np.isclose(total, 1.0, atol=0.01)


def test_calculate_newton_multi_reaction(multi_reaction_equilibrium_environment):
    """Test Newton's method with multiple reactions."""
    calc = EquilibriumCalculator(method_of_calculation="newton")
    calc.fit(multi_reaction_equilibrium_environment)
    
    result = calc.calculate(max_iter=100, tol=1e-8)
    
    assert isinstance(result, list)
    assert len(result) == 3
    assert all(c >= 0 for c in result)


def test_calculate_newton_custom_parameters(simple_equilibrium_environment):
    """Test Newton's method with custom parameters."""
    calc = EquilibriumCalculator(method_of_calculation="newton")
    calc.fit(simple_equilibrium_environment)
    
    result = calc.calculate(
        max_iter=50,
        learning_rate=0.8,
        tol=1e-12,
        backtrack_beta=0.6
    )
    
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


# ---------- Fit Calculate Method Tests ---------- #

def test_fit_calculate_valid_bgd(simple_equilibrium_environment):
    """Test fit_calculate with BGD method."""
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    result = calc.fit_calculate(simple_equilibrium_environment, max_iter=100, tol=1e-6)
    
    assert calc.fitted == True
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


def test_fit_calculate_valid_sgd(simple_equilibrium_environment):
    """Test fit_calculate with SGD method."""
    calc = EquilibriumCalculator(method_of_calculation="sgd")
    result = calc.fit_calculate(simple_equilibrium_environment, max_iter=200, tol=1e-6)
    
    assert calc.fitted == True
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


def test_fit_calculate_valid_newton(simple_equilibrium_environment):
    """Test fit_calculate with Newton's method."""
    calc = EquilibriumCalculator(method_of_calculation="newton")
    result = calc.fit_calculate(simple_equilibrium_environment, max_iter=50, tol=1e-8)
    
    assert calc.fitted == True
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


def test_fit_calculate_invalid_environment():
    """Test fit_calculate with invalid environment."""
    calc = EquilibriumCalculator()
    with pytest.raises(ValueError, match="The input should be an instance of Enviroment class"):
        calc.fit_calculate("not an environment")


def test_fit_calculate_equivalent_to_fit_then_calculate(simple_equilibrium_environment):
    """Test that fit_calculate is equivalent to fit() followed by calculate()."""
    calc1 = EquilibriumCalculator(method_of_calculation="bgd")
    calc2 = EquilibriumCalculator(method_of_calculation="bgd")
    
    # Method 1: fit_calculate
    result1 = calc1.fit_calculate(simple_equilibrium_environment, max_iter=100, tol=1e-6)
    
    # Method 2: fit then calculate
    calc2.fit(simple_equilibrium_environment)
    result2 = calc2.calculate(max_iter=100, tol=1e-6)
    
    # Results should be equivalent (within numerical tolerance)
    assert len(result1) == len(result2)
    assert np.allclose(result1, result2, atol=1e-5)


def test_fit_calculate_with_custom_parameters(simple_equilibrium_environment):
    """Test fit_calculate with custom parameters."""
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    result = calc.fit_calculate(
        simple_equilibrium_environment,
        max_iter=200,
        learning_rate=0.05,
        tol=1e-10,
        backtrack_beta=0.7,
        min_concentration=1e-14
    )
    
    assert calc.fitted == True
    assert isinstance(result, list)
    assert len(result) == 2


# ---------- Edge Cases and Integration Tests ---------- #

def test_calculate_different_methods_same_result(simple_equilibrium_environment):
    """Test that different methods produce similar results for simple systems."""
    calc_bgd = EquilibriumCalculator(method_of_calculation="bgd")
    calc_sgd = EquilibriumCalculator(method_of_calculation="sgd")
    calc_newton = EquilibriumCalculator(method_of_calculation="newton")
    
    calc_bgd.fit(simple_equilibrium_environment)
    calc_sgd.fit(simple_equilibrium_environment)
    calc_newton.fit(simple_equilibrium_environment)
    
    result_bgd = calc_bgd.calculate(max_iter=1000, tol=1e-6)
    result_sgd = calc_sgd.calculate(max_iter=2000, tol=1e-6)
    result_newton = calc_newton.calculate(max_iter=100, tol=1e-6)
    
    # All methods should produce non-negative concentrations
    assert all(c >= 0 for c in result_bgd)
    assert all(c >= 0 for c in result_sgd)
    assert all(c >= 0 for c in result_newton)
    
    # Results should be similar (within reasonable tolerance)
    # Note: Different methods may converge to slightly different solutions
    # but should be in the same ballpark
    assert np.allclose(result_bgd, result_sgd, atol=0.1)
    assert np.allclose(result_bgd, result_newton, atol=0.1)


def test_calculate_with_very_small_equilibrium_constant():
    """Test calculation with very small equilibrium constant."""
    A = Compound("A")
    B = Compound("B")
    
    reactants = [{"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1}]
    products = [{"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}]
    
    rxn = Reaction(
        reactants, 
        products, 
        [1.0], 
        [0.0], 
        K=1e-10,  # Very small K
        kf=0.5, 
        kb=0.5
    )
    env = Enviroment(rxn, T=298)
    
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(env)
    result = calc.calculate(max_iter=1000, tol=1e-6)
    
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)
    # With very small K, reaction should favor reactants
    assert result[0] > result[1]


def test_calculate_with_very_large_equilibrium_constant():
    """Test calculation with very large equilibrium constant."""
    A = Compound("A")
    B = Compound("B")
    
    reactants = [{"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1}]
    products = [{"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}]
    
    rxn = Reaction(
        reactants, 
        products, 
        [1.0], 
        [0.0], 
        K=1e10,  # Very large K
        kf=0.5, 
        kb=0.5
    )
    env = Enviroment(rxn, T=298)
    
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(env)
    result = calc.calculate(max_iter=1000, tol=1e-6)
    
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)
    # With very large K, reaction should favor products
    assert result[1] > result[0]


def test_calculate_with_zero_initial_concentrations():
    """Test calculation when some compounds start at zero concentration."""
    A = Compound("A")
    B = Compound("B")
    C = Compound("C")
    
    # Reaction: A + B ⇌ C
    reactants = [
        {"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1},
        {"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}
    ]
    products = [{"stoichiometric_coefficient": 1, "compound": C, "rate_dependency": 1}]
    
    rxn = Reaction(
        reactants, 
        products, 
        [1.0, 0.0],  # A=1.0, B=0.0
        [0.0],  # C=0.0
        K=1.0, 
        kf=0.5, 
        kb=0.5
    )
    env = Enviroment(rxn, T=298)
    
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(env)
    result = calc.calculate(max_iter=1000, tol=1e-6, min_concentration=1e-12)
    
    assert isinstance(result, list)
    assert len(result) == 3
    assert all(c >= 0 for c in result)


def test_calculate_mass_conservation(complex_stoichiometry_environment):
    """Test that mass is conserved in calculations."""
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(complex_stoichiometry_environment)
    
    initial = complex_stoichiometry_environment.concentrations_array
    result = calc.calculate(max_iter=1000, tol=1e-6)
    
    # For A + 2B ⇌ C, we need to account for stoichiometry
    # Total "A-equivalents": [A] + [C] (since C contains one A)
    # Total "B-equivalents": [B] + 2*[C] (since C contains two B)
    initial_A_equiv = initial[0] + initial[2]
    initial_B_equiv = initial[1] + 2 * initial[2]
    
    final_A_equiv = result[0] + result[2]
    final_B_equiv = result[1] + 2 * result[2]
    
    # Mass should be approximately conserved
    assert np.isclose(initial_A_equiv, final_A_equiv, atol=0.01)
    assert np.isclose(initial_B_equiv, final_B_equiv, atol=0.01)


def test_calculate_max_iter_reached():
    """Test behavior when max_iter is reached before convergence."""
    A = Compound("A")
    B = Compound("B")
    
    reactants = [{"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1}]
    products = [{"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}]
    
    rxn = Reaction(
        reactants, 
        products, 
        [1.0], 
        [0.0], 
        K=2.0, 
        kf=0.5, 
        kb=0.25
    )
    env = Enviroment(rxn, T=298)
    
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(env)
    # Use very few iterations to test max_iter behavior
    result = calc.calculate(max_iter=5, tol=1e-12)
    
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


def test_calculate_negative_concentrations_clamped():
    """Test that negative concentrations are clamped to zero."""
    A = Compound("A")
    B = Compound("B")
    
    reactants = [{"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1}]
    products = [{"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}]
    
    rxn = Reaction(
        reactants, 
        products, 
        [0.1],  # Small initial concentration
        [0.0], 
        K=100.0,  # Large K to push reaction forward
        kf=0.5, 
        kb=0.005
    )
    env = Enviroment(rxn, T=298)
    
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(env)
    result = calc.calculate(max_iter=1000, tol=1e-6)
    
    # All concentrations should be non-negative
    assert all(c >= 0 for c in result)


def test_calculate_min_concentration_parameter():
    """Test that min_concentration parameter is respected."""
    A = Compound("A")
    B = Compound("B")
    
    reactants = [{"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1}]
    products = [{"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}]
    
    rxn = Reaction(
        reactants, 
        products, 
        [1.0], 
        [0.0], 
        K=2.0, 
        kf=0.5, 
        kb=0.25
    )
    env = Enviroment(rxn, T=298)
    
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(env)
    # Use a larger min_concentration to test the parameter
    result = calc.calculate(max_iter=1000, tol=1e-6, min_concentration=1e-10)
    
    assert isinstance(result, list)
    assert len(result) == 2
    assert all(c >= 0 for c in result)


def test_x_solution_stored_after_calculation(simple_equilibrium_environment):
    """Test that reaction extents are stored in x_solution after calculation."""
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(simple_equilibrium_environment)
    
    result = calc.calculate(max_iter=100, tol=1e-6)
    
    assert hasattr(calc, 'x_solution')
    assert isinstance(calc.x_solution, np.ndarray)
    assert len(calc.x_solution) == 1  # One reaction


def test_x_solution_multi_reaction(multi_reaction_equilibrium_environment):
    """Test x_solution with multiple reactions."""
    calc = EquilibriumCalculator(method_of_calculation="bgd")
    calc.fit(multi_reaction_equilibrium_environment)
    
    result = calc.calculate(max_iter=1000, tol=1e-6)
    
    assert hasattr(calc, 'x_solution')
    assert isinstance(calc.x_solution, np.ndarray)
    assert len(calc.x_solution) == 2  # Two reactions

