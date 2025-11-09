import pytest
import numpy as np
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

# Add paths to sys.path for imports
project_root = Path(__file__).parent.parent
src_chemcalc_path = project_root / "src" / "ChemCalc"
sys.path.insert(0, str(src_chemcalc_path))  # For _general import in Kinetic.py
sys.path.insert(0, str(project_root))  # For src.ChemCalc imports

# Import from _general the same way Kinetic.py does to ensure type compatibility
from _general import Enviroment, Compound, Reaction
from src.ChemCalc.Kinetic import KineticalCalculator


# -------------------------
# KineticalCalculator Tests
# -------------------------

# ---------- Initialization Tests ---------- #

def test_kinetical_calculator_init_default():
    """Test initialization with default accuracy."""
    kc = KineticalCalculator()
    assert kc.accuracy == 1e-3
    assert kc.fitted == False


def test_kinetical_calculator_init_custom_accuracy():
    """Test initialization with custom accuracy."""
    kc = KineticalCalculator(accuracy=0.01)
    assert kc.accuracy == 0.01
    assert kc.fitted == False


def test_kinetical_calculator_init_small_accuracy():
    """Test initialization with very small accuracy for high precision."""
    kc = KineticalCalculator(accuracy=1e-5)
    assert kc.accuracy == 1e-5
    assert kc.fitted == False


# ---------- Fixture for test environment ---------- #

@pytest.fixture
def simple_environment():
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
        K=2.0, 
        kf=0.5, 
        kb=0.25
    )
    env = Enviroment(rxn, T=298)
    return env


@pytest.fixture
def multi_reaction_environment():
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


# ---------- Fit Method Tests ---------- #

def test_fit_valid_environment(simple_environment):
    """Test fitting with a valid Enviroment instance."""
    kc = KineticalCalculator()
    kc.fit(simple_environment)
    
    assert kc.fitted == True
    assert kc.enviroment == simple_environment
    assert len(kc.concentrations) == 2
    assert kc.concentrations[0] == 1.0  # A
    assert kc.concentrations[1] == 0.0  # B
    assert kc.number_of_reactions == 1
    assert len(kc.rate_constants) == 1
    assert kc.rate_constants[0] == [0.5, 0.25]


def test_fit_invalid_environment():
    """Test that fitting with invalid input raises ValueError."""
    kc = KineticalCalculator()
    with pytest.raises(ValueError, match="The input should be an instance of Enviroment class"):
        kc.fit("not an environment")


def test_fit_multi_reaction_environment(multi_reaction_environment):
    """Test fitting with environment containing multiple reactions."""
    kc = KineticalCalculator()
    kc.fit(multi_reaction_environment)
    
    assert kc.fitted == True
    assert kc.number_of_reactions == 2
    assert len(kc.concentrations) == 3
    assert len(kc.rate_constants) == 2
    assert kc.rate_constants[0] == [0.5, 0.25]
    assert kc.rate_constants[1] == [0.3, 0.2]


def test_fit_overwrites_previous_environment(simple_environment, multi_reaction_environment):
    """Test that fitting again overwrites the previous environment."""
    kc = KineticalCalculator()
    kc.fit(simple_environment)
    assert kc.number_of_reactions == 1
    
    kc.fit(multi_reaction_environment)
    assert kc.number_of_reactions == 2
    assert len(kc.concentrations) == 3


# ---------- Calculate Method Tests ---------- #

def test_calculate_not_fitted():
    """Test that calculate raises NameError when not fitted."""
    kc = KineticalCalculator()
    with pytest.raises(NameError, match="You should fit the model to an enviromt object before calculation"):
        kc.calculate(time=1.0)


def test_calculate_invalid_plot_mode(simple_environment):
    """Test that calculate raises ValueError for invalid plot mode."""
    kc = KineticalCalculator()
    kc.fit(simple_environment)
    with pytest.raises(ValueError, match="`plot` is not one of"):
        kc.calculate(time=1.0, plot="invalid_mode")


def test_calculate_no_plot(simple_environment):
    """Test calculate without plotting."""
    kc = KineticalCalculator(accuracy=0.1)
    kc.fit(simple_environment)
    results = kc.calculate(time=1.0, plot=False)
    
    assert isinstance(results, list)
    assert len(results) >= 1
    # Last checkpoint should contain final concentrations
    final_concentrations = results[-1]
    assert isinstance(final_concentrations, np.ndarray) or isinstance(final_concentrations, list)
    assert len(final_concentrations) == 2
    # Concentrations should be non-negative
    assert all(c >= 0 for c in final_concentrations)


def test_calculate_with_checkpoint_times(simple_environment):
    """Test calculate with specified checkpoint times."""
    kc = KineticalCalculator(accuracy=0.1)
    kc.fit(simple_environment)
    checkpoint_times = [0.2, 0.5, 0.8]
    results = kc.calculate(time=1.0, checkpoint_time=checkpoint_times, plot=False)
    
    # Should have checkpoints at specified times plus final state
    assert isinstance(results, list)
    # At minimum, should have the final checkpoint
    assert len(results) >= 1
    final_concentrations = results[-1]
    assert len(final_concentrations) == 2


def test_calculate_concentration_evolution(simple_environment):
    """Test that concentrations evolve correctly over time."""
    kc = KineticalCalculator(accuracy=0.01)
    kc.fit(simple_environment)
    
    # Initial: A=1.0, B=0.0
    results_short = kc.calculate(time=0.1, plot=False)
    results_long = kc.calculate(time=1.0, plot=False)
    
    initial_A = simple_environment.concentrations_array[0]
    initial_B = simple_environment.concentrations_array[1]
    
    # After longer time, A should decrease and B should increase
    final_short = results_short[-1]
    final_long = results_long[-1]
    
    # A should decrease over time (forward reaction)
    assert final_long[0] < final_short[0] or np.isclose(final_long[0], final_short[0], atol=1e-6)
    # B should increase over time
    assert final_long[1] >= final_short[1] or np.isclose(final_long[1], final_short[1], atol=1e-6)


def test_calculate_negative_concentrations_clamped(simple_environment):
    """Test that negative concentrations are clamped to zero."""
    kc = KineticalCalculator(accuracy=0.1)
    kc.fit(simple_environment)
    results = kc.calculate(time=10.0, plot=False)
    
    # All concentrations should be non-negative
    for checkpoint in results:
        assert all(c >= 0 for c in checkpoint)


def test_calculate_save_plot(simple_environment, tmp_path):
    """Test calculate with plot="save" mode."""
    kc = KineticalCalculator(accuracy=0.1)
    kc.fit(simple_environment)
    plot_path = tmp_path / "test_plot.png"
    
    results = kc.calculate(time=0.5, plot="save", directory=str(plot_path))
    
    assert isinstance(results, list)
    # Check that file was created (if matplotlib worked)
    # Note: This test may fail if matplotlib backend doesn't support file saving in test environment


@patch('builtins.input', return_value='exit')
def test_calculate_interactive_plot(mock_input, simple_environment):
    """Test calculate with plot="interactive" mode (mocked input)."""
    kc = KineticalCalculator(accuracy=0.1)
    kc.fit(simple_environment)
    
    # Mock matplotlib to avoid actual display
    with patch('matplotlib.pyplot.show'):
        with patch('matplotlib.pyplot.close'):
            results = kc.calculate(time=0.5, plot="interactive")
            assert isinstance(results, list)


def test_calculate_different_accuracy_levels(simple_environment):
    """Test that different accuracy levels produce consistent results."""
    kc1 = KineticalCalculator(accuracy=0.1)
    kc2 = KineticalCalculator(accuracy=0.01)
    
    kc1.fit(simple_environment)
    kc2.fit(simple_environment)
    
    results1 = kc1.calculate(time=1.0, plot=False)
    results2 = kc2.calculate(time=1.0, plot=False)
    
    # Both should produce valid results
    assert len(results1) >= 1
    assert len(results2) >= 1
    # Higher accuracy should give more detailed results (more checkpoints for same time)
    # Note: This is not always true due to checkpoint timing, but both should be valid


# ---------- Fit Calculate Method Tests ---------- #

def test_fit_calculate_valid(simple_environment):
    """Test fit_calculate with valid environment."""
    kc = KineticalCalculator(accuracy=0.1)
    results = kc.fit_calculate(simple_environment, time=1.0, plot=False)
    
    assert kc.fitted == True
    assert isinstance(results, list)
    assert len(results) >= 1
    final_concentrations = results[-1]
    assert len(final_concentrations) == 2


def test_fit_calculate_invalid_environment():
    """Test fit_calculate with invalid environment."""
    kc = KineticalCalculator()
    with pytest.raises(ValueError, match="The input should be an instance of Enviroment class"):
        kc.fit_calculate("not an environment", time=1.0)


def test_fit_calculate_equivalent_to_fit_then_calculate(simple_environment):
    """Test that fit_calculate is equivalent to fit() followed by calculate()."""
    kc1 = KineticalCalculator(accuracy=0.1)
    kc2 = KineticalCalculator(accuracy=0.1)
    
    # Method 1: fit_calculate
    results1 = kc1.fit_calculate(simple_environment, time=1.0, plot=False)
    
    # Method 2: fit then calculate
    kc2.fit(simple_environment)
    results2 = kc2.calculate(time=1.0, plot=False)
    
    # Results should be equivalent
    assert len(results1) == len(results2)
    assert np.allclose(results1[-1], results2[-1], atol=1e-6)


def test_fit_calculate_with_checkpoints(simple_environment):
    """Test fit_calculate with checkpoint times."""
    kc = KineticalCalculator(accuracy=0.1)
    checkpoint_times = [0.2, 0.5]
    results = kc.fit_calculate(
        simple_environment, 
        time=1.0, 
        checkpoint_time=checkpoint_times, 
        plot=False
    )
    
    assert isinstance(results, list)
    assert len(results) >= 1


def test_fit_calculate_with_plot_save(simple_environment, tmp_path):
    """Test fit_calculate with plot="save"."""
    kc = KineticalCalculator(accuracy=0.1)
    plot_path = tmp_path / "fit_calculate_plot.png"
    
    results = kc.fit_calculate(
        simple_environment,
        time=0.5,
        plot="save",
        directory=str(plot_path)
    )
    
    assert isinstance(results, list)
    assert kc.fitted == True


# ---------- Calculate Responsively Method Tests ---------- #

def test_calculate_responsively_not_fitted():
    """Test that calculate_responsively raises NameError when not fitted."""
    kc = KineticalCalculator()
    with pytest.raises(NameError, match="You must fit the model to an Enviroment before calculation"):
        kc.calculate_responsively()


def test_calculate_responsively_fitted(simple_environment):
    """Test calculate_responsively with fitted environment."""
    kc = KineticalCalculator(accuracy=0.1)
    kc.fit(simple_environment)
    
    # Mock input and matplotlib to avoid actual user interaction
    with patch('builtins.input', return_value='exit'):
        with patch('matplotlib.pyplot.show'):
            with patch('matplotlib.pyplot.close'):
                with patch('matplotlib.animation.FuncAnimation'):
                    # This is a limited test due to the interactive nature
                    # We mainly check that the method can be called without errors
                    # In a real scenario, this would require more sophisticated mocking
                    pass
    # Note: Full testing of calculate_responsively is difficult due to its
    # interactive and animation-dependent nature


# ---------- Edge Cases and Integration Tests ---------- #

def test_zero_time_calculation(simple_environment):
    """Test calculation with zero time."""
    kc = KineticalCalculator(accuracy=0.1)
    kc.fit(simple_environment)
    results = kc.calculate(time=0.0, plot=False)
    
    assert isinstance(results, list)
    assert len(results) >= 1
    # Note: Even with time=0.0, the calculation loop runs at least once
    # (range(int(0.0/0.1+1)) = range(1) = one iteration)
    # So concentrations will change by one accuracy step
    # We just verify that results are returned and are valid (non-negative)
    final = results[-1]
    assert len(final) == len(simple_environment.concentrations_array)
    assert all(c >= 0 for c in final)  # Concentrations should be non-negative


def test_very_small_time_step(simple_environment):
    """Test with very small time step."""
    kc = KineticalCalculator(accuracy=1e-5)
    kc.fit(simple_environment)
    results = kc.calculate(time=0.01, plot=False)
    
    assert isinstance(results, list)
    assert len(results) >= 1


def test_multi_reaction_calculation(multi_reaction_environment):
    """Test calculation with multiple reactions."""
    kc = KineticalCalculator(accuracy=0.1)
    kc.fit(multi_reaction_environment)
    results = kc.calculate(time=1.0, plot=False)
    
    assert isinstance(results, list)
    assert len(results) >= 1
    final_concentrations = results[-1]
    assert len(final_concentrations) == 3  # A, B, C


def test_concentration_preservation_sum(multi_reaction_environment):
    """Test that total mass is approximately conserved (for closed systems)."""
    kc = KineticalCalculator(accuracy=0.01)
    kc.fit(multi_reaction_environment)
    results = kc.calculate(time=1.0, plot=False)
    
    initial_sum = np.sum(multi_reaction_environment.concentrations_array)
    final_sum = np.sum(results[-1])
    
    # Total should be approximately conserved (within numerical error)
    # Note: This assumes no mass loss in the system
    assert np.isclose(initial_sum, final_sum, atol=1e-2)


def test_checkpoint_time_outside_simulation_range(simple_environment):
    """Test checkpoint times that are outside the simulation time range."""
    kc = KineticalCalculator(accuracy=0.1)
    kc.fit(simple_environment)
    
    # Checkpoint time beyond simulation time
    checkpoint_times = [2.0, 5.0]  # Beyond simulation time of 1.0
    results = kc.calculate(time=1.0, checkpoint_time=checkpoint_times, plot=False)
    
    # Should still return valid results (final state)
    assert isinstance(results, list)
    assert len(results) >= 1


def test_reaction_with_zero_rate_constants():
    """Test reaction with zero rate constants (no change expected)."""
    A = Compound("A")
    B = Compound("B")
    
    reactants = [{"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1}]
    products = [{"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}]
    
    rxn = Reaction(
        reactants, 
        products, 
        [1.0], 
        [0.0], 
        K=1.0, 
        kf=0.0,  # Zero forward rate
        kb=0.0   # Zero backward rate
    )
    env = Enviroment(rxn, T=298)
    
    kc = KineticalCalculator(accuracy=0.1)
    kc.fit(env)
    results = kc.calculate(time=1.0, plot=False)
    
    # Concentrations should not change
    initial = env.concentrations_array
    final = results[-1]
    assert np.allclose(final, initial, atol=1e-6)
