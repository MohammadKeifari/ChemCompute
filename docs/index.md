# ChemCalc Documentation

Welcome to the ChemCalc documentation! ChemCalc is a comprehensive Python library for chemical reaction simulation, including kinetic modeling and thermodynamic equilibrium calculations.

## Overview

ChemCalc provides a powerful and flexible framework for:

- **Chemical Compound Management**: Create and manage chemical compounds with formulas, phases, and physical properties
- **Reaction Definition**: Define complex chemical reactions with multiple reactants and products
- **Temperature-Dependent Calculations**: Automatic updates of rate constants and equilibrium constants using Arrhenius and van't Hoff equations
- **Thermodynamic Properties**: Support for enthalpy, entropy, and activation energies for realistic temperature-dependent simulations
- **Kinetic Simulation**: Simulate time-dependent concentration changes using numerical integration
- **Equilibrium Calculations**: Calculate equilibrium concentrations using advanced optimization algorithms
- **Phase Handling**: Support for solid, liquid, gas, and aqueous phases with temperature-dependent transitions
- **Visualization**: Interactive and static plotting capabilities for analyzing reaction kinetics

## Quick Navigation

### Getting Started

- [Installation Guide](installation.md) - How to install and set up ChemCalc
- [Quick Start Tutorial](quickstart.md) - Get up and running in minutes
- [Basic Concepts](concepts.md) - Understanding the core concepts

### Core Components

- [Compound Class](compound.md) - Working with chemical compounds
- [Reaction Class](reaction.md) - Defining and managing chemical reactions
- [Enviroment Class](environment.md) - Managing reaction systems
- [KineticalCalculator](kinetic.md) - Simulating reaction kinetics
- [EquilibriumCalculator](equilibrium.md) - Calculating equilibrium states

### Advanced Topics

- [Reaction Syntax Guide](syntax.md) - Understanding reaction string formats
- [Phase Handling](phases.md) - Working with different phases
- [Optimization Methods](optimization.md) - Understanding equilibrium calculation algorithms
- [Visualization](visualization.md) - Creating plots and visualizations
- [Performance Tips](performance.md) - Optimizing your simulations

### Examples

- [Basic Examples](examples/basic.md) - Simple reaction examples
- [Advanced Examples](examples/advanced.md) - Complex multi-reaction systems
- [Real-World Applications](examples/applications.md) - Practical use cases

### API Reference

- [Full API Documentation](api/index.md) - Complete API reference
- [Compound API](api/compound.md) - Compound class methods and properties
- [Reaction API](api/reaction.md) - Reaction class methods and properties
- [Enviroment API](api/environment.md) - Enviroment class methods and properties
- [KineticalCalculator API](api/kinetic.md) - KineticalCalculator methods
- [EquilibriumCalculator API](api/equilibrium.md) - EquilibriumCalculator methods

### Additional Resources

- [Troubleshooting](troubleshooting.md) - Common issues and solutions
- [FAQ](faq.md) - Frequently asked questions
- [Contributing](contributing.md) - How to contribute to ChemCalc
- [Changelog](changelog.md) - Version history and changes

## Key Features

### 🧪 Chemical Compound Representation

Create compounds with formulas, phases, and physical properties:

```python
from ChemCalc import Compound

# Simple compound
water = Compound("H2O")

# With phase information
co2 = Compound("CO2", phase_point_list=[{"phase": "g", "temperature": 298}])

# With melting/boiling points
ethanol = Compound("C2H5OH", mp=-114, bp=78)
```

### ⚗️ Reaction Definition

Define reactions using simple or complex syntax with thermodynamic parameters:

```python
from ChemCalc import Reaction

# Simple syntax with thermodynamic parameters
rxn = Reaction.from_string_simple_syntax(
    "2A + B > 3C",
    concentrations=[1.0, 1.0, 0.0],
    K=10.0,
    kf=0.5,
    kb=0.05,
    enthalpy=-50000,  # J/mol
    entropy=-100,     # J/(mol·K)
    activation_energy_forward=50000,   # J/mol
    activation_energy_backward=100000  # J/mol
)
```

### 🌡️ Temperature-Dependent Calculations

Automatically update rate constants and equilibrium constants with temperature:

```python
# Create reaction with thermodynamic parameters
rxn = Reaction.from_string_simple_syntax(
    "A > B",
    K=2.0,
    kf=0.5,
    kb=0.25,
    enthalpy=-50000,
    activation_energy_forward=50000,
    activation_energy_backward=100000,
    T=298
)

# Change temperature - K, kf, kb automatically update
rxn.T = 350  # Uses Arrhenius and van't Hoff equations
```

### ⏱️ Kinetic Simulation

Simulate reaction kinetics over time:

```python
from ChemCalc.Kinetic import KineticalCalculator

kc = KineticalCalculator(accuracy=1e-3)
kc.fit(env)
results = kc.calculate(time=10.0, plot="interactive")
```

### ⚖️ Equilibrium Calculations

Calculate equilibrium concentrations using multiple algorithms:

```python
from ChemCalc.Thermodynamic import EquilibriumCalculator

eq_calc = EquilibriumCalculator(method_of_calculation="newton")
equilibrium = eq_calc.fit_calculate(env, max_iter=1000, tol=1e-8)
```

## Installation

```bash
pip install numpy matplotlib
```

For detailed installation instructions, see the [Installation Guide](installation.md).

## Quick Example

Here's a complete example to get you started:

```python
from ChemCalc import Compound, Reaction, Enviroment
from ChemCalc.Kinetic import KineticalCalculator
from ChemCalc.Thermodynamic import EquilibriumCalculator

# Create a simple reversible reaction: A ⇌ B
rxn = Reaction.from_string_simple_syntax(
    "A > B",
    concentrations=[1.0, 0.0],
    K=2.0,
    kf=0.5,
    kb=0.25
)

# Create environment
env = Enviroment(rxn, T=298)

# Kinetic simulation
kc = KineticalCalculator(accuracy=0.01)
kc.fit(env)
kinetic_results = kc.calculate(time=10.0, plot=False)

# Equilibrium calculation
eq_calc = EquilibriumCalculator(method_of_calculation="bgd")
equilibrium = eq_calc.fit_calculate(env, max_iter=1000, tol=1e-8)

print(f"Equilibrium concentrations: {equilibrium}")
```

## Documentation Structure

This documentation is organized into several sections:

1. **Getting Started** - Installation and basic usage
2. **Core Components** - Detailed documentation of each class
3. **Advanced Topics** - Advanced features and techniques
4. **Examples** - Practical examples and use cases
5. **API Reference** - Complete API documentation
6. **Additional Resources** - Troubleshooting, FAQ, and more

## Requirements

- Python 3.7 or higher
- NumPy
- Matplotlib (for plotting features)

## Support

For questions, issues, or contributions:

- Check the [FAQ](faq.md) for common questions
- Review [Troubleshooting](troubleshooting.md) for solutions to common issues
- See [Contributing](contributing.md) for how to contribute

## License

This project is licensed under the MIT License - see the [LICENSE](../LICENSE) file for details.

## Version

Current version: 0.1.0

---

**Next Steps**: Start with the [Installation Guide](installation.md) or jump to the [Quick Start Tutorial](quickstart.md) to begin using ChemCalc.
