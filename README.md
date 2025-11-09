# ChemCalc

A Python library for chemical reaction simulation, including kinetic modeling and thermodynamic equilibrium calculations.

## Features

- **Chemical Compound Representation**: Create and manage chemical compounds with formulas, phases, and physical properties
- **Reaction Definition**: Define chemical reactions with reactants, products, stoichiometric coefficients, and rate constants
- **Kinetic Simulation**: Simulate time-dependent concentration changes using numerical integration
- **Equilibrium Calculations**: Calculate equilibrium concentrations using multiple optimization algorithms:
  - Batch Gradient Descent (BGD)
  - Stochastic Gradient Descent (SGD)
  - Newton's Method
- **Phase Support**: Handle different phases (solid, liquid, gas, aqueous) with temperature-dependent phase transitions
- **Unicode Formula Display**: Automatic conversion to Unicode subscripts and superscripts for chemical formulas
- **Visualization**: Interactive and static plotting capabilities for kinetic simulations

## Installation

### Requirements

- Python 3.7+
- NumPy
- Matplotlib (for plotting features)

### Install Dependencies

```bash
pip install numpy matplotlib
```

### Installation from Source

```bash
git clone <repository-url>
cd ChemCalc
pip install -e .
```

**Important:** After installation, you can import the package directly using:

```python
from ChemCalc import Compound, Reaction, Enviroment
from ChemCalc.Kinetic import KineticalCalculator
from ChemCalc.Thermodynamic import EquilibriumCalculator
```

The package must be installed (using `pip install -e .`) for these imports to work. Without installation, you would need to use `from src.ChemCalc import ...` instead.

## Quick Start

### Basic Usage

```python
from ChemCalc import Compound, Reaction, Enviroment
from ChemCalc.Kinetic import KineticalCalculator
from ChemCalc.Thermodynamic import EquilibriumCalculator

# Create compounds
A = Compound("A")
B = Compound("B")

# Create a reaction: A ⇌ B
rxn = Reaction.from_string_simple_syntax(
    "A > B",
    concentrations=[1.0, 0.0],  # [A_initial, B_initial]
    K=2.0,  # Equilibrium constant
    kf=0.5,  # Forward rate constant
    kb=0.25  # Backward rate constant
)

# Create environment
env = Enviroment(rxn, T=298)  # Temperature in Kelvin

# Kinetic simulation
kc = KineticalCalculator(accuracy=1e-3)
kc.fit(env)
results = kc.calculate(time=10.0, plot=False)

# Equilibrium calculation
eq_calc = EquilibriumCalculator(method_of_calculation="bgd")
eq_calc.fit(env)
equilibrium = eq_calc.calculate(max_iter=1000, tol=1e-8)
```

## Core Components

### Compound

Represents a chemical compound with formula, phase information, and physical properties.

```python
# Simple compound
water = Compound("H2O")

# Compound with phase information
co2_gas = Compound("CO2", phase_point_list=[{"phase": "g", "temperature": 298}])

# Compound with melting/boiling points
ethanol = Compound("C2H5OH", mp=-114, bp=78)

# Disable Unicode formatting
simple = Compound("H2O", scription=False)
```

**Attributes:**

- `formula`: Chemical formula string
- `unicode_formula`: Unicode representation with subscripts/superscripts
- `phase_point_list`: List of phase data points
- `mp`: Melting point
- `bp`: Boiling point

**Methods:**

- `phase(temperature)`: Determine phase at given temperature

### Reaction

Represents a chemical reaction with reactants, products, and kinetic/thermodynamic parameters.

#### Creating Reactions

**From Simple Syntax:**

```python
# Simple: A ⇌ B
rxn1 = Reaction.from_string_simple_syntax("A > B", concentrations=[1.0, 0.0])

# With stoichiometry: 2A + B ⇌ 3C
rxn2 = Reaction.from_string_simple_syntax("2A + B > 3C", concentrations=[1.0, 1.0, 0.0])

# With phases: A.g + B.l ⇌ C.aq
rxn3 = Reaction.from_string_simple_syntax("A.g + B.l > C.aq")

# With rate dependencies: A2 + B1 > C1
rxn4 = Reaction.from_string_simple_syntax("A2 + B1 > C1")
```

**From Complex Syntax:**

```python
# Complex syntax supports more flexible compound names
rxn = Reaction.from_string_complex_syntax(
    "2_Fe(CN)6_-3 & Ce+2 > 2_Fe(CN)6_-4 & Ce+3",
    concentrations=[1.0, 1.0, 0.0, 0.0],
    K=1e5,
    kf=0.1,
    kb=1e-6
)
```

**Direct Initialization:**

```python
A = Compound("A")
B = Compound("B")

reactants = [{"stoichiometric_coefficient": 1, "compound": A, "rate_dependency": 1}]
products = [{"stoichiometric_coefficient": 1, "compound": B, "rate_dependency": 1}]

rxn = Reaction(
    reactants,
    products,
    [1.0],  # Reactant concentrations
    [0.0],  # Product concentrations
    K=2.0,
    kf=0.5,
    kb=0.25,
    T=298
)
```

**Parameters:**

- `reactants`: List of reactant dictionaries
- `products`: List of product dictionaries
- `K`: Equilibrium constant
- `kf`: Forward rate constant
- `kb`: Backward rate constant
- `T`: Temperature (Kelvin)

### Enviroment

Manages multiple reactions and compounds in a chemical system.

```python
# Single reaction
env = Enviroment(rxn1, T=298)

# Multiple reactions
env = Enviroment(rxn1, rxn2, rxn3, T=298)

# Set concentrations
env.concentrations = [1.0, 0.5, 0.0, 0.0]

# Add reactions
env.add(new_reaction)
env += another_reaction
```

**Key Properties:**

- `compounds`: List of all unique compounds
- `reactions`: List of reactions
- `concentrations`: Current concentrations
- `stoichiometric_coefficient_array`: Stoichiometric matrix
- `rate_constants_array`: Rate constants matrix

### KineticalCalculator

Simulates chemical reaction kinetics over time.

```python
# Initialize with accuracy (time step)
kc = KineticalCalculator(accuracy=1e-3)

# Fit to environment
kc.fit(env)

# Calculate concentrations over time
results = kc.calculate(
    time=10.0,  # Total simulation time
    checkpoint_time=[1.0, 5.0, 10.0],  # Optional: specific time points
    plot=False  # or "interactive" or "save"
)

# Or fit and calculate in one step
results = kc.fit_calculate(env, time=10.0, plot="interactive")
```

**Plotting Options:**

- `plot=False`: No plotting
- `plot="interactive"`: Interactive matplotlib plot (type 'exit' to close)
- `plot="save"`: Save plot to file (use `directory` parameter)

### EquilibriumCalculator

Calculates equilibrium concentrations using numerical optimization.

```python
# Initialize with method
eq_calc = EquilibriumCalculator(method_of_calculation="bgd")  # or "sgd" or "newton"

# Fit to environment
eq_calc.fit(env)

# Calculate equilibrium
equilibrium = eq_calc.calculate(
    max_iter=5000,
    learning_rate=0.1,
    tol=1e-8,
    backtrack_beta=0.5,
    min_concentration=1e-12
)

# Or fit and calculate in one step
equilibrium = eq_calc.fit_calculate(
    env,
    max_iter=1000,
    tol=1e-8
)
```

**Optimization Methods:**

- `"bgd"`: Batch Gradient Descent (default) - processes all reactions simultaneously
- `"sgd"`: Stochastic Gradient Descent - processes reactions in random order
- `"newton"`: Newton's Method - uses second-order information for faster convergence

**Parameters:**

- `max_iter`: Maximum iterations (default: 5000)
- `learning_rate`: Step size for gradient updates (default: 0.1)
- `tol`: Convergence tolerance (default: 1e-8)
- `backtrack_beta`: Backtracking line search parameter (default: 0.5)
- `min_concentration`: Minimum concentration threshold (default: 1e-12)

## Examples

### Example 1: Simple Reversible Reaction

```python
from ChemCalc import Compound, Reaction, Enviroment
from ChemCalc.Kinetic import KineticalCalculator

# Create reaction: A ⇌ B
rxn = Reaction.from_string_simple_syntax(
    "A > B",
    concentrations=[1.0, 0.0],
    K=2.0,
    kf=0.5,
    kb=0.25
)

env = Enviroment(rxn, T=298)

# Kinetic simulation
kc = KineticalCalculator(accuracy=0.01)
kc.fit(env)
results = kc.calculate(time=10.0, plot="interactive")
```

### Example 2: Multiple Reactions

```python
# Reaction 1: A ⇌ B
rxn1 = Reaction.from_string_simple_syntax("A > B", [1.0, 0.0], K=2.0, kf=0.5, kb=0.25)

# Reaction 2: B ⇌ C
rxn2 = Reaction.from_string_simple_syntax("B > C", [0.0, 0.0], K=1.5, kf=0.3, kb=0.2)

env = Enviroment(rxn1, rxn2, T=298)
env.concentrations = [1.0, 0.0, 0.0]

kc = KineticalCalculator(accuracy=0.01)
kc.fit(env)
results = kc.calculate(time=20.0, checkpoint_time=[5.0, 10.0, 15.0, 20.0])
```

### Example 3: Equilibrium Calculation

```python
from ChemCalc.Thermodynamic import EquilibriumCalculator

# A + 2B ⇌ C
rxn = Reaction.from_string_simple_syntax(
    "A + 2B > C",
    concentrations=[1.0, 2.0, 0.0],
    K=10.0,
    kf=0.5,
    kb=0.05
)

env = Enviroment(rxn, T=298)

# Calculate equilibrium using Newton's method
eq_calc = EquilibriumCalculator(method_of_calculation="newton")
equilibrium = eq_calc.fit_calculate(env, max_iter=100, tol=1e-10)

print(f"Equilibrium concentrations: {equilibrium}")
```

### Example 4: Phase-Dependent Reactions

```python
# Create compounds with specific phases
A_gas = Compound("A", phase_point_list=[{"phase": "g", "temperature": 298}])
B_liquid = Compound("B", phase_point_list=[{"phase": "l", "temperature": 298}])
C_aq = Compound("C", phase_point_list=[{"phase": "aq", "temperature": 298}])

# Reaction with phase annotations
rxn = Reaction.from_string_simple_syntax(
    "A.g + B.l > C.aq",
    concentrations=[1.0, 1.0, 0.0],
    K=5.0
)

env = Enviroment(rxn, T=298)
```

## Testing

Run the test suite using pytest:

```bash
pytest tests/
```

Run specific test files:

```bash
pytest tests/test_general.py
pytest tests/test_kinetic.py
pytest tests/test_thermodynamic.py
```

## Project Structure

```
ChemCalc/
├── src/                          # Source code directory
│   └── ChemCalc/                 # Main package
│       ├── __init__.py           # Package initialization (exports core classes)
│       ├── _general.py           # Core classes: Compound, Reaction, Enviroment
│       ├── Kinetic.py            # KineticalCalculator class for kinetic simulations
│       └── Thermodynamic.py      # EquilibriumCalculator class for equilibrium calculations
│
├── tests/                        # Test suite
│   ├── __init__.py               # Test package initialization
│   ├── test_general.py           # Tests for Compound, Reaction, Enviroment
│   ├── test_kinetic.py           # Tests for KineticalCalculator
│   └── test_thermodynamic.py     # Tests for EquilibriumCalculator
│
├── docs/                         # Documentation
│   └── index.md                  # Documentation index
│
├── LICENSE                       # MIT License
├── README.md                     # This file - project documentation
├── requirements.txt              # Python dependencies
├── setup.py                      # Package setup configuration
├── pyproject.toml                # Modern Python packaging configuration
└── pytest.ini                    # Pytest configuration
```

## Key Features in Detail

### 1. Chemical Compound Representation

The `Compound` class provides a flexible way to represent chemical compounds with various properties:

**Basic Usage:**

```python
# Simple compound
water = Compound("H2O")

# With phase information at specific temperature
co2 = Compound("CO2", phase_point_list=[{"phase": "g", "temperature": 298}])

# With melting and boiling points
ethanol = Compound("C2H5OH", mp=-114, bp=78)

# Multiple phase points
compound = Compound("H2O", phase_point_list=[
    {"phase": "s", "temperature": 273},
    {"phase": "l", "temperature": 298},
    {"phase": "g", "temperature": 373}
])
```

**Phase Determination:**
The library automatically determines the phase of a compound at a given temperature using:

1. Explicit phase points in `phase_point_list`
2. Melting and boiling points (`mp` and `bp`)
3. Priority: phase points override mp/bp logic

```python
water = Compound("H2O", mp=0, bp=100)
print(water.phase(-5))   # "s" (solid)
print(water.phase(50))   # "l" (liquid)
print(water.phase(120))  # "g" (gas)
```

### 2. Unicode Formula Formatting

Chemical formulas are automatically converted to Unicode with proper subscripts and superscripts for beautiful display:

```python
# Subscripts for molecular formulas
compound = Compound("H2O")
print(compound.unicode_formula)  # H₂O

compound = Compound("CO2")
print(compound.unicode_formula)  # CO₂

# Superscripts for ions
ion = Compound("Na+1")
print(ion.unicode_formula)  # Na⁺¹

ion = Compound("SO4-2")
print(ion.unicode_formula)  # SO₄⁻²

# Complex formulas
complex_ion = Compound("Fe(CN)6-3")
print(complex_ion.unicode_formula)  # Fe(CN)₆⁻³
```

**Disable Formatting:**

```python
compound = Compound("H2O", scription=False)
print(compound.unicode_formula)  # H2O (plain text)
```

### 3. Reaction Definition and Syntax

ChemCalc supports multiple ways to define chemical reactions, from simple to complex:

#### Simple Syntax

The simple syntax is intuitive and perfect for most use cases:

**Format:** `"reactants > products"`

**Features:**

- `+` separator for multiple compounds
- Prefix number for stoichiometric coefficient: `2A` means 2 moles of A
- Suffix number for rate dependency: `A2` means rate depends on [A]²
- Phase suffix: `.s`, `.l`, `.g`, or `.aq`
- Simple compound names (alphabetic only)

**Examples:**

```python
# Simple reversible reaction
rxn1 = Reaction.from_string_simple_syntax("A > B")

# With stoichiometry
rxn2 = Reaction.from_string_simple_syntax("2A + B > 3C")

# With phases
rxn3 = Reaction.from_string_simple_syntax("A.g + B.l > C.aq")

# With rate dependencies
rxn4 = Reaction.from_string_simple_syntax("A2 + B1 > C1")

# Combined: stoichiometry, phases, and rate dependencies
rxn5 = Reaction.from_string_simple_syntax("2A.g2 + B.l1 > 3C.aq1")
```

#### Complex Syntax

The complex syntax provides full control and supports complex compound names:

**Format:** `"reactants > products"` with `&` separator

**Structure:** `stoichiometric_coefficient_compound_rate_dependency`

**Features:**

- `&` separator for multiple compounds
- Full control over all parameters
- Supports complex compound names with parentheses, numbers, and charges
- Explicit specification of stoichiometric coefficients and rate dependencies

**Examples:**

```python
# Complex ions and compounds
rxn1 = Reaction.from_string_complex_syntax(
    "2_Fe(CN)6_-3 & Ce+2 > 2_Fe(CN)6_-4 & Ce+3"
)

# Mixed notation
rxn2 = Reaction.from_string_complex_syntax(
    "1_H2O_1 & 1_CO2_1 > 1_H2CO3_1"
)
```

### 4. Phase Handling

The library supports four physical phases with intelligent phase determination:

**Supported Phases:**

- `"s"`: Solid
- `"l"`: Liquid
- `"g"`: Gas
- `"aq"`: Aqueous

**Phase Specification Methods:**

1. **In Reaction Strings:**

   ```python
   rxn = Reaction.from_string_simple_syntax("A.g + B.l > C.aq")
   ```

2. **Via Phase Point List:**

   ```python
   compound = Compound("H2O", phase_point_list=[
       {"phase": "s", "temperature": 273},
       {"phase": "l", "temperature": 298}
   ])
   ```

3. **Using Melting/Boiling Points:**
   ```python
   compound = Compound("H2O", mp=0, bp=100)
   # Automatically determines phase based on temperature
   ```

**Important Note:** In equilibrium calculations, solid and liquid phases are excluded from the mass-action law. Only gas and aqueous phases participate in equilibrium expressions, which is physically correct as pure solids and liquids have unit activity.

### 5. Kinetic Simulation

The `KineticalCalculator` class provides powerful kinetic simulation capabilities:

**Key Features:**

- Numerical integration of reaction kinetics
- Configurable time step (accuracy parameter)
- Automatic concentration clamping (prevents negative values)
- Checkpoint recording at specific times
- Interactive and static plotting

**Usage:**

```python
kc = KineticalCalculator(accuracy=1e-3)  # Smaller = more accurate
kc.fit(env)
results = kc.calculate(
    time=10.0,
    checkpoint_time=[1.0, 5.0, 10.0],  # Record at these times
    plot="interactive"  # or "save" or False
)
```

**Plotting Options:**

- `plot=False`: No plotting, just return results
- `plot="interactive"`: Display interactive matplotlib plot (type 'exit' to close)
- `plot="save"`: Save plot to file (specify path with `directory` parameter)

### 6. Equilibrium Calculations

The `EquilibriumCalculator` class solves for equilibrium concentrations using advanced optimization algorithms:

**Three Optimization Methods:**

1. **Batch Gradient Descent (BGD)** - Default

   - Processes all reactions simultaneously
   - Stable and reliable
   - Good for most systems

2. **Stochastic Gradient Descent (SGD)**

   - Processes reactions in random order
   - Can be faster for large systems
   - Useful when reactions are loosely coupled

3. **Newton's Method**
   - Uses second-order information
   - Fastest convergence when near solution
   - Requires good initial guess

**Advanced Features:**

- Backtracking line search to ensure non-negative concentrations
- Configurable convergence tolerance
- Minimum concentration threshold to prevent numerical issues
- Automatic phase exclusion (solids/liquids excluded from equilibrium)

**Usage:**

```python
eq_calc = EquilibriumCalculator(method_of_calculation="newton")
equilibrium = eq_calc.fit_calculate(
    env,
    max_iter=1000,
    learning_rate=0.1,
    tol=1e-8,
    backtrack_beta=0.5,
    min_concentration=1e-12
)
```

### 7. Multi-Reaction Systems

ChemCalc excels at handling complex systems with multiple reactions:

**Features:**

- Automatic compound aggregation across reactions
- Shared compounds between reactions
- Consistent concentration tracking
- Mass conservation verification

**Example:**

```python
# Reaction 1: A ⇌ B
rxn1 = Reaction.from_string_simple_syntax("A > B", [1.0, 0.0], K=2.0)

# Reaction 2: B ⇌ C (B is shared)
rxn2 = Reaction.from_string_simple_syntax("B > C", [0.0, 0.0], K=1.5)

# Create environment with both reactions
env = Enviroment(rxn1, rxn2, T=298)
env.concentrations = [1.0, 0.0, 0.0]  # [A, B, C]

# Both kinetic and equilibrium calculations work seamlessly
```

### 8. Environment Management

The `Enviroment` class provides a unified interface for managing chemical systems:

**Key Capabilities:**

- Automatic compound deduplication
- Stoichiometric matrix generation
- Rate constant arrays
- Concentration management
- Reaction addition and modification

**Properties:**

- `compounds`: List of all unique compounds
- `reactions`: List of all reactions
- `concentrations`: Current concentrations array
- `stoichiometric_coefficient_array`: Matrix representation
- `rate_constants_array`: Rate constants for all reactions

**Dynamic Reaction Management:**

```python
env = Enviroment(rxn1, T=298)
env.add(rxn2)  # Add reaction
env += rxn3    # Or use += operator
```

## Limitations and Notes

1. **Numerical Stability**: Very small or very large equilibrium constants may require careful tuning of parameters
2. **Convergence**: Some systems may require adjustment of `max_iter`, `learning_rate`, or `tol` for convergence
3. **Phase Exclusion**: Solid and liquid phases are excluded from equilibrium expressions (only gas and aqueous)
4. **Mass Conservation**: The library assumes closed systems; mass conservation should be verified for your specific use case

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Authors

- **Mohammad Keifari** - _Initial work_ - [mohammadkeifari2007@gmail.com](mailto:mohammadkeifari2007@gmail.com)
