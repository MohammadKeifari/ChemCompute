# Reaction

The `Reaction` class defines stoichiometry, equilibrium constant K, rate constants, and optional thermodynamic parameters.

## From string

```python
from ChemCompute import Reaction, XS

rxn = Reaction.from_string(
    "A & B > C",
    concentrations={"A": 1.0, "B": 0.5, "C": 0.0},
    K=10.0,
    kf=0.5,
    kb=0.05,
)
```

See [Reaction syntax](../guide/syntax.md) for the full grammar.

## Concentrations

Pass a `{formula: amount}` dict (missing species default to 0) or a legacy ordered list (reactants then products). Use `XS(amount)` for excess species fixed during the solve.

## Thermodynamics and temperature

```python
rxn = Reaction.from_string(
    "A > B",
    K=2.0,
    kf=0.5,
    kb=0.25,
    enthalpy=-50000,
    entropy=-100,
    activation_energy_forward=50000,
    activation_energy_backward=100000,
    T=298,
)
rxn.T = 350  # K, kf, kb update via van't Hoff and Arrhenius
```

Disable automatic updates in an environment with `adjust_thermodynamics=False`.

## Strong / irreversible reactions

```python
Reaction.from_string("HCl.aq > H+ & Cl-", infinite_K=True)
```

`infinite_K=True` drives the reaction to completion analytically during equilibrium.

## Explicit construction

Build with reactant/product dict lists when string notation is insufficient:

```python
Reaction(
    reactants=[{"stoichiometric_coefficient": 1, "compound": a, "rate_dependency": 1}],
    products=[{"stoichiometric_coefficient": 1, "compound": b, "rate_dependency": 1}],
    reactants_concentration=[1.0],
    products_concentration=[0.0],
    K=2.0,
    kf=0.5,
    kb=0.25,
)
```

Omitted `rate_dependency` defaults to the stoichiometric coefficient.

## Standalone solvers

```python
final = rxn.equilibrium(method="newton", tol=1e-10)
traces = rxn.kinetics(time=5.0, accuracy=1e-3)
```

## Library reactions

```python
from ChemCompute.reactions import water_kw, acetic_acid

assert water_kw().K == 1e-14
env = Enviroment(water_kw(), acetic_acid(), concentrations={"CH3COOH": 0.10})
```

See [Reactions library](../libraries/reactions.md).

## Related

- [Equilibrium](equilibrium.md)
- [Kinetics](kinetics.md)
