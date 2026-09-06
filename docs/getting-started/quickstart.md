# Quick start

Every calculation starts with compounds, reactions, and an environment.

## Single reaction

```python
from ChemCompute import Compound, Reaction, Enviroment

rxn = Reaction.from_string(
    "A > B",
    concentrations=[1.0, 0.0],
    K=2.0,
    kf=0.5,
    kb=0.25,
)

env = Enviroment(rxn, T=298)  # Kelvin
```

## Multiple reactions

Reactions in one environment share compounds automatically:

```python
rxn1 = Reaction.from_string("A > B", [1.0, 0.0], K=2.0, kf=0.5, kb=0.25)
rxn2 = Reaction.from_string("B > C", [0.0, 0.0], K=1.5, kf=0.3, kb=0.2)

env = Enviroment(rxn1, rxn2, T=298)
env.concentrations = [1.0, 0.0, 0.0]  # [A, B, C]
```

## Useful accessors

| Property / method | Purpose |
|-------------------|---------|
| `env.compounds` | Ordered list of unique species |
| `env.concentrations` | Current molar concentrations |
| `env.concentrations_dict` | `{formula: concentration}` mapping |
| `env.compound_labels` | Formula strings aligned with concentration vectors |

## Equilibrium

```python
equilibrium = env.equilibrium(method="newton", tol=1e-10)
print(equilibrium)  # concentration list

result = env.equilibrium(method="newton", tol=1e-10, return_details=True)
print(result.concentrations_dict)
print(result.criterion_met)
```

Or apply concentrations back onto the environment:

```python
env.apply_equilibrium(method="newton", tol=1e-10)
```

## Kinetics

```python
results = env.kinetics(time=10.0, accuracy=1e-3)
final = results[-1]
```

## Standalone reaction solvers

A single `Reaction` can be solved without building an environment manually:

```python
rxn = Reaction.from_string("A > B", [1.0, 0.0], K=2.0, kf=0.5, kb=0.25)
final = rxn.equilibrium(method="newton", tol=1e-10)
traces = rxn.kinetics(time=5.0, accuracy=1e-3)
```

## Next steps

- [Concepts](../concepts.md) — how Compound, Reaction, and Enviroment fit together
- [Reaction syntax](../guide/syntax.md) — string notation
- [Examples](../examples/index.md) — Selenium Pourbaix, titration, env 16, Jungle Model
