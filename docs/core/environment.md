# Environment

The `Enviroment` class (alias `Environment`) is the container for multi-reaction and redox systems.

## Construction

### Reactions only

```python
env = Enviroment(rxn1, rxn2, T=298, volume=1.0)
```

### Compounds only

```python
env = Enviroment.from_compounds({"Na+": 0.1, "Cl-": 0.1}, volume=1.0)
```

### Reactions + concentration overrides

```python
env = Enviroment(
    rxn1, rxn2,
    concentrations={"H+": 0.06, "Ca+2": 0.01},
    volume=0.1,
)
```

Dict values override concentrations summed from reactions.

### Half-reactions

Pass `HalfReaction` instances directly:

```python
env = Enviroment(kw, ka, hr, concentrations={"HSeO4-": 1.0}, buffer=["H+"])
env.set_electrode_potential(Eh=0.44)  # fixed E vs SHE
```

With two or more half-reactions and no fixed `electrode_Eh`, equilibrium solves a shared electrode potential jointly with concentrations.

## Accessors

| Property / method | Purpose |
|-------------------|---------|
| `env.compounds` | Ordered unique species |
| `env.concentrations` | Current molar concentrations |
| `env.concentrations_dict` | `{formula: c}` |
| `env.compound_labels` | Formulas aligned with vectors |
| `env.reactions` | List of `Reaction` objects |
| `env.half_reactions` | Registered half-reactions |
| `env.volume` | Litres |
| `env.last_equilibrium_result` | Last `EquilibriumResult` |

## Copy

`env.copy()` deep-copies the environment for titration workflows. Compounds are rewired by formula so each copy shares one compound object per formula internally.

## Buffering

Hold species fixed during equilibrium or kinetics:

```python
env = Enviroment(
    weak_acid, water_kw,
    concentrations={"H+": 1e-7},
    buffer=["H+"],
)
```

Explicit targets: `buffer={"H+": 1e-7}`. Re-snap after changing concentrations with `env.set_buffer(["H+"])`.

See [Activity and buffering](../guide/activity-and-buffering.md).

## Mixing

```python
envC = envA + envB
envD = 0.5 * envA + 4 * envB
envE = Enviroment.combine((0.5, envA), (4, envB))
```

See [Mixing and titration](../guide/mixing-and-titration.md).

## Solvers

```python
env.equilibrium(method="newton", tol=1e-10, return_details=True)
env.apply_equilibrium(method="newton", tol=1e-10)
env.kinetics(time=10.0, accuracy=1e-3, plot="save", directory="plot.png")
```

## Library environments

Factory functions return a **new** environment wired to shared library objects:

```python
from ChemCompute.environments import phosphoric_acid, silver_chloride

env = phosphoric_acid(concentrations={"H3PO4": 0.10})
```

See [Environments library](../libraries/environments.md).

## Related

- [Equilibrium](equilibrium.md)
- [Kinetics](kinetics.md)
- [Pourbaix guide](../guide/pourbaix.md)
