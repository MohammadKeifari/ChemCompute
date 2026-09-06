# Concepts

ChemCompute organizes calculations around three core types plus a marker for excess species.

## Compound

A **Compound** is a chemical species: formula, phase behaviour, optional melting/boiling points (kelvin), charge, and optional UV-Vis spectrum.

```python
from ChemCompute import Compound

water = Compound("H2O", phase_point_list=[{"phase": "l", "temperature": 298}])
h_plus = Compound("H+", phase_point_list=[{"phase": "aq", "temperature": 298}], charge=1)
```

Library compounds live in `ChemCompute.compounds` as cached getter functions — call `water()`, not `water`.

## Reaction

A **Reaction** holds reactants, products, stoichiometry, equilibrium constant K, and optional rate constants `kf` / `kb`. Build from explicit lists or from a string:

```python
from ChemCompute import Reaction

rxn = Reaction.from_string(
    "A > B",
    concentrations={"A": 1.0, "B": 0.0},
    K=2.0,
    kf=0.5,
    kb=0.25,
)
```

Reactions can carry thermodynamic parameters (`enthalpy`, `entropy`, activation energies) so K and rate constants update when temperature changes.

## Enviroment

An **Enviroment** (alias `Environment`) is the working system: one or more reactions and/or half-reactions, shared species list, concentrations, volume, temperature, optional activity model, buffer specification, and electrode potential for redox coupling.

```python
from ChemCompute import Enviroment

env = Enviroment(rxn1, rxn2, T=298, volume=1.0)
env.equilibrium(method="newton", tol=1e-10)
env.kinetics(time=10.0, accuracy=1e-3)
```

Half-reactions are passed directly to the constructor or registered later; they couple to equilibrium through Nernst electrode potential.

## XS — excess marker

`XS(amount)` marks a species whose concentration stays **fixed** during equilibrium or kinetics (solvent, excess solid). It does not multiply K by bulk molarity. Solids and liquids with phase `.s` / `.l` are also omitted from Q with activity 1.

```python
from ChemCompute import XS

kw = Reaction.from_string(
    "H2O.l > H+ & OH-",
    concentrations={"H2O": XS(55.5), "H+": 1e-7, "OH-": 1e-7},
    K=1e-14,
)
```

## Data flow

```
Compound ──► Reaction ──► Enviroment ──► equilibrium() / kinetics()
                │              │
                │              ├── Titration (mix + equilibrate)
                │              └── Pourbaix (pH–Eh predominance)
                └── HalfReaction (Nernst redox)
```

## Libraries

Tabulated chemistry (K, E°, spectra) lives in subpackages as **getter functions** that return shared singletons:

| Package | Returns | Concentrations |
|---------|---------|----------------|
| `ChemCompute.compounds` | `Compound` | 0 on the object; set on env |
| `ChemCompute.reactions` | `Reaction` | all 0 |
| `ChemCompute.half_reactions` | `HalfReaction` | all 0 |
| `ChemCompute.environments` | new `Enviroment` | 0 unless you pass `concentrations=` |

See [Libraries overview](libraries/overview.md).

## Next steps

- [Compound](core/compound.md)
- [Reaction](core/reaction.md)
- [Environment](core/environment.md)
