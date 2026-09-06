# Libraries overview

ChemCompute ships four related libraries of tabulated chemistry. Each name is a **function that returns a cached singleton** — call it every time you need the object:

```python
from ChemCompute.compounds import water
from ChemCompute.reactions import water_kw
from ChemCompute.half_reactions import iron_iii
from ChemCompute.environments import phosphoric_acid

assert water() is water()
env = phosphoric_acid(concentrations={"H3PO4": 0.10})
```

## Design rules

| Rule | Detail |
|------|--------|
| Concentrations | Start at **0** on library objects; set amounts on the `Enviroment` |
| K / E° | Tabulated at 25 °C unless noted; strong acids use `infinite_K` |
| Sharing | Same object every call — copy before mutating `spectrum`, `mp`, `bp` |
| Wiring | Environments use library reactions/compounds **as-is** (no deep copy) |
| Tokens | Use `.token` in strings to keep library identity and phase |

## Packages

| Package | Returns | Documentation |
|---------|---------|---------------|
| `ChemCompute.compounds` | `Compound` | [Compounds](compounds.md) |
| `ChemCompute.reactions` | `Reaction` | [Reactions](reactions.md) |
| `ChemCompute.half_reactions` | `HalfReaction` | [Half-reactions](half-reactions.md) |
| `ChemCompute.environments` | new `Enviroment` | [Environments](environments.md) |

## Listing everything

```python
from ChemCompute.compounds import all_compounds
from ChemCompute.reactions import all_reactions
from ChemCompute.half_reactions import all_half_reactions
from ChemCompute.environments import all_environments
```

## Formula lookup

```python
from ChemCompute.compounds import get

water = get("H2O")
```

Raises `KeyError` if the formula is missing or ambiguous.

## Related

- [Concepts](../concepts.md)
- [Examples](../examples/index.md)
